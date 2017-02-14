package org.broadinstitute.hellbender.utils.fermi;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

public final class FermiLiteAssembler {
    private static volatile boolean nativeLibLoaded = false;

    public FermiLiteAssembler() {
        loadNativeLibrary();
    }

    public interface BasesAndQuals {
        byte[] getBases();
        byte[] getQuals();
    }

    public FermiLiteAssembly createAssembly( final Iterable<BasesAndQuals> basesAndQuals ) {
        return createAssembly(basesAndQuals, bAndQ -> bAndQ);
    }

    public <T> FermiLiteAssembly createAssembly( final Iterable<T> reads, final Function<T,BasesAndQuals> func ) {
        final ByteBuffer assemblyData = createAssemblyData(makeReadData(reads, func));
        if ( assemblyData == null ) throw new IllegalStateException("Unable to create assembly.");
        try {
            return interpretAssemblyData(assemblyData);
        } finally {
            destroyAssemblyData(assemblyData);
        }
    }

    private static void loadNativeLibrary() {
        if ( !nativeLibLoaded ) {
            synchronized(FermiLiteAssembler.class) {
                if ( !nativeLibLoaded ) {
                    final String libNameOverride = System.getProperty("LIBFML_PATH");
                    if ( libNameOverride != null ) {
                        System.load(libNameOverride);
                    }
                    else {
                        final String osName = System.getProperty("os.name", "unknown").toUpperCase();
                        final String osArch = System.getProperty("os.arch");
                        final String libName;
                        if ( !"x86_64".equals(osArch) && !"amd64".equals(osArch) ) libName = null;
                        else if ( osName.startsWith("MAC") ) libName = "/libfml.Darwin.dylib";
                        else if ( osName.startsWith("LINUX") ) libName = "/libfml.Linux.so";
                        else libName = null;
                        if ( libName == null ) {
                            throw new IllegalStateException("We have a JNI binding for fermi-lite only for x86-64 Linux and Mac.");
                        }
                        try ( final InputStream is = FermiLiteAssembler.class.getResourceAsStream(libName) ) {
                            if ( is == null ) {
                                throw new IllegalStateException("Can't find resource "+libName);
                            }
                            final File tmpFile = File.createTempFile("libfml.",".jnilib");
                            tmpFile.deleteOnExit();
                            Files.copy(is, tmpFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                            System.load(tmpFile.getPath());
                        }
                        catch (IOException ioe ) {
                            throw new IllegalStateException("Misconfiguration: Unable to load fermi-lite native library "+libName, ioe);
                        }
                    }
                    nativeLibLoaded = true;
                }
            }
        }
    }

    // Writes the number of reads, and then seqs and quals for each into a ByteBuffer.
    private static <T> ByteBuffer makeReadData( final Iterable<T> reads, final Function<T,BasesAndQuals> func ) {
        // capacity calculation:
        //   4 bytes to give the number of reads, plus
        //   for each read,
        //     for each base call, we need two bytes (one for the base, one for the qual)
        //     additionally, we need two bytes for the terminating nulls after the calls
        //     i.e., 2*(length+1)
        int nReads = 0;
        int capacity = 4;
        for ( final T read : reads ) {
            nReads += 1;
            capacity += 2*(func.apply(read).getBases().length + 1);
        }
        final ByteBuffer readData = ByteBuffer.allocateDirect(capacity);
        readData.order(ByteOrder.nativeOrder());
        readData.putInt(nReads); // array length
        for ( final T read : reads ) {
            final BasesAndQuals bAndQ = func.apply(read);
            readData.put(bAndQ.getBases()).put((byte)0);
            readData.put(bAndQ.getQuals()).put((byte)0);
        }
        readData.flip();
        return readData;
    }

    // expects a direct ByteBuffer containing:
    //  the number of contigs (4-byte int)
    //  the offset to the beginning of a byte pool containing sequence and per-base coverage bytes (4 byte int)
    //  N.B.: the sequence and per-base coverage bytes are NOT null terminated.
    //  for each contig, an fml_utg_t structure minus the seq and cov pointers:
    //    the length of the sequence (and per-base support) data (4-byte int)
    //    the number of supporting reads (4-byte int)
    //    the number of connections (4-byte int)
    //    a variable number (given by # of connections, above) of fml_ovlp_t's (8 bytes each)
    //  a byte pool containing the seq and cov data
    private static FermiLiteAssembly interpretAssemblyData( final ByteBuffer assemblyData ) {
        assemblyData.order(ByteOrder.nativeOrder()).position(0).limit(assemblyData.capacity());

        // make the contigs
        final int nContigs = assemblyData.getInt();
        int seqOffset = assemblyData.getInt();
        final List<FermiLiteAssembly.Contig> contigs = new ArrayList<>(nContigs);
        for ( int idx = 0; idx != nContigs; ++idx ) {
            final int seqLen = assemblyData.getInt();
            final int nSupportingReads = assemblyData.getInt();
            final int nConnections = assemblyData.getInt();
            final int mark = assemblyData.position()+8*nConnections; // sizeof(fml_ovlp_t) is 8
            assemblyData.position(seqOffset);
            final byte[] seq = new byte[seqLen];
            assemblyData.get(seq);
            final byte[] coverage = new byte[seqLen];
            assemblyData.get(coverage);
            contigs.add(new FermiLiteAssembly.Contig(seq,coverage,nSupportingReads));
            assemblyData.position(mark);
            seqOffset += 2*seqLen;
        }
        // connect the contigs
        assemblyData.position(8); // skip past nContigs and seqOffset
        for ( int idx = 0; idx != nContigs; ++idx ) {
            final FermiLiteAssembly.Contig contig = contigs.get(idx);
            assemblyData.getInt(); // skip seqLen
            assemblyData.getInt(); // skip # of supporting reads
            int nConnections = assemblyData.getInt();
            final List<FermiLiteAssembly.Connection> connections = new ArrayList<>(nConnections);
            while ( nConnections-- > 0 ) {
                int overlapLen = assemblyData.getInt();
                final boolean isRC = overlapLen < 0;
                overlapLen = overlapLen & Integer.MAX_VALUE; // turn off highest bit (which is now in isRC)
                int contigId = assemblyData.getInt();
                final boolean isTargetRC = contigId < 0;
                contigId = contigId & Integer.MAX_VALUE; // turn off highest bit (which is now in isTargetRC)
                connections.add(new FermiLiteAssembly.Connection(contigs.get(contigId), overlapLen, isRC, isTargetRC));
            }
            contig.setConnections(connections);
        }

        return new FermiLiteAssembly(contigs);
    }

    // these should be called in succession by the same thread
    private static native ByteBuffer createAssemblyData( final ByteBuffer readData );
    private static native void destroyAssemblyData( final ByteBuffer assemblyData );
}
