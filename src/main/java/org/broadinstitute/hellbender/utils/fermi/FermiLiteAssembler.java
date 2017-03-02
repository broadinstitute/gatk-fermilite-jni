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

/**
 * Class that allows you to exercise Heng Li's fermi-lite assembler.
 * This is backed by JNI native code, which could fail to load.
 * This class is not thread-safe, but it's very light-weight:  just use a separate instance in each thread.
 */
public final class FermiLiteAssembler implements AutoCloseable {
    private static volatile boolean nativeLibLoaded = false;
    private ByteBuffer opts;

    public FermiLiteAssembler() {
        loadNativeLibrary();
        opts = createDefaultOptions();
        opts.order(ByteOrder.nativeOrder()).position(0).limit(opts.capacity());
    }

    public boolean isOpen() { return opts != null; }

    @Override
    public void close() { checkOpen(); destroyByteBuffer(opts); opts = null; }

    public interface BasesAndQuals {
        byte[] getBases();
        byte[] getQuals();
    }

    /** number of threads; don't use multi-threading for small data sets */
    public int getNThreads() { checkOpen(); return opts.getInt(0); }
    public void setNThreads( final int nThreads ) { checkOpen(); opts.putInt(0,nThreads); }
    /** k-mer length for error correction; 0 for auto estimate */
    public int getECKSize() { checkOpen(); return opts.getInt(4); }
    public void setECKSize( final int kSize ) { checkOpen(); opts.putInt(4,kSize); }

    /** both occ threshold in ec and tip threshold in cleaning lie in [min_cnt,max_cnt] */
    public int getMinCnt() { checkOpen(); return opts.getInt(8); }
    public void setMinCnt( final int minCnt ) { checkOpen(); opts.putInt(8,minCnt); }
    public int getMaxCnt() { checkOpen(); return opts.getInt(12); }
    public void setMaxCnt( final int maxCnt ) { checkOpen(); opts.putInt(12,maxCnt); }
    /** min overlap length during assembly */
    public int getMinAsmOverlap() { checkOpen(); return opts.getInt(16); }
    public void setMinAsmOverlap( final int minAsmOverlap ) { checkOpen(); opts.putInt(16,minAsmOverlap); }
    /** during assembly, don't explicitly merge an overlap if shorter than this value */
    public int getMinMergeLen() { checkOpen(); return opts.getInt(20); }
    public void setMinMergeLen( final int minMergeLen ) { checkOpen(); opts.putInt(20,minMergeLen); }

    /** graph cleaning options -- you'll have to look at the fermi lite code to try to figure out what these do */
    public int getCleaningFlag() { checkOpen(); return opts.getInt(24); }
    public void setCleaningFlag( final int cleaningFlag ) { checkOpen(); opts.putInt(24,cleaningFlag); }
    public int getCleaningMinOverlap() { checkOpen(); return opts.getInt(28); }
    public void setCleaningMinOverlap( final int minOverlap ) { checkOpen(); opts.putInt(28,minOverlap); }
    public int getCleaningELen() { checkOpen(); return opts.getInt(32); }
    public void setCleaningELen( final int eLen ) { checkOpen(); opts.putInt(32,eLen); }
    public int getCleaningMinEnsr() { checkOpen(); return opts.getInt(36); }
    public void setCleaningMinEnsr( final int minEnsr ) { checkOpen(); opts.putInt(36,minEnsr); }
    public int getCleaningMinInsr() { checkOpen(); return opts.getInt(40); }
    public void setCleaningMinInsr( final int minInsr ) { checkOpen(); opts.putInt(40,minInsr); }
    public int getCleaningMaxBDist() { checkOpen(); return opts.getInt(44); }
    public void setCleaningMaxBDist( final int maxBDist ) { checkOpen(); opts.putInt(44,maxBDist); }
    public int getCleaningMaxBDiff() { checkOpen(); return opts.getInt(48); }
    public void setCleaningMaxBDiff( final int maxBDiff ) { checkOpen(); opts.putInt(48,maxBDiff); }
    public int getCleaningMaxBVtx() { checkOpen(); return opts.getInt(52); }
    public void setCleaningMaxBVtx( final int maxBVtx ) { checkOpen(); opts.putInt(52,maxBVtx); }
    public int getCleaningMinMergeLen() { checkOpen(); return opts.getInt(56); }
    public void setCleaningMinMergeLen( final int minMergeLen ) { checkOpen(); opts.putInt(56,minMergeLen); }
    public int getCleaningTrimLen() { checkOpen(); return opts.getInt(60); }
    public void setCleaningTrimLen( final int trimLen ) { checkOpen(); opts.putInt(60,trimLen); }
    public int getCleaningTrimDepth() { checkOpen(); return opts.getInt(64); }
    public void setCleaningTrimDepth( final int trimDepth ) { checkOpen(); opts.putInt(64,trimDepth); }
    public float getCleaningDRatio1() { checkOpen(); return opts.getFloat(68); }
    public void setCleaningDRatio1( final float dRatio1 ) { checkOpen(); opts.putFloat(68,dRatio1); }
    public float getCleaningMaxBCov() { checkOpen(); return opts.getFloat(72); }
    public void setCleaningMaxBCov( final float maxBCov ) { checkOpen(); opts.putFloat(72,maxBCov); }
    public float getCleaningMaxBFrac() { checkOpen(); return opts.getFloat(76); }
    public void setCleaningMaxBFrac( final float maxBFrac ) { checkOpen(); opts.putFloat(76,maxBFrac); }
    public int getExpectedOptsSize() { return 80; }
    public int getOptsSize() { return opts.capacity(); }

    /**
     * Create an assembly from a collection of objects that implement BasesAndQuals.
     */
    public FermiLiteAssembly createAssembly( final Iterable<? extends BasesAndQuals> basesAndQuals ) {
        return createAssembly(basesAndQuals, bAndQ -> bAndQ);
    }

    /**
     * Create an assembly from a collection of objects that can be transformed (with a lambda) into BasesAndQuals.
     */
    public <T> FermiLiteAssembly createAssembly( final Iterable<T> reads, final Function<T,BasesAndQuals> func ) {
        checkOpen();
        final ByteBuffer assemblyData = createAssemblyData(opts,makeReadData(reads, func));
        if ( assemblyData == null ) throw new IllegalStateException("Unable to create assembly.");
        try {
            return interpretAssemblyData(assemblyData);
        } finally {
            destroyByteBuffer(assemblyData);
        }
    }

    public static String getFermiLiteVersion() {
        loadNativeLibrary();
        return getVersion();
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
                        if ( !"x86_64".equals(osArch) && !"amd64".equals(osArch) ) {
                            throw new IllegalStateException(
                                    "We have pre-built fermi-lite binaries only for x86_64 and amd64.  "+
                                    "Your os.arch is "+osArch+"."+
                                    "Set property LIBFML_PATH to point to a native library for your architecture.");
                        }
                        if ( osName.startsWith("MAC") ) libName = "/libfml.Darwin.dylib";
                        else if ( osName.startsWith("LINUX") ) libName = "/libfml.Linux.so";
                        else {
                            throw new IllegalStateException(
                                    "We have pre-built fermi-lite binaries only for Linux and Mac.  "+
                                    "Your os.name is "+osName+"."+
                                    "Set property LIBFML_PATH to point to a native library for your operating system.");
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

    private void checkOpen() {
        if ( opts == null ) {
            throw new IllegalStateException("The assembler has been closed.");
        }
    }

    private static native ByteBuffer createDefaultOptions();
    private static native ByteBuffer createAssemblyData( final ByteBuffer opts, final ByteBuffer readData );
    private static native void destroyByteBuffer( final ByteBuffer byteBuffer );
    private static native String getVersion();
}
