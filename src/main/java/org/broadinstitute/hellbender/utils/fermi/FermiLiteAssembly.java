package org.broadinstitute.hellbender.utils.fermi;

import java.io.*;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/** an assembly is just a collection of contigs */
public final class FermiLiteAssembly {
    private final List<Contig> contigs;

    public FermiLiteAssembly( final List<Contig> contigs ) {
        this.contigs = Collections.unmodifiableList(contigs);
    }

    public int getNContigs() { return contigs.size(); }
    public Contig getContig( final int idx ) { return contigs.get(idx); }
    public List<Contig> getContigs() { return contigs; }

    /** a sequence of bases, coverage data, and connections to other contigs */
    public final static class Contig {
        private final byte[] sequence;
        private final byte[] perBaseCoverage;
        private final int nSupportingReads;
        private List<Connection> connections;

        public Contig( final byte[] sequence, final byte[] perBaseCoverage, final int nSupportingReads ) {
            this.sequence = sequence;
            this.perBaseCoverage = perBaseCoverage;
            this.nSupportingReads = nSupportingReads;
        }

        public byte[] getSequence() { return sequence; }
        public byte[] getPerBaseCoverage() { return perBaseCoverage; }
        public int getNSupportingReads() { return nSupportingReads; }
        public List<Connection> getConnections() { return connections; }

        public void setConnections( final List<Connection> connections ) {
            this.connections = Collections.unmodifiableList(connections);
        }
    }

    /** a connection between contigs */
    public final static class Connection {
        private final Contig target;      // contig that overlaps the one that possesses this connection
        private final int overlapLen;     // bases in common -- negative overlap lengths are legal, and represent gaps
        private final boolean isRC;       // if target is a predecessor (i.e., upstream of the 5' end of this one)
        private final boolean isTargetRC; // if connection is to RC of target contig

        public Connection( final Contig target, final int overlapLen, final boolean isRC, final boolean isTargetRC ) {
            this.target = target;
            this.overlapLen = overlapLen;
            this.isRC = isRC;
            this.isTargetRC = isTargetRC;
        }

        /** contig that overlaps the one that possesses this connection */
        public Contig getTarget() { return target; }
        /** bases in common -- negative overlap lengths are legal, and represent gaps */
        public int getOverlapLen() { return overlapLen; }
        /** if target is a predecessor (i.e., upstream of the 5' end of this one) */
        public boolean isRC() { return isRC; }
        /** if connection is to RC of target contig */
        public boolean isTargetRC() { return isTargetRC; }
    }

    /** Dump a fermi-lite native format description of the assembly. */
    public void writeGFA( final OutputStream os ) throws IOException {
        final HashMap<Contig, Integer> idMap = new HashMap<>((int)((contigs.size()*4L)/3) + 1);
        int id = 0;
        for (final Contig contig : contigs) {
            idMap.put(contig, id++);
        }
        try ( final Writer writer =
                      new OutputStreamWriter(new BufferedOutputStream(os)) ) {
            writer.write("H\tVN:Z:1.0\n");
            for ( final Contig contig : contigs ) {
                final int contigId = idMap.get(contig);
                writer.write("S\t" + contigId + "\t" + new String(contig.getSequence()) +
                        "\tLN:i:" + contig.getSequence().length + "\tRC:i:" + contig.getNSupportingReads() + "\n");
                for ( final Connection connection : contig.getConnections() ) {
                    final int targetId = idMap.get(connection.getTarget());
                    if ( contigId <= targetId ) {
                        final int overlapLen = connection.getOverlapLen();
                        writer.write("L\ttig" + contigId + "\t" + (connection.isRC() ? "-" : "+") +
                                "\ttig" + targetId + "\t" + (connection.isTargetRC() ? "-" : "+") + "\t" +
                                (overlapLen < 0 ? -overlapLen + "H" : overlapLen + "M") + "\n");
                    }
                }
            }
        }
    }
}
