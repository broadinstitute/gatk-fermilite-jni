package org.broadinstitute.hellbender.utils.fermi;

import java.io.*;
import java.util.Arrays;
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

    public int computeN50() {
        final int nContigs = contigs.size();
        if ( nContigs < 1 ) return 0;
        if ( nContigs == 1 ) return contigs.get(0).getSequence().length;

        final int[] lengths = new int[nContigs];
        for ( int idx = 0; idx != nContigs; ++idx ) {
            lengths[idx] = contigs.get(idx).getSequence().length;
        }
        Arrays.sort(lengths);
        final int totalSize = Arrays.stream(lengths).sum();
        int runningTotal = 0;
        for ( int idx = nContigs - 1; idx >= 0; --idx ) {
            runningTotal += 2 * lengths[idx];
            if ( runningTotal >= totalSize ) return lengths[idx];
        }
        throw new ArithmeticException("impossible situation -- sum of array greater than twice the sum of each element");
    }

    /** a sequence of bases, coverage data, and connections to other contigs */
    public static final class Contig {
        private final byte[] sequence;
        private final byte[] perBaseCoverage;
        private final int nSupportingReads;
        private List<Connection> connections;

        public Contig( final byte[] sequence, final byte[] perBaseCoverage, final int nSupportingReads ) {
            this.sequence = sequence;
            this.perBaseCoverage = perBaseCoverage;
            this.nSupportingReads = nSupportingReads;
            this.connections = Collections.emptyList();
        }

        public byte[] getSequence() { return sequence; }
        public byte[] getPerBaseCoverage() { return perBaseCoverage; }
        public int getNSupportingReads() { return nSupportingReads; }
        public List<Connection> getConnections() { return connections; }

        public void setConnections( final List<Connection> connections ) {
            this.connections = Collections.unmodifiableList(connections);
        }

        public Connection getSolePredecessor() {
            return getSingletonConnection(true);
        }

        public Connection getSoleSuccessor() {
            return getSingletonConnection(false);
        }

        public Connection getSingletonConnection( final boolean isRC ) {
            Connection singleton = null;
            for ( Connection conn : connections ) {
                if ( conn.isRC() == isRC ) {
                    if ( singleton != null ) return null; // found multiple connections, return null
                    singleton = conn;
                }
            }
            return singleton;
        }
    }

    /** a connection between contigs */
    public static final class Connection {
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

        /** call with the contig among whose connections this one appears */
        public Connection rcConnection( final Contig contig ) {
            return new Connection(contig, overlapLen, !isTargetRC, !isRC);
        }
    }

    /** Dump a fermi-lite native format description of the assembly. */
    public void writeGFA( final Writer writer ) throws IOException {
        final HashMap<Contig, Integer> idMap = new HashMap<>((int)((contigs.size()*4L)/3) + 1);
        int id = 0;
        for (final Contig contig : contigs) {
            idMap.put(contig, id++);
        }
        writer.write("H\tVN:Z:1.0\n");
        for ( final Contig contig : contigs ) {
            final int contigId = idMap.get(contig);
            writer.write("S\ttig" + contigId + "\t" + new String(contig.getSequence()) +
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
