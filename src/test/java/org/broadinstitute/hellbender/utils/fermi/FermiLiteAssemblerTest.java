package org.broadinstitute.hellbender.utils.fermi;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public final class FermiLiteAssemblerTest {
    private static final String expectedContig =
            "AATTTGCAAAAGGCCTAATAATCGGCAGAGTTGGTGCCTCTGGAGGTGAGTGTGAGGGGGATCTAATAAAAGAAGGTTTA"+
            "ACTGAAGTCTTTTAAGAAACAGGATTTTCACATCTAGTAATGTGACTCTTTTACTGAAATAACTAAAAATGCAGGAATCC"+
            "AGAGAGATAAGAAGAGTAATAAAAACAAGTGTCTATGAAAAGACACCTATAAGAATGTTCATAATAGTTTATTCATAATA"+
            "ACCCCAAACTAGAAGCAACCCAAATATCCACCAAAAAGAGAAGAGATACAAATTGTGGTACAGTCATATGCCACATAATG"+
            "ACATTTTGGTCAATGGCAGACTGCATATATGACAGTGTTCCCATGGGAGCTGAAAAATTCCTATCACCTAGTGACCTTGT"+
            "AGCCAATGTAAAGTCATAGCACAATGCATTACTCATGTGTTTGTGAGGAAGCTGGTGTAAACAAACCTACTACACTGCCA"+
            "GGCATATAAAATAATATAGCATATAGGATTACGTAGAGTACATAATACTTGGTAACATAATAAATGACTATGTTACTGAT"+
            "TTGTGTATTTAGTATACTTCTTATCATTATTGTAGAGTGTACTCCTACTTATTAAAAAAGTTAAGTGTAAAACAGCCTCA"+
            "GGCAGGTCCTTCAGGAGGTATTCCAGGAGAAAGCATTGTTATCACGGGATGACAGCTCCATGCATGTTATTGTCCCTGAA"+
            "GACCTTCCAGTGGGACAAGATATGGAGATAGGAAACAGTGACATTGATGATCCTGATCTTGTGTAGGCCTAGGCTAATGT"+
            "GTGTGTATGTCTTGGGTTTTAATAAAAAAGTTTAAAAAGTAAAAAAATAAAATAAATCTAGAAGTTTTTTAAATAAAAAA"+
            "AGCTTACAGAATAAAGATAAAAATAAAATATTTTTGTACAGCTATAAAGTATGTTTGTATTTTAAGCCAAGTGTTATTAA"+
            "AAGAGTCAAAAAGTTAAAAAGTTTATAATGTAAAAAAGTTACAGTACGCTAAGGTTAATTTATTATTGAAGAAAATTAAA"+
            "AAAAAATTTAGTGTAGCCAAAGTGTACAGTGTTTATAAAGTCCACAGTAGTGTACAGTAATATCCTAGGCCTTCACATTC"+
            "ACTCACCACTCATTCACCACTCACTCAGAGCAACTTCCAGTATTGCAAGCTCCATTCACGGTAAGTACCCTGTATAGGTG"+
            "TACCATTTAAAAAATCTTTTAAACCATAATTTTACTGTATCCTTTTTATCTTTATACATGTTTAGATACACAAATACTTA"+
            "ACCATTGTGCTGCAACTGCCTATAGCACAGTAACATGCTATACAGGTTTGTAGCCTAGGAGCAATGGGCTATACCATCTA"+
            "GGTTTGAGTAAGTACGCTATGACAGGGGTCCCCAACCCCCAGGCTGAAGACCAATACGGGGGTCTACAGCCTGTTAGGAA"+
            "AGGGGCTGCACAGCAGGAGTTGAGCAGCAGGCGAGTGAGCATTACTACTTGAGCTCCGCCTCCTGTCAGATCAACAGTGG"+
            "CATTAAATCCTCACAGGAGCATGAATCCTATTGTAAACTGCGCATGCAAGGGATCTAGGTTGTGTGCTCCTTATGAGAAT"+
            "CTAATTAATGGCTGATGATCTGAGATGGAACAGTTTCATCCTGAAACCATCCCCCATGCCACTGGTCCATATATATATAT"+
            "ATATATATATATATATATATATATATATAT";

    private static List<FakeRead> genReads( final String seq, final int coverage, final int readLen ) {
        final int seqLen = seq.length();
        final int nReads = coverage * seqLen / readLen + 1;
        final List<FakeRead> reads = new ArrayList<>(nReads);
        new Random(0)
                .ints(nReads, 0, seq.length()-readLen)
                .mapToObj(idx -> seq.substring(idx,idx+readLen))
                .map(String::getBytes)
                .map(FakeRead::new)
                .forEach(reads::add);
        // make sure there's adequate coverage at the ends
        reads.add(new FakeRead(seq.substring(0, readLen).getBytes()));
        reads.add(new FakeRead(seq.substring(1, readLen+1).getBytes()));
        reads.add(new FakeRead(seq.substring(seqLen-readLen, seqLen).getBytes()));
        reads.add(new FakeRead(seq.substring(seqLen-readLen-1, seqLen-1).getBytes()));
        return reads;
    }

    @Test
    public void testOptsSize() {
        try ( final FermiLiteAssembler assembler = new FermiLiteAssembler() ) {
            Assert.assertEquals(assembler.getOptsSize(), assembler.getExpectedOptsSize());
        }
    }

    private static final class FakeRead implements FermiLiteAssembler.BasesAndQuals {
        final byte[] seq;
        final byte[] quals;

        public FakeRead( final byte[] seq ) {
            this.seq = seq;
            quals = new byte[seq.length];
            Arrays.fill(quals, (byte)30);
        }

        public byte[] getBases() { return seq; }
        public byte[] getQuals() { return quals; }
    }

    @Test
    public void testSingleContig() {
        final int readLen = 151;
        final FermiLiteAssembly assembly =
                new FermiLiteAssembler().createAssembly(genReads(expectedContig, 30, readLen));
        Assert.assertEquals(assembly.getNContigs(), 1);
        final FermiLiteAssembly.Contig contig = assembly.getContig(0);
        Assert.assertEquals(contig.getSequence(), expectedContig.getBytes());
        final List<FermiLiteAssembly.Connection> connections = contig.getConnections();
        Assert.assertEquals(connections.size(), 1);
        final FermiLiteAssembly.Connection connection = connections.get(0);
        Assert.assertEquals(connection.getTarget(), contig);
        Assert.assertEquals(connection.getOverlapLen(), 42);
        Assert.assertEquals(connection.isRC(), false);
        Assert.assertEquals(connection.isTargetRC(), true);
    }
}
