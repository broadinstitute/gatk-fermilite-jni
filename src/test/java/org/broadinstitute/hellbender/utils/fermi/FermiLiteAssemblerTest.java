package org.broadinstitute.hellbender.utils.fermi;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FermiLiteAssemblerTest {

    private final static class TestRead implements FermiLiteAssembler.BasesAndQuals {
        final byte[] seq;
        TestRead( final String seq ) { this.seq = seq.getBytes(); }
        @Override public byte[] getBases() { return seq; }
        @Override public byte[] getQuals() {
            final byte[] quals = new byte[seq.length];
            Arrays.fill(quals, (byte)30);
            return quals;
        }
    }

    @Test
    public void testFermiLiteAssembly() throws IOException {
        final String expectedContig =
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
                "CTAATTAATGGCTGATGATCTGAGATGGAACAGTTTCATCCTGAAACCATCCCCCATGCCACTGGTCCAT";

        final List<String> seqsList = new ArrayList<>(1200);
        String seq;
        try ( final BufferedReader rdr = new BufferedReader(new FileReader("src/test/resources/test.seq")) ) {
            while ( (seq = rdr.readLine()) != null ) {
                seqsList.add(seq);
            }
        }
        final FermiLiteAssembly assembly = new FermiLiteAssembler().createAssembly(seqsList, TestRead::new);
        Assert.assertEquals(assembly.getNContigs(), 1);
        Assert.assertEquals(assembly.getContig(0).getSequence(), expectedContig.getBytes());
        final byte[] expected = Files.readAllBytes(new File("src/test/resources/test.gfa").toPath());
        try ( final ByteArrayOutputStream os = new ByteArrayOutputStream() ) {
            assembly.writeGFA(os);
            Assert.assertEquals(os.toByteArray(), expected);
        }
    }
}
