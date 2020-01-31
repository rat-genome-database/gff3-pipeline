package edu.mcw.rgd.gff3;

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * @author mtutaj
 * @since 3/8/2017
 * convert fasta file for loading into JBrowse
 */
public class FastaPreprocessor {

    /**
     * convert contig headers into format suitable for loading by JBrowse in RGD
     * input header: (old, with gi)
     * &gt;gi|529417441|ref|NW_004936469.1| Ictidomys tridecemlineatus isolate #75 unplaced genomic scaffold, SpeTri2.0 scaffold00001, whole genome shotgun sequence
     * input header: (new, without gi)
     * &gt;ref|NW_004624730.1| Heterocephalus glaber isolate NMR 29 unplaced genomic scaffold, HetGla_female_1.0 scaffold00001, whole genome shotgun sequence
     * output header:
     * &gt;NW_004936469
     * @param args input file, output file
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {

        String inputFile = args[0];
        String outputFile = args[1];

        long bytesWritten = 0;

        BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile))));
        System.out.println("reading from to "+inputFile);

        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile))));
        String line;
        while( (line=reader.readLine())!=null ) {
            if( line.startsWith(">") ) {
                // try old format with gi
                int pos1 = line.indexOf("|ref|");
                if( pos1<0 ) {

                    // try new format without gi
                    pos1 = line.indexOf(">ref|");
                    if( pos1<0 ) {
                        throw new Exception("unexpected parsing error for line " + line);
                    }
                }
                pos1 += 5;
                int pos2 = line.indexOf('.', pos1);
                if( pos2<0 ) {
                    throw new Exception("unexpected parsing error for line "+line);
                }
                line = ">"+line.substring(pos1, pos2);
            }
            writer.write(line);
            writer.write("\n");

            bytesWritten += line.length()+1;
        }

        reader.close();
        writer.close();

        System.out.println(bytesWritten + " bytes written to "+outputFile);
    }
}
