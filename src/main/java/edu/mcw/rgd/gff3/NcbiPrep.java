package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.datamodel.Chromosome;
import edu.mcw.rgd.process.Utils;

import java.io.*;
import java.util.List;

/**
 * @author mtutaj
 * @since 01/28/2021
 * Prepares a fasta file downloaded from NCBI FTP site to be loaded by JBrowse scripts as an NCBI track
 * <ol>
 *  <li>Prepends all chromosome names with 'Chr'
 *  <li>Splits one file into multiple, one per chromosome
 * </ol>
 */
public class NcbiPrep {

    public static void main(String[] args) throws Exception {

        //prepUthGff3Files();

        int mapKey = 634;
        boolean isScaffoldAssembly = false;

        // for scaffold assemblies, prefix must be an empty string; for chromosome assemblies it must be 'Chr'
        final String chrPrefix = isScaffoldAssembly ? "" : "Chr";

        String fname = "/Users/mtutaj/Downloads/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna.gz";
        String outDir = "/data/ref/fasta/ROS_Cfam_1.0";

        MapDAO dao = new MapDAO();
        List<Chromosome> chromosomes = dao.getChromosomes(mapKey);

        // orig line:
        // >NC_051336.1 Rattus norvegicus strain BN/NHsdMcwi chromosome 1, mRatBN7.2, whole genome shotgun sequence
        // new line
        // Chr1
        String line;
        String chrFileName;
        BufferedWriter out = null;
        BufferedReader in = Utils.openReader(fname);
        while( (line=in.readLine())!=null ) {
            if( line.startsWith(">") ) {
                if( out!=null ) {
                    out.close();
                }
                // parse chr acc
                int spacePos = line.indexOf(' ');
                String chrAcc = line.substring(1, spacePos).trim();
                String chr = getChrName(chrAcc, chromosomes, chrPrefix);
                if( chr==null ) {
                    System.out.println("cannot match chromosome for line: "+line);
                }
                chrFileName = outDir+"/"+chr+".fna";
                System.out.println(chrFileName);
                out = Utils.openWriter(chrFileName);
                out.write(">"+chr+"\n");
            } else {
                out.write(line);
                out.write("\n");
            }
        }
        in.close();
        if( out!=null ) {
            out.close();
        }
        System.out.println("OK");
    }

    static String getChrName(String acc, List<Chromosome> chromosomes, String chrPrefix) {
        for( Chromosome ch: chromosomes ) {
            if( acc.equals(ch.getRefseqId()) ) {
                return chrPrefix+ch.getChromosome();
            }
        }
        return null;
    }

    static void prepUthFastaFiles() throws Exception {

        String fname =   "/Users/mtutaj/Downloads/rat4/shr.fasta";
        String outName = "/Users/mtutaj/Downloads/rat4/shr2.fa";

        // orig lines:
        //   >chrY Rattus norvegicus strain WKY/Bbb RGD_1581635 chromosome Y, whole genome shotgun sequence
        //   >JALPNS010000023.1 Rattus norvegicus strain WKY/Bbb RGD_1581635 unassigned_1, whole genome shotgun sequence
        // new line
        // >ChrY
        // >JAL...  lines are skipped
        BufferedReader in = Utils.openReader(fname);
        BufferedWriter out = Utils.openWriter(outName);
        String line;
        boolean scaffold = false;
        int linesRead = 0;
        int linesWritten = 0;

        while( (line=in.readLine())!=null ) {

            // header line
            linesRead++;
            if( line.startsWith(">") ) {
                if( line.startsWith(">chr") ) {
                    int spacePos = line.indexOf(" ");
                    out.write(">Chr"+line.substring(4, spacePos)+"\n");
                    scaffold = false;
                    linesWritten++;
                } else {
                    scaffold = true;
                }
                continue;
            }

            // data line
            if( !scaffold ) {
                out.write(line);
                out.write("\n");
                linesWritten++;
            }
        }

        in.close();
        out.close();

        System.out.println("lines read: "+linesRead);
        System.out.println("lines written: "+linesWritten);

        System.exit(0);
    }

    static void prepUthGff3Files() throws Exception {

        String fname =   "/Users/mtutaj/Downloads/rat4/wky.gff3.gz";
        String outName = "/Users/mtutaj/Downloads/rat4/wky2.gff3.gz";

        // copy comment lines as is
        // convert chr?  into Chr? lines
        // JAL...  lines are skipped
        BufferedReader in = Utils.openReader(fname);
        BufferedWriter out = Utils.openWriter(outName);
        String line;
        int linesRead = 0;
        int linesWritten = 0;

        while( (line=in.readLine())!=null ) {

            linesRead++;
            if( line.startsWith("#") ) {
                out.write(line);
                out.write("\n");
                linesWritten++;
            }
            else if( line.startsWith("chr") ) {
                out.write("C");
                out.write(line.substring(1));
                out.write("\n");
                linesWritten++;
            }
        }

        in.close();
        out.close();

        System.out.println("lines read: "+linesRead);
        System.out.println("lines written: "+linesWritten);

        System.exit(0);
    }
}
