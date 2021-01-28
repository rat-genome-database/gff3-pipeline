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

        String fname = "/data/ref/fasta/rn7.2/GCF_015227675.2_mRatBN7.2_genomic.fna.gz";
        String outDir = "/data/ref/fasta/rn7.2";

        int mapKey = 372;
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
                String chr = getChrName(chrAcc, chromosomes);
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

    static String getChrName(String acc, List<Chromosome> chromosomes) {
        for( Chromosome ch: chromosomes ) {
            if( acc.equals(ch.getRefseqId()) ) {
                return "Chr"+ch.getChromosome();
            }
        }
        return null;
    }
}
