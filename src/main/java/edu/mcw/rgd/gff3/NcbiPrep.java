package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.datamodel.Chromosome;
import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.process.Utils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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

        String fname = "/chi/GCF_000276665.1_ChiLan1.0_genomic.fna.gz";
        String outDir = "/data/ref/fasta";
        prepScaffolds(fname, outDir);

        //prepUthGff3Files();
        //loadUthPositions();

        int mapKey = 634;
        boolean isScaffoldAssembly = false;

        // for scaffold assemblies, prefix must be an empty string; for chromosome assemblies it must be 'Chr'
        final String chrPrefix = isScaffoldAssembly ? "" : "Chr";

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

    static void prepScaffolds(String fname, String outDir) throws IOException {

        BufferedReader in = Utils.openReader(fname);
        String outName = outDir +"/out.fa.gz";
        BufferedWriter out = Utils.openWriter(outName);

        String line;
        while( (line = in.readLine())!=null ) {

        }
        in.close();
        out.close();
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

        String fname =   "/Users/mtutaj/Downloads/rat4/shrsp.gff3.gz";
        String outName = "/Users/mtutaj/Downloads/rat4/shrsp2.gff3.gz";

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

                //Chr1	BestRefSeq	gene	1761587	1794261	.	-	.	ID=gene-Pcmt1;Dbxref=GeneID:25604,RGD:3268;Name=Pcmt1;description=protein-L-isoaspartate (D-aspartate) O-methyltransferase 1;gbkey=Gene;gene=Pcmt1;gene_biotype=protein_coding;gene_synonym=PCM
                // remove attribute 'description' from the line
                int descPos = line.indexOf(";description=");
                if( descPos>0 ) {
                    int descPos2 = line.indexOf(";", descPos+12);
                    String line2;
                    if( descPos2>descPos ) {
                        line2 = line.substring(0, descPos) + line.substring(descPos2);
                    } else {
                        line2 = line.substring(0, descPos);
                    }
                    line = line2;
                }
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

    static void loadUthPositions() throws Exception {

        //String fname =   "/Users/mtutaj/Downloads/rat4/shr.gff3.gz";
        //int mapKey = 301;
        //String fname =   "/Users/mtutaj/Downloads/rat4/shrsp.gff3.gz";
        //int mapKey = 302;
        String fname =   "/Users/mtutaj/Downloads/rat4/wky.gff3.gz";
        int mapKey = 303;

        int geneLines = 0;
        int allLines = 0;
        List<UthRec> uthGenes = new ArrayList<>();

        BufferedReader in = Utils.openReader(fname);
        String line;
        while( (line=in.readLine())!=null ) {
            allLines++;

            // skip empty and comment lines
            if( line.isEmpty() || line.startsWith("#") ) {
                continue;
            }

            String[] cols = line.split("[\\t]", -1);

            // sample line
            //Chr1	BestRefSeq	gene	1761587	1794261	.	-	.	ID=gene-Pcmt1;Dbxref=GeneID:25604,RGD:3268;Name=Pcmt1;description=protein-L-isoaspartate (D-aspartate) O-methyltransferase 1;gbkey=Gene;gene=Pcmt1;gene_biotype=protein_coding;gene_synonym=PCM
            String attrField = cols[8];
            if( !attrField.startsWith("ID=gene-") ) {
                continue;
            }
            geneLines++;

            UthRec rec = new UthRec();
            rec.mapKey = mapKey;
            rec.chr = cols[0].substring(3).trim();
            rec.startPos = Integer.parseInt(cols[3]);
            rec.stopPos = Integer.parseInt(cols[4]);
            rec.strand = cols[6];

            Map<String,String> attrs = loadAttrs(attrField);
            String geneSymbol = attrs.get("Name");
            if( geneSymbol!=null ) {
                rec.geneSymbol = geneSymbol;
            }
            String dbxref = attrs.get("Dbxref");
            if( dbxref!=null ) {
                String[] xrefs = dbxref.split("[,]");
                for( String xref: xrefs ) {

                    if( xref.startsWith("GeneID:") ) {
                        String egId = xref.substring(7);
                        rec.geneId = egId;
                    } else if( xref.startsWith("RGD:") ) {
                        String rgdId = xref.substring(4);
                        rec.geneRgdId = Integer.parseInt(rgdId);
                    } else {
                        System.out.println("unknown dbxref: "+xref);
                    }
                }
            }
            uthGenes.add(rec);
        }
        in.close();

        RgdGff3Dao dao = new RgdGff3Dao();

        int genesWithoutRgdId = 0;
        for( UthRec rec: uthGenes ) {
            if( rec.geneRgdId==0 ) {
                genesWithoutRgdId++;
            } else {
                MapData md = new MapData();
                md.setChromosome(rec.chr);
                md.setMapKey(rec.mapKey);
                md.setRgdId(rec.geneRgdId);
                md.setSrcPipeline("NCBI");
                md.setStartPos(rec.startPos);
                md.setStopPos(rec.stopPos);
                md.setStrand(rec.strand);
                dao.insertMapData(md);
            }
        }

        System.out.println("read lines: "+allLines);
        System.out.println("gene lines: "+geneLines);
        System.out.println("gene lines, no rgd id: "+genesWithoutRgdId);

        System.exit(0);
    }

    static Map<String,String> loadAttrs(String attrField) {

        Map<String,String> attrs = new HashMap<>();
        String[] fields = attrField.split("[;]");
        for( String field: fields ) {
            int eqPos = field.indexOf("=");
            String attrName = field.substring(0, eqPos);
            String attrValue = field.substring(eqPos+1);
            attrs.put(attrName, attrValue);
        }
        return attrs;
    }
}

class UthRec {
    public int mapKey;
    public String chr;
    public int startPos;
    public int stopPos;
    public String strand;
    public String geneId;
    public String geneSymbol;
    public int geneRgdId;
}
