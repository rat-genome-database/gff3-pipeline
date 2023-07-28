package edu.mcw.rgd.gff3;

import edu.mcw.rgd.process.Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class BamFileConfGenerator {

    public static void main(String[] args) throws IOException {

        new BamFileConfGenerator().run();
    }

    void run() throws IOException {

        String fname = "data/rn7bam.conf";
        Info info = new Info();
        HashMap<String,String> nameToFileMap = loadName2FileMap(fname, info);

        run(nameToFileMap, info);


        // chinchilla();
    }

    void chinchilla() throws IOException {

        String fname = "data/chinchilla_bam.txt";
        fname = "data/otitis_media.txt";

        Info info = new Info();
        HashMap<String,String> nameToFileMap = loadName2FileMap(fname, info);

        info.outFileName = "data/bam.conf";

        info.jbrowseCategory = "Chinchilla Expression Data/Otitis media - response to Streptococcus pneumoniae infection - PRJNA277957";
        info.useBai = true;
        info.trackIdBase = "PRJNA277957";
        info.bamDir = "PRJNA277957/"; // must be '/'-terminated

        run(nameToFileMap, info);
    }

    void run( HashMap<String,String> nameToFileMap, Info info ) throws IOException {

        BufferedWriter out = Utils.openWriter(info.outFileName);

        int trackIndex = 0;
        String outDir = info.bamDir;

        for( Map.Entry<String, String> entry: nameToFileMap.entrySet() ) {

            String bamFileName = entry.getKey();
            String bamTrackName = entry.getValue();

            trackIndex++;
            String trackId = info.trackIdBase + "_" + trackIndex;

            if( bamFileName.endsWith("bam") ) {

                out.write("\n");
                out.write("[tracks." + trackId + "]\n");
                out.write("storeClass = JBrowse/Store/SeqFeature/BAM\n");
                out.write("urlTemplate = " + outDir + bamFileName + "\n");
                if (info.useBai) {
                    out.write("baiUrlTemplate = " + outDir + bamFileName + ".bai\n");
                }
                out.write("category = " + info.jbrowseCategory + "\n");
                out.write("type = JBrowse/View/Track/Alignments2\n");
                out.write("key = " + bamTrackName + "\n");

                // generate coverage tracks

                out.write("\n");
                out.write("[tracks." + trackId + "c]\n");
                out.write("storeClass = JBrowse/Store/SeqFeature/BAM\n");
                out.write("urlTemplate = " + outDir + bamFileName + "\n");
                if (info.useBai) {
                    out.write("baiUrlTemplate = " + outDir + bamFileName + ".bai\n");
                }
                out.write("category = " + info.jbrowseCategory + "\n");
                out.write("type = JBrowse/View/Track/SNPCoverage\n");
                out.write("key = " + bamTrackName + " (Coverage)\n");
            }

            else if(bamFileName.endsWith("vcf.gz") ) {

                out.write("\n");
                out.write("[tracks." + trackId + "]\n");
                out.write("storeClass = JBrowse/Store/SeqFeature/VCFTabix\n");
                out.write("urlTemplate = " + outDir + bamFileName + "\n");
                out.write("category = " + info.jbrowseCategory + "\n");
                out.write("type = JBrowse/View/Track/CanvasVariants\n");
                out.write("key = " + bamTrackName + "\n");
            }
        }

        out.close();
    }

    HashMap<String,String> loadName2FileMap(String fname, Info info) throws IOException {

        HashMap<String,String> nameToFileMap = new HashMap<>();

        BufferedReader in = Utils.openReader(fname);
        String line;
        while( (line=in.readLine())!=null ) {

            if( line.startsWith("#") ) {
                if (line.startsWith("#OUT=")) {
                    info.outFileName = line.substring(5);
                } else if (line.startsWith("#JBCAT=")) {
                    info.jbrowseCategory = line.substring(7);
                } else if (line.startsWith("#BAI=")) {
                    String val = line.substring(5);
                    int i = Integer.parseInt(val);
                    info.useBai = i > 0;
                } else if (line.startsWith("#TRACKIDBASE=")) {
                    info.trackIdBase = line.substring(13);
                } else if (line.startsWith("#BAMDIR=")) {
                    info.bamDir = line.substring(8);
                    if( !info.bamDir.endsWith("/") ) {
                        info.bamDir += "/";
                    }
                }
                continue;
            }

            String s = line.trim();

            // find first blank character
            int blankPos = -1;
            for( int i=0; i<s.length(); i++ ) {
                if( Character.isWhitespace(line.charAt(i)) ) {
                    blankPos = i;
                    break;
                }
            }
            if( blankPos<0 ) {
                continue;
            }

            String fileName = s.substring(0, blankPos).trim();
            String trackName = s.substring(blankPos).trim();
            nameToFileMap.put(fileName, trackName);
        }
        in.close();

        return nameToFileMap;
    }

    class Info {

        public String outFileName;
        public String jbrowseCategory;
        public boolean useBai = true;
        public String trackIdBase;
        public String bamDir;
    }
}
