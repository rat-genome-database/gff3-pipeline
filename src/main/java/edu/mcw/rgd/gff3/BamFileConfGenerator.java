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

        HashMap<String,String> nameToFileMap = loadName2FileMap();

        String outFileName = "data/bam.conf";
        BufferedWriter out = Utils.openWriter(outFileName);

        String jbrowseCategory = "Chinchilla Expression Data/RNA-seq - PRJNA78823";
        boolean useBai = true;
        String trackIdBase = "PRJNA78823";
        int trackIndex = 0;
        String outDir = "chinchilla_RNAseq/";

        for( Map.Entry<String, String> entry: nameToFileMap.entrySet() ) {

            String bamFileName = entry.getKey();
            String bamTrackName = entry.getValue();

            trackIndex++;
            String trackId = trackIdBase + "_" + trackIndex;

            out.write("\n");
            out.write("[tracks."+trackId+"]\n");
            out.write("storeClass = JBrowse/Store/SeqFeature/BAM\n");
            out.write("urlTemplate = "+outDir+bamFileName+"\n");
            if( useBai ) {
                out.write("baiUrlTemplate = "+outDir+bamFileName+".bai\n");
            }
            out.write("category = "+jbrowseCategory+"\n");
            out.write("type = JBrowse/View/Track/Alignments2\n");
            out.write("key = "+bamTrackName+"\n");

            // generate coverage tracks

            out.write("\n");
            out.write("[tracks."+trackId+"c]\n");
            out.write("storeClass = JBrowse/Store/SeqFeature/BAM\n");
            out.write("urlTemplate = "+outDir+bamFileName+"\n");
            if( useBai ) {
                out.write("baiUrlTemplate = "+outDir+bamFileName+".bai\n");
            }
            out.write("category = "+jbrowseCategory+"\n");
            out.write("type = JBrowse/View/Track/SNPCoverage\n");
            out.write("key = "+bamTrackName+" (Coverage)\n");

        }

        out.close();
    }

    HashMap<String,String> loadName2FileMap() throws IOException {

        HashMap<String,String> nameToFileMap = new HashMap<>();


        String fname = "data/chinchilla_bam.txt";
        BufferedReader in = Utils.openReader(fname);
        String line;
        while( (line=in.readLine())!=null ) {

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
}
