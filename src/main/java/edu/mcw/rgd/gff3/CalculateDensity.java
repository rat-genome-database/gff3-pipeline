package edu.mcw.rgd.gff3;

import java.io.PrintWriter;
import java.util.*;

/**
 * @author pjayaraman
 * Date: 11/17/11
 */
public class CalculateDensity {
    public final static long BLOCK_SIZE = 1000000;
    List<RGDInfo> rgdInfoList = new ArrayList<RGDInfo>();
    PrintWriter densityDiseaseWriter;

    public PrintWriter getDensityDiseaseWriter() {
        return densityDiseaseWriter;
    }

    public void setDensityDiseaseWriter(PrintWriter densityDiseaseWriter) {
        this.densityDiseaseWriter = densityDiseaseWriter;
    }

    public List<RGDInfo> getRgdInfoList() {
        return rgdInfoList;
    }

    public void setRgdInfoList(List<RGDInfo> rgdInfoList) {
        this.rgdInfoList = rgdInfoList;
    }

    void runDensityCalculator() throws Exception {

        System.out.print("\nDensity counts:\n  ");

        // map of chromosome block nr to hit count
        // f.e. '1-0' => 5  means on chromosome 1 block 0 located at pos 1..QTL_BLOCK_SIZE has 5 hits
        Map<String, Integer> map1 = new HashMap<String, Integer>();
        for( RGDInfo info: rgdInfoList ) {
            // convert start position into block nr
            int startBlock = (int)(info.startPos / BLOCK_SIZE);
            int stopBlock = (int)(info.stopPos / BLOCK_SIZE);
            if( stopBlock < startBlock ) {
                int tmp = startBlock;
                startBlock = stopBlock;
                stopBlock = tmp;
            }
            for( int block=startBlock; block<=stopBlock; block++ ) {
                String chblock = "chr "+info.chromosome + "_"
                        + (1+block*BLOCK_SIZE) + "_" + ((1+block)*BLOCK_SIZE);
                Integer hitCount = map1.get(chblock);
                if( hitCount==null )
                    hitCount = 1;
                else
                    hitCount++;
                map1.put(chblock, hitCount);
            }
        }


        List<String> chblockKeys = new ArrayList<String>(map1.keySet());
        Collections.sort(chblockKeys);

        String prevChr = "";
        int countPerChr = 0;

        for(String key: chblockKeys){
            //System.out.println(key + " has count: " + map1.get(key));
            int cnt = map1.get(key);
            String[] chrStrtStp = key.split("_");
            String[] chrNum = chrStrtStp[0].split("chr ");
            String chrom = chrNum[1];
            String strt = chrStrtStp[1];
            String stp = chrStrtStp[2];

            if( chrom.equals(prevChr) ) {
                countPerChr += cnt;
            }
            else {
                if( !prevChr.isEmpty() ){
                    System.out.print("c"+prevChr + ":" + countPerChr+", ");
                }
                countPerChr = cnt;
                prevChr = chrom;
            }
            densityDiseaseWriter.println("Chr"+chrom+"\tdensity\tbin\t"+strt+"\t"+stp+"\t"+cnt+"\t+\t.\tbin Chr"+chrom+":density");

        }
        if( !prevChr.isEmpty() ){
            System.out.println("c"+prevChr + ":" + countPerChr);
        }

        densityDiseaseWriter.close();
    }
}
