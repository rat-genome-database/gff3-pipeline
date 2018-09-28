package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;

import java.io.PrintWriter;
import java.util.*;

/**
 * @author pjayaraman
 * Date: Nov 15, 2010
 */
public class CreateGff4SSLP {

    String toFile;
    int mapKey;
    int speciesTypeKey;

    RgdGff3Dao dao = new RgdGff3Dao();

    int countSslps=0;
    int mapSslpCount=0;
    int noMapSslpPos=0;
    int moreThanOneMapSslpPos=0;
    int sslpsAssocGeneRgd=0;
    int sslpsAliasCount=0;

    public void creategff4sslps(boolean compress) throws Exception{

        int mapKey = getMapKey();
        String species = SpeciesType.getCommonName(speciesTypeKey);

        String gffFile = getToFile()+species+"_RGDSSLPS.gff3";
        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gffFile, false, compress);

        // set up writer to create fasta sequences of primer pairs..
        PrintWriter seqPrimerPairsFile = new PrintWriter(getToFile()+species+"_RGDSSLPS.fa");
        FastaWriter fastaWriter = new FastaWriter(seqPrimerPairsFile);
        fastaWriter.setFastaSeqWriter(seqPrimerPairsFile);

        // set up gff file for liftOver
        Gff3ColumnWriter gff4LiftOverCWriter = new Gff3ColumnWriter(getToFile()+species+"_RGD_liftOver.txt", true, compress);

        //initialize RATMINE GFF3Writer
        Gff3ColumnWriter RATMINEgff3Writer = new Gff3ColumnWriter(getToFile()+species+"_RATMINE_RGDSSLPS.gff3", false, compress);
        RATMINEgff3Writer.setRatmineCompatibleFormat(true);

        for( SSLP sslp: dao.getActiveSslps(speciesTypeKey) ) {

            countSslps++;

            int sslpRgdId = sslp.getRgdId();
            String sslpSymbol = sslp.getName();

            String expectedSize="NA";
            if( sslp.getExpectedSize()!=null && sslp.getExpectedSize()>0 ) {
                expectedSize = sslp.getExpectedSize().toString();
            }

            String forwardSeq="", reverseSeq="";
            List<Sequence> seqs = dao.getSslpSequences(sslpRgdId);
            for( Sequence seq: seqs ) {
                if( seq.getSeqTypeKey()==4 ) {
                    // primer pair sequence
                    forwardSeq = seq.getForwardSeq();
                    reverseSeq = seq.getReverseSeq();
                }
            }


            String assocGeneRgdID = "NA";
            Gene gene = dao.getGeneBySslpKey(sslp.getKey());
            if( gene!=null ) {
                sslpsAssocGeneRgd++;
            }


            String aliases = getAliases(sslpRgdId);

            //get maps data position..
            List<MapData> sslpsMapdataList = dao.getMapData(sslpRgdId, mapKey);

            if(sslpsMapdataList.size()==0){
             noMapSslpPos++;
            }else if(sslpsMapdataList.size()>1){
             moreThanOneMapSslpPos++;
            }

            if( forwardSeq==null ) {
                //System.out.println("no forward primer pair sequence for SSLP RGDID:"+sslpRgdId);
            }
            if( reverseSeq==null ) {
                //System.out.println("no reverse primer pair sequence for SSLP RGDID:"+sslpRgdId);
            }
            fastaWriter.writeFastaSequences("forward_"+sslpRgdId, forwardSeq, "reverse_"+sslpRgdId, reverseSeq);

            if(sslpsMapdataList.size()>0){

             mapSslpCount++;

             for(MapData md: sslpsMapdataList){
                String chrom = md.getChromosome();
                Integer start = md.getStartPos();
                Integer stop = md.getStopPos();
                String strand = md.getStrand();
                if(strand==null){
                    strand = ".";
                }

                gff3Writer.writeFirst8Columns(chrom,"RGD","SSLPS",start,stop,".",strand,".");
                RATMINEgff3Writer.writeFirst8Columns(chrom,"RGD","SimpleSequenceLengthVariation",start,stop,".",strand,".");
                gff4LiftOverCWriter.print("Chr"+chrom+"\t"+start+"\t"+stop+"\t"+sslpRgdId+"_"+start+"_"+stop+"\n");

                //initialize hashmap for attributes
                HashMap<String, String> attributesHashMap = new HashMap<>();
                HashMap<String, String> RATMINEattributesHashMap = new HashMap<>();

                attributesHashMap.put("ID", sslpRgdId+"_"+start+"_"+stop);
                attributesHashMap.put("Name", "SSLP:"+sslpSymbol);
                if( aliases!=null )
                    attributesHashMap.put("Alias", aliases+","+"RGD:"+sslpRgdId+","+sslpRgdId+","+sslpSymbol);
                 else
                    attributesHashMap.put("Alias", "RGD:"+sslpRgdId+","+sslpRgdId+","+sslpSymbol);

                attributesHashMap.put("Dbxref", "RGD:"+sslpRgdId);
                attributesHashMap.put("expectedSize",expectedSize);
                attributesHashMap.put("associatedGene", assocGeneRgdID);

                gff3Writer.writeAttributes4Gff3(attributesHashMap);

                RATMINEattributesHashMap.put("ID", String.valueOf(sslpRgdId));
                RATMINEattributesHashMap.put("Name", sslpSymbol);
                 if( aliases!=null )
                    RATMINEattributesHashMap.put("Alias", aliases+","+"RGD:"+sslpRgdId+","+sslpRgdId);
                 else
                    RATMINEattributesHashMap.put("Alias", "RGD:"+sslpRgdId+","+sslpRgdId);
                RATMINEattributesHashMap.put("Dbxref", "RGD:"+sslpRgdId);
                RATMINEattributesHashMap.put("expectedSize",expectedSize);
                RATMINEattributesHashMap.put("associatedGene", assocGeneRgdID);

                RATMINEgff3Writer.writeAttributes4Gff3(RATMINEattributesHashMap);

                gff4LiftOverCWriter.addnewLineInGff3();
            }//end for each map data list loop
            }
        }

        gff3Writer.close();
        RATMINEgff3Writer.close();
        seqPrimerPairsFile.close();
        gff4LiftOverCWriter.close();


        System.out.println("SSLPS processed:" + countSslps);
        System.out.println("Sslps with Map position:" + mapSslpCount);
        System.out.println("Sslps with NO map position:" + noMapSslpPos);
        System.out.println("Sslps with more than one map position:" + moreThanOneMapSslpPos);
        System.out.println("Sslps having aliases:" + sslpsAliasCount);
        System.out.println("Sslps having an associated gene RGDID:" + sslpsAssocGeneRgd);

        System.out.print("\nGFF3 File SUCCESSFUL!");

    }

    String getAliases(int sslpRgdId) throws Exception {
        Set<String> aliases = new TreeSet<>();
        for( Alias a: dao.getAliases(sslpRgdId) ) {
            aliases.add(a.getValue());
        }

        if( aliases.isEmpty() ) {
            return null;
        }

        sslpsAliasCount++;
        return Utils.concatenate(aliases, ",");
    }

    public String getToFile() {
        return toFile;
    }

    public void setToFile(String toFile) {
        this.toFile = toFile;
    }

    public int getMapKey() {
        return mapKey;
    }

    public void setMapKey(int mapKey) {
        this.mapKey = mapKey;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }
}
