package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;

import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.util.*;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * @author pjayaraman
 * Date: Nov 15, 2010
 */
public class CreateGff4SSLP {

    RgdGff3Dao dao = new RgdGff3Dao();

    private List<String> processedAssemblies;

    /**
     * load the species list and assemblies from properties/AppConfigure.xml
     */
    public void run() throws Exception {
        for( String assemblyInfo: processedAssemblies ) {
            CreateInfo info = new CreateInfo();
            info.parseFromString(assemblyInfo);

            createGff4Markers(info);
        }
    }

    public void createGff4Markers(CreateInfo info) throws Exception {

        CounterPool counters = new CounterPool();

        String speciesName = SpeciesType.getCommonName(info.getSpeciesTypeKey());
        String assemblySymbol = info.getAssemblySymbol()!=null ? info.getAssemblySymbol() : Gff3Utils.getAssemblySymbol(info.getMapKey());

        System.out.println("START MARKER GFF3 Generator for  " + speciesName + "  MAP_KEY=" + info.getMapKey() + "  ASSEMBLY " + MapManager.getInstance().getMap(info.getMapKey()).getName());
        System.out.println("========================");


        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(info.getToDir()+"/"+assemblySymbol+"_markers.gff3", false, info.isCompress());
        gff3Writer.print("# RAT GENOME DATABASE (https://rgd.mcw.edu/)\n");
        gff3Writer.print("# Species: "+ speciesName+"\n");
        gff3Writer.print("# Assembly: "+ MapManager.getInstance().getMap(info.getMapKey()).getName()+"\n");
        gff3Writer.print("# Primary Contact: mtutaj@mcw.edu\n");
        gff3Writer.print("# Generated: "+new Date()+"\n");

        Gff3ColumnWriter RATMINEgff3Writer = new Gff3ColumnWriter(info.getToDir()+"/RATMINE_"+assemblySymbol+"_markers.gff3", false, info.isCompress());
        RATMINEgff3Writer.setRatmineCompatibleFormat(true);

        // set up writer to create fasta sequences of primer pairs..
        String fastaFileName = info.getToDir()+"/"+assemblySymbol+"_markers.fa";
        if( info.isCompress() ) {
            fastaFileName += ".gz";
        }
        BufferedWriter fastaWriterFile = Utils.openWriter(fastaFileName);
        FastaWriter fastaWriter = new FastaWriter(fastaWriterFile);

        createGff4Markers(info.getMapKey(), info.getSpeciesTypeKey(), gff3Writer, RATMINEgff3Writer, fastaWriter, counters);

        gff3Writer.close();
        RATMINEgff3Writer.close();
        fastaWriterFile.close();
    }

    public void createGff4Markers(int mapKey, int speciesTypeKey, Gff3ColumnWriter gff3Writer, Gff3ColumnWriter RATMINEgff3Writer, FastaWriter fastaWriter, CounterPool counters) throws Exception{

        int countSslps=0;
        int mapSslpCount=0;
        int noMapSslpPos=0;
        int moreThanOneMapSslpPos=0;
        int sslpsAssocGeneRgd=0;
        int sslpsAliasCount=0;

        for( SSLP sslp: dao.getActiveSslps(speciesTypeKey) ) {

            countSslps++;

            int sslpRgdId = sslp.getRgdId();
            String sslpSymbol = sslp.getName();

            String expectedSize="NA";
            if( sslp.getExpectedSize()!=null && sslp.getExpectedSize()>0 ) {
                expectedSize = sslp.getExpectedSize().toString();
            }

            String assocGeneRgdID = "NA";
            Gene gene = dao.getGeneBySslpKey(sslp.getKey());
            if( gene!=null ) {
                sslpsAssocGeneRgd++;
            }


            String aliases = getAliases(sslpRgdId);
            if( aliases!=null ) {
                sslpsAliasCount++;
            }

            //get maps data position..
            List<MapData> sslpsMapdataList = dao.getMapData(sslpRgdId, mapKey);

            if(sslpsMapdataList.size()==0){
             noMapSslpPos++;
            }else if(sslpsMapdataList.size()>1){
             moreThanOneMapSslpPos++;
            }

            fastaWriter.writeFastaSequences("forward_"+sslpRgdId, sslp.getForwardSeq(), "reverse_"+sslpRgdId, sslp.getReverseSeq());

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
            }//end for each map data list loop
            }
        }

        System.out.println("Markers processed: " + countSslps);
        System.out.println("Markers with Map position: " + mapSslpCount);
        System.out.println("Markers with NO map position: " + noMapSslpPos);
        System.out.println("Markers with more than one map position: " + moreThanOneMapSslpPos);
        System.out.println("Markers having aliases: " + sslpsAliasCount);
        System.out.println("Markers having an associated gene RGDID: " + sslpsAssocGeneRgd);
        System.out.println("");
    }

    String getAliases(int sslpRgdId) throws Exception {
        Set<String> aliases = new TreeSet<>();
        for( Alias a: dao.getAliases(sslpRgdId) ) {
            aliases.add(a.getValue());
        }

        if( aliases.isEmpty() ) {
            return null;
        }

        return Utils.concatenate(aliases, ",");
    }

    public void setProcessedAssemblies(List processedAssemblies) {
        this.processedAssemblies = processedAssemblies;
    }

    public List getProcessedAssemblies() {
        return processedAssemblies;
    }
}
