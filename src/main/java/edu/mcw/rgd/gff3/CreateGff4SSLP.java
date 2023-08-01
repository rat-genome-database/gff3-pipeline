package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;

import java.io.BufferedWriter;
import java.util.*;

/**
 * @author pjayaraman
 * Date: Nov 15, 2010
 */
public class CreateGff4SSLP {

    RgdGff3Dao dao = new RgdGff3Dao();

    private String outDir;
    private List<Integer> processedMapKeys;

    /**
     * load the species list and assemblies from properties/AppConfigure.xml
     */
    public void run() throws Exception {

        for( int mapKey: getProcessedMapKeys() ) {

            String assemblyDir = Manager.getInstance().getAssemblies().get(mapKey);
            if( assemblyDir==null ) {
                break;
            }

            CreateInfo info = new CreateInfo();

            int speciesTypeKey = MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey();

            info.setMapKey( mapKey );
            info.setToDir( assemblyDir + "/" + getOutDir() );
            info.setSpeciesTypeKey( speciesTypeKey );
            info.setCompressMode( Gff3ColumnWriter.COMPRESS_MODE_BGZIP );

            createGff4Markers(info);
        }
    }

    public void createGff4Markers(CreateInfo info) throws Exception {

        CounterPool counters = new CounterPool();

        String speciesName = SpeciesType.getCommonName(info.getSpeciesTypeKey());
        String ucscId = Gff3Utils.getAssemblySymbol(info.getMapKey());
        String refseqId = MapManager.getInstance().getMap(info.getMapKey()).getRefSeqAssemblyName();
        String fileName = info.getToDir() + "/" + speciesName + " " + refseqId+" ("+ucscId+")";

        System.out.println("START MARKER GFF3 Generator for  " + speciesName + "  MAP_KEY=" + info.getMapKey() + "  ASSEMBLY " + MapManager.getInstance().getMap(info.getMapKey()).getName());
        System.out.println("========================");


        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(fileName + " Markers.gff3", false, info.getCompressMode());
        gff3Writer.print("# RAT GENOME DATABASE (https://rgd.mcw.edu/)\n");
        gff3Writer.print("# Species: "+ speciesName+"\n");
        gff3Writer.print("# Assembly: "+ MapManager.getInstance().getMap(info.getMapKey()).getName()+"\n");
        gff3Writer.print("# Primary Contact: mtutaj@mcw.edu\n");
        gff3Writer.print("# Generated: "+new Date()+"\n");

        // set up writer to create fasta sequences of primer pairs.
        String fastaFileName = fileName + " Markers.fa";
        if( info.getCompressMode()!=Gff3ColumnWriter.COMPRESS_MODE_NONE ) {
            fastaFileName += ".gz";
        }
        BufferedWriter fastaWriterFile = Utils.openWriter(fastaFileName);
        FastaWriter fastaWriter = new FastaWriter(fastaWriterFile);

        createGff4Markers(info.getMapKey(), info.getSpeciesTypeKey(), gff3Writer, fastaWriter, counters);

        gff3Writer.close();
        fastaWriterFile.close();

        gff3Writer.sortInMemory();
    }

    public void createGff4Markers(int mapKey, int speciesTypeKey, Gff3ColumnWriter gff3Writer, FastaWriter fastaWriter, CounterPool counters) throws Exception{

        SequenceRegionWatcher sequenceRegionWatcher = new SequenceRegionWatcher(mapKey, gff3Writer, dao);

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

            String expectedSize = null;
            if( sslp.getExpectedSize()!=null && sslp.getExpectedSize()>0 ) {
                expectedSize = sslp.getExpectedSize().toString();
            }

            String assocGeneRgdID = null;
            Gene gene = dao.getGeneBySslpKey(sslp.getKey());
            if( gene!=null ) {
                sslpsAssocGeneRgd++;
                assocGeneRgdID = "RGD:"+gene.getRgdId();
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

                    sequenceRegionWatcher.emit(chrom);

                    gff3Writer.writeFirst8Columns(chrom,"RGD","SSLPS",start,stop,".",strand,".");

                    //initialize hashmap for attributes
                    HashMap<String, String> attributesHashMap = new HashMap<>();

                    attributesHashMap.put("ID", sslpRgdId+"_"+start+"_"+stop);
                    attributesHashMap.put("Name", "SSLP:"+sslpSymbol);
                    if( aliases!=null ) {
                        if( aliases.contains("\t") || aliases.contains("\n") ) {
                            aliases = aliases.replace("\\s+", " ");
                        }
                        attributesHashMap.put("Alias", aliases+","+"RGD:"+sslpRgdId+","+sslpRgdId+","+sslpSymbol);
                    }
                    else
                        attributesHashMap.put("Alias", "RGD:"+sslpRgdId+","+sslpRgdId+","+sslpSymbol);

                    attributesHashMap.put("Dbxref", "RGD:"+sslpRgdId);
                    if( expectedSize!=null ) {
                        attributesHashMap.put("expectedSize", expectedSize);
                    }
                    if( assocGeneRgdID!=null ) {
                        attributesHashMap.put("associatedGene", assocGeneRgdID);
                    }

                    gff3Writer.writeAttributes4Gff3(attributesHashMap);
                }
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

    public String getOutDir() {
        return outDir;
    }

    public void setOutDir(String outDir) {
        this.outDir = outDir;
    }

    public List<Integer> getProcessedMapKeys() {
        return processedMapKeys;
    }

    public void setProcessedMapKeys(List<Integer> processedMapKeys) {
        this.processedMapKeys = processedMapKeys;
    }
}
