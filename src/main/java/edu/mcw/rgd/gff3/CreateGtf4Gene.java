package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;
import java.util.Map;

public class CreateGtf4Gene {

    RgdGff3Dao dao = new RgdGff3Dao();
    Logger log = LogManager.getLogger("gtf");

    private String outDir;
    private List<Integer> processedMapKeys;

    /**
     * load the species list and assemblies from properties/AppConfigure.xml
     */
    public void run() {

        processedMapKeys.stream().forEach( mapKey -> {

            try {
                createGeneGtf(mapKey);
            } catch(Exception e) {
                throw new RuntimeException(e);
            }
        });
    }

    public void createGeneGtf( int mapKey ) throws Exception{

        StringBuffer msgBuf = new StringBuffer();

        long time0 = System.currentTimeMillis();

        CdsUtils utils = new CdsUtils(dao, mapKey);

        int speciesTypeKey = SpeciesType.getSpeciesTypeKeyForMap(mapKey);
        String speciesName = SpeciesType.getCommonName( speciesTypeKey );

        String ucscId = Utils.NVL( Gff3Utils.getAssemblySymbol(mapKey), "" );
        String refseqId = MapManager.getInstance().getMap(mapKey).getRefSeqAssemblyName();

        if( ucscId.isEmpty() ) {
            msgBuf.append("Generate GTF file for " + speciesName + ", MAP_KEY=" + mapKey + "\n");
        } else {
            msgBuf.append("Generate GTF file for " + speciesName + ", MAP_KEY=" + mapKey + " (" + ucscId + ")\n");
        }
        msgBuf.append("    "+dao.getConnectionInfo()+"\n");

        CounterPool counters = new CounterPool();

        String headerInfo =
            "## source: RAT GENOME DATABASE (https://rgd.mcw.edu/)\n"+
            "## species: "+ speciesName+"\n"+
            "## assembly: "+ MapManager.getInstance().getMap(mapKey).getName()+"\n"+
            "## primary_contact: mtutaj@mcw.edu\n"+
            "## generated: "+new Date()+"\n";

        String fileName = getOutDir() + "/" + ucscId;
        if( ucscId.isEmpty() ) {
            fileName = getOutDir() + "/" + refseqId;
        }
        String fullFileName = fileName+".gtf";
        //String fullFileName = fileName+" genes and transcripts.gtf";
        log.info("creating file "+fullFileName);
        System.out.println("creating file "+fullFileName);

        GtfWriter gtfWriter = new GtfWriter(fullFileName, Gff3ColumnWriter.COMPRESS_MODE_BGZIP);
        gtfWriter.print(headerInfo);

        List<Gene> activeGenes = dao.getActiveGenes(speciesTypeKey);
        Collections.shuffle(activeGenes);
        log.info("   genes to be processed "+activeGenes.size());
        System.out.println("   genes to be processed "+activeGenes.size());

        int iGene = 0;
        for( Gene gene: activeGenes ){

            int geneRgdId = gene.getRgdId();
            List<MapData> geneMap = getMapData(geneRgdId, mapKey, counters);
            if( geneMap.isEmpty() ) {
                //System.out.println("no map positions");
                continue;
            }
            counters.increment(" Genes processed");

            List<XdbId> xdbIds = getXdbIds(geneRgdId);

            List<Transcript> geneTrs = dao.getTranscriptsForGene(geneRgdId);

            String geneId = "RGD"+gene.getRgdId();
            String geneType = getGeneType(gene.getType());
            String geneSymbol = gene.getSymbol();
            String geneName = getNameOfGene(gene);
            String dbXref = getDbXref(xdbIds);

            System.out.println((++iGene)+". RGD:"+gene.getRgdId()+" "+geneType+" "+geneSymbol+" "+dbXref+" ["+geneName+"]");

            for(MapData map : geneMap){

                List<Transcript> trsOnMap = getTranscriptsForMap(geneTrs, map, utils);

                Map<String,String> attributesHashMap = new HashMap<>();

                attributesHashMap.put("gene_type", geneType);
                if( geneSymbol != null ) {
                    attributesHashMap.put("gene_symbol", geneSymbol);
                }
                if( geneName != null ) {
                    attributesHashMap.put("gene_name", geneName);
                }
                if( dbXref != null ) {
                    attributesHashMap.put("db_xref", dbXref);
                }

                gtfWriter.writeFirst8Columns(map.getChromosome(),"RGD", "gene", map.getStartPos(), map.getStopPos(),".", map.getStrand(),".");
                gtfWriter.writeAttributes(geneId, "", attributesHashMap);

                for( Transcript tr: trsOnMap ){

                    String transcriptId = tr.getAccId();
                    String transcriptType = CreateGff4GeneAgr.getTrBiotype(gene.getType(), tr);
                    String proteinId = tr.getProteinAccId();

                    for( MapData trMd: tr.getGenomicPositions() ) {
                        if( !CdsUtils.transcriptPositionOverlapsGenePosition(trMd, map) )
                            continue;

                        attributesHashMap.put("gene_type", geneType);
                        if( geneSymbol != null ) {
                            attributesHashMap.put("gene_symbol", geneSymbol);
                        }
                        if( geneName != null ) {
                            attributesHashMap.put("gene_name", geneName);
                        }
                        if( transcriptType != null ) {
                            attributesHashMap.put("transcript_type", transcriptType);
                        }
                        if( proteinId != null ) {
                            attributesHashMap.put("protein_id", proteinId);
                        }

                        gtfWriter.writeFirst8Columns(trMd.getChromosome(), "RGD", "transcript", trMd.getStartPos(),trMd.getStopPos(), ".", trMd.getStrand(), ".");
                        gtfWriter.writeAttributes(geneId, transcriptId, attributesHashMap);


                        List<edu.mcw.rgd.gff3.CodingFeature> cfList = utils.buildCfList(trMd);
                        for(edu.mcw.rgd.gff3.CodingFeature cf: cfList){

                            String featureType = null;

                            if(cf.getFeatureType()== TranscriptFeature.FeatureType.CDS){
                                featureType = "CDS";
                            }
                            else if(cf.getFeatureType()== TranscriptFeature.FeatureType.EXON) {
                                featureType = "exon";
                            }
                            if( featureType == null ) {
                                continue;
                            }

                            attributesHashMap.put("gene_type", geneType);
                            if( geneSymbol != null ) {
                                attributesHashMap.put("gene_symbol", geneSymbol);
                            }
                            if( geneName != null ) {
                                attributesHashMap.put("gene_name", geneName);
                            }
                            if( transcriptType != null ) {
                                attributesHashMap.put("transcript_type", transcriptType);
                            }
                            if( proteinId != null ) {
                                attributesHashMap.put("protein_id", proteinId);
                            }

                            gtfWriter.writeFirst8Columns(trMd.getChromosome(), "RGD", featureType, cf.getStartPos(), cf.getStopPos(), ".", trMd.getStrand(), cf.getCodingPhaseStr());
                            gtfWriter.writeAttributes(geneId, transcriptId, attributesHashMap);
                        }
                    }
                }
            }
        }

        gtfWriter.close();
        gtfWriter.sortInMemory();

        dumpCounters(counters, ucscId.isEmpty() ? refseqId : ucscId, msgBuf);

        msgBuf.append("OK!  elapsed "+Utils.formatElapsedTime(time0, System.currentTimeMillis())+"\n");
        msgBuf.append("========\n");
        synchronized( this.getClass() ) {
            log.info(msgBuf);
            System.out.println(msgBuf);
        }
    }

    // get a subset of transcripts having position on a given map, and overlapping the given gene position
    List<Transcript> getTranscriptsForMap(List<Transcript> geneTrs, MapData md, CdsUtils utils) {

        List<Transcript> trsOnMap = new ArrayList<>();
        for( Transcript tr: geneTrs ) {
            for( MapData trMd: tr.getGenomicPositions() ) {
                // transcript pos must overlap gene locus
                if( CdsUtils.transcriptPositionOverlapsGenePosition(trMd, md) ) {
                    // yes, this transcript is overlapping gene locus
                    trsOnMap.add(tr);
                }
            }
        }
        return trsOnMap;
    }

    List<MapData> getMapData(int geneRgdId, int mapKey, CounterPool counters) throws Exception {

        List<MapData> geneMap = dao.getMapData(geneRgdId, mapKey);

        if( !geneMap.isEmpty() ) {
            if(geneMap.size()>1){
                counters.increment(" Genes with more than one map position");
            }
        }
        return geneMap;
    }

    void dumpCounters(CounterPool counters, String assembly, StringBuffer msgBuf) {

        String msg = counters.dumpAlphabetically();

        // prepend all lines with 'assembly'
        String[] lines = msg.split("[\\n]", -1);
        for( String line: lines ) {
            msgBuf.append(assembly).append("  ").append(line).append("\n");
        }
    }

    String getGeneType(String geneType) throws Exception {

        if( geneType.startsWith("predicted") )
            return "predicted_gene";

        switch( geneType ) {
            case "pseudogene":
            case "pseudo":
                return "pseudogene";
            case "protein-coding":
                return "protein_coding";
            case "ncrna":
                return "ncRNA";
            case "gene":
                return "gene";
            case "rrna":
                return "rRNA";
            case "trna":
                return "tRNA";
            case "snrna":
                return "snRNA";
            case "scrna":
                return "scRNA";
            case "snorna":
                return "snoRNA";
            case "miscrna":
                return "miscRNA";
            default:
                throw new Exception("unsupported gene type " + geneType);
        }
    }

    // get NCBI gene id; if not found: Ensembl gene id; otherwise null
    private String getDbXref(List<XdbId> xdbList) {

        // return NCBI GeneID if available
        for(XdbId externalId : xdbList ){
            if( externalId.getXdbKey()==XdbId.XDB_KEY_NCBI_GENE ) {
                return "GeneID:"+externalId.getAccId();
            }
        }

        // return Ensembl Gene ID if available
        for(XdbId externalId : xdbList ){
            if( externalId.getXdbKey()==XdbId.XDB_KEY_ENSEMBL_GENES ) {
                return "Ensembl:"+externalId.getAccId();
            }
        }

        return null;
    }

    /**
     * returns a comma separated string of key:value pairs that are external identifiers
     * @param geneRgdId gene rgd id
     * @return xdb ids string
     * @throws Exception
     */
    private List<XdbId> getXdbIds(int geneRgdId) throws Exception{

        XdbId extId = new XdbId();
        extId.setRgdId(geneRgdId);

        List<XdbId> xdbList = dao.getXdbIds(extId);
        xdbList.removeIf( externalId ->
            externalId.getXdbKey() != XdbId.XDB_KEY_NCBI_GENE &&
            externalId.getXdbKey() != XdbId.XDB_KEY_ENSEMBL_GENES);
        return xdbList;
    }

    private String getNameOfGene(Gene gene) {

        String nameOfgene = gene.getName();

        if( nameOfgene!=null ){
            if(nameOfgene.contains(",")){
                nameOfgene = nameOfgene.replaceAll(",", "");
            }
            if(nameOfgene.contains(";")){
                nameOfgene = nameOfgene.replaceAll(";", "");
            }
        }
        return nameOfgene;
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
