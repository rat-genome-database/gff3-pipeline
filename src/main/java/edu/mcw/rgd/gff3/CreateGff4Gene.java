package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.atomic.AtomicInteger;

public class CreateGff4Gene {

    RgdGff3Dao dao = new RgdGff3Dao();
    Logger log = LogManager.getLogger("gene");

    private List<String> processedAssemblies;

    /**
     * load the species list and assemblies from properties/AppConfigure.xml
     */
    public void run() {

        processedAssemblies.parallelStream().forEach( assemblyInfo -> {

            CreateInfo info = new CreateInfo();
            try {
                info.parseFromString(assemblyInfo);

                createGeneGff3(info);
            } catch(Exception e) {
                throw new RuntimeException(e);
            }
        });
    }

    public void createGeneGff3(CreateInfo info) throws Exception{

        StringBuffer msgBuf = new StringBuffer();

        long time0 = System.currentTimeMillis();

        CdsUtils utils = new CdsUtils(dao, info.getMapKey());

        String species = SpeciesType.getCommonName(info.getSpeciesTypeKey());
        Map<String, AtomicInteger> idMap = new ConcurrentHashMap<>();

        String assemblySymbol = Gff3Utils.getAssemblySymbol(info.getMapKey());

        msgBuf.append("Generate GFF3 file for "+species+", MAP_KEY="+info.getMapKey()+" ("+assemblySymbol+")\n");
        msgBuf.append("    "+dao.getConnectionInfo()+"\n");

        // March 2 2016: RATMINE gff3 loader does not allow one gene to have multiple loci
        //               so we emit only first loci per gene
        Set<Integer> geneRgdIdsEmitted = new ConcurrentSkipListSet<>();

        CounterPool counters = new CounterPool();

        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(info.getToDir()+"/"+assemblySymbol+"_genes.gff3", false, info.isCompress());
        gff3Writer.print("# RAT GENOME DATABASE (https://rgd.mcw.edu/)\n");
        gff3Writer.print("# Species: "+ species+"\n");
        gff3Writer.print("# Assembly: "+ MapManager.getInstance().getMap(info.getMapKey()).getName()+"\n");
        gff3Writer.print("# Primary Contact: mtutaj@mcw.edu\n");
        gff3Writer.print("# Generated: "+new Date()+"\n");

        Gff3ColumnWriter RATMINEgff3Writer = new Gff3ColumnWriter(info.getToDir()+"/RATMINE_"+assemblySymbol+"_genes.gff3", false, info.isCompress());
        RATMINEgff3Writer.setRatmineCompatibleFormat(true);

        List<Gene> activeGenes = dao.getActiveGenes(info.getSpeciesTypeKey());

        for( Gene gene: activeGenes ){
            int geneRgdId = gene.getRgdId();
            List<MapData> geneMap = getMapData(geneRgdId, info.getMapKey(), counters);
            if( geneMap.isEmpty() ) {
                //System.out.println("no map positions");
                continue;
            }
            counters.increment(" Genes processed");

            String nameOfgene = getNameOfGene(gene, counters);

            List<Transcript> geneTrs = dao.getTranscriptsForGene(geneRgdId);

            String annotDesc = Utils.getGeneDescription(gene);
            if( Utils.isStringEmpty(annotDesc) ){
                annotDesc = null;
            }else{
                annotDesc = annotDesc.replaceAll(";"," AND ").replaceAll(",", " ");
            }


            String nc="N";

            for(MapData map : geneMap){

                boolean writeRATMINE = geneRgdIdsEmitted.add(map.getRgdId());
                if( !writeRATMINE ) {
                    counters.increment(" Gene Loci skipped for RATMINE");
                }

                List<Transcript> trsOnMap = getTranscriptsForMap(geneTrs, map, utils);

                String gType;
                if(gene.getType()==null){
                    gType="gene";
                    counters.increment(" Genes with NULL geneType");

                }else{
                    gType=gene.getType();

                    if(gene.getType().contains("pseudo")){
                        counters.increment(" Pseudo Genes");
                    }
                }

                List<XdbId> xdbIds = getXdbIds(geneRgdId);

                Map<String,String> attributesHashMap = new HashMap<>();
                Map<String,String> RATMINEattributesHashMap = new HashMap<>();

                String uniqueGeneId = getUniqueId("RGD"+gene.getRgdId(), idMap);
                attributesHashMap.put("ID", uniqueGeneId);
                attributesHashMap.put("Name", gene.getSymbol());
                attributesHashMap.put("fullName", nameOfgene);
                attributesHashMap.put("Alias", gene.getSymbol()+","+"RGD"+gene.getRgdId()+","+gene.getRgdId()
                        + "," + nameOfgene + getHgncMgiIds(xdbIds));
                attributesHashMap.put("geneType", gType.replaceAll("\\-","_"));
                attributesHashMap.put("species",species);
                if( gene.getRefSeqStatus()!=null )
                    attributesHashMap.put("refSeqStatus",gene.getRefSeqStatus());
                if( annotDesc!=null )
                    attributesHashMap.put("Note",annotDesc);
                if( !Utils.isStringEmpty(gene.getNcbiAnnotStatus()) ) {
                    attributesHashMap.put("nomenclatureStatus", gene.getNcbiAnnotStatus());
                }

                String RATMINEuniqueGeneId = getUniqueId("RGD:"+gene.getRgdId(), idMap);
                RATMINEattributesHashMap.put("ID", RATMINEuniqueGeneId);
                RATMINEattributesHashMap.put("Name", gene.getSymbol());
                RATMINEattributesHashMap.put("fullName", nameOfgene);
                RATMINEattributesHashMap.put("Alias", "RGD"+gene.getRgdId()+","+gene.getRgdId()+","+nameOfgene);
                if( annotDesc!=null )
                    RATMINEattributesHashMap.put("Note",annotDesc);
                RATMINEattributesHashMap.put("geneType",gene.getType());
                RATMINEattributesHashMap.put("nomenclatureStatus", gene.getNcbiAnnotStatus());


                String extDbString = getXdbString(xdbIds, counters, false);
                if( !extDbString.isEmpty() ){
                    attributesHashMap.put("Dbxref",extDbString);

                    extDbString = getXdbString(xdbIds, counters, true);
                    RATMINEattributesHashMap.put("Dbxref", extDbString);
                }

                gff3Writer.writeFirst8Columns(map.getChromosome(),"RGD", "gene", map.getStartPos(),map.getStopPos(),".",map.getStrand(),".");
                gff3Writer.writeAttributes4Gff3(attributesHashMap);

                if( writeRATMINE ) {
                    if(RATMINEuniqueGeneId.contains("_")) {
                        log.debug("RATMINE error: dash in ID="+RATMINEuniqueGeneId);
                    }
                    RATMINEgff3Writer.writeFirst8Columns(map.getChromosome(), "RGD", "gene", map.getStartPos(), map.getStopPos(), ".", map.getStrand(), ".");
                    RATMINEgff3Writer.writeAttributes4Gff3(RATMINEattributesHashMap);
                }

                if(trsOnMap.size()>1){
                    counters.increment(" Genes with more than one mapped transcript");
                }
                if(trsOnMap.size()>0){
                    counters.increment(" Genes with transcripts");

                    for( Transcript tr: trsOnMap ){

                        if(tr.isNonCoding()){
                            nc="Y";
                            counters.increment(" NonCoding transcripts");
                        }
                        for( MapData trMd: tr.getGenomicPositions() ) {
                            if( !CdsUtils.transcriptPositionOverlapsGenePosition(trMd, map) )
                                continue;

                            String id = getUniqueId("mRNARGD"+tr.getRgdId(), idMap);

                            counters.increment(" Mapped Transcripts");

                            attributesHashMap.put("ID", id);
                            attributesHashMap.put("Name", tr.getAccId());
                            attributesHashMap.put("Parent", uniqueGeneId);
                            attributesHashMap.put("Alias", "RGD:"+tr.getRgdId());
                            if( tr.getRefSeqStatus()!=null )
                                attributesHashMap.put("refSeqStatus", tr.getRefSeqStatus());
                            attributesHashMap.put("isNonCoding", nc);
                            attributesHashMap.put("gene", gene.getSymbol());

                            RATMINEattributesHashMap.put("ID", id);
                            RATMINEattributesHashMap.put("Name", tr.getAccId());
                            RATMINEattributesHashMap.put("Parent", RATMINEuniqueGeneId);
                            RATMINEattributesHashMap.put("Alias", "RGD:"+tr.getRgdId());
                            if( tr.getRefSeqStatus()!=null )
                                RATMINEattributesHashMap.put("RefSeqStatus", tr.getRefSeqStatus());
                            RATMINEattributesHashMap.put("isNon-Coding", nc);
                            RATMINEattributesHashMap.put("gene", gene.getSymbol());

                            gff3Writer.writeFirst8Columns(trMd.getChromosome(), "RGD", "mRNA", trMd.getStartPos(),trMd.getStopPos(), ".", trMd.getStrand(), ".");
                            gff3Writer.writeAttributes4Gff3(attributesHashMap);

                            if( writeRATMINE ) {
                                RATMINEgff3Writer.writeFirst8Columns(trMd.getChromosome(), "RGD", "mRNA", trMd.getStartPos(), trMd.getStopPos(), ".", trMd.getStrand(), ".");
                                RATMINEgff3Writer.writeAttributes4Gff3(RATMINEattributesHashMap);
                            }

                            List<edu.mcw.rgd.gff3.CodingFeature> cfList = utils.buildCfList(trMd);
                            for(edu.mcw.rgd.gff3.CodingFeature cf: cfList){
                                String featureId;

                                if(cf.getFeatureType()== TranscriptFeature.FeatureType.CDS){
                                    featureId = getUniqueId("trFeatureCDS"+cf.getRgdId(), idMap);

                                    counters.increment(" Mapped CDSs");
                                }
                                else {
                                    featureId = getUniqueId("trFeature"+cf.getRgdId(), idMap);

                                    if(cf.getFeatureType()== TranscriptFeature.FeatureType.EXON){
                                        counters.increment(" Mapped Exons");
                                    }else
                                    if(cf.getFeatureType()== TranscriptFeature.FeatureType.UTR5){
                                        counters.increment(" Mapped UTR5");
                                    }else
                                    if(cf.getFeatureType()== TranscriptFeature.FeatureType.UTR3){
                                        counters.increment(" Mapped UTR3");
                                    }
                                }

                                gff3Writer.writeFirst8Columns(cf.getChromosome(), "RGD", cf.getCanonicalName(), cf.getStartPos(), cf.getStopPos(), ".", cf.getStrand(), cf.getCodingPhaseStr());


                                attributesHashMap.put("ID", featureId);
                                attributesHashMap.put("Parent", id);
                                if( cf.getNotes()!=null )
                                    attributesHashMap.put("Note", cf.getNotes());

                                RATMINEattributesHashMap.put("ID", featureId);
                                RATMINEattributesHashMap.put("Parent", id);
                                if( cf.getNotes()!=null )
                                    RATMINEattributesHashMap.put("Note", cf.getNotes());


                                gff3Writer.writeAttributes4Gff3(attributesHashMap);

                                if( writeRATMINE ) {
                                    RATMINEgff3Writer.writeFirst8Columns(cf.getChromosome(), "RGD", String.valueOf(cf.getFeatureType()), cf.getStartPos(), cf.getStopPos(), ".", cf.getStrand(), cf.getCodingPhaseStr());
                                    RATMINEgff3Writer.writeAttributes4Gff3(RATMINEattributesHashMap);
                                }
                            }
                        }
                    }
                }
                else{
                    counters.increment(" Genes with NO transcripts");

                    // generate fake feature for genes without features
                    gff3Writer.writeFirst8Columns(map.getChromosome(), "RGD", getSoFeatureType(gType), map.getStartPos(),map.getStopPos(), ".", map.getStrand(), ".");

                    attributesHashMap.put("ID", getUniqueId("ftRGD"+geneRgdId, idMap));
                    attributesHashMap.put("Parent", uniqueGeneId);
                    gff3Writer.writeAttributes4Gff3(attributesHashMap);
                }
            }//end of map data loop
        }

        gff3Writer.close();
        RATMINEgff3Writer.close();


        dumpCounters(counters, assemblySymbol, msgBuf);

        msgBuf.append("OK!  elapsed "+Utils.formatElapsedTime(time0, System.currentTimeMillis())+"\n");
        msgBuf.append("========\n");
        synchronized( this.getClass() ) {
            log.info(msgBuf);
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

    String getSoFeatureType(String geneType) throws Exception {

        if( geneType.startsWith("predicted") )
            return "predicted_gene";

        switch( geneType ) {
            case "pseudogene":
                return "pseudogene";
            case "pseudo":
                return "pseudogene";
            case "protein-coding":
                return "protein_coding_gene";
            case "ncrna":
                return "ncRNA_gene";
            case "gene":
                return "gene";
            case "allele":
                return "allele";
            case "splice":
                return "alternatively_spliced";
            case "rrna":
                return "rRNA_gene";
            case "trna":
                return "tRNA_gene";
            case "snrna":
                return "snRNA_gene";
            case "scrna":
                return "scRNA_gene";
            case "snorna":
                return "snoRNA_gene";
            case "miscrna":
                return "miscRNA_gene";
            case "biological-region":
                return "region";
            default:
                throw new Exception("unsupported gene type " + geneType);
        }
    }

    /**
     * returns a comma separated string of key:value pairs that are external identifiers
     * @param xdbList list of XdbId-s
     * @return xdb id string
     * @throws Exception
     */
    private String getXdbString(List<XdbId> xdbList, CounterPool counters, boolean ratmineCompatible) throws Exception{

        Set<String> xdbIds = new TreeSet<>();
        for(XdbId externalId : xdbList ){
            if(externalId.getAccId()!=null){
                String xdbEntry;

                if(externalId.getXdbKeyAsString().contains(" ")){
                    xdbEntry = externalId.getXdbKeyAsString().replaceAll(" ","")+":";
                }else{
                    xdbEntry = externalId.getXdbKeyAsString()+":";
                }

                if(externalId.getAccId().contains(":")){
                    if( ratmineCompatible && externalId.getXdbKey()==XdbId.XDB_KEY_MGD ) {
                        xdbEntry += externalId.getAccId();
                    } else {
                        xdbEntry += externalId.getAccId().replaceAll(":", "");
                    }
                }else{
                    xdbEntry += externalId.getAccId();
                }

                xdbIds.add(xdbEntry);
            }
        }
        String extDbString = Utils.concatenate(xdbIds,",");

        if( !extDbString.isEmpty() ){
            if( extDbString.contains("NCBIGene") )
                counters.increment(" Genes with NCBI geneIds");
        }else{
            counters.increment(" Genes with NO XDB Ids");
        }

        return extDbString;
    }

    /**
     * get comma separated and comma prefixed list of MGI or HGNC ids, f.e. ',HGNC:2';
     * this string is to be directly appended to 'Aliases' field
     * @param xdbList list of XdbId-s
     * @return MGI/HGNC ids string
     * @throws Exception
     */
    private String getHgncMgiIds(List<XdbId> xdbList) throws Exception{

        StringBuilder xdbIds = new StringBuilder();
        for(XdbId externalId : xdbList ){
            if(externalId.getXdbKey()==XdbId.XDB_KEY_MGD || externalId.getXdbKey()==XdbId.XDB_KEY_HGNC) {
                xdbIds.append(",").append(externalId.getAccId());
            }
        }
        return xdbIds.toString();
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
        Iterator<XdbId> it = xdbList.iterator();
        while( it.hasNext() ) {
            XdbId externalId = it.next();

            if( externalId.getXdbKey()!=XdbId.XDB_KEY_NCBI_GENE &&
                    externalId.getXdbKey()!=XdbId.XDB_KEY_MGD &&
                    externalId.getXdbKey()!=XdbId.XDB_KEY_OMIM &&
                    externalId.getXdbKey()!=XdbId.XDB_KEY_UNIPROT &&
                    externalId.getXdbKey()!=XdbId.XDB_KEY_HGNC &&
                    externalId.getXdbKey()!=XdbId.XDB_KEY_ENSEMBL_GENES &&
                    externalId.getXdbKey()!=XdbId.XDB_KEY_ENSEMBL_TRANSCRIPT &&
                    true){
                it.remove();
            }
        }
        return xdbList;
    }

    private String getNameOfGene(Gene gene, CounterPool counters) {

        String nameOfgene = gene.getName();

        if(nameOfgene!=null){
            if(nameOfgene.contains(",")){
                nameOfgene = nameOfgene.replaceAll(",", "");
            }
            if(nameOfgene.contains(";")){
                nameOfgene = nameOfgene.replaceAll(";", "");
            }
        }

        if(nameOfgene==null){
            nameOfgene="null";
            counters.increment(" Genes with NULL geneName");
        }
        return nameOfgene;
    }

    String getUniqueId(String idBase, Map<String, AtomicInteger> idMap) {

        AtomicInteger i = idMap.putIfAbsent(idBase, new AtomicInteger(0));
        if( i==null ) {
            i = idMap.get(idBase);
        }
        int cnt = i.incrementAndGet();

        if( cnt==1 )
            return idBase;
        return idBase+"_"+cnt;
    }

    public void setProcessedAssemblies(List<String> processedAssemblies) {
        this.processedAssemblies = processedAssemblies;
    }

    public List<String> getProcessedAssemblies() {
        return processedAssemblies;
    }
}
