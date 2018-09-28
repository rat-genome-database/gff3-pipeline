package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;

import java.util.*;
import java.util.Map;

public class CreateGff4Gene {

    private String newPathGff3;
    RgdGff3Dao dao = new RgdGff3Dao();

    private int speciesTypeKey;
    private int mapKey;
    private List<String> chromosomes;

    public void setNewPathGff3(String newPathGff3) {
        this.newPathGff3 = newPathGff3;
    }

    public String getNewPathGff3() {
        return newPathGff3;
    }

    Map<String, Integer> idMap = new HashMap<>();

    public void createGeneGff3(boolean compress) throws Exception{

        CdsUtils utils = new CdsUtils(dao, mapKey);

        Map<String, Integer> chrMap = dao.getChromosomeSizes(mapKey);
        if( chrMap.isEmpty() ) {
            System.out.println("*** WARNING: no chromosome sizes available for MAP_KEY="+mapKey);
        }
        String species = SpeciesType.getCommonName(speciesTypeKey);
        idMap.clear();

        // March 2 2016: RATMINE gff3 loader does not allow one gene to have multiple loci
        //               so we emit only first loci per gene
        Set<Integer> geneRgdIdsEmitted = new HashSet<>();

        for(String chr: getChromosomes()){

            Counters counters = new Counters();

            Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(newPathGff3+species+"_RGDChr"+chr+".gff3", false, compress);
            Gff3ColumnWriter RATMINEgff3Writer = new Gff3ColumnWriter(newPathGff3+species+"_RATMINE_RGDChr"+chr+".gff3", false, compress);
            RATMINEgff3Writer.setRatmineCompatibleFormat(true);

            Integer chrSize = chrMap.get(chr);
            List<Gene> activeGenes = chr.equals("Scaffold") ? dao.getActiveGenes(speciesTypeKey) : dao.getActiveGenes(chr, 1, chrSize, mapKey);

            //Gene gene = geneDao.getGene(2413);
            //if(gene!=null){
            for( Gene gene: activeGenes ){
                int geneRgdId = gene.getRgdId();

                List<MapData> geneMap = getMapData(geneRgdId, mapKey, chr, counters);
                if( geneMap.isEmpty() ) {
                    //System.out.println("no map positions");
                    continue;
                }
                counters.genesInEachChromosomeCount++;

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
                        counters.geneLociSkippedForRatMine++;
                    }

                    List<Transcript> trsOnMap = getTranscriptsForMap(geneTrs, map, utils);

                    String gType;
                    if(gene.getType()==null){
                        gType="gene";
                        counters.genesGeneTypeNull++;

                    }else{
                        gType=gene.getType();

                        if(gene.getType().contains("pseudo")){
                            counters.pseudoGenesCount++;
                        }
                    }

                    List<XdbId> xdbIds = getXdbIds(geneRgdId);

                    Map<String,String> attributesHashMap = new HashMap<>();
                    Map<String,String> RATMINEattributesHashMap = new HashMap<>();

                    String uniqueGeneId = getUniqueId("RGD"+gene.getRgdId());
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

                    String RATMINEuniqueGeneId = getUniqueId("RGD:"+gene.getRgdId());
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
                            System.out.println("RATMINE error: dash in ID");
                        }
                        RATMINEgff3Writer.writeFirst8Columns(map.getChromosome(), "RGD", "gene", map.getStartPos(), map.getStopPos(), ".", map.getStrand(), ".");
                        RATMINEgff3Writer.writeAttributes4Gff3(RATMINEattributesHashMap);
                    }

                    if(trsOnMap.size()>1){
                        counters.genesMoreThanOneMappedTrans++;
                    }
                    if(trsOnMap.size()>0){
                        counters.transcriptsCount++;
                        counters.genesWithTranscriptsCount++;

                        for( Transcript tr: trsOnMap ){

                            if(tr.isNonCoding()){
                                nc="Y";
                                counters.nonCodingTransCount++;
                            }
                            for( MapData trMd: tr.getGenomicPositions() ) {
                                if( !utils.transcriptPositionOverlapsGenePosition(trMd, map) )
                                    continue;

                                String id = getUniqueId("mRNARGD"+tr.getRgdId());

                                counters.transcriptsMappedCount++;


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
                                        featureId = getUniqueId("trFeatureCDS"+cf.getRgdId());

                                        counters.cdsCount++;
                                    }
                                    else {
                                        featureId = getUniqueId("trFeature"+cf.getRgdId());

                                        if(cf.getFeatureType()== TranscriptFeature.FeatureType.EXON){
                                            counters.exonCount++;
                                        }else
                                        if(cf.getFeatureType()== TranscriptFeature.FeatureType.UTR5){
                                            counters.utr5Count++;
                                        }else
                                        if(cf.getFeatureType()== TranscriptFeature.FeatureType.UTR3){
                                            counters.utr3Count++;
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
                        counters.genesWithNoTranscripts++;

                        // generate fake feature for genes without features
                        gff3Writer.writeFirst8Columns(map.getChromosome(), "RGD", getSoFeatureType(gType), map.getStartPos(),map.getStopPos(), ".", map.getStrand(), ".");

                        attributesHashMap.put("ID", getUniqueId("ftRGD"+geneRgdId));
                        attributesHashMap.put("Parent", uniqueGeneId);
                        gff3Writer.writeAttributes4Gff3(attributesHashMap);
                    }
                }//end of map data loop
            }

            gff3Writer.close();
            RATMINEgff3Writer.close();

            dumpCounters(counters, chr);
        }
    }

    // get a subset of transcripts having position on a given map, and overlapping the given gene position
    List<Transcript> getTranscriptsForMap(List<Transcript> geneTrs, MapData md, CdsUtils utils) {

        List<Transcript> trsOnMap = new ArrayList<>();
        for( Transcript tr: geneTrs ) {
            for( MapData trMd: tr.getGenomicPositions() ) {
                // transcript pos must overlap gene locus
                if( utils.transcriptPositionOverlapsGenePosition(trMd, md) ) {
                    // yes, this transcript is overlapping gene locus
                    trsOnMap.add(tr);
                }
            }
        }
        return trsOnMap;
    }

    List<MapData> getMapData(int geneRgdId, int mapKey, String chr, Counters counters) throws Exception {

        List<MapData> geneMap = dao.getMapData(geneRgdId, mapKey);

        // filter out duplicates (when a gene has loci on multiple chromosomes)
        if( !chr.equals("Scaffold") ) {
            Iterator<MapData> it = geneMap.iterator();
            while (it.hasNext()) {
                if (!it.next().getChromosome().equals(chr))
                    it.remove();
            }
        }

        if( !geneMap.isEmpty() ) {
            if(geneMap.size()>1){
                counters.genesMoreThanOneMapPosCurrentAssembly++;
            }
        }
        return geneMap;
    }

    void dumpCounters(Counters counters, String chr) {

        System.out.println("\nCounts For: Chr" + chr);
        System.out.println(" Genes processed:" + counters.genesInEachChromosomeCount);
        if( counters.genesMoreThanOneMapPosCurrentAssembly>0 )
            System.out.println(" Genes with more than one map position:" + counters.genesMoreThanOneMapPosCurrentAssembly);
        if( counters.genesGeneNameNull>0 )
            System.out.println(" Genes with NULL geneName:"+ counters.genesGeneNameNull);
        if( counters.genesGeneTypeNull>0 )
            System.out.println(" Genes with NULL geneType:" + counters.genesGeneTypeNull);
        if( counters.pseudoGenesCount>0 )
            System.out.println(" Pseudo Genes: " + counters.pseudoGenesCount);
        System.out.println(" Genes with transcripts:" + counters.genesWithTranscriptsCount);
        if( counters.genesWithNoTranscripts>0 )
            System.out.println(" Genes with NO transcripts:" + counters.genesWithNoTranscripts);
        System.out.println(" Mapped Transcripts:" + counters.transcriptsMappedCount);
        System.out.println(" Mapped Exons:" + counters.exonCount);
        System.out.println(" Mapped CDSs:" + counters.cdsCount);
        System.out.println(" Genes with more than one mapped transcript:" + counters.genesMoreThanOneMappedTrans);
        System.out.println(" NonCoding transcripts:" + counters.nonCodingTransCount);
        System.out.println(" Genes with NCBI geneIds:" + counters.genesWithNcbiGeneIds);
        if( counters.genesWithNOXdbIds>0 )
            System.out.println(" Genes with NO XDB Ids:" + counters.genesWithNOXdbIds);
        if( counters.geneLociSkippedForRatMine>0 )
            System.out.println(" Gene Loci skipped for RATMINE:" + counters.geneLociSkippedForRatMine);
        System.out.println("DONE with Chromosome:" + chr+"\n");
        System.out.println("GFF3 File SUCCESSFUL!");
        System.out.println("=============================");
    }

    String getSoFeatureType(String geneType) throws Exception {
        if( geneType.equals("pseudogene") )
            return "pseudogene";
        if( geneType.equals("pseudo") )
            return "pseudogene";
        if( geneType.equals("protein-coding") )
            return "protein_coding_gene";
        if( geneType.equals("ncrna") )
            return "ncRNA_gene";
        if( geneType.equals("gene") )
            return "gene";
        if( geneType.equals("allele") )
            return "allele";
        if( geneType.equals("splice") )
            return "alternatively_spliced";
        if( geneType.startsWith("predicted") )
            return "predicted_gene";
        if( geneType.equals("rrna") )
            return "rRNA_gene";
        if( geneType.equals("trna") )
            return "tRNA_gene";
        if( geneType.equals("snrna") )
            return "snRNA_gene";
        if( geneType.equals("scrna") )
            return "scRNA_gene";
        if( geneType.equals("snorna") )
            return "snoRNA_gene";
        if( geneType.equals("miscrna") )
            return "miscRNA_gene";
        throw new Exception("unsupported gene type "+geneType);
    }

    /**
     * returns a comma separated string of key:value pairs that are external identifiers
     * @param xdbList list of XdbId-s
     * @return xdb id string
     * @throws Exception
     */
    private String getXdbString(List<XdbId> xdbList, Counters counters, boolean ratmineCompatible) throws Exception{

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
                counters.genesWithNcbiGeneIds++;
        }else{
            counters.genesWithNOXdbIds++;
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

        String xdbIds = "";
        for(XdbId externalId : xdbList ){
            if(externalId.getXdbKey()==XdbId.XDB_KEY_MGD || externalId.getXdbKey()==XdbId.XDB_KEY_HGNC) {
                xdbIds += "," + externalId.getAccId();
            }
        }
        return xdbIds;
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
                    externalId.getXdbKey()!=XdbId.XDB_KEY_UNIGENE &&
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

    private String getNameOfGene(Gene gene, Counters counters) {

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
            counters.genesGeneNameNull++;
        }
        return nameOfgene;
    }

    String getUniqueId(String idBase) {

        Integer cnt = idMap.get(idBase);
        if( cnt==null ) {
            cnt = 1;
        }
        else {
            cnt++;
        }
        idMap.put(idBase, cnt);

        if( cnt==1 )
            return idBase;
        return idBase+"_"+cnt;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setMapKey(int mapKey) {
        this.mapKey = mapKey;
    }

    public int getMapKey() {
        return mapKey;
    }

    public List<String> getChromosomes() {
        return chromosomes;
    }

    public void setChromosomes(List<String> chromosomes) {
        this.chromosomes = chromosomes;
    }

    class Counters {
        int genesInEachChromosomeCount=0;
        int genesWithTranscriptsCount=0;
        int genesWithNoTranscripts=0;
        int transcriptsMappedCount=0;
        int genesMoreThanOneMappedTrans=0;
        int pseudoGenesCount=0;
        int nonCodingTransCount=0;
        int genesMoreThanOneMapPosCurrentAssembly=0;
        int genesGeneTypeNull=0;
        int genesGeneNameNull=0;
        int genesWithNcbiGeneIds=0;
        int genesWithNOXdbIds=0;
        int geneLociSkippedForRatMine=0;

        int transcriptsCount=0;
        int exonCount=0;
        int utr5Count=0;
        int utr3Count=0;
        int cdsCount=0;
    }
}
