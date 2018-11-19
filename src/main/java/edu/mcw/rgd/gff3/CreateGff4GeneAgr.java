package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;

import java.util.*;
import java.util.Map;

/**
 * @author mtutaj
 * @since 9/19/2017
 */
public class CreateGff4GeneAgr {

    private String gff3Path;
    RgdGff3Dao dao = new RgdGff3Dao();

    private int speciesTypeKey;
    private int mapKey;

    public void setGff3Path(String gff3Path) {
        this.gff3Path = gff3Path;
    }

    Map<String, Integer> idMap = new HashMap<>();

    public void createGeneGff3(boolean compress) throws Exception{

        CdsUtils utils = new CdsUtils(dao, mapKey);

        String species = SpeciesType.getCommonName(speciesTypeKey);
        idMap.clear();

        Counters counters = new Counters();

        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gff3Path+species+"_RGD_AGR.gff3", false, compress);
        gff3Writer.setAgrCompatibleFormat(true);

        gff3Writer.print("# Genome build: "+ MapManager.getInstance().getMap(mapKey).getName()+"\n");
        gff3Writer.print("# Primary Contact: mtutaj@mcw.edu\n");
        gff3Writer.print("# Tool: AGR GFF3 extractor v.1.1.3 (single locus genes)\n");
        gff3Writer.print("# Generated: "+new Date()+"\n");

        List<Gene> activeGenes = dao.getActiveGenes(speciesTypeKey);
        Collections.sort(activeGenes, new Comparator<Gene>() {
            @Override
            public int compare(Gene o1, Gene o2) {
                return Utils.stringsCompareToIgnoreCase(o1.getSymbol(), o2.getSymbol());
            }
        });
        //Collections.shuffle(activeGenes);

        for( Gene gene: activeGenes ){
            int geneRgdId = gene.getRgdId();

            // process only protein coding genes
            String gType = Utils.NVL(gene.getType(), "gene");

            // 'curie' is unique identifier used by AGR;
            // for rat genes it is 'RGD:xxx', for human genes it is 'HGNC:xxx'
            List<XdbId> xdbIds = getXdbIds(geneRgdId);
            String curie = getCurie(geneRgdId, xdbIds);
            if( curie==null ) {
                continue;
            }

            List<MapData> geneMap = getMapData(geneRgdId, mapKey, counters);
            if( geneMap.isEmpty() || geneMap.size()>1 ) {
                //System.out.println("no map positions");
                continue;
            }
            MapData map = geneMap.get(0);

            counters.genesInEachChromosomeCount++;

            String nameOfgene = getNameOfGene(gene, counters);

            List<Transcript> geneTrs = dao.getTranscriptsForGene(geneRgdId);
            if(geneTrs.size()>1){
                counters.genesMoreThanOneMappedTrans++;
            }

            String annotDesc = Utils.getGeneDescription(gene);
            if( Utils.isStringEmpty(annotDesc) ){
                annotDesc = null;
            }else{
                annotDesc = annotDesc.replaceAll(";"," AND ").replaceAll(",", " ");
            }

            String nc="N";

            Map<String,String> attributes = new HashMap<>();

            String uniqueGeneId = getUniqueId("gene");
            attributes.put("ID", uniqueGeneId);
            attributes.put("Name", gene.getSymbol());
            attributes.put("Note", nameOfgene);
            attributes.put("curie", curie);

            String alias = gene.getSymbol()+","+curie+ "," + nameOfgene;
            if( gene.getSpeciesTypeKey()==SpeciesType.HUMAN ) {
                alias += ",RGD:"+gene.getRgdId();
            }
            attributes.put("Alias", alias);

            attributes.put("geneType", gType.replaceAll("\\-","_"));
            if( gene.getRefSeqStatus()!=null )
                attributes.put("status",gene.getRefSeqStatus());
            if( annotDesc!=null )
                attributes.put("description",annotDesc);
            if( !Utils.isStringEmpty(gene.getNcbiAnnotStatus()) ) {
                attributes.put("nomenclatureStatus", gene.getNcbiAnnotStatus());
            }
            String extDbString = getXdbString(xdbIds, counters);
            if( !extDbString.isEmpty() ){
                attributes.put("Dbxref",extDbString);
            }
            attributes.put("so_term_name", getSoTermNameForGene(gType));

            gff3Writer.writeFirst8Columns(map.getChromosome(),"RGD", "gene", map.getStartPos(),map.getStopPos(),".",map.getStrand(),".");
            gff3Writer.writeAttributes4Gff3(attributes);

            if(geneTrs.size()>0){
                counters.transcriptsCount++;
                counters.genesWithTranscriptsCount++;

                for( Transcript tr: geneTrs ){

                    if(tr.isNonCoding()){
                        nc="Y";
                        counters.nonCodingTransCount++;
                    }
                    for( MapData trMd: tr.getGenomicPositions() ) {
                        if( !utils.transcriptPositionOverlapsGenePosition(trMd, map) )
                            continue;

                        String id = getUniqueId("rna");
                        counters.transcriptsMappedCount++;


                        attributes.put("ID", id);
                        attributes.put("Name", tr.getAccId());
                        attributes.put("Parent", uniqueGeneId);
                        if( tr.getRefSeqStatus()!=null )
                            attributes.put("status", tr.getRefSeqStatus());
                        //attributes.put("isNonCoding", nc);
                        attributes.put("gene", gene.getSymbol());
                        attributes.put("biotype", getTrBiotype(gType, gene.getName(), tr));

                        gff3Writer.writeFirst8Columns(trMd.getChromosome(), "RGD", "mRNA", trMd.getStartPos(),trMd.getStopPos(), ".", trMd.getStrand(), ".");
                        gff3Writer.writeAttributes4Gff3(attributes);

                        String cdsId = null;
                        List<CodingFeature> cfList = utils.buildCfList(trMd);
                        for(CodingFeature cf: cfList){
                            String featureId;

                            if(cf.getFeatureType()== TranscriptFeature.FeatureType.CDS){
                                // one CDS id per multiple fragments of CDS (as it is in NCBI RefSeq gff3 file for rat)
                                featureId = cdsId==null ? getUniqueId("cds") : cdsId;

                                counters.cdsCount++;
                            }
                            else {
                                if(cf.getFeatureType()== TranscriptFeature.FeatureType.EXON){
                                    counters.exonCount++;
                                    featureId = getUniqueId("e");
                                }else
                                if(cf.getFeatureType()== TranscriptFeature.FeatureType.UTR5){
                                    counters.utr5Count++;
                                    featureId = getUniqueId("u");
                                }else
                                if(cf.getFeatureType()== TranscriptFeature.FeatureType.UTR3){
                                    counters.utr3Count++;
                                    featureId = getUniqueId("u");
                                } else {
                                    featureId = getUniqueId("f");
                                }
                            }

                            gff3Writer.writeFirst8Columns(cf.getChromosome(), "RGD", cf.getCanonicalName(), cf.getStartPos(), cf.getStopPos(), ".", cf.getStrand(), cf.getCodingPhaseStr());


                            attributes.put("ID", featureId);
                            attributes.put("Parent", id);
                            if( cf.getNotes()!=null )
                                attributes.put("Note", cf.getNotes());

                            gff3Writer.writeAttributes4Gff3(attributes);
                        }
                    }
                }
            } else {
                counters.genesWithNoTranscripts++;

                // generate fake feature for genes without features
                gff3Writer.writeFirst8Columns(map.getChromosome(), "RGD", "transcript_region", map.getStartPos(),map.getStopPos(), ".", map.getStrand(), ".");

                String regionId = getUniqueId("rna");
                attributes.put("ID", regionId);
                attributes.put("Parent", uniqueGeneId);
                gff3Writer.writeAttributes4Gff3(attributes);

                gff3Writer.writeFirst8Columns(map.getChromosome(), "RGD", "exon", map.getStartPos(),map.getStopPos(), ".", map.getStrand(), ".");

                attributes.put("ID", getUniqueId("e"));
                attributes.put("Parent", regionId);
                gff3Writer.writeAttributes4Gff3(attributes);
            }

        }

        gff3Writer.close();

        dumpCounters(counters);
    }

    String getTrBiotype(String geneType, String geneName, Transcript tr) {
        switch(geneType) {
            case "ncrna":
                if( geneName.startsWith("microRNA") ) {
                    return "miRNA";
                } else {
                    return "misc_RNA";
                }
            case "protein-coding":
                if( tr.getAccId().startsWith("XR_") || tr.getAccId().startsWith("NR_") ) {
                    return "misc_RNA";
                }
                return "protein_coding";
            case "pseudo": return "pseudogene";
            case "rrna": return "rRNA";
            case "snrna": return "snRNA";
            case "gene":
                if( tr.getAccId().startsWith("XR_") ) {
                    return "misc_RNA";
                }
                if( tr.getAccId().startsWith("XM_") ) {
                    return "protein_coding";
                }
                System.out.println("unknown");
            default:
                System.out.println("unknown");
                return "";
        }
    }

    String getCurie(int geneRgdId, List<XdbId> xdbIds) {
        if (speciesTypeKey == SpeciesType.RAT) {
            return "RGD:" + geneRgdId;
        } else if (speciesTypeKey == SpeciesType.HUMAN ) {
            for(XdbId id : xdbIds ){
                if(id.getXdbKey()==XdbId.XDB_KEY_HGNC) {
                    return id.getAccId();
                }
            }
        }
        return null;
    }

    List<MapData> getMapData(int geneRgdId, int mapKey, Counters counters) throws Exception {

        List<MapData> geneMap = dao.getMapData(geneRgdId, mapKey);
        if(geneMap.size()>1){
            counters.genesMoreThanOneMapPosCurrentAssembly++;
        }
        return geneMap;
    }

    void dumpCounters(Counters counters) {

        System.out.println(" Genes processed:" + counters.genesInEachChromosomeCount);
        if( counters.genesMoreThanOneMapPosCurrentAssembly>0 )
            System.out.println(" Genes with more than one map position:" + counters.genesMoreThanOneMapPosCurrentAssembly);
        if( counters.genesGeneNameNull>0 )
            System.out.println(" Genes with NULL geneName:"+ counters.genesGeneNameNull);
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
        System.out.println("GFF3 File SUCCESSFUL!");
        System.out.println("=============================");
    }

    String getSoTermNameForGene(String geneType) throws Exception {
        if( geneType.equals("protein-coding") )
            return "protein_coding_gene";
        if( geneType.equals("pseudogene") )
            return "pseudogene";
        if( geneType.equals("pseudo") )
            return "pseudogene";
        if( geneType.equals("ncrna") )
            return "ncRNA_gene";
        if( geneType.equals("trna") )
            return "tRNA_gene";
        if( geneType.equals("rrna") )
            return "rRNA_gene";
        if( geneType.equals("gene") )
            return "gene";
        if( geneType.equals("snrna") )
            return "snRNA_gene";
        if( geneType.equals("snorna") )
            return "snoRNA_gene";
        if( geneType.equals("scrna") )
            return "scRNA_gene";
        /*
        if( geneType.equals("allele") )
            return "allele";
        if( geneType.equals("splice") )
            return "alternatively_spliced";
        if( geneType.startsWith("predicted") )
            return "predicted_gene";
        if( geneType.equals("miscrna") )
            return "miscRNA_gene";
            */
        throw new Exception("unsupported gene type "+geneType);
    }

    /**
     * returns a comma separated string of key:value pairs that are external identifiers
     * @param xdbList list of XdbId-s
     * @return xdb id string
     * @throws Exception
     */
    private String getXdbString(List<XdbId> xdbList, Counters counters) throws Exception{

        Set<String> xdbIds = new TreeSet<>();
        for(XdbId externalId : xdbList ){
            if(externalId.getAccId()!=null){
                String xdbEntry;

                if(externalId.getAccId().contains(":")){
                    xdbEntry = externalId.getAccId();
                } else {
                    String xdbName = externalId.getXdbKeyAsString();
                    if( xdbName.contains(" ")){
                        xdbEntry = xdbName.replaceAll(" ","");
                    }else{
                        xdbEntry = xdbName;
                    }
                    xdbEntry += ":" + externalId.getAccId();
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

        return idBase+cnt;
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

    class Counters {
        int genesInEachChromosomeCount=0;
        int genesWithTranscriptsCount=0;
        int genesWithNoTranscripts=0;
        int transcriptsMappedCount=0;
        int genesMoreThanOneMappedTrans=0;
        int nonCodingTransCount=0;
        int genesMoreThanOneMapPosCurrentAssembly=0;
        int genesGeneNameNull=0;
        int genesWithNcbiGeneIds=0;
        int genesWithNOXdbIds=0;

        int transcriptsCount=0;
        int exonCount=0;
        int utr5Count=0;
        int utr3Count=0;
        int cdsCount=0;
    }
}
