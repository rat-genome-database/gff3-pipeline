package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;

import java.text.SimpleDateFormat;
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
    private int mapKeyEnsembl;

    public void setGff3Path(String gff3Path) {
        this.gff3Path = gff3Path;
    }

    Map<String, Integer> idMap = new HashMap<>();

    public void createGeneGff3(int compressMode) throws Exception{

        CdsUtils utilsNcbi = new CdsUtils(dao, mapKey);
        CdsUtils utilsEnsembl = new CdsUtils(dao, mapKeyEnsembl);

        String species = SpeciesType.getCommonName(speciesTypeKey);
        idMap.clear();

        Map<Integer, String> dbXrefsForTranscripts = loadDbXrefsForTranscripts(speciesTypeKey);

        Counters counters = new Counters();

        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gff3Path+species+"_RGD_AGR.gff3", false, compressMode);
        gff3Writer.setAgrCompatibleFormat(true);

        // date format as agreed on DQM meeting on May 18, 2021
        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'");

        gff3Writer.print("#!data-source RAT GENOME DATABASE (https://rgd.mcw.edu/)\n");
        gff3Writer.print("#!assembly: "+ MapManager.getInstance().getMap(mapKey).getName()+"\n");
        if( mapKey==372 ) {
            gff3Writer.print("#!annotationSource RefSeq RS_2023_06\n");     // https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_015227675.2/
            gff3Writer.print("#!annotationSource ENSEMBL 110.72\n"); // https://m.ensembl.org/Rattus_norvegicus/Info/Annotation
        } else if( mapKey==38 ) {
            gff3Writer.print("#!annotationSource RefSeq RS_2023_03 (GRCh38.p14)\n");     // https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40
            gff3Writer.print("#!annotationSource ENSEMBL 110.38 (GRCh38.p14)\n"); // https://m.ensembl.org/Homo_sapiens/Info/Annotation
        }

        gff3Writer.print("#!date-produced "+sdt.format(new Date())+"\n");
        gff3Writer.print("#!species "+ species+"\n");
        gff3Writer.print("#!primary-contact mtutaj@mcw.edu\n");
        gff3Writer.print("#!tool AGR GFF3 extractor  v 2023-08-08\n");

        List<Gene> activeGenes = dao.getActiveGenes(speciesTypeKey);
        Collections.sort(activeGenes, new Comparator<Gene>() {
            @Override
            public int compare(Gene o1, Gene o2) {
                return Utils.stringsCompareToIgnoreCase(o1.getSymbol(), o2.getSymbol());
            }
        });

        System.out.println("active genes: "+activeGenes.size());
//Collections.shuffle(activeGenes);
        int i = 0;
        for( Gene gene: activeGenes ){
            i++;
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

            List<MapData> geneMap = getMergedMapData(geneRgdId, mapKey, mapKeyEnsembl, counters);
            if( geneMap.isEmpty() ) {
                continue;
            }

            for( MapData map: geneMap ) {
                counters.genesInEachChromosomeCount++;

                String nameOfgene = getNameOfGene(gene, counters);

                String annotDesc = Utils.getGeneDescription(gene);
                if (Utils.isStringEmpty(annotDesc)) {
                    annotDesc = null;
                } else {
                    annotDesc = annotDesc.replaceAll(";", " AND ").replaceAll(",", " ");
                }

                Map<String, String> attributes = new HashMap<>();

                String uniqueGeneId = getUniqueId2("RGD:" + geneRgdId);
                attributes.put("ID", uniqueGeneId);
                attributes.put("Name", gene.getSymbol());
                if( nameOfgene!=null ) {
                    attributes.put("Note", nameOfgene);
                }
                attributes.put("curie", curie);
                attributes.put("gene_id", curie);

                String alias = gene.getSymbol() + "," + curie;
                if( nameOfgene!=null ) {
                    alias += "," + nameOfgene;
                }
                if (gene.getSpeciesTypeKey() == SpeciesType.HUMAN) {
                    alias += ",RGD:" + gene.getRgdId();
                }
                attributes.put("Alias", alias);

                attributes.put("geneType", gType.replaceAll("\\-", "_"));
                if (gene.getRefSeqStatus() != null)
                    attributes.put("status", gene.getRefSeqStatus());
                if (annotDesc != null)
                    attributes.put("description", annotDesc);
                if (!Utils.isStringEmpty(gene.getNcbiAnnotStatus())) {
                    attributes.put("nomenclatureStatus", gene.getNcbiAnnotStatus());
                }
                String extDbString = getXdbString(xdbIds, counters);
                if (!extDbString.isEmpty()) {
                    attributes.put("Dbxref", extDbString);
                }
                attributes.put("so_term_name", getSoTermNameForGene(gType));
                attributes.put("Ontology_term", getSoTermAccForGene(gType));

                gff3Writer.writeFirst8Columns(map.getChromosome(), "RGD", "gene", map.getStartPos(), map.getStopPos(), ".", map.getStrand(), ".");
                gff3Writer.writeAttributes4Gff3(attributes);


                List<Transcript> geneTrs = dao.getTranscriptsForGene(geneRgdId);
                if (geneTrs.size() > 1) {
                    counters.genesMoreThanOneMappedTrans++;
                }

                if (geneTrs.size() > 0) {
                    counters.transcriptsCount++;
                    counters.genesWithTranscriptsCount++;

                    for (Transcript tr : geneTrs) {

                        if (tr.isNonCoding()) {
                            counters.nonCodingTransCount++;
                        }
                        for (MapData trMd : tr.getGenomicPositions()) {
                            // skip positions from other assemblies
                            if( !(trMd.getMapKey()==mapKey || trMd.getMapKey()==mapKeyEnsembl) ) {
                                continue;
                            }

                            boolean isEnsemblTr = trMd.getMapKey()==mapKeyEnsembl;
                            String dbName = isEnsemblTr ? "ENSEMBL" : "NCBI";
                            CdsUtils utils = isEnsemblTr ? utilsEnsembl : utilsNcbi;

                            if( !positionsOverlap(trMd, map))
                                continue;

                            String id = getUniqueId("rna");
                            counters.transcriptsMappedCount++;

                            attributes.put("ID", id);
                            String trAccVer = dao.getTranscriptVersionInfo(tr.getAccId());
                            if( trAccVer==null ) {
                                System.out.println("no tr acc ver for "+tr.getAccId());
                                trAccVer = tr.getAccId();
                            }
                            String transcriptId = isEnsemblTr ? "ENSEMBL:" : "RefSeq:";
                            transcriptId += trAccVer;
                            attributes.put("transcript_id", transcriptId);

                            // transcript curie: acc id prefixed with either 'RefSeq:' or 'Ensembl:'
                            attributes.put("curie", transcriptId);

                            attributes.put("Name", trAccVer);
                            attributes.put("Parent", uniqueGeneId);
                            if (tr.getRefSeqStatus() != null) {
                                attributes.put("status", tr.getRefSeqStatus());
                            }
                            attributes.put("gene", gene.getSymbol());
                            String proteinId = null;
                            if (!Utils.isStringEmpty(tr.getProteinAccId())) {
                                proteinId = tr.getProteinAccId();
                                attributes.put("protein_id", proteinId);
                            }
                            attributes.put("Ontology_term", "SO:0000673"); // transcript

                            String trDbXrefs = dbXrefsForTranscripts.get(tr.getRgdId());
                            if (trDbXrefs != null) {
                                attributes.put("Dbxref", trDbXrefs);
                            }

                            String trBiotype = getTrBiotype(gType, tr);
                            gff3Writer.writeFirst8Columns(trMd.getChromosome(), dbName, trBiotype, trMd.getStartPos(), trMd.getStopPos(), ".", trMd.getStrand(), ".");
                            gff3Writer.writeAttributes4Gff3(attributes);

                            boolean trIsCoding = trBiotype.equals("mRNA");

                            boolean hasCds = false;
                            List<CodingFeature> cfList = utils.buildCfList(trMd);
                            for (CodingFeature cf : cfList) {
                                String featureId;

                                if (cf.getFeatureType() == TranscriptFeature.FeatureType.CDS) {

                                    if( !trIsCoding ) {
                                        // CdsUtils always tries to generate CDS; therefore it is important to know if the transcript is coding
                                        continue;
                                    }

                                    // one CDS id per multiple fragments of CDS (as it is in NCBI RefSeq gff3 file for rat)
                                    featureId = getUniqueId("cds");

                                    if( proteinId!=null ) {
                                        attributes.put("protein_id", proteinId);
                                    }

                                    counters.cdsCount++;
                                    hasCds = true;
                                } else {
                                    if (cf.getFeatureType() == TranscriptFeature.FeatureType.EXON) {
                                        counters.exonCount++;
                                        featureId = getUniqueId("e");
                                    } else if (cf.getFeatureType() == TranscriptFeature.FeatureType.UTR5) {
                                        counters.utr5Count++;
                                        featureId = getUniqueId("u");
                                    } else if (cf.getFeatureType() == TranscriptFeature.FeatureType.UTR3) {
                                        counters.utr3Count++;
                                        featureId = getUniqueId("u");
                                    } else {
                                        featureId = getUniqueId("f");
                                    }
                                }

                                gff3Writer.writeFirst8Columns(cf.getChromosome(), dbName, cf.getCanonicalName(), cf.getStartPos(), cf.getStopPos(), ".", cf.getStrand(), cf.getCodingPhaseStr());


                                attributes.put("ID", featureId);
                                attributes.put("Parent", id);
                                if (cf.getNotes() != null)
                                    attributes.put("Note", cf.getNotes());

                                gff3Writer.writeAttributes4Gff3(attributes);
                            }

                            if( !hasCds && trIsCoding ) {
                                System.out.println("hasNoCDS && type="+trBiotype+"    gene nr="+i);
                            }
                        }
                    }
                } else {
                    counters.genesWithNoTranscripts++;

                    // generate fake feature for genes without features
                    gff3Writer.writeFirst8Columns(map.getChromosome(), "RGD", "transcript_region", map.getStartPos(), map.getStopPos(), ".", map.getStrand(), ".");

                    String regionId = getUniqueId("rna");
                    attributes.put("ID", regionId);
                    attributes.put("Parent", uniqueGeneId);
                    gff3Writer.writeAttributes4Gff3(attributes);

                    gff3Writer.writeFirst8Columns(map.getChromosome(), "RGD", "exon", map.getStartPos(), map.getStopPos(), ".", map.getStrand(), ".");

                    attributes.put("ID", getUniqueId("e"));
                    attributes.put("Parent", regionId);
                    gff3Writer.writeAttributes4Gff3(attributes);
                }
            }
        }

        gff3Writer.close();

        dumpCounters(counters);
    }

    String getTrBiotype(String geneType, Transcript tr) throws Exception {
        if( tr.getType()==null ) {
            switch (geneType) {
                case "protein-coding":
                    if (tr.getAccId().startsWith("XR_") || tr.getAccId().startsWith("NR_")) {
                        return "transcript";
                    }
                    return "mRNA";
                case "pseudo":
                case "pseudogene":
                case "unprocessed_pseudogene":
                    return "pseudogenic_transcript";

                case "lincrna":
                    return "lincRNA";
                case "mirna":
                    return "miRNA";
                case "misc_rna":
                    return "transcript";
                case "ncrna":
                    return "ncRNA";
                case "rrna":
                    return "rRNA";
                case "scrna":
                    return "scRNA";
                case "snorna":
                    return "snoRNA";
                case "snrna":
                    return "snRNA";

                case "gene":
                    if (tr.getAccId().startsWith("XR_")) {
                        return "transcript";
                    }
                    if (tr.getAccId().startsWith("XM_") || tr.getAccId().startsWith("NM_")) {
                        return "mRNA";
                    }
                    return "transcript";
            }
        } else {
            switch(tr.getType()) {
                case "protein_coding":
                case "non_stop_decay":
                case "IG_C_gene":
                case "IG_D_gene":
                case "IG_J_gene":
                case "IG_V_gene":
                case "TR_C_gene":
                case "TR_D_gene":
                case "TR_J_gene":
                case "TR_V_gene":
                    return "mRNA";
                case "processed_transcript":
                    return "processed_transcript";
                case "pseudogene":
                case "polymorphic_pseudogene":
                case "protein_coding_LoF": // replaced 'polymorphic_pseudogene'; per Stan's recommendation, use 'pseudogenic_transcript' mapping
                case "processed_pseudogene":
                case "unprocessed_pseudogene":
                case "translated_processed_pseudogene":
                case "translated_unprocessed_pseudogene":
                case "transcribed_processed_pseudogene":
                case "transcribed_unprocessed_pseudogene":
                case "transcribed_unitary_pseudogene":
                case "unitary_pseudogene":
                case "IG_C_pseudogene":
                case "IG_J_pseudogene":
                case "IG_V_pseudogene":
                case "rRNA_pseudogene":
                case "TR_J_pseudogene":
                case "TR_V_pseudogene":
                    return "pseudogenic_transcript";
                case "TEC":
                    return "unconfirmed_transcript";
                case "lincRNA":
                    return "lincRNA";
                case "sense_intronic":
                case "lncRNA":
                    return "lnc_RNA";
                case "miRNA":
                    return "miRNA";
                case "Mt_rRNA":
                    return "mt_rRNA";
                case "Mt_tRNA":
                    return "mt_tRNA";
                case "misc_RNA":
                case "rRNA":
                    return "rRNA";
                case "ribozyme":
                    return "enzymatic_RNA";
                case "scRNA":
                    return "scRNA";
                case "scaRNA":
                    return "scaRNA";
                case "snRNA":
                    return "snRNA";
                case "snoRNA":
                    return "snoRNA";
                case "vault_RNA":
                case "Y_RNA":
                    return "ncRNA";
                case "retained_intron": // no suitable mappings
                case "antisense":
                case "nonsense_mediated_decay":
                case "protein_coding_CDS_not_defined":
                    return "transcript";
            }
        }
        throw new Exception("tr biotype: unsupported for gene type ["+geneType+"], tr type ["+tr.getType()+"]    "+tr.getAccId());
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

        switch(geneType) {
            case "protein-coding":
            case "protein_coding":
                return "protein_coding_gene";
            case "pseudogene":
            case "pseudo":
                return "pseudogene";
            case "ncrna":
                return "ncRNA_gene";
            case "trna":
                return "tRNA_gene";
            case "rrna":
                return "rRNA_gene";
            case "gene":
            case "misc_rna":
                return "gene";
            case "snrna":
                return "snRNA_gene";
            case "snorna":
                return "snoRNA_gene";
            case "scrna":
                return "scRNA_gene";
            case "lincrna":
                return "lincRNA_gene";
            case "lncrna":
                return "lncRNA_gene";
            case "mirna":
                return "miRNA_gene";
            case "biological-region":
                return "biological_region";
            case "processed_transcript":
                return "processed_transcript";
            case "processed_pseudogene":
                return "processed_pseudogene";
            case "unprocessed_pseudogene":
                return "non_processed_pseudogene";
            case "transcribed_unprocessed_pseudogene":
                return "transcribed_unprocessed_pseudogene";
            case "transcribed_processed_pseudogene":
                return "transcribed_processed_pseudogene";
            case "sense_intronic":
                return "sense_intronic_ncRNA";
            case "antisense":
                return "antisense";
            case "tec":
                return "unconfirmed_transcript";
            case "scarna":
                return "scaRNA";
            case "mt_rrna":
                return "mt_rRNA";
            case "mt_trna":
                return "mt_tRNA";
            case "ribozyme":
                return "ribozyme_gene";
            case "ig_v_gene":
                return "IG_V_gene";
            case "tr_c_gene":
                return "TR_C_Gene";
            case "tr_j_gene":
                return "TR_J_Gene";
            case "tr_v_gene":
                return "TR_V_Gene";
            case "y_rna":
                return "Y_RNA_gene";
            default:
                throw new Exception("unsupported gene type "+geneType);
        }
    }

    String getSoTermAccForGene(String geneType) throws Exception {

        switch(geneType) {
            case "protein-coding":
            case "protein_coding":
                return "SO:0001217";
            case "pseudogene":
            case "pseudo":
                return "SO:0000336";
            case "ncrna":
                return "SO:0001263";
            case "trna":
                return "SO:0001272";
            case "rrna":
                return "SO:0001637";
            case "gene":
            case "misc_rna":
                return "SO:0000704";
            case "snrna":
                return "SO:0001268";
            case "snorna":
                return "SO:0001267";
            case "scrna":
                return "SO:0001266";
            case "lincrna":
                return "SO:0001641";
            case "lncrna":
                return "SO:0002127";
            case "biological-region":
                return "SO:0001411";
            case "processed_transcript":
                return "SO:0001503";
            case "mirna":
                return "SO:0001265";
            case "processed_pseudogene":
                return "SO:0000043";
            case "unprocessed_pseudogene":
                return "SO:0001760";
            case "transcribed_unprocessed_pseudogene":
                return "SO:0002107";
            case "transcribed_processed_pseudogene":
                return "SO:0002109";
            case "sense_intronic":
                return "SO:0002131";
            case "antisense":
                return "SO:0000077";
            case "tec":
                return "SO:0002139";
            case "scarna":
                return "SO:0002095";
            case "mt_rrna":
                return "SO:0002128";
            case "mt_trna":
                return "SO:0002129";
            case "ribozyme":
                return "SO:0002181";
            case "ig_v_gene":
                return "SO:0002126";
            case "tr_c_gene":
                return "SO:0002134";
            case "tr_j_gene":
                return "SO:0002136";
            case "tr_v_gene":
                return "SO:0002137";
            case "y_rna":
                return "SO:0002359";
            default:
                throw new Exception("unsupported gene type "+geneType);
        }
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

        if( nameOfgene!=null ){
            if(nameOfgene.contains(",")){
                nameOfgene = nameOfgene.replaceAll(",", "");
            }
            if(nameOfgene.contains(";")){
                nameOfgene = nameOfgene.replaceAll(";", "");
            }
        } else {
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

    String getUniqueId2(String idBase) {

        Integer cnt = idMap.get(idBase);
        if( cnt==null ) {
            cnt = 1;
        }
        else {
            cnt++;
        }
        idMap.put(idBase, cnt);

        if( cnt==1 ) {
            return idBase;
        }
        return idBase+"_"+cnt;
    }

    // currently we have only UniProtKB accessions as DbXrefs
    Map<Integer, String> loadDbXrefsForTranscripts(int speciesTypeKey) throws Exception {
        XdbId filter = new XdbId();
        filter.setXdbKey(XdbId.XDB_KEY_UNIPROT);
        List<XdbId> ids = dao.getXdbIds(filter, speciesTypeKey, RgdId.OBJECT_KEY_TRANSCRIPTS);
        Map<Integer, String> results = new HashMap<>();
        for( XdbId xdbId: ids ) {
            String uniprotAcc = xdbId.getAccId();
            String str = results.get(xdbId.getRgdId());
            if( str==null ) {
                str = "UniProtKB:"+uniprotAcc;
            } else if( !str.contains(uniprotAcc) ) {
                str += ",UniProtKB:"+uniprotAcc;
            }
            results.put(xdbId.getRgdId(), str);
        }
        return results;
    }

    List<MapData> getMergedMapData(int geneRgdId, int mapKey1, int mapKey2, Counters counters) throws Exception {
        List<MapData> geneMap = getMapData(geneRgdId, mapKey1, counters);
        geneMap.addAll(getMapData(geneRgdId, mapKey2, counters));

        List<MapData> results = new ArrayList<>();
        while( !geneMap.isEmpty() ) {
            MapData md = geneMap.remove(0);

            // merge it with other overlapping gene positions
            boolean wasMerged = false;
            for( int i=0; i<geneMap.size(); i++ ) {
                MapData md2 = geneMap.get(i);
                if( positionsOverlap(md, md2) ) {
                    // merge positions
                    geneMap.remove(i);
                    if( md2.getStartPos()<md.getStartPos() ) {
                        md.setStartPos(md2.getStartPos());
                    }
                    if( md2.getStopPos()>md.getStopPos() ) {
                        md.setStopPos(md2.getStopPos());
                    }
                    geneMap.add(md);
                    wasMerged = true;
                    break;
                }
            }

            if( !wasMerged ) {
                results.add(md);
            }
        }
        return results;
    }

    boolean positionsOverlap(MapData md1, MapData md2) {
        // chromosomes must match
        if( !Utils.stringsAreEqualIgnoreCase(md1.getChromosome(), md2.getChromosome()) )
            return false;
        // positions must overlap
        if( md1.getStopPos() < md2.getStartPos() )
            return false;
        if( md2.getStopPos() < md1.getStartPos() )
            return false;
        return true;
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

    public int getMapKeyEnsembl() {
        return mapKeyEnsembl;
    }

    public void setMapKeyEnsembl(int mapKeyEnsembl) {
        this.mapKeyEnsembl = mapKeyEnsembl;
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
