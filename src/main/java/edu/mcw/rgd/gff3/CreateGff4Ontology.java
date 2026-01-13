package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.datamodel.ontologyx.*;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.net.URLEncoder;
import java.util.*;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class CreateGff4Ontology {
    AssociationDAO assDao = new AssociationDAO();
    RgdGff3Dao dao = new RgdGff3Dao();

    OntologyXDAO ontXdao = new OntologyXDAO();
    QTLDAO qtldao = new QTLDAO();

    /// new properties
    private String outDirForDiseases;
    private String outDirForChebi;
    private List<Integer> processedMapKeys;
    private Map<String, String> termTrackNames;


    public void runDiseaseOntology() throws Exception {

        Logger log = LogManager.getLogger("disease");

        Collection<String> doTermAccs = getTermsInJBrowseSlim();

        AtomicInteger mapKeysDone = new AtomicInteger(0);
        getProcessedMapKeys().stream().parallel().forEach( processedMapKey -> {
            try {
                runOntology(doTermAccs, getOutDirForDiseases(), "D", log, processedMapKey, mapKeysDone);
            } catch( Exception e ) {
                Utils.printStackTrace(e, log);
                throw new RuntimeException(e);
            }
        });
    }

    public void runChebiOntology() throws Exception {

        Logger log = LogManager.getLogger("chebi");

        final String termAcc = "CHEBI:24432"; // CHEBI term 'biological_role'
        Collection<String> termAccs = dao.getTermDescendants(termAcc).keySet();

        AtomicInteger mapKeysDone = new AtomicInteger(0);
        for( Integer processedMapKey: getProcessedMapKeys() ) {
            runOntology(termAccs, getOutDirForChebi(), "E", log, processedMapKey, mapKeysDone);
        }
    }

    public void runOntology(Collection<String> doTermAccs, String outDirName, String ontAspect, Logger log, int mapKey, AtomicInteger mapKeysDone) throws Exception {
        long t0 = System.currentTimeMillis();

        int compressMode = Gff3ColumnWriter.COMPRESS_MODE_BGZIP;
        int speciesTypeKey = MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey();

        String assemblyDir = Manager.getInstance().getAssemblies().get(mapKey);
        if( assemblyDir==null ) {
            return;
        }
        assemblyDir += "/" + Gff3Utils.getAssemblyDirStandardized(mapKey);

        String mainDir = assemblyDir + "/" + outDirName;

        log.info(mainDir+" ...");

        SequenceRegionWatcher sequenceRegionWatcher = new SequenceRegionWatcher(0, null, null);

        ArrayList<String> doTermAccessions = new ArrayList<>(doTermAccs);
        Collections.shuffle(doTermAccessions);
        for( String termAcc: doTermAccessions ) {

            log.debug("  "+mainDir+"  "+termAcc+" ...");

            String trackName = getTermTrackNames().get(termAcc);
            if (trackName == null) {
                log.error("**** ERROR *** TERM " + termAcc + " DOES NOT HAVE A MAPPING IN CONFIG FILE -- SKIPPED!");
                continue;
            }

            if (!isTermAnnotated(termAcc, speciesTypeKey))
                continue;

            String outDir = mainDir + "/" + trackName;
            new File(outDir).mkdirs(); // make sure the directory is always created

            // write as comment: date and time it was generated, species, assembly and ontology term
            String gffHeader = "";
            gffHeader += ("#generated on " + new Date() + "\n");
            gffHeader += ("#" + SpeciesType.getCommonName(speciesTypeKey) + " assembly " + MapManager.getInstance().getMap(mapKey).getName() + "\n");
            gffHeader += ("#disease ontology track for term '" + dao.getTerm(termAcc).getTerm() + "' (" + termAcc + ")\n");

            Map<String, Term> mapAccAnnList = getAnnotatedChildTerms(termAcc, speciesTypeKey);

            String gffFile = outDir + "/" + trackName + " Related Genes.gff3";
            Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gffFile, compressMode);
            gff3Writer.print(gffHeader);
            sequenceRegionWatcher.init(mapKey, gff3Writer, dao);

            List<RGDInfo> rgdInfoList = new ArrayList<>();
            int counter = processGenes("*", termAcc, gff3Writer, mapAccAnnList, rgdInfoList, sequenceRegionWatcher, ontAspect, mapKey, speciesTypeKey);

            gff3Writer.close();
            if (counter > 0) {
                gff3Writer.sortInMemory();
            } else {
                // no data -- delete it
                new File(gff3Writer.getOutFileName()).delete();
            }

            //Term rootTerm = dao.getTerm(termAcc);
            String summaryMsg = "  "+mainDir+"  "+termAcc + " written genes: "+counter;


            //// QTLS

            counter = 0;

            gffFile = outDir + "/" + trackName + " Related QTLs.gff3";
            gff3Writer = new Gff3ColumnWriter(gffFile, compressMode);
            gff3Writer.print(gffHeader);
            sequenceRegionWatcher.init(mapKey, gff3Writer, dao);

            for (MapData md : dao.getMapDataByMapKeyChr("*", mapKey, RgdId.OBJECT_KEY_QTLS)) {
                if (md.getStartPos() != null && md.getStopPos() != null) {

                    Gff3Entry entry = new Gff3Entry(RgdId.OBJECT_KEY_QTLS, md);
                    processQtl(entry, speciesTypeKey);

                    counter += processAnnotations(mapAccAnnList, entry, ontAspect, speciesTypeKey);

                    if( !Utils.isStringEmpty(entry.anns) ) {
                        createRGdInfo(entry, rgdInfoList);
                        writeGff3Line(gff3Writer, entry, termAcc, speciesTypeKey);
                        sequenceRegionWatcher.emit(md.getChromosome());
                    }
                }
            }

            gff3Writer.close();
            if (counter > 0) {
                gff3Writer.sortInMemory();
            } else {
                // no data -- delete it
                new File(gff3Writer.getOutFileName()).delete();
            }
            summaryMsg += ", qtls: "+counter;


            //// RAT STRAINS

            if (speciesTypeKey == SpeciesType.RAT) {

                counter = 0;

                gffFile = outDir + "/" + trackName + " Related Strains.gff3";
                gff3Writer = new Gff3ColumnWriter(gffFile, compressMode);
                gff3Writer.print(gffHeader);
                sequenceRegionWatcher.init(mapKey, gff3Writer, dao);

                for (MapData md : dao.getMapDataByMapKeyChr("*", mapKey, RgdId.OBJECT_KEY_STRAINS)) {
                    if (md.getStartPos() != null && md.getStopPos() != null) {

                        Gff3Entry entry = new Gff3Entry(RgdId.OBJECT_KEY_STRAINS, md);
                        processStrain(entry);

                        counter += processAnnotations(mapAccAnnList, entry, ontAspect, speciesTypeKey);

                        if( !Utils.isStringEmpty(entry.anns) ) {
                            createRGdInfo(entry, rgdInfoList);
                            writeGff3Line(gff3Writer, entry, termAcc, speciesTypeKey);
                            sequenceRegionWatcher.emit(md.getChromosome());
                        }
                    }
                }

                gff3Writer.close();
                if (counter > 0) {
                    gff3Writer.sortInMemory();
                } else {
                    // no data -- delete it
                    new File(gff3Writer.getOutFileName()).delete();
                }
                summaryMsg += ", strains: "+counter;
            }
            log.info(summaryMsg);
        }

        long t1 = System.currentTimeMillis();
        int jobsDone = mapKeysDone.incrementAndGet();
        int jobCount = getProcessedMapKeys().size();
        log.info("===");
        log.info("=== "+jobsDone+"/"+jobCount+".  "+mainDir+" DONE! Time elapsed " + Utils.formatElapsedTime(t0, t1));
        log.info("===");
    }

    Map<String, Term> getAnnotatedChildTerms(String termAcc, int speciesTypeKey) throws Exception {
        List<TermDagEdge> termWithStatsList = ontXdao.getAllChildEdges(termAcc);

        Map<String, Term> mapAccAnnList = new HashMap<>();
        for (TermDagEdge tws : termWithStatsList) {
            TermWithStats t = dao.getTerm(tws.getChildTermAcc());
            if( isTermAnnotated(tws.getChildTermAcc(), speciesTypeKey) ) {
                mapAccAnnList.put(t.getAccId(), t);
            }
        }
        return mapAccAnnList;
    }

    void processGene(Gff3Entry entry) throws Exception {

        Gene geneObj = dao.getGeneWithDescription(entry.rgdId);
        entry.name = Utils.defaultString(geneObj.getName()); // ensure gene name is not null
        entry.symbol = geneObj.getSymbol();


        if( !Utils.isStringEmpty(geneObj.getDescription()) ) {
            entry.note = geneObj.getDescription().replace(";", " AND ");
        }

        if( geneObj.getType() != null ) {
            entry.gType = geneObj.getType();
        }
    }

    void processQtl(Gff3Entry entry, int speciesTypeKey) throws Exception {

        QTL qtlObj = qtldao.getQTL(entry.rgdId);
        entry.name = qtlObj.getName();
        entry.symbol = qtlObj.getSymbol();

        if( qtlObj.getLod() != null ) {
            entry.lod = qtlObj.getLod().toString();
        }
        if (qtlObj.getPValue() != null) {
            entry.pValue = qtlObj.getPValue().toString();
        }
        if (qtlObj.getNotes() != null) {
            entry.note = qtlObj.getNotes().replaceAll(";", " ");
            entry.note = URLEncoder.encode(entry.note, "UTF-8");
        }


        //get related strains
        if (speciesTypeKey == SpeciesType.RAT) {
            List<Strain> relatedStrainsList = assDao.getStrainAssociationsForQTL(entry.rgdId);

            if (relatedStrainsList.size() > 0) {

                for (Strain st : relatedStrainsList) {
                    entry.relStrain += st.getSymbol() + ":" + st.getRgdId() + ",";
                }
            }
            if (entry.relStrain.endsWith(",")) {
                entry.relStrain = entry.relStrain.substring(0, entry.relStrain.length() - 1);
            }
        }

        //get related genes.
        List<Gene> relatedGenesList = assDao.getGeneAssociationsByQTL(entry.rgdId);
        if (relatedGenesList.size() > 0) {
            for (Gene g : relatedGenesList) {
                entry.relGenes += g.getSymbol() + ":" + g.getRgdId() + ",";
            }
        }
        if (entry.relGenes.endsWith(",")) {
            entry.relGenes = entry.relGenes.substring(0, entry.relGenes.length() - 1);
        }

        //get qtl - qtl associations.
        Map<Integer, String> relatedQtlMap = assDao.getQtlToQtlAssociations(qtlObj.getKey());
        if (relatedQtlMap != null) {
            if (relatedQtlMap.size() > 0) {
                Set<Integer> relQtlKeys = relatedQtlMap.keySet();

                for (int keys : relQtlKeys) {
                    entry.relQtls += keys + ":" + relatedQtlMap.get(keys) + ",";
                }
            }
            if (entry.relQtls.endsWith(",")) {
                entry.relQtls = entry.relQtls.substring(0, entry.relQtls.length() - 1);
            }
        }
    }

    void processStrain(Gff3Entry entry) throws Exception {

        Strain strObject = dao.getStrain(entry.rgdId);
        entry.name = strObject.getStrain();
        entry.symbol = strObject.getSymbol();

        entry.note = strObject.getOrigination();
        if (entry.note != null) {
            entry.note = URLEncoder.encode(entry.note, "UTF-8");
        }

        entry.parentStr = strObject.getStrain();
        if (entry.parentStr != null) {
            entry.parentStr = URLEncoder.encode(entry.parentStr, "UTF-8");
        }

        entry.subStr = strObject.getSubstrain();
        if (entry.subStr != null) {
            entry.subStr = URLEncoder.encode(entry.subStr, "UTF-8");
        }

        entry.src = strObject.getSource();
        if (entry.src != null) {
            entry.src = URLEncoder.encode(entry.src, "UTF-8");
        }
    }

    int processAnnotations(Map<String, Term> mapAccAnnList, Gff3Entry entry, String ontAspect, int speciesTypeKey) throws Exception {

        Map<String,String> annotMap = new TreeMap<>();
        for( Annotation ann: dao.getAnnotationsByAspect(entry.rgdId, ontAspect, speciesTypeKey) ) {

            if (mapAccAnnList.containsKey(ann.getTermAcc())) {

                String qualifier;
                if( Utils.isStringEmpty(ann.getQualifier()) ) {
                    qualifier = "";
                } else {
                    qualifier = "||"+ann.getQualifier();
                    if (qualifier.contains(";")) {
                        qualifier = qualifier.replaceAll(";", "");
                    }
                    if (qualifier.contains(",")) {
                        qualifier = qualifier.replaceAll(",", "");
                    }
                }

                String fullterm = ann.getTerm();
                if (fullterm.contains(";")) {
                    fullterm = fullterm.replaceAll(";", "");
                }
                if (fullterm.contains(",")) {
                    fullterm = fullterm.replaceAll(",", "");
                }

                String key = ann.getTermAcc() + "_" + fullterm + qualifier;
                String value = "Disease:" + fullterm;
                annotMap.put(key, value);
            }
        }

        if( !annotMap.isEmpty() ) {
            entry.anns = Utils.concatenate(annotMap.keySet(), ",");
            entry.aliasTerm = Utils.concatenate(annotMap.values(), ",");
        }

        return annotMap.size();
    }

    synchronized void createRGdInfo(Gff3Entry entry, List<RGDInfo> rgdInfoList) {
        RGDInfo rgdInfo = new RGDInfo();
        rgdInfo.setChromosome(entry.chrom);
        rgdInfo.setRgdId(entry.rgdId);
        rgdInfo.setStartPos(entry.start);
        rgdInfo.setStopPos(entry.stop);
        rgdInfoList.add(rgdInfo);
    }

    void writeGff3Line(Gff3ColumnWriter gff3Writer, Gff3Entry entry, String termAcc, int speciesTypeKey) {

        gff3Writer.print( prepGff3Line(gff3Writer, entry, termAcc, speciesTypeKey) );
    }

    String prepGff3Line(Gff3ColumnWriter gff3Writer, Gff3Entry entry, String termAcc, int speciesTypeKey) {

        StringBuffer buf = new StringBuffer(
            gff3Writer.prepFirst8Columns(entry.chrom,
                "RGD_Ontology_" + termAcc.replaceAll(":", ""),
                RgdId.getObjectTypeName(entry.objectKey).toLowerCase() + "_ont",
                entry.start,
                entry.stop,
                ".",
                entry.strand,
                ".")
        );

        Map<String, String> attributesHashMap = new HashMap<>();
        attributesHashMap.put("Dbxref", "RGD:" + entry.rgdId);
        attributesHashMap.put("ID", "RGD" + entry.rgdId);
        if( !Utils.isStringEmpty(entry.name) ) {
            attributesHashMap.put("geneName", entry.name.replaceAll(",", "").replaceAll(";", "|"));
        }
        attributesHashMap.put("Name", entry.symbol);
        if( !Utils.isStringEmpty(entry.note) ) {
            attributesHashMap.put("info", entry.note.replaceAll(",", ""));
        }
        attributesHashMap.put("Alias", entry.aliasTerm);
        if (entry.objectKey == 1) {
            if( !Utils.isStringEmpty(entry.gType) ) {
                attributesHashMap.put("geneType", entry.gType);
            }
            attributesHashMap.put("objectTypeName", "Gene");
        } else if( speciesTypeKey == 3 && entry.objectKey == 5 ) {
            if( !Utils.isStringEmpty(entry.parentStr) ) {
                attributesHashMap.put("parentStrain", entry.parentStr);
            }
            if( !Utils.isStringEmpty(entry.subStr) ) {
                attributesHashMap.put("subStrain", entry.subStr);
            }
            if( !Utils.isStringEmpty(entry.src) ) {
                attributesHashMap.put("strainSource", entry.src);
            }
            attributesHashMap.put("objectTypeName", "Strain");
        } else if( entry.objectKey == 6 ) {
            if( !Utils.isStringEmpty(entry.lod) ) {
                attributesHashMap.put("lod", entry.lod);
            }
            if( !Utils.isStringEmpty(entry.pValue) ) {
                attributesHashMap.put("pValue", entry.pValue);
            }
            if( !Utils.isStringEmpty(entry.relQtls) ) {
                attributesHashMap.put("relatedQTLs", entry.relQtls);
            }
            if( !Utils.isStringEmpty(entry.relStrain) ) {
                attributesHashMap.put("relatedStrains", entry.relStrain);
            }
            if( !Utils.isStringEmpty(entry.relGenes) ) {
                attributesHashMap.put("relatedGenes", entry.relGenes);
            }
            attributesHashMap.put("objectTypeName", "QTL");
        }

        if( !Utils.isStringEmpty(entry.anns) ) {
            attributesHashMap.put("Ontology_term", entry.anns);
        }

        buf.append(gff3Writer.prepAttributes4Gff3(attributesHashMap));
        return buf.toString();
    }

    ///////
    //// highly parallel code

    int processGenes(String chr, String termAcc, Gff3ColumnWriter gff3Writer, Map<String, Term> mapAccAnnList,
                     List<RGDInfo> rgdInfoList, SequenceRegionWatcher sequenceRegionWatcher, String ontAspect,
                     int mapKey, int speciesTypeKey) throws Exception {

        AtomicInteger annotCount = new AtomicInteger(0);
        StringBuffer gff3Lines = new StringBuffer();

        dao.getMapDataByMapKeyChr(chr, mapKey, RgdId.OBJECT_KEY_GENES).parallelStream().forEach(md -> {

            if (md.getStartPos() != null && md.getStopPos() != null) {

                sequenceRegionWatcher.emit(md.getChromosome());

                try {
                    Gff3Entry entry = new Gff3Entry(RgdId.OBJECT_KEY_GENES, md);
                    processGene(entry);

                    annotCount.addAndGet(processAnnotations(mapAccAnnList, entry, ontAspect, speciesTypeKey));

                    if( !Utils.isStringEmpty(entry.anns) ) {
                        createRGdInfo(entry, rgdInfoList);
                        gff3Lines.append(prepGff3Line(gff3Writer, entry, termAcc, speciesTypeKey));
                    }
                } catch(Exception e) {
                    throw new RuntimeException(e);
                }
            }
        });

        gff3Writer.print(gff3Lines.toString());
        return annotCount.get();
    }

    boolean isTermAnnotated(String termAcc, int speciesTypeKey) throws Exception {
        TermWithStats t = dao.getTerm(termAcc);
        return t.getAnnotObjectCountForSpecies(speciesTypeKey, true)>0;
    }

    Collection<String> getTermsInJBrowseSlim() throws Exception {
        Set<String> termAccs = new HashSet<>();

        List<TermSynonym> synonyms = ontXdao.getActiveSynonymsByName("RDO", "RGD_JBrowse_slim");
        for( TermSynonym synonym: synonyms ) {
            termAccs.add(synonym.getTermAcc());
        }
        return termAccs;
    }

    class Gff3Entry {
        int rgdId;
        int objectKey;
        String chrom;
        Integer start;
        Integer stop;
        String name = "";
        String symbol = "";
        String anns = "";
        String aliasTerm = "";
        String note = "";
        String gType = "";
        String lod = "";
        String pValue = "";
        String relStrain = "";
        String relGenes = "";
        String relQtls = "";
        String parentStr = "";
        String subStr = "";
        String src = "";
        String strand;

        public Gff3Entry(int objectKey, MapData md) {
            this.rgdId = md.getRgdId();
            this.objectKey = objectKey;
            this.chrom = md.getChromosome();
            this.start = md.getStartPos();
            this.stop = md.getStopPos();
            setStrand(md.getStrand());
        }

        public void setStrand(String strand) {
            this.strand = strand==null ? "." : strand;
        }
    }

    public String getOutDirForDiseases() {
        return outDirForDiseases;
    }

    public void setOutDirForDiseases(String outDirForDiseases) {
        this.outDirForDiseases = outDirForDiseases;
    }

    public String getOutDirForChebi() {
        return outDirForChebi;
    }

    public void setOutDirForChebi(String outDirForChebi) {
        this.outDirForChebi = outDirForChebi;
    }

    public List<Integer> getProcessedMapKeys() {
        return processedMapKeys;
    }

    public void setProcessedMapKeys(List<Integer> processedMapKeys) {
        this.processedMapKeys = processedMapKeys;
    }

    public Map<String, String> getTermTrackNames() {
        return termTrackNames;
    }

    public void setTermTrackNames(Map<String, String> termTrackNames) {
        this.termTrackNames = termTrackNames;
    }
}

