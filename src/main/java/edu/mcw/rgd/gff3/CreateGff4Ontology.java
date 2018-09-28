package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.datamodel.ontologyx.Ontology;
import edu.mcw.rgd.datamodel.ontologyx.Term;
import edu.mcw.rgd.datamodel.ontologyx.TermDagEdge;
import edu.mcw.rgd.datamodel.ontologyx.TermWithStats;
import edu.mcw.rgd.process.Utils;

import java.io.File;
import java.io.PrintWriter;
import java.net.URLEncoder;
import java.util.*;
import java.util.Map;

public class CreateGff4Ontology {
    AssociationDAO assDao = new AssociationDAO();
    RgdGff3Dao dao = new RgdGff3Dao();

    OntologyXDAO ontXdao = new OntologyXDAO();
    QTLDAO qtldao = new QTLDAO();

    List<String> chromosomes;
    int speciesTypeKey=0;
    int mapKey=0;
    String toFile;
    String ontAspect;


    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }

    public int getMapKey() {
        return mapKey;
    }

    public void setMapKey(int mapKey) {
        this.mapKey = mapKey;
    }

    public String getToFile() {
        return toFile;
    }

    public void setToFile(String toFile) {
        this.toFile = toFile;
    }

    public String getOntAspect() {
        return ontAspect;
    }

    public void setOntAspect(String ontAspect) {
        this.ontAspect = ontAspect;
    }

    public List<String> getChromosomes() {
        return chromosomes;
    }

    public void setChromosomes(List<String> chromosomes) {
        this.chromosomes = chromosomes;
    }

    public void run(boolean compress) throws Exception{

        switch( getOntAspect() ) {
            case "D":
                Ontology ont = ontXdao.getOntologyFromAspect(getOntAspect());
                String rootTermAcc = ontXdao.getRootTerm(ont.getId());
                run(compress, rootTermAcc);
                break;

            case "E":
                final String termAcc = "CHEBI:24432"; // CHEBI term 'biological_role'
                run(compress, termAcc);
                break;

            default:
                throw new Exception("Unsupported aspect "+getOntAspect());
        }
    }

    public void run(boolean compress, String rootTermAcc) throws Exception{

        long t0 = System.currentTimeMillis();

        String fname = getToFile();
        new File(fname).mkdirs(); // make sure the directory is always created

        String gffFile;
        Gff3ColumnWriter gff3Writer;

        List<RGDInfo> rgdInfoList = new ArrayList<>();

        Ontology ont = ontXdao.getOntologyFromAspect(getOntAspect());

        FileGuard fileGuard = new FileGuard();
        fileGuard.init(fname, ont.getId());

        for( String termAcc: dao.getTermDescendants(rootTermAcc).keySet() ) {
            if( !isTermAnnotated(termAcc) )
                continue;

            Map<String, Term> mapAccAnnList = getAnnotatedChildTerms(termAcc);

            Term rootTerm = dao.getTerm(termAcc);
            System.out.println("Number of annotations for " + rootTerm.getAccId() + " " +rootTerm.getTerm());

            for (String chr : getChromosomes()) {

                int counter = 0;

                gffFile = fname + termAcc.replaceAll(":", "") + "_Ontology_" + SpeciesType.getCommonName(speciesTypeKey) + "_RGDChr" + chr + ".gff3";
                gff3Writer = new Gff3ColumnWriter(gffFile, false, compress);

                for (MapData md : dao.getMapDataByMapKeyChr(chr, mapKey, RgdId.OBJECT_KEY_GENES)) {
                    if( md.getStartPos() != null && md.getStopPos() != null ) {

                        Gff3Entry entry = new Gff3Entry(RgdId.OBJECT_KEY_GENES, md);
                        processGene(entry);

                        counter += processAnnotations(mapAccAnnList, entry);

                        if (!entry.anns.equals("NA")) {
                            createRGdInfo(entry, rgdInfoList);
                            writeGff3Line(gff3Writer, entry, termAcc);
                        }
                    }
                }

                for (MapData md : dao.getMapDataByMapKeyChr(chr, mapKey, RgdId.OBJECT_KEY_QTLS)) {
                    if( md.getStartPos() != null && md.getStopPos() != null ) {

                        Gff3Entry entry = new Gff3Entry(RgdId.OBJECT_KEY_QTLS, md);
                        processQtl(entry);

                        counter += processAnnotations(mapAccAnnList, entry);

                        if (!entry.anns.equals("NA")) {
                            createRGdInfo(entry, rgdInfoList);
                            writeGff3Line(gff3Writer, entry, termAcc);
                        }
                    }
                }

                if( speciesTypeKey==SpeciesType.RAT )
                for (MapData md : dao.getMapDataByMapKeyChr(chr, mapKey, RgdId.OBJECT_KEY_STRAINS)) {
                    if( md.getStartPos() != null && md.getStopPos() != null ) {

                        Gff3Entry entry = new Gff3Entry(RgdId.OBJECT_KEY_STRAINS, md);
                        processStrain(entry);

                        counter += processAnnotations(mapAccAnnList, entry);

                        if (!entry.anns.equals("NA")) {
                            createRGdInfo(entry, rgdInfoList);
                            writeGff3Line(gff3Writer, entry, termAcc);
                        }
                    }
                }

                gff3Writer.close();

                System.out.print("c"+chr + ":" + counter+", ");
            }
            System.out.println();
        }

		calculateDensity(rgdInfoList);

        fileGuard.check(mapKey);

        long t1 = System.currentTimeMillis();
        System.out.println("============");
        System.out.println("\nTotal time elapsed " + Utils.formatElapsedTime(t0, t1));
    }

    Map<String, Term> getAnnotatedChildTerms(String termAcc) throws Exception {
        List<TermDagEdge> termWithStatsList = ontXdao.getAllChildEdges(termAcc);

        Map<String, Term> mapAccAnnList = new HashMap<>();
        for (TermDagEdge tws : termWithStatsList) {
            TermWithStats t = dao.getTerm(tws.getChildTermAcc());
            if( isTermAnnotated(tws.getChildTermAcc()) ) {
                mapAccAnnList.put(t.getAccId(), t);
            }
        }
        return mapAccAnnList;
    }

    void processGene(Gff3Entry entry) throws Exception {

        Gene geneObj = dao.getGeneWithDescription(entry.rgdId);
        entry.name = geneObj.getName();
        entry.symbol = geneObj.getSymbol();


        if (Utils.isStringEmpty(geneObj.getDescription())) {
            entry.note = "null";
        } else {
            entry.note = geneObj.getDescription().replace(";", " AND ");
        }


        if (geneObj.getType() == null) {
            entry.gType = "null-null";
        } else {
            entry.gType = geneObj.getType();
        }
    }

    void processQtl(Gff3Entry entry) throws Exception {

        QTL qtlObj = qtldao.getQTL(entry.rgdId);
        entry.name = qtlObj.getName();
        entry.symbol = qtlObj.getSymbol();

        if (qtlObj.getLod() != null) {
            entry.lod = String.valueOf(qtlObj.getLod());
        } else {
            entry.lod = "null";
        }
        if (qtlObj.getPValue() != null) {
            entry.pValue = String.valueOf(qtlObj.getPValue());
        } else {
            entry.pValue = "null";
        }
        if (qtlObj.getNotes() != null) {
            entry.note = qtlObj.getNotes().replaceAll(";", " ");
            entry.note = URLEncoder.encode(entry.note, "UTF-8");
        } else {
            entry.note = "null";
        }


        //get related strains
        if (speciesTypeKey == SpeciesType.RAT) {
            List<Strain> relatedStrainsList = assDao.getStrainAssociationsForQTL(entry.rgdId);

            if (relatedStrainsList.size() > 0) {

                for (Strain st : relatedStrainsList) {
                    entry.relStrain += st.getSymbol() + ":" + st.getRgdId() + ",";
                }
            } else {

                entry.relStrain = "NA";
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
        } else {
            entry.relGenes = "NA";
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
            } else if (relatedQtlMap.size() == 0) {


                entry.relQtls = "NA";
            }
            if (entry.relQtls.endsWith(",")) {
                entry.relQtls = entry.relQtls.substring(0, entry.relQtls.length() - 1);
            }
        } else {
            entry.relQtls = "NA";
        }
    }

    void processStrain(Gff3Entry entry) throws Exception {

        Strain strObject = dao.getStrain(entry.rgdId);
        entry.name = strObject.getStrain();
        entry.symbol = strObject.getSymbol();

        entry.note = strObject.getOrigin();
        if (entry.note != null) {
            entry.note = URLEncoder.encode(entry.note, "UTF-8");
        } else {
            entry.note = "NA";
        }

        entry.parentStr = strObject.getStrain();
        if (entry.parentStr != null) {
            entry.parentStr = URLEncoder.encode(entry.parentStr, "UTF-8");
        } else {
            entry.parentStr = "NA";
        }

        entry.subStr = strObject.getSubstrain();
        if (entry.subStr != null) {
            entry.subStr = URLEncoder.encode(entry.subStr, "UTF-8");
        } else {
            entry.subStr = "NA";
        }

        entry.src = strObject.getSource();
        if (entry.src != null) {
            entry.src = URLEncoder.encode(entry.src, "UTF-8");
        } else {
            entry.src = "NA";
        }
    }

    int processAnnotations(Map<String, Term> mapAccAnnList, Gff3Entry entry) throws Exception {

        Map<String,String> annotMap = new TreeMap<>();
        for( Annotation ann: dao.getAnnotationsByAspect(entry.rgdId, getOntAspect()) ) {

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

        if( annotMap.isEmpty() ) {
            entry.anns = "NA";
        } else {
            entry.anns = Utils.concatenate(annotMap.keySet(), ",");
            entry.aliasTerm = Utils.concatenate(annotMap.values(), ",");
        }

        return annotMap.size();
    }

    void createRGdInfo(Gff3Entry entry, List<RGDInfo> rgdInfoList) {
        RGDInfo rgdInfo = new RGDInfo();
        rgdInfo.setChromosome(entry.chrom);
        rgdInfo.setRgdId(entry.rgdId);
        rgdInfo.setStartPos(entry.start);
        rgdInfo.setStopPos(entry.stop);
        rgdInfoList.add(rgdInfo);
    }

    void writeGff3Line(Gff3ColumnWriter gff3Writer, Gff3Entry entry, String termAcc) throws Exception {

        gff3Writer.writeFirst8Columns(entry.chrom,
                "RGD_Ontology_" + termAcc.replaceAll(":", ""),
                RgdId.getObjectTypeName(entry.objectKey).toLowerCase() + "_ont",
                entry.start,
                entry.stop,
                ".",
                entry.strand,
                ".");

        Map<String, String> attributesHashMap = new HashMap<>();
        attributesHashMap.put("Dbxref", "RGD:" + entry.rgdId);
        attributesHashMap.put("ID", "RGD" + entry.rgdId);
        attributesHashMap.put("Name", entry.name.replaceAll(",", "").replaceAll(";","|"));
        attributesHashMap.put("symbol", entry.symbol);
        attributesHashMap.put("note", entry.note.replaceAll(",", ""));
        attributesHashMap.put("Alias", entry.aliasTerm);
        if (entry.objectKey == 1) {
            attributesHashMap.put("geneType", entry.gType);
            attributesHashMap.put("objectTypeName", "Gene");
        } else if( speciesTypeKey == 3 && entry.objectKey == 5 ) {
            attributesHashMap.put("parentStr", entry.parentStr);
            attributesHashMap.put("subStr", entry.subStr);
            attributesHashMap.put("src", entry.src);
            attributesHashMap.put("objectTypeName", "Strain");
        } else if( entry.objectKey == 6 ) {
            attributesHashMap.put("LOD", entry.lod);
            attributesHashMap.put("pValue", entry.pValue);
            attributesHashMap.put("relatedQTLs", entry.relQtls);
            if (speciesTypeKey == 3) {
                attributesHashMap.put("relatedStrains", entry.relStrain);
            }
            attributesHashMap.put("relatedGenes", entry.relGenes);
            attributesHashMap.put("objectTypeName", "QTL");
        }

        attributesHashMap.put("Ontology_term", entry.anns);

        gff3Writer.writeAttributes4Gff3(attributesHashMap);
    }

    void calculateDensity(List<RGDInfo> rgdInfoList) throws Exception {

        HashMap<Integer, RGDInfo> rgdinfoMap = new HashMap<>();
        for(RGDInfo rgdinf:rgdInfoList){
            rgdinfoMap.put(rgdinf.getRgdId(),rgdinf);
        }

        List<RGDInfo> listOfRgdInfos= new ArrayList<>(rgdinfoMap.values());
        CalculateDensity calcDensity = new CalculateDensity();
        calcDensity.setRgdInfoList(listOfRgdInfos);
        PrintWriter densityWriter = new PrintWriter(getToFile() + "density.gff");
        calcDensity.setDensityDiseaseWriter(densityWriter);
        calcDensity.runDensityCalculator();
    }

    boolean isTermAnnotated(String termAcc) throws Exception {
        TermWithStats t = dao.getTerm(termAcc);
        return t.getAnnotObjectCountForSpecies(getSpeciesTypeKey(), true)>0;
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
}

