package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.util.ArrayList;
import java.util.List;

/**
* @author BBakir
* @since Jun 23, 2008
* @author mtutaj
*/
public class Manager {

    private int speciesTypekey;
    private int mapKey;
    private int mapKey2;
    private int objectTypeKey;
    private List<String> chromosomes;
    private String toFile;
    private String toDir;
    int sampleID;
    int patientID;
    String build;
    String source="RGD_CNRat";
    int compressMode = Gff3ColumnWriter.COMPRESS_MODE_NONE;
    String flavor;
    private String version;

    private java.util.Map<Integer, String> assemblies;

    Logger log = LogManager.getLogger("status");

    public static final int OBJECT_KEY_DB_SNP = -1;


    static private Manager _instance;

    static public Manager getInstance() {
        return _instance;
    }

    static public void main(String[] args) throws Exception {
        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        Manager creator = (Manager) (bf.getBean("manager"));
        _instance = creator;

        System.out.println(creator.getVersion());

        try{
            creator.doMain(args, bf);
        }catch (Exception e){
            Utils.printStackTrace(e, creator.log);
        }
    }


    public void doMain(String[] args, DefaultListableBeanFactory bf) throws Exception {
        if( parseArgs(args, bf) ) {
            return;
        }

        if(objectTypeKey!=0){
            handleObjects();
        }
        else if(sampleID!=0 || patientID!=0){
            if( source==null || toFile==null  ){
                throw new ArgumentsException("This script requires '-source:' '-toFile:' as parameters");
            }

            CreateGff4CarpeNovo createGffCarpe = new CreateGff4CarpeNovo();
            createGffCarpe.setFileSource(source);
            createGffCarpe.setToFile(toFile);
            //createGffCarpe.setChromosomes(getChromosomes());
            if( sampleID!=0 )
                createGffCarpe.createGff3ForSample(sampleID);
            if( patientID!=0 ) {
                createGffCarpe.createGff3ForPatient(patientID);
            }

        }
        else if( flavor.equals("ensembl_prep") ) {
            EnsemblPrep ep = (EnsemblPrep) bf.getBean("ensemblPrep");
            ep.run();

        } else {
            throw new ArgumentsException("This script requires '-object:' or '-sampleID' as a parameter:\n" +
                    "-sampleID:309/329/330\t-object:gene/qtl/sslp/strain");
        }
    }

    void handleObjects() throws Exception {
        switch(objectTypeKey) {
            case RgdId.OBJECT_KEY_GENES:
                if( mapKey>0 && toDir!=null && speciesTypekey!=0 ){
                    if( flavor==null ) {
                        CreateGff4Gene createGff = new CreateGff4Gene();
                        CreateInfo info = new CreateInfo();
                        info.setMapKey(mapKey);
                        info.setToDir(toDir);
                        info.setSpeciesTypeKey(speciesTypekey);
                        info.setCompressMode(compressMode);
                        createGff.createGeneGff3(info);
                    } else if( flavor.equals("AGR") ){
                        CreateGff4GeneAgr createGff = new CreateGff4GeneAgr();
                        createGff.setMapKey(mapKey);
                        createGff.setMapKeyEnsembl(mapKey2);
                        createGff.setGff3Path(toDir);
                        createGff.setSpeciesTypeKey(speciesTypekey);
                        createGff.createGeneGff3(compressMode);
                    } else {
                        throw new ArgumentsException("unknown gene flavor: "+flavor);
                    }
                }else{
                    throw new ArgumentsException("This Script requires '-mapKey:','-toDir:','-species:' as parameters");
                }
                break;

            case RgdId.OBJECT_KEY_QTLS:
                if(speciesTypekey!=0){
                    CreateGff4QTL createGff = new CreateGff4QTL();
                    if( mapKey>0 && toDir!=null ){
                        CreateInfo info = new CreateInfo();
                        info.setMapKey(mapKey);
                        info.setToDir(toDir);
                        info.setSpeciesTypeKey(speciesTypekey);
                        info.setCompressMode(compressMode);
                        createGff.creategff4QTL(info);

                    }else{
                        throw new ArgumentsException("This Script requires '-mapKey: and -toDir:' as parameters");
                    }
                }else{
                    throw new ArgumentsException("This Script requires '-species:' as a parameter:\n" +
                            "-species:RAT/MOUSE/HUMAN");
                }
                break;

            case RgdId.OBJECT_KEY_SSLPS:
                CreateGff4SSLP create4Sslp = new CreateGff4SSLP();
                if( mapKey>0 && toDir!=null ){
                    CreateInfo info = new CreateInfo();
                    info.setMapKey(mapKey);
                    info.setToDir(toDir);
                    info.setSpeciesTypeKey(speciesTypekey);
                    info.setCompressMode(compressMode);
                    create4Sslp.createGff4Markers(info);

                }else{
                    throw new ArgumentsException("This Script requires '-mapKey: -fromFile: -toFile:' as parameters");
                }
                break;

            case RgdId.OBJECT_KEY_PROMOTERS:
                CreatePromoters4Gene createPromoters4Gene = new CreatePromoters4Gene();

                if( toDir!=null ){
                    createPromoters4Gene.setToDir(toDir);
                    createPromoters4Gene.createGenomicElements(compressMode);
                }else{
                    throw new ArgumentsException("This Script requires -toDir:' as parameters");
                }
                break;

            case OBJECT_KEY_DB_SNP:
                CreateGff4DbSnp createGff4DbSnp = new CreateGff4DbSnp();
                if( mapKey>0 && toFile!=null && speciesTypekey!=0 ){
                    createGff4DbSnp.setMapKey(mapKey);
                    createGff4DbSnp.setToFile(toFile);
                    createGff4DbSnp.setSpeciesTypeKey(speciesTypekey);
                    createGff4DbSnp.setBuild(build);
                    createGff4DbSnp.run(compressMode);
                }else{
                    throw new ArgumentsException("This Script requires '-mapKey: -species: -toFile:' as parameters");
                }
                break;

            default:
                throw new ArgumentsException("Invalid Object Type Key found!");
        }
    }

    boolean parseArgs(String[] args, DefaultListableBeanFactory bf) throws Exception {

        if(args.length>0){
            for(String obj: args ){
                String argArr[] = obj.split(":");
                if(obj.startsWith("-object:")){

                    switch (argArr[1]) {
                        // all species // assemblies version
                        case "genes":
                            CreateGff4Gene gm = (CreateGff4Gene) (bf.getBean("geneManager"));
                            gm.run();
                            return true;
                        case "ensembl":
                            CreateGff4Ensembl ensm = (CreateGff4Ensembl) (bf.getBean("ensemblManager"));
                            ensm.run();
                            return true;
                        case "jb2_eva":
                            JBrowse2PrepEva ee = (JBrowse2PrepEva) (bf.getBean("jbrowse2PrepEva"));
                            ee.run();
                            return true;
                        case "variants":
                            CreateGff4Variants varm = (CreateGff4Variants) (bf.getBean("variantManager"));
                            varm.run();
                            return true;
                        case "phenotypic_variants":
                            CreateGff4PhenotypicVariants pvarm = (CreateGff4PhenotypicVariants) (bf.getBean("phenotypicVariantManager"));
                            pvarm.run();
                            return true;
                        case "qtls":
                            CreateGff4QTL qm = (CreateGff4QTL) (bf.getBean("qtlManager"));
                            qm.run();
                            return true;
                        case "markers":
                            CreateGff4SSLP qs = (CreateGff4SSLP) (bf.getBean("markerManager"));
                            qs.run();
                            return true;
                        case "proteinDomains":
                            CreateGff4ProteinDomains pdm = (CreateGff4ProteinDomains) (bf.getBean("proteinDomainManager"));
                            pdm.run();
                            return true;
                        case "strains":
                            CreateGff4CongenicStrains s = (CreateGff4CongenicStrains) (bf.getBean("strainManager"));
                            s.run();
                            return true;
                        case "Eva":
                            CreateGff4Eva em = (CreateGff4Eva) (bf.getBean("evaManager"));
                            em.run();
                            return true;
                        case "diseases":
                            CreateGff4Ontology pdo = (CreateGff4Ontology) (bf.getBean("ontologyManager"));
                            pdo.runDiseaseOntology();
                            return true;
                        case "chebi":
                            CreateGff4Ontology pdw = (CreateGff4Ontology) (bf.getBean("ontologyManager"));
                            pdw.runChebiOntology();
                            return true;
                        case "aliases":
                            JBrowse2Aliases jba = new JBrowse2Aliases();
                            jba.run();
                            return true;
                        case "vep":
                            CreateGff4Vep vep = (CreateGff4Vep) (bf.getBean("vepManager"));
                            vep.run();
                            return true;
                        case "clinvar": {
                            CreateGff4ClinVar clinvar = (CreateGff4ClinVar) (bf.getBean("clinvarManager"));
                            clinvar.run( 17, "data/Variants/ClinVar/GRCh37.p13", Gff3ColumnWriter.COMPRESS_MODE_BGZIP);
                            clinvar.run( 38, "data/Variants/ClinVar/GRCh38.p14", Gff3ColumnWriter.COMPRESS_MODE_BGZIP);
                            return true;
                        }
                        case "gtfGenes": {
                            CreateGtf4Gene gtf = (CreateGtf4Gene) (bf.getBean("geneGtfWriter"));
                            gtf.run();
                            return true;
                        }

                        case "gene":
                            objectTypeKey = RgdId.OBJECT_KEY_GENES;
                            break;
                        case "qtl":
                            objectTypeKey = RgdId.OBJECT_KEY_QTLS;
                            break;
                        case "sslp":
                            objectTypeKey = RgdId.OBJECT_KEY_SSLPS;
                            break;
                        case "chromosome":
                            objectTypeKey = RgdId.OBJECT_KEY_CHROMOSOME;
                            break;
                        case "promoter":
                            objectTypeKey = RgdId.OBJECT_KEY_PROMOTERS;
                            break;
                        case "variant":
                            objectTypeKey = RgdId.OBJECT_KEY_VARIANTS;
                            break;
                        case "dbSnp":
                            objectTypeKey = OBJECT_KEY_DB_SNP;
                            break;
                    }
                }else
                if(obj.startsWith("-species:")){
                    speciesTypekey = SpeciesType.parse(argArr[1]);
                }else
                if(obj.startsWith("-mapKey:")){
                    mapKey = Integer.parseInt(argArr[1]);
                }else
                if(obj.startsWith("-mapKey2:")){
                    mapKey2 = Integer.parseInt(argArr[1]);
                }else
                if(obj.startsWith("-toFile:")){
                    toFile = argArr[1];
                }else
                if(obj.startsWith("-toDir:")){
                    toDir = argArr[1];
                }else
                if(obj.startsWith("-chr:")){
                    setChromosomes(expandChromosomeString(argArr[1]));
                }else
                if(obj.startsWith("-sampleID:")){
                    sampleID = Integer.parseInt(argArr[1]);
                }else
                if(obj.startsWith("-patientID:")){
                    patientID = Integer.parseInt(argArr[1]);
                }else
                if(obj.startsWith("-source:")){
                    source = argArr[1];
                }else
                if(obj.startsWith("-build:")){
                    build = argArr[1];
                }else
                if(obj.startsWith("-compress")){
                    compressMode = Gff3ColumnWriter.COMPRESS_MODE_BGZIP;
                }else
                if(obj.startsWith("-flavor:")){
                    flavor = argArr[1];
                }
            }
        }else{
            throw new ArgumentsException("Arguments not found..this script needs parameters");
        }

        // species post-processing: in case the species type is not parsable
        if( speciesTypekey<=0 && (mapKey>0 || mapKey2>0) ) {
            if( mapKey>0 ) {
                Map map = MapManager.getInstance().getMap(mapKey);
                speciesTypekey = map.getSpeciesTypeKey();
            }

            if( mapKey2>0 ) {

                Map map = MapManager.getInstance().getMap(mapKey2);
                int speciesTypeKey2 = map.getSpeciesTypeKey();
                if( speciesTypekey<=0 ) {
                    speciesTypekey = speciesTypeKey2;
                } else {
                    if( speciesTypekey!=speciesTypeKey2 ) {
                        throw new Exception("ERROR: 2 map keys specified, mapped to different species");
                    }
                }
            }
        }
        return false;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public class ArgumentsException extends Exception{
        public ArgumentsException(String msg){
            super(msg);
        }
    }

    public List<String> getChromosomes() {
        return chromosomes;
    }

    public void setChromosomes(List<String> chromosomes) {
        this.chromosomes = chromosomes;
    }

    // expand chromosome string into list of chromosomes
    // f.e. 1-19,X,Y      becomes      1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y
    List<String> expandChromosomeString(String chr) {

        List<String> chromosomes = new ArrayList<>();
        for( String range: chr.split("[,]") ) {
            int dashPos = range.indexOf("-");
            if( dashPos>0 ) {
                int start = Integer.parseInt(range.substring(0, dashPos));
                int stop = Integer.parseInt(range.substring(dashPos + 1));
                for( int i=start; i<=stop; i++ ) {
                    chromosomes.add(Integer.toString(i));
                }
            }
            else {
                chromosomes.add(range);
            }
        }
        return chromosomes;
    }

    public int getCompressMode() {
        return compressMode;
    }

    public void setCompressMode(int compressMode) {
        this.compressMode = compressMode;
    }

    public java.util.Map<Integer, String> getAssemblies() {
        return assemblies;
    }

    public void setAssemblies(java.util.Map<Integer, String> assemblies) {
        this.assemblies = assemblies;
    }
}
