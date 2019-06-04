package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.impl.TranscriptDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.SeqUtils;
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
    private String mapKey;
    private int objectTypeKey;
    private List<String> chromosomes;
    private String toFile;
    private String fromFile;
    int sampleID;
    int patientID;
    String build;
    String source="RGD_CNRat";
    String ont_aspect;
    boolean compress = false;
    String flavor;
    private String version;

    public static final int OBJECT_KEY_DB_SNP = -1;

    public String getUsage() {
        return usage;
    }

    public void setUsage(String usage) {
        this.usage = usage;
    }

    private String usage;


    static public void main(String[] args) throws Exception {
        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        Manager creator = (Manager) (bf.getBean("manager"));

        System.out.println(creator.getVersion());

        try{
            System.out.println("starting gff3 pipeline.. ");

            String codeUsage =
                    "USAGE for RGD OBjects:\n" +
                        "-object:gene/qtl/sslp/chromosome\n" +
                        "-species:RAT/HUMAN/MOUSE\n Other Parameters may include:\n" +
                        "-mapKey:360/70/60/18/13/35/17\n" +
                        "-chr:<chromosome-range>\n"+
                        "in Rat: chromosomes = 1-20,X,Y\n" +
                        "in Mouse: chromosomes = 1-19,X,Y\n"+
                        "in Human: chromosomes = 1-22,X,Y\n" +
                        "-toFile:../log/RGDGFF3/Output/\n"+
                        "-fromFile:/home/rgddata/data/RGDGFF3/Input/\n" +
                        "-compress   (compress output files with gzip) \n" +
                        "\n---------------------------------------------------\n" +
                    "USAGE for CarpeNovo Variants:\n" +
                        "-sampleID:<sampleId>\n" +
                        "-patientID:<patientId>,  f.e. -patientID:500\n" +
                        "-source:(this is freeText..byDefault it is CN\n" +
                        "-toFile:../log/RGDGFF3/Output/\n"+
                        "-chr:<chromosome-range>\n"+
                        "in Rat: chromosomes = 1-20,X\n" +
                        "in Mouse: chromosomes = 1-19,X,Y\n"+
                        "in Human: chromosomes = 1-22,X,Y\n" +
                        "Note: either sampleID or patientID is mandatory\n" +
                        "\n-------------------------------------------------------\n";
            creator.setUsage(codeUsage);
            creator.doMain(args);
        }catch (Exception createNewGffException){
            createNewGffException.printStackTrace();
        }
    }

    public void doMain(String[] args) throws Exception {
        parseArgs(args);

        if(objectTypeKey!=0){
            handleObjects();
        }
        else if(sampleID!=0 || patientID!=0){
            if( source==null || toFile==null || getChromosomes()==null ){
                throw new ArgumentsException("This script requires '-source:' '-toFile:' '-chr:' as parameters:\n" +
                           getUsage());
            }

            CreateGff4CarpeNovo createGffCarpe = new CreateGff4CarpeNovo();
            createGffCarpe.setFileSource(source);
            createGffCarpe.setToFile(toFile);
            createGffCarpe.setChromosomes(getChromosomes());
            if( sampleID!=0 )
                createGffCarpe.createGff3ForSample(sampleID);
            if( patientID!=0 )
                createGffCarpe.createGff3ForPatient(patientID);

        }
        else if( ont_aspect!=null ){
            if( speciesTypekey!=0 && mapKey!=null && toFile!=null && getChromosomes()!=null ){
                CreateGff4Ontology createGff4Ont = new CreateGff4Ontology();
                createGff4Ont.setSpeciesTypeKey(speciesTypekey);
                createGff4Ont.setToFile(toFile);
                createGff4Ont.setMapKey(Integer.parseInt(mapKey));
                createGff4Ont.setOntAspect(ont_aspect);
                createGff4Ont.setChromosomes(getChromosomes());
                createGff4Ont.run(compress);
            }else{
                throw new ArgumentsException("This script requires '-speciesTypeKey:','-mapKey:','-toFile:','-chr:' as parameters:\n" +
                        getUsage());
            }

        }
        else if( flavor.equals("ensembl_prep") ) {
            EnsemblPrep.run(fromFile);

        } else {
            throw new ArgumentsException("This script requires '-object:' or '-sampleID' or '-ontAspect' as a parameter:\n" +
                    "-sampleID:309/329/330\t-object:gene/qtl/sslp/strain\t-ontAspect:D/M/P/W/..\n" + getUsage());
        }
    }

    void handleObjects() throws Exception {
        switch(objectTypeKey) {
            case RgdId.OBJECT_KEY_GENES:
                if((mapKey!=null)&&(getChromosomes()!=null)&&(toFile!=null)&&(speciesTypekey!=0)){
                    if( flavor==null ) {
                        CreateGff4Gene createGff = new CreateGff4Gene();
                        createGff.setMapKey(Integer.parseInt(mapKey));
                        createGff.setChromosomes(getChromosomes());
                        createGff.setNewPathGff3(toFile);
                        createGff.setSpeciesTypeKey(speciesTypekey);
                        createGff.createGeneGff3(compress);
                    } else if( flavor.equals("AGR") ){
                        CreateGff4GeneAgr createGff = new CreateGff4GeneAgr();
                        createGff.setMapKey(Integer.parseInt(mapKey));
                        createGff.setGff3Path(toFile);
                        createGff.setSpeciesTypeKey(speciesTypekey);
                        createGff.createGeneGff3(compress);
                    } else {
                        throw new ArgumentsException("unknown gene flavor: "+flavor);
                    }
                }else{
                    throw new ArgumentsException("This Script requires '-mapKey:','-toFile:'," +
                            "'-chr:','-species:' as parameters:\n" + getUsage());
                }
                break;

            case RgdId.OBJECT_KEY_QTLS:
                if(speciesTypekey!=0){
                    CreateGff4QTL createGff = new CreateGff4QTL();
                    if((mapKey!=null)&&(toFile!=null)){
                        createGff.setMapKey(Integer.parseInt(mapKey));
                        createGff.setToFile(toFile);
                        createGff.setSpeciesTypeKey(speciesTypekey);
                        createGff.creategff4QTL(compress);

                    }else{
                        throw new ArgumentsException("This Script requires '-mapKey: and -toFile:' " +
                                "as parameters:\n" + getUsage());
                    }
                }else{
                    throw new ArgumentsException("This Script requires '-species:' as a parameter:\n" +
                            "-species:RAT/MOUSE/HUMAN\n" + getUsage());
                }
                break;

            case RgdId.OBJECT_KEY_SSLPS:
                CreateGff4SSLP create4Sslp = new CreateGff4SSLP();
                if( mapKey!=null && toFile!=null ){
                    create4Sslp.setToFile(toFile);
                    create4Sslp.setMapKey(Integer.parseInt(mapKey));
                    create4Sslp.setSpeciesTypeKey(speciesTypekey);
                    create4Sslp.creategff4sslps(compress);

                }else{
                    throw new ArgumentsException("This Script requires '-mapKey: -fromFile: -toFile:' " +
                            "as parameters:\n" + getUsage());
                }
                break;

            case RgdId.OBJECT_KEY_STRAINS:
                CreateGff4CongenicStrains create4Strains = new CreateGff4CongenicStrains();

                if( mapKey!=null && toFile!=null ){
                    create4Strains.setMap_key(Integer.parseInt(mapKey));
                    create4Strains.setToFile(toFile);
                    create4Strains.creategff4CongenicStrains(compress);

                }else{
                    throw new ArgumentsException("This Script requires '-mapKey: -toFile:' as parameter:\n" +
                            getUsage());
                }
                break;

            case RgdId.OBJECT_KEY_PROMOTERS:
                CreatePromoters4Gene createPromoters4Gene = new CreatePromoters4Gene();

                if( mapKey!=null && toFile!=null && speciesTypekey!=0 ){
                    createPromoters4Gene.setMapKey(Integer.parseInt(mapKey));
                    createPromoters4Gene.setObjectKey(objectTypeKey);
                    createPromoters4Gene.setToFile(toFile);
                    createPromoters4Gene.setSpeciesTypekey(speciesTypekey);
                    createPromoters4Gene.createGenomicElements(compress);
                }else{
                    throw new ArgumentsException("This Script requires '-mapKey: -species: -toFile:' " +
                            "as parameters:\n" + getUsage());
                }
                break;

            case RgdId.OBJECT_KEY_VARIANTS:
                CreateGff4ClinVar creator = new CreateGff4ClinVar();
                if( mapKey!=null && toFile!=null && speciesTypekey!=0 ){
                    creator.setMapKey(Integer.parseInt(mapKey));
                    creator.setToFile(toFile);
                    creator.setSpeciesTypeKey(speciesTypekey);
                    creator.run(compress);
                }else{
                    throw new ArgumentsException("This Script requires '-mapKey: -species: -toFile:' " +
                            "as parameters:\n" + getUsage());
                }
                break;


            case OBJECT_KEY_DB_SNP:
                CreateGff4DbSnp createGff4DbSnp = new CreateGff4DbSnp();
                if( mapKey!=null && toFile!=null && speciesTypekey!=0 ){
                    createGff4DbSnp.setMapKey(Integer.parseInt(mapKey));
                    createGff4DbSnp.setToFile(toFile);
                    createGff4DbSnp.setSpeciesTypeKey(speciesTypekey);
                    createGff4DbSnp.setBuild(build);
                    createGff4DbSnp.run(compress);
                }else{
                    throw new ArgumentsException("This Script requires '-mapKey: -species: -toFile:' " +
                            "as parameters:\n" + getUsage());
                }
                break;

            default:
                throw new ArgumentsException("Invalid Object Type Key found!:\n" + getUsage());
        }
    }

    void parseArgs(String[] args) throws Exception {

        if(args.length>0){
            for(String obj: args ){
                String argArr[] = obj.split(":");
                if(obj.startsWith("-object:")){

                    switch (argArr[1]) {
                        case "gene":
                            objectTypeKey = RgdId.OBJECT_KEY_GENES;
                            break;
                        case "qtl":
                            objectTypeKey = RgdId.OBJECT_KEY_QTLS;
                            break;
                        case "strain":
                            objectTypeKey = RgdId.OBJECT_KEY_STRAINS;
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
                    mapKey = argArr[1];
                }else
                if(obj.startsWith("-toFile:")){
                    toFile = argArr[1];
                }else
                if(obj.startsWith("-fromFile:")){
                    fromFile=argArr[1];
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
                if(obj.startsWith("-ontAspect:")){
                    ont_aspect = argArr[1];
                }else
                if(obj.startsWith("-build:")){
                    build = argArr[1];
                }else
                if(obj.startsWith("-compress")){
                    compress = true;
                }else
                if(obj.startsWith("-flavor:")){
                    flavor = argArr[1];
                }
            }
        }else{
            throw new ArgumentsException("Arguments not found..this script needs parameters.\n" + getUsage());
        }
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

    public boolean isCompress() {
        return compress;
    }

    public void setCompress(boolean compress) {
        this.compress = compress;
    }
}


