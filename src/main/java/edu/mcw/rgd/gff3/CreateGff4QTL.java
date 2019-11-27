package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.mapping.MapManager;

import java.io.PrintWriter;
import java.net.URLEncoder;
import java.util.*;
import java.util.Map;

/**
 * @author BBakir
 * Date: Aug 11, 2008
 */
public class CreateGff4QTL {

    private RgdGff3Dao dao = new RgdGff3Dao();
    private List<String> processedAssemblies;

    // counts
    int activeQtlCount;
    int qtlsWithMapPos;
    int qtlsMoreThanOneMapPos;
    int qtlsNoMapPos;
    int qtlsFlanking;
    int qtlsFlankingPeak;
    int qtlsPeakOnly;
    int qtlsSingleFlanking;
    int qtlsPeakWithSizeAdjusted;
    int qtlsImportedFromExternal;
    int noMapsPosMethodId;
    int qtlsWithRelStrains;
    int qtlsWithNoRelStrains;
    int qtlsWithRelGenes;
    int qtlswithNoRelGenes;
    int qtlsWithRelQtls;
    int qtlswithNoRelQtls;

    void clearCounts() {
        activeQtlCount=0;
        qtlsWithMapPos=0;
        qtlsMoreThanOneMapPos=0;
        qtlsNoMapPos=0;
        qtlsFlanking=0;
        qtlsFlankingPeak=0;
        qtlsPeakOnly=0;
        qtlsSingleFlanking=0;
        qtlsPeakWithSizeAdjusted=0;
        qtlsImportedFromExternal=0;
        noMapsPosMethodId=0;
        qtlsWithRelStrains=0;
        qtlsWithNoRelStrains=0;
        qtlsWithRelGenes=0;
        qtlswithNoRelGenes=0;
        qtlsWithRelQtls=0;
        qtlswithNoRelQtls=0;
    }

    /**
     * load the species list and assemblies from properties/AppConfigure.xml
     */
    public void run() throws Exception {
        for( String assemblyInfo: processedAssemblies ) {
            CreateInfo info = new CreateInfo();
            info.parseFromString(assemblyInfo);

            creategff4QTL(info);
        }
    }

    /**
     * create Gff3 for Rat, Human and Mouse.. only for Rat it will print out related strains
     * @throws Exception
     */
    public void creategff4QTL(CreateInfo info) throws Exception {

        clearCounts();

        String speciesName = SpeciesType.getCommonName(info.getSpeciesTypeKey());

        System.out.println("START QTL GFF3 Generator for  "+speciesName+"  MAP_KEY="+info.getMapKey()+"  ASSEMBLY "+ MapManager.getInstance().getMap(info.getMapKey()).getName());
        System.out.println("========================");

        String gffFile = info.getToDir()+speciesName+"_RGDQTLS.gff3";
        String RATMINEGffFile = info.getToDir()+speciesName+"_RATMINE_RGDQTLS.gff3";

        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gffFile, false, info.isCompress());
        Gff3ColumnWriter RATMINEgff3Writer = new Gff3ColumnWriter(RATMINEGffFile, false, info.isCompress());
        RATMINEgff3Writer.setRatmineCompatibleFormat(true);

        PrintWriter densityFile = new PrintWriter(info.getToDir()+speciesName+"_density.gff");

        List<QTL> qtlList = dao.getActiveQTLs(info.getSpeciesTypeKey());
        List<RGDInfo> rgdInfoList = new ArrayList<RGDInfo>();

        for(QTL qtlObj: qtlList){
            createGffFromQtlObject(qtlObj, gff3Writer, RATMINEgff3Writer, rgdInfoList, info.getMapKey(), info.getSpeciesTypeKey());
        }

        System.out.println("Number of Active Qtls processed:"+ activeQtlCount);
        System.out.println("Number of Qtls with Map positions:"+ qtlsWithMapPos);
        System.out.println("Qtls with NO Map positions:"+ qtlsNoMapPos);
        System.out.println("Qtls with more than one Map position:"+ qtlsMoreThanOneMapPos);
        System.out.println("Qtl Mapping Method stats:");
        System.out.println("-Positioned by flanking markers:"+qtlsFlanking);
        System.out.println("-Positioned by flanking marker and peak marker:"+qtlsFlankingPeak);
        System.out.println("-Positioned by peak marker only:"+qtlsPeakOnly);
        System.out.println("-Positioned by single flanking marker only:"+qtlsSingleFlanking);
        System.out.println("-Positioned by peak marker with size adjusted to avg qtl size for species:"+qtlsPeakWithSizeAdjusted);
        System.out.println("-Position imported from external source:"+qtlsImportedFromExternal);
        System.out.println("-NO Map Position Method ID:"+noMapsPosMethodId);
        System.out.println("\nQtls related to Strains:" + qtlsWithRelStrains);
        System.out.println("Qtls NOT related to Strains:" + qtlsWithNoRelStrains);
        System.out.println("Qtls related to Genes:" + qtlsWithRelGenes);
        System.out.println("Qtls NOT related to Genes:" + qtlswithNoRelGenes);
        System.out.println("Qtls related to Qtls:" + qtlsWithRelQtls);
        System.out.println("Qtls NOT related to Qtls:" + qtlswithNoRelQtls);
        System.out.println("\nGFF3 File SUCCESSFUL!\n");

         //close file
        gff3Writer.close();
        RATMINEgff3Writer.close();

        CalculateDensity calcDensity = new CalculateDensity();
        calcDensity.setRgdInfoList(rgdInfoList);
        calcDensity.setDensityWriter(densityFile);

        boolean verbose = false;
        calcDensity.runDensityCalculator(verbose);
    }

    public void createGffFromQtlObject(QTL qtlObject, Gff3ColumnWriter gff3Writer, Gff3ColumnWriter RATMINEgff3Writer,
                                       List<RGDInfo> rgdInfoList, int mapKey, int speciesTypeKey) throws Exception {

        int qtlKey = qtlObject.getKey();
        int rgdId = qtlObject.getRgdId();
        String type="QTL";
        String source="RGD";
        String full_name = qtlObject.getName();
        String symbol = qtlObject.getSymbol();
        String mappingMethod;
        String chrom;
        Integer start;
        Integer stop;
        String strand=".";
        String lod;
        String pValue;
        String notes;

        activeQtlCount++;

        List<MapData> mdList = dao.getMapData(rgdId,mapKey);


        if(qtlObject.getLod()!=null){
            lod=String.valueOf(qtlObject.getLod());
        }else{
            lod = "null";
        }
        if(qtlObject.getPValue()!=null){
            pValue=String.valueOf(qtlObject.getPValue());
        }else{
            pValue="null";
        }
        if(qtlObject.getNotes()!=null){
            notes=qtlObject.getNotes().replaceAll(";"," ");
            //notes = URLEncoder.encode(notes, "UTF-8");
        }else{
            notes = "null";
        }


        if(mdList.size()>1){
            qtlsMoreThanOneMapPos++;
        }else if(mdList.size()==0){
            qtlsNoMapPos++;
        }

        if(mdList.size()>0){

            qtlsWithMapPos++;

            for(MapData md: mdList){
                start = md.getStartPos();
                stop = md.getStopPos();
                chrom = md.getChromosome();
                RGDInfo rgdInf = new RGDInfo();
                rgdInf.setRgdId(rgdId);
                rgdInf.setChromosome(chrom);
                rgdInf.setStartPos(start);
                rgdInf.setStopPos(stop);

                //add object to rgdInfoList
                rgdInfoList.add(rgdInf);

                if(!(md.getStrand()==null)){
                    strand = md.getStrand();
                }

                if(md.getMapPositionMethod().contains("-")){
                    String mapMethArr[] = md.getMapPositionMethod().split(" - ");

                    mappingMethod = mapMethArr[1];
                }else{
                    mappingMethod = md.getMapPositionMethod();
                }


                switch(md.getMapsDataPositionMethodId()){
                    case 1 :
                        //Statements
                        qtlsFlanking++;
                        break; //optional
                    case 2 :
                        //Statements
                        qtlsFlankingPeak++;
                        break; //optional
                    case 3 :
                        //Statements
                        qtlsPeakOnly++;
                        break; //optional
                    case 4 :
                        //Statements
                        qtlsSingleFlanking++;
                        break; //optional
                    case 5 :
                        //Statements
                        qtlsPeakWithSizeAdjusted++;
                        break; //optional
                    case 6 :
                        //Statements
                        qtlsImportedFromExternal++;
                        break; //optional
                    //You can have any number of case statements.
                    default : //Optional
                        noMapsPosMethodId++;
                       //Statements
                }

                //write first 8 columns
                gff3Writer.writeFirst8Columns(chrom, source, type, start, stop, ".",strand,".");
                RATMINEgff3Writer.writeFirst8Columns(chrom, source, type, start, stop, ".",strand,".");


                //initialize hashmap for attributes.
                Map<String, String> attributesHashMap = new HashMap<String, String>();
                Map<String, String> RATMINEattributesHashMap = new HashMap<String, String>();


                //get related strains
                String relStrain="";
                if(speciesTypeKey==3){
                    List<Strain> relatedStrainsList = dao.getStrainAssociationsForQtl(rgdId);

                    if(relatedStrainsList.size()>0){

                        qtlsWithRelStrains++;

                        for(Strain st : relatedStrainsList){
                            String smbl = st.getSymbol();
                            smbl = smbl.replaceAll(":","_");
                            relStrain += URLEncoder.encode(smbl,"UTF-8")+":"+st.getRgdId()+",";
                        }
                    }else{

                        qtlsWithNoRelStrains++;

                        relStrain = "NA";
                    }
                    if(relStrain.endsWith(",")){
                        relStrain = relStrain.substring(0,relStrain.length()-1);
                    }
                }


                //get related genes.
                String relGenes="";
                List<Gene> relatedGenesList = dao.getGeneAssociationsForQtl(rgdId);
                if(relatedGenesList.size()>0){

                    qtlsWithRelGenes++;

                    for(Gene g: relatedGenesList){
                        relGenes += g.getSymbol()+":"+g.getRgdId()+",";
                    }
                }else{

                    qtlswithNoRelGenes++;

                    relGenes="NA";
                }
                if(relGenes.endsWith(",")){
                    relGenes = relGenes.substring(0,relGenes.length()-1);
                }


                //get qtl - qtl associations.
                Map<Integer, String> relatedQtlMap = dao.getQtlToQtlAssociations(qtlKey);

                String relQtls="";
                if(relatedQtlMap!=null){
                    if(relatedQtlMap.size()>0){
                        Set<Integer> relQtlKeys = relatedQtlMap.keySet();
                        qtlsWithRelQtls++;

                        for(int keys: relQtlKeys){
                           relQtls += keys+":"+relatedQtlMap.get(keys)+",";
                        }
                    }else
                    if(relatedQtlMap.size()==0){

                        qtlswithNoRelQtls++;

                        relQtls = "NA";
                    }
                    if(relQtls.endsWith(",")){
                        relQtls = relQtls.substring(0,relQtls.length()-1);
                    }
                }else{
                    relQtls = "NA";
                }

                attributesHashMap.put("ID",String.valueOf(rgdId)+"_"+start+"_"+stop);
                attributesHashMap.put("Name", "QTL:"+symbol);
                attributesHashMap.put("fullName", full_name);
                attributesHashMap.put("Alias", "RGD:"+rgdId+", QTL:"+full_name+","+symbol);
                attributesHashMap.put("LOD",lod);
                attributesHashMap.put("pValue",pValue);
                attributesHashMap.put("Note",notes);
                attributesHashMap.put("Dbxref","RGD:"+rgdId);
                attributesHashMap.put("mappingMethod",mappingMethod);
                attributesHashMap.put("relatedQTLs",relQtls);
                if(speciesTypeKey==3){
                    attributesHashMap.put("relatedStrains", relStrain);
                }
                attributesHashMap.put("relatedGenes",relGenes);
                attributesHashMap.put("Index","1");

                gff3Writer.writeAttributes4Gff3(attributesHashMap);


                RATMINEattributesHashMap.put("ID",String.valueOf(rgdId));
                RATMINEattributesHashMap.put("Name", symbol);
                RATMINEattributesHashMap.put("fullName", full_name);
                RATMINEattributesHashMap.put("Alias", "RGD:"+rgdId+", QTL:"+full_name);
                RATMINEattributesHashMap.put("LOD",lod);
                RATMINEattributesHashMap.put("pValue",pValue);
                RATMINEattributesHashMap.put("Note",notes);
                RATMINEattributesHashMap.put("Dbxref","RGD:"+rgdId);
                RATMINEattributesHashMap.put("mappingMethod",mappingMethod);
                RATMINEattributesHashMap.put("relatedQTLs",relQtls);
                if(speciesTypeKey==3){
                    RATMINEattributesHashMap.put("relatedStrains", relStrain);
                }
                RATMINEattributesHashMap.put("relatedGenes",relGenes);
                RATMINEattributesHashMap.put("Index","1");

                RATMINEgff3Writer.writeAttributes4Gff3(RATMINEattributesHashMap);
            }
        }
    }

    public void setProcessedAssemblies(List processedAssemblies) {
        this.processedAssemblies = processedAssemblies;
    }

    public List getProcessedAssemblies() {
        return processedAssemblies;
    }
}
