package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
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

        CounterPool counters = new CounterPool();

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
            createGffFromQtlObject(qtlObj, gff3Writer, RATMINEgff3Writer, rgdInfoList, info.getMapKey(), info.getSpeciesTypeKey(), counters);
        }

        System.out.println("Number of Active Qtls processed:"+ counters.get("activeQtlCount"));
        System.out.println("Number of Qtls with Map positions:"+ counters.get("qtlsWithMapPos"));
        System.out.println("Qtls with NO Map positions:"+ counters.get("qtlsNoMapPos"));
        System.out.println("Qtls with more than one Map position:"+ counters.get("qtlsMoreThanOneMapPos"));
        System.out.println("Qtl Mapping Method stats:");
        System.out.println("-Positioned by flanking markers:"+counters.get("qtlsFlanking"));
        System.out.println("-Positioned by flanking marker and peak marker:"+counters.get("qtlsFlankingPeak"));
        System.out.println("-Positioned by peak marker only:"+counters.get("qtlsPeakOnly"));
        System.out.println("-Positioned by single flanking marker only:"+counters.get("qtlsSingleFlanking"));
        System.out.println("-Positioned by peak marker with size adjusted to avg qtl size for species:"+counters.get("qtlsPeakWithSizeAdjusted"));
        System.out.println("-Position imported from external source:"+counters.get("qtlsImportedFromExternal"));
        System.out.println("-NO Map Position Method ID:"+counters.get("noMapsPosMethodId"));
        System.out.println("\nQtls related to Strains:" + counters.get("qtlsWithRelStrains"));
        System.out.println("Qtls NOT related to Strains:" + counters.get("qtlsWithNoRelStrains"));
        System.out.println("Qtls related to Genes:" + counters.get("qtlsWithRelGenes"));
        System.out.println("Qtls NOT related to Genes:" + counters.get("qtlswithNoRelGenes"));
        System.out.println("Qtls related to Qtls:" + counters.get("qtlsWithRelQtls"));
        System.out.println("Qtls NOT related to Qtls:" + counters.get("qtlswithNoRelQtls"));
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
                                       List<RGDInfo> rgdInfoList, int mapKey, int speciesTypeKey, CounterPool counters) throws Exception {

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

        counters.increment("activeQtlCount");

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
            counters.increment("qtlsMoreThanOneMapPos");
        }else if(mdList.size()==0){
            counters.increment("qtlsNoMapPos");
        }

        if(mdList.size()>0){

            counters.increment("qtlsWithMapPos");

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
                    String[] mapMethArr = md.getMapPositionMethod().split(" - ");

                    mappingMethod = mapMethArr[1];
                }else{
                    mappingMethod = md.getMapPositionMethod();
                }


                switch(md.getMapsDataPositionMethodId()){
                    case 1 :
                        counters.increment("qtlsFlanking");
                        break;
                    case 2 :
                        counters.increment("qtlsFlankingPeak");
                        break;
                    case 3 :
                        counters.increment("qtlsPeakOnly");
                        break;
                    case 4 :
                        counters.increment("qtlsSingleFlanking");
                        break;
                    case 5 :
                        counters.increment("qtlsPeakWithSizeAdjusted");
                        break;
                    case 6 :
                        counters.increment("qtlsImportedFromExternal");
                        break;
                    default:
                        counters.increment("noMapsPosMethodId");
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

                        counters.increment("qtlsWithRelStrains");

                        for(Strain st : relatedStrainsList){
                            String smbl = st.getSymbol();
                            smbl = smbl.replaceAll(":","_");
                            relStrain += URLEncoder.encode(smbl,"UTF-8")+":"+st.getRgdId()+",";
                        }
                    }else{

                        counters.increment("qtlsWithNoRelStrains");

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

                    counters.increment("qtlsWithRelGenes");

                    for(Gene g: relatedGenesList){
                        relGenes += g.getSymbol()+":"+g.getRgdId()+",";
                    }
                }else{

                    counters.increment("qtlswithNoRelGenes");

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
                        counters.increment("qtlsWithRelQtls");

                        for(int keys: relQtlKeys){
                           relQtls += keys+":"+relatedQtlMap.get(keys)+",";
                        }
                    }else
                    if(relatedQtlMap.size()==0){

                        counters.increment("qtlswithNoRelQtls");

                        relQtls = "NA";
                    }
                    if(relQtls.endsWith(",")){
                        relQtls = relQtls.substring(0,relQtls.length()-1);
                    }
                }else{
                    relQtls = "NA";
                }

                attributesHashMap.put("ID", rgdId +"_"+start+"_"+stop);
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
