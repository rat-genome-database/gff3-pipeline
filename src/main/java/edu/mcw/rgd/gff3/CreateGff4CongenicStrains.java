package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.Strain;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;

import java.util.*;

/**
 * @author pjayaraman
 * Date: Mar 23, 2011
 */
public class CreateGff4CongenicStrains {

    private RgdGff3Dao dao = new RgdGff3Dao();
    private List<String> processedAssemblies;

    final String source = "RGD_Rat_Strain";
    
    /**
     * load the species list and assemblies from properties/AppConfigure.xml
     */
    public void run() throws Exception {
        for( String assemblyInfo: processedAssemblies ) {
            CreateInfo info = new CreateInfo();
            info.parseFromString(assemblyInfo);

            creategff4CongenicStrains(info);
        }
    }

    public void creategff4CongenicStrains(CreateInfo info) throws Exception{
        CounterPool counters = new CounterPool();

        String speciesName = SpeciesType.getCommonName(info.getSpeciesTypeKey());
        String assemblySymbol = info.getAssemblySymbol()!=null ? info.getAssemblySymbol() : Gff3Utils.getAssemblySymbol(info.getMapKey());

        System.out.println("START Strain GFF3 Generator for  "+speciesName+"  MAP_KEY="+info.getMapKey()+"  ASSEMBLY "+ MapManager.getInstance().getMap(info.getMapKey()).getName());
        System.out.println("========================");

        String gffFileCongenic = info.getToDir()+"/"+assemblySymbol+"_CongenicStrains.gff3";
        String RATMINEGffFile = info.getToDir()+"/"+assemblySymbol+"_Strains_RATMINE.gff3";
        String gffFileMutant = info.getToDir()+"/"+assemblySymbol+"_MutantStrains.gff3";

        //initialization of the gff3writer
        Gff3ColumnWriter gff3WriterCongenic = new Gff3ColumnWriter(gffFileCongenic, false, info.isCompress());
        Gff3ColumnWriter gff3WriterRATMINE = new Gff3ColumnWriter(RATMINEGffFile, false, info.isCompress());
        Gff3ColumnWriter gff3WriterMutant = new Gff3ColumnWriter(gffFileMutant, false, info.isCompress());

        String header = "# RAT GENOME DATABASE (https://rgd.mcw.edu/)\n";
        header += "# Species: "+ speciesName+"\n";
        header += "# Assembly: "+ MapManager.getInstance().getMap(info.getMapKey()).getName()+"\n";
        header += "# Primary Contact: mtutaj@mcw.edu\n";
        header += "# Generated: "+new Date()+"\n";

        gff3WriterCongenic.print(header);
        gff3WriterMutant.print(header);


        // get all active strains having positions on the given assembly
        List<Strain> strainList = dao.getMappedStrains(info.getMapKey());

        //for each strain get the annotations from full annot table
        // and get the map data positions
        for(Strain strain: strainList){

            // we only have map positions for congenic and mutant strains
            String strainTypeLc = Utils.defaultString(strain.getStrainTypeName()).toLowerCase();
            boolean isCongenic = strainTypeLc.equals("congenic");
            boolean isMutant = strainTypeLc.equals("mutant")
                    || strainTypeLc.equals("mutant_knockin")
                    || strainTypeLc.equals("mutant_knockout")
                    || strainTypeLc.equals("transgenic");
            if( !(isCongenic || isMutant) ){
                continue;
            }
            if( isCongenic )
                counters.increment("congenicStrains");
            if( isMutant )
                counters.increment("mutantStrains");


            String note = strain.getOrigin();
            if(note!=null){
                note = note.replace(';','-')
                           .replace(","," AND ")
                           .replace('\"','\'')  // per Aurash's request
                           .replace('\n',' ');
            }else{
                note = "NA";
            }

            String parentStr = Utils.NVL(strain.getStrain(), "NA");
            String subStr = Utils.NVL(strain.getSubstrain(), "NA");
            String src = Utils.NVL(strain.getSource(), "NA");

            Map<String,String> annotAttributes = processAnnotations(strain.getRgdId(), counters);


            List<MapData> newmdList = dao.getMapData(strain.getRgdId(), info.getMapKey());

            if(newmdList.size()>1){
                counters.increment("moreThanOneMapPosCount");
                if( isCongenic ) {
                    strainTypeLc = "multicongenic";
                }
            }else if(newmdList.size()==1){
                counters.increment("mappedStrainCount");
            }else if(newmdList.size()==0){
                counters.increment("noMapPos");
            }

            //for each new map data list
            for(int m=0; m<newmdList.size(); m++){
                MapData md = newmdList.get(m);

                String strand = Utils.NVL(md.getStrand(), ".");

                //write first each columns for this line
                if( isCongenic ) {
                    gff3WriterCongenic.writeFirst8Columns(md.getChromosome(), source, strainTypeLc,
                            md.getStartPos(),md.getStopPos(),".",strand,".");
                }
                if( isMutant ) {
                    gff3WriterMutant.writeFirst8Columns(md.getChromosome(), source, strainTypeLc,
                            md.getStartPos(),md.getStopPos(),".",strand,".");
                }
                gff3WriterRATMINE.writeFirst8Columns(md.getChromosome(), source, strainTypeLc,
                        md.getStartPos(),md.getStopPos(),".",strand,".");

                //create attributes hash map for this line
                Map<String, String> attributesHashMap = new HashMap<>();
                Map<String, String> RATMINE_attributesHashMap = new HashMap<>();

                emitAttributesForMultiCongenics(attributesHashMap, newmdList, md);

                attributesHashMap.put("ID",strain.getRgdId()+"_"+md.getChromosome()+"_"+m);
                attributesHashMap.put("Name", strain.getSymbol());
                attributesHashMap.put("Alias", "RGD:"+strain.getRgdId()+","+strain.getRgdId()+","+
                        "Strain:"+parentStr+","+subStr+","+strain.getSymbol());
                attributesHashMap.put("Note", note);
                attributesHashMap.put("origin", src);
                attributesHashMap.put("Index","1");

                RATMINE_attributesHashMap.put("ID",String.valueOf(strain.getRgdId()));
                RATMINE_attributesHashMap.put("Name", strain.getSymbol());
                RATMINE_attributesHashMap.put("Alias", "RGD"+strain.getRgdId()+","+strain.getRgdId()+","+
                        "Strain:"+parentStr+","+subStr);
                RATMINE_attributesHashMap.put("Note", note);
                RATMINE_attributesHashMap.put("source", src);
                RATMINE_attributesHashMap.put("Index","1");

                for( Map.Entry<String,String> entry: annotAttributes.entrySet() ) {
                    attributesHashMap.put(entry.getKey(), entry.getValue());
                    RATMINE_attributesHashMap.put(entry.getKey(), entry.getValue());
                }

                //write attributes into gff3 file
                if( isCongenic ) {
                    gff3WriterCongenic.writeAttributes4Gff3(attributesHashMap);
                }
                if( isMutant ) {
                    gff3WriterMutant.writeAttributes4Gff3(attributesHashMap);
                }
                gff3WriterRATMINE.writeAttributes4Gff3(RATMINE_attributesHashMap);
            }
        }

        gff3WriterCongenic.close();
        gff3WriterRATMINE.close();
        gff3WriterMutant.close();

        System.out.println("\nCongenic Strain records processed: "+ counters.get("congenicStrains"));
        System.out.println("Mutant Strain records processed: "+ counters.get("mutantStrains"));
        System.out.println("Total count of Strain having only one Map position: "+ counters.get("mappedStrainCount"));
        System.out.println("Number of Strains with more than one Map position: "+ counters.get("moreThanOneMapPosCount"));
        System.out.println("Number of Strains having NO Map positions: "+ counters.get("noMapPos"));
        System.out.println("Number of Strains with Annotation: "+ counters.get("strainsWithAnn"));
        System.out.println("Number of Strains with NO Annotations at all: "+ counters.get("strainsNoAnn"));
        System.out.println("Number of Strain Records with Disease annotations: "+ counters.get("strainsWithDiseaseAnn"));
        System.out.println("Number of Strain Records with NO Disease annotations: "+ counters.get("strainsNoDisAnn"));
        System.out.println("Number of Strain Records with Phenotype annotations: "+ counters.get("strainsWithPhenoAnn"));
        System.out.println("Number of Strain Records with NO Phenotype annotations: "+ counters.get("strainsNoPhenoAnn"));
        System.out.println("\nGFF3 File SUCCESSFUL!");
    }

    Map<String,String> processAnnotations(int strainRgdId, CounterPool counters) throws Exception {

        List<Annotation> disAnn = new ArrayList<>();
        List<Annotation> pheAnn = new ArrayList<>();

        List<Annotation> annList = dao.getAnnotations(strainRgdId);
        if( annList.isEmpty() ){
            counters.increment("strainsNoAnn");
            return Collections.emptyMap();
        }

        counters.increment("strainsWithAnn");
        Map<String,String> annotAttributes = new HashMap<>();

        for(Annotation annObj:annList){
            if(annObj.getTermAcc().startsWith("RDO")){
                disAnn.add(annObj);
            }else if(annObj.getTermAcc().startsWith("MP")){
                pheAnn.add(annObj);
            }
        }

        if(disAnn.size()!=0){
            counters.increment("strainsWithDiseaseAnn");

            String disOntValues="";
            for( Annotation disease: disAnn ){
                disOntValues += disease.getTermAcc()+":"+disease.getTerm()+",";
            }
            //remove last "," from string
            disOntValues = disOntValues.substring(0,disOntValues.length()-1);

            annotAttributes.put("DiseaseOntologyAssociation",disOntValues);
        }else {
            counters.increment("strainsNoDisAnn");
        }

        if(pheAnn.size()!=0){
            counters.increment("strainsWithPhenoAnn");

            String pheOntValues="";

            for( Annotation phenotype: pheAnn ){
                pheOntValues += phenotype.getTermAcc()+":"+phenotype.getTerm()+",";
            }
            //remove last "," from this line
            pheOntValues = pheOntValues.substring(0,pheOntValues.length()-1);

            annotAttributes.put("PhenotypeOntologyAssociation",pheOntValues);
        }else {
            counters.increment("strainsNoPhenoAnn");
        }

        return annotAttributes;
    }

    void emitAttributesForMultiCongenics(Map<String,String> attributes, List<MapData> mds, MapData md) {

        // emit positional information for additional loci
        int locus = 2;
        for( MapData md1: mds ) {
            if( md1==md )
                continue;
            attributes.put("loc"+locus, "Chr"+md1.getChromosome()+":"+md1.getStartPos()+".."+md1.getStopPos());
            locus++;
        }
    }

    public void setProcessedAssemblies(List processedAssemblies) {
        this.processedAssemblies = processedAssemblies;
    }

    public List getProcessedAssemblies() {
        return processedAssemblies;
    }
}
