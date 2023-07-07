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
    private List<Integer> processedMapKeys;
    private String outDir;
    private String congenicStrains;
    private String mutantStrains;

    final String source = "RGD_Rat_Strain";
    
    /**
     * load the species list and assemblies from properties/AppConfigure.xml
     */
    public void run() throws Exception {

        getProcessedMapKeys().parallelStream().forEach( mapKey -> {

            String assemblyDir = Manager.getInstance().getAssemblies().get(mapKey);
            if( assemblyDir==null ) {
                return;
            }

            try {
                CreateInfo info = new CreateInfo();
                int speciesTypeKey = MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey();

                info.setMapKey(mapKey);
                info.setToDir(assemblyDir + "/" + getOutDir());
                info.setSpeciesTypeKey(speciesTypeKey);
                info.setCompressMode(Gff3ColumnWriter.COMPRESS_MODE_BGZIP);

                run4Strains(info, false);
                run4Strains(info, true);

            } catch( Exception e ) {
                throw new RuntimeException(e);
            }
        });
    }

    public void run4Strains(CreateInfo info, boolean processCongenicStrains) throws Exception{
        CounterPool counters = new CounterPool();

        String speciesName = SpeciesType.getCommonName(info.getSpeciesTypeKey());

        String ucscId = Gff3Utils.getAssemblySymbol(info.getMapKey());
        String refseqId = MapManager.getInstance().getMap(info.getMapKey()).getRefSeqAssemblyName();
        String fileName = info.getToDir() + "/" + speciesName + " " + refseqId+" ("+ucscId+") ";

        System.out.println("START "+(processCongenicStrains?"Congenic":"Mutant") + " Strain GFF3 Generator for  "+speciesName
                +"  MAP_KEY="+info.getMapKey()+"  ASSEMBLY "+ MapManager.getInstance().getMap(info.getMapKey()).getName());
        System.out.println("========================");

        String gffFile = fileName + (processCongenicStrains ? getCongenicStrains() : getMutantStrains()) + ".gff3";

        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gffFile, false, info.getCompressMode());

        String header = "# RAT GENOME DATABASE (https://rgd.mcw.edu/)\n";
        header += "# Species: "+ speciesName+"\n";
        header += "# Assembly: "+ refseqId+"\n";
        header += "# Primary Contact: mtutaj@mcw.edu\n";
        header += "# Generated: "+new Date()+"\n";

        gff3Writer.print(header);


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

            if( processCongenicStrains && !isCongenic ) {
                continue;
            }
            if( !processCongenicStrains && !isMutant ) {
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
                gff3Writer.writeFirst8Columns(md.getChromosome(), source, strainTypeLc,
                        md.getStartPos(),md.getStopPos(),".",strand,".");

                //create attributes hash map for this line
                Map<String, String> attributesHashMap = new HashMap<>();

                emitAttributesForMultiCongenics(attributesHashMap, newmdList, md);

                attributesHashMap.put("ID",strain.getRgdId()+"_"+md.getChromosome()+"_"+m);
                attributesHashMap.put("Name", strain.getSymbol());
                attributesHashMap.put("Alias", "RGD:"+strain.getRgdId()+","+strain.getRgdId()+","+
                        "Strain:"+parentStr+","+subStr+","+strain.getSymbol());
                attributesHashMap.put("Note", note);
                attributesHashMap.put("origin", src);
                attributesHashMap.put("Index","1");

                for( Map.Entry<String,String> entry: annotAttributes.entrySet() ) {
                    attributesHashMap.put(entry.getKey(), entry.getValue());
                }

                //write attributes into gff3 file
                gff3Writer.writeAttributes4Gff3(attributesHashMap);
            }
        }

        gff3Writer.close();
        gff3Writer.sortInMemory();

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
        System.out.println("\nOK!");
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

    public List<Integer> getProcessedMapKeys() {
        return processedMapKeys;
    }

    public void setProcessedMapKeys(List<Integer> processedMapKeys) {
        this.processedMapKeys = processedMapKeys;
    }

    public String getOutDir() {
        return outDir;
    }

    public void setOutDir(String outDir) {
        this.outDir = outDir;
    }

    public String getCongenicStrains() {
        return congenicStrains;
    }

    public void setCongenicStrains(String congenicStrains) {
        this.congenicStrains = congenicStrains;
    }

    public String getMutantStrains() {
        return mutantStrains;
    }

    public void setMutantStrains(String mutantStrains) {
        this.mutantStrains = mutantStrains;
    }
}
