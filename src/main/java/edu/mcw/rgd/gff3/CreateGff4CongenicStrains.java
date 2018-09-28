package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.Strain;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.process.Utils;

import java.util.*;

/**
 * @author pjayaraman
 * Date: Mar 23, 2011
 */
public class CreateGff4CongenicStrains {
    String toFile;
    RgdGff3Dao dao = new RgdGff3Dao();
    int map_key;
    String source = "RGD_Rat_Strain";
    
    int strainsWithDiseaseAnn;
    int strainsNoDisAnn;
    int strainsWithPhenoAnn;
    int strainsNoPhenoAnn;
    int strainsNoAnn;
    int strainsWithAnn;

    public void creategff4CongenicStrains(boolean compress) throws Exception{

        //create new file to write and add first line ##gff-version-3 is mandatory..!!!!
        String gffFileCongenic = getToFile()+"RatCongenicStrains_RGD.gff3";
        String RATMINEGffFile = getToFile()+"RatStrains_RATMINE.gff3";
        String gffFileMutant = getToFile()+"RatMutantStrains_RGD.gff3";

        //initialization of the gff3writer
        Gff3ColumnWriter gff3WriterCongenic = new Gff3ColumnWriter(gffFileCongenic, false, compress);
        Gff3ColumnWriter gff3WriterRATMINE = new Gff3ColumnWriter(RATMINEGffFile, false, compress);
        Gff3ColumnWriter gff3WriterMutant = new Gff3ColumnWriter(gffFileMutant, false, compress);

        // get all active strains having positions on the given assembly
        List<Strain> strainList = dao.getMappedStrains(map_key);

        int congenicStrains=0;
        int mutantStrains=0;
        int mappedStrainCount=0;
        int noMapPos=0;
        int moreThanOneMapPosCount=0;

        strainsWithDiseaseAnn=0;
        strainsNoDisAnn=0;
        strainsWithPhenoAnn=0;
        strainsNoPhenoAnn=0;
        strainsNoAnn=0;
        strainsWithAnn=0;


        //for each strain rgdid... get the annotations from the full annot table..
        // and get the map data positions.
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
                congenicStrains++;
            if( isMutant )
                mutantStrains++;


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

            Map<String,String> annotAttributes = processAnnotations(strain.getRgdId());


            List<MapData> newmdList = dao.getMapData(strain.getRgdId(), map_key);

            if(newmdList.size()>1){
                moreThanOneMapPosCount++;
                if( isCongenic ) {
                    strainTypeLc = "multicongenic";
                }
            }else if(newmdList.size()==1){
                mappedStrainCount++;
            }else if(newmdList.size()==0){
                noMapPos++;
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

        System.out.println("\nCongenic Strain records processed: "+ congenicStrains);
        System.out.println("Mutant Strain records processed: "+ mutantStrains);
        System.out.println("Total count of Strain having only one Map position: "+ mappedStrainCount);
        System.out.println("Number of Strains with more than one Map position: "+ moreThanOneMapPosCount);
        System.out.println("Number of Strains having NO Map positions: "+ noMapPos);
        System.out.println("Number of Strains with Annotation: "+ strainsWithAnn);
        System.out.println("Number of Strains with NO Annotations at all: "+ strainsNoAnn);
        System.out.println("Number of Strain Records with Disease annotations: "+ strainsWithDiseaseAnn);
        System.out.println("Number of Strain Records with NO Disease annotations: "+ strainsNoDisAnn);
        System.out.println("Number of Strain Records with Phenotype annotations: "+ strainsWithPhenoAnn);
        System.out.println("Number of Strain Records with NO Phenotype annotations: "+ strainsNoPhenoAnn);
        System.out.println("\nGFF3 File SUCCESSFUL!");

    }

    Map<String,String> processAnnotations(int strainRgdId) throws Exception {

        List<Annotation> disAnn = new ArrayList<>();
        List<Annotation> pheAnn = new ArrayList<>();

        List<Annotation> annList = dao.getAnnotations(strainRgdId);
        if( annList.isEmpty() ){
            strainsNoAnn++;
            return Collections.emptyMap();
        }

        strainsWithAnn++;
        Map<String,String> annotAttributes = new HashMap<>();

        for(Annotation annObj:annList){
            if(annObj.getTermAcc().startsWith("RDO")){
                disAnn.add(annObj);
            }else if(annObj.getTermAcc().startsWith("MP")){
                pheAnn.add(annObj);
            }
        }

        if(disAnn.size()!=0){
            strainsWithDiseaseAnn++;

            String disOntValues="";
            for( Annotation disease: disAnn ){
                disOntValues += disease.getTermAcc()+":"+disease.getTerm()+",";
            }
            //remove last "," from string
            disOntValues = disOntValues.substring(0,disOntValues.length()-1);

            annotAttributes.put("DiseaseOntologyAssociation",disOntValues);
        }else {
            strainsNoDisAnn++;
        }

        if(pheAnn.size()!=0){
            strainsWithPhenoAnn++;

            String pheOntValues="";

            for( Annotation phenotype: pheAnn ){
                pheOntValues += phenotype.getTermAcc()+":"+phenotype.getTerm()+",";
            }
            //remove last "," from this line
            pheOntValues = pheOntValues.substring(0,pheOntValues.length()-1);

            annotAttributes.put("PhenotypeOntologyAssociation",pheOntValues);
        }else {
            strainsNoPhenoAnn++;
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

    public String getToFile() {
        return toFile;
    }

    public void setToFile(String toFile) {
        this.toFile = toFile;
    }
    public int getMap_key() {
        return map_key;
    }

    public void setMap_key(int map_key) {
        this.map_key = map_key;
    }
}
