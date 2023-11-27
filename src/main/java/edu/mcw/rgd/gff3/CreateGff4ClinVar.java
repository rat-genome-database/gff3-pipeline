package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.VariantInfo;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by mtutaj on 11/17/2016.
 */
public class CreateGff4ClinVar {

    private int mapKey;
    private int speciesTypeKey;
    private String toFile;
    private int compressMode;
    private Map<String,String> allowedVarTypes;
    private RgdGff3Dao dao = new RgdGff3Dao();

    public CreateGff4ClinVar() {

        allowedVarTypes = new HashMap<>();
        allowedVarTypes.put("single nucleotide variant","SNV");
        allowedVarTypes.put("deletion","deletion");
        allowedVarTypes.put("duplication","duplication");
        allowedVarTypes.put("insertion","insertion");
        allowedVarTypes.put("indel","indel");
        allowedVarTypes.put("variation","sequence alteration");
    }

    /**
     * generate gff3 file for ClinVar
     */
    public void run(int compressMode) throws Exception {
        this.compressMode = compressMode;

        String gffFile = getToFile();
        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gffFile, compressMode);

        int variantsBadType = 0;
        int variantsWithPos = 0;
        int variantsWithoutPos = 0;
        int variantsWithMultiPos = 0;

        int variantsWithClinicalSignificance = 0;
        int variantsWithMethodType = 0;
        int variantsWithMolecularConsequence = 0;
        int variantsWithAgeOfOnset = 0;
        int variantsWithPrevalence = 0;
        int variantsWithSubmitter = 0;
        int variantsWithTraitName = 0;

        List<VariantInfo> variants = dao.getClinVarVariants();
        for( VariantInfo var: variants ) {

            // only small variants are exported
            String varType = allowedVarTypes.get(var.getObjectType());
            if( varType==null ) {
                variantsBadType++;
                continue;
            }

            // these variants must have a position on given assembly
            List<MapData> mds = dao.getMapData(var.getRgdId(), getMapKey());
            for( MapData md: mds ) {
                gff3Writer.writeFirst8Columns(md.getChromosome(), "ClinVar", varType, md.getStartPos(), md.getStopPos(), ".", ".", ".");
                HashMap<String, String> attributes = new HashMap<>();

                attributes.put("ID", var.getSymbol());
                attributes.put("Name", var.getName());
                attributes.put("Dbxref", "RGD:" + var.getRgdId());
                attributes.put("Symbol", var.getSymbol());

                if( var.getClinicalSignificance()!=null ) {
                    attributes.put("clinicalSignificance", var.getClinicalSignificance());
                    variantsWithClinicalSignificance++;
                }
                if( var.getMethodType()!=null ) {
                    attributes.put("methodType", var.getMethodType());
                    variantsWithMethodType++;
                }
                if( var.getMolecularConsequence()!=null ) {
                    attributes.put("molecularConsequence", var.getMolecularConsequence());
                    variantsWithMolecularConsequence++;
                }
                if( var.getAgeOfOnset()!=null ) {
                    attributes.put("ageOfOnset", var.getAgeOfOnset());
                    variantsWithAgeOfOnset++;
                }
                if( var.getPrevalence()!=null ) {
                    attributes.put("prevalence", var.getPrevalence());
                    variantsWithPrevalence++;
                }
                if( var.getSubmitter()!=null ) {
                    attributes.put("submitter", var.getSubmitter());
                    variantsWithSubmitter++;
                }
                if( var.getTraitName()!=null ) {
                    attributes.put("traitName", var.getTraitName());
                    variantsWithTraitName++;
                }
                gff3Writer.writeAttributes4Gff3(attributes);
            }

            if( mds.isEmpty() ) {
                variantsWithoutPos++;
            }
            else if( mds.size()==1 ) {
                variantsWithPos++;
            }
            else {
                variantsWithMultiPos++;
            }
        }

        gff3Writer.close();


        System.out.println("ClinVar variants processed:" + variants.size());
        System.out.println(" variants skipped: type not in(snv,del,dup,ins):" + variantsBadType);
        System.out.println(" variants with Map position:" + variantsWithPos);
        System.out.println(" variants with NO map position:" + variantsWithoutPos);
        System.out.println(" variants with more than one map position:" + variantsWithMultiPos);

        System.out.println(" variants with clinical significance :" + variantsWithClinicalSignificance);
        System.out.println(" variants with method type           :" + variantsWithMethodType);
        System.out.println(" variants with molecular consequence :" + variantsWithMolecularConsequence);
        System.out.println(" variants with age of onset          :" + variantsWithAgeOfOnset);
        System.out.println(" variants with prevalence            :" + variantsWithPrevalence);
        System.out.println(" variants with submitter             :" + variantsWithSubmitter);
        System.out.println(" variants with trait name            :" + variantsWithTraitName);

        //System.out.println("ClinVar variants having aliases:" + sslpsAliasCount);
        //System.out.println("ClinVar variants having an associated gene RGDID:" + sslpsAssocGeneRgd);

        System.out.print("\nGFF3 File SUCCESSFUL!");
    }

    public int getMapKey() {
        return mapKey;
    }

    public void setMapKey(int mapKey) {
        this.mapKey = mapKey;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }

    public String getToFile() {
        return toFile;
    }

    public void setToFile(String toFile) {
        this.toFile = toFile;
    }

    public int getCompressMode() {
        return compressMode;
    }
}
