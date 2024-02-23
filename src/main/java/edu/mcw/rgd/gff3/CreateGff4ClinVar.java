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

    private Map<String, Map<String,String>> tracks;
    private RgdGff3Dao dao = new RgdGff3Dao();

    /**
     * generate gff3 file for ClinVar  Variants/ClinVar/GRCh37.p13
     *                                 Variants/ClinVar/GRCh38.p14
     */
    public void run(int mapKey, String outDir, int compressMode) throws Exception {

        String gffDir = outDir;
        System.out.println("  "+outDir);

        Gff3ColumnWriter pathogenicWriter = null;
        Gff3ColumnWriter snpsIndelsWriter = null;

        // create all gff3 writers
        Map<String, Gff3ColumnWriter> gffWriters = new HashMap<>();
        for( String trackName: tracks.keySet() ) {
            String gffName = gffDir + "/" + trackName + ".gff3";
            Gff3ColumnWriter gffWriter = new Gff3ColumnWriter(gffName, compressMode);
            gffWriters.put(trackName, gffWriter);

            if( trackName.equals("ClinVar - SNPs and Indels") ) {
                snpsIndelsWriter = gffWriter;
            }
            if( trackName.equals("ClinVar - SNPs and Indels (Pathogenic)") ) {
                pathogenicWriter = gffWriter;
            }
        }

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
        System.out.println("  variants retrieved: "+variants.size());

        for( VariantInfo var: variants ) {

            String trackName = null;
            String soType = null;
            Gff3ColumnWriter gff3Writer = null;
            for( Map.Entry<String, Map<String,String>> entry: tracks.entrySet() ) {
                Map<String,String> typeMap = entry.getValue();
                soType = typeMap.get(var.getObjectType());
                if( soType!=null ) {
                    trackName = entry.getKey();
                    gff3Writer = gffWriters.get(trackName);
                    break;
                }
            }

            if( soType==null ) {
                variantsBadType++;
                continue;
            }

            // these variants must have a position on given assembly
            List<MapData> mds = dao.getMapData(var.getRgdId(), mapKey);
            for( MapData md: mds ) {
                gff3Writer.writeFirst8Columns(md.getChromosome(), "ClinVar", soType, md.getStartPos(), md.getStopPos(), ".", ".", ".");
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

                if( gff3Writer == snpsIndelsWriter && var.getClinicalSignificance().contains("pathogenic") ) {

                    pathogenicWriter.writeFirst8Columns(md.getChromosome(), "ClinVar", soType, md.getStartPos(), md.getStopPos(), ".", ".", ".");
                    pathogenicWriter.writeAttributes4Gff3(attributes);
                }
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


        // close all writers
        for( Gff3ColumnWriter writer: gffWriters.values() ) {
            writer.close();
            writer.sortInMemory();
        }


        System.out.println("ClinVar variants processed: " + variants.size());
        System.out.println(" variants skipped: unsupported type: " + variantsBadType);
        System.out.println(" variants with a map position: " + variantsWithPos);
        System.out.println(" variants with nO map position: " + variantsWithoutPos);
        System.out.println(" variants with more than one map position: " + variantsWithMultiPos);

        System.out.println(" variants with clinical significance :" + variantsWithClinicalSignificance);
        System.out.println(" variants with method type           :" + variantsWithMethodType);
        System.out.println(" variants with molecular consequence :" + variantsWithMolecularConsequence);
        System.out.println(" variants with age of onset          :" + variantsWithAgeOfOnset);
        System.out.println(" variants with prevalence            :" + variantsWithPrevalence);
        System.out.println(" variants with submitter             :" + variantsWithSubmitter);
        System.out.println(" variants with trait name            :" + variantsWithTraitName);

        System.out.print("\nGFF3 File SUCCESSFUL!");
    }

    public Map<String, Map<String,String>> getTracks() {
        return tracks;
    }

    public void setTracks( Map<String, Map<String,String>> tracks ) {
        this.tracks = tracks;
    }
}

