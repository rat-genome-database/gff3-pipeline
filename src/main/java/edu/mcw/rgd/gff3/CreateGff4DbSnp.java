package edu.mcw.rgd.gff3;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author mtutaj
 * @since July 21, 2017
 */
public class CreateGff4DbSnp {

    private int mapKey;
    private int speciesTypeKey;
    private String toFile;
    private String build;
    private int compressMode;
    private Map<String,String> allowedVarTypes;
    private RgdGff3Dao dao = new RgdGff3Dao();

    public CreateGff4DbSnp() {

        allowedVarTypes = new HashMap<>();
        allowedVarTypes.put("snp","SNV");
    }

    /**
     * generate gff3 file for DB_SNP
     * @param compressMode
     */
    public void run(int compressMode) throws Exception {
        this.compressMode = compressMode;

        String gffFile = getToFile();
        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gffFile, compressMode);

        int snpsBadType = 0;

        int snpsWithClinicalSignificance = 0;
        int snpsWithFunctionClass = 0;

        List<RgdGff3Dao.DbSnp> dbSnps = dao.getDbSnps(getMapKey(), getBuild());
        for( RgdGff3Dao.DbSnp dbSnp: dbSnps ) {

            // only small variants are exported
            String snpType = allowedVarTypes.get(dbSnp.snpClass);
            if( snpType==null ) {
                snpsBadType++;
                continue;
            }

            gff3Writer.writeFirst8Columns(dbSnp.chr, "dbSNP", snpType, dbSnp.pos, dbSnp.pos, ".", ".", ".");
            HashMap<String, String> attributes = new HashMap<>();

            attributes.put("ID", dbSnp.snpName);
            attributes.put("Name", dbSnp.snpName);
            attributes.put("Symbol", dbSnp.snpName);
            attributes.put("allele", dbSnp.refAllele+">"+dbSnp.allele);

            if( dbSnp.clinicalSignificance!=null ) {
                attributes.put("clinicalSignificance", dbSnp.clinicalSignificance.toLowerCase());
                snpsWithClinicalSignificance++;
            }
            if( dbSnp.functionClass!=null ) {
                attributes.put("functionClass", dbSnp.functionClass);
                snpsWithFunctionClass++;
            }
            gff3Writer.writeAttributes4Gff3(attributes);
        }

        gff3Writer.close();


        System.out.println("DbSnp objects processed:" + dbSnps.size());
        System.out.println(" dbSnps skipped: type not in(snv)  :" + snpsBadType);
        System.out.println(" dbSnps with clinical significance :" + snpsWithClinicalSignificance);
        System.out.println(" dbSnps with functional class      :" + snpsWithFunctionClass);

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

    public String getBuild() {
        return build;
    }

    public void setBuild(String build) {
        this.build = build;
    }
}
