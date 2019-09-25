package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author mtutaj
 * @since Sep 23, 2019
 */
public class CreateGff4ProteinDomains {

    private int mapKey;
    private int speciesTypeKey;
    private String toFile;
    private boolean compress;
    private RgdGff3Dao dao = new RgdGff3Dao();
    Logger log = Logger.getLogger("domains");

    /**
     * generate gff3 file for protein domains
     * @param compress
     */
    public void run(boolean compress) throws Exception {
        this.compress = compress;

        long time0 = System.currentTimeMillis();

        String species = SpeciesType.getCommonName(speciesTypeKey);

        String gffFile = getToFile();
        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gffFile, false, compress);

        List<String> chromosomes = getChromosomes();

        Map<Integer, String> domainNamesMap = dao.getProteinDomainNames();

        log.info("Generate GFF3 file for "+species+", MAP_KEY="+mapKey+" ("+ MapManager.getInstance().getMap(mapKey).getName()+")");
        log.info("    "+dao.getConnectionInfo());

        int dataLinesWritten = 0;
        for( String chr: chromosomes ) {

            List<MapData> mds = dao.getMapDataByMapKeyChr(chr, mapKey, RgdId.OBJECT_KEY_PROTEIN_DOMAINS);

            log.debug("  data lines written for chr "+chr+":  "+mds.size());

            for (MapData md : mds) {
                Map<String, String> attributesHashMap = new HashMap<>();

                String domainName = domainNamesMap.get(md.getRgdId());
                attributesHashMap.put("ID", "RGD:" + md.getRgdId());
                attributesHashMap.put("Name", domainName);
                attributesHashMap.put("Alias", domainName + "," + md.getRgdId());
                if (md.getNotes() != null) {
                    // in notes there is information about proteins, like this:
                    //   Q924C6 part 1; E9Q600 part 1; Q8CEU1 part 1
                    // convert it to
                    //   Q924C6, E9Q600, Q8CEU1
                    String[] words = md.getNotes().split("[\\;] ");
                    String proteins = "";
                    for( String word: words ) {
                        int spacePos = word.indexOf(' ');
                        String protein;
                        if( spacePos>0 ) {
                            protein = word.substring(0, spacePos);
                        } else {
                            protein = word;
                        }
                        if( proteins.isEmpty() ) {
                            proteins = protein;
                        } else {
                            proteins += ", "+protein;
                        }
                    }
                    attributesHashMap.put("Note", proteins);
                }

                gff3Writer.writeFirst8Columns(md.getChromosome(), "RGD", "sequence feature", md.getStartPos(), md.getStopPos(), ".", md.getStrand(), ".");
                gff3Writer.writeAttributes4Gff3(attributesHashMap);
                dataLinesWritten++;
            }
        }
        gff3Writer.close();

        log.info("  data lines written:  " + dataLinesWritten);
        log.info("OK!  elapsed " + Utils.formatElapsedTime(time0, System.currentTimeMillis()));
        log.info("========");
    }

    List<String> getChromosomes() throws Exception {

        // truncate version numbers from scaffold accessions
        List<String> result = new ArrayList<>();
        List<Chromosome> chromosomes = dao.getChromosomes(mapKey);
        for( Chromosome chr: chromosomes ) {
            String c = chr.getChromosome();
            if( c.startsWith("NW_") ) {
                int dotPos = c.indexOf(".");
                if( dotPos>0 ) {
                    result.add(c.substring(0, dotPos));
                    continue;
                }
            }
            result.add(c);
        }
        return result;
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

    public boolean isCompress() {
        return compress;
    }
}
