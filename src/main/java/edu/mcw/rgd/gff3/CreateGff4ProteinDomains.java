package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.util.*;
import java.util.Map;

/**
 * @author mtutaj
 * @since Sep 23, 2019
 */
public class CreateGff4ProteinDomains {

    private RgdGff3Dao dao = new RgdGff3Dao();
    Logger log = LogManager.getLogger("domains");

    private List<Integer> processedMapKeys;
    private String outDir;
    private String trackName;

    /**
     * load the species list and assemblies from properties/AppConfigure.xml
     */
    public void run() throws Exception {

        try {
            log.info("PROTEIN DOMAIN GFF3 GENERATOR -- all species");
            log.info(dao.getConnectionInfo());
            log.info("");

            long time0 = System.currentTimeMillis();

            getProcessedMapKeys().parallelStream().forEach( mapKey -> {
                try {
                    String assemblyDir = Manager.getInstance().getAssemblies().get(mapKey);
                    if( assemblyDir==null ) {
                        return;
                    }
                    assemblyDir += "/" + Gff3Utils.getAssemblyDirStandardized(mapKey);

                    CreateInfo info = new CreateInfo();

                    int speciesTypeKey = MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey();

                    info.setMapKey( mapKey );
                    info.setToDir( assemblyDir + "/" + getOutDir() );
                    info.setSpeciesTypeKey( speciesTypeKey );
                    info.setCompressMode( Gff3ColumnWriter.COMPRESS_MODE_BGZIP );

                    run(info);
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            });

            log.info("");
            log.info("OK  elapsed " + Utils.formatElapsedTime(time0, System.currentTimeMillis()));
            log.info("");

        } catch( Exception e ) {
            Utils.printStackTrace(e, log);
            throw e;
        }
    }

    /**
     * generate gff3 file for protein domains
     */
    public void run(CreateInfo info) throws Exception {

        String speciesName = SpeciesType.getCommonName(info.getSpeciesTypeKey());

        String ucscId = Gff3Utils.getAssemblySymbol(info.getMapKey());
        String refseqId = MapManager.getInstance().getMap(info.getMapKey()).getRefSeqAssemblyName();
        String fileName = info.getToDir() + "/" + speciesName + " " + refseqId+" ("+ucscId+") ";

        String gffFile = fileName + getTrackName() + ".gff3";
        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gffFile, info.getCompressMode());

        SequenceRegionWatcher sequenceRegionWatcher = new SequenceRegionWatcher(info.getMapKey(), gff3Writer, dao);

        List<String> chromosomes = getChromosomes(info.getMapKey());

        Map<Integer, String> domainNamesMap = dao.getProteinDomainNames();

        int dataLinesWritten = 0;
        for( String chr: chromosomes ) {

            List<MapData> mds = dao.getMapDataByMapKeyChr(chr, info.getMapKey(), RgdId.OBJECT_KEY_PROTEIN_DOMAINS);

            log.debug("  "+ucscId+": data lines written for chr "+chr+":  "+mds.size());

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
                    attributesHashMap.put("proteins", proteins);
                }

                // lazy create of gff3 writer
                if( gff3Writer==null ) {

                    gff3Writer.print("# RAT GENOME DATABASE (https://rgd.mcw.edu/)\n");
                    gff3Writer.print("# Species: "+ speciesName+"\n");
                    gff3Writer.print("# Assembly: "+ refseqId+"\n");
                    gff3Writer.print("# Primary Contact: mtutaj@mcw.edu\n");
                    gff3Writer.print("# Generated: "+new Date()+"\n");
                }

                gff3Writer.writeFirst8Columns(md.getChromosome(), "RGD", "sequence feature", md.getStartPos(), md.getStopPos(), ".", md.getStrand(), ".");
                gff3Writer.writeAttributes4Gff3(attributesHashMap);
                dataLinesWritten++;

                sequenceRegionWatcher.emit(chr);
            }
        }

        gff3Writer.close();
        if( dataLinesWritten>0 ) {
            gff3Writer.sortInMemory();
        } else {
            // no data in file: don't generate the file
            File f = new File(gff3Writer.getOutFileName());
            f.deleteOnExit();
        }

        synchronized( this.getClass() ) {
            log.info(speciesName+", MAP_KEY="+info.getMapKey()+" ("+ ucscId+")   -- data lines: "+Utils.formatThousands(dataLinesWritten));
        }
    }

    List<String> getChromosomes(int mapKey) throws Exception {

        // strip version numbers from scaffold accessions
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

    public String getTrackName() {
        return trackName;
    }

    public void setTrackName(String trackName) {
        this.trackName = trackName;
    }
}
