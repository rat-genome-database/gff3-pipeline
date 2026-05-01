package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;
import java.util.Map;


/**
 * Created by llamers on 3/5/2020.
 */
public class CreateGff4Eva {

    private RgdGff3Dao dao = new RgdGff3Dao();
    Logger log = LogManager.getLogger("evas");
    private List<String> processedAssemblies;
    private String evaRelease;
    private String jbrowse2OutDir;

    private String getMemoryUsage() {
        Runtime rt = Runtime.getRuntime();
        long usedMB = (rt.totalMemory() - rt.freeMemory()) / (1024 * 1024);
        long totalMB = rt.totalMemory() / (1024 * 1024);
        long maxMB = rt.maxMemory() / (1024 * 1024);
        return "Memory: used=" + usedMB + "MB, allocated=" + totalMB + "MB, max=" + maxMB + "MB";
    }

    public void run() throws Exception{
        try {
            log.info("EVA GFF3 GENERATOR -- all species");
            log.info(dao.getConnectionInfo()+"\n");
            log.info(getMemoryUsage());

            long timeStart = System.currentTimeMillis();

            getProcessedAssemblies().parallelStream().forEach(assemblyInfo ->{
                try{
                    CreateInfo info = new CreateInfo();
                    info.parseFromString(assemblyInfo);
                    run(info);
                }
                catch(Exception e){
                    throw new RuntimeException(e);
                }
            });

            log.info("");
            log.info("OK elapsed "+Utils.formatElapsedTime(timeStart, System.currentTimeMillis()));
            log.info(getMemoryUsage());
            log.info("");
        }
        catch (Exception e) {
            Utils.printStackTrace(e, log);
            throw e;
        }
    }

    public void run(CreateInfo info) throws Exception {

        int mapKey = info.getMapKey();
        String fileName = info.getCanonicalFileName()+" "+getEvaRelease()+".gff3";
        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(fileName, info.getCompressMode());

        List<String> chromosomes = getChromosomes(mapKey);
        SequenceRegionWatcher sequenceRegionWatcher = new SequenceRegionWatcher(mapKey, gff3Writer, dao);
        int[] dataLinesWritten = {0};

        for(String chr : chromosomes) {
            sequenceRegionWatcher.emit(chr);
            log.debug("  MAP_KEY="+mapKey+" chr "+chr+" start -- "+getMemoryUsage());

            List<Eva> currentGroup = new ArrayList<>();

            dao.streamEvaObjectsByKeyAndChrom(mapKey, chr, rs -> {
                Eva eva = mapResultSetToEva(rs);

                // if rsId changed, flush the previous group
                if(!currentGroup.isEmpty() && !currentGroup.get(0).getRsId().equals(eva.getRsId())) {
                    try {
                        dataLinesWritten[0] += writeGroup(currentGroup, gff3Writer);
                    } catch(Exception e) {
                        throw new RuntimeException(e);
                    }
                    currentGroup.clear();
                }

                currentGroup.add(eva);
            });

            // flush remaining group for this chromosome
            if(!currentGroup.isEmpty()) {
                dataLinesWritten[0] += writeGroup(currentGroup, gff3Writer);
                currentGroup.clear();
            }
            log.debug("  MAP_KEY="+mapKey+" chr "+chr+" done -- "+getMemoryUsage());
        }

        gff3Writer.close();
        log.info("MAP_KEY="+mapKey+" before sortInMemory -- "+getMemoryUsage());
        gff3Writer.sortInMemory();
        log.info("MAP_KEY="+mapKey+" after sortInMemory -- "+getMemoryUsage());

        synchronized( this.getClass() ) {
            log.info(info.speciesName+", MAP_KEY="+info.getMapKey()+" ("+ info.refseqId+")   -- data lines: "+Utils.formatThousands(dataLinesWritten[0]));
        }

        copyToJBrowse2OutDir( gff3Writer.getOutFileName(), info );
    }

    private Eva mapResultSetToEva(ResultSet rs) throws SQLException {
        Eva eva = new Eva();
        eva.setEvaid(rs.getInt("EVA_ID"));
        eva.setChromosome(rs.getString("CHROMOSOME"));
        eva.setPos(rs.getInt("POS"));
        eva.setRsid(rs.getString("RS_ID"));
        eva.setRefnuc(rs.getString("REF_NUC"));
        eva.setVarnuc(rs.getString("VAR_NUC"));
        eva.setPadBase(rs.getString("PADDING_BASE"));
        eva.setSoterm(rs.getString("SO_TERM_ACC"));
        eva.setMapkey(rs.getInt("MAP_KEY"));
        return eva;
    }

    private String mapSoTerm(String evaSoTerm) {
        return switch (evaSoTerm) {
            case "SO:0002007", "0002007" -> "MNP";
            case "SO:0000159", "0000159" -> "DELETION";
            case "SO:0000667", "0000667" -> "INSERTION";
            case "SO:1000032", "1000032" -> "DELIN";
            case "SO:0000705", "0000705" -> "TANDEM_REPEAT";
            default -> "SNP";
        };
    }

    private int writeGroup(List<Eva> group, Gff3ColumnWriter gff3Writer) throws Exception {
        if(group.isEmpty()) return 0;

        String rsId = group.get(0).getRsId();

        // Sub-group by (chr, pos, refNuc, soTerm). Each sub-group represents one distinct variant
        // sharing the rsId. Within a sub-group, rows may carry different varNuc values — those are
        // genuine multi-allelic alternatives at the same reference span and collapse into a single
        // GFF3 line with allele = refNuc/alt1/alt2/...
        // A LinkedHashMap keeps insertion order stable; we sort keys explicitly below for determinism.
        Map<String, SubGroup> subGroups = new LinkedHashMap<>();
        for(Eva eva : group) {
            String soTerm = mapSoTerm(eva.getSoTerm());
            String refNuc = Utils.NVL(eva.getRefNuc(), "-");
            String varNuc = Utils.NVL(eva.getVarNuc(), "-");
            String key = eva.getChromosome()+"\u0001"+eva.getPos()+"\u0001"+refNuc+"\u0001"+soTerm;
            SubGroup sg = subGroups.computeIfAbsent(key, k -> new SubGroup(eva, soTerm, refNuc));
            if(!sg.varNucs.contains(varNuc)) {
                sg.varNucs.add(varNuc);
            }
        }

        // Deterministic ordering of sub-groups: by pos, then refNuc, then soTerm
        List<SubGroup> ordered = new ArrayList<>(subGroups.values());
        ordered.sort(Comparator
                .comparingInt((SubGroup sg) -> sg.representative.getPos())
                .thenComparing(sg -> sg.refNuc)
                .thenComparing(sg -> sg.soTerm));

        int suffix = 0;
        for(SubGroup sg : ordered) {
            Eva rep = sg.representative;
            int offset = rep.getRefNuc() == null ? 0 : rep.getRefNuc().length() - 1;

            gff3Writer.writeFirst8Columns(rep.getChromosome(), "EVA", sg.soTerm, rep.getPos(), rep.getPos() + offset, ".", ".", ".");

            HashMap<String, String> attributes = new HashMap<>();
            attributes.put("ID", suffix == 0 ? rsId : rsId + "_" + suffix);
            attributes.put("Alias", rsId);

            Collections.sort(sg.varNucs);
            StringBuilder allele = new StringBuilder(sg.refNuc);
            for(String v : sg.varNucs) {
                allele.append("/").append(v);
            }
            attributes.put("allele", allele.toString());

            gff3Writer.writeAttributes4Gff3(attributes);
            suffix++;
        }
        return ordered.size();
    }

    private static class SubGroup {
        final Eva representative;
        final String soTerm;
        final String refNuc;
        final List<String> varNucs = new ArrayList<>();
        SubGroup(Eva representative, String soTerm, String refNuc) {
            this.representative = representative;
            this.soTerm = soTerm;
            this.refNuc = refNuc;
        }
    }

    void copyToJBrowse2OutDir( String fullFileName, CreateInfo info ) throws IOException {

        int sepPos = fullFileName.lastIndexOf('/');
        if( sepPos<0 ) {
            sepPos = fullFileName.lastIndexOf('\\');
        }
        String fileName = fullFileName.substring(sepPos+1);

        String outFileDir = getJbrowse2OutDir()
                .replace("{SPECIES}", info.speciesName)
                .replace("{ASSEMBLY}", info.refseqId);
        File outDir = new File(outFileDir);
        if( !outDir.exists() ) {
            outDir.mkdirs();
        }
        String outFileName = outFileDir + "/"+fileName;

        File inFile = new File(fullFileName);
        File outFile = new File(outFileName);
        log.info("copying "+inFile.getAbsolutePath()+" to "+outFile.getAbsolutePath());
        Files.copy( inFile.toPath(), outFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
    }

    public void setProcessedAssemblies(List<String> processedAssemblies) {
        this.processedAssemblies = processedAssemblies;
    }

    public List<String> getProcessedAssemblies() {
        return processedAssemblies;
    }

    List<String> getChromosomes(int mapKey) throws Exception {

        // truncate version numbers from scaffold accessions
        List<String> result = new ArrayList<>();
        List<Chromosome> chromosomes = dao.getChromosomes(mapKey);
        for( Chromosome chr: chromosomes ) {
            String c = chr.getChromosome();
            result.add(c);
        }
        return result;
    }

    public String getEvaRelease() {
        return evaRelease;
    }

    public void setEvaRelease(String evaRelease) {
        this.evaRelease = evaRelease;
    }

    public String getJbrowse2OutDir() {
        return jbrowse2OutDir;
    }

    public void setJbrowse2OutDir(String jbrowse2OutDir) {
        this.jbrowse2OutDir = jbrowse2OutDir;
    }
}
