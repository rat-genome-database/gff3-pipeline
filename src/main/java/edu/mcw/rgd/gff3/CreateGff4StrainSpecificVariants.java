package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.dao.impl.SampleDAO;
import edu.mcw.rgd.datamodel.Sample;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.sql.Connection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Strain-specific variants GFF3 generator for jbrowse2.
 * Patient ids are configured in AppConfigure.xml; for each patient, every sample is processed
 * in parallel. Output paths are derived from {@code sample.getMapKey()} so each assembly's
 * files land under {@code .../Rat/<assemblyDir>/Variants/{Strain Specific Variants,Damaging Variants}/}.
 */
public class CreateGff4StrainSpecificVariants {

    Logger log = LogManager.getLogger("status");

    private List<Integer> patientIds;
    private Map<Integer, String> assembliesVariants;
    private String fileSource = "CN";
    private int parallelSamples = 10;

    /** existing carpe-novo helper; reused for the actual GFF3 generation per sample */
    private final CreateGff4CarpeNovo helper = new CreateGff4CarpeNovo();

    public void run() throws Exception {
        long time0 = System.currentTimeMillis();
        helper.setFileSource(fileSource);

        log.info("=== Strain Specific Variants GFF3 generation ===");
        log.info("   patient ids: " + patientIds);
        log.info("   parallel samples: " + parallelSamples);

        for (int patientId : patientIds) {
            processPatient(patientId);
        }
        log.info("=== DONE === elapsed: " + Utils.formatElapsedTime(time0, System.currentTimeMillis()));
    }

    void processPatient(int patientId) throws Exception {
        long time0 = System.currentTimeMillis();

        SampleDAO sampleDAO = new SampleDAO();
        sampleDAO.setDataSource(DataSourceFactory.getInstance().getCarpeNovoDataSource());
        List<Sample> samples = sampleDAO.getSamples(patientId);
        log.info("PATIENT " + patientId + ": " + samples.size() + " samples");

        AtomicInteger remaining = new AtomicInteger(samples.size());
        ForkJoinPool pool = new ForkJoinPool(parallelSamples);
        try {
            pool.submit(() -> samples.parallelStream().forEach(sample -> {
                try {
                    processSample(sample);
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
                int left = remaining.decrementAndGet();
                log.info("PATIENT " + patientId + ": " + left + "/" + samples.size() + " samples remaining");
            })).get();
        } finally {
            pool.shutdown();
        }
        log.info("PATIENT " + patientId + " OK -- elapsed: " + Utils.formatElapsedTime(time0, System.currentTimeMillis()));
    }

    void processSample(Sample sample) throws Exception {
        long time0 = System.currentTimeMillis();
        int sampleId = sample.getId();
        int mapKey = sample.getMapKey();
        String sampleName = sample.getAnalysisName().replace('/', '_');

        String prefix = assembliesVariants.get(mapKey);
        if (prefix == null) {
            log.warn("SAMPLE " + sampleId + " (" + sampleName + ", mapKey=" + mapKey
                    + "): skipped -- no assembliesVariants entry for this mapKey");
            return;
        }
        String assemblyDir = Gff3Utils.getAssemblyDirStandardized(mapKey);

        String ssvDir = prefix + "/" + assemblyDir + "/Variants/Strain Specific Variants";
        String dmgDir = prefix + "/" + assemblyDir + "/Variants/Damaging Variants";
        new File(ssvDir).mkdirs();
        new File(dmgDir).mkdirs();

        String gffPath = ssvDir + "/" + sampleName + ".gff3";
        String dmgPath = dmgDir + "/" + sampleName + "_damaging.gff3";

        Gff3ColumnWriter gffWriter = new Gff3ColumnWriter(gffPath, Gff3ColumnWriter.COMPRESS_MODE_BGZIP);
        Gff3ColumnWriter dmgWriter = new Gff3ColumnWriter(dmgPath, Gff3ColumnWriter.COMPRESS_MODE_BGZIP);
        try (gffWriter; dmgWriter) {

            String assemblyName = MapManager.getInstance().getMap(mapKey).getRefSeqAssemblyName();
            String header = "#RGD SAMPLE ID: " + sampleId + "\n#ASSEMBLY: " + assemblyName + "\n";
            gffWriter.print(header);
            dmgWriter.print(header);

            RgdGff3Dao dao = new RgdGff3Dao();
            SequenceRegionWatcher srw1 = new SequenceRegionWatcher(mapKey, gffWriter, dao);
            SequenceRegionWatcher srw2 = new SequenceRegionWatcher(mapKey, dmgWriter, dao);
            String[] chromosomes = {"1","2","3","4","5","6","7","8","9","10",
                    "11","12","13","14","15","16","17","18","19","20","X","Y","MT"};
            for (String c : chromosomes) {
                srw1.emit(c);
                srw2.emit(c);
            }

            try (Connection conn = DataSourceFactory.getInstance().getCarpeNovoDataSource().getConnection()) {
                helper.process(sampleId, gffWriter, dmgWriter, conn);
            }
        }
        gffWriter.sortInMemory();
        dmgWriter.sortInMemory();

        log.info("SAMPLE " + sampleId + " (" + sampleName + ", mapKey=" + mapKey + ") OK -- elapsed: "
                + Utils.formatElapsedTime(time0, System.currentTimeMillis()));
    }

    public void setPatientIds(List<Integer> patientIds) { this.patientIds = patientIds; }
    public List<Integer> getPatientIds() { return patientIds; }

    public void setAssembliesVariants(Map<Integer, String> assembliesVariants) { this.assembliesVariants = assembliesVariants; }
    public Map<Integer, String> getAssembliesVariants() { return assembliesVariants; }

    public void setFileSource(String fileSource) { this.fileSource = fileSource; }
    public String getFileSource() { return fileSource; }

    public void setParallelSamples(int parallelSamples) { this.parallelSamples = parallelSamples; }
    public int getParallelSamples() { return parallelSamples; }
}
