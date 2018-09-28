package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.dao.impl.PhenominerDAO;
import edu.mcw.rgd.datamodel.pheno.Condition;
import edu.mcw.rgd.datamodel.pheno.IndividualRecord;
import edu.mcw.rgd.datamodel.pheno.Record;
import edu.mcw.rgd.datamodel.pheno.Study;
import edu.mcw.rgd.process.Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.List;

/**
 * Created by mtutaj on 5/18/2017.
 */
public class PhenoLoader {

    public static void main(String[] args) throws Exception {

        PhenominerDAO dao = new PhenominerDAO();

        List<Study> studies = dao.getStudies();
        for( Study  study: studies ) {
            int cnt = dao.getRecordCountForStudy(study.getId());
            System.out.println("STUDY_ID="+study.getId()+", REC_CNT="+cnt);
            List<Integer> cnts = dao.getRecordStatusCountForStudy(study.getId());
            for( int c: cnts ) {
                System.out.println("  STUDY_ID="+study.getId()+", STATUS_CNT="+c);
            }
        }

        String fname = "/tmp/experiment2.txt"; // CMO:0000257
        BufferedReader reader = new BufferedReader(new FileReader(fname));
        String line = reader.readLine(); // skip header
        while( (line=reader.readLine())!=null ) {
            String[] cols = line.split("[\\t]", -1);
            String animalNr = cols[0]; // f.e. 1, 2, etc
            String indAnimalId = cols[1]; // f.e. 6060_B
            String sampleNotes = cols[2];
            String sampleAgeDaysFromDobLowBound = cols[3];
            String sampleAgeDaysFromDobHighBound = cols[4];
            String sampleSex = cols[5];
            String sampleStrainOntId = cols[6];
            String measurementValue = cols[7];
            String measurementUnits = cols[8];
            String indMeasurementValue = cols[9];

            Record record = new Record();
            record.setExperimentId(21064); // exp1
            record.setStudyId(6060);

            if( !Utils.isStringEmpty(measurementValue) ) {
                record.setMeasurementValue(measurementValue);
            }
            if( !Utils.isStringEmpty(measurementUnits) ) {
                record.setMeasurementUnits(measurementUnits);
            }

            if( !Utils.isStringEmpty(sampleNotes) ) {
                record.getSample().setNotes(sampleNotes);
            }
            if( !Utils.isStringEmpty(sampleAgeDaysFromDobLowBound)) {
                record.getSample().setAgeDaysFromLowBound(Integer.parseInt(sampleAgeDaysFromDobLowBound));
            }
            if( !Utils.isStringEmpty(sampleAgeDaysFromDobHighBound)) {
                record.getSample().setAgeDaysFromHighBound(Integer.parseInt(sampleAgeDaysFromDobHighBound));
            }
            if( !Utils.isStringEmpty(sampleSex) ) {
                record.getSample().setSex(sampleSex);
            }
            if( !Utils.isStringEmpty(sampleStrainOntId) ) {
                record.getSample().setStrainAccId(sampleStrainOntId);
            }
            record.getSample().setNumberOfAnimals(1);

            record.setHasIndividualRecord(true);

            // placeholder for measurement method
            record.getMeasurementMethod().setAccId("MMO:0000000");

            // placeholder for clinical measurement
            record.getClinicalMeasurement().setAccId("CMO:0000625");

            // placeholder for condition
            Condition cond = new Condition();
            cond.setOntologyId("XCO:0000000");
            cond.setOrdinality(1);
            record.getConditions().add(cond);

            dao.insertRecord(record);


            IndividualRecord irec = new IndividualRecord();
            if( !Utils.isStringEmpty(indMeasurementValue) ) {
                irec.setMeasurementValue(indMeasurementValue);
            }
            irec.setAnimalId(indAnimalId);
            irec.setRecordId(record.getId());
            dao.insertIndividualRecord(irec);
        }
        reader.close();
    }
}
