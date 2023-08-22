package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.dao.impl.SampleDAO;
import edu.mcw.rgd.datamodel.Sample;
import edu.mcw.rgd.gff3.dataModel.Polyphen;
import edu.mcw.rgd.gff3.dataModel.VarTranscript;
import edu.mcw.rgd.gff3.dataModel.Variant;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author pjayaraman
 * Date: Dec 9, 2010
 * Time: 12:36:06 PM
 */
public class CreateGff4CarpeNovo2 {

    String toFile;
    String fileSource;

    public void createGff3ForPatient(int patientId) throws Exception{

        SampleDAO sampleDAO = new SampleDAO();
        sampleDAO.setDataSource(DataSourceFactory.getInstance().getCarpeNovoDataSource());
        List<Sample> samples = sampleDAO.getSamples(patientId);

        for (Sample sample : samples) {
            createGff3ForSample(sample.getId());
        }
    }

    public void createGff3ForSample(int sampleId) throws Exception{

        SampleDAO sampleDAO = new SampleDAO();
        sampleDAO.setDataSource(DataSourceFactory.getInstance().getCarpeNovoDataSource());
        Sample sample = sampleDAO.getSample(sampleId);
        String sampleName = sample.getAnalysisName().replace('/','_');
        String gffFile = getToFile()+sampleName+".gff3";
        String gffDamagingFile = getToFile()+sampleName+"_damaging.gff3";
        Connection conn = getConnection();
        creategff4CarpeNovo(gffFile, gffDamagingFile, sample.getMapKey(), conn, sampleId);

        conn.close();
    }

    void creategff4CarpeNovo(String gffFile, String gffDamagingFile, int mapKey, Connection conn, int sampleId) throws Exception{

        Gff3ColumnWriter gffWriter = new Gff3ColumnWriter(gffFile, false, Gff3ColumnWriter.COMPRESS_MODE_ZIP);
        Gff3ColumnWriter gffDmgVariantWriter = new Gff3ColumnWriter(gffDamagingFile, false, Gff3ColumnWriter.COMPRESS_MODE_ZIP);

        RgdGff3Dao dao = new RgdGff3Dao();
        SequenceRegionWatcher sequenceRegionWatcher1 = new SequenceRegionWatcher(mapKey, gffWriter, dao);
        SequenceRegionWatcher sequenceRegionWatcher2 = new SequenceRegionWatcher(mapKey, gffDmgVariantWriter, dao);

        // all chromosomes
        String[] chromosomes = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X", "Y", "MT"};
        for( String c: chromosomes ) {
            sequenceRegionWatcher1.emit(c);
            sequenceRegionWatcher2.emit(c);
        }

        process(sampleId, gffWriter, gffDmgVariantWriter, conn);

        gffWriter.close();
        gffDmgVariantWriter.close();

        gffWriter.sortInMemory();
        gffDmgVariantWriter.sortInMemory();
    }

    void process( int sampleId, Gff3ColumnWriter gffWriter, Gff3ColumnWriter gffDmgVariantWriter, Connection conn) throws Exception {

        CounterPool counters = new CounterPool();

        String getvarTranscript = "select vt.SYN_STATUS, vt.TRANSCRIPT_RGD_ID, vt.REF_AA, " +
                "vt.VAR_AA, g.RGD_ID, vt.LOCATION_NAME, vt.NEAR_SPLICE_SITE, g.GENE_SYMBOL " +
                "from (variant_transcript vt " +
                "inner join TRANSCRIPTS t on t.TRANSCRIPT_RGD_ID=vt.TRANSCRIPT_RGD_ID) " +
                "inner join GENES g on g.RGD_ID=t.GENE_RGD_ID where vt.VARIANT_RGD_ID=?";

        PreparedStatement findVarTran = conn.prepareStatement(getvarTranscript);


        String getpph = "select p.UNIPROT_ACC, p.PREDICTION, p.PROTEIN_ID, p.POSITION, p.AA1, p.AA2 " +
                "from polyphen p where p.VARIANT_RGD_ID=? and p.transcript_rgd_id=?";

        PreparedStatement findPph = conn.prepareStatement(getpph);


        String sql = "SELECT v.RGD_ID, m.CHROMOSOME, m.START_POS, m.END_POS, v.REF_NUC, v.VAR_NUC," +
                "s.TOTAL_DEPTH, s.VAR_FREQ, s.ZYGOSITY_STATUS, m.GENIC_STATUS " +
                "FROM variant v ,variant_sample_detail s,variant_map_data m " +
                "WHERE v.rgd_id=s.rgd_id AND m.rgd_id=v.rgd_id AND s.SAMPLE_ID=?";

        PreparedStatement findVar = conn.prepareStatement(sql);
        findVar.setInt(1, sampleId);

        ResultSet rsvar = findVar.executeQuery();
        while (rsvar.next()) {
            Variant v = new Variant();

            v.setVariantRgdId(rsvar.getInt(1));
            v.setChr(rsvar.getString(2).trim());
            v.setStart(rsvar.getInt(3));
            v.setStop(rsvar.getInt(4));
            v.setRef(rsvar.getString(5));
            v.setVar(rsvar.getString(6));
            v.setDepth(rsvar.getInt(7));
            v.setFreq(rsvar.getInt(8));
            v.setZygosity(rsvar.getString(9));
            v.setGenicStat(rsvar.getString(10));

            int cnt = counters.increment("varNumCount");
            if( cnt%10000 == 0 ) {
                System.out.println("  -- processing "+Utils.formatThousands(cnt)+" sampleId="+sampleId);
            }

            String chr = v.getChr();

            if (v.getDepth() > 10) {
                counters.increment("varDepthMoreThan10");
            } else if (v.getDepth() < 5) {
                counters.increment("varDepthLessThan5");
            }

            if (v.getGenicStat() != null) {
                if (v.getGenicStat().equalsIgnoreCase("GENIC")) {
                    counters.increment("varGenic");
                } else if (v.getGenicStat().equalsIgnoreCase("INTERGENIC")) {
                    counters.increment("varInterGenic");
                }
            }

            HashMap<Integer, VarTranscript> varTransHashObj = getVarGeneTransInfo(v.getVariantRgdId(), findVarTran, findPph);

            if (varTransHashObj.size() == 0) {

                counters.increment("varNotInTranscript");
                CreateGff4CarpeNovo.printNonTranscriptGFF3(chr, gffWriter, v, sampleId, getFileSource());

            } else {

                counters.increment("varInTranscript");
                CreateGff4CarpeNovo.printTranscriptGff3(chr, gffWriter, v, varTransHashObj, counters, sampleId, getFileSource());
                CreateGff4CarpeNovo.printDamagingTranscriptGff3(chr, gffDmgVariantWriter, v, varTransHashObj, sampleId, getFileSource());
            }
        }
        findVarTran.close();
        findPph.close();
        findVar.close();



        CreateGff4CarpeNovo.printStats(counters);
    }




    HashMap<Integer, VarTranscript> getVarGeneTransInfo(int varRgdId, PreparedStatement findVarTran, PreparedStatement findPph) throws Exception {

        HashMap<Integer, VarTranscript> newvarTransHash = new HashMap<>();

        findVarTran.setInt(1, varRgdId);
        ResultSet rsvarTran= findVarTran.executeQuery();

        while (rsvarTran.next()){
            int trRgdId = rsvarTran.getInt(2);

            VarTranscript newvarTransObj = new VarTranscript();
            newvarTransObj.setVariantRgdId(varRgdId);
            newvarTransObj.setTranscriptRgdId(trRgdId);

            findPph.setInt(1, varRgdId);
            findPph.setInt(2, trRgdId);
            ResultSet rspph= findPph.executeQuery();

            while(rspph.next()){
                Polyphen ppO = new Polyphen();

                String protid;
                if(!(rspph.getString(3)==null)){
                    protid= rspph.getString(3);
                }else{
                    protid = null;
                }
                ppO.setProteinID(protid);

                int pos = rspph.getInt(4);
                ppO.setPosition(pos);

                String aa1 = rspph.getString(5);
                ppO.setAA1(aa1);

                String aa2 = rspph.getString(6);
                ppO.setAA2(aa2);

                String prediction = rspph.getString(2);
                ppO.setPrediction(prediction);

                String uprotAcc = rspph.getString(1);
                ppO.setUniprotACC(uprotAcc);

                newvarTransObj.setPpObj(ppO);
            }
            rspph.close();

            String refaa = rsvarTran.getString(3);
            if(refaa!=null){
                refaa = refaa.trim();
            }else{
                refaa = null;
            }
            newvarTransObj.setRefAA(refaa);

            String varaa = rsvarTran.getString(4);
            if (null == varaa) {
                varaa = null;
            } else {
                varaa = varaa.trim();
            }
            newvarTransObj.setVarAA(varaa);

            String locname = rsvarTran.getString(6);
            if( locname!=null ){
                locname = locname.trim().replaceAll("\'","");
            }else{
                locname = null;
            }
            newvarTransObj.setLocName(locname);

            String nearss = rsvarTran.getString(7);
            if( nearss!=null ){
                nearss = nearss.trim();
            }else{
                nearss = null;
            }
            newvarTransObj.setNearSpliceSite(nearss);

            String gsym = rsvarTran.getString(8);
            if( gsym!=null ){
                gsym = gsym.trim();
            }else{
                gsym = null;
            }
            newvarTransObj.setAssociatedGeneSym(gsym);

            String grgd = rsvarTran.getString(5);
            if( grgd!=null ){
                grgd = grgd.trim();
            }else{
                grgd = null;
            }
            newvarTransObj.setAssociatedGeneRGD(grgd);

            String synStat = rsvarTran.getString(1);
            if( synStat!=null ){
                synStat = synStat.trim();
            }
            else{
                synStat = null;
            }
            newvarTransObj.setSynStat(synStat);
            newvarTransHash.put(trRgdId, newvarTransObj);
        }
        rsvarTran.close();

        return newvarTransHash;
    }

    public void setFileSource(String fileSource) {
        this.fileSource = fileSource;
    }

    public String getFileSource() {
        return fileSource;
    }



    public void setToFile(String toFile) {
        this.toFile = toFile;
    }

    public String getToFile() {
        return toFile;
    }

    Connection getConnection() throws Exception {
        return DataSourceFactory.getInstance().getCarpeNovoDataSource().getConnection();
    }

}