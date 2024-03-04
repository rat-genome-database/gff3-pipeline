package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.dao.impl.SampleDAO;
import edu.mcw.rgd.datamodel.Sample;
import edu.mcw.rgd.gff3.dataModel.Polyphen;
import edu.mcw.rgd.gff3.dataModel.VarTranscript;
import edu.mcw.rgd.gff3.dataModel.Variant;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.*;
import java.util.concurrent.ForkJoinPool;

/**
 * @author pjayaraman
 * Date: Dec 9, 2010
 * Time: 12:36:06 PM
 */
public class CreateGff4CarpeNovo {

    String toFile;
    String fileSource;

    public void createGff3ForPatient(int patientId) throws Exception{

        SampleDAO sampleDAO = new SampleDAO();
        sampleDAO.setDataSource(DataSourceFactory.getInstance().getCarpeNovoDataSource());
        List<Sample> samples = sampleDAO.getSamples(patientId);

        int maxThreadCount = 10;
        {
            ForkJoinPool customThreadPool = new ForkJoinPool(maxThreadCount);
            customThreadPool.submit(() -> samples.parallelStream().forEach(sample -> {
                try {
                    createGff3ForSample(sample.getId());
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            })).get();
            customThreadPool.shutdown();
        }
    }

    public void createGff3ForSample(int sampleId) throws Exception{

        long time0 = System.currentTimeMillis();

        SampleDAO sampleDAO = new SampleDAO();
        sampleDAO.setDataSource(DataSourceFactory.getInstance().getCarpeNovoDataSource());
        Sample sample = sampleDAO.getSample(sampleId);
        String sampleName = sample.getAnalysisName().replace('/','_');
        String gffFile = getToFile()+sampleName+".gff3";
        String gffDamagingFile = getToFile()+sampleName+"_damaging.gff3";
        Connection conn = getConnection();
        creategff4CarpeNovo(gffFile, gffDamagingFile, sample.getMapKey(), conn, sampleId);

        conn.close();

        System.out.println("### SAMPLE "+sampleId+" OK!    elapsed: "+Utils.formatElapsedTime(time0, System.currentTimeMillis()));
    }

    void creategff4CarpeNovo(String gffFile, String gffDamagingFile, int mapKey, Connection conn, int sampleId) throws Exception{

        Gff3ColumnWriter gffWriter = new Gff3ColumnWriter(gffFile, Gff3ColumnWriter.COMPRESS_MODE_ZIP);
        Gff3ColumnWriter gffDmgVariantWriter = new Gff3ColumnWriter(gffDamagingFile, Gff3ColumnWriter.COMPRESS_MODE_ZIP);

        String msg = "#RGD SAMPLE ID: "+sampleId+"\n";
        gffWriter.print(msg);
        gffDmgVariantWriter.print(msg);

        String assemblyName = MapManager.getInstance().getMap(mapKey).getRefSeqAssemblyName();
        msg = "#ASSEMBLY: "+assemblyName+"\n";
        gffWriter.print(msg);
        gffDmgVariantWriter.print(msg);


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
                "s.TOTAL_DEPTH, s.VAR_FREQ, s.ZYGOSITY_STATUS, m.GENIC_STATUS, v.variant_type " +
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
            v.setVarType(rsvar.getString(11));

            // positions for SNVs: END_POS must be always the same as START_POS
            if( v.getVarType().equalsIgnoreCase("snv") ) {
                v.setStop( v.getStart() );
            }
            // MNVs: subtract 1 from stop pos
            if( v.getVarType().equalsIgnoreCase("mnv") ) {
                int newStopPos = v.getStop()-1;
                if( newStopPos >= v.getStart() ) {
                    v.setStop(newStopPos);
                }
            }
            // deletions
            if( v.getVarType().equalsIgnoreCase("del") || v.getVarType().equalsIgnoreCase("deletion") ) {
                int newStopPos = v.getStart() + v.getRef().length() - 1;
                if( newStopPos != v.getStop() ) {
                    v.setStop(newStopPos);
                }
            }

            int cnt = counters.increment("varNumCount");
            if( cnt%100000 == 0 ) {
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
                printNonTranscriptGFF3(chr, gffWriter, v, getFileSource());

            } else {

                counters.increment("varInTranscript");
                printTranscriptGff3(chr, gffWriter, v, varTransHashObj, counters, getFileSource());
                printDamagingTranscriptGff3(chr, gffDmgVariantWriter, v, varTransHashObj, getFileSource());
            }
        }
        findVarTran.close();
        findPph.close();
        findVar.close();

        printStats(counters, sampleId);
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

    synchronized static void printStats( CounterPool counters, int sampleId ) {
        System.out.println("SAMPLE "+sampleId+": Variants processed:" + counters.get("varNumCount"));
        System.out.println("SAMPLE "+sampleId+": Genic Variants processed:" + counters.get("varGenic"));
        System.out.println("SAMPLE "+sampleId+": Intergenic Variants processed:" + counters.get("varInterGenic"));
        System.out.println("SAMPLE "+sampleId+": Variants near Splice Sites:" + counters.get("varNearSpliceSite"));
        System.out.println("SAMPLE "+sampleId+": Variants with depth > 10:" + counters.get("varDepthMoreThan10"));
        System.out.println("SAMPLE "+sampleId+": Variants with depth < 5:" + counters.get("varDepthLessThan5"));
        System.out.println("SAMPLE "+sampleId+": Variants NOT associated with ANY transcript:" + counters.get("varNotInTranscript"));
        System.out.println("SAMPLE "+sampleId+": Variants associated with one or more Transcripts:" + counters.get("varInTranscript"));
        System.out.println("SAMPLE "+sampleId+": Genic Variants with a defined LocationName:" + counters.get("varTransWithLocName"));
        System.out.println("SAMPLE "+sampleId+": Genic Variants with an associated Gene Symbol:" + counters.get("varWithAssocGene"));
        System.out.println("SAMPLE "+sampleId+": Variants in INTRON:" + counters.get("varInINTRON"));
        System.out.println("SAMPLE "+sampleId+": Variants in 3UTRS:" + counters.get("varIn3UTR"));
        System.out.println("SAMPLE "+sampleId+": Variants in 5UTRS:" + counters.get("varIn5UTR"));
        System.out.println("SAMPLE "+sampleId+": Variants in EXON:" + counters.get("varInEXON"));
        System.out.println("SAMPLE "+sampleId+": Genic Variants with a synonymous/nonsynonymous status:" + counters.get("varWithSynStat"));
        System.out.println("SAMPLE "+sampleId+": Genic Variants that are synonymous:" + counters.get("varSyn"));
        System.out.println("SAMPLE "+sampleId+": Genic Variants that are NONsynonymous:" + counters.get("varNonSyn"));
        System.out.println("SAMPLE "+sampleId+": Genic Variants with PolyPhen Predictions:" + counters.get("varWithPPPred"));
        System.out.println("SAMPLE "+sampleId+": BENIGN PolyPhen Predictions:" + counters.get("varPPBenign"));
        System.out.println("SAMPLE "+sampleId+": POSSIBLY DAMAGING PolyPhen Predictions:" + counters.get("varPPPossibly"));
        System.out.println("SAMPLE "+sampleId+": PROBABLY DAMAGING PolyPhen Predictions:" + counters.get("varPPProbably"));
        System.out.println("SAMPLE "+sampleId+": =====");
    }

    static void printTranscriptGff3(String chr, Gff3ColumnWriter gffWriter, Variant varOB, HashMap<Integer, VarTranscript> varTransHashObj,
                                    CounterPool counters, String source) throws Exception {

        String transcriptRgdId="";
        String refAA="";
        String varAA="";
        String associatedGene="";
        String varLoc = "";
        String nearSpliceSite="";
        String transcriptPolyPred="";
        String synStatus="";
        HashMap<String,String> attributesHashMap = new HashMap<String, String>();

        attributesHashMap.put("ID", String.valueOf(varOB.getVariantRgdId()));
        attributesHashMap.put("Name", String.valueOf(varOB.getVariantRgdId()));
        attributesHashMap.put("reference", Utils.NVL(varOB.getRef(), "-"));
        attributesHashMap.put("variant", Utils.NVL(varOB.getVar(), "-"));
        attributesHashMap.put("depth", String.valueOf(varOB.getDepth()));
        attributesHashMap.put("frequency", String.valueOf(varOB.getFreq()));

        for(int transcriptID : varTransHashObj.keySet() ){
            VarTranscript vt = varTransHashObj.get(transcriptID);
            if( vt.getLocName()!=null ){
                transcriptRgdId+=vt.getTranscriptRgdId()+",";
                refAA+=vt.getRefAA()+",";
                varAA+=vt.getVarAA()+",";
                associatedGene+=vt.getAssociatedGeneSym()+":"+vt.getAssociatedGeneRGD()+",";
                varLoc+=vt.getTranscriptRgdId()+":"+vt.getLocName().replaceAll(",", "||")+",";
                synStatus+=vt.getSynStat()+",";
                nearSpliceSite+=vt.getNearSpliceSite()+",";

                if(vt.getPpObj()!=null){
                    counters.increment("varWithPPPred");
                    transcriptPolyPred+=vt.getTranscriptRgdId()+":"+
                            vt.getPpObj().getProteinID()+"||"+
                            vt.getPpObj().getPosition()+"||"+
                            vt.getPpObj().getAA1()+"||"+
                            vt.getPpObj().getAA2()+"||"+
                            vt.getPpObj().getPrediction()+"||"+
                            vt.getPpObj().getUniprotACC()+",";


                    if(vt.getPpObj().getPrediction().equalsIgnoreCase("benign")){
                        counters.increment("varPPBenign");
                    }else
                    if(vt.getPpObj().getPrediction().equalsIgnoreCase("possibly damaging")){
                        counters.increment("varPPPossibly");
                    }else
                    if(vt.getPpObj().getPrediction().equalsIgnoreCase("probably damaging")){
                        counters.increment("varPPProbably");
                    }
                }else{
                    transcriptPolyPred+="NA,";
                }

                if(vt.getLocName().contains("INTRON")){
                    counters.increment("varInINTRON");
                }
                if(vt.getLocName().contains("EXON")){
                    counters.increment("varInEXON");
                }
                if(vt.getLocName().contains("5UTRS")){
                    counters.increment("varIn5UTR");
                }
                if(vt.getLocName().contains("3UTRS")){
                    counters.increment("varIn3UTR");
                }

                if( vt.getAssociatedGeneSym()!=null ){
                    counters.increment("varWithAssocGene");
                }

                if( vt.getSynStat()!=null ){
                    counters.increment("varWithSynStat");
                    if(vt.getSynStat().equalsIgnoreCase("synonymous")){
                        counters.increment("varSyn");
                    }else
                    if(vt.getSynStat().equalsIgnoreCase("nonsynonymous")){
                        counters.increment("varNonSyn");
                    }
                }

                if( vt.getNearSpliceSite()!=null ){
                    counters.increment("varNearSpliceSite");
                }
            }
        }

        attributesHashMap.put("transcriptRGDID", transcriptRgdId.substring(0, transcriptRgdId.length()-1));
        attributesHashMap.put("refAA", refAA.substring(0, refAA.length()-1));
        attributesHashMap.put("varAA", varAA.substring(0, varAA.length()-1));
        attributesHashMap.put("associatedGene", associatedGene.substring(0, associatedGene.length()-1));
        attributesHashMap.put("synStatus", synStatus.substring(0, synStatus.length()-1));
        attributesHashMap.put("locName", varLoc.substring(0, varLoc.length()-1));
        attributesHashMap.put("nearSpliceSite", nearSpliceSite.substring(0, nearSpliceSite.length()-1));
        attributesHashMap.put("polyPred", transcriptPolyPred.substring(0, transcriptPolyPred.length()-1));


        gffWriter.writeFirst8Columns(chr, source,getSoVarType(varOB.getVarType()), varOB.getStart(), varOB.getStop(),".",".",".");
        gffWriter.writeAttributes4Gff3(attributesHashMap);
    }

    static void printDamagingTranscriptGff3(String chr, Gff3ColumnWriter gffDmgVariantWriter, Variant varOB,
                                            HashMap<Integer, VarTranscript> varTransHashObj, String source) throws Exception {

        String transcriptRgdId="";
        String refAA="";
        String varAA="";
        String associatedGene="";
        String varLoc = "";
        String nearSpliceSite="";
        String transcriptPolyPred="";
        String synStatus="";
        boolean damaging= false;
        for(int transcriptID : varTransHashObj.keySet() ) {
            VarTranscript vt = varTransHashObj.get(transcriptID);

            if( vt.getLocName()!=null ) {
                if (vt.getPpObj() != null && (vt.getPpObj().getPrediction().equalsIgnoreCase("possibly damaging") ||
                        vt.getPpObj().getPrediction().equalsIgnoreCase("probably damaging"))) {
                    damaging = true;
                    transcriptRgdId += vt.getTranscriptRgdId() + ",";
                    refAA += vt.getRefAA() + ",";
                    varAA += vt.getVarAA() + ",";
                    associatedGene += vt.getAssociatedGeneSym() + ":" +
                            vt.getAssociatedGeneRGD() + ",";
                    varLoc += vt.getTranscriptRgdId() + ":" +
                            vt.getLocName().replaceAll(",", "||") + ",";
                    synStatus += vt.getSynStat() + ",";
                    nearSpliceSite += vt.getNearSpliceSite() + ",";
                    transcriptPolyPred += vt.getTranscriptRgdId() + ":" +
                            vt.getPpObj().getProteinID() + "||" +
                            vt.getPpObj().getPosition() + "||" +
                            vt.getPpObj().getAA1() + "||" +
                            vt.getPpObj().getAA2() + "||" +
                            vt.getPpObj().getPrediction() + "||" +
                            vt.getPpObj().getUniprotACC() + ",";
                }
            }
        }
        if(damaging) {
            Map<String,String> attributesHashMap = new HashMap<String, String>();

            attributesHashMap.put("ID", String.valueOf(varOB.getVariantRgdId()));
            attributesHashMap.put("Name", String.valueOf(varOB.getVariantRgdId()));
            attributesHashMap.put("reference", Utils.NVL(varOB.getRef(), "-"));
            attributesHashMap.put("variant", Utils.NVL(varOB.getVar(), "-"));
            attributesHashMap.put("depth", String.valueOf(varOB.getDepth()));
            attributesHashMap.put("frequency", String.valueOf(varOB.getFreq()));
            attributesHashMap.put("transcriptRGDID", transcriptRgdId.substring(0, transcriptRgdId.length() - 1));
            attributesHashMap.put("refAA", refAA.substring(0, refAA.length() - 1));
            attributesHashMap.put("varAA", varAA.substring(0, varAA.length() - 1));
            attributesHashMap.put("associatedGene", associatedGene.substring(0, associatedGene.length() - 1));
            attributesHashMap.put("synStatus", synStatus.substring(0, synStatus.length() - 1));
            attributesHashMap.put("locName", varLoc.substring(0, varLoc.length() - 1));
            attributesHashMap.put("nearSpliceSite", nearSpliceSite.substring(0, nearSpliceSite.length() - 1));
            attributesHashMap.put("polyPred", transcriptPolyPred.substring(0, transcriptPolyPred.length() - 1));

            gffDmgVariantWriter.writeFirst8Columns(chr, source, getSoVarType(varOB.getVarType()), varOB.getStart(), varOB.getStop(), ".", ".", ".");
            gffDmgVariantWriter.writeAttributes4Gff3(attributesHashMap);
        }
    }

    static void printNonTranscriptGFF3(String chr, Gff3ColumnWriter gffWriter, Variant varOB, String source) throws Exception {
        Map<String, String> attributeHashMap = new HashMap<>();

        attributeHashMap.put("ID", String.valueOf(varOB.getVariantRgdId()));
        attributeHashMap.put("Name", String.valueOf(varOB.getVariantRgdId()));
        attributeHashMap.put("reference", Utils.NVL(varOB.getRef(), "-"));
        attributeHashMap.put("variant", Utils.NVL(varOB.getVar(), "-"));
        attributeHashMap.put("depth", String.valueOf(varOB.getDepth()));
        attributeHashMap.put("frequency", String.valueOf(varOB.getFreq()));

        gffWriter.writeGff3AllColumns(chr, source, getSoVarType(varOB.getVarType()), varOB.getStart(), varOB.getStop(), ".", ".", ".", attributeHashMap);
    }

    static String getSoVarType( String varTypeInRgd ) {

        return switch (varTypeInRgd) {
            case "snv" -> "SNV";
            case "ins", "insertion" -> "insertion";
            case "del", "deletion" -> "deletion";
            case "delins" -> "delins";
            case "mnv" -> "MNV";
            case "tandem_repeat" -> "tandem_repeat";
            case "sequence_alteration" -> "sequence_alteration";
            default -> varTypeInRgd;
        };
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