package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.dao.impl.SampleDAO;
import edu.mcw.rgd.datamodel.Sample;

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
public class CreateGff4CarpeNovo {

    String toFile;
    List<String> chromosomes;
    String fileSource;
    int sampleID;

    int varNumCount=0;
    int varGenic=0;
    int varInterGenic=0;
    int varInTranscript=0;
    int varNotInTranscript=0;
    int varDepthMoreThan10=0;
    int varDepthLessThan5=0;
    int varTransWithLocName=0;
    int varIn5UTR=0;
    int varIn3UTR=0;
    int varInINTRON=0;
    int varInEXON=0;
    int varWithAssocGene=0;
    int varWithSynStat=0;
    int varSyn=0;
    int varNonSyn=0;
    int varNearSpliceSite=0;
    int varWithPPPred=0;
    int varPPBenign=0;
    int varPPPossibly=0;
    int varPPProbably=0;


    public void createGff3ForPatient(int patientId) throws Exception{

        SampleDAO sampleDAO = new SampleDAO();
        sampleDAO.setDataSource(DataSourceFactory.getInstance().getCarpeNovoDataSource());

        List<Sample> samples = sampleDAO.getSamples(patientId);
        samples.parallelStream().forEach(sample-> {
            try{
                createGff3ForSample(sample.getId());
            }catch (Exception e){
                e.printStackTrace();
            }
        });
    }

    public void createGff3ForSample(int sampleId) throws Exception{
        this.sampleID = sampleId;

        Connection conn = getConnection();

        SampleDAO sampleDAO = new SampleDAO();
        sampleDAO.setDataSource(DataSourceFactory.getInstance().getCarpeNovoDataSource());
        Sample sample = sampleDAO.getSampleBySampleId(sampleId);
        String sampleName = sample.getAnalysisName().replace('/','_');

        String gffFile = getToFile()+sampleName+".gff3";
        creategff4CarpeNovo(gffFile, conn);

        conn.close();
    }

    void creategff4CarpeNovo(String gffFile, Connection conn) throws Exception{

        Gff3ColumnWriter gffWriter = new Gff3ColumnWriter(gffFile, false, true);

        for(String chr: getChromosomes()){

            varNumCount=0;
            varGenic=0;
            varInterGenic=0;
            varInTranscript=0;
            varNotInTranscript=0;
            varDepthMoreThan10=0;
            varDepthLessThan5=0;
            varTransWithLocName=0;
            varIn5UTR=0;
            varIn3UTR=0;
            varInINTRON=0;
            varInEXON=0;
            varWithAssocGene=0;
            varWithSynStat=0;
            varSyn=0;
            varNonSyn=0;
            varNearSpliceSite=0;
            varWithPPPred=0;
            varPPBenign=0;
            varPPPossibly=0;
            varPPProbably=0;


            HashMap<Integer, varTranscript> varTransHashObj;
            HashMap<Integer, variant> varhashOB = getVariants(chr);

            for(int variantid: varhashOB.keySet()){
                varNumCount++;

                variant v = varhashOB.get(variantid);

                if(v.getDepth()>10){
                    varDepthMoreThan10++;
                }else if(v.getDepth()<5){
                    varDepthLessThan5++;
                }

                if(v.getGenicStat()!=null){
                    if(v.getGenicStat().equalsIgnoreCase("GENIC")){
                        varGenic++;
                    }else
                    if(v.getGenicStat().equalsIgnoreCase("INTERGENIC")){
                        varInterGenic++;
                    }
                }

                varTransHashObj = getVarGeneTransInfo(variantid, conn);

                if(varTransHashObj.size()==0){

                    varNotInTranscript++;
                    printNonTranscriptGFF3(chr, gffWriter, v);

                }else{

                    varInTranscript++;
                    printTranscriptGff3(chr, gffWriter, v, varTransHashObj);
                }

            }

            System.out.println("\nChromosome:" + chr + " stats-");
            System.out.println("Variants processed:" + varNumCount);
            System.out.println("Genic Variants processed:" + varGenic);
            System.out.println("Intergenic Variants processed:" + varInterGenic);
            System.out.println("Variants near Splice Sites:" + varNearSpliceSite);
            System.out.println("Variants with depth > 10:" + varDepthMoreThan10);
            System.out.println("Variants with depth < 5:" + varDepthLessThan5);
            System.out.println("Variants NOT associated with ANY transcript:" + varNotInTranscript);
            System.out.println("Variants associated with one or more Transcripts:" + varInTranscript);
            System.out.println("Genic Variants with a defined LocationName:" + varTransWithLocName);
            System.out.println("Genic Variants with an associated Gene Symbol:" + varWithAssocGene);
            System.out.println("Variants in INTRON:" + varInINTRON);
            System.out.println("Variants in 3UTRS:" + varIn3UTR);
            System.out.println("Variants in 5UTRS:" + varIn5UTR);
            System.out.println("Variants in EXON:" + varInEXON);
            System.out.println("Genic Variants with a synonymous/nonsynonymous status:" + varWithSynStat);
            System.out.println("Genic Variants that are synonymous:" + varSyn);
            System.out.println("Genic Variants that are NONsynonymous:" + varNonSyn);
            System.out.println("Genic Variants with PolyPhen Predictions:" + varWithPPPred);
            System.out.println("BENIGN PolyPhen Predictions:" + varPPBenign);
            System.out.println("POSSIBLY DAMAGING PolyPhen Predictions:" + varPPPossibly);
            System.out.println("PROBABLY DAMAGING PolyPhen Predictions:" + varPPProbably);
            System.out.println("DONE with Chromosome:" + chr);
            System.out.println("GFF3 File SUCCESSFUL!");
        }

        gffWriter.close();
    }



    void printTranscriptGff3(String chr, Gff3ColumnWriter gffWriter, variant varOB, HashMap<Integer, varTranscript> varTransHashObj) throws Exception {

        String transcriptId="";
        String transcriptRgdId="";
        String refAA="";
        String varAA="";
        String associatedGene="";
        String varLoc = "";
        String nearSpliceSite="";
        String transcriptPolyPred="";
        String synStatus="";

        for(int transcriptID : varTransHashObj.keySet() ){
            varTranscript vt = varTransHashObj.get(transcriptID);
            if(!vt.getLocName().equals("NA")){
                transcriptId += transcriptID+",";
                transcriptRgdId+=vt.getTranscriptRgdId()+",";
                refAA+=vt.getRefAA()+",";
                varAA+=vt.getVarAA()+",";
                associatedGene+=vt.getAssociatedGeneSym()+":"+
                        vt.getAssociatedGeneRGD()+",";
                varLoc+=vt.getTranscriptRgdId()+":"+
                        vt.getLocName().replaceAll(",", "||")+",";
                synStatus+=vt.getSynStat()+",";
                nearSpliceSite+=vt.getNearSpliceSite()+",";

                if(vt.getPpObj()!=null){
                    transcriptPolyPred+=vt.getTranscriptRgdId()+":"+
                            vt.getPpObj().getProteinID()+"||"+
                            vt.getPpObj().getPosition()+"||"+
                            vt.getPpObj().getAA1()+"||"+
                            vt.getPpObj().getAA2()+"||"+
                            vt.getPpObj().getPrediction()+"||"+
                            vt.getPpObj().getUniprotACC()+",";


                    if(vt.getPpObj().getPrediction().equalsIgnoreCase("benign")){
                        varPPBenign++;
                    }else
                    if(vt.getPpObj().getPrediction().equalsIgnoreCase("possibly damaging")){
                        varPPPossibly++;
                    }else
                    if(vt.getPpObj().getPrediction().equalsIgnoreCase("probably damaging")){
                        varPPProbably++;
                    }
                }else{
                    transcriptPolyPred+="NA,";
                }

                if(vt.getLocName().contains("INTRON")){
                    varInINTRON++;
                }
                if(vt.getLocName().contains("EXON")){
                    varInEXON++;
                }
                if(vt.getLocName().contains("5UTRS")){
                    varIn5UTR++;
                }
                if(vt.getLocName().contains("3UTRS")){
                    varIn3UTR++;
                }

                if(!(vt.getAssociatedGeneSym().equals("NA"))){
                    varWithAssocGene++;
                }

                if(!(vt.getSynStat().equals("NA"))){
                    varWithSynStat++;
                    if(vt.getSynStat().equalsIgnoreCase("synonymous")){
                        varSyn++;
                    }else
                    if(vt.getSynStat().equalsIgnoreCase("nonsynonymous")){
                        varNonSyn++;
                    }
                }

                if(!(vt.getNearSpliceSite().equals("NA"))){
                    varNearSpliceSite++;
                }

                if(vt.getPpObj()!=null){
                    varWithPPPred++;
                }
            }
        }

        Map<String,String> attributesHashMap = new HashMap<String, String>();
        attributesHashMap.put("ID", String.valueOf(varOB.getVariantID()));
        attributesHashMap.put("Name", String.valueOf(varOB.getVariantID()));
        attributesHashMap.put("Alias", String.valueOf(varOB.getVariantID()));
        attributesHashMap.put("Reference", varOB.getRef());
        attributesHashMap.put("Variant", varOB.getVar());
        attributesHashMap.put("Depth", String.valueOf(varOB.getDepth()));
        attributesHashMap.put("Frequency", String.valueOf(varOB.getFreq()));
        attributesHashMap.put("VariantTranscriptID", transcriptId.substring(0, transcriptId.length()-1));
        attributesHashMap.put("transcriptRGDID", transcriptRgdId.substring(0, transcriptRgdId.length()-1));
        attributesHashMap.put("refAA", refAA.substring(0, refAA.length()-1));
        attributesHashMap.put("varAA", varAA.substring(0, varAA.length()-1));
        attributesHashMap.put("associatedGene", associatedGene.substring(0, associatedGene.length()-1));
        attributesHashMap.put("synStatus", synStatus.substring(0, synStatus.length()-1));
        attributesHashMap.put("locName", varLoc.substring(0, varLoc.length()-1));
        attributesHashMap.put("nearSpliceSite", nearSpliceSite.substring(0, nearSpliceSite.length()-1));
        attributesHashMap.put("polyPred", transcriptPolyPred.substring(0, transcriptPolyPred.length()-1));


        gffWriter.writeFirst8Columns(chr,getFileSource(),"SNV_"+this.sampleID,varOB.getStart(),varOB.getStart(),".",".",".");
        gffWriter.writeAttributes4Gff3(attributesHashMap);
    }

    void printNonTranscriptGFF3(String chr, Gff3ColumnWriter gffWriter, variant varOB) throws Exception {
        Map<String, String> attributeHashMap = new HashMap<String, String>();
        attributeHashMap.put("ID", String.valueOf(varOB.getVariantID()));
        attributeHashMap.put("Name", String.valueOf(varOB.getVariantID()));
        attributeHashMap.put("Alias", String.valueOf(varOB.getVariantID()));
        attributeHashMap.put("Reference", varOB.getRef());
        attributeHashMap.put("Variant", varOB.getVar());
        attributeHashMap.put("Depth", String.valueOf(varOB.getDepth()));
        attributeHashMap.put("Frequency", String.valueOf(varOB.getFreq()));
        attributeHashMap.put("nearSpliceSite","NA");
        attributeHashMap.put("synStatus", "NA");
        attributeHashMap.put("Index","1");

        gffWriter.writeGff3AllColumns(chr, getFileSource(), "SNV_" + this.sampleID, varOB.getStart(), varOB.getStart(), ".", ".", ".", attributeHashMap);
    }


    HashMap<Integer, variant> getVariants(String chrNum) throws Exception{

        HashMap<Integer, variant> newvarHash = new HashMap<Integer, variant>();

        Connection connection = getConnection();
        // System.out.println("here is the connection:" + connection);
        String findAllVariants = "SELECT v.VARIANT_ID, v.CHROMOSOME, v.START_POS, v.END_POS, v.REF_NUC, v.VAR_NUC," +
                "TOTAL_DEPTH, v.VAR_FREQ, v.ZYGOSITY_STATUS, v.GENIC_STATUS FROM VARIANT v " +
                "WHERE v.SAMPLE_ID = ? and v.CHROMOSOME = ? ";

        PreparedStatement findVar = connection.prepareStatement(findAllVariants);
        findVar.setInt(1, sampleID);
        findVar.setString(2, chrNum);

        ResultSet rsvar= findVar.executeQuery();
        while (rsvar.next()) {
            variant newVar = new variant();

            newVar.setVariantID(rsvar.getInt(1));
            newVar.setChr(rsvar.getString(2).trim());
            newVar.setStart(rsvar.getInt(3));
            newVar.setStop(rsvar.getInt(4));
            newVar.setRef(rsvar.getString(5));
            newVar.setVar(rsvar.getString(6));
            newVar.setDepth(rsvar.getInt(7));
            newVar.setFreq(rsvar.getInt(8));
            newVar.setZygosity(rsvar.getString(9));
            newVar.setGenicStat(rsvar.getString(10));

            newvarHash.put(newVar.getVariantID(), newVar);
        }

        connection.close();
        return newvarHash;
    }


    HashMap<Integer, varTranscript> getVarGeneTransInfo(int varObj, Connection conn) throws Exception {

        varTranscript newvarTransObj = new varTranscript();
        HashMap<Integer, varTranscript> newvarTransHash = new HashMap<Integer, varTranscript>();

        String getvarTranscript = "select vt.VARIANT_TRANSCRIPT_ID, vt.SYN_STATUS, vt.TRANSCRIPT_RGD_ID, vt.REF_AA, " +
                "vt.VAR_AA, g.RGD_ID, vt.LOCATION_NAME, vt.NEAR_SPLICE_SITE, g.GENE_SYMBOL " +
                "from (VARIANT_TRANSCRIPT vt " +
                "inner join TRANSCRIPTS t on t.TRANSCRIPT_RGD_ID=vt.TRANSCRIPT_RGD_ID) " +
                "inner join GENES g on g.RGD_ID=t.GENE_RGD_ID where vt.VARIANT_ID=?";

        PreparedStatement findVarTran = conn.prepareStatement(getvarTranscript);


        String getpph = "select p.UNIPROT_ACC, p.PREDICTION, p.PROTEIN_ID, p.POSITION, p.AA1, p.AA2 " +
                "from POLYPHEN p where p.VARIANT_TRANSCRIPT_ID=?";

        PreparedStatement findPph = conn.prepareStatement(getpph);

        findVarTran.setInt(1, varObj);
        ResultSet rsvarTran= findVarTran.executeQuery();

        while (rsvarTran.next()){
            int vartranid = rsvarTran.getInt(1);
            newvarTransObj.setVariantTranscriptID(vartranid);

            findPph.setInt(1, vartranid);
            ResultSet rspph= findPph.executeQuery();

            while(rspph.next()){
                polyphenObj ppO = new polyphenObj();

                String protid;
                if(!(rspph.getString(3)==null)){
                    protid= rspph.getString(3);
                }else{
                    protid = "NA";
                }
                ppO.setProteinID(protid);

                int pos = rspph.getInt(4);
                ppO.setPosition(pos);

                String aa1 = rspph.getString(5);
                ppO.setAA1(aa1);

                String aa2 = rspph.getString(6);
                ppO.setAA2(aa2);

                String prediction;
                if(rspph.getString(2)==null){
                    prediction="null";
                }else{
                    prediction=rspph.getString(2);
                }
                ppO.setPrediction(prediction);

                String uprotAcc = rspph.getString(1);
                ppO.setUniprotACC(uprotAcc);

                newvarTransObj.setPpObj(ppO);
            }

            int tranRgdid = rsvarTran.getInt(3);
            newvarTransObj.setTranscriptRgdId(tranRgdid);

            String refaa;
            if(!(rsvarTran.getString(4)==null)){
                refaa = rsvarTran.getString(4).trim();
            }else{
                refaa = "NA";
            }
            newvarTransObj.setRefAA(refaa);

            String varaa;
            if(!(rsvarTran.getString(5)==null)){
                varaa = rsvarTran.getString(5).trim();
            }else{
                varaa = "NA";
            }
            newvarTransObj.setVarAA(varaa);

            String locname;
            if(!(rsvarTran.getString(7)==null)){
                locname = rsvarTran.getString(7).trim().replaceAll("\'","");
            }else{
                locname = "NA";
            }
            newvarTransObj.setLocName(locname);

            String nearss;
            if(!(rsvarTran.getString(8)==null)){
                nearss = rsvarTran.getString(8).trim();
            }else{
                nearss = "NA";
            }
            newvarTransObj.setNearSpliceSite(nearss);

            String gsym;
            if(!(rsvarTran.getString(9)==null)){
                gsym = rsvarTran.getString(9).trim();
            }else{
                gsym = "NA";
            }
            newvarTransObj.setAssociatedGeneSym(gsym);

            String grgd;
            if(!(rsvarTran.getString(6)==null)){
                grgd = rsvarTran.getString(6).trim();
            }else{
                grgd = "NA";
            }
            newvarTransObj.setAssociatedGeneRGD(grgd);

            String synStat;
            if(!(rsvarTran.getString(2)==null)){
                synStat = rsvarTran.getString(2).trim();
            }
            else{
                synStat = "NA";
            }
            newvarTransObj.setSynStat(synStat);
            newvarTransHash.put(newvarTransObj.getVariantTranscriptID(), newvarTransObj);
        }


        findVarTran.close();
        findPph.close();
        return newvarTransHash;
    }

    public void setFileSource(String fileSource) {
        this.fileSource = fileSource;
    }

    public String getFileSource() {
        return fileSource;
    }


    class polyphenObj{
        String proteinID;
        int position=0;
        String AA1;
        String AA2;
        String prediction;
        String UniprotACC;

        public String getProteinID() {
            return proteinID;
        }

        public void setProteinID(String proteinID) {
            this.proteinID = proteinID;
        }

        public int getPosition() {
            return position;
        }

        public void setPosition(int position) {
            this.position = position;
        }

        public String getAA1() {
            return AA1;
        }

        public void setAA1(String AA1) {
            this.AA1 = AA1;
        }

        public String getAA2() {
            return AA2;
        }

        public void setAA2(String AA2) {
            this.AA2 = AA2;
        }

        public String getPrediction() {
            return prediction;
        }

        public void setPrediction(String prediction) {
            this.prediction = prediction;
        }

        public String getUniprotACC() {
            return UniprotACC;
        }

        public void setUniprotACC(String uniprotACC) {
            UniprotACC = uniprotACC;
        }
    }

    class varTranscript{
        private int variantTranscriptID=0;
        private int transcriptRgdId=0;
        private String refAA;
        private String varAA;
        private String SynStat;
        private String LocName;
        private String NearSpliceSite;
        private String associatedGeneRGD;
        private String associatedGeneSym;
        private polyphenObj ppObj;

        public polyphenObj getPpObj() {
            return ppObj;
        }

        public void setPpObj(polyphenObj ppObj) {
            this.ppObj = ppObj;
        }

        public int getVariantTranscriptID() {
            return variantTranscriptID;
        }

        public void setVariantTranscriptID(int variantTranscriptID) {
            this.variantTranscriptID = variantTranscriptID;
        }

        public int getTranscriptRgdId() {
            return transcriptRgdId;
        }

        public void setTranscriptRgdId(int transcriptRgdId) {
            this.transcriptRgdId = transcriptRgdId;
        }

        public String getRefAA() {
            return refAA;
        }

        public void setRefAA(String refAA) {
            this.refAA = refAA;
        }

        public String getVarAA() {
            return varAA;
        }

        public void setVarAA(String varAA) {
            this.varAA = varAA;
        }

        public String getSynStat() {
            return SynStat;
        }

        public void setSynStat(String synStat) {
            SynStat = synStat;
        }

        public String getLocName() {
            return LocName;
        }

        public void setLocName(String locName) {
            LocName = locName;
        }

        public String getNearSpliceSite() {
            return NearSpliceSite;
        }

        public void setNearSpliceSite(String nearSpliceSite) {
            NearSpliceSite = nearSpliceSite;
        }

        public String getAssociatedGeneRGD() {
            return associatedGeneRGD;
        }

        public void setAssociatedGeneRGD(String associatedGeneRGD) {
            this.associatedGeneRGD = associatedGeneRGD;
        }

        public String getAssociatedGeneSym() {
            return associatedGeneSym;
        }

        public void setAssociatedGeneSym(String associatedGeneSym) {
            this.associatedGeneSym = associatedGeneSym;
        }
    }


    class variant {

        private int variantID=0;
        private String chr;
        private int start=0;
        private int stop=0;
        private String ref;
        private String var;
        private int depth=0;
        private int freq=0;
        private String genicStat;
        private String zygosity;

        public String getGenicStat() {
            return genicStat;
        }

        public void setGenicStat(String genicStat) {
            this.genicStat = genicStat;
        }

        public String getZygosity() {
            return zygosity;
        }

        public void setZygosity(String zygosity) {
            this.zygosity = zygosity;
        }

        public int getVariantID() {
            return variantID;
        }

        public void setVariantID(int variantID) {
            this.variantID = variantID;
        }

        public String getChr() {
            return chr;
        }

        public void setChr(String chr) {
            this.chr = chr;
        }

        public int getStart() {
            return start;
        }

        public void setStart(int start) {
            this.start = start;
        }

        public int getStop() {
            return stop;
        }

        public void setStop(int stop) {
            this.stop = stop;
        }

        public String getRef() {
            return ref;
        }

        public void setRef(String ref) {
            this.ref = ref;
        }

        public String getVar() {
            return var;
        }

        public void setVar(String var) {
            this.var = var;
        }

        public int getDepth() {
            return depth;
        }

        public void setDepth(int depth) {
            this.depth = depth;
        }

        public int getFreq() {
            return freq;
        }

        public void setFreq(int freq) {
            this.freq = freq;
        }
    }

    public void setToFile(String toFile) {
        this.toFile = toFile;
    }

    public String getToFile() {
        return toFile;
    }

    public List<String> getChromosomes() {
        return chromosomes;
    }

    public void setChromosomes(List<String> chromosomes) {
        this.chromosomes = chromosomes;
    }

    Connection getConnection() throws Exception {
        return DataSourceFactory.getInstance().getCarpeNovoDataSource().getConnection();
    }

}
