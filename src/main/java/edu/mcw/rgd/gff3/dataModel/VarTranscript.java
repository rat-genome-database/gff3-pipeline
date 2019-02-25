package edu.mcw.rgd.gff3.dataModel;

/**
 * Created by hsnalabolu on 2/22/2019.
 */
public class VarTranscript{
    private int variantTranscriptID=0;
    private int transcriptRgdId=0;
    private String refAA;
    private String varAA;
    private String SynStat;
    private String LocName;
    private String NearSpliceSite;
    private String associatedGeneRGD;
    private String associatedGeneSym;
    private Polyphen ppObj;

    public Polyphen getPpObj() {
        return ppObj;
    }

    public void setPpObj(Polyphen ppObj) {
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

