package edu.mcw.rgd.gff3.dataModel;

/**
 * Created by hsnalabolu on 2/22/2019.
 */
public class Polyphen{
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

