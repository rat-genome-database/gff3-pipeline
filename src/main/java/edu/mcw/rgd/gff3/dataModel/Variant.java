package edu.mcw.rgd.gff3.dataModel;

/**
 * Created by hsnalabolu on 2/22/2019.
 */
public class Variant {

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

