package edu.mcw.rgd.gff3;

/**
 * Created by IntelliJ IDEA.
 * User: pjayaraman
 * Date: 11/17/11
 * Time: 5:45 PM
 */
public class RGDInfo {
    public int getRgdId() {
        return rgdId;
    }

    public void setRgdId(int rgdId) {
        this.rgdId = rgdId;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public long getStartPos() {
        return startPos;
    }

    public void setStartPos(long startPos) {
        this.startPos = startPos;
    }

    public long getStopPos() {
        return stopPos;
    }

    public void setStopPos(long stopPos) {
        this.stopPos = stopPos;
    }

    public int rgdId;
    public String chromosome;
    public long startPos;
    public long stopPos;
}
