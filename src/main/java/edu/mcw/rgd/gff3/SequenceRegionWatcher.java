package edu.mcw.rgd.gff3;

import java.util.HashSet;
import java.util.Set;

public class SequenceRegionWatcher {

    private Set<String> encounteredChromosomes = new HashSet<>();

    private int mapKey;
    private Gff3ColumnWriter gff3Writer;
    private RgdGff3Dao dao;

    private String chromosomePrefix = "Chr";

    public SequenceRegionWatcher( int mapKey, Gff3ColumnWriter gff3Writer, RgdGff3Dao dao ) {
        this.mapKey = mapKey;
        this.gff3Writer = gff3Writer;
        this.dao = dao;
    }

    public SequenceRegionWatcher( int mapKey, Gff3ColumnWriter gff3Writer, RgdGff3Dao dao, String chrPrefix ) {
        this.mapKey = mapKey;
        this.gff3Writer = gff3Writer;
        this.dao = dao;
        this.chromosomePrefix = chrPrefix;
    }

    public void init( int mapKey, Gff3ColumnWriter gff3Writer, RgdGff3Dao dao ) {
        this.mapKey = mapKey;
        this.gff3Writer = gff3Writer;
        this.dao = dao;

        encounteredChromosomes.clear();
    }

    synchronized public boolean emit( String chr ) {

        if( !encounteredChromosomes.contains(chr) ) {
            encounteredChromosomes.add(chr);
            int chrSize = 0;

            try {
                chrSize = dao.getChromosomeSize(mapKey, chr);
            } catch( Exception e ) {
                throw new RuntimeException(e);
            }

            if( chrSize>0 ) {
                gff3Writer.print("##sequence-region " + (chr.length()>2 ? chr : chromosomePrefix+chr) + " 1 " + chrSize + "\n");
                return true;
            }
        }
        return false;
    }
}
