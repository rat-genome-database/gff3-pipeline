package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.SpeciesType;

/**
 * configuration info used to run an gf33 create job
 */
public class CreateInfo {
    private String toDir;
    private int mapKey;
    private int speciesTypeKey;
    private boolean compress;

    public void parseFromString(String s) throws Exception {
        // sample string:
        //    species:HUMAN mapKey:13 toDir:data/Qtl/Human/human36/ compress:yes
        String[] fields = s.split("[\\s]+");
        for( String field: fields ) {
            // every field has name and value separated by ':'
            int colonPos = field.indexOf(':');
            if( colonPos>0 ) {
                String name = field.substring(0, colonPos).trim();
                String value = field.substring(colonPos+1).trim();
                switch(name) {
                    case "species": speciesTypeKey = SpeciesType.parse(value); break;
                    case "mapKey": mapKey = Integer.parseInt(value); break;
                    case "toDir": toDir = value; break;
                    case "compress": compress = (value.equals("yes") ? true : false); break;
                    default: throw new Exception("unknown field: "+name);
                }
            }
        }
    }

    public String getToDir() {
        return toDir;
    }

    public void setToDir(String toDir) {
        this.toDir = toDir;
    }

    public int getMapKey() {
        return mapKey;
    }

    public void setMapKey(int mapKey) {
        this.mapKey = mapKey;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }

    public boolean isCompress() {
        return compress;
    }

    public void setCompress(boolean compress) {
        this.compress = compress;
    }

    public void dump() {
        System.out.println("toDir = "+toDir);
        System.out.println("mapKey = "+mapKey);
        System.out.println("speciesTypeKey = "+speciesTypeKey);
        System.out.println("compress = "+compress);
    }
}
