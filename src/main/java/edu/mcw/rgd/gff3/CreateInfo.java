package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;



/**
 * configuration info used to run a gf3 create job
 */
public class CreateInfo {
    private String toDir;
    private int mapKey;
    private int speciesTypeKey;
    private int compressMode;
    private String assemblySymbol;


    // properties available after calling 'getCanonicalFileName()'
    public String ucscId;
    public String refseqId;
    public String speciesName = SpeciesType.getCommonName(getSpeciesTypeKey());

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
                    case "species": break; // ignore it: we will derive speciesTypeKey from mapKey
                    case "mapKey": mapKey = Integer.parseInt(value); break;
                    case "toDir": toDir = value; break;
                    case "assemblySymbol": assemblySymbol = value; break;
                    case "compress": {
                        if( value.equals("bgzip") || value.equals("bgz") || value.equals("bgzf") ) {
                            compressMode = Gff3ColumnWriter.COMPRESS_MODE_BGZIP;
                        }
                        else if( value.equals("gzip") || value.equals("gz") || value.equals("yes") ) {
                            compressMode = Gff3ColumnWriter.COMPRESS_MODE_ZIP;
                        }
                        else {
                            compressMode = Gff3ColumnWriter.COMPRESS_MODE_NONE;
                        }
                    } break;
                    default: throw new Exception("unknown field: "+name);
                }
            }
        }

        if( mapKey>0 ) {
            speciesTypeKey = MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey();
        }
    }

    // creates a file name based on mapKey, RefSeq assembly name and/or UCSC assembly name
    public String getCanonicalFileName() throws Exception {

        ucscId = Utils.NVL( Gff3Utils.getAssemblySymbol(getMapKey()), "" );
        refseqId = MapManager.getInstance().getMap(getMapKey()).getRefSeqAssemblyName();
        speciesName = SpeciesType.getCommonName(getSpeciesTypeKey());

        String fileName = getToDir() + "/" + speciesName + " " + refseqId;
        if( !ucscId.isEmpty()  &&  !ucscId.equalsIgnoreCase(refseqId) ) {
            fileName += " ("+ucscId+")";
        }
        return fileName;
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

    public int getCompressMode() {
        return compressMode;
    }

    public void setCompressMode(int compressMode) {
        this.compressMode = compressMode;
    }

    public String getAssemblySymbol() {
        return assemblySymbol;
    }

    public void setAssemblySymbol(String assemblySymbol) {
        this.assemblySymbol = assemblySymbol;
    }

    public void dump() {
        String compressModeStr = "none";
        if( compressMode==Gff3ColumnWriter.COMPRESS_MODE_ZIP ) {
            compressModeStr = "zip";
        }
        else if( compressMode==Gff3ColumnWriter.COMPRESS_MODE_BGZIP ) {
            compressModeStr = "bgzip";
        }

        System.out.println("toDir = "+toDir);
        System.out.println("mapKey = "+mapKey);
        System.out.println("speciesTypeKey = "+speciesTypeKey);
        System.out.println("compressMode = "+compressModeStr);
        if( assemblySymbol!=null ) {
            System.out.println("assemblySymbol = " + assemblySymbol);
        }
    }
}
