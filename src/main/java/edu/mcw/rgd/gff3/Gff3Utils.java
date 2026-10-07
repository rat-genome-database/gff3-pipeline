package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.Map;
import edu.mcw.rgd.process.mapping.MapManager;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author mtutaj
 * @since 2019-10-24
 */
public class Gff3Utils {

    /**
     * @param mapKeys configured map keys, in their configured order
     * @param speciesTypeKey 0 keeps every map key; otherwise only map keys of assemblies of this species are kept
     */
    public static List<Integer> filterMapKeysBySpecies(Collection<Integer> mapKeys, int speciesTypeKey) throws Exception {
        List<Integer> result = new ArrayList<>();
        for( int mapKey: mapKeys ) {
            if( speciesTypeKey<=0 || MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey()==speciesTypeKey ) {
                result.add(mapKey);
            }
        }
        return result;
    }

    /// return human friendly assembly symbol
    synchronized static public String getAssemblySymbol(int mapKey) throws Exception {

        // first, return UCSC assembly symbol, if available
        Map map = MapManager.getInstance().getMap(mapKey);
        if( map == null ) {
            return "";
        }
        String symbol = map.getUcscAssemblyId();
        if( symbol != null ) {
            return symbol;
        }

        // if not, create UCSC-like assembly symbol, by taking the assembly name,
        // making lowercase the first letter and getting rid of version
        // f.e. 'ChiLan1.0' becomes 'chiLan1'
        int dotPos = map.getName().indexOf('.');
        if( dotPos > 0 ) {
            symbol = Character.toLowerCase(map.getName().charAt(0)) + map.getName().substring(1, dotPos);
        }
        return symbol;
    }

    synchronized static public String getAssemblyDirStandardized( int mapKey ) throws Exception {
        String assemblyName = MapManager.getInstance().getMap(mapKey).getRefSeqAssemblyName();
        String stdAssemblyName = assemblyName.replace(" ", "");
        return stdAssemblyName;
    }
}
