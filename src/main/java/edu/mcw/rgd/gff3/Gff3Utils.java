package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.Map;
import edu.mcw.rgd.process.mapping.MapManager;

/**
 * @author mtutaj
 * @since 2019-10-24
 */
public class Gff3Utils {

    /// return human friendly assembly symbol
    synchronized static public String getAssemblySymbol(int mapKey) throws Exception {

        // first, return UCSC assembly symbol, if available
        Map map = MapManager.getInstance().getMap(mapKey);
        String symbol = map.getUcscAssemblyId();
        if( symbol!=null ) {
            return symbol;
        }

        // if not, create UCSC-like assembly symbol, by taking the assembly name,
        // making lowercase the first letter and getting rid of version
        // f.e. 'ChiLan1.0' becomes 'chiLan1'
        int dotPos = map.getName().indexOf('.');
        symbol = Character.toLowerCase(map.getName().charAt(0)) + map.getName().substring(1, dotPos);
        return symbol;
    }
}
