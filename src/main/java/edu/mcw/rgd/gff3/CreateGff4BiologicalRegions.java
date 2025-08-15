package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class CreateGff4BiologicalRegions {


    RgdGff3Dao dao = new RgdGff3Dao();
    Logger log = LogManager.getLogger("gene");

    private String outDir;
    private List<Integer> processedMapKeys;

    /**
     * load the species list and assemblies from properties/AppConfigure.xml
     */
    public void run() {

        processedMapKeys.parallelStream().forEach( mapKey -> {

            String assemblyDir = Manager.getInstance().getAssemblies().get(mapKey);
            if( assemblyDir==null ) {
                return;
            }

            int speciesTypeKey = 0;
            try {
                assemblyDir += "/" + Gff3Utils.getAssemblyDirStandardized(mapKey);
                speciesTypeKey = MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey();
            } catch( Exception e ) {
            }
            if( speciesTypeKey==0 ) {
                return;
            }

            CreateInfo info = new CreateInfo();
            info.setMapKey( mapKey );
            info.setToDir( assemblyDir + "/" + getOutDir() );
            info.setSpeciesTypeKey( speciesTypeKey );
            info.setCompressMode( Gff3ColumnWriter.COMPRESS_MODE_BGZIP );

            try {
                createGff3(info);
            } catch(Exception e) {
                throw new RuntimeException(e);
            }
        });
    }

    public void createGff3(CreateInfo info) throws Exception{

        StringBuffer msgBuf = new StringBuffer();

        long time0 = System.currentTimeMillis();

        String speciesName = SpeciesType.getCommonName(info.getSpeciesTypeKey());
        Map<String, AtomicInteger> idMap = new ConcurrentHashMap<>();

        String ucscId = Utils.NVL( Gff3Utils.getAssemblySymbol(info.getMapKey()), "" );
        String refseqId = MapManager.getInstance().getMap(info.getMapKey()).getRefSeqAssemblyName();

        if( ucscId.isEmpty() ) {
            msgBuf.append("Generate GFF3 file for " + speciesName + ", MAP_KEY=" + info.getMapKey() + "\n");
        } else {
            msgBuf.append("Generate GFF3 file for " + speciesName + ", MAP_KEY=" + info.getMapKey() + " (" + ucscId + ")\n");
        }
        msgBuf.append("    "+dao.getConnectionInfo()+"\n");

        CounterPool counters = new CounterPool();

        String headerInfo =
            "# RAT GENOME DATABASE (https://rgd.mcw.edu/)\n"+
            "# Species: "+ speciesName+"\n"+
            "# Assembly: "+ MapManager.getInstance().getMap(info.getMapKey()).getName()+"\n"+
            "# Primary Contact: mtutaj@mcw.edu\n"+
            "# Generated: "+new Date()+"\n";

        String fileName = info.getToDir() + "/" + speciesName + " " + refseqId;
        if( !ucscId.isEmpty()  &&  !ucscId.equalsIgnoreCase(refseqId) ) {
            fileName += " ("+ucscId+")";
        }

        Gff3ColumnWriter gff3BiologicalRegions = new Gff3ColumnWriter(fileName + " Biological Regions.gff3", info.getCompressMode());
        gff3BiologicalRegions.print(headerInfo);
        SequenceRegionWatcher sequenceRegionWatcher = new SequenceRegionWatcher(info.getMapKey(), gff3BiologicalRegions, dao);

        List<GenomicElement> activeElements = dao.getActiveBiologicalRegions(info.getSpeciesTypeKey());

        for( GenomicElement ge: activeElements ){

            int geRgdId = ge.getRgdId();

            List<MapData> geMap = getMapData(geRgdId, info.getMapKey(), counters);
            if( geMap.isEmpty() ) {
                //System.out.println("no map positions");
                continue;
            }
            counters.increment(" Biological Regions processed");

            String name = getName(ge);

            String annotDesc = ge.getDescription();
            if( Utils.isStringEmpty(annotDesc) ){
                annotDesc = null;
            }else{
                annotDesc = annotDesc.replaceAll(";"," AND ").replaceAll(",", " ");
            }

            for( MapData map : geMap ){

                String gType = ge.getObjectType();

                List<XdbId> xdbIds = dao.getXdbIds(ge.getRgdId());

                String aliasesStr = ge.getSymbol()+","+"RGD"+ge.getRgdId()+","+ge.getRgdId() + getHgncMgiIds(xdbIds);
                if( name != null ) {
                    aliasesStr += "," + name;
                }

                Map<String,String> attributesHashMap = new HashMap<>();

                String uniqueGeId = getUniqueId("RGD"+ge.getRgdId(), idMap);
                attributesHashMap.put("ID", uniqueGeId);
                attributesHashMap.put("Name", ge.getSymbol());
                if( name!=null ) {
                    attributesHashMap.put("fullName", name);
                }
                attributesHashMap.put("Alias", aliasesStr);
                attributesHashMap.put("regionType", gType.replaceAll("\\-","_"));
                attributesHashMap.put("species", speciesName);
                if( annotDesc!=null )
                    attributesHashMap.put("info", annotDesc);

                String extDbString = getXdbString(xdbIds, counters);
                if( !extDbString.isEmpty() ){
                    attributesHashMap.put("Dbxref",extDbString);
                }

                String attrStr = Gff3ColumnWriter.prepAttributes4Gff3(attributesHashMap);

                sequenceRegionWatcher.emit(map.getChromosome());

                gff3BiologicalRegions.writeFirst8Columns(map.getChromosome(), "RGD", "biological_region", map.getStartPos(), map.getStopPos(), ".", map.getStrand(), ".");
                gff3BiologicalRegions.print(attrStr);
            }
        }

        gff3BiologicalRegions.close();
        gff3BiologicalRegions.sortInMemory();

        dumpCounters(counters, ucscId.isEmpty() ? refseqId : ucscId, msgBuf);

        msgBuf.append("OK!  elapsed "+Utils.formatElapsedTime(time0, System.currentTimeMillis())+"\n");
        msgBuf.append("========\n");
        synchronized( this.getClass() ) {
            log.info(msgBuf);
        }
    }

    List<MapData> getMapData(int geneRgdId, int mapKey, CounterPool counters) throws Exception {

        List<MapData> geneMap = dao.getMapData(geneRgdId, mapKey);

        if( !geneMap.isEmpty() ) {
            if(geneMap.size()>1){
                counters.increment(" Genes with more than one map position");
            }
        }
        return geneMap;
    }

    void dumpCounters(CounterPool counters, String assembly, StringBuffer msgBuf) {

        String msg = counters.dumpAlphabetically();

        // prepend all lines with 'assembly'
        String[] lines = msg.split("[\\n]", -1);
        for( String line: lines ) {
            msgBuf.append(assembly).append("  ").append(line).append("\n");
        }
    }

    /**
     * returns a comma separated string of key:value pairs that are external identifiers
     * @param xdbList list of XdbId-s
     * @return xdb id string
     * @throws Exception
     */
    private String getXdbString(List<XdbId> xdbList, CounterPool counters) {

        Set<String> xdbIds = new TreeSet<>();
        for(XdbId externalId : xdbList ){
            if(externalId.getAccId()!=null){
                String xdbEntry;

                if(externalId.getXdbKeyAsString().contains(" ")){
                    xdbEntry = externalId.getXdbKeyAsString().replaceAll(" ","")+":";
                }else{
                    xdbEntry = externalId.getXdbKeyAsString()+":";
                }

                if(externalId.getAccId().contains(":")){
                    xdbEntry += externalId.getAccId().replaceAll(":", "");
                }else{
                    xdbEntry += externalId.getAccId();
                }

                xdbIds.add(xdbEntry);
            }
        }
        String extDbString = Utils.concatenate(xdbIds,",");

        if( !extDbString.isEmpty() ){
            if( extDbString.contains("NCBIGene") )
                counters.increment(" Biological Regions with NCBI geneIds");
        }else{
            counters.increment(" Biological Regions with NO XDB Ids");
        }

        return extDbString;
    }

    /**
     * get comma separated and comma prefixed list of MGI or HGNC ids, f.e. ',HGNC:2';
     * this string is to be directly appended to 'Aliases' field
     * @param xdbList list of XdbId-s
     * @return MGI/HGNC ids string
     * @throws Exception
     */
    private String getHgncMgiIds(List<XdbId> xdbList) {

        StringBuilder xdbIds = new StringBuilder();
        for(XdbId externalId : xdbList ){
            if(externalId.getXdbKey()==XdbId.XDB_KEY_MGD || externalId.getXdbKey()==XdbId.XDB_KEY_HGNC) {
                xdbIds.append(",").append(externalId.getAccId());
            }
        }
        return xdbIds.toString();
    }

    private String getName(GenomicElement ge) {

        String name = ge.getName();

        if( name!=null ){
            if(name.contains(",")){
                name = name.replaceAll(",", "");
            }
            if(name.contains(";")){
                name = name.replaceAll(";", "");
            }
        }
        return name;
    }

    String getUniqueId(String idBase, Map<String, AtomicInteger> idMap) {

        AtomicInteger i = idMap.putIfAbsent(idBase, new AtomicInteger(0));
        if( i==null ) {
            i = idMap.get(idBase);
        }
        int cnt = i.incrementAndGet();

        if( cnt==1 )
            return idBase;
        return idBase+"_"+cnt;
    }

    public String getOutDir() {
        return outDir;
    }

    public void setOutDir(String outDir) {
        this.outDir = outDir;
    }

    public List<Integer> getProcessedMapKeys() {
        return processedMapKeys;
    }

    public void setProcessedMapKeys(List<Integer> processedMapKeys) {
        this.processedMapKeys = processedMapKeys;
    }
}
