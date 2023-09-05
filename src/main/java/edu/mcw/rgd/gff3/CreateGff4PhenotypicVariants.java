package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.impl.RgdVariantDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontologyx.TermWithStats;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;

import java.net.URLEncoder;
import java.util.Map;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class CreateGff4PhenotypicVariants {

    private RgdGff3Dao dao = new RgdGff3Dao();
    private String outDir;
    private List<Integer> processedMapKeys;

    /**
     * load the species list and assemblies from properties/AppConfigure.xml
     */
    public void run() throws Exception {

        for( int mapKey: getProcessedMapKeys() ) {

            String assemblyDir = Manager.getInstance().getAssemblies().get(mapKey);
            if( assemblyDir==null ) {
                break;
            }

            CreateInfo info = new CreateInfo();

            int speciesTypeKey = MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey();

            info.setMapKey( mapKey );
            info.setToDir( assemblyDir + "/" + getOutDir() );
            info.setSpeciesTypeKey( speciesTypeKey );
            info.setCompressMode( Gff3ColumnWriter.COMPRESS_MODE_BGZIP );

            creategff4QTL(info);
        }
    }

    /**
     * create Gff3 for Rat, Human and Mouse.. only for Rat it will print out related strains
     * @throws Exception
     */
    public void creategff4QTL(CreateInfo info) throws Exception {

        String speciesName = SpeciesType.getCommonName(info.getSpeciesTypeKey());

        String ucscId = Gff3Utils.getAssemblySymbol(info.getMapKey());
        String refseqId = MapManager.getInstance().getMap(info.getMapKey()).getRefSeqAssemblyName();
        String fileName = info.getToDir() + "/" + speciesName + " " + refseqId+" ("+ucscId+")";

        System.out.println("START Phenotypic Alleles Variant GFF3 Generator for  "+speciesName+"  MAP_KEY="+info.getMapKey()+"  ASSEMBLY "+ MapManager.getInstance().getMap(info.getMapKey()).getName());
        System.out.println("========================");

        String gffFile = fileName + " Phenotypic Alleles Variants.gff3";

        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(gffFile, false, info.getCompressMode());

        SequenceRegionWatcher sequenceRegionWatcher = new SequenceRegionWatcher(info.getMapKey(), gff3Writer, dao);

        RgdVariantDAO vdao = new RgdVariantDAO();
        List<RgdVariant> variants = vdao.getVariantsForSpecies(info.getSpeciesTypeKey());

        int dataLines = 0;

        for( RgdVariant var: variants ){

            TermWithStats soTerm = dao.getTerm(var.getType());

            List<Alias> aliases = dao.getAliases(var.getRgdId());
            String aliasesStr = Utils.concatenate(",", aliases, "getValue");

            String varNuc = Utils.NVL(var.getVarNuc(), "-");
            String refNuc = Utils.NVL(var.getRefNuc(), "-");

            List<MapData> mds = dao.getMapData(var.getRgdId(), info.getMapKey());

            for( MapData md: mds ) {

                sequenceRegionWatcher.emit(md.getChromosome());

                gff3Writer.writeFirst8Columns(md.getChromosome(), "RGD", soTerm.getTerm(), md.getStartPos(), md.getStopPos(), ".", md.getStrand(),".");

                Map<String, String> attributesHashMap = new HashMap<String, String>();

                String id = getUniqueId(var.getRgdId());
                attributesHashMap.put("ID", id);
                attributesHashMap.put("Name", var.getName());

                String aliasField = "RGD:"+var.getRgdId();
                if( !Utils.isStringEmpty(aliasesStr) ) {
                    aliasField += "," + aliasesStr;
                }
                attributesHashMap.put("Alias", aliasField);
                attributesHashMap.put("Dbxref","RGD:"+var.getRgdId());
                attributesHashMap.put("Ontology_term", soTerm.getAccId());

                StringBuffer associatedStrains = new StringBuffer();
                String associatedGenes = getAssociatedGenes(var.getRgdId(), associatedStrains);
                if( associatedGenes!=null ) {
                    attributesHashMap.put("associatedGenes", associatedGenes);
                }
                if( associatedStrains.length()>0 ) {
                    attributesHashMap.put("associatedStrains", associatedStrains.toString());
                }

                attributesHashMap.put("refNuc", refNuc);
                attributesHashMap.put("varNuc", varNuc);

                gff3Writer.writeAttributes4Gff3(attributesHashMap);

                dataLines++;
            }
        }

        System.out.println("Number of phenotypic allele variants processed:"+ variants.size());
        System.out.println("Number of gff data lines written:"+ dataLines);

        gff3Writer.close();
        gff3Writer.sortInMemory();
    }

    String getUniqueId(Integer id) {

        AtomicInteger i = idMap.putIfAbsent(id, new AtomicInteger(0));
        if( i==null ) {
            i = idMap.get(id);
        }
        int cnt = i.incrementAndGet();

        if( cnt==1 )
            return id+"";
        return id+"_"+cnt;
    }
    Map<Integer, AtomicInteger> idMap = new HashMap<>();

    String getAssociatedGenes(int rgdId, StringBuffer associatedStrains) throws Exception {

        List<Gene> geneList = dao.getAssociatedGenes(rgdId);
        if( geneList.isEmpty() ) {
            return null;
        }
        List<String> geneStrList = new ArrayList<>();
        for( Gene g: geneList ) {

            List<Strain> strains = dao.getAssociatedStrains(g.getRgdId());
            if( !strains.isEmpty() ) {

                List<String> strainStrList = new ArrayList<>();
                for( Strain s: strains ) {
                    strainStrList.add(s.getSymbol()+" (RGD:"+s.getRgdId()+")");
                }
                associatedStrains.append( Utils.concatenate(strainStrList, ",") );
            }

            String geneStr = g.getSymbol()+" (RGD:"+g.getRgdId()+")";
            geneStrList.add(geneStr);
        }
        return Utils.concatenate(geneStrList, ",");
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
