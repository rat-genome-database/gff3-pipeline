package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.impl.AssociationDAO;
import edu.mcw.rgd.dao.impl.GeneDAO;
import edu.mcw.rgd.dao.impl.GenomicElementDAO;
import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.text.DecimalFormat;
import java.util.*;
import java.util.Map;

/**
 * @author pjayaraman
 * @since 8/26/12
 */
public class CreatePromoters4Gene {
    Logger log = LogManager.getLogger("promoters");

    GenomicElementDAO genEleDao = new GenomicElementDAO();
    MapDAO mdao = new MapDAO();
    private String toDir="";

    public String getToDir() {
        return toDir;
    }

    public void setToDir(String toDir) {
        this.toDir = toDir;
    }

    public void createGenomicElements(boolean compress) throws Exception{

        long time0 = System.currentTimeMillis();

        Map<Integer, Integer> rowsWritten = new HashMap<>();

        AssociationDAO assDao = new AssociationDAO();
        GeneDAO gDao = new GeneDAO();

        List<GenomicElement> listOfGenElements = genEleDao.getActiveElements(RgdId.OBJECT_KEY_PROMOTERS);
        Collections.shuffle(listOfGenElements);

        listOfGenElements.parallelStream().forEach( ge -> {

            try {
                int rgdid = ge.getRgdId();
                List<MapData> mapObjects = mdao.getMapData(rgdid);
                for (MapData mdObject : mapObjects) {

                    Gff3ColumnWriter gff3Writer = getWriter(mdObject.getMapKey(), compress);

                    String strand = Utils.NVL(mdObject.getStrand(), ".");

                    Map<String, String> attributeHashMap = new HashMap<String, String>();
                    attributeHashMap.put("ID", ge.getRgdId() + "_" + mdObject.getStartPos() + "_" + mdObject.getStopPos());

                    attributeHashMap.put("Name", Utils.NVL(ge.getName(), ge.getSymbol()));

                    List<Association> rgdIdsList = assDao.getAssociationsForMasterRgdId(rgdid, "promoter_to_gene");
                    String promoterGeneAssoc = "";
                    for (Association assObj : rgdIdsList) {
                        Gene geneObj = gDao.getGene(assObj.getDetailRgdId());
                        promoterGeneAssoc += geneObj.getSymbol() + "_" + geneObj.getRgdId() + ",";
                    }
                    if (promoterGeneAssoc.length() > 0) {
                        String geneAssoc = promoterGeneAssoc.substring(0, (promoterGeneAssoc.length() - 1));
                        attributeHashMap.put("promoterGeneRlt", geneAssoc);
                    } else {
                        attributeHashMap.put("promoterGeneRlt", "NA");
                    }

                    String geName = "NA";
                    String geDesc = "NA";
                    if ((ge.getName() != null) && (!(ge.getName().equals("")))) {
                        geName = ge.getName().replaceAll(";", " ");
                    }
                    if ((ge.getDescription() != null) && (!(ge.getDescription().equals("")))) {
                        geDesc = ge.getDescription().replaceAll(";", "-");
                    }
                    attributeHashMap.put("Alias", +ge.getRgdId() + "," + geName + "," + geDesc);

                    String geNotes = Utils.defaultString(ge.getNotes());
                    if( !Utils.isStringEmpty(geNotes) ) {
                        String notes = geNotes.replaceAll(";", ",").replaceAll("=", "->");
                        attributeHashMap.put("Note", notes);
                    }

                    String geObjectType = "NA";
                    if ((ge.getObjectType() != null) && (!(ge.getObjectType().equals("")))) {
                        geObjectType = ge.getObjectType();
                    }
                    attributeHashMap.put("objectType", geObjectType);


                    Set<String> tissueArr = new TreeSet<String>();
                    Set<String> transcriptSet = new TreeSet<String>();

                    String chipSeqDensity = null, transcripts = null, expMethods = null, regul = null, tissues = null;

                    List<ExpressionData> promExprData = genEleDao.getExpressionData(rgdid);

                    for (ExpressionData expd : promExprData) {

                        String tissue = Utils.defaultString(expd.getTissue());
                        if( !Utils.isStringEmpty(tissue) ) {
                            if( !tissueArr.contains(tissue) ) {
                                tissueArr.add(tissue);
                            }
                        }

                        String transcriptList = Utils.defaultString(expd.getTranscripts());
                        String[] transcriptArr = transcriptList.split("[\\,]");
                        for( String tr: transcriptArr ) {
                            if( !Utils.isStringEmpty(tr) ) {
                                if( !transcriptSet.contains(tr) ) {
                                    transcriptSet.add(tr);
                                }
                            }
                        }

                        DecimalFormat df = new DecimalFormat("#.####");
                        String modChipSeqReads = df.format(expd.getChipSeqReadDensity());

                        if( !(modChipSeqReads.equals("0.0") || modChipSeqReads.equals("0")) ) {
                            if( chipSeqDensity==null ) {
                                chipSeqDensity = modChipSeqReads;
                            } else {
                                chipSeqDensity += "," + modChipSeqReads;
                            }
                        }


                        if( !Utils.isStringEmpty(expd.getExperimentMethods()) ) {

                            String expMeStr = expd.getExperimentMethods().replaceAll(";", " ");
                            if( expMethods==null ) {
                                expMethods = expMeStr;
                            } else {
                                if( !expMethods.contains(expMeStr) ) {
                                    expMethods += "," + expMeStr;
                                }
                            }
                        }


                        if( !Utils.isStringEmpty(expd.getRegulation()) ) {
                            String regulStr = expd.getRegulation().replaceAll(";", " AND ");
                            if( regul==null ) {
                                regul = regulStr;
                            } else {
                                if( !regul.contains(regulStr) ) {
                                    regul += "," + regulStr;
                                }
                            }
                        }
                    }

                    if (tissueArr.size() > 0) {
                        tissues = Utils.concatenate(tissueArr, ",");
                    }
                    if( !transcriptSet.isEmpty() ) {
                        transcripts = Utils.concatenate(transcriptSet, ",");
                    }


                    if( tissues!=null ) {
                        attributeHashMap.put("tissues", tissues);
                    }
                    if( transcripts!=null ) {
                        attributeHashMap.put("transcripts", transcripts);
                    }
                    if( chipSeqDensity!=null ) {
                        attributeHashMap.put("chipSeqDensityValues", chipSeqDensity);
                    }
                    if( expMethods!=null ) {
                        attributeHashMap.put("experimentMethods", expMethods);
                    }
                    if( regul!=null ) {
                        attributeHashMap.put("regulation", regul);
                    }

                    synchronized(gff3Writer) {
                        gff3Writer.writeFirst8Columns(mdObject.getChromosome(), ge.getSource() + "_RGD", ge.getSoAccId(), mdObject.getStartPos(),
                                mdObject.getStopPos(), ".", strand, ".");

                        gff3Writer.writeAttributes4Gff3(attributeHashMap);
                    }

                    incrementRowsWritten(mdObject.getMapKey(), rowsWritten);
                }
            } catch(Exception e) {
                throw new RuntimeException(e);
            }
        });


        log.info("PROMOTERS EXPORTED TO GFF3 FILES IN DIRECTORY "+getToDir()+"    elapsed "+Utils.formatElapsedTime(time0, System.currentTimeMillis()));
        log.info("===");
        for( int mapKey: rowsWritten.keySet() ) {
            Gff3ColumnWriter writer = getWriter(mapKey, compress);
            log.info("   "+writer.getFileName()+":   "+rowsWritten.get(mapKey)+" rows");
        }
        log.info("===");

        closeGff3Writers();
    }

    synchronized void incrementRowsWritten( int mapKey, Map<Integer, Integer> rowsWrittenMap ) {
        Integer rowsWritten = rowsWrittenMap.get(mapKey);
        if( rowsWritten==null ) {
            rowsWritten = 1;
        } else {
            rowsWritten++;
        }
        rowsWrittenMap.put(mapKey, rowsWritten);
    }

    synchronized Gff3ColumnWriter getWriter(int mapKey, boolean compress) throws Exception {

        Gff3ColumnWriter gff3Writer = gff3Writers.get(mapKey);
        if( gff3Writer==null ) {
            //String assemblySymbol = Gff3Utils.getAssemblySymbol(mapKey);
            edu.mcw.rgd.datamodel.Map map = MapManager.getInstance().getMap(mapKey);
            int speciesTypeKey = map.getSpeciesTypeKey();
            String species = SpeciesType.getCommonName(speciesTypeKey);

            gff3Writer = new Gff3ColumnWriter(getToDir()+species+"/"+map.getName()+"_promoters.gff3", false, compress);
            gff3Writer.print("# RAT GENOME DATABASE (https://rgd.mcw.edu/)\n");
            gff3Writer.print("# Species: "+ species+"\n");
            gff3Writer.print("# Assembly: "+ map.getName()+"\n");
            gff3Writer.print("# Primary Contact: mtutaj@mcw.edu\n");
            gff3Writer.print("# Generated: "+new Date()+"\n");
            gff3Writers.put(mapKey, gff3Writer);
        }
        return gff3Writer;
    }

    void closeGff3Writers() {
        for( Gff3ColumnWriter gff3Writer: gff3Writers.values() ) {
            gff3Writer.close();
        }
    }
    Map<Integer, Gff3ColumnWriter> gff3Writers = new HashMap<>();
}
