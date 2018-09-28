package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.impl.AssociationDAO;
import edu.mcw.rgd.dao.impl.GeneDAO;
import edu.mcw.rgd.dao.impl.GenomicElementDAO;
import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.datamodel.*;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: pjayaraman
 * Date: 8/26/12
 * Time: 4:41 PM
 */
public class CreatePromoters4Gene {
    GenomicElementDAO genEleDao = new GenomicElementDAO();
    MapDAO mdao = new MapDAO();
    private int objectKey;
    int mapKey;
    String toFile="";
    int speciesTypekey;

    public int getMapKey() {
        return mapKey;
    }

    public void setMapKey(int mapKey) {
        this.mapKey = mapKey;
    }

    public String getToFile() {
        return toFile;
    }

    public void setToFile(String toFile) {
        this.toFile = toFile;
    }

    public int getObjectKey() {
        return objectKey;
    }

    public void setObjectKey(int objectKey) {
        this.objectKey = objectKey;
    }

    public int getSpeciesTypekey() {
        return speciesTypekey;
    }

    public void setSpeciesTypekey(int speciesTypekey) {
        this.speciesTypekey = speciesTypekey;
    }

    public void createGenomicElements(boolean compress) throws Exception{

        Gff3ColumnWriter gff3Writer = new Gff3ColumnWriter(getToFile()+ SpeciesType.getCommonName(speciesTypekey)
                +"_RGD"+ RgdId.getObjectTypeName(getObjectKey())+".gff3", false, compress);

        AssociationDAO assDao = new AssociationDAO();
        GeneDAO gDao = new GeneDAO();

        List<GenomicElement> listOfGenElements = genEleDao.getActiveElements(getObjectKey());
        for(GenomicElement ge:listOfGenElements){
            if(ge.getSpeciesTypeKey()==speciesTypekey){
                int rgdid = ge.getRgdId();
                List<MapData> mapObjects = mdao.getMapData(rgdid, mapKey);
                for(MapData mdObject : mapObjects){
                    String strand = ".";
                    if(mdObject.getStrand()!=null){
                        strand = mdObject.getStrand();
                    }else{
                        strand = ".";
                    }

                    gff3Writer.writeFirst8Columns(mdObject.getChromosome(), ge.getSource()+"_RGD", ge.getSoAccId(), mdObject.getStartPos(),
                            mdObject.getStopPos(), ".", strand, ".");

                    Map<String, String> attributeHashMap = new HashMap<String, String>();
                    attributeHashMap.put("ID", ge.getRgdId()+"_"+mdObject.getStartPos()+"_"+mdObject.getStopPos());
                    attributeHashMap.put("Name", ge.getSymbol());

                    List<Association> rgdIdsList = assDao.getAssociationsForMasterRgdId(rgdid, "promoter_to_gene");
                    String promoterGeneAssoc = "";
                    for(Association assObj: rgdIdsList){
                        Gene geneObj = gDao.getGene(assObj.getDetailRgdId());
                        promoterGeneAssoc+=geneObj.getSymbol()+"_"+geneObj.getRgdId()+",";
                    }
                    if(promoterGeneAssoc.length()>0){
                        String geneAssoc = promoterGeneAssoc.substring(0, (promoterGeneAssoc.length()-1));
                        attributeHashMap.put("promoterGeneRlt", geneAssoc);
                    }else{
                        attributeHashMap.put("promoterGeneRlt", "NA");
                    }

                    String geName="NA"; String geDesc="NA";
                    if((ge.getName()!=null) && (!(ge.getName().equals("")))){
                        geName=ge.getName().replaceAll(";"," ");
                    }
                    if((ge.getDescription()!=null) && (!(ge.getDescription().equals("")))){
                        geDesc=ge.getDescription().replaceAll(";","-");
                    }
                    attributeHashMap.put("Alias", +ge.getRgdId()+","+geName+","+geDesc);

                    String geNotes="NA";
                    if((ge.getNotes()!=null) && (!(ge.getNotes().equals("")))){
                        String notes = ge.getNotes().replaceAll(";",",");
                        geNotes = notes.replaceAll("=", "->");
                    }
                    attributeHashMap.put("Note", geNotes);

                    String geObjectType="NA";
                    if((ge.getObjectType()!=null) && (!(ge.getObjectType().equals("")))){
                        geObjectType=ge.getObjectType();
                    }
                    attributeHashMap.put("objectType", geObjectType);


                    ArrayList<String> tissueArr = new ArrayList<String>();

                    String tissues="", chipSeqDensity="", transcripts="", expMethods="";
                    String regul="", series="";
                    String modTissues="NA", modTranscripts="NA", modChipSeqDensity="0.0";
                    String modExpMethods="NA", modRegul="NA";

                    List<ExpressionData> promExprData = genEleDao.getExpressionData(rgdid);


                    for(ExpressionData expd : promExprData){

                        if((expd.getTissue()!=null)&&(!expd.getTissue().equals(""))){
                            if(!(tissueArr.contains(expd.getTissue()))){
                                //tissues+=expd.getTissue()+",";
                                tissueArr.add(expd.getTissue());
                            }
                        }else{
                            if(!(tissues.contains("NA"))){
                                //tissues+="NA,";
                                tissueArr.add("NA");
                            }
                        }


                        if((expd.getTranscripts()!=null)&&(!expd.getTranscripts().equals(""))){
                            if(!(transcripts.contains(expd.getTranscripts()))){
                                transcripts+=expd.getTranscripts()+",";
                            }
                        }else{
                            if(!(transcripts.contains("NA"))){
                                transcripts+="NA,";
                            }
                        }

                        DecimalFormat df = new DecimalFormat("#.####");
                        String modChipSeqReads = df.format(expd.getChipSeqReadDensity());

                        if((modChipSeqReads!=null)&&(!modChipSeqReads.equals(""))){
                            if(!(chipSeqDensity.contains(modChipSeqReads))){
                                chipSeqDensity+=modChipSeqReads+",";
                            }

                            if(!(series.contains(modChipSeqReads))){
                                series=modChipSeqReads;
                            }
                        }else{
                            if(!(chipSeqDensity.contains("0.0"))){
                                chipSeqDensity+="0.0,";
                            }

                            if(!(series.contains("0.0"))){
                                series="0.0";
                            }
                        }


                        if((expd.getExperimentMethods()!=null)&&(!expd.getExperimentMethods().equals(""))){
                            if(!(expMethods.contains(expd.getExperimentMethods()))){
                                expMethods+=expd.getExperimentMethods().replaceAll(";", " AND ")+",";
                            }
                        }else{
                            if(!(expMethods.contains("NA"))){
                                expMethods+="NA,";
                            }
                        }


                        if((expd.getRegulation()!=null)&&(!expd.getRegulation().equals(""))){
                            if(!(regul.contains(expd.getRegulation()))){
                                regul+=expd.getRegulation().replaceAll(";", " AND ")+",";
                            }
                        }else{
                            if(!(regul.contains("NA"))){
                                regul+="NA,";
                            }
                        }

                    }

                    if(tissueArr.size()>0){
                        for(String tiss : tissueArr){
                            tissues+=tiss+",";
                        }
                        if(tissues.endsWith(",")){
                            modTissues = tissues.substring(0, (tissues.length()-1));
                        }
                    }else{
                        modTissues = "NA";
                    }


                    if(transcripts.endsWith(",")){
                        modTranscripts = transcripts.substring(0, (transcripts.length()-1));
                    }
                    if(chipSeqDensity.endsWith(",")){
                        modChipSeqDensity = chipSeqDensity.substring(0, (chipSeqDensity.length()-1));
                    }
                    if(expMethods.endsWith(",")){
                        modExpMethods = expMethods.substring(0, (expMethods.length()-1));
                    }
                    if(regul.endsWith(",")){
                        modRegul = regul.substring(0, (regul.length()-1));
                    }


                    attributeHashMap.put("tissues", modTissues);
                    attributeHashMap.put("transcripts", modTranscripts);
                    attributeHashMap.put("chipSeqDensityValues", modChipSeqDensity);
                    attributeHashMap.put("experimentMethods", modExpMethods);
                    attributeHashMap.put("regulation", modRegul);
                    if( !series.isEmpty() )
                        attributeHashMap.put("series", series);

                    gff3Writer.writeAttributes4Gff3(attributeHashMap);
                }
            }
        }
        gff3Writer.close();
    }
}
