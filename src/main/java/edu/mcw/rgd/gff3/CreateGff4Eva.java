package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


/**
 * Created by llamers on 3/5/2020.
 */
public class CreateGff4Eva {

    private RgdGff3Dao dao = new RgdGff3Dao();
    Logger log = Logger.getLogger("evas");
    private List<String> processedAssemblies;

    public void run() throws Exception{
        try {
            log.info("EVA GFF3 GENERATOR -- all species");
            log.info(dao.getConnectionInfo()+"\n");

            long timeStart = System.currentTimeMillis();

            getProcessedAssemblies().parallelStream().forEach(assemblyInfo ->{
                try{
                    CreateInfo info = new CreateInfo();
                    info.parseFromString(assemblyInfo);
                    run(info);
                }
                catch(Exception e){
                    throw new RuntimeException(e);
                }
            });

            log.info("");
            log.info("OK elapsed "+Utils.formatElapsedTime(timeStart, System.currentTimeMillis()));
            log.info("");
        }
        catch (Exception e) {
            Utils.printStackTrace(e, log);
            throw e;
        }
    }

    public void run(CreateInfo info) throws Exception {
        String species = SpeciesType.getCommonName(info.getSpeciesTypeKey());
        int mapKey = info.getMapKey();
        String assemblyName = MapManager.getInstance().getMap(mapKey).getName();
        Gff3ColumnWriter gff3Writer = null;
        if(gff3Writer==null) {
            String gffFile = info.getToDir() + "EVA_" + assemblyName+ ".gff3";
            gff3Writer = new Gff3ColumnWriter(gffFile, false, info.isCompress());
        }
        List<String> chromosomes = getChromosomes(info.getMapKey());

        int dataLinesWritten = 0;
        for(String chr : chromosomes) {
            List<Eva> data = dao.getEvaObjectsbyKeyandChrom(mapKey,chr);
            log.debug(" "+assemblyName+": data lines for Eva in chrom "+chr+": "+data.size());

            if(data.size()==0)
                continue;

            String prevVarNuc = "";

            for(int i = 0; i < data.size(); i++)
            {
                if(i+1==data.size()) {
                    gff3Writer.writeFirst8Columns(data.get(i).getChromosome(), "EVA", "SNP", data.get(i).getPos(), data.get(i).getPos(), ".", ".", ".");
                    HashMap<String, String> attributes = new HashMap<>();
                    attributes.put("ID", Integer.toString(data.get(i).getEvaId()));
                    attributes.put("Name", data.get(i).getRsId());
                    attributes.put("Alias", data.get(i).getRsId());
                    prevVarNuc = "/"+data.get(i).getVarNuc();
                    attributes.put("allele", data.get(i).getRefNuc() + prevVarNuc);
                    gff3Writer.writeAttributes4Gff3(attributes);
                    dataLinesWritten++;
                }
                else { // check +1 rs id, if different add current data to file else add data to allele
                    String nextRsId = data.get(i+1).getRsId();
                    if(data.get(i).getRsId().equals(nextRsId)) {
                        prevVarNuc += "/" + data.get(i).getVarNuc();
                    }
                    else {

                        gff3Writer.writeFirst8Columns(data.get(i).getChromosome(), "EVA", "SNP", data.get(i).getPos(), data.get(i).getPos(), ".", ".", ".");
                        HashMap<String, String> attributes = new HashMap<>();
                        attributes.put("ID", Integer.toString(data.get(i).getEvaId()));
                        attributes.put("Name", data.get(i).getRsId());
                        attributes.put("Alias", data.get(i).getRsId());
                        String notHere = "-";
                        if(data.get(i).getRefNuc().equals(null)) {
                            prevVarNuc += "/" + data.get(i).getVarNuc();
                            attributes.put("allele", notHere+prevVarNuc);
                        }
                        else if(data.get(i).getVarNuc().equals(null)) {
                            attributes.put("allele",data.get(i).getRefNuc()+"/"+notHere);
                        }
                        else {
                            prevVarNuc += "/" + data.get(i).getVarNuc();
                            attributes.put("allele", data.get(i).getRefNuc() + prevVarNuc);
                        }
                        prevVarNuc = "";
                        gff3Writer.writeAttributes4Gff3(attributes);
                        dataLinesWritten++;
                    }
                }
            }
        }
        if(gff3Writer!=null)
            gff3Writer.close();

        synchronized( this.getClass() ) {
            log.info(species+", MAP_KEY="+info.getMapKey()+" ("+ assemblyName+")   -- data lines: "+Utils.formatThousands(dataLinesWritten));
        }
    }

    public void setProcessedAssemblies(List<String> processedAssemblies) {
        this.processedAssemblies = processedAssemblies;
    }

    public List<String> getProcessedAssemblies() {
        return processedAssemblies;
    }

    List<String> getChromosomes(int mapKey) throws Exception {

        // truncate version numbers from scaffold accessions
        List<String> result = new ArrayList<>();
        List<Chromosome> chromosomes = dao.getChromosomes(mapKey);
        for( Chromosome chr: chromosomes ) {
            String c = chr.getChromosome();
            result.add(c);
        }
        return result;
    }
}
