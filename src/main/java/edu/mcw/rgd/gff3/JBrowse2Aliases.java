package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.Chromosome;
import edu.mcw.rgd.process.Utils;

import java.io.BufferedWriter;
import java.util.List;
import java.util.Map;

public class JBrowse2Aliases {

    public void run() throws Exception {

        RgdGff3Dao dao = new RgdGff3Dao();

        Map<Integer, String> assemblies = Manager.getInstance().getAssemblies();

        for( Map.Entry<Integer,String> entry: assemblies.entrySet() ) {

            int mapKey = entry.getKey();
            String assemblyPath = entry.getValue();
            int slashPos = assemblyPath.lastIndexOf("/");
            String assemblyName = assemblyPath.substring(slashPos+1);

            String fileName = "aliases_"+assemblyName+".txt";

            System.out.println(fileName);

            BufferedWriter out = Utils.openWriter(fileName);

            List<Chromosome> chromosomes = dao.getChromosomes(mapKey);
            for( Chromosome c: chromosomes ) {

                out.write(c.getChromosome());

                // chromosome specific: add prefix 'chr' and 'Chr'
                if( !c.getChromosome().startsWith("N") ) { // if not a scaffold

                    out.write("\tchr"+c.getChromosome());
                    out.write("\tChr"+c.getChromosome());
                }

                String refSeqAcc = c.getRefseqId();
                if( !Utils.isStringEmpty(refSeqAcc) ) {
                    out.write("\t"+refSeqAcc);
                }
                String genbankAcc = c.getGenbankId();
                if( !Utils.isStringEmpty(genbankAcc) ) {
                    out.write("\t"+genbankAcc);
                }

                out.write("\n");
            }

            out.close();
        }
    }
}
