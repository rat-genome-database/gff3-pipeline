package edu.mcw.rgd.gff3;

import edu.mcw.rgd.datamodel.Chromosome;
import edu.mcw.rgd.process.Utils;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class JBrowse2Aliases {

    public void run() throws Exception {

        RgdGff3Dao dao = new RgdGff3Dao();

        Map<Integer, String> assemblies = Manager.getInstance().getAssemblies();

        for( int mapKey: assemblies.keySet() ) {

            String assemblyName = Gff3Utils.getAssemblyDirStandardized(mapKey);

            String fileName = "aliases_"+assemblyName+".txt";

            System.out.println(fileName);

            BufferedWriter out = Utils.openWriter(fileName);

            List<Chromosome> chromosomes = dao.getChromosomes(mapKey);
            for( Chromosome c: chromosomes ) {

                ArrayList<String> aliases = new ArrayList<>();

                // chromosome specific: add prefix 'chr' and 'Chr'
                if( !c.getChromosome().startsWith("N") ) { // if not a scaffold

                    aliases.add("Chr"+c.getChromosome());
                    aliases.add("chr"+c.getChromosome());
                }

                aliases.add(c.getChromosome());

                String refSeqAcc = c.getRefseqId();
                if( !Utils.isStringEmpty(refSeqAcc) ) {
                    aliases.add(refSeqAcc);
                }
                String genbankAcc = c.getGenbankId();
                if( !Utils.isStringEmpty(genbankAcc) ) {
                    aliases.add(genbankAcc);
                }

                out.write( Utils.concatenate(aliases, "\t") );
                out.write("\n");
            }

            out.close();
        }
    }
}
