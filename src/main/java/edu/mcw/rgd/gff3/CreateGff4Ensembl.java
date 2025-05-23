package edu.mcw.rgd.gff3;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Map;

public class CreateGff4Ensembl {

    private String outDir;
    private Map<String,String> mappings;
    private Map<String,Integer> mappings2;

    public void run() throws Exception {

        for( Map.Entry<String, String> entry: getMappings().entrySet() ) {

            String srcFileName = "data/Ensembl/"+entry.getKey();

            int mapKey = getMappings2().get(entry.getKey());
            String assemblyDir = Manager.getInstance().getAssemblies().get(mapKey);
            String outDir2 = assemblyDir+"/"+Gff3Utils.getAssemblyDirStandardized(mapKey)+"/"+getOutDir();
            new File(outDir2).mkdirs();
            String outPath = outDir2+"/"+entry.getValue();

            File inFile = new File(srcFileName);
            File outFile = new File(outPath);
            Files.copy( inFile.toPath(), outFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
        }
    }

    public String getOutDir() {
        return outDir;
    }

    public void setOutDir(String outDir) {
        this.outDir = outDir;
    }

    public Map<String, String> getMappings() {
        return mappings;
    }

    public void setMappings(Map<String, String> mappings) {
        this.mappings = mappings;
    }

    public Map<String, Integer> getMappings2() {
        return mappings2;
    }

    public void setMappings2(Map<String, Integer> mappings2) {
        this.mappings2 = mappings2;
    }
}
