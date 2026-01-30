package edu.mcw.rgd.gff3;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Map;

public class JBrowse2PrepEva {

    private String outDir;
    private String outPrefix;
    private Map<Integer, String> srcFiles;

    public void run() throws IOException {

        for( Map.Entry<Integer, String> entry: getSrcFiles().entrySet() ) {

            int mapKey = entry.getKey();
            String srcFileName = entry.getValue();

            String assemblyDir = Manager.getInstance().getAssemblies().get(mapKey);
            String outDir2 = assemblyDir+"/"+getOutDir()+"/";
            new File(outDir2).mkdirs();
            String outPath = outDir2+"/"+getOutPrefix()+".gff3.gz";

            File inFile = new File(srcFileName);
            File outFile = new File(outPath);
            System.out.println("copying "+inFile.getAbsolutePath()+" to "+outFile.getAbsolutePath());
            Files.copy( inFile.toPath(), outFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

            Gff3ColumnWriter.sortInMemory(outPath, Gff3ColumnWriter.COMPRESS_MODE_BGZIP);
        }
    }

    public String getOutDir() {
        return outDir;
    }

    public void setOutDir(String outDir) {
        this.outDir = outDir;
    }

    public String getOutPrefix() {
        return outPrefix;
    }

    public void setOutPrefix(String outPrefix) {
        this.outPrefix = outPrefix;
    }

    public Map<Integer, String> getSrcFiles() {
        return srcFiles;
    }

    public void setSrcFiles(Map<Integer, String> srcFiles) {
        this.srcFiles = srcFiles;
    }
}
