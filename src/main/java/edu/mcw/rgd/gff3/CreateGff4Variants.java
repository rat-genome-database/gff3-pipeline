package edu.mcw.rgd.gff3;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Map;

public class CreateGff4Variants {

    private String outDir1;
    private String outDir2;
    private String suffix1;
    private String suffix2;

    private Map<Integer,String> srcDataDir;

    public void run() throws IOException {

        // copy files for all assemblies
        for( Map.Entry<Integer, String> entry: getSrcDataDir().entrySet() ) {

            int mapKey = entry.getKey();
            String srcDir = entry.getValue();

            String assemblyDir = Manager.getInstance().getAssemblies().get(mapKey);

            String outDir1Name = assemblyDir+"/"+getOutDir1();
            new File(outDir1Name).mkdirs();

            String outDir2Name = assemblyDir+"/"+getOutDir2();
            new File(outDir2Name).mkdirs();

            File[] files = new File(srcDir).listFiles();
            if( files!=null )
            for( File f: files ) {
                if( !f.isFile() ) {
                    continue;
                }
                if( f.getName().endsWith( getSuffix1() )) {

                    File inFile = f;
                    File outFile = new File(outDir1Name+"/"+f.getName());
                    Files.copy( inFile.toPath(), outFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                }

                if( f.getName().endsWith( getSuffix2() )) {

                    File inFile = f;
                    String newName = f.getName().replace("_damaging","");
                    File outFile = new File(outDir2Name+"/"+newName);
                    Files.copy( inFile.toPath(), outFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                }
            }
        }
    }

    public String getOutDir1() {
        return outDir1;
    }

    public void setOutDir1(String outDir1) {
        this.outDir1 = outDir1;
    }

    public String getOutDir2() {
        return outDir2;
    }

    public void setOutDir2(String outDir2) {
        this.outDir2 = outDir2;
    }

    public String getSuffix1() {
        return suffix1;
    }

    public void setSuffix1(String suffix1) {
        this.suffix1 = suffix1;
    }

    public String getSuffix2() {
        return suffix2;
    }

    public void setSuffix2(String suffix2) {
        this.suffix2 = suffix2;
    }

    public Map<Integer, String> getSrcDataDir() {
        return srcDataDir;
    }

    public void setSrcDataDir(Map<Integer, String> srcDataDir) {
        this.srcDataDir = srcDataDir;
    }
}
