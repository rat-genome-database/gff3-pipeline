package edu.mcw.rgd.gff3;

import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 * @author mtutaj
 * @since 11/20/2018
 * Prepares a gff3 file downloaded from Ensembl FTP site to be loaded by JBrowse scripts as an Ensembl track
 * <ol>
 *  <li>for supercontigs (scaffolds), emits contig RefSeq acc as the chromosome (f.e. NW_004956882)
 *  <li>Skips 'chromosome' lines
 *  <li>'biological_region' lines are written to 'features' file; everything else to 'model' file
 *  <li>'description=' attribute is changed to 'notes=' attribute;
 *      purpose: to turn off automatic display of descriptions, which could be lengthy and confusing
 *     (f.e. MGI source information is shown for hundreds of rat genes)
 * </ol>
 */
public class EnsemblPrep {
    private Map<Integer, String> ensemblGff;
    private Map<Integer, String> ensemblJBrowseDataDirs;
    private String outDir;
    Logger log = LogManager.getLogger("ensembl");
    //SimpleDateFormat sdt = new SimpleDateFormat("yyyyMMdd");

    public void run() throws Exception {

        // ensure output directory does exist
        File outputDir = new File(getOutDir());
        outputDir.mkdirs();

        // open input and output files
        Set<Integer> mapKeys = ensemblGff.keySet();
        for(Integer mapKey : mapKeys) {
            String inputFile = downloadEnsemblGffFile(ensemblGff.get(mapKey));
            int pos = inputFile.indexOf(".gff3");
            log.info("Downloaded gff file from Ensembl: "+inputFile);

            String modelFileName = inputFile.substring(0, pos) + "-model" + inputFile.substring(pos);
            String featureFileName = inputFile.substring(0, pos) + "-feature" + inputFile.substring(pos);

            //int lastSlashPos = inputFile.lastIndexOf("/");
            //int dotPos = inputFile.indexOf(".",lastSlashPos);
            //String assemblyName = inputFile.substring(dotPos+1, pos); // f.e. GRCh38.106

            BufferedReader in = Utils.openReader(inputFile);
            BufferedWriter modelFile = Utils.openWriter(modelFileName);
            BufferedWriter featureFile = Utils.openWriter(featureFileName);

            log.info("opened file " + inputFile);
            log.info("writing file " + modelFileName);
            log.info("writing file " + featureFileName);

            int commentLines = 0;
            int notesLines = 0;
            int modelFileLines = 0;
            int featureFileLines = 0;
            int skippedHashLines = 0;
            int badChrLines = 0;

            Map<String, String> genbankToRefseqAccMap = parseSupercontigs(inputFile);
            Set<String> chromosomes = new HashSet<>();

            String line;
            while ((line = in.readLine()) != null) {
                // copy comment lines to both output files
                if (line.startsWith("#")) {
                    if( line.equals("###") ) {
                        skippedHashLines++;
                        continue;
                    }
                    modelFile.write(line + "\n");
                    featureFile.write(line + "\n");
                    commentLines++;
                    continue;
                }
                // skip chromosome lines
                if (line.contains("ID=chromosome:") || line.contains("ID=region:") ) {
                    continue;
                }
                // process supercontig lines
                if (line.contains("ID=supercontig:")) {
                    continue;
                }
                // convert description attr into notes
                pos = line.indexOf(";description=");
                if (pos > 0) {
                    line = line.substring(0, pos) + ";notes=" + line.substring(pos + 13);
                    notesLines++;
                }
                // prepend line with 'Chr'
                if (!line.startsWith("Chr")) {
                    int chrLen = line.indexOf('\t');
                    if (chrLen <= 3) {
                        ;
                    } else {
                        // replace genbank acc ids with refseq acc ids for supercontigs (scaffolds)
                        String genbankAcc = line.substring(0, chrLen);
                        String refseqAcc = genbankToRefseqAccMap.get(genbankAcc);
                        if (refseqAcc == null) {
                            log.debug("*** WARN: null scaffold acc "+genbankAcc);
                            badChrLines++;
                            continue;
                        }
                        else if (chromosomes.add(refseqAcc)) {
                            log.info(refseqAcc);
                        }
                        line = refseqAcc + line.substring(chrLen);
                    }
                }

                // write out the line to the proper file
                if (line.contains("biological_region")) {
                    featureFile.write(line);
                    featureFile.write("\n");
                    featureFileLines++;
                } else {
                    modelFile.write(line);
                    modelFile.write("\n");
                    modelFileLines++;
                }
            }

            // close the files
            in.close();
            modelFile.close();
            featureFile.close();

            log.info("");
            log.info("comment lines: " + commentLines);
            log.info("converted notes lines: " + notesLines);
            log.info("skipped ### lines: " + skippedHashLines);
            log.info("invalid chromosome lines: " + badChrLines);
            log.info("data lines written to model file: " + modelFileLines);
            log.info("data lines written to feature file: " + featureFileLines);
            log.info("***********************\n\n");

            Gff3ColumnWriter.sortInMemory(modelFileName, Gff3ColumnWriter.COMPRESS_MODE_BGZIP);
            Gff3ColumnWriter.sortInMemory(featureFileName, Gff3ColumnWriter.COMPRESS_MODE_BGZIP);
        }
    }

    static Map<String,String> parseSupercontigs(String inputFile) throws IOException {
        BufferedReader in = Utils.openReader(inputFile);
        System.out.println("opened file "+inputFile);

        int skippedChrLines = 0;

        Map<String,String> genbankToRefseqAccMap = new HashMap<>();

        String line;
        while( (line=in.readLine())!=null ) {
            // copy comment lines to both output files
            if( line.startsWith("#") ) {
                continue;
            }
            // skip chromosome lines
            if( line.contains("ID=chromosome:") || line.contains("ID=region:") ) {
                skippedChrLines++;
                continue;
            }
            // process supercontig lines
            if( line.contains("ID=supercontig:") ) {
                //AGCD01080321.1	Ensembl	supercontig	1	2322	.	.	.	ID=supercontig:AGCD01080321.1;Alias=NW_004956926.1
                int pos2 = line.indexOf("ID=supercontig:") + "ID=supercontig:".length();
                int pos3 = line.indexOf(";Alias=", pos2);
                int pos4 = line.indexOf(".", pos3);
                String genbankAcc = line.substring(pos2, pos3);
                String refseqAcc;
                if( pos4>pos3 ) {
                    refseqAcc = line.substring(pos3 + ";Alias=".length(), pos4);
                } else {
                    refseqAcc = line.substring(pos3 + ";Alias=".length()).trim();
                }
                int pos5 = refseqAcc.indexOf(",");
                if( pos5>0 ) {
                    refseqAcc = refseqAcc.substring(pos5+1);
                }
                genbankToRefseqAccMap.put(genbankAcc, refseqAcc);
                skippedChrLines++;
                continue;
            }
        }

        // close the files
        in.close();

        System.out.println("");
        System.out.println("skipped chromosome lines: "+skippedChrLines);

        return genbankToRefseqAccMap;
    }

    String downloadEnsemblGffFile(String file) throws Exception{

        int lastSlashPos = file.lastIndexOf('/');
        String localFileName = getOutDir()+file.substring(lastSlashPos+1);
        FileDownloader downloader = new FileDownloader();
        downloader.setExternalFile(file);
        downloader.setLocalFile(localFileName);
        downloader.setUseCompression(true);
        return downloader.downloadNew();
    }

    public void setEnsemblGff(Map<Integer, String> ensemblGff) {
        this.ensemblGff = ensemblGff;
    }

    public Map<Integer, String> getEnsemblGff() {
        return ensemblGff;
    }

    public String getOutDir() {
        return outDir;
    }

    public void setOutDir(String outDir) {
        this.outDir = outDir;
    }

    public Map<Integer, String> getEnsemblJBrowseDataDirs() {
        return ensemblJBrowseDataDirs;
    }

    public void setEnsemblJBrowseDataDirs(Map<Integer, String> ensemblJBrowseDataDirs) {
        this.ensemblJBrowseDataDirs = ensemblJBrowseDataDirs;
    }
}
