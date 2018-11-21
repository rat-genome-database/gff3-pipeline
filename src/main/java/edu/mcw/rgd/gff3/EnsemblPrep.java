package edu.mcw.rgd.gff3;

import edu.mcw.rgd.process.Utils;

import java.io.*;
import java.util.zip.GZIPOutputStream;

/**
 * @author mtutaj
 * @since 11/20/2018
 * Prepares a gff3 file downloaded from Ensembl FTP site to be loaded by JBrowse scripts as an Ensembl track
 * <ol>
 *  <li>Prepends all chromosome names with 'Chr'
 *  <li>Skips 'chromosome' lines
 *  <li>'biological_region' lines are written to 'features' file; everything else to 'model' file
 *  <li>'description=' attribute is changed to 'notes=' attribute;
 *      purpose: to turn off automatic display of descriptions, which could be lengthy and confusing
 *     (f.e. MGI source information is shown for hundreds of rat genes)
 * </ol>
 */
public class EnsemblPrep {

    public static void run(String inputFile) throws Exception {

        // open input and output files
        int pos = inputFile.indexOf(".gff3");
        String modelFileName = inputFile.substring(0, pos)+"-model"+inputFile.substring(pos);
        String featureFileName = inputFile.substring(0, pos)+"-feature"+inputFile.substring(pos);

        BufferedReader in = Utils.openReader(inputFile);
        BufferedWriter modelFile = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(modelFileName))));
        BufferedWriter featureFile = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(featureFileName))));

        System.out.println("opened file "+inputFile);
        System.out.println("writing file "+modelFileName);
        System.out.println("writing file "+featureFileName);

        int commentLines = 0;
        int skippedChrLines = 0;
        int notesLines = 0;
        int modelFileLines = 0;
        int featureFileLines = 0;

        String line;
        while( (line=in.readLine())!=null ) {
            // copy comment lines to both output files
            if( line.startsWith("#") ) {
                modelFile.write(line+"\n");
                featureFile.write(line+"\n");
                commentLines++;
                continue;
            }
            // skip chromosome lines
            if( line.contains("ID=chromosome:") ) {
                skippedChrLines++;
                continue;
            }
            // convert description attr into notes
            pos = line.indexOf(";description=");
            if( pos>0 ) {
                line = line.substring(0, pos)+";notes="+line.substring(pos+13);
                notesLines++;
            }
            // write out the line to the proper file
            if( line.contains("biological_region") ) {
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

        System.out.println("");
        System.out.println("comment lines: "+commentLines);
        System.out.println("skipped chromosome lines: "+skippedChrLines);
        System.out.println("converted notes lines: "+notesLines);
        System.out.println("data lines written to model file: "+modelFileLines);
        System.out.println("data lines written to feature file: "+featureFileLines);

    }
}
