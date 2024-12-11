package edu.mcw.rgd.gff3;

import edu.mcw.rgd.process.Utils;
import htsjdk.samtools.util.BlockCompressedOutputStream;

import java.io.*;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

public class GtfWriter {
    public static int COMPRESS_MODE_NONE = 0;
    public static int COMPRESS_MODE_ZIP = 1;
    public static int COMPRESS_MODE_BGZIP = 2;

    private PrintWriter writer;
    private int compressMode;
    private String gtfFileName;
    private String outFileName;

    // do not compress output
    public GtfWriter(String fileName) throws IOException {
        init(fileName, COMPRESS_MODE_NONE);
    }

    // do not compress output
    public GtfWriter(String fileName, int compressMode) throws IOException {
        init(fileName, compressMode);
    }

    private void init(String fileName, int compressMode) throws IOException {
        // ensure the directory is created
        int lastSlashPos = fileName.lastIndexOf('/');
        if( lastSlashPos < 0 )
            lastSlashPos = fileName.lastIndexOf('\\');
        if( lastSlashPos>0 ) {
            new File(fileName.substring(0, lastSlashPos)).mkdirs();
        }

        outFileName = fileName;

        this.compressMode = compressMode;
        if( compressMode != COMPRESS_MODE_NONE ) {
            if( !outFileName.endsWith(".gz") ) {
                outFileName += ".gz";
            }
            if( compressMode == COMPRESS_MODE_ZIP ) {
                writer = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFileName))));
            } else if( compressMode == COMPRESS_MODE_BGZIP ) {
                writer = new PrintWriter(new OutputStreamWriter(new BlockCompressedOutputStream(outFileName)));
            }
        }
        else {
            writer = new PrintWriter(outFileName);
        }
        gtfFileName = fileName;

        writer.println("#gtf-version 2.2");
    }

    public String getFileName() {
        return gtfFileName;
    }
    public String getOutFileName() {
        return outFileName;
    }

    public void close() {
        writer.close();
    }

    public void print(String str) {
        writer.print(str);
    }

    /**
     * write the first 8 columns for each file; has to be called for each line that gets printed into the gtf file
     * @return the contents written
     */
    public String writeFirst8Columns(String chrNum, String source, String type, Integer start, Integer stop, String score, String strand, String phase) {

        String text = prepFirst8Columns(chrNum, source, type, start, stop, score, strand, phase);
        writer.print(text);
        return text;
    }

    public String prepFirst8Columns(String chrNum, String source, String type, Integer start, Integer stop, String score, String strand, String phase) {

        String chr = chrNum;
        switch(type) {
            case "gene": break;
            case "transcript": break;
            case "exon": break;
            case "CDS": break;
            default:
                System.out.println("WARNING! unexpected object type: "+type);
        }

        // ensure strand is not null
        strand = Utils.NVL(strand, "");

        if( start>stop ) {
            System.out.println("WARNING: reverse start pos > stop pos: " + getFileName() + " " +
                chr + "\t" + source + "\t" + type + "\t" + start + "\t" + stop + "\t" + score + "\t" + strand + "\t" + phase);

            int tmp = start;
            start = stop;
            stop = tmp;
        }

        String text = chr + "\t" + source + "\t" + type + "\t" + start + "\t" + stop + "\t" + score + "\t" + strand + "\t" + phase + "\t";
        return text;
    }

    /**
     * writes the attributes, then writes the new line; in the end drops all attributes from attr map;
     * note: this method automatically percent-encodes TAB,CR,NL,'%',';','=', per GFF3 spec
     * @param attributesMap map of attributes to be written
     * @return return count of attributes written
     * @throws Exception
     */
    public int writeAttributes( String geneId, String transcriptId, Map<String, String> attributesMap ) {

        String str = prepAttributes( geneId, transcriptId, attributesMap );
        writer.print(str);
        int attrCount = attributesMap.size();
        return attrCount;
    }

    /**
     * writes the attributes, then writes the new line; in the end drops all attributes from attr map;
     * note: this method automatically percent-encodes TAB,CR,NL,'%',';','=', per GFF3 spec
     * @return gtf content
     */
    public String prepAttributes( String geneId, String transcriptId, Map<String, String> attrMap ) {

        StringBuilder buf = new StringBuilder();

        // gene_id is mandatory
        add( "gene_id", geneId, buf);

        // transcript_id is mandatory
        add( "transcript_id", transcriptId, buf);

        for( Map.Entry<String, String>entry: attrMap.entrySet() ) {
            add( entry.getKey(), entry.getValue(), buf);
        }

        buf.append("\n");
        attrMap.clear();
        return buf.toString();
    }

    void add( String attrName, String attrValue, StringBuilder buf ) {

        if( buf.length()>0 )
            buf.append(" ");

        buf.append(attrName+ " \"" + Gff3ColumnWriter.encodeAttrValueForGff3(attrValue) +"\";");
    }

    /**
     * sorts the given file by 1 column (chromosome) and then by start pos (4th col) and stop pos (5th col)
     * goal: tabix-ready file
     */
    public void sortInMemory() throws IOException {

        String fname = getFileName();
        Gff3ColumnWriter.sortInMemory(fname, compressMode);
    }

    public int getCompressMode() {
        return compressMode;
    }

    public void setCompressMode(int compressMode) {
        this.compressMode = compressMode;
    }
}
