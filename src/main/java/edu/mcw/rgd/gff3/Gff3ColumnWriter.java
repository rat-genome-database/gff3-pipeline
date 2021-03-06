package edu.mcw.rgd.gff3;

import java.io.*;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

/**
 * @author pjayaraman
 * @since 9/9/11
 * @version %I%, %G%
 * <p>
 * convenience class to handle writing of gff or gff3 files
 */
public class Gff3ColumnWriter {

    private PrintWriter gff3Writer;
    private boolean ratmineCompatibleFormat;
    private boolean agrCompatibleFormat;
    private String gff3FileName;

    // init gff writer in gff3 format; do not compress output
    public Gff3ColumnWriter(String fileName) throws IOException {
        init(fileName, false, false);
    }

    // init gff writer in gff or gff3 format; do not compress output
    public Gff3ColumnWriter(String fileName, boolean useGffFormat) throws IOException {
        init(fileName, useGffFormat, false);
    }

    // init gff writer in gff or gff3 format; do not compress output
    public Gff3ColumnWriter(String fileName, boolean useGffFormat, boolean compress) throws IOException {
        init(fileName, useGffFormat, compress);
    }

    private void init(String fileName, boolean useGffFormat, boolean compress) throws IOException {
        // ensure the directory is created
        int lastSlashPos = fileName.lastIndexOf('/');
        if( lastSlashPos < 0 )
            lastSlashPos = fileName.lastIndexOf('\\');
        if( lastSlashPos>0 ) {
            new File(fileName.substring(0, lastSlashPos)).mkdirs();
        }

        String outFileName = fileName;
        if( compress ) {
            if( !outFileName.endsWith(".gz") ) {
                outFileName += ".gz";
            }
            gff3Writer = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFileName))));
        }
        else {
            gff3Writer = new PrintWriter(outFileName);
        }
        gff3FileName = fileName;

        if( !useGffFormat )
            gff3Writer.println("##gff-version 3");
    }

    public String getFileName() {
        return gff3FileName;
    }

    public void close() {
        gff3Writer.close();
    }

    public void print(String str) {
        gff3Writer.print(str);
    }

    //write Gff3 file for each line.
    // This is a combined method that can be called for each line.
    public void writeGff3AllColumns(String chrNum, String source, String type, Integer start, Integer stop, String score,
                                   String strand, String phase, Map<String, String> attributesMap) throws Exception{
        writeFirst8Columns(chrNum,source, type,start,stop,score,strand,phase);
        writeAttributesForGff3(attributesMap);
    }

    /**
     * write the first 8 columns for each file; has to be called for each line that gets printed into the gff3 file
     * @param chrNum
     * @param source
     * @param type
     * @param start
     * @param stop
     * @param score
     * @param strand
     * @param phase
     * @return the contents written
     * @throws Exception
     */
    public String writeFirst8Columns(String chrNum, String source, String type, Integer start, Integer stop, String score,
                                   String strand, String phase) throws Exception{

        if( start>stop ) {
            System.out.println("WARNING: reverse start pos > stop pos");
        }

        String chr = chrNum;
        if( isRatmineCompatibleFormat() ) {
            switch(type) {
                case "EXON": type = "Exon"; break;
                case "UTR5": type = "FivePrimeUTR"; break;
                case "UTR3": type = "ThreePrimeUTR"; break;
            }
        } else if( !isAgrCompatibleFormat() ) {
            if( chrNum.length() <= 2 ) {
                chr = "Chr" + chrNum;
            }
        }

        String text = chr + "\t" + source + "\t" + type + "\t" + start + "\t" + stop + "\t" + score + "\t" + strand + "\t" + phase + "\t";
        gff3Writer.print(text);
        return text;
    }

    public void writeAttributes4Gff(Map<String, String> attributesMap) throws Exception{
        for (Map.Entry<String, String> entry: attributesMap.entrySet()) {
            gff3Writer.print(entry.getKey() + " \"" + entry.getValue() + "\"; ");
        }
    }

    //add new line.
    public void addnewLineInGff3(){
        gff3Writer.print("\n");
    }

    /**
     * writes the attributes, then writes the new line; in the end drops all attributes from attr map
     * @param attributesMap map of attributes to be written
     * @return return count of attributes written
     * @throws Exception
     * @deprecated
     */
    public int writeAttributesForGff3(Map<String, String> attributesMap) throws Exception{
        int attrCount = 0;
        for( Map.Entry<String, String> entry: attributesMap.entrySet() ) {
            if( attrCount>0 )
                gff3Writer.print(';');
            gff3Writer.print(entry.getKey() + "=" + entry.getValue());
            attrCount++;
        }

        addnewLineInGff3();
        attributesMap.clear();
        return attrCount;
    }

    /**
     * writes the attributes, then writes the new line; in the end drops all attributes from attr map;
     * note: this method automatically percent-encodes TAB,CR,NL,'%',';','=', per GFF3 spec
     * @param attributesMap map of attributes to be written
     * @return return count of attributes written
     * @throws Exception
     */
    public int writeAttributes4Gff3(Map<String, String> attributesMap) throws Exception{
        int attrCount = 0;
        for( Map.Entry<String, String> entry: attributesMap.entrySet() ) {
            if( attrCount>0 )
                gff3Writer.print(';');
            gff3Writer.print(entry.getKey() + "=" + encodeAttrValueForGff3(entry.getValue()));
            attrCount++;
        }

        addnewLineInGff3();
        attributesMap.clear();
        return attrCount;
    }

    /** per GFF3 specification:<p>
     * Literal use of tab, newline, carriage return, the percent (%) sign, and control characters must be encoded
     * using RFC 3986 Percent-Encoding; no other characters may be encoded.
     * <p>
     * We also encode ;=
     */
    public String encodeAttrValueForGff3(String value) {
        if( value==null ) {
            return "null";
        }

        StringBuilder buf = new StringBuilder();
        for( int i=0; i<value.length(); i++ ) {
            char c = value.charAt(i);
            switch(c) {
                case '\t': buf.append("%09"); break;
                case '\r': buf.append("%0D"); break;
                case '\n': buf.append("%0A"); break;
                case ';': buf.append("%3B"); break;
                case '=': buf.append("%3D"); break;
                case '%':
                    // sometimes the value contains an url, that is already percent-encoded
                    // if we percent-encode percent char for that url, the url will be broken in JBrowse!
                    // so we must detect the already percent-encoded urls in the source data to prevent double encoding
                    int pos1 = value.lastIndexOf("href=", i);
                    int pos2 = value.indexOf("</a>");
                    if( pos1>0 && pos2>i ) {
                        // already percent-encoded, only emit '%'
                        buf.append(c);
                    } else {
                        buf.append("%25");
                    }
                    break;

                default: buf.append(c); break;
            }
        }
        return buf.toString();
    }

    public boolean isRatmineCompatibleFormat() {
        return ratmineCompatibleFormat;
    }

    public void setRatmineCompatibleFormat(boolean ratmineCompatibleFormat) {
        this.ratmineCompatibleFormat = ratmineCompatibleFormat;
    }

    public boolean isAgrCompatibleFormat() {
        return agrCompatibleFormat;
    }

    public void setAgrCompatibleFormat(boolean agrCompatibleFormat) {
        this.agrCompatibleFormat = agrCompatibleFormat;
    }
}
