package edu.mcw.rgd.gff3;

import edu.mcw.rgd.process.Utils;
import htsjdk.samtools.util.BlockCompressedOutputStream;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Map;
import java.util.TreeSet;
import java.util.zip.GZIPOutputStream;

/**
 * @author pjayaraman
 * @since 9/9/11
 * convenience class to handle writing of gff3 files
 */
public class Gff3ColumnWriter {

    public static int COMPRESS_MODE_NONE = 0;
    public static int COMPRESS_MODE_ZIP = 1;
    public static int COMPRESS_MODE_BGZIP = 2;

    private PrintWriter gff3Writer;
    private int compressMode;
    private boolean agrCompatibleFormat;
    private String gff3FileName;
    private String outFileName;

    // init gff writer in gff3 format; do not compress output
    public Gff3ColumnWriter(String fileName) throws IOException {
        init(fileName, COMPRESS_MODE_NONE);
    }

    // init gff writer in gff3 format; do not compress output
    public Gff3ColumnWriter(String fileName, int compressMode) throws IOException {
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
                gff3Writer = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFileName))));
            } else if( compressMode == COMPRESS_MODE_BGZIP ) {
                gff3Writer = new PrintWriter(new OutputStreamWriter(new BlockCompressedOutputStream(outFileName)));
            }
        }
        else {
            gff3Writer = new PrintWriter(outFileName);
        }
        gff3FileName = fileName;

        gff3Writer.println("##gff-version 3");
    }

    public String getFileName() {
        return gff3FileName;
    }
    public String getOutFileName() {
        return outFileName;
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
                                   String strand, String phase, Map<String, String> attributesMap) {
        writeFirst8Columns(chrNum,source, type,start,stop,score,strand,phase);
        writeAttributesForGff3(attributesMap);
    }

    /**
     * write the first 8 columns for each file; has to be called for each line that gets printed into the gff3 file
     * @return the contents written
     */
    public String writeFirst8Columns(String chrNum, String source, String type, Integer start, Integer stop, String score, String strand, String phase) {

        String text = prepFirst8Columns(chrNum, source, type, start, stop, score, strand, phase);
        gff3Writer.print(text);
        return text;
    }

    public String prepFirst8Columns(String chrNum, String source, String type, Integer start, Integer stop, String score, String strand, String phase) {

        String chr = chrNum;
        if( false ) {
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

        // ensure strand is not null
        strand = Utils.NVL(strand, "");

        String text = chr + "\t" + source + "\t" + type + "\t" + start + "\t" + stop + "\t" + score + "\t" + strand + "\t" + phase + "\t";

        if( start>stop ) {
            System.out.println("WARNING: reverse start pos > stop pos: "+getFileName()+" "+text);

            int tmp = start;
            start = stop;
            stop = tmp;
            text = chr + "\t" + source + "\t" + type + "\t" + start + "\t" + stop + "\t" + score + "\t" + strand + "\t" + phase + "\t";
        }
        return text;
    }

    //add new line.
    public void addnewLineInGff3(){
        gff3Writer.print("\n");
    }

    /**
     * writes the attributes, then writes the new line; in the end drops all attributes from attr map
     * @param attributesMap map of attributes to be written
     * @return return count of attributes written
     * @deprecated
     */
    @Deprecated
    public int writeAttributesForGff3(Map<String, String> attributesMap) {
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

    /**
     * writes the attributes, then writes the new line; in the end drops all attributes from attr map;
     * note: this method automatically percent-encodes TAB,CR,NL,'%',';','=', per GFF3 spec
     * @param attributesMap map of attributes to be written
     * @return gff3 content
     */
    static public String prepAttributes4Gff3(Map<String, String> attributesMap) {
        StringBuilder buf = new StringBuilder();
        int attrCount = 0;
        for( Map.Entry<String, String> entry: attributesMap.entrySet() ) {
            if( attrCount>0 )
                buf.append(";");
            buf.append(entry.getKey() + "=" + encodeAttrValueForGff3(entry.getValue()));
            attrCount++;
        }

        buf.append("\n");
        attributesMap.clear();
        return buf.toString();
    }

    /** per GFF3 specification:<p>
     * Literal use of tab, newline, carriage return, the percent (%) sign, and control characters must be encoded
     * using RFC 3986 Percent-Encoding; no other characters may be encoded.
     * <p>
     * We also encode ;=
     */
    static public String encodeAttrValueForGff3(String value) {
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

    public void sortAndRemoveDuplicates() throws IOException {

        String fname = getFileName();
        sortAndRemoveDuplicates(fname, compressMode);
    }

    static public void sortAndRemoveDuplicates( String fname, int compressMode ) throws IOException {

        if( compressMode!=Gff3ColumnWriter.COMPRESS_MODE_NONE ) {
            if( !fname.endsWith(".gz") ) {
                fname += ".gz";
            }
        }

        int linesRead = 0;
        int linesWritten = 0;
        BufferedReader in = Utils.openReader(fname);
        String line;
        TreeSet<String> lines = new TreeSet<>(new Gff3Comparator());
        while( (line=in.readLine())!=null ) {
            linesRead++;
            lines.add(line);
        }
        in.close();

        BufferedWriter out;
        if( compressMode == COMPRESS_MODE_BGZIP ) {
            out = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(fname)));
        } else {
            out = Utils.openWriter(fname);
        }

        for( String l: lines ) {
            out.write(l);
            out.write("\n");
            linesWritten++;
        }
        out.close();

        System.out.println("sortAndRemoveDuplicates for "+fname+":   lines read: "
                +Utils.formatThousands(linesRead)+",    lines written: "+Utils.formatThousands(linesWritten));
    }


    /**
     * sorts the given gff3 file by 1 column (chromosome) and then by start pos (4th col) and stop pos (5th col)
     * goal: tabix-ready gff3 file
     */
    public void sortInMemory() throws IOException {

        String fname = getFileName();
        sortInMemory(fname, compressMode);
    }

    /**
     * sorts the given gff3 file by 1 column (chromosome) and then by start pos (4th col) and stop pos (5th col)
     * goal: tabix-ready gff3 file
     */
    static public void sortInMemory( String fname, int compressMode ) throws IOException {

        if( compressMode!=Gff3ColumnWriter.COMPRESS_MODE_NONE ) {
            if( !fname.endsWith(".gz") ) {
                fname += ".gz";
            }
        }

        BufferedReader in = Utils.openReader(fname);
        String line;
        ArrayList<String> lines = new ArrayList<>();
        while( (line=in.readLine())!=null ) {
            lines.add(line);
        }
        in.close();

        lines.sort(new Gff3Comparator());

        BufferedWriter out;
        if( compressMode == COMPRESS_MODE_BGZIP ) {
            out = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(fname)));
        } else {
            out = Utils.openWriter(fname);
        }

        for( String l: lines ) {
            out.write(l);
            out.write("\n");
        }
        out.close();
    }

    static class Gff3Comparator implements Comparator<String> {
        public int compare(String o1, String o2) {
            String[] cols1 = o1.split("[\\t]", -1);
            String[] cols2 = o2.split("[\\t]", -1);
            if (cols1.length >= 9 && cols2.length >= 9) {
                // compare chromosomes
                int r = cols1[0].compareTo(cols2[0]);
                if (r != 0) {
                    return r;
                }
                // compare start positions
                int pos1 = Integer.parseInt(cols1[3]);
                int pos2 = Integer.parseInt(cols2[3]);
                if (pos1 != pos2) {
                    return pos1 - pos2;
                }
                // compare end positions
                pos1 = Integer.parseInt(cols1[4]);
                pos2 = Integer.parseInt(cols2[4]);
                if (pos1 != pos2) {
                    return pos1 - pos2;
                }
                // chr, start & stop positions the same: just compare the lines
                return o1.compareTo(o2);
            } else {
                // '##' lines have absolute priority before the rest of lines
                int v1 = o1.startsWith("##") ? 0 : 1;
                int v2 = o2.startsWith("##") ? 0 : 1;
                if (v1 != v2) {
                    return v1 - v2;
                }
                return o1.compareTo(o2);
            }
        }
    }


    public int getCompressMode() {
        return compressMode;
    }

    public void setCompressMode(int compressMode) {
        this.compressMode = compressMode;
    }

    public boolean isAgrCompatibleFormat() {
        return agrCompatibleFormat;
    }

    public void setAgrCompatibleFormat(boolean agrCompatibleFormat) {
        this.agrCompatibleFormat = agrCompatibleFormat;
    }
}
