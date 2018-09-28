package edu.mcw.rgd.gff3;

import java.io.PrintWriter;

/**
 * Created by IntelliJ IDEA.
 * User: pjayaraman
 * Date: 9/9/11
 * Time: 8:50 AM
 */
public class FastaWriter {

    PrintWriter fastaSeqWriter;

    public FastaWriter(PrintWriter fastaSeqWriter) {
          this.fastaSeqWriter = fastaSeqWriter;
    }

    public PrintWriter getFastaSeqWriter() {
        return fastaSeqWriter;
    }

    public void setFastaSeqWriter(PrintWriter fastaSeqWriter) {
        this.fastaSeqWriter = fastaSeqWriter;
    }

    //write both forward and reverse sequences with header.
    // This method has to be called for each line that gets printed into the fasta file.
    public void writeFastaSequences(String fwdSeqHeader, String fwdSeq, String revSeqHeader, String revSeq) throws Exception{
        writeFastaSeq(fwdSeqHeader, fwdSeq);
        writeFastaSeq(revSeqHeader, revSeq);
    }

    public void writeFastaSeq(String seqheader, String seq) throws Exception{
        if(seq!=null && seq.length()>0){
            fastaSeqWriter.println(">"+seqheader);
            fastaSeqWriter.println(seq);
        }
    }

    //add new line.
    public void addnewLineInFasta(){
        fastaSeqWriter.print("\n");
    }



}
