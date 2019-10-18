package edu.mcw.rgd.gff3;

import edu.mcw.rgd.process.Utils;

import java.io.File;
import java.net.*;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;

/**
 * @author mtutaj
 * @since 1/28/15
 */
public class FileGuard {

    private Date time0; // time when the processing started
    private String dirName; // name of dir to guard
    private String ontId; // f.e. DOID, CHEBI
    private Set<String> originalAccIds;

    void init(String dirName, String ontId) {
        this.time0 = new Date();
        this.dirName = dirName;
        this.ontId = ontId;

        originalAccIds = loadAccIds();
    }

    void check(int mapKey) throws UnknownHostException {
        // keep track of the accession ids that are new, current or stale
        Set<String> newAccIds = new HashSet<String>();
        Set<String> staleAccIds = new HashSet<String>();

        File dir = new File(dirName);
        for( File file: dir.listFiles() ) {
            String fileName = file.getName();
            if( fileName.startsWith(ontId) && fileName.endsWith("gff3.gz")) {
                int pos = fileName.indexOf("_Ontology");
                if( pos>0 ) {
                    String accId = fileName.substring(0, pos);
                    // check the file
                    long timeLastModified = file.lastModified();
                    if( timeLastModified>time0.getTime() ) {
                        if( !originalAccIds.contains(accId) ) {
                            newAccIds.add(accId);
                        }
                    } else {
                        staleAccIds.add(accId);
                    }
                }
            }
        }

        String msg = "";
        if( !newAccIds.isEmpty() ) {
            msg += "NEW "+ontId+" TRACKS for acc ids: "+ Utils.concatenate(newAccIds,",")+" \n";
        }
        if( !staleAccIds.isEmpty() ) {
            msg += "DISCONTINUED "+ontId+" TRACKS for acc ids: "+ Utils.concatenate(staleAccIds,",")+" \n";
        }
        if( !msg.isEmpty() ) {
            InetAddress addr = InetAddress.getLocalHost();
            String hostname = addr.getHostName();

            String mailServer = "localhost";
            String mailFrom = "rgddata@"+hostname;
            String[] recipients = {"mtutaj@mcw.edu","jrsmith@mcw.edu"};
            String title = "["+hostname.toUpperCase()+"] JBROWSE "+ontId+" TRACKS WARNING for map_key="+mapKey;
            Utils.sendMail(mailServer, mailFrom, recipients, title, msg);
        }
    }

    Set<String> loadAccIds() {
        Set<String> accIds = new HashSet<String>();
        File dir = new File(dirName);
        for( String fileName: dir.list() ) {
            if( fileName.startsWith(ontId) && fileName.endsWith("gff3.gz") ) {
                int pos = fileName.indexOf("_Ontology");
                if( pos>0 )
                    accIds.add(fileName.substring(0, pos));
            }
        }
        return accIds;
    }
}
