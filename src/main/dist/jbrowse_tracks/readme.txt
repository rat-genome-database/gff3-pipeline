
These scripts run under user 'rgdpub' in the directory '/data/data/gff3' (or /rgd/data/gff3)

The master script LoadGff3IntoJBrowse_Reed.sh must be placed in /rgd/scripts directory
and it is beeing run automatically by cron on weekly basis.

The directory /jbrowse is a remote mount to a directory on HOSHI server:
  - on HANSEN       /jbrowse points to /data/jbrowsedev on HOSHI
  - on REED         /jbrowse points to /data/jbrowsepipeline on HOSHI
  - on HANCOCK/OWEN /jbrowse points to /data/jbrowseprod on HOSHI
  
