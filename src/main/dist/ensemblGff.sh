# 1) download gff3 files from Ensembl for species configured in properties/AppConfigure.xml
# 2) split those files into model file and feature file
# 3) generate bash script to load those files into JBrowse as JBrowse tracks
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

LOGDIR=$APP_HOME/logs
DATA_RELEASE_DIR=/home/rgddata/data_release/GFF3


SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "REED" ]; then
  EMAILLIST=mtutaj@mcw.edu,llamers@mcw.edu
else
  EMAILLIST=mtutaj@mcw.edu,llamers@mcw.edu
fi

RUNLOAD="$APP_HOME/run.sh"

$RUNLOAD -flavor:ensembl_prep > ensembl.log

mailx -s "[$SERVER]Pipeline to download Ensembl Gff3 files ran" $EMAILLIST < logs/ensembl.log

echo "copy generated gff3 files to data_release directory:"
echo "  rsync -avd $LOGDIR/Ensembl/ $DATA_RELEASE_DIR/Ensembl/"
rsync -avd $LOGDIR/Ensembl/ $DATA_RELEASE_DIR/Ensembl/
echo "  rsync OK!"
