# generate gff3 files for human dbSnp variants
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

LOGDIR=$APP_HOME/data
DATA_RELEASE_DIR=/home/rgddata/data_release/GFF3


SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "REED" ]; then
  EMAILLIST=mtutaj@mcw.edu
else
  EMAILLIST=mtutaj@mcw.edu
fi

RUNLOAD="$APP_HOME/run.sh"

$RUNLOAD -object:dbSnp -species:HUMAN -mapKey:38 -build:dbSnp150 -toFile:$LOGDIR/DbSnp/DbSnp_hg38.gff3 -compress &> dbSnp_hg38.log
mailx -s "[$SERVER] Pipeline to create DbSnp Gff3 data for human assembly 38 ran" $EMAILLIST<dbSnp_hg38.log

$RUNLOAD -object:dbSnp -species:HUMAN -mapKey:17 -build:dbSnp150 -toFile:$LOGDIR/DbSnp/DbSnp_hg19.gff3 -compress &> dbSnp_hg19.log
mailx -s "[$SERVER] Pipeline to create DbSnp Gff3 data for human assembly 37 ran" $EMAILLIST<dbSnp_hg19.log


echo "copy generated gff3 files to data_release directory:"
echo "  rsync -avd $LOGDIR/DbSnp/ $DATA_RELEASE_DIR/DbSnp/"
rsync -avd $LOGDIR/DbSnp/ $DATA_RELEASE_DIR/DbSnp/
echo "  rsync OK!"
