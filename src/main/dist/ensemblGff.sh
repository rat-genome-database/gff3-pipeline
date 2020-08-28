# generate gff3 files for Eva Variants
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

LOGDIR=$APP_HOME/logs
DATA_RELEASE_DIR=/home/rgddata/data_release/GFF3


SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "REED" ]; then
  EMAILLIST=jrsmith@mcw.edu,RGD.Developers@mcw.edu
#  EMAILLIST=mtutaj@mcw.edu,llamers@mcw.edu
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
