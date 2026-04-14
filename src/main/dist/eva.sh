# generate gff3 files for Eva Variants
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

LOGDIR=$APP_HOME/data
DATA_RELEASE_DIR=/home/rgddata/data_release/GFF3


SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "REED" ]; then
  EMAILLIST="jrsmith@mcw.edu rgd.devops@mcw.edu"
#  EMAILLIST="mtutaj@mcw.edu llamers@mcw.edu"
else
  EMAILLIST="mtutaj@mcw.edu llamers@mcw.edu"
fi

RUNLOAD="$APP_HOME/run.sh"

$RUNLOAD -object:Eva > Eva.log

mailx -s "[$SERVER]Pipeline to create Gff3 data for Eva Variants ran" $EMAILLIST < logs/evas.log

echo "copy generated gff3 files to data_release directory:"
echo "  rsync -avd $LOGDIR/Eva/ $DATA_RELEASE_DIR/Eva/"
rsync -avd $LOGDIR/Eva/ $DATA_RELEASE_DIR/Eva/
echo "  rsync OK!"
