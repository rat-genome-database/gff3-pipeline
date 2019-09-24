# script env setup
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline/
OUTDIR=$APP_HOME/data/strain_specific_variants
#EMAILLIST=jrsmith@mcw.edu,RGD.Developers@mcw.edu
EMAILLIST=hsnalabolu@mcw.edu

CHR_DOG="1-38,X,MT"

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

RUNLOAD="$APP_HOME/run.sh"

cd $APP_HOME


LOGFILE=cron_dog.log
rm $LOGFILE
$RUNLOAD -patientID:3 -source:CN -toFile:$OUTDIR/dog/ -chr:$CHR_DOG -compress  2>&1 | tee -a $LOGFILE
mailx -s "[$SERVER]Pipeline to create Gff3 data for dog ran" $EMAILLIST<$LOGFILE

