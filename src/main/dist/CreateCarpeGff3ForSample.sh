# script env setup
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline/
OUTDIR=$APP_HOME/data/strain_specific_variants
EMAILLIST=mtutaj@mcw.edu

CHR_RAT="1-20,X,Y,MT"

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

RUNLOAD="$APP_HOME/run.sh"

cd $APP_HOME

SAMPLE=$1

if [ $SAMPLE -ge 500 ] && [ $SAMPLE -le 699 ]; then
  RATDIR=rat34
elif [ $SAMPLE -ge 700 ] && [ $SAMPLE -le 899 ]; then
  RATDIR=rat50
elif [ $SAMPLE -ge 900 ] && [ $SAMPLE -le 999 ]; then
  RATDIR=rat60
else
  RATDIR="."
fi

LOGFILE="cron_ratcn${SAMPLE}.log"
$RUNLOAD -sampleID:$SAMPLE -source:CN -toFile:$OUTDIR/$RATDIR/ -chr:$CHR_RAT -compress  2>&1 | tee -a $LOGFILE
mailx -s "[$SERVER]Pipeline to create Gff3 data for sample $SAMPLE ran" $EMAILLIST<$LOGFILE
