# script env setup
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline/
OUTDIR=$APP_HOME/data/strain_specific_variants
#EMAILLIST=jrsmith@mcw.edu,RGD.Developers@mcw.edu
EMAILLIST=mtutaj@mcw.edu

CHR_RAT="1-20,X,Y,MT"

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

RUNLOAD="$APP_HOME/run.sh"

cd $APP_HOME

# Rnor6.0 samples
LOGFILE=cron_ratcn60.log
rm $LOGFILE
$RUNLOAD -patientID:600 -source:CN -toFile:$OUTDIR/rat6/ -chr:$CHR_RAT -compress  2>&1 | tee -a $LOGFILE
mailx -s "[$SERVER]Pipeline to create Gff3 data for ratcn 6.0 ran" $EMAILLIST<$LOGFILE

# Rnor5.0 samples
LOGFILE=cron_ratcn50.log
rm $LOGFILE
$RUNLOAD -patientID:500 -source:CN -toFile:$OUTDIR/rat5/ -chr:$CHR_RAT -compress  2>&1 | tee -a $LOGFILE
mailx -s "[$SERVER]Pipeline to create Gff3 data for ratcn 5.0 ran" $EMAILLIST<$LOGFILE

# RGSC3.4 samples
LOGFILE=cron_ratcn34.log
rm $LOGFILE
$RUNLOAD -patientID:180 -source:CN -toFile:$OUTDIR/rat3_4/ -chr:$CHR_RAT -compress  2>&1 | tee -a $LOGFILE
mailx -s "[$SERVER]Pipeline to create Gff3 data for ratcn 3.4 ran" $EMAILLIST<$LOGFILE
