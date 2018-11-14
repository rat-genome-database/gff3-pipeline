#  create files for AGR (alliance of genome resources)
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline/dist
cd $APP_HOME

LOGDIR=$APP_HOME/log/RGDGFF3/Output
DATA_RELEASE_DIR_AGR=/home/rgddata/data_release/agr


SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "REED" ]; then
  EMAILLIST=jrsmith@mcw.edu,RGD.Developers@mcw.edu
else
  EMAILLIST=mtutaj@mcw.edu
fi

RUNLOAD="$APP_HOME/run.sh"

$RUNLOAD -object:gene -species:RAT -mapKey:360 -toFile:$LOGDIR/AGR/ -chr:* -flavor:AGR -compress  &> cron_ratAGR.log
mailx -s "[$SERVER]Pipeline to create AGR Gff3 data for Rat Gene assembly 6.0 ran" $EMAILLIST<cron_ratAGR.log

$RUNLOAD -object:gene -species:HUMAN -mapKey:38 -toFile:$LOGDIR/AGR/ -chr:* -flavor:AGR -compress  &> cron_humanAGR.log
mailx -s "[$SERVER]Pipeline to create AGR Gff3 data for Human Gene assembly 38 ran" $EMAILLIST<cron_humanAGR.log

cp -p $LOGDIR/AGR/Rat_RGD_AGR.gff3.gz $DATA_RELEASE_DIR_AGR/genes_10116.gff3.gz
cp -p $LOGDIR/AGR/Human_RGD_AGR.gff3.gz $DATA_RELEASE_DIR_AGR/genes_9606.gff3.gz
