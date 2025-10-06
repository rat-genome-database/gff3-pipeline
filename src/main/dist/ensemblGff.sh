# 1) download gff3 files from Ensembl for species configured in properties/AppConfigure.xml
# 2) split those files into model file and feature file
# 3) split files are stored in data/Ensembl directory
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

LOGDIR=$APP_HOME/data


SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "REED" ]; then
  EMAILLIST="mtutaj@mcw.edu llamers@mcw.edu"
else
  EMAILLIST=mtutaj@mcw.edu
fi

RUNLOAD="$APP_HOME/run.sh"

$RUNLOAD -flavor:ensembl_prep > ensembl.log

mailx -s "[$SERVER]Pipeline to download Ensembl Gff3 files ran" $EMAILLIST < logs/ensembl.log


