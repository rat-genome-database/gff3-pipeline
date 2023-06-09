# script env setup
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
RUNLOAD="$APP_HOME/run.sh"

$RUNLOAD -object:diseases > diseases.log
$RUNLOAD -object:chebi > chebi.log
$RUNLOAD -object:genes > genes.log
$RUNLOAD -object:qtls > qtls.log

