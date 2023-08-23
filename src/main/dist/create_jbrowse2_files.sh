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
$RUNLOAD -object:markers > markers.log
$RUNLOAD -object:strains > strains.log
$RUNLOAD -object:proteinDomains > domains.log
$RUNLOAD -object:ensembl > ensembl.log
$RUNLOAD -object:variants > variants.log
$RUNLOAD -object:jb2_eva > jb2_eva.log

