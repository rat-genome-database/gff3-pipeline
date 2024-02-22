# generate gff3 files for human clinvar variants, assemblies 37 and 38
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

DATA_RELEASE_DIR="/home/rgddata/data_release/GFF3/ClinVar/"

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "REED" ]; then
  EMAILLIST=mtutaj@mcw.edu
else
  EMAILLIST=mtutaj@mcw.edu
fi

RUNLOAD="$APP_HOME/run.sh"

$RUNLOAD -object:clinvar > clinvar.log

echo "copy generated gff3 files to data_release directory"
echo "  rsync -avd data/Variants/ClinVar/ $DATA_RELEASE_DIR"
rsync -avd data/Variants/ClinVar/ $DATA_RELEASE_DIR
echo "  rsync OK!"
