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
# data/jbrowse2/Human/<assembly>/Variants/ClinVar/  -->  $DATA_RELEASE_DIR<assembly>/
for CLINVAR_DIR in data/jbrowse2/Human/*/Variants/ClinVar; do
  ASSEMBLY=`basename $(dirname $(dirname $CLINVAR_DIR))`
  echo "  rsync -avd $CLINVAR_DIR/ $DATA_RELEASE_DIR$ASSEMBLY/"
  rsync -avd $CLINVAR_DIR/ $DATA_RELEASE_DIR$ASSEMBLY/
done
echo "  rsync OK!"
