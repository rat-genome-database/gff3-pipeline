# generate gff3 files for human clinvar variants
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

LOGDIR=$APP_HOME/data
DATA_RELEASE_DIR=/home/rgddata/data_release/GFF3


SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "REED" ]; then
  #EMAILLIST=jrsmith@mcw.edu,RGD.Developers@mcw.edu
  EMAILLIST=mtutaj@mcw.edu
else
  EMAILLIST=mtutaj@mcw.edu
fi

RUNLOAD="$APP_HOME/run.sh"

$RUNLOAD -object:variant -species:HUMAN -mapKey:38 -toFile:$LOGDIR/ClinVar/ClinVar_hg38.gff3 -compress &> clinvar_hg38.log
mailx -s "[$SERVER] Pipeline to create ClinVar Gff3 data for human assembly 38 ran" $EMAILLIST<clinvar_hg38.log

$RUNLOAD -object:variant -species:HUMAN -mapKey:17 -toFile:$LOGDIR/ClinVar/ClinVar_hg19.gff3 -compress &> clinvar_hg19.log
mailx -s "[$SERVER] Pipeline to create ClinVar Gff3 data for human assembly 37 ran" $EMAILLIST<clinvar_hg19.log

$RUNLOAD -object:variant -species:HUMAN -mapKey:13 -toFile:$LOGDIR/ClinVar/ClinVar_hg18.gff3 -compress &> clinvar_hg18.log
mailx -s "[$SERVER] Pipeline to create ClinVar Gff3 data for human assembly 36 ran" $EMAILLIST<clinvar_hg18.log


echo "copy generated gff3 files to data_release directory:"
echo "  rsync -avd $LOGDIR/ClinVar/ $DATA_RELEASE_DIR/ClinVar/"
rsync -avd $LOGDIR/ClinVar/ $DATA_RELEASE_DIR/ClinVar/
echo "  rsync OK!"
