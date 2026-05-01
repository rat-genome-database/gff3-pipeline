# NOTE: this script REPLACES CreateRatStrainSpecificVariantGff3.sh.
#       The old script made one run.sh call per patient and wrote into
#       data/strain_specific_variants/ratN/, then a separate -object:variants
#       step copied files into the final tree. This single script does both
#       in one pass.
#
# Generate per-sample strain-specific variant GFF3 files for jbrowse2.
# Patient list and per-mapKey output prefixes live in properties/AppConfigure.xml under the
# 'gff3VariantsManager' bean. Output goes to data/jbrowse2_variants/Rat/<assembly>/Variants/...
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline/
EMAILLIST="mtutaj@mcw.edu llamers@mcw.edu"
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

cd $APP_HOME

LOGFILE=cron_jbrowse2_gff3_variants.log
rm -f $LOGFILE

$APP_HOME/run.sh -object:gff3Variants 2>&1 | tee -a $LOGFILE
mailx -s "[$SERVER] jbrowse2 strain-specific variants GFF3 pipeline ran" $EMAILLIST < $LOGFILE
