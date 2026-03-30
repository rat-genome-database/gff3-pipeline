#!/usr/bin/env bash
#
# Run all GFF3 generation scripts sequentially:
#   1. create_jbrowse2_files
#   2. clinvar
#   3. agr
#   4. eva
#   5. ensemblGff
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=rgd.devops@mcw.edu
fi

echo "=== runAll.sh started at $(date) ==="

echo "--- create_jbrowse2_files ---"
bash $APP_HOME/create_jbrowse2_files.sh

echo "--- clinvar ---"
bash $APP_HOME/clinvar.sh

echo "--- agr ---"
bash $APP_HOME/agr.sh

echo "--- eva ---"
bash $APP_HOME/eva.sh

echo "--- ensemblGff ---"
bash $APP_HOME/ensemblGff.sh

echo "=== runAll.sh finished at $(date) ==="
