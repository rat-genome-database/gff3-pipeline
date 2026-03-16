#!/usr/bin/env bash
# weekly script to generate all gff3 files
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

LOGFILE=$APP_HOME/run_weekly.log

log() {
  echo "$(date) $1" | tee -a $LOGFILE
}

log "starting weekly gff3 pipeline run"

log "running create_jbrowse2_files.sh"
bash $APP_HOME/create_jbrowse2_files.sh

log "running clinvar.sh"
bash $APP_HOME/clinvar.sh

log "running agr.sh"
bash $APP_HOME/agr.sh

log "running eva.sh"
bash $APP_HOME/eva.sh

log "running ensemblGff.sh"
bash $APP_HOME/ensemblGff.sh

log "weekly gff3 pipeline run complete"
