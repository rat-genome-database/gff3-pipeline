#!/usr/bin/env bash
# shell script to run RGDGff3Pipeline
. /etc/profile

APPNAME=RGDGff3Pipeline
APPDIR=/home/rgddata/pipelines/$APPNAME/dist

cd $APPDIR
pwd
DB_OPTS="-Dspring.config=$APPDIR/../../properties/default_db.xml"
LOG4J_OPTS="-Dlog4j.configuration=file://$APPDIR/properties/log4j.properties"
export RGD_GFF3_PIPELINE_OPTS="$DB_OPTS $LOG4J_OPTS"

bin/$APPNAME "$@"
