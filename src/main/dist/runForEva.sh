#!/usr/bin/env bash
# shell script to run RGDGff3Pipeline
. /etc/profile

APPNAME=RGDGff3Pipeline
APPDIR=/home/rgddata/pipelines/$APPNAME

cd $APPDIR
java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -Xmx100g -jar lib/$APPNAME.jar "$@"
