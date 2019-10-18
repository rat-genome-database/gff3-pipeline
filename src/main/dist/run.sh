#!/usr/bin/env bash
# shell script to run RGDGff3Pipeline
. /etc/profile

APPNAME=RGDGff3Pipeline
APPDIR=/home/rgddata/pipelines/$APPNAME

cd $APPDIR
java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configuration=file://$APPDIR/properties/log4j.properties \
    -jar lib/$APPNAME.jar "$@"
