#!/usr/bin/env bash
# shell script to run RGDGff3Pipeline
. /etc/profile

APPNAME=RGDGff3Pipeline
APPDIR=/home/rgddata/pipelines/$APPNAME/dist

cd $APPDIR
java -Dspring.config=$APPDIR/../../properties/default_db.xml \
    -Dlog4j.configuration=file://$APPDIR/properties/log4j.properties \
    -jar $APPNAME.jar "$@"
