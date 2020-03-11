#!/bin/bash
# load EVA tracks

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAILLIST=mtutaj@mcw.edu,llamers@mcw.edu
JBROWSE_HOME="/rgd/JBrowse-1.12.3/"
GFF3_LOC="/rgd/data/gff3/Eva"
#"/home/rgddata/pipelines/RGDGff3Pipeline/data/Eva"

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "HANSEN" ]; then
    cd /rgd/data/gff3/Eva
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Eva/EVA_Rnor_6.0.gff3.gz .
    cd ../../../../..
fi

cd $JBROWSE_HOME/bin

echo
echo "RAT 6.0"

TMP_INPUT_FILE=$GFF3_LOC/EVA_Rnor_6.0.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/rat60_EVA.gff3

    ./remove-track.pl --dir /jbrowse/data_rgd6 --trackLabel EVA --delete

    ./flatfile-to-json.pl \
      --gff /tmp/rat60_EVA.gff3 \
      --trackLabel EVA \
      --key "RGD Rat (rn6) EVA" \
      --out /jbrowse/data_rgd6 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Variants\DbSNPs\EVA\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "=== OK ==="