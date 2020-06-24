#!/bin/bash
# load EVA tracks

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAILLIST=mtutaj@mcw.edu,llamers@mcw.edu
JBROWSE_HOME="/rgd/JBrowse-1.12.3/"

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "HANSEN" ]; then
    GFF3_LOC="/rgd/data/gff3/Eva"
    cd /rgd/data/gff3/Eva
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Eva/EVA_Rnor_6.0.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Eva/EVA_Rnor_5.0.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Eva/EVA_GRCm38.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Eva/EVA_CanFam3.1.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Eva/EVA_Sscrofa11.1.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Eva/EVA_Sscrofa10.2.gff3.gz .
    cd ../../../../..
else
    GFF3_LOC="/home/rgddata/pipelines/RGDGff3Pipeline/data/Eva"
fi

cd $JBROWSE_HOME/bin

echo
#echo "RAT 6.0"
#
#TMP_INPUT_FILE=$GFF3_LOC/EVA_Rnor_6.0.gff3.gz
#if [ -f $TMP_INPUT_FILE ]; then
#    gunzip -c $TMP_INPUT_FILE > /tmp/rat60_EVA.gff3
#
#    ./remove-track.pl --dir /jbrowse/data_rgd6 --trackLabel EVA --delete
#
#    ./flatfile-to-json.pl \
#      --gff /tmp/rat60_EVA.gff3 \
#      --trackLabel EVA \
#      --key "EVA" \
#      --out /jbrowse/data_rgd6 \
#      --trackType JBrowse/View/Track/CanvasFeatures \
#      --config "{ \"category\" : \"Variants/DbSNPs\" }"
#else
#    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
#fi
#
#
#echo
#echo "RAT 5.0"
#
#TMP_INPUT_FILE=$GFF3_LOC/EVA_Rnor_5.0.gff3.gz
#if [ -f $TMP_INPUT_FILE ]; then
#    gunzip -c $TMP_INPUT_FILE > /tmp/rat50_EVA.gff3
#
#    ./remove-track.pl --dir /jbrowse/data_rgd5 --trackLabel EVA --delete
#
#    ./flatfile-to-json.pl \
#      --gff /tmp/rat50_EVA.gff3 \
#      --trackLabel EVA \
#      --key "EVA" \
#      --out /jbrowse/data_rgd5 \
#      --trackType JBrowse/View/Track/CanvasFeatures \
#      --config "{ \"category\" : \"Variants/DbSNPs\" }"
#else
#    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
#fi
#
#
#echo
#echo "MOUSE 38"
#
#TMP_INPUT_FILE=$GFF3_LOC/EVA_GRCm38.gff3.gz
#if [ -f $TMP_INPUT_FILE ]; then
#    gunzip -c $TMP_INPUT_FILE > /tmp/mouse38_EVA.gff3
#
#    ./remove-track.pl --dir /jbrowse/data_mm38 --trackLabel EVA --delete
#
#    ./flatfile-to-json.pl \
#      --gff /tmp/mouse38_EVA.gff3 \
#      --trackLabel EVA \
#      --key "EVA" \
#      --out /jbrowse/data_mm38 \
#      --trackType JBrowse/View/Track/CanvasFeatures \
#      --config "{ \"category\" : \"Variants/DbSNPs\" }"
#else
#    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
#fi


echo
echo "DOG"

TMP_INPUT_FILE=$GFF3_LOC/EVA_CanFam3.1.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/dog31_EVA.gff3

    ./remove-track.pl --dir /jbrowse/data_dog3_1 --trackLabel EVA --delete

    ./flatfile-to-json.pl | tee -a EvaOutput.log \
      --gff /tmp/dog31_EVA.gff3 \
      --trackLabel EVA \
      --key "EVA" \
      --out /jbrowse/data_dog3_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Variants/DbSNPs\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi
#
#
#echo
#echo "PIG 11.1"
#
#TMP_INPUT_FILE=$GFF3_LOC/EVA_Sscrofa11.1.gff3.gz
#if [ -f $TMP_INPUT_FILE ]; then
#    gunzip -c $TMP_INPUT_FILE > /tmp/pig11_EVA.gff3
#
#    ./remove-track.pl --dir /jbrowse/data_pig11_1 --trackLabel EVA --delete
#
#    ./flatfile-to-json.pl \
#      --gff /tmp/pig11_EVA.gff3 \
#      --trackLabel EVA \
#      --key "EVA" \
#      --out /jbrowse/data_pig11_1 \
#      --trackType JBrowse/View/Track/CanvasFeatures \
#      --config "{ \"category\" : \"Variants/DbSNPs\" }"
#else
#    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
#fi


echo
echo "PIG 10.2"

TMP_INPUT_FILE=$GFF3_LOC/EVA_Sscrofa10.2.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/pig10_EVA.gff3

    ./remove-track.pl --dir /jbrowse/data_pig10_2 --trackLabel EVA --delete

    ./flatfile-to-json.pl | tee -a EvaOutput.log \
      --gff /tmp/pig10_EVA.gff3 \
      --trackLabel EVA \
      --key "EVA" \
      --out /jbrowse/data_pig10_2 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Variants/DbSNPs\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "=== OK ==="