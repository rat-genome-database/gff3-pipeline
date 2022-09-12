#!/bin/bash
# load Promoters tracks

JBROWSE_HOME="/rgd/JBrowse-1.16.11/"
GFF3_LOC="/home/rgddata/pipelines/RGDGff3Pipeline/data/Promoter"

cd $JBROWSE_HOME/bin

echo
echo "RAT 6.0"

./remove-track.pl --dir /jbrowse/data_rgd6 --trackLabel Promoters --delete

if [ -f $GFF3_LOC/Rat/Rnor_6.0_promoters.gff3.gz ]; then
    gunzip -c $GFF3_LOC/Rat/Rnor_6.0_promoters.gff3.gz > /tmp/rat60_promoters.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/rat60_promoters.gff3 \
      --trackLabel Promoters \
      --key "EPD (rn6) Promoters" \
      --out /jbrowse/data_rgd6 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Promoters\" }"
fi


echo "RAT 3.4"

./remove-track.pl --dir /jbrowse/data_rgd3_4 --trackLabel Promoters --delete

if [ -f $GFF3_LOC/Rat/RGSC_v3.4_promoters.gff3.gz ]; then
    gunzip -c $GFF3_LOC/Rat/RGSC_v3.4_promoters.gff3.gz > /tmp/rat34_promoters.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/rat34_promoters.gff3 \
      --trackLabel Promoters \
      --key "EPD (rn3.4) Promoters" \
      --out /jbrowse/data_rgd3_4 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Promoters\" }"
fi


echo "HUMAN 38"

./remove-track.pl --dir /jbrowse/data_hg38 --trackLabel Promoters --delete

if [ -f $GFF3_LOC/Human/GRCh38_promoters.gff3.gz ]; then
    gunzip -c $GFF3_LOC/Human/GRCh38_promoters.gff3.gz > /tmp/human38_promoters.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/human38_promoters.gff3 \
      --trackLabel Promoters \
      --key "EPD (hg38) Promoters" \
      --out /jbrowse/data_hg38 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Promoters\" }"
fi


echo "HUMAN 36"

./remove-track.pl --dir /jbrowse/data_hg18 --trackLabel Promoters --delete

if [ -f "${GFF3_LOC}/Human/Build 36_promoters.gff3.gz" ]; then
    gunzip -c "${GFF3_LOC}/Human/Build 36_promoters.gff3.gz" > /tmp/human36_promoters.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/human36_promoters.gff3 \
      --trackLabel Promoters \
      --key "EPD (hg18) Promoters" \
      --out /jbrowse/data_hg18 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Promoters\" }"
fi


echo "MOUSE 38"

./remove-track.pl --dir /jbrowse/data_mm38 --trackLabel Promoters --delete

if [ -f $GFF3_LOC/Mouse/GRCm38_promoters.gff3.gz ]; then
    gunzip -c $GFF3_LOC/Mouse/GRCm38_promoters.gff3.gz > /tmp/mouse38_promoters.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/mouse38_promoters.gff3 \
      --trackLabel Promoters \
      --key "EPD (mm38) Promoters" \
      --out /jbrowse/data_mm38 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Promoters\" }"
fi


echo "MOUSE 37"

./remove-track.pl --dir /jbrowse/data_mm37 --trackLabel Promoters --delete

if [ -f $GFF3_LOC/Mouse/MGSCv37_promoters.gff3.gz ]; then
    gunzip -c $GFF3_LOC/Mouse/MGSCv37_promoters.gff3.gz > /tmp/mouse37_promoters.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/mouse37_promoters.gff3 \
      --trackLabel Promoters \
      --key "EPD (mm37) Promoters" \
      --out /jbrowse/data_mm37 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Promoters\" }"
fi


echo "DOG"

./remove-track.pl --dir /jbrowse/data_dog3_1 --trackLabel Promoters --delete

if [ -f $GFF3_LOC/Dog/CanFam3.1_promoters.gff3.gz ]; then
    gunzip -c $GFF3_LOC/Dog/CanFam3.1_promoters.gff3.gz > /tmp/dog31_promoters.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/dog31_promoters.gff3 \
      --trackLabel Promoters \
      --key "EPD (dog3.1) Promoters" \
      --out /jbrowse/data_dog3_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Promoters\" }"
fi


echo "=== OK ==="
