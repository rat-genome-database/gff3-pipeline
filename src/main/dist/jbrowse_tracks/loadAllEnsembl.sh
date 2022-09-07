#!/bin/bash
# load Ensembl tracks

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAILLIST=mtutaj@mcw.edu,llamers@mcw.edu

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "HANSEN" ]; then
    JBROWSE_HOME="/rgd/JBrowse-1.12.3/"
    GFF3_LOC="/rgd/data/gff3/Ensembl"
    cd /rgd/data/gff3/Ensembl
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/Rnor_6.0_Ensembl-model.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/Rnor_6.0_Ensembl-feature.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/GRCh38.p13_Ensembl-model.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/GRCh38.p13_Ensembl-feature.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/GRCm38.p6_Ensembl-model.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/GRCm38.p6_Ensembl-feature.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/Sscrofa11.1_Ensembl-model.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/Sscrofa11.1_Ensembl-feature.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/PanPan1.1_Ensembl-model.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/PanPan1.1_Ensembl-feature.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/CanFam3.1_Ensembl-model.gff3.gz .
    scp -p rgddata@travis.rgd.mcw.edu:/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl/CanFam3.1_Ensembl-feature.gff3.gz .
    cd ../../../../..
else
    JBROWSE_HOME="/rgd/JBrowse-1.16.11/"
    GFF3_LOC="/home/rgddata/pipelines/RGDGff3Pipeline/data/Ensembl"
fi

cd $JBROWSE_HOME/bin


set -e
TMP_INPUT_FILE=$GFF3_LOC/Rnor_6.0_Ensembl-model.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/Rnor_6.0_Ensembl-model.gff3
    echo "Rnor 6.0"
    echo
    echo "STEP 1: DELETE TRACK Ensembl_genes"
    echo


    ./remove-track.pl \
      --dir /jbrowse/data_rgd6 \
      --trackLabel Ensembl_genes \
      --delete

    echo
    echo "STEP 2: LOAD TRACK Ensembl_genes"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/Rnor_6.0_Ensembl-model.gff3 \
      --trackLabel Ensembl_genes \
      --key "Ensembl (rn6) Genes and Transcripts" \
      --out /jbrowse/data_rgd6 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

TMP_INPUT_FILE=$GFF3_LOC/Rnor_6.0_Ensembl-feature.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/Rnor_6.0_Ensembl-feature.gff3
    echo
    echo "STEP 3: DELETE TRACK Ensembl_features"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_rgd6 \
      --trackLabel Ensembl_features \
      --delete

    echo
    echo "STEP 4: LOAD TRACK Ensembl_features"
    echo

    ./flatfile-to-json.pl \
      --gff /ref/gff3_ensembl/Ensembl_rn6_features.gff3 \
      --trackLabel Ensembl_features \
      --key "Ensembl (rn6) Features" \
      --out /jbrowse/data_rgd6 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"

    echo "STEP 5: === DONE ! ==="
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

echo "Human"

TMP_INPUT_FILE=$GFF3_LOC/GRCh38.p13_Ensembl-model.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/GRCh38.p13_Ensembl-model.gff3
    echo
    echo "STEP 1: DELETE TRACK Ensembl_genes"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_hg38 \
      --trackLabel Ensembl_genes \
      --delete

    echo
    echo "STEP 2: LOAD TRACK Ensembl_genes"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/GRCh38.p13_Ensembl-model.gff3 \
      --trackLabel Ensembl_genes \
      --key "Ensembl (hg38) Genes and Transcripts" \
      --out /jbrowse/data_hg38 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

TMP_INPUT_FILE=$GFF3_LOC/GRCh38.p13_Ensembl-feature.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/GRCh38.p13_Ensembl-feature.gff3
    echo
    echo "STEP 3: DELETE TRACK Ensembl_features"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_hg38 \
      --trackLabel Ensembl_features \
      --delete

    echo
    echo "STEP 4: LOAD TRACK Ensembl_features"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/GRCh38.p13_Ensembl-feature.gff3 \
      --trackLabel Ensembl_features \
      --key "Ensembl (hg38) Features" \
      --out /jbrowse/data_hg38 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"

    echo "STEP 5: === DONE ! ==="

else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

echo "Mouse"

TMP_INPUT_FILE=$GFF3_LOC/GRCm38.p6_Ensembl-model.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/GRCm38.p6_Ensembl-model.gff3
    echo
    echo "STEP 1: DELETE TRACK Ensembl_genes"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_mm38 \
      --trackLabel Ensembl_genes \
      --delete

    echo
    echo "STEP 2: LOAD TRACK Ensembl_genes"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/GRCm38.p6_Ensembl-model.gff3 \
      --trackLabel Ensembl_genes \
      --key "Ensembl (mm38) Genes and Transcripts" \
      --out /jbrowse/data_mm38 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"

else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

TMP_INPUT_FILE=$GFF3_LOC/GRCm38.p6_Ensembl-feature.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/GRCm38.p6_Ensembl-feature.gff3
    echo
    echo "STEP 3: DELETE TRACK Ensembl_features"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_mm38 \
      --trackLabel Ensembl_features \
      --delete

    echo
    echo "STEP 4: LOAD TRACK Ensembl_features"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/GRCm38.p6_Ensembl-feature.gff3 \
      --trackLabel Ensembl_features \
      --key "Ensembl (mm38) Features" \
      --out /jbrowse/data_mm38 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"

    echo "STEP 5: === DONE ! ==="
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

echo "Pig"

TMP_INPUT_FILE=$GFF3_LOC/Sscrofa11.1_Ensembl-model.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/Sscrofa11.1_Ensembl-model.gff3
    echo
    echo "STEP 1: DELETE TRACK Ensembl_genes"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_pig11_1 \
      --trackLabel Ensembl_genes \
      --delete

    echo
    echo "STEP 2: LOAD TRACK Ensembl_genes"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/Sscrofa11.1_Ensembl-model.gff3 \
      --trackLabel Ensembl_genes \
      --key "Ensembl (Sscrofa11.1) Genes and Transcripts" \
      --out /jbrowse/data_pig11_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"

else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

TMP_INPUT_FILE=$GFF3_LOC/Sscrofa11.1_Ensembl-feature.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/Sscrofa11.1_Ensembl-feature.gff3
    echo
    echo "STEP 3: DELETE TRACK Ensembl_features"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_pig11_1 \
      --trackLabel Ensembl_features \
      --delete

    echo
    echo "STEP 4: LOAD TRACK Ensembl_features"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/Sscrofa11.1_Ensembl-feature.gff3 \
      --trackLabel Ensembl_features \
      --key "Ensembl (Sscrofa11.1) Features" \
      --out /jbrowse/data_pig11_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"

    echo "STEP 5: === DONE ! ==="

else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

echo "Bonobo"

TMP_INPUT_FILE=$GFF3_LOC/PanPan1.1_Ensembl-model.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/PanPan1.1_Ensembl-model.gff3
    echo
    echo "STEP 1: DELETE TRACK Ensembl_genes"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_bonobo1_1 \
      --trackLabel Ensembl_genes \
      --delete

    echo
    echo "STEP 2: LOAD TRACK Ensembl_genes"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/PanPan1.1_Ensembl-model.gff3 \
      --trackLabel Ensembl_genes \
      --key "Ensembl (panpan1.1) Genes and Transcripts" \
      --out /jbrowse/data_bonobo1_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"

else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

TMP_INPUT_FILE=$GFF3_LOC/PanPan1.1_Ensembl-feature.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/PanPan1.1_Ensembl-feature.gff3
    echo
    echo "STEP 3: DELETE TRACK Ensembl_features"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_bonobo1_1 \
      --trackLabel Ensembl_features \
      --delete

    echo
    echo "STEP 4: LOAD TRACK Ensembl_features"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/PanPan1.1_Ensembl-feature.gff3 \
      --trackLabel Ensembl_features \
      --key "Ensembl (panpan1.1) Features" \
      --out /jbrowse/data_bonobo1_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"

    echo "STEP 5: === DONE ! ==="

else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

echo "Dog"

TMP_INPUT_FILE=$GFF3_LOC/CanFam3.1_Ensembl-model.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/CanFam3.1_Ensembl-model.gff3
    echo
    echo "STEP 1: DELETE TRACK Ensembl_genes"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_dog3_1 \
      --trackLabel Ensembl_genes \
      --delete

    echo
    echo "STEP 2: LOAD TRACK Ensembl_genes"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/CanFam3.1_Ensembl-model.gff3 \
      --trackLabel Ensembl_genes \
      --key "Ensembl (dog3.1) Genes and Transcripts" \
      --out /jbrowse/data_dog3_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"

else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

TMP_INPUT_FILE=$GFF3_LOC/CanFam3.1_Ensembl-feature.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/CanFam3.1_Ensembl-feature.gff3
    echo
    echo "STEP 3: DELETE TRACK Ensembl_features"
    echo

    ./remove-track.pl \
      --dir /jbrowse/data_dog3_1 \
      --trackLabel Ensembl_features \
      --delete

    echo
    echo "STEP 4: LOAD TRACK Ensembl_features"
    echo

    ./flatfile-to-json.pl \
      --gff /tmp/CanFam3.1_Ensembl-feature.gff3 \
      --trackLabel Ensembl_features \
      --key "Ensembl (dog3.1) Features" \
      --out /jbrowse/data_dog3_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Gene Models/Ensembl Gene Features\" }"

    echo "STEP 5: === DONE ! ==="
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" #|  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi