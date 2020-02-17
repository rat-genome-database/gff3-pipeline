#!/bin/bash
# load ProteinDomain tracks

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAILLIST=mtutaj@mcw.edu
JBROWSE_HOME="/rgd/JBrowse-1.16.3/"
GFF3_LOC="/home/rgddata/pipelines/RGDGff3Pipeline/data/ProteinDomain"

cd $JBROWSE_HOME/bin

echo
echo "RAT 6.0"

TMP_INPUT_FILE=$GFF3_LOC/Rnor_6.0_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/rat60_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_rgd6 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/rat60_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Rat (rn6) Protein Domains" \
      --out /jbrowse/data_rgd6 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "RAT 5.0"

TMP_INPUT_FILE=$GFF3_LOC/Rnor_5.0_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/rat50_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_rgd5 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/rat50_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Rat (rn5) Protein Domains" \
      --out /jbrowse/data_rgd5 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "RAT 3.4"

TMP_INPUT_FILE=$GFF3_LOC/RGSC_v3.4_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/rat34_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_rgd3_4 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/rat34_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Rat (rn3.4) Protein Domains" \
      --out /jbrowse/data_rgd3_4 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "HUMAN 38"

TMP_INPUT_FILE=$GFF3_LOC/GRCh38_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/human38_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_hg38 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/human38_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Human (hg38) Protein Domains" \
      --out /jbrowse/data_hg38 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "HUMAN 37"

TMP_INPUT_FILE=$GFF3_LOC/GRCh37_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/human37_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_hg19 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/human37_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Human (hg19) Protein Domains" \
      --out /jbrowse/data_hg19 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "MOUSE 38"

TMP_INPUT_FILE=$GFF3_LOC/GRCm38_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/mouse38_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_mm38 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/mouse38_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Mouse (mm38) Protein Domains" \
      --out /jbrowse/data_mm38 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "MOUSE 37"

TMP_INPUT_FILE=$GFF3_LOC/MGSCv37_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/mouse37_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_mm37 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/mouse37_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Mouse (mm37) Protein Domains" \
      --out /jbrowse/data_mm37 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "DOG"

TMP_INPUT_FILE=$GFF3_LOC/CanFam3.1_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/dog31_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_dog3_1 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/dog31_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Dog (dog3.1) Protein Domains" \
      --out /jbrowse/data_dog3_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "BONOBO"

TMP_INPUT_FILE=$GFF3_LOC/PanPan1.1_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/bonobo11_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_bonobo1_1 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/bonobo11_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Bonobo (panpan1.1) Protein Domains" \
      --out /jbrowse/data_bonobo1_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "SQUIRREL"

TMP_INPUT_FILE=$GFF3_LOC/SpeTri2.0_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/squirrel20_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_squirrel2_0 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/squirrel20_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Squirrel (SpeTri2.0) Protein Domains" \
      --out /jbrowse/data_squirrel2_0 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "CHINCHILLA"

TMP_INPUT_FILE=$GFF3_LOC/ChiLan1.0_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/chinchilla10_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_cl1_0 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/chinchilla10_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Chinchilla (ChiLan1.0) Protein Domains" \
      --out /jbrowse/data_cl1_0 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi


echo "PIG Sscrofa11.1"

TMP_INPUT_FILE=$GFF3_LOC/Sscrofa11.1_domains.gff3.gz
if [ -f $GFF3_LOC/TMP_INPUT_FILE ]; then
    gunzip -c $GFF3_LOC/TMP_INPUT_FILE > /tmp/pig11_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_pig11_1 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/pig11_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Pig (Sscrofa11.1) Protein Domains" \
      --out /jbrowse/data_pig11_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

### SO FAR, NO DATA FOR PIG Sscrofa10.2, Vervet 1.1 and HetGla 1.0
echo "=== OK ==="
exit 0



echo "PIG Sscrofa10.2"

TMP_INPUT_FILE=$GFF3_LOC/Sscrofa10.2_domains.gff3.gz
if [ -f $TMP_INPUT_FILE ]; then
    gunzip -c $TMP_INPUT_FILE > /tmp/pig10_domains.gff3

    ./remove-track.pl --dir /jbrowse/data_pig10_2 --trackLabel ProteinDomain --delete

    ./flatfile-to-json.pl \
      --gff /tmp/pig10_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Pig (Sscrofa10.2) Protein Domains" \
      --out /jbrowse/data_pig10_2 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
else
    echo "ERROR: File not found: $TMP_INPUT_FILE" |  mailx -s "[$SERVER] GFF3 JBrowse Loader: missing file $TMP_INPUT_FILE" $EMAILLIST
fi

echo "=== OK ==="
