#!/bin/bash
# load ProteinDomain tracks

JBROWSE_HOME="/rgd/JBrowse-1.16.3/"
GFF3_LOC="/home/rgddata/pipelines/RGDGff3Pipeline/data/ProteinDomain"

cd $JBROWSE_HOME/bin

echo
echo "RAT 6.0"

./remove-track.pl --dir /jbrowse/data_rgd6 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/rat60_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/rat60_domains.gff3.gz > /tmp/rat60_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/rat60_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Rat (rn6) Protein Domains" \
      --out /jbrowse/data_rgd6 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "RAT 5.0"

./remove-track.pl --dir /jbrowse/data_rgd5 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/rat50_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/rat50_domains.gff3.gz > /tmp/rat50_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/rat50_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Rat (rn5) Protein Domains" \
      --out /jbrowse/data_rgd5 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "RAT 3.4"

./remove-track.pl --dir /jbrowse/data_rgd3_4 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/rat34_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/rat34_domains.gff3.gz > /tmp/rat34_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/rat34_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Rat (rn3.4) Protein Domains" \
      --out /jbrowse/data_rgd3_4 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "HUMAN 38"

./remove-track.pl --dir /jbrowse/data_hg38 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/human38_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/human38_domains.gff3.gz > /tmp/human38_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/human38_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Human (hg38) Protein Domains" \
      --out /jbrowse/data_hg38 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "HUMAN 37"

./remove-track.pl --dir /jbrowse/data_hg19 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/human37_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/human37_domains.gff3.gz > /tmp/human37_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/human37_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Human (hg19) Protein Domains" \
      --out /jbrowse/data_hg19 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "MOUSE 38"

./remove-track.pl --dir /jbrowse/data_mm38 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/mouse38_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/mouse38_domains.gff3.gz > /tmp/mouse38_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/mouse38_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Mouse (mm38) Protein Domains" \
      --out /jbrowse/data_mm38 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "MOUSE 37"

./remove-track.pl --dir /jbrowse/data_mm37 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/mouse37_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/mouse37_domains.gff3.gz > /tmp/mouse37_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/mouse37_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Mouse (mm37) Protein Domains" \
      --out /jbrowse/data_mm37 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "DOG"

./remove-track.pl --dir /jbrowse/data_dog3_1 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/dog31_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/dog31_domains.gff3.gz > /tmp/dog31_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/dog31_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Dog (dog3.1) Protein Domains" \
      --out /jbrowse/data_dog3_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "BONOBO"

./remove-track.pl --dir /jbrowse/data_bonobo1_1 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/bonobo11_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/bonobo11_domains.gff3.gz > /tmp/bonobo11_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/bonobo11_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Bonobo (panpan1.1) Protein Domains" \
      --out /jbrowse/data_bonobo1_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "SQUIRREL"

./remove-track.pl --dir /jbrowse/data_squirrel2_0 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/squirrel20_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/squirrel20_domains.gff3.gz > /tmp/squirrel20_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/squirrel20_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Squirrel (SpeTri2.0) Protein Domains" \
      --out /jbrowse/data_squirrel2_0 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi

echo "CHINCHILLA"

./remove-track.pl --dir /jbrowse/data_cl1_0 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/chinchilla10_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/chinchilla10_domains.gff3.gz > /tmp/chinchilla10_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/chinchilla10_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Chinchilla (ChiLan1.0) Protein Domains" \
      --out /jbrowse/data_cl1_0 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "PIG 11"

./remove-track.pl --dir /jbrowse/data_pig11_1 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/pig11_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/pig11_domains.gff3.gz > /tmp/pig11_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/pig11_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Pig (Sscrofa11.1) Protein Domains" \
      --out /jbrowse/data_pig11_1 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi


echo "PIG 10"

./remove-track.pl --dir /jbrowse/data_pig10_2 --trackLabel ProteinDomain --delete

if [ -f $GFF3_LOC/pig10_domains.gff3.gz ]; then
    gunzip -c $GFF3_LOC/pig10_domains.gff3.gz > /tmp/pig10_domains.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/pig10_domains.gff3 \
      --trackLabel ProteinDomain \
      --key "RGD Pig (Sscrofa10.2) Protein Domains" \
      --out /jbrowse/data_pig10_2 \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --config "{ \"category\" : \"Protein Domains\" }"
fi

echo "=== OK ==="
