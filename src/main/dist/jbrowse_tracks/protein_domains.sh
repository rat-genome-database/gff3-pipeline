#!/bin/bash
# load ProteinDomain tracks

JBROWSE_HOME="/rgd/JBrowse-1.16.3/"
GFF3_LOC="/home/rgddata/pipelines/RGDGff3Pipeline/data/ProteinDomain"

cd $JBROWSE_HOME/bin

echo
echo "RAT"

./remove-track.pl --dir /jbrowse/data_rgd6 --trackLabel ProteinDomain --delete

gunzip -c $GFF3_LOC/rat60_domains.gff3.gz > /tmp/rat60_domains.gff3

./flatfile-to-json.pl \
  --gff /tmp/rat60_domains.gff3 \
  --trackLabel ProteinDomain \
  --key "RGD Rat (rn6) Protein Domains" \
  --out /jbrowse/data_rgd6 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }"


echo "HUMAN"

./remove-track.pl --dir /jbrowse/data_hg38 --trackLabel ProteinDomain --delete

gunzip -c $GFF3_LOC/human38_domains.gff3.gz > /tmp/human38_domains.gff3

./flatfile-to-json.pl \
  --gff /tmp/human38_domains.gff3 \
  --trackLabel ProteinDomain \
  --key "RGD Human (hg38) Protein Domains" \
  --out /jbrowse/data_hg38 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }"


echo "MOUSE"

./remove-track.pl --dir /jbrowse/data_mm38 --trackLabel ProteinDomain --delete

gunzip -c $GFF3_LOC/mouse38_domains.gff3.gz > /tmp/mouse38_domains.gff3

./flatfile-to-json.pl \
  --gff /tmp/mouse38_domains.gff3 \
  --trackLabel ProteinDomain \
  --key "RGD Mouse (mm38) Protein Domains" \
  --out /jbrowse/data_mm38 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }"


echo "DOG"

./remove-track.pl --dir /jbrowse/data_dog3_1 --trackLabel ProteinDomain --delete

gunzip -c $GFF3_LOC/dog31_domains.gff3.gz > /tmp/dog31_domains.gff3

./flatfile-to-json.pl \
  --gff /tmp/dog31_domains.gff3 \
  --trackLabel ProteinDomain \
  --key "RGD Dog (dog3.1) Protein Domains" \
  --out /jbrowse/data_dog3_1 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }"


echo "BONOBO"

./remove-track.pl --dir /jbrowse/data_bonobo1_1 --trackLabel ProteinDomain --delete

gunzip -c $GFF3_LOC/bonobo11_domains.gff3.gz > /tmp/bonobo11_domains.gff3

./flatfile-to-json.pl \
  --gff /tmp/bonobo11_domains.gff3 \
  --trackLabel ProteinDomain \
  --key "RGD Bonobo (panpan1.1) Protein Domains" \
  --out /jbrowse/data_bonobo1_1 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }"


echo "SQUIRREL"

./remove-track.pl --dir /jbrowse/data_squirrel2_0 --trackLabel ProteinDomain --delete

gunzip -c $GFF3_LOC/squirrel20_domains.gff3.gz > /tmp/squirrel20_domains.gff3

./flatfile-to-json.pl \
  --gff /tmp/squirrel20_domains.gff3 \
  --trackLabel ProteinDomain \
  --key "RGD Squirrel (SpeTri2.0) Protein Domains" \
  --out /jbrowse/data_squirrel2_0 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }"


echo "CHINCHILLA"

./remove-track.pl --dir /jbrowse/data_cl1_0 --trackLabel ProteinDomain --delete

gunzip -c $GFF3_LOC/chinchilla10_domains.gff3.gz > /tmp/chinchilla10_domains.gff3

./flatfile-to-json.pl \
  --gff /tmp/chinchilla10_domains.gff3 \
  --trackLabel ProteinDomain \
  --key "RGD Chinchilla (ChiLan1.0) Protein Domains" \
  --out /jbrowse/data_cl1_0 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }"


echo "PIG"

./remove-track.pl --dir /jbrowse/data_pig11_1 --trackLabel ProteinDomain --delete

gunzip -c $GFF3_LOC/pig11_domains.gff3.gz > /tmp/pig11_domains.gff3

./flatfile-to-json.pl \
  --gff /tmp/pig11_domains.gff3 \
  --trackLabel ProteinDomain \
  --key "RGD Pig (Sscrofa11.1) Protein Domains" \
  --out /jbrowse/data_pig11_1 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }"


echo "=== OK ==="
