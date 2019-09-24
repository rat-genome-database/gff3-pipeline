#!/bin/bash
# load ProteinDomain tracks

JBROWSE_HOME="/rgd/JBrowse-1.16.3/"
GFF3_LOC="/home/rgddata/pipelines/RGDGff3Pipeline/data/ProteinDomain"

cd $JBROWSE_HOME

echo
echo "RAT"
echo

./remove-track.pl --dir /jbrowse/data_rgd6 --trackLabel ProteinDomain --delete

gunzip -c $GFF3_LOC/rat60_domains.gff3.gz > /tmp/rat60_domains.gff3

./flatfile-to-json.pl \
  --gff /tmp/rat60_domains.gff3 \
  --trackLabel ProteinDomain \
  --key "RGD Rat (rn6) Protein Domains" \
  --out /jbrowse/data_rgd6 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }"



