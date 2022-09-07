#!/bin/bash
# load marker tracks
#
# NOTE: this script is not finished! it is a start ...

JBROWSE_HOME="/rgd/JBrowse-1.16.11/"
GFF3_LOC="/home/rgddata/pipelines/RGDGff3Pipeline/data/Marker"

cd $JBROWSE_HOME/bin

echo
echo "RAT 7.2"

DATADIR=/jbrowse/data_rn7_2

./remove-track.pl --dir $DATADIR --trackLabel SSLP --delete

if [ -f $GFF3_LOC/Rat/rn7_markers.gff3.gz ]; then
    gunzip -c $GFF3_LOC/Rat/rn7_markers.gff3.gz > /tmp/rn7_markers.gff3

    ./flatfile-to-json.pl \
      --gff /tmp/rn7_markers.gff3 \
      --trackLabel SSLP \
      --key "RGD mRatBN7.2 Markers" \
      --out $DATADIR \
      --trackType JBrowse/View/Track/CanvasFeatures \
      --clientConfig "{ \"color\" : \"#378D15\", \"featureScale\" : 0.0001, \"label\" : \"symbol,name\", \"description\" : \"\" }" \
      --config "{ \"category\" : \"Markers/RGD mRatBN7.2 Markers\", \"histograms\" : { \"color\" : \"#378D15\", \"binsPerBlock\" : 25 }, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"

			./bin/generate-names.pl --out $DATADIR --tracks "SSLP" --verbose --mem 4096000000

fi
echo "=== OK ==="
