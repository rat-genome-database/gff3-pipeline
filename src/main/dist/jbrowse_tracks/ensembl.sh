# Load Ensembl tracks -- this script is run only manually
#
#RAT TRACKS FROM ENSEMBL
#
#1) Download the gff3 file from Ensembl
#2) Reference sequence in RGD JBrowse has the chromosomes prefixed with 'Chr', f.e. 'Chr1';
#   thus in column one make sure the chromosomes are prefixed with 'Chr'
#3) Remove the chromosome objects (they are not needed in JBrowse)
#4) Split the original file into 2 files:
#   - genes and transcripts model  (grep -v biological_region <file> > Ensembl_rn6_model.gff3)
#   - other features               (grep biological_region <file> > Ensembl_rn6_features.gff3)
#5) Run the script 'ensembl.sh'

cd /rgd/JBrowse-1.12.3/bin

set -e

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
  --gff /rgd/data/gff3/Ensembl/Ensembl_rn6_model.gff3 \
  --trackLabel Ensembl_genes \
  --key "Ensembl (rn6) Genes" \
  --out /jbrowse/data_rgd6 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }" 

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
  --gff /rgd/data/gff3/Ensembl/Ensembl_rn6_features.gff3 \
  --trackLabel Ensembl_features \
  --key "Ensembl (rn6) Features" \
  --out /jbrowse/data_rgd6 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RGD Gene Features\" }"

echo "STEP 5: === DONE ! ==="
