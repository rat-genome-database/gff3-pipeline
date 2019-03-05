# called from master script 'rnacentral.sh' to do the actual refresh of RNAcentral track for given species
#
# parameters:
# 1 - name of file to be downloaded from RNAcentral (without .gz extension)
# 2 - JBrowse data directory for species

RNACENTRAL_URL=ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/gff3
JBROWSE_BIN=/rgd/JBrowse-1.12.3/bin

FILEURL=${RNACENTRAL_URL}/$1.gz
TRACKNAME=RNAcentral_features

echo "==="
echo "processing [$1]       [$2]"

cd /data/data/gff3/rnacentral
wget $FILEURL
gunzip $1.gz

#gff3 files from RNAcentral have chromosomes as '1','2' -- for RGD JBrowse installation, we need to prefix them with 'Chr', f.e. 'Chr1',...
cp -p $1 original_$1
grep "^#" $1 > c_$1
grep -v "^#" $1 > d_$1
sed  -e "s/^/Chr/" d_$1 > f_$1
cat c_$1 f_$1 > $1
rm c_$1 d_$1 f_$1

echo
echo "STEP 1: DELETE TRACK $TRACKNAME"
echo

$JBROWSE_BIN/remove-track.pl \
  --dir $2 \
  --trackLabel $TRACKNAME \
  --delete

echo
echo "STEP 2: LOAD TRACK $TRACKNAME"
echo

$JBROWSE_BIN/flatfile-to-json.pl \
  --gff $1 \
  --trackLabel $TRACKNAME \
  --key "RNAcentral Features" \
  --out $2 \
  --trackType JBrowse/View/Track/CanvasFeatures \
  --config "{ \"category\" : \"Gene Models/RNAcentral Gene Features\" }"

echo "=== OK! ==="
