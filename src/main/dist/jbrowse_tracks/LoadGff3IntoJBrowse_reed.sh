# script env setup
. /etc/profile

# variables setup
JBROWSE="JBrowse_RGDtracks"

pushd /rgd/data/gff3
bash ${JBROWSE}_updater5.sh \
  ${JBROWSE}_updaterInput.conf \
  ${JBROWSE}_updaterOntInput.conf &> ${JBROWSE}_updater.log

#generate ClinVar tracks for human hg38, hg19 and hg18
./clinvar.sh 2>&1 > clinvar.log

popd

