# script env setup
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
RUNLOAD="$APP_HOME/run.sh"

$RUNLOAD -object:phenotypic_variants > phenotypic_variants.log
$RUNLOAD -object:diseases > diseases.log
$RUNLOAD -object:chebi > chebi.log
$RUNLOAD -object:genes > genes.log
$RUNLOAD -object:qtls > qtls.log
$RUNLOAD -object:markers > markers.log
$RUNLOAD -object:strains > strains.log
$RUNLOAD -object:proteinDomains > domains.log
$RUNLOAD -object:ensembl > ensembl.log
$RUNLOAD -object:variants > variants.log
$RUNLOAD -object:jb2_eva > jb2_eva.log

# uth files
mkdir "data/jbrowse2/Rat/UTH_Rnor_SHR_Utx/Gene Models"
scp -p /ref/gff3/uth_shr2.gff3.gz "data/jbrowse2/Rat/UTH_Rnor_SHR_Utx/Gene Models/UTH_Rnor_SHR_Utx Genes and Transcripts.gff3.gz"

mkdir "data/jbrowse2/Rat/UTH_Rnor_SHRSP_BbbUtx_1.0/Gene Models"
scp -p /ref/gff3/uth_shrsp2.gff3.gz "data/jbrowse2/Rat/UTH_Rnor_SHRSP_BbbUtx_1.0/Gene Models/UTH_Rnor_SHRSP_BbbUtx_1.0 Genes and Transcripts.gff3.gz"

mkdir "data/jbrowse2/Rat/UTH_Rnor_WKY_Bbb_1.0/Gene Models"
scp -p /ref/gff3/uth_wky2.gff3.gz "data/jbrowse2/Rat/UTH_Rnor_WKY_Bbb_1.0/Gene Models/UTH_Rnor_WKY_Bbb_1.0 Genes and Transcripts.gff3.gz"