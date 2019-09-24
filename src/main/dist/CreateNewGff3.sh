# script env setup
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

CHR_RAT="1-20,X,Y,MT"
CHR_MOUSE="1-19,X,Y,MT"
CHR_HUMAN="1-22,X,Y,MT"
CHR_DOG="1-38,X,MT"
CHR_PIG="1-18,X,Y,MT"
CHR_BONOBO="1,2A,2B,3-22,X,MT"
CHR_SQUIRREL="Scaffold"
CHR_CHINCHILLA="Scaffold"

LOGDIR=$APP_HOME/data
DATA_RELEASE_DIR_GFF3=/home/rgddata/data_release/GFF3
DATA_RELEASE_DIR_GFF=/home/rgddata/data_release/GFF


SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" == "REED" ]; then
#  EMAILLIST=jrsmith@mcw.edu,RGD.Developers@mcw.edu
  EMAILLIST=mtutaj@mcw.edu
else
  EMAILLIST=mtutaj@mcw.edu
fi

RUNLOAD="$APP_HOME/run.sh"

##### PROTEIN DOMAINS -- currently there is data only for primary assembly

$RUNLOAD -object:proteinDomain -species:PIG -mapKey:911 -toFile:$LOGDIR/ProteinDomain/pig11_domains.gff3 -compress  &> cron_pig11Domains.log
$RUNLOAD -object:proteinDomain -species:CHINCHILLA -mapKey:44 -toFile:$LOGDIR/ProteinDomain/chinchilla10_domains.gff3 -compress  &> cron_chinchilla10Domains.log
$RUNLOAD -object:proteinDomain -species:SQUIRREL -mapKey:720 -toFile:$LOGDIR/ProteinDomain/squirrel20_domains.gff3 -compress  &> cron_squirrel20Domains.log
$RUNLOAD -object:proteinDomain -species:BONOBO -mapKey:511 -toFile:$LOGDIR/ProteinDomain/bonobo11_domains.gff3 -compress  &> cron_bonobo31Domains.log
$RUNLOAD -object:proteinDomain -species:DOG -mapKey:631 -toFile:$LOGDIR/ProteinDomain/dog31_domains.gff3 -compress  &> cron_dog31Domains.log
$RUNLOAD -object:proteinDomain -species:RAT -mapKey:360 -toFile:$LOGDIR/ProteinDomain/rat60_domains.gff3 -compress  &> cron_rat60Domains.log
$RUNLOAD -object:proteinDomain -species:HUMAN -mapKey:38 -toFile:$LOGDIR/ProteinDomain/human38_domains.gff3 -compress  &> cron_human38Domains.log
$RUNLOAD -object:proteinDomain -species:MOUSE -mapKey:35 -toFile:$LOGDIR/ProteinDomain/mouse38_domains.gff3 -compress  &> cron_mouse38Domains.log

mailx -s "[$SERVER]Pipeline to create Gff3 data for protein domains ran" $EMAILLIST < logs/domains.log


##### GENES

$RUNLOAD -object:gene -species:PIG -mapKey:911 -toFile:$LOGDIR/Gene/Pig/pig11/ -chr:$CHR_PIG -compress  &> cron_pigGene11.log
$RUNLOAD -object:gene -species:PIG -mapKey:910 -toFile:$LOGDIR/Gene/Pig/pig10/ -chr:$CHR_PIG -compress  &> cron_pigGene10.log
$RUNLOAD -object:gene -species:CHINCHILLA -mapKey:44 -toFile:$LOGDIR/Gene/Chinchilla/chinchilla10/ -chr:$CHR_CHINCHILLA -compress  &> cron_chinchillaGene10.log
$RUNLOAD -object:gene -species:SQUIRREL -mapKey:720 -toFile:$LOGDIR/Gene/Squirrel/squirrel20/ -chr:$CHR_SQUIRREL -compress  &> cron_squirrelGene20.log
$RUNLOAD -object:gene -species:BONOBO -mapKey:511 -toFile:$LOGDIR/Gene/Bonobo/bonobo11/ -chr:$CHR_BONOBO -compress  &> cron_bonoboGene31.log
$RUNLOAD -object:gene -species:DOG -mapKey:631 -toFile:$LOGDIR/Gene/Dog/dog31/ -chr:$CHR_DOG -compress  &> cron_dogGene31.log

$RUNLOAD -object:gene -species:RAT -mapKey:60 -toFile:$LOGDIR/Gene/Rat/rat34/ -chr:$CHR_RAT -compress  &> cron_ratGene34.log
$RUNLOAD -object:gene -species:RAT -mapKey:70 -toFile:$LOGDIR/Gene/Rat/rat50/ -chr:$CHR_RAT -compress  &> cron_ratGene50.log
$RUNLOAD -object:gene -species:RAT -mapKey:360 -toFile:$LOGDIR/Gene/Rat/rat60/ -chr:$CHR_RAT -compress  &> cron_ratGene60.log

$RUNLOAD -object:gene -species:HUMAN -mapKey:13 -toFile:$LOGDIR/Gene/Human/human36/ -chr:$CHR_HUMAN -compress  &> cron_humanGene36.log
$RUNLOAD -object:gene -species:HUMAN -mapKey:17 -toFile:$LOGDIR/Gene/Human/human37/ -chr:$CHR_HUMAN -compress  &> cron_humanGene37.log
$RUNLOAD -object:gene -species:HUMAN -mapKey:38 -toFile:$LOGDIR/Gene/Human/human38/ -chr:$CHR_HUMAN -compress  &> cron_humanGene38.log

$RUNLOAD -object:gene -species:MOUSE -mapKey:18 -toFile:$LOGDIR/Gene/Mouse/mouse37/ -chr:$CHR_MOUSE -compress  &> cron_mouseGene37.log
$RUNLOAD -object:gene -species:MOUSE -mapKey:35 -toFile:$LOGDIR/Gene/Mouse/mouse38/ -chr:$CHR_MOUSE -compress  &> cron_mouseGene38.log

mailx -s "[$SERVER]Pipeline to create Gff3 data for genes ran" $EMAILLIST < logs/gene.log



$RUNLOAD -object:qtl -species:RAT -mapKey:60 -toFile:$LOGDIR/Qtl/Rat/rat34/ -compress  &> cron_ratQtl34.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Qtl assembly 3.4 ran" $EMAILLIST<cron_ratQtl34.log

$RUNLOAD -object:strain -species:RAT -mapKey:60 -toFile:$LOGDIR/Strain/Rat/rat34/ -compress  &> cron_ratStrain34.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Strain assembly 3.4 ran" $EMAILLIST<cron_ratStrain34.log

$RUNLOAD -object:sslp -species:RAT -mapKey:60 -toFile:$LOGDIR/Sslp/Rat/rat34/ -compress  &> cron_ratSslp34.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Sslp assembly 3.4 ran" $EMAILLIST<cron_ratSslp34.log


$RUNLOAD -object:qtl -species:RAT -mapKey:70 -toFile:$LOGDIR/Qtl/Rat/rat50/ -compress &> cron_ratQtl50.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Qtl assembly 5.0 ran" $EMAILLIST<cron_ratQtl50.log

$RUNLOAD -object:strain -species:RAT -mapKey:70 -toFile:$LOGDIR/Strain/Rat/rat50/ -compress &> cron_ratStrain50.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Strain assembly 5.0 ran" $EMAILLIST<cron_ratStrain50.log

$RUNLOAD -object:sslp -species:RAT -mapKey:70 -toFile:$LOGDIR/Sslp/Rat/rat50/ -compress &> cron_ratSslp50.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Sslp assembly 5.0 ran" $EMAILLIST<cron_ratSslp50.log


$RUNLOAD -object:qtl -species:RAT -mapKey:360 -toFile:$LOGDIR/Qtl/Rat/rat60/ -compress &> cron_ratQtl60.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Qtl assembly 6.0 ran" $EMAILLIST<cron_ratQtl60.log

$RUNLOAD -object:strain -species:RAT -mapKey:360 -toFile:$LOGDIR/Strain/Rat/rat60/ -compress &> cron_ratStrain60.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Strain assembly 6.0 ran" $EMAILLIST<cron_ratStrain60.log

$RUNLOAD -object:sslp -species:RAT -mapKey:360 -toFile:$LOGDIR/Sslp/Rat/rat60/ -compress &> cron_ratSslp60.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Sslp assembly 6.0 ran" $EMAILLIST<cron_ratSslp60.log


$RUNLOAD -object:qtl -species:HUMAN -mapKey:13 -toFile:$LOGDIR/Qtl/Human/human36/ -compress &> cron_humanQtl36.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Human Qtl Build 36 ran" $EMAILLIST<cron_humanQtl36.log

$RUNLOAD -object:qtl -species:HUMAN -mapKey:17 -toFile:$LOGDIR/Qtl/Human/human37/ -compress &> cron_humanQtl37.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Human Qtl Build 37 ran" $EMAILLIST<cron_humanQtl37.log

$RUNLOAD -object:qtl -species:HUMAN -mapKey:38 -toFile:$LOGDIR/Qtl/Human/human38/ -compress &> cron_humanQtl38.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Human Qtl Build 38 ran" $EMAILLIST<cron_humanQtl38.log


$RUNLOAD -object:sslp -species:HUMAN -mapKey:13 -toFile:$LOGDIR/Sslp/Human/human36/ -compress &> cron_humanSslp36.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Human Sslp assembly 36 ran" $EMAILLIST<cron_humanSslp36.log

$RUNLOAD -object:sslp -species:HUMAN -mapKey:17 -toFile:$LOGDIR/Sslp/Human/human37/ -compress &> cron_humanSslp37.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Human Sslp assembly 37 ran" $EMAILLIST<cron_humanSslp37.log


$RUNLOAD -object:sslp -species:MOUSE -mapKey:18 -toFile:$LOGDIR/Sslp/Mouse/mouse37/ -compress &> cron_mouseSslp37.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Mouse Sslp assembly 37 ran" $EMAILLIST<cron_mouseSslp37.log

$RUNLOAD -object:qtl -species:MOUSE -mapKey:18 -toFile:$LOGDIR/Qtl/Mouse/mouse37/ -compress &> cron_mouseQtl37.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Mouse Qtl assembly 37 ran" $EMAILLIST<cron_mouseQtl37.log


mailx -s "[$SERVER]Pipeline to create Gff3 data for Mouse Gene assembly 38 ran" $EMAILLIST < logs/gene.log

$RUNLOAD -object:sslp -species:MOUSE -mapKey:35 -toFile:$LOGDIR/Sslp/Mouse/mouse38/ -compress &> cron_mouseSslp38.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Mouse Sslp assembly 38 ran" $EMAILLIST<cron_mouseSslp38.log

$RUNLOAD -object:qtl -species:MOUSE -mapKey:35 -toFile:$LOGDIR/Qtl/Mouse/mouse38/ -compress &> cron_mouseQtl38.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Mouse Qtl assembly 38 ran" $EMAILLIST<cron_mouseQtl38.log


$RUNLOAD -ontAspect:D -species:RAT -mapKey:60 -toFile:$LOGDIR/Ont/Rat/rat34/ -chr:$CHR_RAT -compress &> cron_ratDisease34.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 3.4 Disease related elements ran" $EMAILLIST<cron_ratDisease34.log

$RUNLOAD -ontAspect:D -species:RAT -mapKey:70 -toFile:$LOGDIR/Ont/Rat/rat50/ -chr:$CHR_RAT -compress &> cron_ratDisease50.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 5 Disease related elements ran" $EMAILLIST<cron_ratDisease50.log

$RUNLOAD -ontAspect:D -species:RAT -mapKey:360 -toFile:$LOGDIR/Ont/Rat/rat60/ -chr:$CHR_RAT -compress &> cron_ratDisease60.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 6 Disease related elements ran" $EMAILLIST<cron_ratDisease60.log


$RUNLOAD -ontAspect:D -species:HUMAN -mapKey:13 -toFile:$LOGDIR/Ont/Human/human36/ -chr:$CHR_HUMAN -compress &> cron_humanDisease36.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human Disease assembly 36 related elements ran" $EMAILLIST<cron_humanDisease36.log

$RUNLOAD -ontAspect:D -species:HUMAN -mapKey:17 -toFile:$LOGDIR/Ont/Human/human37/ -chr:$CHR_HUMAN -compress &> cron_humanDisease37.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human Disease assembly 37 related elements ran" $EMAILLIST<cron_humanDisease37.log

$RUNLOAD -ontAspect:D -species:HUMAN -mapKey:38 -toFile:$LOGDIR/Ont/Human/human38/ -chr:$CHR_HUMAN -compress &> cron_humanDisease38.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human Disease assembly 38 related elements ran" $EMAILLIST<cron_humanDisease38.log


$RUNLOAD -ontAspect:D -species:MOUSE -mapKey:18 -toFile:$LOGDIR/Ont/Mouse/mouse37/ -chr:$CHR_MOUSE -compress &> cron_mouseDisease37.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for mouse Disease assembly 37 related elements ran" $EMAILLIST<cron_mouseDisease37.log

$RUNLOAD -ontAspect:D -species:MOUSE -mapKey:35 -toFile:$LOGDIR/Ont/Mouse/mouse38/ -chr:$CHR_MOUSE -compress &> cron_mouseDisease38.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for mouse Disease assembly 38 related elements ran" $EMAILLIST<cron_mouseDisease38.log



$RUNLOAD -ontAspect:E -species:RAT -mapKey:60 -toFile:$LOGDIR/Ont/Rat/rat34/ -chr:$CHR_RAT -compress &> cron_ratDrugGene34.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 3.4 DrugGene related elements ran" $EMAILLIST<cron_ratDrugGene34.log

$RUNLOAD -ontAspect:E -species:RAT -mapKey:70 -toFile:$LOGDIR/Ont/Rat/rat50/ -chr:$CHR_RAT -compress &> cron_ratDrugGene50.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 5 DrugGene related elements ran" $EMAILLIST<cron_ratDrugGene50.log

$RUNLOAD -ontAspect:E -species:RAT -mapKey:360 -toFile:$LOGDIR/Ont/Rat/rat60/ -chr:$CHR_RAT -compress &> cron_ratDrugGene60.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 6 DrugGene related elements ran" $EMAILLIST<cron_ratDrugGene60.log


$RUNLOAD -ontAspect:E -species:HUMAN -mapKey:13 -toFile:$LOGDIR/Ont/Human/human36/ -chr:$CHR_HUMAN -compress &> cron_humanDrugGene36.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human DrugGene assembly 36 related elements ran" $EMAILLIST<cron_humanDrugGene36.log

$RUNLOAD -ontAspect:E -species:HUMAN -mapKey:17 -toFile:$LOGDIR/Ont/Human/human37/ -chr:$CHR_HUMAN -compress &> cron_humanDrugGene37.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human DrugGene assembly 37 related elements ran" $EMAILLIST<cron_humanDrugGene37.log

$RUNLOAD -ontAspect:E -species:HUMAN -mapKey:38 -toFile:$LOGDIR/Ont/Human/human38/ -chr:$CHR_HUMAN -compress &> cron_humanDrugGene38.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human DrugGene assembly 38 related elements ran" $EMAILLIST<cron_humanDrugGene38.log


$RUNLOAD -ontAspect:E -species:MOUSE -mapKey:18 -toFile:$LOGDIR/Ont/Mouse/mouse37/ -chr:$CHR_MOUSE -compress &> cron_mouseDrugGene37.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for mouse DrugGene assembly 37 related elements ran" $EMAILLIST<cron_mouseDrugGene37.log

$RUNLOAD -ontAspect:E -species:MOUSE -mapKey:35 -toFile:$LOGDIR/Ont/Mouse/mouse38/ -chr:$CHR_MOUSE -compress &> cron_mouseDrugGene38.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for mouse DrugGene assembly 38 related elements ran" $EMAILLIST<cron_mouseDrugGene38.log


$RUNLOAD -object:promoter -species:RAT -mapKey:60 -toFile:$LOGDIR/Promoter/Rat/rat34/ -compress &> cron_ratPromoter34.log
mailx -s "[$SERVER]Pipeline to create Promoter Gff3 data for rat assembly 3.4 ran" $EMAILLIST<cron_ratPromoter34.log

$RUNLOAD -object:promoter -species:HUMAN -mapKey:13 -toFile:$LOGDIR/Promoter/Human/human36/ -compress &> cron_humanPromoter34.log
mailx -s "[$SERVER]Pipeline to create Promoter Gff3 data for human assembly 36 ran" $EMAILLIST<cron_humanPromoter34.log


echo "copy generated gff3 files to data_release directory:"
echo "  rsync -zarv --prune-empty-dirs --include='*/' --include='*.gff3.gz' --exclude='*' $LOGDIR/ $DATA_RELEASE_DIR_GFF3/"
rsync -zarv --prune-empty-dirs --include='*/' --include='*.gff3.gz' --exclude='*' $LOGDIR/ $DATA_RELEASE_DIR_GFF3/

echo "copy generated gff files to data_release directory:"
echo "  rsync -zarv --prune-empty-dirs --include='*/' --include='*.gff.gz' --exclude='*' $LOGDIR/ $DATA_RELEASE_DIR_GFF/"
rsync -zarv --prune-empty-dirs --include='*/' --include='*.gff.gz' --exclude='*' $LOGDIR/ $DATA_RELEASE_DIR_GFF/

echo "  rsync OK!"



#$APP_HOME/LoadIntoJBrowse.sh
#echo "OK: gff3 files loaded into JBrowse!" | mailx -s "[$SERVER] Gff3 files loaded into JBrowse!" $EMAILLIST
