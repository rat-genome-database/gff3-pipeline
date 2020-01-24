# script env setup
#
. /etc/profile

APP_HOME=/home/rgddata/pipelines/RGDGff3Pipeline
cd $APP_HOME

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

$RUNLOAD -object:proteinDomain -species:PIG -mapKey:911 -toFile:$LOGDIR/ProteinDomain/pig11_domains.gff3 -compress  > pig11Domains.log &
$RUNLOAD -object:proteinDomain -species:PIG -mapKey:910 -toFile:$LOGDIR/ProteinDomain/pig10_domains.gff3 -compress  > pig10Domains.log &
$RUNLOAD -object:proteinDomain -species:CHINCHILLA -mapKey:44 -toFile:$LOGDIR/ProteinDomain/chinchilla10_domains.gff3 -compress  > chinchilla10Domains.log &
$RUNLOAD -object:proteinDomain -species:SQUIRREL -mapKey:720 -toFile:$LOGDIR/ProteinDomain/squirrel20_domains.gff3 -compress  > squirrel20Domains.log &
$RUNLOAD -object:proteinDomain -species:BONOBO -mapKey:511 -toFile:$LOGDIR/ProteinDomain/bonobo11_domains.gff3 -compress  > bonobo31Domains.log &
$RUNLOAD -object:proteinDomain -species:DOG -mapKey:631 -toFile:$LOGDIR/ProteinDomain/dog31_domains.gff3 -compress  > dog31Domains.log &
$RUNLOAD -object:proteinDomain -species:RAT -mapKey:360 -toFile:$LOGDIR/ProteinDomain/rat60_domains.gff3 -compress  > rat60Domains.log &
$RUNLOAD -object:proteinDomain -species:RAT -mapKey:70 -toFile:$LOGDIR/ProteinDomain/rat50_domains.gff3 -compress  > rat50Domains.log &
$RUNLOAD -object:proteinDomain -species:RAT -mapKey:60 -toFile:$LOGDIR/ProteinDomain/rat34_domains.gff3 -compress  > rat34Domains.log &
$RUNLOAD -object:proteinDomain -species:HUMAN -mapKey:38 -toFile:$LOGDIR/ProteinDomain/human38_domains.gff3 -compress  > human38Domains.log &
$RUNLOAD -object:proteinDomain -species:HUMAN -mapKey:17 -toFile:$LOGDIR/ProteinDomain/human37_domains.gff3 -compress  > human37Domains.log &
$RUNLOAD -object:proteinDomain -species:MOUSE -mapKey:35 -toFile:$LOGDIR/ProteinDomain/mouse38_domains.gff3 -compress  > mouse38Domains.log &
$RUNLOAD -object:proteinDomain -species:MOUSE -mapKey:18 -toFile:$LOGDIR/ProteinDomain/mouse37_domains.gff3 -compress  > mouse37Domains.log &

wait
mailx -s "[$SERVER]Pipeline to create Gff3 data for protein domains ran" $EMAILLIST < logs/domains.log


##### GENES

$RUNLOAD -object:gene -species:VERVET -mapKey:1311 -toDir:$LOGDIR/Gene/Vervet/ -chr:* -compress  > vervetGene11.log &
$RUNLOAD -object:gene -species:PIG -mapKey:911 -toDir:$LOGDIR/Gene/Pig/ -chr:* -compress  > pigGene11.log &
$RUNLOAD -object:gene -species:PIG -mapKey:910 -toDir:$LOGDIR/Gene/Pig/ -chr:* -compress  > pigGene10.log &
$RUNLOAD -object:gene -species:CHINCHILLA -mapKey:44 -toDir:$LOGDIR/Gene/Chinchilla/ -chr:* -compress  > chinchillaGene10.log &
$RUNLOAD -object:gene -species:SQUIRREL -mapKey:720 -toDir:$LOGDIR/Gene/Squirrel/ -chr:* -compress  > squirrelGene20.log &
$RUNLOAD -object:gene -species:BONOBO -mapKey:511 -toDir:$LOGDIR/Gene/Bonobo/ -chr:* -compress  > bonoboGene31.log &
$RUNLOAD -object:gene -species:DOG -mapKey:631 -toDir:$LOGDIR/Gene/Dog/ -chr:* -compress  > dogGene31.log &

$RUNLOAD -object:gene -species:MOUSE -mapKey:18 -toDir:$LOGDIR/Gene/Mouse/ -chr:* -compress  > mouseGene37.log &
$RUNLOAD -object:gene -species:MOUSE -mapKey:35 -toDir:$LOGDIR/Gene/Mouse/ -chr:* -compress  > mouseGene38.log &
wait

$RUNLOAD -object:gene -species:RAT -mapKey:60 -toDir:$LOGDIR/Gene/Rat/ -chr:* -compress  > ratGene34.log &
$RUNLOAD -object:gene -species:RAT -mapKey:70 -toDir:$LOGDIR/Gene/Rat/ -chr:* -compress  > ratGene50.log &
$RUNLOAD -object:gene -species:RAT -mapKey:360 -toDir:$LOGDIR/Gene/Rat/ -chr:* -compress  > ratGene60.log &

$RUNLOAD -object:gene -species:HUMAN -mapKey:13 -toDir:$LOGDIR/Gene/Human/ -chr:* -compress  > humanGene36.log &
$RUNLOAD -object:gene -species:HUMAN -mapKey:17 -toDir:$LOGDIR/Gene/Human/ -chr:* -compress  > humanGene37.log &
$RUNLOAD -object:gene -species:HUMAN -mapKey:38 -toDir:$LOGDIR/Gene/Human/ -chr:* -compress  > humanGene38.log &
wait

mailx -s "[$SERVER]Pipeline to create Gff3 data for genes ran" $EMAILLIST < logs/gene.log


##### QTLS

# assemblies, species info and output dirs are read from properties file AppConfigure.xml
$RUNLOAD -object:qtls > qtls.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for QTLS ran" $EMAILLIST < qtls.log


##### STRAINS, SSLPS

$RUNLOAD -object:strain -species:RAT -mapKey:60 -toFile:$LOGDIR/Strain/Rat/rat34/ -compress  > ratStrain34.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Strain assembly 3.4 ran" $EMAILLIST < ratStrain34.log

$RUNLOAD -object:sslp -species:RAT -mapKey:60 -toFile:$LOGDIR/Sslp/Rat/rat34/ -compress  > ratSslp34.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Sslp assembly 3.4 ran" $EMAILLIST < ratSslp34.log


$RUNLOAD -object:strain -species:RAT -mapKey:70 -toFile:$LOGDIR/Strain/Rat/rat50/ -compress > ratStrain50.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Strain assembly 5.0 ran" $EMAILLIST < ratStrain50.log

$RUNLOAD -object:sslp -species:RAT -mapKey:70 -toFile:$LOGDIR/Sslp/Rat/rat50/ -compress > ratSslp50.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Sslp assembly 5.0 ran" $EMAILLIST < ratSslp50.log


$RUNLOAD -object:strain -species:RAT -mapKey:360 -toFile:$LOGDIR/Strain/Rat/rat60/ -compress > ratStrain60.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Strain assembly 6.0 ran" $EMAILLIST < ratStrain60.log

$RUNLOAD -object:sslp -species:RAT -mapKey:360 -toFile:$LOGDIR/Sslp/Rat/rat60/ -compress > ratSslp60.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Sslp assembly 6.0 ran" $EMAILLIST < ratSslp60.log


$RUNLOAD -object:sslp -species:HUMAN -mapKey:13 -toFile:$LOGDIR/Sslp/Human/human36/ -compress > humanSslp36.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Human Sslp assembly 36 ran" $EMAILLIST < humanSslp36.log

$RUNLOAD -object:sslp -species:HUMAN -mapKey:17 -toFile:$LOGDIR/Sslp/Human/human37/ -compress > humanSslp37.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Human Sslp assembly 37 ran" $EMAILLIST < humanSslp37.log


$RUNLOAD -object:sslp -species:MOUSE -mapKey:18 -toFile:$LOGDIR/Sslp/Mouse/mouse37/ -compress > mouseSslp37.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Mouse Sslp assembly 37 ran" $EMAILLIST < mouseSslp37.log

$RUNLOAD -object:sslp -species:MOUSE -mapKey:35 -toFile:$LOGDIR/Sslp/Mouse/mouse38/ -compress > mouseSslp38.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Mouse Sslp assembly 38 ran" $EMAILLIST < mouseSslp38.log


##### DISEASE ONTOLOGY

$RUNLOAD -ontAspect:D -species:RAT -mapKey:60 -toFile:$LOGDIR/Ont/Rat/rat34/ -chr:* -compress &> cron_ratDisease34.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 3.4 Disease related elements ran" $EMAILLIST<cron_ratDisease34.log

$RUNLOAD -ontAspect:D -species:RAT -mapKey:70 -toFile:$LOGDIR/Ont/Rat/rat50/ -chr:* -compress &> cron_ratDisease50.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 5 Disease related elements ran" $EMAILLIST<cron_ratDisease50.log

$RUNLOAD -ontAspect:D -species:RAT -mapKey:360 -toFile:$LOGDIR/Ont/Rat/rat60/ -chr:* -compress &> cron_ratDisease60.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 6 Disease related elements ran" $EMAILLIST<cron_ratDisease60.log


$RUNLOAD -ontAspect:D -species:HUMAN -mapKey:13 -toFile:$LOGDIR/Ont/Human/human36/ -chr:* -compress &> cron_humanDisease36.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human Disease assembly 36 related elements ran" $EMAILLIST<cron_humanDisease36.log

$RUNLOAD -ontAspect:D -species:HUMAN -mapKey:17 -toFile:$LOGDIR/Ont/Human/human37/ -chr:* -compress &> cron_humanDisease37.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human Disease assembly 37 related elements ran" $EMAILLIST<cron_humanDisease37.log

$RUNLOAD -ontAspect:D -species:HUMAN -mapKey:38 -toFile:$LOGDIR/Ont/Human/human38/ -chr:* -compress &> cron_humanDisease38.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human Disease assembly 38 related elements ran" $EMAILLIST<cron_humanDisease38.log


$RUNLOAD -ontAspect:D -species:MOUSE -mapKey:18 -toFile:$LOGDIR/Ont/Mouse/mouse37/ -chr:* -compress &> cron_mouseDisease37.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for mouse Disease assembly 37 related elements ran" $EMAILLIST<cron_mouseDisease37.log

$RUNLOAD -ontAspect:D -species:MOUSE -mapKey:35 -toFile:$LOGDIR/Ont/Mouse/mouse38/ -chr:* -compress &> cron_mouseDisease38.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for mouse Disease assembly 38 related elements ran" $EMAILLIST<cron_mouseDisease38.log


$RUNLOAD -ontAspect:D -species:CHINCHILLA -mapKey:44 -toFile:$LOGDIR/Ont/Chinchilla/chinchilla10/ -chr:* -compress &> cron_chinchillaDisease10.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for chinchilla Disease assembly 1.0 related elements ran" $EMAILLIST<cron_chinchillaDisease10.log

$RUNLOAD -ontAspect:D -species:BONOBO -mapKey:511 -toFile:$LOGDIR/Ont/Bonobo/bonobo11/ -chr:* -compress &> cron_bonoboDisease11.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for bonobo Disease assembly 1.1 related elements ran" $EMAILLIST<cron_bonoboDisease11.log

$RUNLOAD -ontAspect:D -species:SQUIRREL -mapKey:720 -toFile:$LOGDIR/Ont/Squirrel/squirrel20/ -chr:* -compress &> cron_squirrelDisease20.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for squirrel Disease assembly 2.0 related elements ran" $EMAILLIST<cron_squirrelDisease20.log

$RUNLOAD -ontAspect:D -species:DOG -mapKey:631 -toFile:$LOGDIR/Ont/Dog/dog31/ -chr:* -compress &> cron_dogDisease31.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for dog Disease assembly 3.1 related elements ran" $EMAILLIST<cron_dogDisease31.log

$RUNLOAD -ontAspect:D -species:PIG -mapKey:910 -toFile:$LOGDIR/Ont/Pig/pig10/ -chr:* -compress &> cron_pigDisease10.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for pig Disease assembly 10 related elements ran" $EMAILLIST<cron_pigDisease10.log

$RUNLOAD -ontAspect:D -species:PIG -mapKey:911 -toFile:$LOGDIR/Ont/Pig/pig11/ -chr:* -compress &> cron_pigDisease11.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for pig Disease assembly 11 related elements ran" $EMAILLIST<cron_pigDisease11.log


##### CHEBI ONTOLOGY

$RUNLOAD -ontAspect:E -species:RAT -mapKey:60 -toFile:$LOGDIR/Ont/Rat/rat34/ -chr:* -compress &> cron_ratDrugGene34.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 3.4 DrugGene related elements ran" $EMAILLIST<cron_ratDrugGene34.log

$RUNLOAD -ontAspect:E -species:RAT -mapKey:70 -toFile:$LOGDIR/Ont/Rat/rat50/ -chr:* -compress &> cron_ratDrugGene50.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 5 DrugGene related elements ran" $EMAILLIST<cron_ratDrugGene50.log

$RUNLOAD -ontAspect:E -species:RAT -mapKey:360 -toFile:$LOGDIR/Ont/Rat/rat60/ -chr:* -compress &> cron_ratDrugGene60.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 6 DrugGene related elements ran" $EMAILLIST<cron_ratDrugGene60.log


$RUNLOAD -ontAspect:E -species:HUMAN -mapKey:13 -toFile:$LOGDIR/Ont/Human/human36/ -chr:* -compress &> cron_humanDrugGene36.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human DrugGene assembly 36 related elements ran" $EMAILLIST<cron_humanDrugGene36.log

$RUNLOAD -ontAspect:E -species:HUMAN -mapKey:17 -toFile:$LOGDIR/Ont/Human/human37/ -chr:* -compress &> cron_humanDrugGene37.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human DrugGene assembly 37 related elements ran" $EMAILLIST<cron_humanDrugGene37.log

$RUNLOAD -ontAspect:E -species:HUMAN -mapKey:38 -toFile:$LOGDIR/Ont/Human/human38/ -chr:* -compress &> cron_humanDrugGene38.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for human DrugGene assembly 38 related elements ran" $EMAILLIST<cron_humanDrugGene38.log


$RUNLOAD -ontAspect:E -species:MOUSE -mapKey:18 -toFile:$LOGDIR/Ont/Mouse/mouse37/ -chr:* -compress &> cron_mouseDrugGene37.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for mouse DrugGene assembly 37 related elements ran" $EMAILLIST<cron_mouseDrugGene37.log

$RUNLOAD -ontAspect:E -species:MOUSE -mapKey:35 -toFile:$LOGDIR/Ont/Mouse/mouse38/ -chr:* -compress &> cron_mouseDrugGene38.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for mouse DrugGene assembly 38 related elements ran" $EMAILLIST<cron_mouseDrugGene38.log


$RUNLOAD -ontAspect:E -species:CHINCHILLA -mapKey:44 -toFile:$LOGDIR/Ont/Chinchilla/chinchilla10/ -chr:* -compress &> cron_chinchillaDrugGene10.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for chinchilla DrugGene assembly 1.0 related elements ran" $EMAILLIST<cron_chinchillaDrugGene10.log

$RUNLOAD -ontAspect:E -species:BONOBO -mapKey:511 -toFile:$LOGDIR/Ont/Bonobo/bonobo11/ -chr:* -compress &> cron_bonoboDrugGene11.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for bonobo DrugGene assembly 1.1 related elements ran" $EMAILLIST<cron_bonoboDrugGene11.log

$RUNLOAD -ontAspect:E -species:SQUIRREL -mapKey:720 -toFile:$LOGDIR/Ont/Squirrel/squirrel20/ -chr:* -compress &> cron_squirrelDrugGene20.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for squirrel DrugGene assembly 2.0 related elements ran" $EMAILLIST<cron_squirrelDrugGene20.log

$RUNLOAD -ontAspect:E -species:DOG -mapKey:631 -toFile:$LOGDIR/Ont/Dog/dog31/ -chr:* -compress &> cron_dogDrugGene31.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for dog DrugGene assembly 3.1 related elements ran" $EMAILLIST<cron_dogDrugGene31.log

$RUNLOAD -ontAspect:E -species:PIG -mapKey:910 -toFile:$LOGDIR/Ont/Pig/pig10/ -chr:* -compress &> cron_pigDrugGene10.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for pig DrugGene assembly 10 related elements ran" $EMAILLIST<cron_pigDrugGene10.log

$RUNLOAD -ontAspect:E -species:PIG -mapKey:911 -toFile:$LOGDIR/Ont/Pig/pig11/ -chr:* -compress &> cron_pigDrugGene11.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for pig DrugGene assembly 11 related elements ran" $EMAILLIST<cron_pigDrugGene11.log


##### PROMOTERS
#     all species: rat, mouse, human, dog

$RUNLOAD -object:promoter -toDir:$LOGDIR/Promoter/ -compress &> promoters.log
mailx -s "[$SERVER]Pipeline to create promoter gff3 files ran" $EMAILLIST<promoters.log



##### FINAL WRAPUP

echo "copy generated gff3 files to data_release directory:"
echo "  rsync -zarv --prune-empty-dirs --include='*/' --include='*.gff3.gz' --exclude='*' $LOGDIR/ $DATA_RELEASE_DIR_GFF3/"
rsync -zarv --prune-empty-dirs --include='*/' --include='*.gff3.gz' --exclude='*' $LOGDIR/ $DATA_RELEASE_DIR_GFF3/

echo "copy generated gff files to data_release directory:"
echo "  rsync -zarv --prune-empty-dirs --include='*/' --include='*.gff.gz' --exclude='*' $LOGDIR/ $DATA_RELEASE_DIR_GFF/"
rsync -zarv --prune-empty-dirs --include='*/' --include='*.gff.gz' --exclude='*' $LOGDIR/ $DATA_RELEASE_DIR_GFF/

echo "  rsync OK!"



#$APP_HOME/LoadIntoJBrowse.sh
#echo "OK: gff3 files loaded into JBrowse!" | mailx -s "[$SERVER] Gff3 files loaded into JBrowse!" $EMAILLIST
