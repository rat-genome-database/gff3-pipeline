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

$RUNLOAD -object:proteinDomains > proteinDomains.log

mailx -s "[$SERVER]Pipeline to create Gff3 data for protein domains ran" $EMAILLIST < logs/domains.log


##### GENES

# assemblies, species info and output dirs are read from properties file AppConfigure.xml
$RUNLOAD -object:genes > gene.log

mailx -s "[$SERVER]Pipeline to create Gff3 data for GENES ran" $EMAILLIST < logs/gene.log


##### QTLS
# assemblies, species info and output dirs are read from properties file AppConfigure.xml
$RUNLOAD -object:qtls > qtls.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for QTLS ran" $EMAILLIST < qtls.log


##### MARKERS
# assemblies, species info and output dirs are read from properties file AppConfigure.xml
$RUNLOAD -object:markers > markers.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for MARKERS ran" $EMAILLIST < markers.log





##### STRAINS

$RUNLOAD -object:strain -species:RAT -mapKey:60 -toFile:$LOGDIR/Strain/Rat/rat34/ -compress  > ratStrain34.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Strain assembly 3.4 ran" $EMAILLIST < ratStrain34.log

$RUNLOAD -object:strain -species:RAT -mapKey:70 -toFile:$LOGDIR/Strain/Rat/rat50/ -compress > ratStrain50.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Strain assembly 5.0 ran" $EMAILLIST < ratStrain50.log

$RUNLOAD -object:strain -species:RAT -mapKey:360 -toFile:$LOGDIR/Strain/Rat/rat60/ -compress > ratStrain60.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Strain assembly 6.0 ran" $EMAILLIST < ratStrain60.log

$RUNLOAD -object:strain -species:RAT -mapKey:372 -toFile:$LOGDIR/Strain/Rat/rat72/ -compress > ratStrain72.log
mailx -s "[$SERVER]Pipeline to create Gff3 data for Rat Strain assembly 7.2 ran" $EMAILLIST < ratStrain72.log


##### DISEASE ONTOLOGY

$RUNLOAD -ontAspect:D -species:RAT -mapKey:60 -toFile:$LOGDIR/Ont/Rat/rat34/ -chr:* -compress &> cron_ratDisease34.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 3.4 Disease related elements ran" $EMAILLIST<cron_ratDisease34.log

$RUNLOAD -ontAspect:D -species:RAT -mapKey:70 -toFile:$LOGDIR/Ont/Rat/rat50/ -chr:* -compress &> cron_ratDisease50.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 5 Disease related elements ran" $EMAILLIST<cron_ratDisease50.log

$RUNLOAD -ontAspect:D -species:RAT -mapKey:360 -toFile:$LOGDIR/Ont/Rat/rat60/ -chr:* -compress &> cron_ratDisease60.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 6 Disease related elements ran" $EMAILLIST<cron_ratDisease60.log

$RUNLOAD -ontAspect:D -species:RAT -mapKey:372 -toFile:$LOGDIR/Ont/Rat/rat72/ -chr:* -compress &> cron_ratDisease72.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 7 Disease related elements ran" $EMAILLIST<cron_ratDisease72.log


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

$RUNLOAD -ontAspect:D -species:MOUSE -mapKey:239 -toFile:$LOGDIR/Ont/Mouse/mouse39/ -chr:* -compress &> cron_mouseDisease39.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for mouse Disease assembly 39 related elements ran" $EMAILLIST<cron_mouseDisease39.log


$RUNLOAD -ontAspect:D -species:CHINCHILLA -mapKey:44 -toFile:$LOGDIR/Ont/Chinchilla/chinchilla10/ -chr:* -compress &> cron_chinchillaDisease10.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for chinchilla Disease assembly 1.0 ran" $EMAILLIST<cron_chinchillaDisease10.log

$RUNLOAD -ontAspect:D -species:BONOBO -mapKey:511 -toFile:$LOGDIR/Ont/Bonobo/bonobo11/ -chr:* -compress &> cron_bonoboDisease11.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for bonobo Disease assembly 1.1 ran" $EMAILLIST<cron_bonoboDisease11.log

$RUNLOAD -ontAspect:D -species:BONOBO -mapKey:513 -toFile:$LOGDIR/Ont/Bonobo/bonobo2/ -chr:* -compress &> cron_bonoboDisease2.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for bonobo Disease assembly Mhudiblu_PPA_v0 ran" $EMAILLIST<cron_bonoboDisease2.log

$RUNLOAD -ontAspect:D -species:SQUIRREL -mapKey:720 -toFile:$LOGDIR/Ont/Squirrel/squirrel20/ -chr:* -compress &> cron_squirrelDisease20.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for squirrel Disease assembly 2.0 ran" $EMAILLIST<cron_squirrelDisease20.log

$RUNLOAD -ontAspect:D -species:DOG -mapKey:631 -toFile:$LOGDIR/Ont/Dog/dog31/ -chr:* -compress &> cron_dogDisease31.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for dog Disease assembly 3.1 ran" $EMAILLIST<cron_dogDisease31.log

$RUNLOAD -ontAspect:D -species:PIG -mapKey:910 -toFile:$LOGDIR/Ont/Pig/pig10/ -chr:* -compress &> cron_pigDisease10.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for pig Disease assembly 10 ran" $EMAILLIST<cron_pigDisease10.log

$RUNLOAD -ontAspect:D -species:PIG -mapKey:911 -toFile:$LOGDIR/Ont/Pig/pig11/ -chr:* -compress &> cron_pigDisease11.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for pig Disease assembly 11 ran" $EMAILLIST<cron_pigDisease11.log

$RUNLOAD -ontAspect:D -species:SQUIRREL -mapKey:720 -toFile:$LOGDIR/Ont/Squirrel/squirrel20/ -chr:* -compress &> cron_squirrelDisease20.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for squirrel Disease assembly 2.0 ran" $EMAILLIST<cron_squirrelDisease20.log

$RUNLOAD -ontAspect:D -species:VERVET -mapKey:1311 -toFile:$LOGDIR/Ont/Vervet/chlSab2/ -chr:* -compress &> cron_vervetDisease2.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for vervet ChlSab1.1 assembly ran" $EMAILLIST<cron_vervetDisease2.log

$RUNLOAD -ontAspect:D -species:VERVET -mapKey:1313 -toFile:$LOGDIR/Ont/Vervet/Vero_WHO_p1/ -chr:* -compress &> cron_vervetDiseaseVero.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for vervet Vero_WHO_p1.0 assembly ran" $EMAILLIST<cron_vervetDiseaseVero.log

$RUNLOAD -ontAspect:D -species:MOLERAT -mapKey:1410 -toFile:$LOGDIR/Ont/Molerat/hetGla2/ -chr:* -compress &> cron_moleratDisease2.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for molerat hetGla2 assembly ran" $EMAILLIST<cron_moleratDisease2.log

##### CHEBI ONTOLOGY

$RUNLOAD -ontAspect:E -species:RAT -mapKey:60 -toFile:$LOGDIR/Ont/Rat/rat34/ -chr:* -compress &> cron_ratDrugGene34.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 3.4 DrugGene related elements ran" $EMAILLIST<cron_ratDrugGene34.log

$RUNLOAD -ontAspect:E -species:RAT -mapKey:70 -toFile:$LOGDIR/Ont/Rat/rat50/ -chr:* -compress &> cron_ratDrugGene50.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 5 DrugGene related elements ran" $EMAILLIST<cron_ratDrugGene50.log

$RUNLOAD -ontAspect:E -species:RAT -mapKey:360 -toFile:$LOGDIR/Ont/Rat/rat60/ -chr:* -compress &> cron_ratDrugGene60.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 6 DrugGene related elements ran" $EMAILLIST<cron_ratDrugGene60.log

$RUNLOAD -ontAspect:E -species:RAT -mapKey:372 -toFile:$LOGDIR/Ont/Rat/rat72/ -chr:* -compress &> cron_ratDrugGene72.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for rat 7.2 DrugGene related elements ran" $EMAILLIST<cron_ratDrugGene72.log


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

$RUNLOAD -ontAspect:E -species:MOUSE -mapKey:239 -toFile:$LOGDIR/Ont/Mouse/mouse39/ -chr:* -compress &> cron_mouseDrugGene39.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for mouse DrugGene assembly 39 related elements ran" $EMAILLIST<cron_mouseDrugGene39.log


$RUNLOAD -ontAspect:E -species:CHINCHILLA -mapKey:44 -toFile:$LOGDIR/Ont/Chinchilla/chinchilla10/ -chr:* -compress &> cron_chinchillaDrugGene10.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for chinchilla DrugGene assembly 1.0 ran" $EMAILLIST<cron_chinchillaDrugGene10.log

$RUNLOAD -ontAspect:E -species:BONOBO -mapKey:511 -toFile:$LOGDIR/Ont/Bonobo/bonobo11/ -chr:* -compress &> cron_bonoboDrugGene11.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for bonobo DrugGene assembly 1.1 ran" $EMAILLIST<cron_bonoboDrugGene11.log

$RUNLOAD -ontAspect:E -species:BONOBO -mapKey:513 -toFile:$LOGDIR/Ont/Bonobo/bonobo2/ -chr:* -compress &> cron_bonoboDrugGene2.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for bonobo DrugGene assembly Mhudiblu_PPA_v0 ran" $EMAILLIST<cron_bonoboDrugGene2.log

$RUNLOAD -ontAspect:E -species:SQUIRREL -mapKey:720 -toFile:$LOGDIR/Ont/Squirrel/squirrel20/ -chr:* -compress &> cron_squirrelDrugGene20.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for squirrel DrugGene assembly 2.0 related elements ran" $EMAILLIST<cron_squirrelDrugGene20.log

$RUNLOAD -ontAspect:E -species:DOG -mapKey:631 -toFile:$LOGDIR/Ont/Dog/dog31/ -chr:* -compress &> cron_dogDrugGene31.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for dog DrugGene assembly 3.1 related elements ran" $EMAILLIST<cron_dogDrugGene31.log

$RUNLOAD -ontAspect:E -species:PIG -mapKey:910 -toFile:$LOGDIR/Ont/Pig/pig10/ -chr:* -compress &> cron_pigDrugGene10.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for pig DrugGene assembly 10 related elements ran" $EMAILLIST<cron_pigDrugGene10.log

$RUNLOAD -ontAspect:E -species:PIG -mapKey:911 -toFile:$LOGDIR/Ont/Pig/pig11/ -chr:* -compress &> cron_pigDrugGene11.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for pig DrugGene assembly 11 related elements ran" $EMAILLIST<cron_pigDrugGene11.log

$RUNLOAD -ontAspect:E -species:VERVET -mapKey:1311 -toFile:$LOGDIR/Ont/Vervet/chlSab2/ -chr:* -compress &> cron_vervetDrugGene2.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for vervet ChlSab1.1 assembly ran" $EMAILLIST<cron_vervetDrugGene2.log

$RUNLOAD -ontAspect:E -species:VERVET -mapKey:1313 -toFile:$LOGDIR/Ont/Vervet/Vero_WHO_p1/ -chr:* -compress &> cron_vervetDrugGeneVero.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for vervet Vero_WHO_p1.0 assembly ran" $EMAILLIST<cron_vervetDrugGeneVero.log

$RUNLOAD -ontAspect:E -species:MOLERAT -mapKey:1410 -toFile:$LOGDIR/Ont/Molerat/hetGla2/ -chr:* -compress &> cron_moleratDrugGene2.log
mailx -s "[$SERVER]Pipeline to create Ontology Gff3 data for molerat hetGla2 assembly ran" $EMAILLIST<cron_moleratDrugGene2.log

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
