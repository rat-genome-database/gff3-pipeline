#!/bin/bash
# load ClinVar tracks into human hg38, hg19 and hg18 assemblies

###### VARIABLE DECLARATIONS ######
jbrowseHome="/rgd/JBrowse-1.16.11/"

mainGFF3home="/rgd/data/gff3"
workingDir="pipelineOutput"

reedPulldown="scp -q rgdpub@reed.rgd.mcw.edu:"

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" = "REED" ]; then
  reedPulldown="cp -f "
fi

mainPipelineOutput="/home/rgddata/pipelines/RGDGff3Pipeline/data"

contentBlock="/rgdweb/jbrowse/contextMenu.html?source={source}&type={type}&seq_id={seq_id}&start={start}&end={end}&name={name}&id={id}&dbxRef={dbxref}&alias={alias}&symbol={symbol}&clinicalSignificance={clinicalSignificance}&methodType={methodType}&molecularConsequence={molecularConsequence}&ageOfOnset={ageOfOnset}&prevalence={prevalence}&submitter={submitter}&traitName={traitName}"
###################################

assemblies=( "38" "19" "18" )

for ASSEMBLY in "${assemblies[@]}"
do

	mainOrgDir="Human/hg${ASSEMBLY}"
	datasetDir="data_hg${ASSEMBLY}"
	organism=Human
	assembly="hg${ASSEMBLY}"


	printf "***********************************\n"
	printf "Beginning processing of ClinVar for assembly ${ASSEMBLY} ...\n"
	date
	printf "***********************************\n\n"
	
	## Initial path setting
	startingHome="${PWD}"
	currentOrgPath="${mainGFF3home}/${mainOrgDir}"
	currentWorkingPath="${currentOrgPath}/${workingDir}"
	printf "startingHome: ${startingHome}\n"
	printf "currentOrgPath: ${currentOrgPath}\n"
	printf "currentWorkingPath: ${currentWorkingPath}\n"
	
	## Initial directory setup
	printf "Initiating directory setup ... "
	if [[ ! -d ${currentWorkingPath} ]]; then
		mkdir -p ${currentWorkingPath}
	fi
	printf 'done!\n\n'
	
	## Initial clean-up
	printf "Initiating clean-up ... "
	rm -f  ${currentWorkingPath}/*.gff3* ${currentOrgPath}/*.gff3*
	printf 'done!\n\n'

	#pulldown from REED
	echo PULLDOWN from reed "${reedPulldown}${mainPipelineOutput}/ClinVar/ClinVar_hg${ASSEMBLY}.gff3.gz ${currentOrgPath}"
	echo
	${reedPulldown}${mainPipelineOutput}/ClinVar/ClinVar_hg${ASSEMBLY}.gff3.gz ${currentOrgPath}
	echo "PULLDOWN ok"


	## Unzip imported GFF3 files
	printf "Unzipping GFF3 files ... "
	for file in ${currentOrgPath}/*.gff3.gz ${currentWorkingPath}/*.gff3.gz 
	do
		gzip -d ${file}
	done
	printf 'done!\n\n'

	## JBrowse handling - where the magic happens!
	pushd ${jbrowseHome}
	printf "\nMoved to JBrowse home directory: ${PWD}\n\n"

	# Purging old track
	echo "remove ClinVar track"
	./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "ClinVar" --delete

	echo "insert clinvar track"
	./bin/flatfile-to-json.pl --gff "${currentOrgPath}/ClinVar_hg${ASSEMBLY}.gff3" --trackLabel "ClinVar" --key "RGD ${organism} (${assembly}) ClinVar" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"color\" : \"blue\", \"label\" : \"name,symbol\", \"description\" : \"\" }" --config "{ \"category\" : \"Variants/ClinVar\", \"maxHeight\" : 2400, \"hooks\" : { \"modify\" : \"function( track, f, fdiv ) { if( f.get('Clinicalsignificance').indexOf('uncertain')>=0 ) { fdiv.style.background = '#FF0000'; } }\" }, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"RGD Feature Data for {symbol}\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
	echo "+++ inserted track ClinVar"


	# Generation of names index
	echo "=== generation of names index ==="
	./bin/generate-names.pl --incremental --out "${datasetDir}" --tracks "ClinVar" --verbose --mem 4096000000

	popd
	printf "COMPLETE!\nMoved back to: ${PWD}\n\n"
	## End of JBrowse handling

done

