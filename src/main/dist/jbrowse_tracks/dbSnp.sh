#!/bin/bash
# load dbSnp tracks into human hg38 and hg19

###### VARIABLE DECLARATIONS ######
jbrowseHome="/rgd/JBrowse-1.12.3/"

mainGFF3home="/rgd/data/gff3"
workingDir="pipelineOutput"

reedPulldown="scp -q rgdpub@reed.rgd.mcw.edu:"

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" = "REED" ]; then
  reedPulldown="cp -f "
fi

mainPipelineOutput="/home/rgddata/pipelines/RGDGff3Pipeline/dist/log/RGDGFF3/Output"

contentBlock="/rgdweb/jbrowse/contextMenu.html?source={source}&type={type}&seq_id={seq_id}&start={start}&end={end}&name={name}&id={id}&dbxRef={dbxref}&alias={alias}&symbol={symbol}&clinicalSignificance={clinicalSignificance}&allele={allele}&functionClass={functionClass}"
###################################

assemblies=( "38" "19" )

for ASSEMBLY in "${assemblies[@]}"
do

	mainOrgDir="Human/hg${ASSEMBLY}"
	datasetDir="data_hg${ASSEMBLY}"
	organism=Human
	assembly="hg${ASSEMBLY}"


	printf "***********************************\n"
	printf "Beginning processing of dbSnp for assembly ${ASSEMBLY} ...\n"
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
	echo PULLDOWN from reed "${reedPulldown}${mainPipelineOutput}/DbSnp/DbSnp_hg${ASSEMBLY}.gff3.gz ${currentOrgPath}"
	echo
	${reedPulldown}${mainPipelineOutput}/DbSnp/DbSnp_hg${ASSEMBLY}.gff3.gz ${currentOrgPath}
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
	echo "remove dbSNP track"
	./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "dbSNP" --delete

	echo "insert dbSNP track"
	./bin/flatfile-to-json.pl --gff "${currentOrgPath}/DbSnp_hg${ASSEMBLY}.gff3" --trackLabel "dbSNP" --key "RGD ${organism} (${assembly}) dbSNP build 150" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"color\" : \"green\", \"label\" : \"name,symbol\", \"description\" : \"\" }" --config "{ \"category\" : \"Variants/dbSNP\", \"maxHeight\" : 2400, \"hooks\" : { \"modify\" : \"function( track, f, fdiv ) { }\" }, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"RGD Feature Data for {symbol}\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
	echo "+++ inserted track dbSNP"


	# Generation of names index
	echo "=== generation of names index ==="
	./bin/generate-names.pl --incremental --out "${datasetDir}" --tracks "dbSNP" --verbose --mem 4096000000

	popd
	printf "COMPLETE!\nMoved back to: ${PWD}\n\n"
	## End of JBrowse handling

done

