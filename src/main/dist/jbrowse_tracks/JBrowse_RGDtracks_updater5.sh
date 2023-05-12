#!/bin/bash

############### DOCUMENTATION ###############
# This script is the weekly tracks updater pipeline for the RGD JBrowse installation.
# 
# It *MUST* be stored and executed from:
# /rgd/data/gff3
# 
# This pipeline uses two detailed input files for its execution, the first of which
# relates to general assembly information (e.g. Rat5, Rat3.4, and so on), while the
# second pertains to specific ontology track-related information.
# 
# An example run of this pipeline:
# 
# cd /rgd/data/gff3
# bash JBrowse_RGDtracks_updater5.sh JBrowse_RGDtracks_updaterInput.conf JBrowse_RGDtracks_updaterOntInput.conf | tee JBrowse_RGDtracks_updater.log
# 
# (The tee command allows for simultaneous re-direction of standard output to both
# the terminal and a log file for review at a later time.)
# 
# Normally, this pipeline script is run automatically as part of the weekly data release.
# 
# The wrapper script used for this weekly pipeline can be found in:
# /rgd/scripts/LoadGff3IntoJBrowse_Reed.sh
# 
# Aurash
#############################################

jbrowseHome="/rgd/JBrowse-1.16.11/"

mainGFF3home="/rgd/data/gff3"
workingDir="pipelineOutput"
workingOntDir="ontology"

reedPulldown="scp -q rgdpub@reed.rgd.mcw.edu:"

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" = "REED" ]; then
  reedPulldown="cp -f "
fi

mainPipelineOutput="/home/rgddata/pipelines/RGDGff3Pipeline/data"

contentBlock="https://rgd.mcw.edu/rgdweb/jbrowse/contextMenu.html?source={source}&type={type}&seq_id={seq_id}&start={start}&end={end}&name={name}&fullName={fullname}&geneType={genetype}&refSeqStatus={refseqstatus}&species={species}&note={note}&id={id}&mappingMethod={mappingmethod}&LOD={lod}&pValue={pvalue}&dbxRef={dbxref}&relatedQTLs={relatedqtls}&relatedStrains={relatedStrains}&relatedGenes={relatedgenes}&alias={alias}&diseaseOntologyAssociation={diseaseontologyassociation}&phenotypeOntologyAssociation={phenotypeontologyassociation}&synBlockLength={synblocklength}&synStart={synstart}&synEnd={synend}&synChrNum={synchrnum}&gene={gene}&isNonCoding={isnoncoding}&mapWeight={mapweight}&allele={allele}&strain={strain}&expectedSize={expectedsize}&symbol={symbol}&ontologyTerms={ontology_term}&product={product}&transcriptID={transcript_id}&proteinID={protein_id}&codonStart={codon_start}&promoterGeneRlt={promotergenerlt}&experimentMethods={experimentmethods}&regulation={regulation}&objectType={objecttype}&tissues={tissues}&score={score}&method={method}&reference={reference}&loc1={loc1}&loc2={loc2}&loc3={loc3}&loc4={loc4}&loc5={loc5}&loc6={loc6}&origin={origin}"
###################################

# Main while-loop execution block
while IFS=$'\t' read mainOrgDir pipelineGenes pipelineQTL pipelineSSLP pipelineStrains pipelineOnt datasetDir organism assembly pipelineGenesOnly
do
	printf "***********************************\n"
	printf "Beginning processing of ${assembly} for species ${organism} ...\n"
	date
	printf "***********************************\n\n"
	
	## Initial path setting
	startingHome="${PWD}"
	currentOrgPath="${mainGFF3home}/${mainOrgDir}"
	currentWorkingPath="${currentOrgPath}/${workingDir}"
	currentWorkingOntPath="${currentWorkingPath}/${workingOntDir}"
	printf "startingHome: ${startingHome}\n"
	printf "currentOrgPath: ${currentOrgPath}\n"
	printf "currentWorkingPath: ${currentWorkingPath}\n"
	printf "currentWorkingOntPath: ${currentWorkingOntPath}\n\n"
	
	## Initial directory setup
	printf "Initiating directory setup ... "
	if [[ ! -d ${currentWorkingPath} ]]; then
		mkdir -p ${currentWorkingPath}
	fi
	if [[ ! -d ${currentWorkingOntPath} ]]; then
		mkdir -p ${currentWorkingOntPath}
	fi
	printf 'done!\n\n'
	
	## Import pipeline GFF3 files
	printf "Importing pipeline GFF3 files ... "
	  ${reedPulldown}${mainPipelineOutput}/${pipelineGenes} ${currentOrgPath}/ARGD_curated_genes.gff3.gz
	  ${reedPulldown}${mainPipelineOutput}/${pipelineGenesOnly} ${currentOrgPath}/ARGD_curated_genesOnly.gff3.gz
	printf "  GENES OK"
	if [[ "${pipelineQTL}" != "null" ]]; then
		${reedPulldown}${mainPipelineOutput}/${pipelineQTL} ${currentOrgPath}/AQTLS.gff3.gz
		printf "  QTLS OK"
	fi
	if [[ "${pipelineSSLP}" != "null" ]]; then
		${reedPulldown}${mainPipelineOutput}/${pipelineSSLP} ${currentOrgPath}/markers.gff3.gz
		printf "  SSLPS OK"
	fi
	if [[ "${pipelineStrains}" != "null" ]]; then
		${reedPulldown}${mainPipelineOutput}/${pipelineStrains} ${currentOrgPath}
		printf "  STRAINS OK"
	fi
	
	# Pull ontology-related tracks
	if [[ "${pipelineOnt}" != "null" ]]; then
		${reedPulldown}${mainPipelineOutput}/${pipelineOnt} ${currentWorkingOntPath}
		printf "  ONTOLOGIES OK"
	fi
	printf ' done!\n\n'
	## End of pipeline importing
	
	
	# And here come the ontology tracks!
	if [[ "${pipelineOnt}" != "null" ]]; then
	while IFS=$'\t' read ontAccNum ontLabel ontName ontCategory ; do
		if [[ ${ontAccNum} =~ ^DOID.* ]] ; then
			# DO+ Genes track
			# Every organism gets an ontology genes track
			printf "Building ${ontLabel}Genes ... "
			cat $(ls -v ${currentWorkingOntPath}/${ontAccNum}*.gff3) | sed '/^##.*/d' | awk '{ FS=OFS="\t" } { if( $3 == "gene_ont" ) { print; } }' > ${currentWorkingOntPath}/${ontLabel}Genes.gff3
			printf 'done!\n'
			
			# DO+ QTLs track
			# Rat and Human get ontology QTL tracks
			if [[ ${organism} =~ ^[rR]at.* || ${organism} =~ ^[hH]uman.* ]] ; then
				printf "Building ${ontLabel}QTLs ... "
				cat $(ls -v ${currentWorkingOntPath}/${ontAccNum}*.gff3) | sed '/^##.*/d' | awk '{ FS=OFS="\t" } { if( $3 == "qtl_ont" ) { print; } }' > ${currentWorkingOntPath}/${ontLabel}QTLs.gff3
				printf 'done!\n'
			fi
			
			# DO+ Strains track
			# Only Rat gets ontology strain tracks
			if [[ ${organism} =~ ^[rR]at.* ]] ; then
				printf "Building ${ontLabel}Strains ... "
				cat $(ls -v ${currentWorkingOntPath}/${ontAccNum}*.gff3) | sed '/^##.*/d' | awk '{ FS=OFS="\t" } { if( $3 == "strain_ont" ) { print; } }' > ${currentWorkingOntPath}/${ontLabel}Strains.gff3
				printf 'done!\n'
			fi
			
			# Clean-up
			rm -f ${currentWorkingOntPath}/${ontAccNum}*.gff3
		elif [[ ${ontAccNum} =~ ^CHEBI.* ]] ; then
			# CHEBI Genes track
			printf "Building ${ontLabel} ... "
			cat $(ls -v ${currentWorkingOntPath}/${ontAccNum}*.gff3) | sed '/^##.*/d' > ${currentWorkingOntPath}/${ontLabel}.gff3
			printf 'done!\n'
			
			# Clean-up
			rm -f ${currentWorkingOntPath}/${ontAccNum}*.gff3
		fi
	done < ${2}
	fi

	## JBrowse handling - where the magic happens!
	pushd ${jbrowseHome}
	printf "\nMoved to JBrowse home directory: ${PWD}\n\n"
	
	# Purging old Gene Models, QTLs, and SSLP tracks
	if [[ "${pipelineQTL}" != null ]]; then
        ./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "AQTLS" --delete
	fi
	if [[ "${pipelineSSLP}" != null ]]; then
        	./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "SSLP" --delete
	fi

	# Removal of congenic and mutant strains tracks for Rat assemblies
	if [[ ${organism} =~ ^[rR]at.* ]]; then
		./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "CongenicStrains" --delete
		./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "MutantStrains" --delete
	fi
	
	# (One-time) Removal of out-of-date ontology tracks
	if [[ -e ${currentOrgPath}/testParser.log ]]
	then
		while IFS=$'\t' read trackLabel trackName trackCategory
		do
			if [[ ${trackCategory} =~ ^Disease\ Related.* ]]
			then
				./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "${trackLabel}" --delete
			elif [[ ${trackCategory} =~ ^Gene-Chemical.* ]]
			then
				./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "${trackLabel}" --delete
			fi
		done < ${currentOrgPath}/testParser.log
	fi
	
	# Removal of up-to-date ontology tracks:
	# Attempt to remove every possible ontology track,
	# even if it does not exist within the current JBrowse assembly repository
	if [[ "${pipelineOnt}" != "null" ]]; then
	while IFS=$'\t' read ontAccNum ontLabel ontName ontCategory
	do
		if [[ ${ontAccNum} =~ ^DOID.* ]]
		then
			# Remove 3 DO+ tracks per accession number
			./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "${ontLabel}Genes" --delete
			./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "${ontLabel}QTLs" --delete
			./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "${ontLabel}Strains" --delete
		elif [[ ${ontAccNum} =~ ^CHEBI.* ]]
		then
			# Remove 1 CHEBI track per accession number
			./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "${ontLabel}" --delete
		fi
	done < ${startingHome}/${2}
	fi
	
	###
	# GENES

	# obsolete track ARGD_curated_transcripts
  ./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "ARGD_curated_transcripts" --delete

	# ARGD_curated_genes insertion
  zcat "${currentOrgPath}/ARGD_curated_genes.gff3.gz" > "${currentOrgPath}/ARGD_curated_genes.gff3"
	./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "ARGD_curated_genes" --delete
	./bin/flatfile-to-json.pl --gff "${currentOrgPath}/ARGD_curated_genes.gff3" --trackLabel "ARGD_curated_genes" --key "RGD ${organism} (${assembly}) Genes and Transcripts" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"color\" : \"#DF255A\", \"featureScale\" : 0.0001, \"label\" : \"symbol,name\", \"description\" : \"\" }" --config "{ \"category\" : \"Gene Models/RGD Gene Features\", \"histograms\" : { \"color\" : \"#DF255A\", \"binsPerBlock\" : 25 }, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
  rm "${currentOrgPath}/ARGD_curated_genes.gff3"
	echo "+++ inserted track ARGD_curated_genes"
	
	# ARGD_curated_genesOnly insertion
  zcat "${currentOrgPath}/ARGD_curated_genesOnly.gff3.gz" > "${currentOrgPath}/ARGD_curated_genesOnly.gff3"
  ./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "ARGD_curated_genesOnly" --delete
	./bin/flatfile-to-json.pl --gff "${currentOrgPath}/ARGD_curated_genesOnly.gff3" --trackLabel "ARGD_curated_genesOnly" --key "RGD ${organism} (${assembly}) Genes" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"color\" : \"#DF255A\", \"featureScale\" : 0.0001, \"label\" : \"symbol,name\", \"description\" : \"\" }" --config "{ \"category\" : \"Gene Models/RGD Gene Features\", \"glyph\" : \"JBrowse/View/FeatureGlyph/Box\", \"histograms\" : { \"color\" : \"#DF255A\", \"binsPerBlock\" : 25 }, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
	rm "${currentOrgPath}/ARGD_curated_genesOnly.gff3"
	echo "+++ inserted track ARGD_curated_genesOnly"


  ###
	# AQTLS insertion
	if [[ "${pipelineQTL}" != "null" ]]; then
	./bin/flatfile-to-json.pl --gff "${currentOrgPath}/AQTLS.gff3" --trackLabel "AQTLS" --key "RGD ${organism} (${assembly}) QTLs" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"color\" : \"#255ADF\", \"featureScale\" : 0.00002, \"label\" : \"symbol,name\", \"description\" : \"\" }" --config "{ \"category\" : \"QTLs/${organism} QTLs\", \"maxHeight\" : 2400, \"glyph\" : \"JBrowse/View/FeatureGlyph/Alignment\", \"histograms\" : { \"color\" : \"#255ADF\", \"binsPerBlock\" : 25 }, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
	echo "+++ inserted track AQTLS"
	fi
	
	# SSLP insertion
	if [[ "${pipelineSSLP}" != "null" ]]; then
	./bin/flatfile-to-json.pl --gff "${currentOrgPath}/markers.gff3" --trackLabel "SSLP" --key "RGD ${organism} (${assembly}) Micro Satellite Markers" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"color\" : \"#378D15\", \"featureScale\" : 0.0001, \"label\" : \"symbol,name\", \"description\" : \"\" }" --config "{ \"category\" : \"Variants/Micro Satellite Markers\", \"histograms\" : { \"color\" : \"#378D15\", \"binsPerBlock\" : 25 }, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
		echo "+++ inserted track SSLP"
	fi

	# Congenic and Mutant Strains track insertion for Rat assemblies
	if [[ ${organism} =~ ^[rR]at.* ]]
	then
		# Congenic Strains insertion
		./bin/flatfile-to-json.pl --gff "${currentOrgPath}/RatCongenicStrains_RGD.gff3" --trackLabel "CongenicStrains" --key "RGD ${organism} (${assembly}) Congenic Strains Assembly" --out "${datasetDir}" --trackType JBrowse/View/Track/HTMLFeatures --clientConfig "{ \"className\" : \"feature3\", \"featureCss\" : \"background-color: #5ADF25;\", \"label\" : \"symbol,name\", \"description\" : \"\" }" --config "{ \"category\" : \"Strains/Congenic Strains\", \"maxHeight\" : 2400, \"hooks\" : { \"modify\" : \"function( track, f, fdiv ) { if( f.get( \\\"type\\\" ).match( /^multicongenic/i ) != null ) { fdiv.style.background = \\\"#AA25DF\\\"; } }\" }, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
		echo "+++ inserted track CongenicStrains"

		# Mutant Strains insertion
		./bin/flatfile-to-json.pl --gff "${currentOrgPath}/RatMutantStrains_RGD.gff3" --trackLabel "MutantStrains" --key "RGD ${organism} (${assembly}) Mutant Strains Assembly" --out "${datasetDir}" --trackType JBrowse/View/Track/HTMLFeatures --clientConfig "{ \"className\" : \"feature3\", \"featureCss\" : \"background-color: #5ADF25;\", \"label\" : \"symbol,name\", \"description\" : \"\" }" --config "{ \"category\" : \"Strains/Mutant Strains\", \"maxHeight\" : 2400, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
		echo "+++ inserted track MutantStrains"
	fi
	
	# Ontology track insertion: here goes the big fish!
	if [[ "${pipelineOnt}" != "null" ]]; then
	while IFS=$'\t' read ontAccNum ontLabel ontName ontCategory
	do
		if [[ ${ontAccNum} =~ ^DOID.* ]]; then
			# Insertion of DO+ tracks
			# 
			# DO+ Genes track insertion:
			# this occurs for every assembly
			./bin/flatfile-to-json.pl --gff "${currentWorkingOntPath}/${ontLabel}Genes.gff3" --trackLabel "${ontLabel}Genes" --key "${ontName} Related Genes" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"color\" : \"#DF255A\", \"label\" : \"symbol,name\", \"description\" : \"\" }" --config "{ \"category\" : \"${ontCategory}\", \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
			echo "  +++ inserted track ${ontLabel}Genes"

			# DO+ QTLs track insertion:
			# only occurs for Rat and Human
			if [[ ${organism} =~ ^[rR]at.* || ${organism} =~ ^[hH]uman.* ]]
			then
				./bin/flatfile-to-json.pl --gff "${currentWorkingOntPath}/${ontLabel}QTLs.gff3" --trackLabel "${ontLabel}QTLs" --key "${ontName} Related QTLs" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"color\" : \"#255ADF\", \"label\" : \"symbol,name\", \"description\" : \"\" }" --config "{ \"category\" : \"${ontCategory}\", \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
				echo "  +++ inserted track ${ontLabel}QTLs"
			fi
			
			# DO+ Strains track insertion:
			# only occurs for Rat
			if [[ ${organism} =~ ^[rR]at.* ]]
			then
				./bin/flatfile-to-json.pl --gff "${currentWorkingOntPath}/${ontLabel}Strains.gff3" --trackLabel "${ontLabel}Strains" --key "${ontName} Related Strains" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"color\" : \"#5ADF25\", \"label\" : \"symbol,name\", \"description\" : \"\" }" --config "{ \"category\" : \"${ontCategory}\", \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
				echo "  +++ inserted track ${ontLabel}Strains"
			fi
		
		elif [[ ${ontAccNum} =~ ^CHEBI.* ]]
		then

			# only rat, mouse and human have CHEBI tracks -- there is no (or very little) data for other species
			if [[ ${organism} =~ ^[rR]at.* || ${organism} =~ ^[hH]uman.* || ${organism} =~ ^[M]ouse.* ]] ; then

			  # Insert 1 CHEBI track per accession number
			  # CHEBI Genes track insertion
			  ./bin/flatfile-to-json.pl --gff "${currentWorkingOntPath}/${ontLabel}.gff3" --trackLabel "${ontLabel}" --key "${ontName} Related Genes" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"color\" : \"#DF255A\", \"label\" : \"symbol,name\", \"description\" : \"\" }" --config "{ \"category\" : \"${ontCategory}\", \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, \"menuTemplate\" : [ { \"label\" : \"View RGD details\", \"iconClass\" : \"dijitIconDatabase\", \"action\" : \"contentDialog\", \"content\" : \"<iframe src=\\\"${contentBlock}\\\" frameborder='0' marginheight='0' width='425' height='240' id='id_IFrame'></iframe>\", \"title\" : \"<center>RGD Feature Data for {name}</center>\" }, { \"label\" : \"Highlight this feature\" } ] }"
			  echo "  +++ inserted track ${ontLabel}"
			fi
		fi
	done < ${startingHome}/${2}
	fi
	
	# Generation of names index
	echo "=== generation of names index ==="
	if [[ "${pipelineSSLP}" != "null" ]]; then
		if [[ "${pipelineQTL}" != "null" ]]; then
			./bin/generate-names.pl --out "${datasetDir}" --tracks "ARGD_curated_genes,ARGD_curated_genesOnly,AQTLS,SSLP" --verbose --mem 4096000000
		else
			./bin/generate-names.pl --out "${datasetDir}" --tracks "ARGD_curated_genes,ARGD_curated_genesOnly,SSLP" --verbose --mem 4096000000
		fi
	else
		if [[ "${pipelineQTL}" != "null" ]]; then
			./bin/generate-names.pl --out "${datasetDir}" --tracks "ARGD_curated_genes,ARGD_curated_genesOnly,AQTLS" --verbose --mem 4096000000
		else
			./bin/generate-names.pl --out "${datasetDir}" --tracks "ARGD_curated_genes,ARGD_curated_genesOnly" --verbose --mem 4096000000
		fi
	fi

	popd
	printf "COMPLETE!\nMoved back to: ${PWD}\n\n"
	## End of JBrowse handling


done < ${1}

printf "***********************************\n"
printf 'Pipeline complete!\n'
date
printf "***********************************\n\n"

