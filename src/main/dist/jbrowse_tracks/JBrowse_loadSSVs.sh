#!/bin/bash
# load Strain Specific Variant tracks into rgd6, rgd5 and rgd3_4 assemblies

###### VARIABLE DECLARATIONS ######
jbrowseHome="/rgd/JBrowse-1.16.3/"


workingDir="pipelineOutput"
gff3Dir="strain_specific_variants"
reedPulldown="scp -q rgdpub@reed.rgd.mcw.edu:"

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
if [ "$SERVER" = "REED" ]; then
  reedPulldown="cp -f "
fi

mainGFF3home="/rgd/data/gff3"
mainPipelineOutput="/home/rgddata/pipelines/RGDGff3Pipeline/data/strain_specific_variants"

###################################

assemblies=( "6" "5" "3_4" )

for ASSEMBLY in "${assemblies[@]}"
do

        mainOrgDir="Rat/rn${ASSEMBLY}"
        datasetDir="data_rgd${ASSEMBLY}"
        organism=Rat
        assembly="data_rgd${ASSEMBLY}"


        printf "***********************************\n"
        printf "Beginning processing of Strain Specific Variants for assembly ${ASSEMBLY} ...\n"
        date
        printf "***********************************\n\n"

        ## Initial path setting
        startingHome="${PWD}"
        currentOrgPath="${mainGFF3home}/${mainOrgDir}"
        currentWorkingPath="${currentOrgPath}/${workingDir}"
        currentgff3Path="${currentOrgPath}/${gff3Dir}"
        printf "startingHome: ${startingHome}\n"
        printf "currentOrgPath: ${currentOrgPath}\n"
        printf "currentWorkingPath: ${currentWorkingPath}\n"
        printf "currentgff3Path: ${currentgff3Path}\n"

        ## Initial directory setup
       printf "Initiating directory setup ... "
      if [[ ! -d ${currentWorkingPath} ]]; then
               mkdir -p ${currentWorkingPath}
       fi
       if [[ ! -d ${currentgff3Path} ]]; then
               mkdir -p ${currentgff3Path}
       fi
       printf 'done!\n\n'


        ## Initial clean-up
       printf "Initiating clean-up ... "
       rm -f  ${currentWorkingPath}/*.gff3* ${currentgff3Path}/*.gff3*
       printf 'done!\n\n'

        #pulldown from REED
       echo PULLDOWN from reed "${reedPulldown}${mainPipelineOutput}/rat${ASSEMBLY}/*.gff3.gz ${currentgff3Path}"
       echo
       ${reedPulldown}${mainPipelineOutput}/rat${ASSEMBLY}/*.gff3.gz ${currentgff3Path}
       echo "PULLDOWN ok"


        ## Unzip imported GFF3 files
        printf "Unzipping GFF3 files ... "
        for file in ${currentgff3Path}/*.gff3.gz
        do
                gzip -d "${file}"
        printf "File: ${file}\n"

        trackName="$(printf "${file##*/}" | sed 's/_/\//' | sed 's/\.gff3.*//')"
        trackLabel="SSV_$(printf "${file##*/}" | sed 's/\.gff3.*//')"
        fileName="$(printf "${file##*/}" | sed 's/\.gz//' )"
        printf "trackLabel: ${trackLabel}\n"
        printf "fileName: ${fileName}\n"
        ## JBrowse handling - where the magic happens!
        pushd ${jbrowseHome}
        printf "\nMoved to JBrowse home directory: ${PWD}\n\n"

        # Purging old track
        echo "remove ${trackLabel} track"
        ./bin/remove-track.pl --dir "${datasetDir}" --trackLabel "${trackLabel}" --delete

        echo "insert Variants track"
        ./bin/flatfile-to-json.pl --gff "${currentgff3Path}/${fileName}" --trackLabel "${trackLabel}" --key "${trackName}" --out "${datasetDir}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"className\" : \"feature\", \"color\" : \"function( feature ) { var pred = feature.get( \\\"polypred\\\" ); if( pred == null ) pred = \\\"NA\\\"; else pred = decodeURIComponent( pred.toString() ).split( \\\",\\\" )[ 0 ].split( \\\"||\\\" )[ 4 ]; switch( pred ) { case \\\"probably damaging\\\": return \\\"#DF255A\\\"; case \\\"possibly damaging\\\": return \\\"#C4C400\\\"; case \\\"benign\\\": return \\\"#378D15\\\"; case \\\"NA\\\": default: return \\\"#000000\\\"; } }\", \"label\" : \"symbol,name\", \"featureScale\" : 0.0001, \"description\" : \"\" }" --config "{ \"category\" : \"Variants/Strain Specific Variants\", \"histograms\" : { \"binsPerBlock\" : 25, \"color\" : \"#6600CC\" }, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"content\" : \"<iframe src=\\\"/rgdweb/front/detail_gbrowse.html?vid={name}&mapKey=${mapKey}\\\" frameborder='0' marginheight='0' scrolling='yes' width='680' height='315'></iframe>\", \"action\" : \"contentDialog\", \"title\" : \"<center>RGD Strain-Specific Data for Variant #{name}</center>\" }, \"menuTemplate\" : [ { \"iconClass\" : \"dijitIconDatabase\", \"content\" : \"<iframe src=\\\"/rgdweb/front/detail_gbrowse.html?vid={name}&mapKey=${mapKey}\\\" frameborder='0' marginheight='0' scrolling='yes' width='680' height='315'></iframe>\", \"action\" : \"contentDialog\", \"title\" : \"<center>RGD Strain-Specific Data for Variant #{name}</center>\", \"label\" : \"View RGD Details\" }, { \"label\" : \"Highlight this feature\" } ] }"
        echo "+++ inserted track Variants"


        popd
        printf "COMPLETE!\nMoved back to: ${PWD}\n\n"
        ## End of JBrowse handling
        done
done
