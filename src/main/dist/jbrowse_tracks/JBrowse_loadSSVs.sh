# Assignment of the selected assembly's
# GFF3 directory
GFF3home="${PWD}"

# Main execution block
for file in *.gff3*
do
        # Variable assignments
        fileLocation="${GFF3home}/${file}"
        trackName="$(printf "${file}" | sed 's/_/\//' | sed 's/\.gff3.*//')"
        trackLabel="SSV_$(printf "${file}" | sed 's/\.gff3.*//')"

        # DEBUG - sanity check
        # printf "${fileLocation}\n${trackName}\n${trackLabel}\n\n"

        # Check if the track is GZipped--if it is, unzip the file
        if [[ ${file} =~ ".*\.gff3\.gz$" ]]
        then
                printf "Unzipping ${file} ... "
                gzip -d "${file}"
                printf 'done!\n'
        fi

        # Move to the JBrowse home directory
        pushd ${jbrowseHome}

        # Remove any pre-existing SSV track with an identical track label
        ./bin/remove-track.pl --dir "${dataset}" --trackLabel "${trackLabel}" --delete

        # Insert the new SSV track
        ./bin/flatfile-to-json.pl --gff "${fileLocation%.gz}" --trackLabel "${trackLabel}" --key "${trackName}" --out "${dataset}" --trackType JBrowse/View/Track/CanvasFeatures --clientConfig "{ \"className\" : \"feature\", \"color\" : \"function( feature ) { var pred = feature.get( \\\"polypred\\\" ); if( pred == null ) pred = \\\"NA\\\"; else pred = decodeURIComponent( pred.toString() ).split( \\\",\\\" )[ 0 ].split( \\\"||\\\" )[ 4 ]; switch( pred ) { case \\\"probably damaging\\\": return \\\"#DF255A\\\"; case \\\"possibly damaging\\\": return \\\"#C4C400\\\"; case \\\"benign\\\": return \\\"#378D15\\\"; case \\\"NA\\\": default: return \\\"#000000\\\"; } }\", \"label\" : \"symbol,name\", \"featureScale\" : 0.0001, \"description\" : \"\" }" --config "{ \"category\" : \"Variants/Strain Specific Variants\", \"histograms\" : { \"binsPerBlock\" : 25, \"color\" : \"#6600CC\" }, \"onClick\" : { \"iconClass\" : \"dijitIconDatabase\", \"content\" : \"<iframe src=\\\"/rgdweb/front/detail_gbrowse.html?vid={name}&mapKey=${mapKey}\\\" frameborder='0' marginheight='0' scrolling='yes' width='680' height='315'></iframe>\", \"action\" : \"contentDialog\", \"title\" : \"<center>RGD Strain-Specific Data for Variant #{name}</center>\" }, \"menuTemplate\" : [ { \"iconClass\" : \"dijitIconDatabase\", \"content\" : \"<iframe src=\\\"/rgdweb/front/detail_gbrowse.html?vid={name}&mapKey=${mapKey}\\\" frameborder='0' marginheight='0' scrolling='yes' width='680' height='315'></iframe>\", \"action\" : \"contentDialog\", \"title\" : \"<center>RGD Strain-Specific Data for Variant #{name}</center>\", \"label\" : \"View RGD Details\" }, { \"label\" : \"Highlight this feature\" } ] }"

        # Move back to the specific assembly
        # directory housing this set of SSVs
        popd

        # Re-zip the SSV *.gff3 file
        #
        # This should probably be re-factored
        # into an entirely new variable and used
        # both here and in the track-insertion
        # clause above!
        printf "GZipping ${file%.gz} ... "
        gzip "${file%.gz}"
        printf 'done!\n'
done

# Final popd to return to the *general*
# SSVs housing directory
popd

