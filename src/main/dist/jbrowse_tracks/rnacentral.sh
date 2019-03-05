# MASTER SCRIPT to refresh all JBrowse tracks from RNAcentral
#
cd /data/data/gff3
mkdir rnacentral
rm rnacentral/*.gff3

./rnacentral_load.sh pan_paniscus.panpan1.1.gff3       '/jbrowse/data_bonobo1_1'
./rnacentral_load.sh sus_scrofa.Sscrofa11.1.gff3       '/jbrowse/data_pig11_1'
./rnacentral_load.sh rattus_norvegicus.Rnor_6.0.gff3   '/jbrowse/data_rgd6'
./rnacentral_load.sh mus_musculus.GRCm38.gff3          '/jbrowse/data_mm38'
./rnacentral_load.sh homo_sapiens.GRCh38.gff3          '/jbrowse/data_hg38'
./rnacentral_load.sh canis_familiaris.CanFam3.1.gff3   '/jbrowse/data_dog3_1'
