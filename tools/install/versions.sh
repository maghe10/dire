# use with  &> filename.log
set -x #echo on

echo "environment"
python --version
conda --version

echo ""
echo "raw reads quality programs"
confindr --version
/root/sequencing/trimgalore/TrimGalore-0.6.10/trim_galore --version
cutadapt --version
fastqc --version

echo ""
echo "assembly"
python /root/sequencing/spades/SPAdes-3.15.4/bin/spades.py --version
spades.py --version
unicycler --version
seqkit version

echo ""
echo "post assembly tools"
quast  --version
multiqc --version
checkm | head -n 3
amrfinder --version


