conda create -n dire python=3.11.5
conda activate dire

#trimgalore => not in conda
#spades => not in conda

#unicycler, install with

conda install -c conda-forge python
conda install -c bioconda seqkit
conda install -c bioconda ncbi-amrfinderplus


pip3 install /root/sequencing/unicycler/Unicycler-0.5.0

# The following does not install under python 3.11+ 
#conda install -c conda-forge gsl
#conda install -c bioconda confindr
#conda install -c bioconda quast

