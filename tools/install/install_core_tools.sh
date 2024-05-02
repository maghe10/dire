#install spades 3.15.4
# inspired by https://singularityhub.github.io/singularityhub-archive/containers/TomHarrop-assemblers-spades_3.14.1/
apt install cmake
apt install g++
apt install libbz2-dev
apt install zlib1g-dev

mkdir -p /root/sequencing/spades
cd /root/sequencing/spades
wget http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4.tar.gz
tar -xzf SPAdes-3.15.4.tar.gz
cd SPAdes-3.15.4
./spades_compile.sh

#wget http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5.tar.gz
#tar -xzf SPAdes-3.15.5.tar.gz
#cd SPAdes-3.15.5
#./spades_compile.sh

#trimgalore requires cutadapt and fastqc
apt install cutadapt
apt install fastqc
mkdir -p /root/sequencing/trimgalore
cd /root/sequencing/trimgalore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz

# checkm requires tools
# please use conda instead
#apt install hmmer
#apt install pplacer
#apt install prodigal


# for quast
apt install zlib1g-dev

#unicycler
mkdir -p /root/sequencing/unicycler
cd /root/sequencing/unicycler
wget https://github.com/rrwick/Unicycler/archive/refs/tags/v0.5.0.tar.gz
tar -xzf v0.5.0.tar.gz
#cd Unicycler-0.5.0/

