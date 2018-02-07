
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod ugo+x Miniconda3-latest-Linux-x86_64.sh 
bash ./Miniconda3-latest-Linux-x86_64.sh -b
source ~/.bashrc
export PATH=~/miniconda3/bin:$PATH 
conda update conda -y
pip install -U psutil
pip install -U snakemake
pip install -U PyYAML
pip install matplotlib
pip install pandas
pip install ipython
pip install numpy
conda config --add channels bioconda
conda install graphviz -y
#conda install retrying
mkdir ~/data_ref
echo " " >> ~/.bashrc
echo "PATH=~/miniconda3/bin:$PATH"  >> ~/.bashrc

