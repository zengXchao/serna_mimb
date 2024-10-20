curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
if [ -f Miniforge3-Linux-x86_64.sh ]; then
	rm -rf Miniforge3-Linux-x86_64.sh
fi
mamba install -c bioconda samtools star stringtie subread bedtools
mamba install r
chmod +x ${HOME}/miniforge3/bin/*
chmod +x ${HOME}/miniforge3/lib/R/bin/*
