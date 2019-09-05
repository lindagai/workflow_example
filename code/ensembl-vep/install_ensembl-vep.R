################################################################################

#1_installing EnsemblVEP

################################################################################

#NOTE: All code below should be sent to the Terminal, not R.

################################################################################

#0. Sign into JHPCE

ssh -X lgai@jhpce01.jhsph.edu
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=16G

#################################################################################

#1. Set up conda

module load conda

#Add
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#Check version
conda -V

#Output:
#conda 4.6.14

#################################################################################

#2. Set up perl environment

#You will need to follow the directions here before running the code below.
#https://jhpce.jhu.edu/knowledge-base/community/perl/

#Load packages that allow you to install perl modules locally
#You do not have permission to install perl modules globally!
module load JHPCE_CENTOS7_DEFAULT_ENV
module load gcc/4.4.7
module load .perl/5.10.1
module load myperl

#################################################################################

#3. Check if required perl libraries are installed by looking for filepaths

module load perl
perl -V

#Output:
#/jhpce/shared/jhpce/core/perl/5.24.4/lib/5.24.4

#Check if required perl modules are installed by looking for filepaths
#If using perl 5.10.1 on JHPCE, they should already be installed
#THEY DO NOT QORK
pmpath Archive::Zip
pmpath DBI
pmpath DBD::mysql

#Load perl modules needed to install ensembl-VEP
module load htslib

#Install perl modules needed to install ensembl-VEP locally
cpanm --local-lib $MYPERL_INSTALL_BASE Bio::DB::HTS --force

#################################################################################

#4. Install ensembl-VEP

#A. Clone conda into your directory

# You do not have the necessary permissions to install packages
# into the install area '/jhpce/shared/jhpce/core/conda/miniconda-3-4.6.14'.
# So you will have to clone this environment into your home directory and
# then make changes to it.
# This may be done using the command:
conda create -n my_root --clone="/jhpce/shared/jhpce/core/conda/miniconda3-4.6.14"

#Output
# NotWritableError: The current user does not have write permissions to a required path.
  # path: /jhpce/shared/jhpce/core/conda/miniconda-3/envs/.conda_envs_dir_test
  # uid: 41693
  # gid: 100

# If you feel that permissions on this path are set incorrectly, you can manually
# change them by executing

  # $ sudo chown 41693:100 /jhpce/shared/jhpce/core/conda/miniconda-3/envs/.conda_envs_dir_test

# In general, it's not advisable to use 'sudo conda'.

#To activate this environment, use:
source activate my_root

#B. Install ensembl-VEP and update it
conda install ensembl-vep
conda update ensembl-vep

#To deactivate an active environment, use:
#source deactivate

#################################################################################

#5. Install ensembl-VEP

conda install ensembl-vep

#Install ensembl-vep and install

#Do not try to update conda even if you get a warning that it's out of date
#You do not have write-permission to update conda

#################################################################################

#6. Install hg38 VEP library

#Fix this!

vep_install -a cf -s homo_sapiens -y GRCh38 -c /output/path/to/GRCh38/vep --CONVERT

# Version check reports a newer release of ensembl-vep is available (installed: 92, available: 97)

# You should exit this installer and re-download ensembl-vep if you wish to update

# Do you wish to exit so you can get updates (y) or continue (n): n
# OK, bye!

# NB: Remember to re-run INSTALL.pl after updating to check for API updates

vep_install -a cf -s homo_sapiens -y GRCh38 -c /output/path/to/GRCh38/vep --CONVERT

# Version check reports a newer release of ensembl-vep is available (installed: 92, available: 97)

# You should exit this installer and re-download ensembl-vep if you wish to update

# Do you wish to exit so you can get updates (y) or continue (n): n
# ERROR: Could not create directory /output/path/to/GRCh38/vep
