#$ -S /bin/sh
#$ -cwd			        #run in current working directory
#$ -j y
#$ -N radCHARRAD   		# Name of job
##$ -pe mpich 10 		# Number of CPUs
#$ -P OzelGroup		    #  
##$ -l h=compute-0-8    #

# Libraries
source ~/.bashrc

# Source Anaconda
export PATH="/home/ao57/anaconda3/bin:$PATH"

# Case folder
caseFolder=$PWD
#caseFolderName=$(basename $caseFolder)
caseFolderName=$(echo $caseFolder | awk -F "/" '{print $(NF-1)""$NF}')

# Go to tmp folder
cd /tmp
if [ ! -d $USER ]; then
  mkdir $USER
fi
cd $USER
if [ -d $caseFolderName ]; then
  rm -rf $caseFolderName
fi

mkdir $caseFolderName
rsync -a $caseFolder/ ./$caseFolderName/

# Go to case folder
cd $caseFolderName

# Run case
python /home/ao57/projects/pyLBM/src/pyLBM.py > logRun

cd ..
rsync -uavz ./$caseFolderName/ $caseFolder/
rm -rf $caseFolderName

