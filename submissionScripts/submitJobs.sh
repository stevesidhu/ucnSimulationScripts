#!/bin/sh

## This is a script to automate the process of submitting jobs, and allows additional functionality
## in naming the output directories, allowing for more simulations to be ran at once without confusion

## First prompts the user for what they have renamed the study, and how many different geometries are within it

read -p $'\n\nWhat have you re-named the study? ' parameter

path=${parameter}Run
pathTop=${parameter}TopRun
pathBot=${parameter}BottomRun

read -p $'\n\nIs this an emptying/storage simulation? (y/n) ' emptying

if [ $emptying = 'y' ]
then
	numGeometries=$(ls -d -1q $pathTop* | wc -l)
else
	numGeometries=$(ls -d -1q $path* | wc -l)
fi

## This loops through and creates a new batch file for each separameter variation of the parameter in the study
## It will set a doubled walltime if the simulation is emptying

if [ $emptying = 'y' ]
then
	for i in $(seq 1 $numGeometries)
	do
	
		## First creates a batch file for the top cell emptying

		cat > ${parameter}TopBatch$i.sh <<-EOF
		#!/bin/sh
				
		#SBATCH --time=200
		#SBATCH --mem=4G
		#SBATCH --array=1-50
		#SBATCH --output=/home/stewmo/slurmOutput/$pathTop$i.%A_%a.out
		#SBATCH --error=/home/stewmo/slurmErrors/$pathTop$i.%A_%a.out
		#PBS -l walltime=12:00:00
		
		ID=\$SLURM_ARRAY_TASK_ID
		JOB=\$SLURM_ARRAY_JOB_ID
		
		/home/stewmo/PENTrack/PENTrack \$JOB\$ID /home/stewmo/$pathTop$i/config.in /home/stewmo/scratch/$pathTop$i
		EOF

		## Prior to submitting them, new folders in the scratch drive need to be created to accomodate the output files
	
		mkdir -p scratch/$pathTop$i


		## Repeat the same process, but for the bottom cell emptying

		cat > ${parameter}BottomBatch$i.sh <<-EOF
		#!/bin/sh
				
		#SBATCH --time=200
		#SBATCH --mem=4G
		#SBATCH --array=1-50
		#SBATCH --output=/home/stewmo/slurmOutput/$pathBot$i.%A_%a.out
		#SBATCH --error=/home/stewmo/slurmErrors/$pathBot$i.%A_%a.out
		#PBS -l walltime=12:00:00
		
		ID=\$SLURM_ARRAY_TASK_ID
		JOB=\$SLURM_ARRAY_JOB_ID
		
		/home/stewmo/PENTrack/PENTrack \$JOB\$ID /home/stewmo/$pathBot$i/config.in /home/stewmo/scratch/$pathBot$i
		EOF

		## Prior to submitting them, new folders in the scratch drive need to be created to accomodate the output files
	
		mkdir -p scratch/$pathBot$i
	
	done
else
	for i in $(seq 1 $numGeometries)
	do

		cat > ${parameter}Batch$i.sh <<-EOF
		#!/bin/sh
				
		#SBATCH --time=200
		#SBATCH --mem=4G
		#SBATCH --array=1-50
		#SBATCH --output=/home/stewmo/slurmOutput/$path$i.%A_%a.out
		#SBATCH --error=/home/stewmo/slurmErrors/$path$i.%A_%a.out
		#PBS -l walltime=12:00:00
	
		ID=\$SLURM_ARRAY_TASK_ID
		JOB=\$SLURM_ARRAY_JOB_ID
	
		/home/stewmo/PENTrack/PENTrack \$JOB\$ID /home/stewmo/$path$i/config.in /home/stewmo/scratch/$path$i
		EOF

		## Prior to submitting them, new folders in the scratch drive need to be created to accomodate the output files

		mkdir -p scratch/$path$i

	done

fi


## Once there is a batch file for each, it loops through all of them and submits the jobs, while writing
## the jobs IDs from Slurm to a text file to make it easier to keep track of them

for i in $(seq 1 $numGeometries)
do

	if [ $emptying = "y" ]
	then
		sbatch --account=def-rpicker ${parameter}TopBatch$i.sh >> temp${parameter}Top.txt
		sbatch --account=def-rpicker ${parameter}BottomBatch$i.sh >> temp${parameter}Bottom.txt
		
		## Brief bit to keep only batch IDs in text files
	
		grep -Eo '[0-9]{1,20}' temp${parameter}Top.txt >> ${parameter}Top.txt
		grep -Eo '[0-9]{1,20}' temp${parameter}Bottom.txt >> ${parameter}Bottom.txt
		rm temp${parameter}Top.txt
		rm temp${parameter}Bottom.txt
	
	else
		sbatch --account=def-rpicker ${parameter}Batch$i.sh >> temp$parameter.txt 
		grep -Eo '[0-9]{1,20}' temp$parameter.txt >> ${parameter}.txt
		rm temp$parameter.txt
	fi
done


