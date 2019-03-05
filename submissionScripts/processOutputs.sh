#!/bin/sh



## This script combines several smaller scripts in order to automate the entire process from when
## PENTrack outputs several .out files to the checking, renaming, and moving to Angerona of root files



## Makes sure the user doesn't have anything in the toMerge folder they don't want to overwrite

printf "\n\nWARNING: This will overwrite existing run folders in the toMerge folder! Type 'y' (no quotes)
to proceed, or type anything if you don't want this to happen and the script will exit! "

read goAhead

if [ $goAhead != "y" ]
then
	exit
fi

## Firstly, the user is asked for the runs they want to process

## Asks the user the name of the study

read -p $'\n\nEnter the name of the study: ' parameter

read -p $'\n\nWhat is the first run you want to process? ' firstRun
read -p $'\n\nWhat is the last run you want to process? ' lastRun

## This moves the run files to a separate sub-folder which handles the merging processing of 
## completed simulations, enabling a separate simulation to be ran concurrently with merging

## First asks the user if the files have already been copied:

read -p $'\n\nHave the files already been moved to the toMerge folder? If you want to move them, type "y" (no quotes).
Otherwise, type anything if they are already located there. ' alreadyDone

if [ $alreadyDone = 'y' ]
then
	printf "\n\nMoving files to the toMerge folder!\n\n"

	for i in $(seq $firstRun $lastRun)
		do
			rm -rf toMerge/${parameter}run$i;
			mkdir -p toMerge/${parameter}run$i;
			mv ${parameter}run$i toMerge;
		done
fi


## Firstly, this loops through all specified runs in order to combine the .out files into a .root file

## This will also look to check if there is already an out.root file in the run folder. If so,
## it will prompt the user before overwriting it

## Confirm that user does in fact want to merge files

read -p $'\n\nDo you want to merge .out files to a .root file? If so, type "y" (no quotes). 
If they are already merged, and you want to overwrite them all at once, type "w" (no quotes). 
If they are already merged, and you don\'t want to overwrite them - or overwrite them separately, type anything: ' answer

if [ $answer = "y" ] || [ $answer = "w" ]
then
	for i in $(seq $fistRun $lastRun);
		do
			cd toMerge/${parameter}run$i
	
			if [ -e out.root ] && [ $answer != "w" ]
			then
				read -p "There is already an out.root file for Run $i. Type 'w' (no quotes) to overwrite, or type anything to leave current out.root file: " overwrite
				if [ $overwrite = 'w' ]
				then
					rm out.root
				else
					continue
				fi
			fi
			
			root -q -l ../../cScripts/merge_all.c

			printf "\n\n\nCompleted merging $(basename $i)!\n\n\n"
			cd ../..;
		done
fi


## This part automates the renaming of the root files


## Asks the user the name of the study (e.g. KinkHeightVariation)


## Loops through each run, asking the user the value of the parameter in that run

read -p $'\n\nHave the root files been renamed yet? If you want to rename them then, type "y" (no quotes). If not, type anything: ' proceed

printf "\n\nThe current study is $parameter!\n\n\n"

if [ $proceed = 'y' ]
then
	for i in $(seq $firstRun $lastRun);	
	do
		read -p "Enter the value of the parameter in Run $i: " val;
		mv toMerge/${parameter}run$i/out.root ${parameter}${val}.root;
	done
fi




## This part quickly copies all root files to Angerona, and then moves the root files to a backup folder to reduce clutter

## Asks the user if they want to copy to Angerona right now

read -p $'\n\nIf you want to copy over to UCN cluster now, type "y" (no quotes). If not, type anything. ' proceed

if [ $proceed = 'y' ]
then
	scp ${parameter}*.root smorawetz@angerona.triumf.ca:/ucn/orithyia_data/ssidhu2/cryoSimulations
fi

## This part is a short script called when ssh-ing into Angerona to check the plots.
## This must be done on Angerona because Cedar has issues with plotting in root.

read -p $'\n\nWould you like to view the zx neutronend plot? Type "y" (no quotes) to do so, or type anything to skip this. ' proceed

if [ $proceed = 'y' ]
then
	ssh -Y smorawetz@angerona.triumf.ca 'bash -s' < checkPlots.sh $parameter
fi
