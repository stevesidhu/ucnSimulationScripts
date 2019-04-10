## This script combines several smaller scripts in order to automate the entire process from when
## PENTrack outputs several .out files to the checking, renaming, and moving to Angerona of root files


#!/bin/sh


## Firstly, the user is asked for the Runs they want to process

## Asks the user the name of the study

read -p $'\n\nEnter the name of the study: ' parameter

read -p $'\n\nYou must move the files to the toMerge folder to proceed. If you have not done that yet, would you like to now?
WARNING: This will overwrite ALL existing .out files of the same study name in the toMerge folder!! (y/n) ' moveToMerge

read -p $'\n\nIs this an emptying/storage simulation? (y/n) ' emptying

if [ $emptying = 'y' ]
then
		pathTop=${parameter}TopRun
		pathBot=${parameter}BottomRun
		
		if [ $moveToMerge = 'y' ]
		then
	        	numValues=$(ls -d -1q $pathTop* | wc -l)
		else
			numValues=$(ls -d -1q toMerge/$pathTop* | wc -l)
		fi
	else
		path=${parameter}Run

		if [ $moveToMerge = 'y'  ]
		then
	        	numValues=$(ls -d -1q $path* | wc -l)
		else
			numValues=$(ls -d -1q toMerge/$path* | wc -l)
		fi
fi

## This moves the Run files to a separate sub-folder which handles the merging processing of 
## completed simulations, enabling a separate simulation to be ran concurrently with merging

## First asks the user if the files have already been copied:

if [ $moveToMerge = 'y' ]
then
	printf "\n\nMoving files to the toMerge folder!\n\n"

	if [ $emptying = 'y' ]
	then
		for i in $(seq 1 $numValues)
		do
			rm -rf toMerge/$pathTop$i;
			rm -rf toMerge/$pathBot$i;
			mkdir -p toMerge/$pathTop$i;
			mkdir -p toMerge/$pathBot$i;
			mv $pathTop$i toMerge;
			mv $pathBot$i toMerge;
		done
	else
		for i in $(seq 1 $numValues)
		do
			rm -rf toMerge/$path$i;
			mkdir -p toMerge/$path$i;
			mv $path$i toMerge;
		done

	fi
fi

## Firstly, this loops through all specified Runs in order to combine the .out files into a .root file

## This will also look to check if there is already an out.root file in the Run folder. If so,
## it will prompt the user before overwriting it

## Confirm that user does in fact want to merge files

read -p $'\n\nDo you wish to merge .out files to a .root file? (y/n) ' answer

overwritePrompt=

if [ $answer = "y" ]
then

	## If there are already root files in the run folders, asks if user wants to overwrite them all

	if [ $emptying = "y" ]
	then
		for i in $(seq 1 $numValues)
		do
			cd toMerge/$pathTop$i

			if [ -e out.root ]
			then
				read -p $'\n\nSome runs already have out.root files. Would you like to overwrite all existing out.root files? (y/n) ' overwrite
				break
			fi

			cd ../..

			cd toMerge/$pathBot$i

			if [ -e out.root ]
			then
				read -p $'\n\nSome runs already have out.root files. Would you like to overwrite all existing out.root files? (y/n) ' overwrite
				break
			fi

			cd ../..

		done
	else
		for i in $(seq 1 $numValues)
		do
			cd toMerge/$path$i

			if [ -e out.root ]
			then
				read -p $'\n\nSome runs already have out.root files. Would you like to overwrite all existing out.root files? (y/n) ' overwrite
				break
			fi

			cd ../..

		done

	fi

	## Now do the merging, removing existing out.root files as necessary

	if [ $emptying = "y" ]
	then
		for i in $(seq 1 $numValues)
		do
			cd toMerge/$pathTop$i

			if [ -e out.root ] && [ $overwrite = "y" ]
			then
				rm out.root
			fi

			root -q -l ../../cScripts/merge_all.c

			printf "\n\n\nCompleted merging Top Run $i!\n\n\n"
			
			cd ../..
			
			cd toMerge/$pathBot$i

			if [ -e out.root ] && [ $overwrite = "y" ]
			then
				rm out.root
			fi

			root -q -l ../../cScripts/merge_all.c

			printf "\n\n\nCompleted merging Bottom Run $i!\n\n\n"
			
			cd ../..
		done

	else
		for i in $(seq 1 $numValues)
		do
			cd toMerge/$path$i
	
			if [ -e out.root ] && [ $overwrite = "y" ]
				then
					rm out.root
				else
					continue
			fi
				
			root -q -l ../../cScripts/merge_all.c
	
			printf "\n\n\nCompleted merging Run $i!\n\n\n"
			cd ../..
		done
	fi
fi

## This part automates the renaming of the root files


## Asks the user the name of the study (e.g. KinkHeightVariation)


## Loops through each Run, asking the user the value of the parameter in that Run

read -p $'\n\nWould you like to rename the root files now? (y/n) ' rename

if [ $rename = "y" ]
then
	if [ $emptying = "y" ]
	then
		for i in $(seq 1 $numValues);	
		do
				read -p "Enter the value of the parameter in Run $i: " val
				mv toMerge/$pathTop$i/out.root ${parameter}Top${val}.root
				mv toMerge/$pathBot$i/out.root ${parameter}Bottom${val}.root
		done
	else
		for i in $(seq 1 $numValues);	
		do
				read -p "Enter the value of the parameter in Run $i: " val
				mv toMerge/$path$i/out.root ${parameter}${val}.root
		done
	fi
fi

## This part quickly copies all root files to Angerona, and then moves the root files to a backup folder to reduce clutter

## Asks the user if they want to copy to Angerona right now

read -p $'\n\nWould you like to copy over to UCN cluster now? (y/n) ' proceed

if [ $proceed = 'y' ]
then
	scp ${parameter}*.root smorawetz@angerona.triumf.ca:/ucn/orithyia_data/ssidhu2/cryoSimulations
	
	mv ${parameter}*.root rootFileBackups

fi

## This part speeds up the sanity check of verifying the plots are reasonable


## The function then loops through all Runs specified by the user, checking the
## number of entries, in addition to drawing the plot for verification

ssh -Y smorawetz@angerona.triumf.ca 'bash -s' < checkPlots.sh
