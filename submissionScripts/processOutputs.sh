## This script combines several smaller scripts in order to automate the entire process from when
## PENTrack outputs several .out files to the checking, renaming, and moving to Angerona of root files


#!/bin/sh



## Makes sure the user doesn't have anything in the toMerge folder they don't want to overwrite

printf "\n\nWARNING: This will overwrite existing Run folders in the toMerge folder! Type 'y' (no quotes)
to proceed, or type anything if you don't want this to happen and the script will exit! "

read goAhead

if [ $goAhead != "y" ]
then
	exit
fi

## Firstly, the user is asked for the Runs they want to process

## Asks the user the name of the study

read -p $'\n\nEnter the name of the study: ' parameter

read -p $'\n\nWhat is the first Run you want to process? ' firstRun
read -p $'\n\nWhat is the last Run you want to process? ' lastRun

## This moves the Run files to a separate sub-folder which handles the merging processing of 
## completed simulations, enabling a separate simulation to be ran concurrently with merging

## First asks the user if the files have already been copied:

read -p $'\n\nHave the files already been moved to the toMerge folder? If not, type "n" (no quotes).
Otherwise, type anything if they are already located there. ' alreadyDone

if [ $alreadyDone = 'n' ]
then
	printf "\n\nMoving files to the toMerge folder!\n\n"

	for i in $(seq $firstRun $lastRun)
		do
			rm -rf toMerge/${parameter}Run$i;
			mkdir -p toMerge/${parameter}Run$i;
			mv ${parameter}Run$i toMerge;
		done
fi


## Firstly, this loops through all specified Runs in order to combine the .out files into a .root file

## This will also look to check if there is already an out.root file in the Run folder. If so,
## it will prompt the user before overwriting it

## Confirm that user does in fact want to merge files

read -p $'\n\nDo you want to merge .out files to a .root file? If so, type "y" (no quotes). 
If they are already merged, and you want to overwrite them all at once, type "w" (no quotes). 
If they are already merged, and you don\'t want to overwrite them - or overwrite them separately, type anything: ' answer

if [ $answer = "y" ] || [ $answer = "w" ]
then
	for i in $(seq $fistRun $lastRun);
		do
			cd toMerge/${parameter}Run$i
	
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


## Loops through each Run, asking the user the value of the parameter in that Run

for i in $(seq $firstRun $lastRun);	
do
		read -p "Enter the value of the parameter in Run $i: " val;
		mv toMerge/${parameter}Run$i/out.root ${parameter}${val}.root;
done



## This part deletes the .out files as they are unnecessary after they are merged into a .root file


## Asks the user prior to looping through and deletes all the Run folders in the toMerge

read -p $'\n\nNow that you have the root files, are you ready to delete the Run folders with the .out files? Type "y" (no quotes) if so. ' proceed

if [ $proceed = 'y' ]
then
	for i in $(seq $firstRun $lastRun)
	do
		rm -r toMerge/${parameter}Run$i;
	done
fi


## This part quickly copies all root files to Angerona, and then moves the root files to a backup folder to reduce clutter

## Asks the user if they want to copy to Angerona right now

read -p $'\n\nIf you want to copy over to UCN cluster now, type "y" (no quotes). If not, type anything. ' proceed

if [ $proceed = 'y' ]
then
	scp *.root smorawetz@angerona.triumf.ca:/ucn/orithyia_data/ssidhu2/cryoSimulations
	
	mv *.root rootFileBackups

fi

## This part deletes the .out files as they are unnecessary after they are merged into a .root file


## Asks the user prior to looping through and deletes all the Run folders in the toMerge

read -p $'\n\nNow that you have the root files, are you ready to delete the Run folders with the .out files? Type "y" (no quotes) if so. ' proceed

if [ $proceed = 'y' ]
then
	for i in $(seq $firstRun $lastRun)
fi


## Finally, this deletes the Run folders containing STL files on the home drive, and the batch files

read -p $'\n\nNow that the .root files are on the UCN cluster and have been verified, are you ready to delete the STL,
output and batch files? Type "y" (no quotes) if so. ' proceed

if [ $proceed = 'y' ]
then
	cd ..

	rm -r ${parameter}Run*

	rm ${parameter}batch*.sh
fi

## This part speeds up the sanity check of verifying the plots are reasonable


## The function then loops through all Runs specified by the user, checking the
## number of entries, in addition to drawing the plot for verification

exit

ssh -Y smorawetz@angerona.triumf.ca

cd ../../ucn/orithyia_data/ssidhu2/cryoSimulations

for i in $(ls ${parameter}*.root); 
do
	echo "Currently checking Run $i"	
	root -q -l "cScripts/outputCheck.C(\"$i\")"
	root -l "cScripts/plotCheck.C(\"$i\")" <<-EOF
	std::cout << "\n\nClose the plot to proceed\n\n" << endl;
	gPad->WaitPrimitive();
	.q
	EOF
done

## ## The following is no longer necessary as verifying the plots is the last step
## 
## read -p '\n\n\nIMPORTANT: To proceed, type "y" (no quotes). Otherwise, if any Runs are
## suspicious, type anything and this program will exit - you can repeat this process again
## after the problems are fixed! ' goForward
## 
## if [ $goForward != "y" ]
## then
## 	exit
## fi
