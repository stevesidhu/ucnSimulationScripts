## This is a shorter script which deletes all of the stuff (e.g. batch files for submission, .out files, etc.)
## once the simulations are complete and the plots have been inspected to ensure correctness



## Asks the user the name of the study

read -p $'\n\nEnter the name of the study: ' parameter

read -p $'\n\nWhat is the first Run you want to delete? ' firstRun
read -p $'\n\nWhat is the last Run you want to delete? ' lastRun



## Finally, this deletes the Run folders containing STL files on the home drive, and the batch and .out files

read -p $'\n\nThis will delete the .out files, the batch files used to submit jobs to cedar, and the text files which record job IDs.
Are you sure you are ready to delete all this? Type "y" (no quotes) to delete it all. ' proceed

if [ $proceed = 'y' ]
then
	for i in $(seq $firstRun $lastRun)
	do
		rm -r toMerge/${parameter}Run$i;
	done
	
	mv $parameter*.root rootFileBackups

	cd ..
	
	for i in $(seq $firstRun $lastRun)
	do
		rm -r ${parameter}Run$i;
		rm ${parameter}Batch$i.sh;
	done
	
	rm ${parameter}.txt
fi
