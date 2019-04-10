## This is a shorter script which deletes all of the stuff (e.g. batch files for submission, .out files, etc.)
## once the simulations are complete and the plots have been inspected to ensure correctness



## Asks the user the name of the study

read -p $'\n\nEnter the name of the study: ' parameter

read -p $'\n\nIs this an emptying/storage study? (y/n) ' emptying

## Finally, this deletes the Run folders containing STL files on the home drive, and the batch and .out files

read -p $'\n\nThis will delete the .out files, the batch files used to submit jobs to cedar, and the text files which record job IDs.
Are you sure you are ready to delete all this? (y/n) ' proceed

path=${parameter}Run
pathTop=${parameter}TopRun
pathBot=${parameter}BottomRun


if [ $emptying = 'y' ]
then
	numValues=$(ls -d -1q toMerge/$pathTop* | wc -l)
else
	numValues=$(ls -d -1q toMerge/$path* | wc -l)
fi

if [ $proceed = 'y' ]
then
	for i in $(seq 1 $numValues)
	do
		rm -rf toMerge/$path$i
		rm -rf $path$i
		rm -rf toMerge/$pathTop$i
		rm -rf $pathTop$i
		rm -rf toMerge/$pathBot$i
		rm -rf $pathBot$i
	done
	
	mv $parameter*.root rootFileBackups

	cd ..
	
	for i in $(seq 1 $numValues)
	do
		rm -rf $path$i
		rm -rf $pathTop$i
		rm -rf $pathBot$i
		rm -rf ${parameter}Batch$i.sh
		rm -rf ${parameter}TopBatch$i.sh
		rm -rf ${parameter}BottomBatch$i.sh
	done
	
	rm ${parameter}.txt
fi
