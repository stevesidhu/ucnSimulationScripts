## Informs user about folder-naming restrictions

echo "Ensure that all folders are within the study are named beginning with the value of their parameter. For example, 35cmKinkHeight."

## Asks the user for the name of the study

read -p $'\n\nWhat is the name of the study? ' studyName

read -p $'\n\nWhat would you like to rename the study? ' parameter

read -p $'\n\nIs this an emptying or storage study? If so, type "y" (no quotes). ' emptying

if [ $emptying = "y" ]
then
	printf "\n\nTwo separate folders, ${study}Top and ${study}Bottom will be created.\n\n\n"
fi

## Loop through all of the folders within current study, removing
## spaces in the file names

cd $studyName

for j in $(ls -d *);
	do
		cd $j;
		for f in *\ *; do mv "$f" "${f// /}"; done
		cd ..;
	done

## Loop through all of the different folders of STL files which contain the
## geometries for different values of the parameter being varied.

## Requires an index i to keep track of the run folders

i=1

## Move into chosen study to retrieve a list of all the folders to loop through

for j in $(ls -d -v *);
do
	if [ $emptying != 'y' ]
	then
		rm -rf ../${parameter}Run$i/STL;
		mkdir -p ../${parameter}Run$i/STL;
		echo "Currently copying $j"
		cp -r $j/* ../${parameter}Run$i/STL;
	else
		rm -rf ../${parameter}*Run$i/STL;
		mkdir -p ../${parameter}TopRun$i/STL;
		mkdir -p ../${parameter}BottomRun$i/STL;
		echo "Currently copying $j"
		cp -r $j/* ../${parameter}TopRun$i/STL;
		cp -r $j/* ../${parameter}BottomRun$i/STL;
	fi
	i=$((i+1));
done
