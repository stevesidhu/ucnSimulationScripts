## Make sure the user has the desired file to copy in this directory

echo "Make sure you have the file you want to copy in this directory!"


## Asks the user how many different geometries - i.e., values of the
## parameter in the study - there are, and asks what the study is called

read -p $'\n\nWhat name have you assigned the study? ' parameter

## Asks the user if this is an emptying simulation, if so, if they'd like to copy the emptying config files

read -p $'\n\nIs this an emptying simulation? Type "y" (no quotes) if so. ' emptying

## Determines the number of different geometries from $parameter so that it knows how many times to loop

if [ $emptying = "y"]
then
	numGeometries=$(ls -d -1q ${parameter}Toprun* | wc -l)
else
	numGeometries=$(ls -d -1q ${parameter}run* | wc -l)
fi

if [ $emptying = "y" ]
then
	read -p $'\n\nIf you would like to copy the emptying config files, type "y" (no quotes). ' emptConfig
	
	if [ $emptConfig = "y" ]
	then
		read -p $'\n\nWhat is the name of the config file for the top cell emptying? ' topEmpt
		read -p $'\n\nWhat is the name of the config file for the bottom cell emptying? ' botEmpt
	fi

fi

## Asks the user what file they want to copy

fileToCopy="placeholder"

while [ $fileToCopy != "done" ]
do

	read -p $'\n\nWhat file would you like to copy now? If finished, type 'done' (no quotes) ' fileToCopy
	
	if [ $fileToCopy != "done" ]
	then
		for i in $(seq 1 $numGeometries);
		do
		
			## Makes a new directory for each run without overwriting pre-existing
			## directory, and copies user-specified file to it
			
			if [ $emptying != 'y' ]
			then
				mkdir -p ${parameter}run$i;
				cp $fileToCopy ${parameter}run$i;
			else
			
				## Copies the earlier specified config files into the run folders in the case of an emptying simulation
				
				if [ $emptConfig = "y" ]
				then
					mkdir -p ${parameter}Toprun$i;
					cp $topEmpt ${parameter}Toprun$i/$topEmpt;
					cd ${parameter}Toprun$i
					mv $topEmpt config.in
					cd ..

					mkdir -p ${parameter}Bottomrun$i;
					cp $botEmpt ${parameter}Bottomrun$i/$botEmpt;
					cd ${parameter}Bottomrun$i
					mv $botEmpt config.in
					cd ..
				fi

				## As there are two separate run folders created for emptying simulations - Top and Bottom - must copy to each

				mkdir -p ${parameter}Toprun$i
				mkdir -p ${parameter}Bottomrun$i
				cp $fileToCopy ${parameter}Toprun$i/$fileToCopy
				cp $fileToCopy ${parameter}Bottomrun$i/$fileToCopy
			fi
		done
	fi
done
