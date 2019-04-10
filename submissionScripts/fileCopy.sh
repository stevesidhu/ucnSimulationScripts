## Make sure the user has the desired file to copy in this directory

printf "\n\nMake sure you have the file you want to copy in this directory!\n"


## Asks the user how many different geometries - i.e., values of the
## parameter in the study - there are, and asks what the study is called

read -p $'\n\nWhat name have you assigned the study? ' parameter

path=${parameter}Run
pathTop=${parameter}TopRun
pathBot=${parameter}BottomRun

## Asks the user if this is an emptying simulation, if so, if they'd like to copy the emptying config files

read -p $'\n\nIs this an emptying/storage simulation? (y/n) ' emptying

## Determines the number of different geometries from $parameter so that it knows how many times to loop

if [ $emptying = "y" ]
then
	numGeometries=$(ls -d -1q $pathTop* | wc -l)
else
	numGeometries=$(ls -d -1q $path* | wc -l)
fi

if [ $emptying = "y" ]
then
	read -p $'\n\nWould you like to copy the emptying config files? (y/n) ' emptConfig
	
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

	read -p $'\n\nWhat file would you like to copy now? If finished, type done: ' fileToCopy
	
	if [ $fileToCopy != "done" ]
	then
		for i in $(seq 1 $numGeometries)
		do
		
			## Makes a new directory for each Run without overwriting pre-existing
			## directory, and copies user-specified file to it
			
			if [ $emptying != 'y' ]
			then
				mkdir -p $path$i
				cp $fileToCopy $path$i
			else
			
				## Copies the earlier specified config files into the Run folders in the case of an emptying simulation
				
				if [ $emptConfig = "y" ]
				then
					mkdir -p $pathTop$i
					cp $topEmpt $pathTop$i/$topEmpt
					cd $pathTop$i
					mv $topEmpt config.in
					cd ..

					mkdir -p $pathBot$i
					cp $botEmpt $pathBot$i/$botEmpt
					cd $pathBot$i
					mv $botEmpt config.in
					cd ..
				fi

				## As there are two separate Run folders created for emptying simulations - Top and Bottom - must copy to each

				mkdir -p $pathTop$i
				mkdir -p $pathBot$i
				cp $fileToCopy $pathTop$i/$fileToCopy
				cp $fileToCopy $pathBot$i/$fileToCopy
			fi
		done
	fi
done
