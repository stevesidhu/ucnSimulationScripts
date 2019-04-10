## This script briefly tests PENTrack locally in order ensure that the previous steps have been done correctly.

## First, prompt the user for the number of geometries, and the name of the study

read -p $'\n\nWhat have you re-named the study? ' parameter

path=${parameter}Run
pathTop=${parameter}TopRun
pathBot=${parameter}BottomRun

## Counts the number of different values of the parameter being simulated

read -p $'\n\nIs this an emptying/storage simulation? (y/n) ' emptying

if [ $emptying = 'y' ]
then
	numValues=$(ls -d -1q $pathTop* | wc -l)
else
	numValues=$(ls -d -1q $path* | wc -l)
fi

for i in $(seq 1 $numValues);
do
	if [ $emptying = "y" ]
	then
		echo "Checking Run$i Top!"
		until ./PENTrack/PENTrack 99 $pathTop$i . | grep -m 1 "Particle no.: 1"; do : ; done
		echo "Run $i Top is good to go!"
		
		echo "Checking Run$i Bottom!"
		until ./PENTrack/PENTrack 99 $pathBot$i . | grep -m 1 "Particle no.: 1"; do : ; done
		echo "Run $i Bottom is good to go!"
	else
		echo "Checking Run$i!"
		until ./PENTrack/PENTrack 99 $path$i . | grep -m 1 "Particle no.: 1"; do : ; done
		echo "Run $i is good to go!"
	fi
done
