# ucnSimulationScripts

There are six scrips within the submissionScripts folder, and they are used to considerably speed up the process of submitting simulation jobs to Cedar, the cluster at Compute Canada which we use to run simulations.

**IMPORTANT NOTE 1**: All of these scripts will need to be changed so that they point to the correct paths. Currently they will all point to my home directory on Cedar and Angerona (stewmo, smorawetz), and must be changed to yours if they are to work correctly.

**IMPORTANT NOTE 2**: The naming scheme of the subfolders of the master folder for a study - e.g. the folders within a study 'KinkHeight' - MUST be named with the value first, e.g. 35cmKinkHeight. This will ensure that the runs are named correctly. For example, a sample structure would be KinkHeight -> 35cmKinkHeight, 45cmKinkHeight, 55cmKinkHeight -> STL files.


Firstly, the scripts that are used to prepare files for submission to Compute Canada.

## Home folder scripts

### stlCopy.sh

This script is typically placed in one's home folder on Cedar. This script automates the process of creating separate folders for the STL files created in Solidworks for the simulation. For this script to work properly, it must be located within the same folder as the 'master' folder of a study. For example, if you want to process the 'KinkHeight' study, the folder 'KinkHeight' should be organized as follows (example in brackets)

Master Folder (KinkHeight) -> Parameter Variation Folders (35cmKinkHeight, 45cmKinkHeight, 55cmKinkHeight) -> STL files

This script will ask for three inputs - firstly, the name of the Master Folder (e.g. KinkHeight). It will then ask for you to rename the study to something new. While you could repeat the original name, it is often more convenient to rename it something shorter, e.g. kinkHt. Finally, you are asked if this is an emptying or storage simulation. Type 'y' (no quotes) if it is. This will then create separate run folders for the Top and Bottom emptying or storage simulations. IMPORTANT: The order of the runs will be from smallest to largest value - e.g., if the values of KinkHeight in the study are 35cm, 45cm, and 55cm, run1 will correspond to 35cm, run2 to 45cm, and so on.

### fileCopy.sh

This script automates the process of repeatedly moving several files to all the different run folders. The stlCopy.sh script will create a new folder, and the copied files will be added to the STL folder in the following structure.

Run folder (kinkHtrun1) -> STL ->  STL files

This script must be located in the same folder as stlCopy.sh. The script will ask for two specific inputs, and then which files you would like to copy. It firsts asks the name of the study - this is the name you have chosen, NOT the original name (e.g. kinkHt, not KinkHeight). Secondly, it will ask if this is an emptying or storage simulation. If so, it will ask explicity what are the config files for the top and bottom cells. If it is not an emptying or storage simulation, or you when those two files have been specified (with extensions!), you are asked which files you want to copy (specify file extensions!). Additionally, when copying any file which begins with "config", it will automatically be renamed to config.in when it is copied to the run folders. THis enables you to keep multiple config files for different geometries, without needing to manually rename them each time. A copy of these will be placed in each of the Run folders. Type 'done' (no quotes) when finished. 


### testLocally.sh

This script must also be located in the same folder as stlCopy.sh and fileCopy.sh. This script enables you to test whether PENTrack is able to get started - NOT whether the results are believable, but just to ensure that PENTrack is able to start. This asks for two inputs, firstly, the new name of the study that you have chosen. Secondly, it asks if this is an emptying or storage simulation. It will print out success messages if PENTrack is able to start, and error messages if it is not. If all runs are able to start, then you are ready to submit to Cedar.


### submitJobs.sh

This script belongs in the same folder as the previous three. This is a script used to submit jobs to Cedar. This also asks two questions - firstly, what the study has been renamed, and secondly if it is an emptying or storage simulation. This simulation creates two things. Firstly, it creates a batch file for each run folder. This is necessary in order to submit jobs to Cedar. Secondly, it creates a .txt file for each study. This will have a new Job ID on each line, corresponding to the Job ID of one run. This first line is for the first run, etc. This can be easily coped into an Excel spreadsheet to keep track of what simulations correspond to which job IDs.

Once this script has run, you can confirm the jobs have been submitted to cedar by: squeue -u (your username). It should list all the jobs you have queued, among which should be the Job IDs of the jobs you just submitted.

Now, you wait. You can check the progress of your jobs by typing: sacct. Once all tasks have been marked 'COMPLETED', you are ready to move on to processing the results. PENTracks .out output files will be placed in copied to run folders of the same name, this time in the scratch folder. As the scratch folder does not have space limitations, it makes a suitable destination for all of the outputs, in addition to backups of all files. 

## Scratch folder scripts


### processResults.sh


This script is where a lot of the important stuff happens. This will take you from the raw outputs from PENTrack to ROOT files, renaming them, and copying them to our local cluster. This is done in several steps, and the script is quite modular, so you can use to just rename the files, or just to move them, etc.

First, this script will prompt you to enter the name you have assigned the study. Once this is done, it will ask if you would like to move the the files to the toMerge folder. If you have not done this yet, you must do so in order for the subsequent parts of the script to work. If you have not done this yet, is okay to select 'n'. Be warned, if you already have folders in the toMerge folder of the same name, they WILL be overwritten. Once this is done, you are given the opporunity to merge the output files into a root file, which will be done for all runs and cells consecutively. This can take a while. 

Once that is done, you are given the opportunity to rename the root files. You are prompted to enter the relevant parameter for each run. The files withh be renamed as follows: (study name)(input parameter).root. If, for example, you assigned the name kinkHt to the study, and entered the value of the parameter as 45cm in Run 2, the root file for that run will be named kinkHt45cm.root. After the files have been renamed, you are given the opportunity to copy them over to Angerona. This can be done all at once.

**NOTE:** For copying the files to the local cluster, and running the plot-checking script, you will require your own login to access that cluster. Additionally, the location where the root files are stored may change over time, and so it may be necessary to alter the path.

Once this is done, the final step is to quickly check the x-z plots of where neutrons die. This is done automatically, provided you have specified the correct path in this file to where the .root files are stored on the UCN cluster. When examining these, ensure that the outline looks similar to the outline of the geometry in Solidworks in the x-z plane. Possible problems include an incorrect coordinate system, and holes in the geometry.


### deleteOutputs.sh

Finally, this brief script automates the process of deleting all of the files produced along the way. This only requires you to specify the name of the study you assigned, whether it is an emptying or storage simulation, and a confirmation. This will delete all of the temporary batch files produced, the text file which stores the job IDs given by Compute Canada, and the .out files produced by PENTrack. It will not delete ROOT files in the rootFileBackups folder, however.

**NOTE:** Be sure only to do this once you are certain that the results are correct and the .out files have been merged correctly. This author has wasted far too much time by accidentally deleting .out files when they are still needed.
