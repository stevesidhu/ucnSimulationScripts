There are six scrips within this folder, and they are used to considerably speed up the process of submitting simulation jobs to Cedar, the cluster at Compute Canada which we use to run simulations.


IMPORTANT NOTE: All of these scripts will need to be changed so that they point to the correct paths. Currently they will all point to my home directory on Cedar and Angerona (stewmo, smorawetz), and must be changed to yours if they are to work correctly.

ALSO IMPORTANT NOTE: The naming scheme of the subfolders of the Master folder for a study, e.g. the folders within a study 'KinkHeight', must be named with the value first, e.g. 35cmKinkHeight. This will ensure that the runs are named correctly. For example, a sample structure would be KinkHeight -> 35cmKinkHeight, 45cmKinkHeight, 55cmKinkHeight -> STL files


Firstly, the scripts that are used to prepare files for submission to Cedar.


stlCopy.sh

This script is typically placed in one's home folder on Cedar. This script automates the process of creating separate folders for the STL files created in Solidworks for the simulation. For this script to work properly, it must be located within the same folder as the 'master' folder of a study. For example, if you want to process the 'KinkHeight' study, the folder 'KinkHeight' should be organized as follows (example in brackets)

Master Folder (KinkHeight) -> Parameter Variation Folders (35cmKinkHeight, 45cmKinkHeight, 55cmKinkHeight) -> STL files

This script will ask for three inputs - firstly, the name of the Master Folder (e.g. KinkHeight). It will then ask for you to rename the study to something new. While you could repeat the original name, it is often more convenient to rename it something shorter, e.g. kinkHt. Finally, you are asked if this is an emptying simulation. Type 'y' (no quotes) if it is. This will then create separate run folders for the Top and Bottom emptying simulations. IMPORTANT: The order of the runs will be from smallest to largest value - e.g., if the values of KinkHeight in the study are 35cm, 45cm, and 55cm, run1 will correspond to 35cm, run2 to 45cm, and so on.


fileCopy.sh

This script automates the process of repeatedly moving several files to all the different run folders. The stlCopy.sh script will create a new folder, and the copied files will be added to the STL folder in the following structure.

Run folder (kinkHtrun1) -> STL ->  STL files

This script must be located in the same folder as stlCopy.sh. The script will ask for two specific inputs, and then which files you would like to copy. It firsts asks the name of the study - this is the name you have chosen, NOT the original name (e.g. kinkHt, not KinkHeight). Secondly, it will ask if this is an emptying simulation. If so, it will ask explicity what are the config files for the top and bottom cells. If it is not an emptying simulation, or you when those two files have been specified (with extensions!), you are asked which files you want to copy (specify file extensions!). A copy of these will be placed in each of the Run folders. Type 'done' (no quotes) when finished.


testLocally.sh

This script must also be located in the same folder as stlCopy.sh and fileCopy.sh. This script enables you to test whether PENTrack is able to get started - NOT whether the results are believable, but just to ensure that PENTrack is able to start. This asks for two inputs, firstly, the new name of the study that you have chosen. Secondly, it asks if this is an emptying simulation. It will print out success messages if PENTrack is able to start, and error messages if it is not. If all runs are able to start, then you are ready to submit to Cedar.


submitJobs.sh

This script belongs in the same folder as the previous three. This is a script used to submit jobs to Cedar. This also asks two questions - firstly, what the study has been renamed, and secondly if it is an emptying simulation. This simulation creates two things. Firstly, it creates a batch file for each run folder. This is necessary in order to submit jobs to Cedar. Secondly, it creates a .txt file for each study. This will have a new Job ID on each line, corresponding to the Job ID of one run. This first line is for the first run, etc. This can be easily coped into an Excel spreadsheet to keep track of what simulations correspond to which job IDs.

Once this script has run, you can confirm the jobs have been submitted to cedar by: squeue -u (your username). It should list all the jobs you have queued, among which should be the Job IDs of the jobs you just submitted.


Now, you wait. You can check the progress of your jobs by typing: sacct. Once all tasks have been marked 'COMPLETED', you are ready to move on to processing the results. PENTracks .out output files will be placed in copied to run folders of the same name, this time in the scratch folder. As the scratch folder does not have space limitations, it makes a suitable destination for all of the outputs, in addition to backups of all files. 


processResults.sh


This script is where a lot of the important stuff happens. This will take you from the raw outputs from PENTrack to ROOT files, renaming them, and copying them to our local cluster. This is done in several steps, and the script is quite modular, so you can use to just rename the files, or just to move them, etc.

This script requires the use of a sub-folder within the scratch folder, called toMerge. After confirming that you do in fact want to overwrite any existing simulation (of the SAME name, not any simulation) in the toMerge folder by typing 'y', you will be prompted for three pieces of information. Firstly, the name of study you have chosen, e.g. kinkHt.

IMPORTANT NOTE: Whereas the previous four scripts will prompt you to ask whether the simulation is an emptying simulation, and will automatically create/process the Toprun and Bottomrun folders, you are now required to handle them separately, i.e. run this script once for {study}Top and {study}Bottom.

Secondly, you are asked what are the first and last runs you want to process. This allows you to select a smaller subset of all runs if you wish to do so. Once you have done this, you can begin processing the files. The script will ask you to explicity confirm what things you would like to do, during which you can type 'y' to do so, or anything to skip it. If you have .out files fresh from PENTrack, generally you want to do all of these things. In order, they are:

Move files to the toMerge folder:
This just prompts you to move the run folders into the toMerge folder, to reduce clutter. You MUST do this if you want to use this script to merge the .out files to a root file.

Merge .out files to a .root file:
As the analysis script works with .root files, the raw .out files created by PENTrack must be merged together into one larger .root file for analysis. This questions presents you with three options. Type 'y' to merge the .out files to a .root file, which gives you the option to overwrite individual .root files already created if they exist. Secondly, type 'w' to create new .root files and overwrite all existing .root files, without prompting you. Thirdly, type anything to skip the merging step.

Rename the .root files:
By default, all .root files are given the name {study}.root. In order to differentiate the simulations in which the studied value is different, they must be renamed. Provided you type 'y' to rename them, this part of the script prompts you to manually enter the value of the parameter in each study, in order of the runs. For example, if run 1 corresponded to a KinkHeight value of 35cm, when prompted "Enter the value of the parameter in Run 1: " you would enter 35cm, or 45cm in the case of run 2, etc.

Copying to the UCN cluster:
This part prompts you to copy the files over to our local cluster. As the files are currently stored on Cedar with limited space, we want the .root files created by the simulations to be stored on a drive that can be accessed by multiple people, and with more space. Thus, the .root files are copied to a shared folder on the local cluster. At the time of writing, this folder is @angerona:/ucn/orithyia_data/ssidhu2/cryoSimulations. However, this may need to be changed if that drive runs out of space.

IMPORT NOTE: You will want to change the values in here from the login of the previous person (currently smorawetz) to whatever you username is.

Viewing the neutronend XZ plots:
This final step is used to inspect the profile of the locations where UCN 'die'. This is useful because it will enable you to confirm that there are no holes in the geometry, that the coordinate system was properly assigned when exporting from Solidworks, etc. This step would be much better if it were done immediatly after merging the .out files, however there is a bug (feature?) where Cedar is very inconsistent about whether ROOT plots can be viewed on that drive. To make it simpler, this is just done on our local cluster where there are no problems.

Once this is done, you will now have the .root files ready for analysis, and confirmed whether the simulation was successful. Now all that is left is cleanup.


checkPlots.sh


This is a brief script called by processResults.sh which is run on Angerona in order to check the .root plots. Having this as a separate script seems to enable viewing .ROOT files on an ssh connection, whereas putting the contents of this script in processResults.sh does not. So, although it seems like it could be merged with processResults.sh, it is best kept separate.


deleteOutputs.sh


Finally, this script will delete all of the intermediate requirements and products used in this process - the temporary batch files, the text files listing the job IDs, the .out files, etc. This will prompt you for the name you have assigned the study, e.g. kinkHt, and the first and last run you wish to delete. Finally, it will ask you to confirm that you do in fact want to delete this stuff. Hint: don't do so until you are sure that either the results are reliable, or you are certain that the error was done in Solidworks and thus the .out files are inherently incorrect.





## Feel free to add to this if any additional scripts are created or the existing ones are modified. -SM
