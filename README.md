# Bacterial Genomics

This repository includes example scripts that model possible steps in a genomics workflow, 
designed to be implemented on virtual machines available through 
[XSEDE's Jetstream environment](https://portal.xsede.org/jetstream).

To access your cloud instance, you will use the following software:

* Windows: GitBash
* Mac: Terminal 

Log on to your instance using secure shell:

* `ssh USERNAME@IP.ADDRESS`
* `USERNAME` is `hertweck` for all cloud instances in this class
* `IP.ADDRESS` has been included in an individual email to you personally
* The first time you log in, you will be asked 
`Are you sure you want to continue connecting (yes/no)?`; 
you should type "yes" then hit enter.
* You will be prompted for your password, 
which was also included in an email from your instructor. 
Remember that you will not see anything appear on the screen as you type your password 
(this is ok!).

Running project scripts:

* Obtain the scripts from GitHub using: 
`git clone https://github.com/k8hertweck/bacterialGenomics.git`
* You should now see a `bacterialGenomics` directory
* Change your working directory to `bacterialGenomics` prior to executing all scripts, 
following usage instructions described within each script's comments, and in the following order:
	* `projectSetup.sh`
	* `bacterialAssembly.sh`
	* `bacterialReadMapping.sh`
* Please look through the comments in each script carefully to discern the processes 
being implemented. You should also consult the project documentation and rubric on Canvas 
for more information on tasks to be completed.

Transferring data from Jetstream to your local computer:
* Start in a shell on your local computer
* Use secure copy to transfer files to home directory of cloud instance: 
`scp USERNAME@IP.ADDRESS:PATH/TO/FILENAME .`
* `USERNAME` and `IP.ADDRESS` are the same as used for login
* `PATH/TO/FILENAME` is the location and filename for what you would like to download
* If you don't know where the file is located, go to the directory in the cloud instance 
and type `pwd`
* You will be required to enter your password each time you execute an `scp` command
* If you are attempting to download an entire directory, it is recommended to 
archive and compress the directory as follows: 
`tar -cvzf DIRECTORY.tar.gz DIRECTORY/`
* `ARCHIVE` is the name of the directory you want to compress
	
Other notes:
* To log out of your cloud instance, type `exit` or simply close the shell window.
* You can have multiple bash shells open at once; it is common to have one open for 
your cloud instance and another to access your local computer
