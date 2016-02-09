How to use MC_EventGen to generator neutron beta decay events.

Author: Xuan Sun
December 10, 2015


Step 0: if you need to compile, edit the name of the makefiles. Change the current makefile name to something else.
		And edit 'makefile_MC_EventGen' to just 'makefile'. Then enter 'make' in the command line.
		This makefile should work with MC_EventGen.cc.

Step 1: After compiling, run './MC_EventGen'.

Step 2: At the first prompt, enter 'run'. 
		I don't know why but it seems to say so in the GeantSimManager.py file.

Step 3: At the next prompt, the generator prompt, enter 'n1'. 
		The program looks in the directory set in the code.
		In this particular case, it starts in the current directory and looks in NuclearDecaysGenerators/ for the n1.txt file.
		Changing the generator file as needed for whichever type you want to generate events for.

Step 4: At the output prompt, enter where you'd like the output files to go.

Step 5: At the vertex prompt, enter 'f'. 
		This sets the position of the generated particle randomly inside the fiducial volume.
		There are obviously other options - read the code to find out which ones.

Step 6: At the random source prompt, enter 'n'.
		This chooses the Niederretier method of producing quasi-random numbers for all random number generation.

Step 7: Events per TTree [10000]. Enter '10000'.
		This puts 10000 events per tree/file in the print out. I'm not sure if it has to be 10000.
		But to be safe just make files of 10000.

Step 8: N. TTrees [100]. Enter any number N.
		This produces N files of 10000 events titled by Evts_N.root.
