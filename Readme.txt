Xuan Sun
October 2015

My attempt to reconfigure the B1Example example from the standard GEANT4 10.0 distribution into Michael Mendenhall's simulation.

This attempts to massage the code from B1 to the UCNA simulation written by Mendenhall.

In addition, my simulation will be designed to compile with cmake.

It attempts to separate the GEANT4 simulation code from the separate ROOT/Python analysis.

As such, a lot of Mendenhall's code will be compressed into more complicated, singular classes.

This is not to improve the robustness of his simulation (indeed, it makes it far worse) but rather to make the code easier to understand.

The project itself has too many classes for me to comfortably understand so I'm compiling them to larger files.

In essence, the UCNA simulation will attempt to move all classes into the following "categories":

	1. Total geometry (including EM field)
	2. Physics list
	3. Primary generator action
	4. End of run/event data processing and management.
	5. Track and step level cuts and data recording.


Update: November 30th, 2015.

Completed the basic functions of the simulation. Geometry and magnetic fields reproduced.
Physics list matches Michael's and the particle generation is my own.
Information I/O is done via simple C++ I/O, not compiling with ROOT.


Update: December 22, 2015.

Completed the addition of Particle Generation using M. Mendenhall's tools.
Stand alone code is used to generate event initial conditions, then piped to .txt files.
Simple C++ I/O is used except for conversion between .root files to .txt files.
