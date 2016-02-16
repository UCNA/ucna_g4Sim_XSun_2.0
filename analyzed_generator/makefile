#####################################
# Makefile template for ROOT programs
#####################################

#didn't use CC or LDFLAGS. Don't understand why not or what's happening
CC		= g++
CXX		= `root-config --cxx`
CXXFLAGS	= `root-config --cflags`
LDFLAGS		= `root-config --ldflags`
LDLIBS		= `root-config --glibs`

CFLAGS 		= $(CXX) $(CXXFLAGS) -W -Wall -o $@ $^ $(LDLIBS)

#These variables need to be changed to the name of the executable/source
SOURCE		= ttree.cc
EXECUTABLE 	= ttree

all:	$(EXECUTABLE)

clean:
	rm $(EXECUTABLE)

$(EXECUTABLE): $(SOURCE)
	$(CFLAGS) -o $(EXECUTABLE)

.cc.o:
	$(CFLAGS) -c $*.cc

