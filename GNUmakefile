##**************************************************
## \file GNUmakefile
## \brief: GNUmakefile of G4HidraSim project
## \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
## \start date: 7 July 2021
##**************************************************

name := HidraSim
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

##**************************************************
