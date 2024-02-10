//**************************************************
// \file HidraSimSignalHelper.hh
// \brief: Definition of HidraSimSignalHelper class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 1 September 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef HidraSimSignalHelper_h
#define HidraSimSignalHelper_h

//Includers from Geant4
//
#include "globals.hh"

class HidraSimSignalHelper {

    private:

        static HidraSimSignalHelper* instance;

	//Private constructor (singleton)
        //
	HidraSimSignalHelper();

    public:

    	static HidraSimSignalHelper* Instance();

    	G4double ApplyBirks( const G4double& de, const G4double& steplength );

	G4int SmearSSignal( const G4double& de );

    	G4int SmearCSignal( );

};

#endif

//**************************************************
