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
#include "G4Step.hh"

class HidraSimSignalHelper {

    private:

        static HidraSimSignalHelper* instance;

	const G4double fk_B = 0.126; //Birks constant

	const G4double fSAttenuationLength = 191.6*CLHEP::cm; // from test beam data
	const G4double fCAttenuationLength = 388.9*CLHEP::cm; // from test beam data

	//Private constructor (singleton)
        //
	HidraSimSignalHelper();

    public:

    	static HidraSimSignalHelper* Instance();

    	G4double ApplyBirks( const G4double& de, const G4double& steplength );

	G4int SmearSSignal( const G4double& de );

    	G4int SmearCSignal( );

    	G4double GetDistanceToSiPM(const G4Step* step);

	G4int AttenuateHelper(const G4int& signal, const G4double& distance, const G4double& attenuation_length);

	G4int AttenuateSSignal(const G4int& signal, const G4double& distance);

	G4int AttenuateCSignal(const G4int& signal, const G4double& distance);

};

#endif

//**************************************************