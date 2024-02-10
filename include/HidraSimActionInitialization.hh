//**************************************************
// \file HidraSimActionInitialization.hh
// \brief: Definition of HidraSimActionInitialization class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef HidraSimActionInitialization_h
#define HidraSimActionInitialization_h 1

//Includers from Geant4
//
#include "G4VUserActionInitialization.hh"
#include "G4Types.hh"

//Includers from C++
//
#include <chrono>
#include <random>

//Forward declaration
//
class HidraSimDetectorConstruction;

class HidraSimActionInitialization : public G4VUserActionInitialization {
    
    public:
        //Constructor
        //
        HidraSimActionInitialization(HidraSimDetectorConstruction*, const G4bool FullOptic );
        virtual ~HidraSimActionInitialization();

        virtual void BuildForMaster() const;
        virtual void Build() const;

    private:

        G4bool fFullOptic;

	HidraSimDetectorConstruction* fDetConstruction;

};

#endif

//**************************************************
