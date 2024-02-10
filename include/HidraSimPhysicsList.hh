//**************************************************
// \file HidraSimPhysicsList.hh
// \brief: Definition of HidraSimPhysicsList class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef HidraSimPhysicsList_h
#define HidraSimPhysicsList_h 1

//Includers from Geant4
//
#include "G4VModularPhysicsList.hh"

//Includers from project files
//
#include "HidraSimOpticalPhysics.hh"

class HidraSimPhysicsList : public G4VModularPhysicsList{
    
    public:
        //Constructor
        //
        HidraSimPhysicsList(G4String, const G4bool FullOptic );
        //De-constructor
        //
        virtual ~HidraSimPhysicsList();
    
        HidraSimOpticalPhysics* OpPhysics;
    
        G4bool AbsorptionOn;
    
    private:
        
        G4bool fFullOptic;   

};

#endif

//**************************************************


