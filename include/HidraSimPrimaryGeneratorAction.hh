//**************************************************
// \file HidraSimPrimaryGeneratorAction.hh
// \brief: Definition of HidraSimPrimaryGeneratorAction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including headers multiple times
//
#ifndef HidraSimPrimaryGeneratorAction_h
#define HidraSimPrimaryGeneratorAction_h 1

//Includers from Geant4
//
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
//class G4ParticleGun;         //you can switch to G4ParticleGun
class G4Event;

class HidraSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    
    public:
        //Constructor()
        //
        HidraSimPrimaryGeneratorAction();    
        //De-constructor()
        //
        virtual ~HidraSimPrimaryGeneratorAction();

        virtual void GeneratePrimaries(G4Event* event);
  
    private:
        G4GeneralParticleSource* fGeneralParticleSource;
        //G4ParticleGun*  fParticleGun; // G4ParticleGun

};

#endif

//**************************************************
