//**************************************************
// \file HidraSimActionInitialization.cc
// \brief: Implementation of HidraSimActionInitialization class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "HidraSimActionInitialization.hh"
#include "HidraSimPrimaryGeneratorAction.hh"
#include "HidraSimRunAction.hh"
#include "HidraSimEventAction.hh"
#include "HidraSimSteppingAction.hh"

//Constructor
//
HidraSimActionInitialization::HidraSimActionInitialization( HidraSimDetectorConstruction* detConstruction, const G4bool FullOptic )
    : G4VUserActionInitialization(),
    fFullOptic( FullOptic ),
    fDetConstruction( detConstruction )		
{}

//De-constructor
//
HidraSimActionInitialization::~HidraSimActionInitialization() {}

//BuildForMaster() method
//
void HidraSimActionInitialization::BuildForMaster() const {
    
    auto eventAction = new HidraSimEventAction;
    SetUserAction( new HidraSimRunAction( eventAction ) );

}

//Build() method
//
void HidraSimActionInitialization::Build() const {
  
    SetUserAction(new HidraSimPrimaryGeneratorAction);
    auto eventAction = new HidraSimEventAction;
    SetUserAction(new HidraSimRunAction( eventAction ));
    SetUserAction(eventAction);
    SetUserAction(new HidraSimSteppingAction(eventAction, fDetConstruction, fFullOptic));

}  

//**************************************************
