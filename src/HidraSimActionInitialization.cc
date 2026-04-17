//**************************************************
// \file HidraSimActionInitialization.cc
// \brief: Implementation of HidraSimActionInitialization class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "HidraSimActionInitialization.hh"
#include "HidraSimDetectorConstruction.hh"
#include "HidraSimEventAction.hh"
#include "HidraSimGeneratorConfig.hh"
#include "HidraSimPrimaryGeneratorAction.hh"
#include "HidraSimRunAction.hh"
#include "HidraSimSteppingAction.hh"

#include "G4GenericMessenger.hh"

HidraSimActionInitialization::HidraSimActionInitialization(
    HidraSimDetectorConstruction* detConstruction)
    : G4VUserActionInitialization(),
      fDetConstruction(detConstruction) {
  fMessenger = new G4GenericMessenger(this, "/hidra/gen/", "Primary generator controls.");

  auto& modeCmd = fMessenger->DeclareProperty(
      "mode",
      HidraSimGeneratorConfig::Mode(),
      "Select generator mode: gps, gun, pythia");

  modeCmd.SetParameterName("mode", false);
  modeCmd.SetCandidates("gps gun pythia");
}

HidraSimActionInitialization::~HidraSimActionInitialization() {
  delete fMessenger;
}

void HidraSimActionInitialization::BuildForMaster() const {
  auto eventAction = new HidraSimEventAction;
  SetUserAction(new HidraSimRunAction(eventAction));
}

void HidraSimActionInitialization::Build() const {
  SetUserAction(new HidraSimPrimaryGeneratorAction);

  auto eventAction = new HidraSimEventAction;
  SetUserAction(new HidraSimRunAction(eventAction));
  SetUserAction(eventAction);
  SetUserAction(new HidraSimSteppingAction(eventAction, fDetConstruction));
}


//**************************************************
