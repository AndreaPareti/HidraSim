#include "HidraSimPrimaryGeneratorMessenger.hh"
#include "HidraSimPrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

HidraSimPrimaryGeneratorMessenger::HidraSimPrimaryGeneratorMessenger(
    HidraSimPrimaryGeneratorAction* gen)
  : G4UImessenger(), fGen(gen)
{
  fDir = new G4UIdirectory("/hidra/gen/");
  fDir->SetGuidance("Primary generator controls.");

  fModeCmd = new G4UIcmdWithAString("/hidra/gen/mode", this);
  fModeCmd->SetGuidance("Select generator mode: gps, gun, or pythia.");
  fModeCmd->SetParameterName("mode", false);
  fModeCmd->SetCandidates("gps gun pythia");
}

HidraSimPrimaryGeneratorMessenger::~HidraSimPrimaryGeneratorMessenger()
{
  delete fModeCmd;
  delete fDir;
}

void HidraSimPrimaryGeneratorMessenger::SetNewValue(
    G4UIcommand* cmd, G4String newValue)
{
  if (cmd == fModeCmd && fGen != nullptr) {
    fGen->SetMode(newValue);
  }
}