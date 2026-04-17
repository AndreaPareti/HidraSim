//**************************************************
// \file HidraSimPrimaryGeneratorAction.hh
// \brief: Definition of HidraSimPrimaryGeneratorAction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including headers multiple times
//
// include/HidraSimPrimaryGeneratorAction.hh

#ifndef HidraSimPrimaryGeneratorAction_h
#define HidraSimPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <memory>

class G4Event;
class G4PrimaryVertex;
class G4GeneralParticleSource;
class G4ParticleGun;

class HidraSimPrimaryGeneratorMessenger;

namespace Pythia8 {
  class Pythia;
  class Event;
}

class HidraSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  HidraSimPrimaryGeneratorAction();
  ~HidraSimPrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

  void SetMode(const G4String& mode) { fMode = mode; }
  const G4String& GetMode() const { return fMode; }

private:
  void InitializePythia();
  void GeneratePythiaPrimaries(G4Event* event);
  void StoreFinalStateParticles(const Pythia8::Event& event);

  G4ThreeVector RotateFromAxisToZ(const G4ThreeVector& momentum,
                                  const G4ThreeVector& axis) const;

  int FillSelectedHemisphere(const Pythia8::Event& event,
                             G4PrimaryVertex* vertex) const;

private:
  G4GeneralParticleSource* fGeneralParticleSource = nullptr;
  G4ParticleGun*           fParticleGun = nullptr;

  G4String fMode = "gps";

  HidraSimPrimaryGeneratorMessenger* fMessenger = nullptr;

  std::unique_ptr<Pythia8::Pythia> fPythia;
  G4bool fPythiaInitialized = false;
};

#endif