//**************************************************
// \file HidraSimPrimaryGeneratorAction.cc
// \brief: Implementation of HidraSimPrimaryGeneratorAction class
//**************************************************

// Project
#include "HidraSimPrimaryGeneratorAction.hh"
#include "HidraSimGeneratorConfig.hh"

#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"

#include "Pythia8/Analysis.h"
#include "Pythia8/Pythia.h"
#include <memory>

#include <cmath>

HidraSimPrimaryGeneratorAction::HidraSimPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction() {
  fGeneralParticleSource = new G4GeneralParticleSource();
  fParticleGun = new G4ParticleGun(1);

  auto* particleDefinition =
      G4ParticleTable::GetParticleTable()->FindParticle("e-");

  fGeneralParticleSource->SetParticleDefinition(particleDefinition);
  fGeneralParticleSource->SetParticlePosition(G4ThreeVector(0., 0., 0.));

  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fParticleGun->SetParticleEnergy(10. * GeV);
}

HidraSimPrimaryGeneratorAction::~HidraSimPrimaryGeneratorAction() {
  delete fGeneralParticleSource;
  delete fParticleGun;
  delete fMessenger;

}

void HidraSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  const G4String& mode = HidraSimGeneratorConfig::Mode();

  if (mode == "gps" || mode == "GPS") {
    fGeneralParticleSource->GeneratePrimaryVertex(anEvent);
    return;
  }

  if (mode == "gun" || mode == "GUN") {
    fParticleGun->GeneratePrimaryVertex(anEvent);
    return;
  }

  if (mode == "pythia" || mode == "PYTHIA") {
    GeneratePythiaPrimaries(anEvent);
    return;
  }

  G4String msg = "Unknown generator mode: '" + mode +
                 "'. Allowed: gps, gun, pythia.";
  G4Exception("HidraSimPrimaryGeneratorAction::GeneratePrimaries",
              "HIDRA_BAD_MODE",
              FatalException,
              msg);
}

void HidraSimPrimaryGeneratorAction::InitializePythia() {
  if (fPythiaInitialized) {
    return;
  }

  fPythia = std::make_unique<Pythia8::Pythia>();

  fPythia->readString("Beams:idA = -11");
  fPythia->readString("Beams:idB = 11");
  fPythia->readString("Beams:eCM = 91.1876");
  fPythia->readString("PDF:lepton = off");

  fPythia->readString("WeakSingleBoson:ffbar2gmZ = on");
  fPythia->readString("WeakZ0:gmZmode = 2");
  fPythia->readString("23:onMode = off");
  //fPythia->readString("23:onIfAny = 1 2 3 4 5");
  fPythia->readString("23:onIfAny = 1 2"); // restrict to u and d quarks only


  fPythia->readString("Next:numberShowInfo = 0");
  fPythia->readString("Next:numberShowProcess = 0");
  fPythia->readString("Next:numberShowEvent = 0");

  if (!fPythia->init()) {
    G4Exception("HidraSimPrimaryGeneratorAction::InitializePythia",
                "HIDRA_PYTHIA_INIT_FAILED",
                FatalException,
                "Pythia8 failed to initialize for e+e- -> Z -> qq.");
  }

  fPythiaInitialized = true;
}

G4ThreeVector HidraSimPrimaryGeneratorAction::RotateFromAxisToZ(
    const G4ThreeVector& momentum,
    const G4ThreeVector& axis) const {
  G4ThreeVector from = axis.unit();
  const G4ThreeVector to(0., 0., 1.);

  const G4double cosTheta = from.dot(to);

  if (cosTheta > 1. - 1.e-12) {
    return momentum;
  }

  G4ThreeVector rotated = momentum;

  if (cosTheta < -1. + 1.e-12) {
    rotated.rotate(CLHEP::pi, G4ThreeVector(0., 1., 0.));
    return rotated;
  }

  G4ThreeVector rotAxis = from.cross(to).unit();
  const G4double angle = std::acos(cosTheta);
  rotated.rotate(angle, rotAxis);
  return rotated;
}

int HidraSimPrimaryGeneratorAction::FillSelectedHemisphere(
    const Pythia8::Event& event,
    G4PrimaryVertex* vertex) const {
  Pythia8::Thrust thrust;
  if (!thrust.analyze(event)) {
    return 0;
  }

  auto thrustAxis4 = thrust.eventAxis(1);
  G4ThreeVector thrustAxis(thrustAxis4.px(), thrustAxis4.py(), thrustAxis4.pz());

  if (thrustAxis.mag2() <= 0.) {
    return 0;
  }

  if (thrustAxis.z() < 0.) {
    thrustAxis = -thrustAxis;
  }

  auto* particleTable = G4ParticleTable::GetParticleTable();
  int nAdded = 0;

  for (int i = 0; i < event.size(); ++i) {
    const auto& p = event[i];

    if (!p.isFinal()) {
      continue;
    }

    G4ThreeVector p3(p.px() * GeV, p.py() * GeV, p.pz() * GeV);

    if (p3.mag2() <= 0.) {
      continue;
    }

    if (p3.dot(thrustAxis) <= 0.) {
      continue;
    }

    G4ParticleDefinition* definition = particleTable->FindParticle(p.id());
    if (!definition) {
      continue;
    }

    G4ThreeVector rotatedMomentum = RotateFromAxisToZ(p3, thrustAxis);

    auto* primary = new G4PrimaryParticle(
        definition,
        rotatedMomentum.x(),
        rotatedMomentum.y(),
        rotatedMomentum.z());

    vertex->SetPrimary(primary);
    ++nAdded;
  }

  return nAdded;
}

void HidraSimPrimaryGeneratorAction::StoreFinalStateParticles(
    const Pythia8::Event& event) {
  // Clear previous event data
  HidraSimGeneratorConfig::ClearFinalStateData();

  // Count and store all final state particles
  int numFinalState = 0;
  for (int i = 0; i < event.size(); ++i) {
    const auto& p = event[i];

    if (!p.isFinal()) {
      continue;
    }

    // Store particle information
    HidraSimGeneratorConfig::FinalStatePDGID().push_back(p.id());
    HidraSimGeneratorConfig::FinalStateEnergy().push_back(p.e() * GeV);
    HidraSimGeneratorConfig::FinalStatePx().push_back(p.px() * GeV);
    HidraSimGeneratorConfig::FinalStatePy().push_back(p.py() * GeV);
    HidraSimGeneratorConfig::FinalStatePz().push_back(p.pz() * GeV);

    ++numFinalState;
  }

  HidraSimGeneratorConfig::NumFinalStateParticles() = numFinalState;
}

void HidraSimPrimaryGeneratorAction::GeneratePythiaPrimaries(G4Event* anEvent) {
  InitializePythia();

  constexpr int kMaxAttempts = 1000;

  for (int attempt = 0; attempt < kMaxAttempts; ++attempt) {
    if (!fPythia->next()) {
      continue;
    }

    // Store all final state particles from the Pythia event
    StoreFinalStateParticles(fPythia->event);

    //auto* vertex = new G4PrimaryVertex(G4ThreeVector(0., 0., 0.), 0.);
    //auto* vertex = new G4PrimaryVertex(G4ThreeVector(0., 0., -3.125), 0.); // move to negative z, since particles are shot towards +z
    auto* vertex = new G4PrimaryVertex(G4ThreeVector(0., 0., -3.125*m), 0.); // move to negative z, since particles are shot towards +z

    const int nAdded = FillSelectedHemisphere(fPythia->event, vertex);
    if (nAdded > 0) {
      anEvent->AddPrimaryVertex(vertex);
      return;
    }

    delete vertex;
  }

  G4Exception("HidraSimPrimaryGeneratorAction::GeneratePythiaPrimaries",
              "HIDRA_PYTHIA_NO_ACCEPTED_EVENT",
              FatalException,
              "Failed to generate a non-empty +z hemisphere after 1000 Pythia attempts.");
}
