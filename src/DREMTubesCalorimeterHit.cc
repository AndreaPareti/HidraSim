/// \file DREMTubesCalorimeterHit.cc
/// \brief Implementation of the DREMTubesCalorimeterHit class


#include "DREMTubesCalorimeterHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iomanip>

G4ThreadLocal G4Allocator<DREMTubesCalorimeterHit>* DREMTubesCalorimeterHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DREMTubesCalorimeterHit::DREMTubesCalorimeterHit()
 : G4VHit(),
   fSiPMID(0),
   fZcoord(0.),
   fPhe(0.),
   fEdep(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DREMTubesCalorimeterHit::~DREMTubesCalorimeterHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DREMTubesCalorimeterHit::DREMTubesCalorimeterHit(const DREMTubesCalorimeterHit& fiberHit)
  : G4VHit()
{
  fSiPMID = fiberHit.fSiPMID;
  fTowerID = fiberHit.fTowerID;
  fZcoord = fiberHit.fZcoord;
  fPhe = fiberHit.fPhe;
  fEdep = fiberHit.fEdep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const DREMTubesCalorimeterHit& DREMTubesCalorimeterHit::operator=(const DREMTubesCalorimeterHit& right)
{
  fSiPMID = right.fSiPMID;
  fZcoord = right.fZcoord;
  fPhe = right.fPhe;
  fEdep = right.fEdep;
  fTowerID = right.fTowerID;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DREMTubesCalorimeterHit::operator==(const DREMTubesCalorimeterHit& right) const
{
  return ( this == &right ) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DREMTubesCalorimeterHit::Print()
{
  if(fPhe>0.)
  {
  G4cout
     << "Fiber number: " << fSiPMID 
     << "Tower: " << fTowerID
     << "\tDeposited photoelectrons: " 
     << std::setw(7) << fPhe
     << "\tZ coordinate: " 
     << std::setw(7) << G4BestUnit( fZcoord,"Length")
     //<< " SiPM ID: "
     //<<  fSiPMID
     << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
