/// \file HidraSimCalorimeterHit.cc
/// \brief Implementation of the HidraSimCalorimeterHit class


#include "HidraSimCalorimeterHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iomanip>

G4ThreadLocal G4Allocator<HidraSimCalorimeterHit>* HidraSimCalorimeterHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HidraSimCalorimeterHit::HidraSimCalorimeterHit()
 : G4VHit(),
   fSiPMID(0),
   fZcoord(0.),
   fPhe(0.),
   fEdep(0.),
   fPheVec(),
   fZVec()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HidraSimCalorimeterHit::~HidraSimCalorimeterHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HidraSimCalorimeterHit::HidraSimCalorimeterHit(const HidraSimCalorimeterHit& fiberHit)
  : G4VHit()
{
  fSiPMID = fiberHit.fSiPMID;
  fTowerID = fiberHit.fTowerID;
  fZcoord = fiberHit.fZcoord;
  fPhe = fiberHit.fPhe;
  fEdep = fiberHit.fEdep;
  fPheVec = fiberHit.fPheVec;
  fZVec = fiberHit.fZVec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const HidraSimCalorimeterHit& HidraSimCalorimeterHit::operator=(const HidraSimCalorimeterHit& right)
{
  fSiPMID = right.fSiPMID;
  fZcoord = right.fZcoord;
  fPhe = right.fPhe;
  fEdep = right.fEdep;
  fTowerID = right.fTowerID;
  fPheVec = right.fPheVec;
  fZVec = right.fZVec;  

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HidraSimCalorimeterHit::operator==(const HidraSimCalorimeterHit& right) const
{
  return ( this == &right ) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HidraSimCalorimeterHit::Print()
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
