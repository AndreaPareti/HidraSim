//**************************************************
// \file HidraSimSignalHelper.cc
// \brief: Implementation of HidraSimSignalHelper
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 1 September 2021
//**************************************************

//Includers from project files
//
#include "HidraSimSignalHelper.hh"

//Includers from Geant4
#include "G4Poisson.hh"

HidraSimSignalHelper* HidraSimSignalHelper::instance = 0;

//Define (private) constructor (singleton)
//
HidraSimSignalHelper::HidraSimSignalHelper(){}

//Define Instance() method
//
HidraSimSignalHelper* HidraSimSignalHelper::Instance(){
    if (instance==0){
        instance = new HidraSimSignalHelper;
    }
    return HidraSimSignalHelper::instance;
}

//Define ApplyBirks() method
//
G4double HidraSimSignalHelper::ApplyBirks( const G4double& de, const G4double& steplength ) {
		
    const G4double k_B = 0.126; //Birks constant
    return (de/steplength) / ( 1+k_B*(de/steplength) ) * steplength;

}

//Define SmearSSignal() method
//
G4int HidraSimSignalHelper::SmearSSignal( const G4double& satde ) {
		
    return G4Poisson(satde*9.5);
		
}

//Define SmearCSignal() method
//
G4int HidraSimSignalHelper::SmearCSignal( ){
		
    return G4Poisson(0.153);

}

//**************************************************
