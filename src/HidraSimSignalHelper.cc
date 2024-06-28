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
#include "G4Tubs.hh"
#include "G4NavigationHistory.hh"
#include "TGraph.h"

G4double wavelengths[32] = {609.63618486, 599.6131528,  589.6338564,  579.71014493, 569.59118052,
 559.566787,   549.64539007, 539.59965187, 529.68816745, 519.69823973,
 509.65885738, 499.79846836, 489.7314376,  479.69052224, 469.6969697,
 459.77011494, 449.76423649, 439.71631206, 429.80935875, 419.76980366,
 409.78189028, 399.74210187, 389.81452373, 379.7856049, 369.81807337,
 359.83749275, 349.78843441, 339.81912853, 329.78723404, 319.83492391,
 309.84507746, 299.8065764};


G4double att_length[32] = {0., 0., 0., 1275., 1330., 
    1240., 1160., 1100., 1075., 1010., 
    950., 850., 725., 630., 540., 
    425., 260., 60., 0., 0., 
    0., 0., 0., 0., 0., 
    0., 0., 0., 0., 0., 
    0., 0.};

G4double att_length_C[32] = { 3474.35585523,  7896.26330733,  9650.98848674, 10857.36204758,
 12408.41376866, 10857.36204758,  7896.26330733,  7552.94751136,
  8685.88963807,  9650.98848674,  9143.04172428,  8685.88963807,
  8272.27584578,  8685.88963807,  8685.88963807,  8272.27584578,
  7896.26330733,  7238.24136505,  6681.45356774,  6204.20688433,
  5428.68102379,  4571.52086214,  3776.47375568,  3216.99616225,
  2714.3405119,  2171.47240952,  1737.17792761,  1447.64827301,
  1240.84137687,  1085.73620476,   965.09884867,   868.58896381};


G4double Spmt_PDEs[32] = {0.0141, 0.02,  0.0262, 0.0335, 0.0412, 0.0491, 0.0584, 0.0703, 0.0898, 0.1156,
        0.1334, 0.1433, 0.1522, 0.1629, 0.1744, 0.1874, 0.1966, 0.2059, 0.2138, 0.2208,
        0.2257, 0.2296, 0.2331, 0.2398, 0.2379, 0.2379, 0.2401, 0.2388, 0.2351, 0.2307,
        0.2224, 0.2101};    

G4double Cpmt_PDEs[32] = {0.0345, 0.0428, 0.0517, 0.0609, 0.0703, 0.0805, 0.0908, 0.1024, 0.1178, 0.1439,
    0.1841, 0.218,  0.2337, 0.2474, 0.2633, 0.2845, 0.3068, 0.3259, 0.3406, 0.3516,
    0.3623, 0.3714, 0.3764, 0.3807, 0.381,  0.3824, 0.3863, 0.3862, 0.3811, 0.3729,
    0.3553, 0.3274};    


G4double Ssipm_PDEs[32] = {0.13, 0.135, 0.14, 0.145, 0.15,
		0.15, 0.16, 0.165, 0.17, 0.175,
		0.177,   0.178, 0.179, 0.18, 0.179, 
		0.178, 0.175, 0.173, 0.170, 0.168,
		0.160, 0.145, 0.14, 0.135, 0.125, 
		0.12, 0.12, 0.11, 0.095, 0.08, 
		0.07, 0.05};

G4double Csipm_PDEs[32] = {0.22, 0.23, 0.235, 0.24, 0.25, 
		0.26, 0.27, 0.28, 0.29, 0.3, 
		0.305, 0.31, 0.32, 0.32,  0.32,
		0.318, 0.315, 0.31, 0.3, 0.29, 
		0.275, 0.265, 0.26, 0.25, 0.23, 
		0.22, 0.21, 0.2, 0.17, 0.16, 
		0.14, 0.1};  


TGraph* attenuation_Sfibres = new TGraph(32, wavelengths, att_length);
TGraph* attenuation_Cfibres = new TGraph(32, wavelengths, att_length_C);
TGraph* pmtSpde = new TGraph(32, wavelengths, Spmt_PDEs);
TGraph* pmtCpde = new TGraph(32, wavelengths, Cpmt_PDEs);
TGraph* sipmSpde = new TGraph(32, wavelengths, Ssipm_PDEs);
TGraph* sipmCpde = new TGraph(32, wavelengths, Csipm_PDEs);


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
    return G4Poisson(satde*9.5);        // Original
    //return G4Poisson(satde*21.32);		// TB2023 
}

//Define SmearCSignal() method
//
G4int HidraSimSignalHelper::SmearCSignal( ){
   // return G4Poisson(0.153);            // Original
    return G4Poisson(0.243);            // TB2023

}

G4int HidraSimSignalHelper::SmearSSignal_new( ) {
    return G4Poisson(0.165);        // try to recover #phe from tb
    //return G4Poisson(satde*21.32);		// TB2023 
}

G4int HidraSimSignalHelper::SmearCSignal_new( ) {
    return G4Poisson(0.850);        // try to recover #phe from tb
}


//Define GetDistanceToSiPM() method
//
G4double HidraSimSignalHelper::GetDistanceToSiPM(const G4Step* step) {

    // Get the pre-step point
    const G4StepPoint* preStepPoint = step->GetPreStepPoint();
    // Get the global position of the pre-step point
    G4ThreeVector globalPos = preStepPoint->GetPosition();
    // Get the local position of the pre-step point in the current volume's coordinate system
    G4ThreeVector localPos = preStepPoint->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(globalPos);
    // G4cout << "Local Position (X,Y,Z): (" << localPos.x()/CLHEP::mm << ", " << localPos.y()/CLHEP::mm << ", " << localPos.z()/CLHEP::mm << ") mm" << G4endl;

    // Get the logical volume of the current step
    G4LogicalVolume* currentVolume = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    // Get the solid associated with the logical volume
    G4Tubs* solid = dynamic_cast<G4Tubs*>(currentVolume->GetSolid());
    // Get the dimensions of the solid (size of the volume)
    G4double size = solid->GetZHalfLength();

    G4double distance_to_sipm = size - localPos.z();
    return distance_to_sipm;

}


//Define AttenuateHelper() method
G4int HidraSimSignalHelper::AttenuateHelper(const G4int& signal, const G4double& distance, const G4double& attenuation_length) {
    double probability_of_survival = exp(-distance/attenuation_length);

    G4int survived_photons = 0;
    for (int i=0; i<signal; i++)
    {
        // Simulate drawing between 0 and 1 with probability x of getting 1
        if (G4UniformRand() <= probability_of_survival) survived_photons++;
    }

    return survived_photons;

}

//Define AttenuateSSignal() method
//
G4int HidraSimSignalHelper::AttenuateSSignal(const G4int& signal, const G4double& distance) {
	//const G4double SAttenuationLength = 191.6*CLHEP::cm; // from test beam data
	const G4double SAttenuationLength = 600.0*CLHEP::cm; // from Bedeschi Datasheet
	//const G4double SAttenuationLength = attenuation->Eval(att_length[8])*CLHEP::cm; // from Bedeschi Datasheet


	//const G4double SAttenuationLength = 1.*CLHEP::km; // 

    return AttenuateHelper(signal, distance, SAttenuationLength);    

}

//Define AttenuateCSignal() method
//
G4int HidraSimSignalHelper::AttenuateCSignal(const G4int& signal, const G4double& distance) {
	//const G4double CAttenuationLength = 388.9*CLHEP::cm; // from test beam data
	//const G4double CAttenuationLength = 2000.0*CLHEP::cm; // from test beam data
	const G4double CAttenuationLength = 600.0*CLHEP::cm; // from test beam data

	//const G4double CAttenuationLength = 1.*CLHEP::km; // 

    return AttenuateHelper(signal, distance, CAttenuationLength);    
    
}

//Define AttenuateSSignal() method
//
G4int HidraSimSignalHelper::AttenuateSSignal_WL(const G4int& signal, const G4double& distance, const G4double& wavelength) {
	const G4double SAttenuationLength = attenuation_Sfibres->Eval(wavelength)*CLHEP::cm; // from Bedeschi Datasheet
    return AttenuateHelper(signal, distance, SAttenuationLength);    
}

G4int HidraSimSignalHelper::AttenuateCSignal_WL(const G4int& signal, const G4double& distance, const G4double& wavelength) {
	const G4double CAttenuationLength = attenuation_Cfibres->Eval(wavelength)*CLHEP::cm; // from ESKA Datasheet
    return AttenuateHelper(signal, distance, CAttenuationLength);    
}

G4double HidraSimSignalHelper::GetSPMTpde(const G4double& wavelength){
    return pmtSpde->Eval(wavelength);
}

G4double HidraSimSignalHelper::GetCPMTpde(const G4double& wavelength){
    return pmtCpde->Eval(wavelength);
}

G4double HidraSimSignalHelper::GetSsipmpde(const G4double& wavelength){
    return sipmSpde->Eval(wavelength);
}

G4double HidraSimSignalHelper::GetCsipmpde(const G4double& wavelength){
    return sipmCpde->Eval(wavelength);
}

//**************************************************
