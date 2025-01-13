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
#include "HidraSimGeoPar.hh"

//Includers from Geant4
#include "G4Poisson.hh"
#include "G4Tubs.hh"
#include "G4NavigationHistory.hh"



// Includers for interpolation
#include<vector>
#include<algorithm>

G4double linearInterpolate(const std::vector<G4double>& x, const std::vector<G4double>& y, G4double xi) {

    auto it = std::lower_bound(x.begin(), x.end(), xi);

    if (it == x.begin()) {
        return y.front();
    }
    if (it == x.end()) {
        return y.back();
    }

    size_t idx = it - x.begin();
    double x1 = x[idx - 1], x2 = x[idx];
    double y1 = y[idx - 1], y2 = y[idx];

    return y1 + (xi - x1) * (y2 - y1) / (x2 - x1);
}




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
    return G4Poisson(0.153);            // Original
    //return G4Poisson(0.243);            // TB2023

}

// Dummy functions to recover #phe from TB 
G4int HidraSimSignalHelper::SmearSSignalOpticalPhoton( ) {
    return G4Poisson(0.165);        
    //return G4Poisson(satde*21.32);		// TB2023 
}

G4int HidraSimSignalHelper::SmearCSignalOpticalPhoton( ) {
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
	//const G4double SAttenuationLength = 700.0*CLHEP::cm; // from Bedeschi Datasheet
	const G4double SAttenuationLength = 350.0*CLHEP::cm; // from Bedeschi Datasheet
	//const G4double SAttenuationLength = attenuation->Eval(att_length[8])*CLHEP::cm; // from Bedeschi Datasheet


	//const G4double SAttenuationLength = 1.*CLHEP::km; // 

    return AttenuateHelper(signal, distance, SAttenuationLength);    

}

//Define AttenuateCSignal() method
//
G4int HidraSimSignalHelper::AttenuateCSignal(const G4int& signal, const G4double& distance) {
	//const G4double CAttenuationLength = 388.9*CLHEP::cm; // from test beam data
	//const G4double CAttenuationLength = 2000.0*CLHEP::cm; // from test beam data
	//const G4double CAttenuationLength = 700.0*CLHEP::cm; // from test beam data
	const G4double CAttenuationLength = 350.0*CLHEP::cm; // from test beam data

	//const G4double CAttenuationLength = 1.*CLHEP::km; // 

    return AttenuateHelper(signal, distance, CAttenuationLength);    
    
}



//Define ApplyPMTdishomogeneity() method
G4int HidraSimSignalHelper::ApplyPMTdishomogeneity(const G4int& signal, const G4int& SiPMID) {
    double probability_of_survival = 0.7;
    G4int colID = static_cast<G4int>(SiPMID/(NofFibersrow/2));
    G4int rowID = 2*static_cast<G4int>(SiPMID%(NofFibersrow/2));     

    G4int survived_photons = 0;

    //if(SiPMID <10){G4cout << "SiPMID: " << SiPMID << "\tRow: " << rowID << "\tCol: " << colID << G4endl;}
    if(rowID<=2 || rowID >= (NofFibersrow-2) || colID <=2 || colID >= (NofFiberscolumn-2) ){

        for (int i=0; i<signal; i++)
        {
            // Simulate drawing between 0 and 1 with probability x of getting 1
            if (G4UniformRand() <= probability_of_survival) survived_photons++;
        }
    }

    else{survived_photons = signal;}  

    return survived_photons;

}



/*******************************************************************************************************/
// Not used in current simulation



/*
// Considered wavelengths for optical photons
std::vector<G4double> wavelengths = {299.8065764 , 309.84507746, 319.83492391, 329.78723404,
                                    339.81912853, 349.78843441, 359.83749275, 369.81807337,
                                    379.7856049 , 389.81452373, 399.74210187, 409.78189028,
                                    419.76980366, 429.80935875, 439.71631206, 449.76423649,
                                    459.77011494, 469.6969697 , 479.69052224, 489.7314376 ,
                                    499.79846836, 509.65885738, 519.69823973, 529.68816745,
                                    539.59965187, 549.64539007, 559.566787  , 569.59118052,
                                    579.71014493, 589.6338564 , 599.6131528 , 609.63618486};

// Attenuation lenght in mm, as measured in lab (in cm, S fibres) 
std::vector<G4double> att_length = { 0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
                                     0.,    0.,    0.,    0.,    0.,   60.,  260.,  425.,  540.,
                                     630.,  725.,  850.,  950., 1010., 1075., 1100., 1160., 1240.,
                                     1330., 1275.,    0.,    0.,    0.};


// Attenuation lenght in mm, as measured in lab (in cm, C fibres) 
std::vector<G4double> att_length_C = {  868.58896381,   965.09884867,  1085.73620476,  1240.84137687,
                                        1447.64827301,  1737.17792761,  2171.47240952,  2714.3405119 ,
                                        3216.99616225,  3776.47375568,  4571.52086214,  5428.68102379,
                                        6204.20688433,  6681.45356774,  7238.24136505,  7896.26330733,
                                        8272.27584578,  8685.88963807,  8685.88963807,  8272.27584578,
                                        8685.88963807,  9143.04172428,  9650.98848674,  8685.88963807,
                                        7552.94751136,  7896.26330733, 10857.36204758, 12408.41376866,
                                        10857.36204758,  9650.98848674,  7896.26330733,  3474.35585523};


  
std::vector<G4double> Spmt_PDEs = {0.2101, 0.2224, 0.2307, 0.2351, 0.2388, 0.2401, 0.2379, 0.2379,
                                    0.2398, 0.2331, 0.2296, 0.2257, 0.2208, 0.2138, 0.2059, 0.1966,
                                    0.1874, 0.1744, 0.1629, 0.1522, 0.1433, 0.1334, 0.1156, 0.0898,
                                    0.0703, 0.0584, 0.0491, 0.0412, 0.0335, 0.0262, 0.02  , 0.0141};

    
std::vector<G4double> Cpmt_PDEs = {0.3274, 0.3553, 0.3729, 0.3811, 0.3862, 0.3863, 0.3824, 0.381 ,
                                    0.3807, 0.3764, 0.3714, 0.3623, 0.3516, 0.3406, 0.3259, 0.3068,
                                    0.2845, 0.2633, 0.2474, 0.2337, 0.218 , 0.1841, 0.1439, 0.1178,
                                    0.1024, 0.0908, 0.0805, 0.0703, 0.0609, 0.0517, 0.0428, 0.0345};



std::vector<G4double> Ssipm_PDEs = {0.05 , 0.07 , 0.08 , 0.095, 0.11 , 0.12 , 0.12 , 0.125, 0.135,
                                    0.14 , 0.145, 0.16 , 0.168, 0.17 , 0.173, 0.175, 0.178, 0.179,
                                    0.18 , 0.179, 0.178, 0.177, 0.175, 0.17 , 0.165, 0.16 , 0.15 ,
                                    0.15 , 0.145, 0.14 , 0.135, 0.13 };                            

  
std::vector<G4double> Csipm_PDEs = {0.1  , 0.14 , 0.16 , 0.17 , 0.2  , 0.21 , 0.22 , 0.23 , 0.25 ,
                                    0.26 , 0.265, 0.275, 0.29 , 0.3  , 0.31 , 0.315, 0.318, 0.32 ,
                                    0.32 , 0.32 , 0.31 , 0.305, 0.3  , 0.29 , 0.28 , 0.27 , 0.26 ,
                                    0.25 , 0.24 , 0.235, 0.23 , 0.22 };                                       



std::vector<G4double> pmtS_correction = {0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
                                        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
                                        0.       , 0.       , 0.       , 0.008847 , 0.120873 , 0.13952  ,
                                        0.14661  , 0.137741 , 0.130403 , 0.122061 , 0.106352 , 0.0828854,
                                        0.0649572, 0.05402  , 0.045663 , 0.038316 , 0.031155 , 0.024366 ,
                                        0.0186   , 0.013113 };                                


std::vector<G4double> pmtC_correction = {0.3274, 0.3553, 0.3729, 0.3811, 0.3862, 0.3863, 0.3824, 0.381 ,
                                        0.3807, 0.3764, 0.3714, 0.3623, 0.3516, 0.3406, 0.3259, 0.3068,
                                        0.2845, 0.2633, 0.2474, 0.2337, 0.218 , 0.1841, 0.1439, 0.1178,
                                        0.1024, 0.0908, 0.0805, 0.0703, 0.0609, 0.0517, 0.0428, 0.0345};                                


std::vector<G4double> sipmS_correction = {0.      , 0.      , 0.      , 0.      , 0.      , 0.      ,
                                        0.      , 0.      , 0.      , 0.      , 0.      , 0.      ,
                                        0.      , 0.      , 0.      , 0.007875, 0.11481 , 0.1432  ,
                                        0.162   , 0.161995, 0.16198 , 0.161955, 0.161   , 0.15691 ,
                                        0.15246 , 0.148   , 0.1395  , 0.1395  , 0.13485 , 0.1302  ,
                                        0.12555 , 0.1209  };           


                                 
std::vector<G4double> sipmC_correction = {0.1  , 0.14 , 0.16 , 0.17 , 0.2  , 0.21 , 0.22 , 0.23 , 0.25 ,
                                            0.26 , 0.265, 0.275, 0.29 , 0.3  , 0.31 , 0.315, 0.318, 0.32 ,
                                            0.32 , 0.32 , 0.31 , 0.305, 0.3  , 0.29 , 0.28 , 0.27 , 0.26 ,
                                            0.25 , 0.24 , 0.235, 0.23 , 0.22};

*/



/*
//Define AttenuateSSignal() method
//
G4int HidraSimSignalHelper::AttenuateSSignalOverWL(const G4int& signal, const G4double& distance, const G4double& wavelength) {
	//const G4double SAttenuationLength = attenuation_Sfibres->Eval(wavelength)*CLHEP::cm; // from Bedeschi Datasheet
	const G4double SAttenuationLength = 700*CLHEP::cm; // for now, to be changed with more precise measurements
    return AttenuateHelper(signal, distance, SAttenuationLength);    
}

G4int HidraSimSignalHelper::AttenuateCSignalOverWL(const G4int& signal, const G4double& distance, const G4double& wavelength) {
	//const G4double CAttenuationLength = attenuation_Cfibres->Eval(wavelength)*CLHEP::cm; // from ESKA Datasheet
	const G4double CAttenuationLength = 700*CLHEP::cm; // for now, to be changed with more precise measurements
    return AttenuateHelper(signal, distance, CAttenuationLength);    
}



G4double HidraSimSignalHelper::GetSPMTpde(const G4double& wavelength){
    return linearInterpolate(wavelengths, Spmt_PDEs, wavelength);
}

G4double HidraSimSignalHelper::GetCPMTpde(const G4double& wavelength){
    //return pmtCpde->Eval(wavelength);
    return linearInterpolate(wavelengths, Cpmt_PDEs, wavelength);
}

G4double HidraSimSignalHelper::GetSsipmpde(const G4double& wavelength){
    return linearInterpolate(wavelengths, Ssipm_PDEs, wavelength);
}

G4double HidraSimSignalHelper::GetCsipmpde(const G4double& wavelength){
    return linearInterpolate(wavelengths, Csipm_PDEs, wavelength);
}


G4double HidraSimSignalHelper::GetSpmtCorrection(const G4double& wavelength){
    return linearInterpolate(wavelengths, pmtS_correction, wavelength);
}

G4double HidraSimSignalHelper::GetCpmtCorrection(const G4double& wavelength){
    return linearInterpolate(wavelengths, pmtC_correction, wavelength);
}

G4double HidraSimSignalHelper::GetSsipmCorrection(const G4double& wavelength){
    return linearInterpolate(wavelengths, sipmS_correction, wavelength);
}

G4double HidraSimSignalHelper::GetCsipmCorrection(const G4double& wavelength){
    return linearInterpolate(wavelengths, sipmC_correction, wavelength);
}
*/



