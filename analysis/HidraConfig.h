#include<cmath>
#include <array>
#include <string>
#include<TMath.h>
#include<TRandom.h>

/*****************************/
//     Analysis Parameters   //
/*****************************/    

    // Input directory
    //string infolder = "../build/output_electrons_0p70BorderEfficiency_3p5mAttenuation/";
    string infolder = "../build/";

        
    // Brass calibration parameters (80 module geo)
    double SciPheGeV_Brass = 218.73299999999998;   
    double CerPheGeV_Brass = 53.669528571428565;
    double brass_elcontainment = 1.0035;
    double brass_picontainment[9] = {1.08605833, 1.08210717, 1.07979589, 1.07815601, 1.07584473, 1.07420486, 1.07293287, 1.06898171, 1.06375857};
    double chi_brass = 0.3109;

    // Steel calibration parameters (80 module geo)
    //double SciPheGeV_Steel = 255.58200000000002;    //TB2023
    //double CerPheGeV_Steel = 49.32668571428572;     //TB2023 
    //double SciPheGeV_Steel = 510.40757142857143;    // 1Km Attenuation Length
    //double CerPheGeV_Steel = 88.42804285714286;     // 1Km Attenuation Length

    // Original HidraSim values
    //double SciPheGeV_Steel = 225.96657142857143;
    //double CerPheGeV_Steel = 55.18902857142858;
    // With original sim parameters, 7m attenuation length
    double SciPheGeV_Steel = 162.5;
    double CerPheGeV_Steel = 39.9;


    //double CerPheGeV_Steel = 4700;    
    double steel_elcontainment = 1.0006985714285714;
    //double chi_steel = 0.157; //datasheet
    double chi_steel = 0.33;    // 1km
    //double chi_steel = -0.15; 
    double steel_picontainment[9] = {1.08782751, 1.08371817, 1.08131437, 1.07960884, 1.07720503, 1.0754995, 1.07417659, 1., 1.};    



    // PDG numbers for em objects
    std::vector<int> emPIDs {-11, 11, 22, 111, -111};


    /**********************************************/
    /* Geometry Parameters used in the simulation */
    /**********************************************/

    // TB 24 (DRAGO -> 36 modules)
    const int NofmodulesX = 3;
    const int NofmodulesY = 12;
    const int nPMT = NofmodulesX*NofmodulesY;
    // This one is DRAGO (if irot is true)
    const int modflag[36]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                             12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                             24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
    
    const int NofSiPMTowersX = 3;
    const int NofSiPMTowersY=12;
    const int NofModulesSiPM=NofSiPMTowersX*NofSiPMTowersY;
    const int NofActiveModules = 36;
    const int SiPMMod[NofModulesSiPM]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                             12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                             24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
    const bool irot=false;
    //const int NofModulesSiPM = 10;
    //std::vector<int> SiPMtowers = {36, 37, 38, 39, 40, 41, 42, 43, 44, 45};




    /*
    // 80 modules, rotated
    const int NofmodulesX = 20;
    const int NofmodulesY = 5;
    const int NofActiveModules = 80;
    const int nPMT = NofmodulesX*NofmodulesY;
    // HiDRa Geometry
    const int modflag[nPMT]={-1,-1,-1,-1,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-1,-1,-1,-1,-1,
                              10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,
                              30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,
                              50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,
                              -1,-1,-1,-1,-1,70,71,72,73,74,75,76,77,78,79,-1,-1,-1,-1,-1};
    */                          


    
    
 

    // Old 84 Modules
    /*
    const int NofmodulesX = 24;
    const int NofmodulesY = 5;
    const int NofActiveModules = 84;    
    const int nPMT = NofmodulesX*NofmodulesY;    
    const int modflag[nPMT]={-1,-1,-1,-1,-1,-1,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-1,-1,-1,-1,-1,-1,-1,
                       -1,-1,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,-1,-1,
                       30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,
                       -1,-1,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,-1,-1,
                       -1,-1,-1,-1,-1,-1,-1,74,75,76,77,78,79,80,81,82,83,-1,-1,-1,-1,-1,-1,-1};
    */



    // Tower geometry
    const double sq3=1.733;
    const double sq3m1=sq3/3.;
    const double tuberadius = 1.0;
    const int NofFiberscolumn = 64;
    const int NofFibersrow = 16;
    const double moduleZ = 2500.;

    const int NofSciSiPM=(NofFibersrow*NofFiberscolumn*NofModulesSiPM)/2;
    const int NofCerSiPM=(NofFibersrow*NofFiberscolumn*NofModulesSiPM)/2;

    double dtubeY=sq3*tuberadius;                                       // dtubeY given in simulation before rotation
    double dtubeX=2.*tuberadius;                                        // For the reconstruction they must be inverted
    double moduleX = (2*NofFiberscolumn+1)*tuberadius; 
    double moduleY = (NofFibersrow-1.)*dtubeY+4.*tuberadius*sq3m1;




    /*********************/
    /* Utility functions */
    /*********************/

// Scintillator on even rows (start from 0), Cerenkov on uneven rows
// row ID input ranges from 0 to 31 
// need to know fiber type to distinguish global row

void GetSiPMcoordinate(int TowID, int rowID, int colID_original, double &SiPM_X, double &SiPM_Y, std::string fiber, unsigned int grouping)
{
    double TowerOffset = ( (NofModulesSiPM-1)*moduleY )/2 - TowID*moduleY;
    // N.B.: in HidraSim calorimeter if rotated by 90 degrees --> associate sipmX to Y direction of modules and vice versa
    //std::cout << "Col: " << colID << "\tchannel " << colID%grouping << "\tY: " <<  << std::endl; 
    int myCol = colID_original;
    // after grouping, there are nfibercolumns/grouping channels
    unsigned int channel = static_cast<unsigned int>(colID_original/grouping);
    // colID should be such that the corresponding coordinate is in the middle of the channel,
    // between 0 and grouping-1 -> For 8 fibre grouping should be 3.5 (cast it to a double)
    double colID = (static_cast<double>(grouping)-1)/2 + channel*grouping;
    //std::cout << "Pre: " << myCol << "\tChannel: " << channel << "\tPost: " << colID << std::endl;

    if(fiber == "S"){                                                         
        SiPM_Y = +moduleX/2 - tuberadius - (tuberadius*2)*colID -tuberadius;
        SiPM_X = -(-TowerOffset -moduleY/2 + tuberadius + (sq3*tuberadius)*(rowID)+tuberadius*(sq3m1-2) );
        //std::cout << "S X: " << SiPM_X << " Y: " << SiPM_Y << std::endl;
    }    

    if(fiber == "C"){                                                         
        SiPM_Y =  moduleX/2 - tuberadius - tuberadius - (tuberadius*2)*colID -tuberadius;
        SiPM_X = -( -TowerOffset -moduleY/2 + tuberadius + (sq3*tuberadius)*rowID+tuberadius*(2.*sq3m1-1.) );
        //std::cout << "C X: " << SiPM_X << " Y: " << SiPM_Y << std::endl;
    }
}




    // parameters for SiPM grouping and smearing(timing-like)
    // Used for Spatial Resolution studies 
    const unsigned int grouping = 8;
    double smearingZ = 50; // in millimetres
    double RadLenght = 26.5753;   // estimated radiation lenght for HiDRa



   // Mean Z barycenter for electromagnetic showers at 40 GeV
    double meanZbarS_ele = 227.718;  // in mm
    double ZbarVecS_ele[7] = {207.904, 215.358, 222.621, 227.718, 237.043, 243.6, 248.814};
    double meanZbarC_ele = 228.643;
    double ZbarVecC_ele[7] = {209.283, 216.613, 223.671, 228.643, 237.887, 244.488, 249.495};
    // Mean Z barycenter for had showers at 40 GeV
    double meanZbarS_had = 590.164;
    double ZbarVecS_had[7] = {526.466, 545.984, 570.074, 590.164, 612.029, 633.78, 657.165};
    double meanZbarC_had = 575.515;
    double ZbarVecC_had[7] = {561.77, 545.831, 560.253, 575.515, 591.021, 609.396, 630.135};
    // attenuation length - TB2023
    //double att_length_S = 1916.; 
    //double att_length_C = 3889.;
    // attenuation length - 7m
    double att_length_S = 7000.; 
    double att_length_C = 7000.;

    double S_attenuation_correction = (TMath::Exp(-(2500-meanZbarS_had)/att_length_S) ) / (TMath::Exp(-(2500-meanZbarS_ele)/att_length_S) );
    double C_attenuation_correction = (TMath::Exp(-(2500-meanZbarC_had)/att_length_C) ) / (TMath::Exp(-(2500-meanZbarC_ele)/att_length_C) );



// Smearing function to emulate SiPM timing information
void SmearZ(double &recoZ)
{
    double trueZ = recoZ;
    TRandom2 *random = new TRandom2(0); 
    recoZ = random->Gaus(recoZ, smearingZ);
    //std::cout << "TrueZ: " << trueZ << "\tReco: " << recoZ << std::endl;
}


double EnergyContainmentConverter(double energy, double picont, double *array_cont)
{
    if(energy<11. and energy>9.){picont = array_cont[0];}
    if(energy<21. and energy>19.){picont = array_cont[1];}
    if(energy<31. and energy>29.){picont = array_cont[2];}
    if(energy<41. and energy>39.){picont = array_cont[3];}
    if(energy<61. and energy>59.){picont = array_cont[4];}
    if(energy<81. and energy>79.){picont = array_cont[5];}
    if(energy<101. and energy>99.){picont = array_cont[6];}
    if(energy<201. and energy>199.){picont = array_cont[7];}
    if(energy<501. and energy>499.){picont = array_cont[8];}
    return picont;
}



void AngleCorrection(double &X_orig, double &Y_orig, double &X_corr, double &Y_corr, double ang_horiz, double ang_vert, double Z)
{
    double ang_radX = (M_PI*ang_horiz)/180.;
    double ang_radY = (M_PI*ang_vert)/180.;

    X_orig-=8.;
    Y_orig+=8.;
    X_corr = X_orig + Z*TMath::Tan(ang_radX);
    Y_corr = Y_orig - Z*TMath::Tan(ang_radY);

    //X_corr = X_orig*TMath::Cos(ang_radX);
    //std::cout << "Input angle: " << ang_horiz << "\t radiants: " << ang_radX << std::endl;
    //std::cout << "Cosine value: " << TMath::Cos(ang_radX) << std::endl;
    //std::cout << "Sin Value: " << TMath::Sin(ang_radX) << std::endl;
    //std::cout << "Tan Value " << TMath::Tan(ang_radX) << std::endl;

    //std::cout << "Input X: " << X_orig << "\t Corrected X: " << X_corr << std::endl;
    //std::cout << "Correction: " << Z*TMath::Tan(ang_radX) << std::endl;
    //std::cout << "Input Y: " << Y_orig << "\t Corrected Y: " << Y_corr << std::endl;
    //std::cout << "Correction: " << Z*TMath::Tan(ang_radY) << std::endl;

}

