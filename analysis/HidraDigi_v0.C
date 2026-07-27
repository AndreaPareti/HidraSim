// ROOT macro that digitises the binned calorimeter hits written by HidraSim.
//
// Example (attenuation is disabled by default):
//   root -l -b -q 'HidraDigi_v0.C("../DREMTubesout_Run0.root")'
//
// Enable the two attenuation models independently (lengths are in mm):
//   root -l -b -q 'HidraDigi_v0.C("input.root","digi.root",true,6500.,true,9000.)'

#include <TFile.h>
#include <TParameter.h>
#include <TRandom3.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <tuple>
#include <vector>

#include "HidraGeo.h"

namespace {

double DistanceFromBin(int bin, double start, double width, int nBins)
{
  // HidraSim uses bin 0 and nBins+1 for underflow and overflow.
  if (bin <= 0) return std::max(0.0, start);
  if (bin > nBins) return std::max(0.0, start + nBins * width);
  return std::max(0.0, start + (bin - 0.5) * width);
}

int DetectPhotoelectrons(TRandom3 &random, double mean, bool attenuate,
                         double distance_mm, double attenuationLength_mm)
{
  if (!(mean > 0.0)) return 0;

  // This is the same ordering as DREMTubesSignalHelper: Poisson conversion,
  // followed by independent survival of every p.e. when attenuation is used.
  const double poissonDraw = random.PoissonD(mean);
  const auto maximum = static_cast<double>(std::numeric_limits<int>::max());
  const int photoelectrons = static_cast<int>(std::min(poissonDraw, maximum));
  if (!attenuate) return photoelectrons;

  const double survival = std::exp(-distance_mm / attenuationLength_mm);
  return random.Binomial(photoelectrons, std::clamp(survival, 0.0, 1.0));
}

bool HasBranch(TTree *tree, const char *name)
{
  if (tree->GetBranch(name)) return true;
  std::cerr << "HidraDigi_v0: required branch '" << name << "' is missing.\n";
  return false;
}

double TimeFromBin(int bin, double start_ns, double width_ns, int nBins)
{
  // Represent underflow and overflow at the centres of the adjacent bins.
  if (bin <= 0) return start_ns - 0.5 * width_ns;
  if (bin > nBins) return start_ns + (nBins + 0.5) * width_ns;
  return start_ns + (bin - 0.5) * width_ns;
}

bool IsSiPMTower(int towerID)
{
  return std::find(std::begin(SiPMMod), std::end(SiPMMod), towerID) !=
         std::end(SiPMMod);
}

bool BuildModulePosition(std::vector<int> &moduleColumn,
                         std::vector<int> &moduleRow)
{
  moduleColumn.assign(NoModulesActive, -1);
  moduleRow.assign(NoModulesActive, -1);
  for (int slot = 0; slot < NofmodulesX * NofmodulesY; ++slot) {
    const int module = modflag[slot];
    if (module < 0) continue;
    if (module >= NoModulesActive) return false;
    moduleColumn[module] = slot % NofmodulesX;
    moduleRow[module] = slot / NofmodulesX;
  }
  return std::find(moduleColumn.begin(), moduleColumn.end(), -1) ==
             moduleColumn.end() &&
         std::find(moduleRow.begin(), moduleRow.end(), -1) == moduleRow.end();
}

bool FiberPosition(int towerID, int fiberID, bool isCherenkov,
                   const std::vector<int> &moduleColumn,
                   const std::vector<int> &moduleRow, double &x_mm,
                   double &y_mm)
{
  const int rowsPerType = NofFibersrow / 2;
  if (towerID < 0 || towerID >= NoModulesActive || fiberID < 0 ||
      fiberID >= NofFiberscolumn * rowsPerType)
    return false;

  const int column = fiberID / rowsPerType;
  const int row = 2 * (fiberID % rowsPerType) + (isCherenkov ? 1 : 0);
  const double moduleCenterX =
      dtubeX * NofFiberscolumn *
      ((NofmodulesX - 1.0) / 2.0 - moduleColumn[towerID]);
  const double moduleCenterY =
      dtubeY * NofFibersrow *
      (moduleRow[towerID] - (NofmodulesY - 1.0) / 2.0);

  x_mm = moduleCenterX + moduleX / 2.0 - tuberadius -
         2.0 * tuberadius * column - (isCherenkov ? tuberadius : 0.0);
  y_mm = moduleCenterY - moduleY / 2.0 + tuberadius + sq3 * tuberadius * row +
         tuberadius * (2.0 * sq3m1 - 1.0);
  return true;
}

struct FiberSignal {
  int totalPE{0};
  std::vector<double> arrivalTimes_ns;
};

struct PMTSignal {
  Long64_t totalPE{0};
  std::map<int, Long64_t> photoelectronsByTimeBin;

  void Add(int photoelectrons, int timeBin)
  {
    if (photoelectrons <= 0) return;
    totalPE += photoelectrons;
    photoelectronsByTimeBin[timeBin] += photoelectrons;
  }
};

double PMTTimeOfArrival(const PMTSignal &signal, int thresholdPE,
                        double timeStart_ns, double timeBinWidth_ns,
                        int timeNBins)
{
  Long64_t accumulatedPE = 0;
  for (const auto &[timeBin, photoelectrons] :
       signal.photoelectronsByTimeBin) {
    accumulatedPE += photoelectrons;
    if (accumulatedPE >= thresholdPE)
      return TimeFromBin(timeBin, timeStart_ns, timeBinWidth_ns, timeNBins);
  }
  return std::numeric_limits<double>::quiet_NaN();
}

struct FiberOutput {
  std::vector<int> towerID;
  std::vector<int> fiberID;
  std::vector<double> x_mm;
  std::vector<double> y_mm;
  std::vector<int> totalPE;
  std::vector<std::vector<double>> timeOfArrival_ns;

  void Clear()
  {
    towerID.clear();
    fiberID.clear();
    x_mm.clear();
    y_mm.clear();
    totalPE.clear();
    timeOfArrival_ns.clear();
  }
};

struct PMTOutput {
  std::vector<int> towerID;
  std::vector<Long64_t> scintillationTotalPE;
  std::vector<Long64_t> cherenkovTotalPE;
  std::vector<double> scintillationTimeOfArrival_ns;
  std::vector<double> cherenkovTimeOfArrival_ns;

  void Clear()
  {
    towerID.clear();
    scintillationTotalPE.clear();
    cherenkovTotalPE.clear();
    scintillationTimeOfArrival_ns.clear();
    cherenkovTimeOfArrival_ns.clear();
  }
};

} // namespace

void HidraDigi_v0(const char *inputFile = "DREMTubesout_Run0.root",
                  const char *outputFile = "HidraDigi_v0.root",
                  bool attenuateScintillation = false,
                  double scintillationAttenuationLength_mm = 6500.0,
                  bool attenuateCherenkov = false,
                  double cherenkovAttenuationLength_mm = 9000.0,
                  ULong_t randomSeed = 0,
                  int pmtThresholdPE = 1)
{
  constexpr double kScintillationPEPerMeV = 9.5;
  constexpr double kCherenkovPEPerPhoton = 0.153;

  if ((attenuateScintillation && scintillationAttenuationLength_mm <= 0.0) ||
      (attenuateCherenkov && cherenkovAttenuationLength_mm <= 0.0)) {
    std::cerr << "HidraDigi_v0: enabled attenuation lengths must be positive.\n";
    return;
  }
  if (pmtThresholdPE <= 0) {
    std::cerr << "HidraDigi_v0: the PMT timing threshold must be positive.\n";
    return;
  }

  TFile input(inputFile, "READ");
  if (input.IsZombie()) {
    std::cerr << "HidraDigi_v0: cannot open input file '" << inputFile << "'.\n";
    return;
  }
  auto *simulation = dynamic_cast<TTree *>(input.Get("DREMTubesout"));
  if (!simulation) {
    std::cerr << "HidraDigi_v0: tree 'DREMTubesout' was not found.\n";
    return;
  }

  const char *required[] = {
      "CherenkovTimeStart_ns", "CherenkovTimeBinWidth_ns",
      "CherenkovTimeNBins", "CherenkovDistanceStart_mm",
      "CherenkovDistanceBinWidth_mm", "CherenkovDistanceNBins",
      "CherenkovTimeTowerID", "CherenkovTimeFiberID", "CherenkovTimeBin",
      "CherenkovDistanceBin", "CherenkovTimeBinCount",
      "ScintillationTowerID", "ScintillationFiberID", "ScintillationTimeBin",
      "ScintillationDistanceBin", "ScintillationVisibleEnergy_MeV"};
  for (const char *name : required)
    if (!HasBranch(simulation, name)) return;

  double distanceStart_mm = 0.0;
  double distanceBinWidth_mm = 0.0;
  int distanceNBins = 0;
  double timeStart_ns = 0.0;
  double timeBinWidth_ns = 0.0;
  int timeNBins = 0;
  std::vector<int> *cherenkovTowerID = nullptr;
  std::vector<int> *cherenkovFiberID = nullptr;
  std::vector<int> *cherenkovTimeBin = nullptr;
  std::vector<int> *cherenkovDistanceBin = nullptr;
  std::vector<int> *cherenkovPhotonCount = nullptr;
  std::vector<int> *scintillationTowerID = nullptr;
  std::vector<int> *scintillationFiberID = nullptr;
  std::vector<int> *scintillationTimeBin = nullptr;
  std::vector<int> *scintillationDistanceBin = nullptr;
  std::vector<double> *scintillationVisibleEnergy_MeV = nullptr;

  simulation->SetBranchAddress("CherenkovTimeStart_ns", &timeStart_ns);
  simulation->SetBranchAddress("CherenkovTimeBinWidth_ns", &timeBinWidth_ns);
  simulation->SetBranchAddress("CherenkovTimeNBins", &timeNBins);
  simulation->SetBranchAddress("CherenkovDistanceStart_mm", &distanceStart_mm);
  simulation->SetBranchAddress("CherenkovDistanceBinWidth_mm", &distanceBinWidth_mm);
  simulation->SetBranchAddress("CherenkovDistanceNBins", &distanceNBins);
  simulation->SetBranchAddress("CherenkovTimeTowerID", &cherenkovTowerID);
  simulation->SetBranchAddress("CherenkovTimeFiberID", &cherenkovFiberID);
  simulation->SetBranchAddress("CherenkovTimeBin", &cherenkovTimeBin);
  simulation->SetBranchAddress("CherenkovDistanceBin", &cherenkovDistanceBin);
  simulation->SetBranchAddress("CherenkovTimeBinCount", &cherenkovPhotonCount);
  simulation->SetBranchAddress("ScintillationTowerID", &scintillationTowerID);
  simulation->SetBranchAddress("ScintillationFiberID", &scintillationFiberID);
  simulation->SetBranchAddress("ScintillationTimeBin", &scintillationTimeBin);
  simulation->SetBranchAddress("ScintillationDistanceBin", &scintillationDistanceBin);
  simulation->SetBranchAddress("ScintillationVisibleEnergy_MeV",
                               &scintillationVisibleEnergy_MeV);

  TFile output(outputFile, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "HidraDigi_v0: cannot create output file '" << outputFile << "'.\n";
    return;
  }
  output.cd();
  auto *digitised = simulation->CloneTree(0);
  digitised->SetName("HidraDigi");
  digitised->SetTitle("HidraSim calorimeter hits with digitised signals");

  std::vector<int> scintillationPE;
  std::vector<int> cherenkovPE;
  Long64_t totalScintillationPE = 0;
  Long64_t totalCherenkovPE = 0;
  Long64_t totalScintillationSiPMPE = 0;
  Long64_t totalCherenkovSiPMPE = 0;
  Long64_t totalScintillationPMTPE = 0;
  Long64_t totalCherenkovPMTPE = 0;
  FiberOutput scintillationFiberOutput;
  FiberOutput cherenkovFiberOutput;
  PMTOutput pmtOutput;
  digitised->Branch("DigiScintillationPE", &scintillationPE);
  digitised->Branch("DigiCherenkovPE", &cherenkovPE);
  digitised->Branch("DigiTotalScintillationPE", &totalScintillationPE);
  digitised->Branch("DigiTotalCherenkovPE", &totalCherenkovPE);
  digitised->Branch("DigiTotalScintillationSiPMPE",
                    &totalScintillationSiPMPE);
  digitised->Branch("DigiTotalCherenkovSiPMPE",
                    &totalCherenkovSiPMPE);
  digitised->Branch("DigiTotalScintillationPMTPE",
                    &totalScintillationPMTPE);
  digitised->Branch("DigiTotalCherenkovPMTPE",
                    &totalCherenkovPMTPE);
  digitised->Branch("DigiScintillationFiberTowerID",
                    &scintillationFiberOutput.towerID);
  digitised->Branch("DigiScintillationFiberID",
                    &scintillationFiberOutput.fiberID);
  digitised->Branch("DigiScintillationFiberX_mm",
                    &scintillationFiberOutput.x_mm);
  digitised->Branch("DigiScintillationFiberY_mm",
                    &scintillationFiberOutput.y_mm);
  digitised->Branch("DigiScintillationFiberTotalPE",
                    &scintillationFiberOutput.totalPE);
  digitised->Branch("DigiScintillationSiPMTotalPE",
                    &scintillationFiberOutput.totalPE);
  digitised->Branch("DigiScintillationFiberTimeOfArrival_ns",
                    &scintillationFiberOutput.timeOfArrival_ns);
  digitised->Branch("DigiCherenkovFiberTowerID",
                    &cherenkovFiberOutput.towerID);
  digitised->Branch("DigiCherenkovFiberID",
                    &cherenkovFiberOutput.fiberID);
  digitised->Branch("DigiCherenkovFiberX_mm",
                    &cherenkovFiberOutput.x_mm);
  digitised->Branch("DigiCherenkovFiberY_mm",
                    &cherenkovFiberOutput.y_mm);
  digitised->Branch("DigiCherenkovFiberTotalPE",
                    &cherenkovFiberOutput.totalPE);
  digitised->Branch("DigiCherenkovSiPMTotalPE",
                    &cherenkovFiberOutput.totalPE);
  digitised->Branch("DigiCherenkovFiberTimeOfArrival_ns",
                    &cherenkovFiberOutput.timeOfArrival_ns);
  digitised->Branch("DigiPMTTowerID", &pmtOutput.towerID);
  digitised->Branch("DigiScintillationPMTTotalPE",
                    &pmtOutput.scintillationTotalPE);
  digitised->Branch("DigiCherenkovPMTTotalPE",
                    &pmtOutput.cherenkovTotalPE);
  digitised->Branch("DigiScintillationPMTTimeOfArrival_ns",
                    &pmtOutput.scintillationTimeOfArrival_ns);
  digitised->Branch("DigiCherenkovPMTTimeOfArrival_ns",
                    &pmtOutput.cherenkovTimeOfArrival_ns);

  std::vector<int> moduleColumn;
  std::vector<int> moduleRow;
  if (!BuildModulePosition(moduleColumn, moduleRow)) {
    std::cerr << "HidraDigi_v0: invalid module map in HidraGeo.h.\n";
    output.Close();
    return;
  }

  TRandom3 random(randomSeed);
  const Long64_t entries = simulation->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    simulation->GetEntry(entry);
    scintillationPE.clear();
    cherenkovPE.clear();
    totalScintillationPE = 0;
    totalCherenkovPE = 0;
    totalScintillationSiPMPE = 0;
    totalCherenkovSiPMPE = 0;
    totalScintillationPMTPE = 0;
    totalCherenkovPMTPE = 0;

    std::map<std::tuple<bool, int, int>, FiberSignal> fiberSignals;
    std::vector<PMTSignal> scintillationPMTSignals(NoModulesActive);
    std::vector<PMTSignal> cherenkovPMTSignals(NoModulesActive);

    const bool invalidScintillation =
        !scintillationTowerID || !scintillationFiberID ||
        !scintillationTimeBin || !scintillationDistanceBin ||
        !scintillationVisibleEnergy_MeV ||
        scintillationTowerID->size() != scintillationVisibleEnergy_MeV->size() ||
        scintillationFiberID->size() != scintillationVisibleEnergy_MeV->size() ||
        scintillationTimeBin->size() != scintillationVisibleEnergy_MeV->size() ||
        scintillationDistanceBin->size() != scintillationVisibleEnergy_MeV->size();
    const bool invalidCherenkov =
        !cherenkovTowerID || !cherenkovFiberID || !cherenkovTimeBin ||
        !cherenkovDistanceBin || !cherenkovPhotonCount ||
        cherenkovTowerID->size() != cherenkovPhotonCount->size() ||
        cherenkovFiberID->size() != cherenkovPhotonCount->size() ||
        cherenkovTimeBin->size() != cherenkovPhotonCount->size() ||
        cherenkovDistanceBin->size() != cherenkovPhotonCount->size();
    if (invalidScintillation || invalidCherenkov) {
      std::cerr << "HidraDigi_v0: inconsistent hit-vector sizes at entry "
                << entry << ". Output was not written.\n";
      output.Close();
      return;
    }

    scintillationPE.reserve(scintillationVisibleEnergy_MeV->size());
    for (std::size_t i = 0; i < scintillationVisibleEnergy_MeV->size(); ++i) {
      const double distance = DistanceFromBin(scintillationDistanceBin->at(i),
                                              distanceStart_mm,
                                              distanceBinWidth_mm, distanceNBins);
      const int pe = DetectPhotoelectrons(
          random, scintillationVisibleEnergy_MeV->at(i) * kScintillationPEPerMeV,
          attenuateScintillation, distance, scintillationAttenuationLength_mm);
      scintillationPE.push_back(pe);
      totalScintillationPE += pe;
      const int tower = scintillationTowerID->at(i);
      if (tower < 0 || tower >= NoModulesActive) {
        std::cerr << "HidraDigi_v0: invalid scintillation tower ID at entry "
                  << entry << ". Output was not written.\n";
        output.Close();
        return;
      }
      if (IsSiPMTower(tower)) {
        totalScintillationSiPMPE += pe;
        auto &fiber =
            fiberSignals[{false, tower, scintillationFiberID->at(i)}];
        fiber.totalPE += pe;
        const double arrival_ns = TimeFromBin(scintillationTimeBin->at(i),
                                              timeStart_ns, timeBinWidth_ns,
                                              timeNBins);
        fiber.arrivalTimes_ns.insert(fiber.arrivalTimes_ns.end(), pe,
                                     arrival_ns);
      } else {
        totalScintillationPMTPE += pe;
        scintillationPMTSignals[tower].Add(pe,
                                           scintillationTimeBin->at(i));
      }
    }

    cherenkovPE.reserve(cherenkovPhotonCount->size());
    for (std::size_t i = 0; i < cherenkovPhotonCount->size(); ++i) {
      const double distance = DistanceFromBin(cherenkovDistanceBin->at(i),
                                              distanceStart_mm,
                                              distanceBinWidth_mm, distanceNBins);
      const int pe = DetectPhotoelectrons(
          random, cherenkovPhotonCount->at(i) * kCherenkovPEPerPhoton,
          attenuateCherenkov, distance, cherenkovAttenuationLength_mm);
      cherenkovPE.push_back(pe);
      totalCherenkovPE += pe;
      const int tower = cherenkovTowerID->at(i);
      if (tower < 0 || tower >= NoModulesActive) {
        std::cerr << "HidraDigi_v0: invalid Cherenkov tower ID at entry "
                  << entry << ". Output was not written.\n";
        output.Close();
        return;
      }
      if (IsSiPMTower(tower)) {
        totalCherenkovSiPMPE += pe;
        auto &fiber = fiberSignals[{true, tower, cherenkovFiberID->at(i)}];
        fiber.totalPE += pe;
        const double arrival_ns =
            TimeFromBin(cherenkovTimeBin->at(i), timeStart_ns,
                        timeBinWidth_ns, timeNBins);
        fiber.arrivalTimes_ns.insert(fiber.arrivalTimes_ns.end(), pe,
                                     arrival_ns);
      } else {
        totalCherenkovPMTPE += pe;
        cherenkovPMTSignals[tower].Add(pe, cherenkovTimeBin->at(i));
      }
    }

    scintillationFiberOutput.Clear();
    cherenkovFiberOutput.Clear();
    for (auto &[key, signal] : fiberSignals) {
      const auto [isCherenkov, tower, channel] = key;
      double x_mm = 0.0;
      double y_mm = 0.0;
      if (!FiberPosition(tower, channel, isCherenkov, moduleColumn, moduleRow,
                         x_mm, y_mm)) {
        std::cerr << "HidraDigi_v0: invalid tower/fibre ID at entry " << entry
                  << ". Output was not written.\n";
        output.Close();
        return;
      }
      auto &fiberOutput =
          isCherenkov ? cherenkovFiberOutput : scintillationFiberOutput;
      fiberOutput.towerID.push_back(tower);
      fiberOutput.fiberID.push_back(channel);
      fiberOutput.x_mm.push_back(x_mm);
      fiberOutput.y_mm.push_back(y_mm);
      fiberOutput.totalPE.push_back(signal.totalPE);
      std::sort(signal.arrivalTimes_ns.begin(), signal.arrivalTimes_ns.end());
      fiberOutput.timeOfArrival_ns.push_back(std::move(signal.arrivalTimes_ns));
    }

    pmtOutput.Clear();
    for (int tower = 0; tower < NoModulesActive; ++tower) {
      if (IsSiPMTower(tower)) continue;
      const auto &scintillationSignal = scintillationPMTSignals[tower];
      const auto &cherenkovSignal = cherenkovPMTSignals[tower];
      pmtOutput.towerID.push_back(tower);
      pmtOutput.scintillationTotalPE.push_back(
          scintillationSignal.totalPE);
      pmtOutput.cherenkovTotalPE.push_back(cherenkovSignal.totalPE);
      pmtOutput.scintillationTimeOfArrival_ns.push_back(PMTTimeOfArrival(
          scintillationSignal, pmtThresholdPE, timeStart_ns,
          timeBinWidth_ns, timeNBins));
      pmtOutput.cherenkovTimeOfArrival_ns.push_back(PMTTimeOfArrival(
          cherenkovSignal, pmtThresholdPE, timeStart_ns, timeBinWidth_ns,
          timeNBins));
    }
    digitised->Fill();
  }

  TParameter<double>("ScintillationPEPerMeV", kScintillationPEPerMeV).Write();
  TParameter<double>("CherenkovPEPerPhoton", kCherenkovPEPerPhoton).Write();
  TParameter<int>("ScintillationAttenuationEnabled", attenuateScintillation).Write();
  TParameter<double>("ScintillationAttenuationLength_mm",
                     scintillationAttenuationLength_mm).Write();
  TParameter<int>("CherenkovAttenuationEnabled", attenuateCherenkov).Write();
  TParameter<double>("CherenkovAttenuationLength_mm",
                     cherenkovAttenuationLength_mm).Write();
  // TRandom3 consumes a 32-bit seed; storing it as double avoids relying on a
  // TParameter<unsigned long> dictionary that is absent in some ROOT builds.
  TParameter<double>("DigitisationRandomSeed", static_cast<double>(randomSeed)).Write();
  TParameter<int>("PMTTimeThresholdPE", pmtThresholdPE).Write();
  digitised->Write();
  output.Close();

  std::cout << "HidraDigi_v0: wrote " << entries << " events to '"
            << outputFile << "'.\n";
}
