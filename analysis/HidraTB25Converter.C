#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <iostream>
#include <stdint.h>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>
#include <nlohmann/json.hpp>
#include "HidraGeo.h"

using json = nlohmann::json;
const bool printMapCheck = false;

// -----------------------------------------------------------------------------
// Global constants and parameters
// -----------------------------------------------------------------------------

const double sciPheGeV = 119.001;
const double cerPheGeV = 29.4;

const unsigned int grouping = 8;
const unsigned int ChannelsPerFers = 64;
const unsigned int TotalFers = 16;
const unsigned int TotalCalibChannels = ChannelsPerFers * TotalFers;
const unsigned int FersMultiplicity = 2;

const unsigned int NoiseRandomSeed = 12345;
const bool ClampNegativeAdcToZero = true;

std::map<unsigned int, double> fers_to_thr_map = {
    { 1, 0.17 },
    { 2, 0.20 },
    { 3, 0.30 },
    { 4, 0.28 },
    { 5, 0.105},
    { 6, 0.095},
    { 7, 0.08 },
    { 8, 0.08 },
    { 9, 0.05 },
    {10, 0.07 },
    {11, 0.05 },
    {12, 0.025},
    {13, 0.06 },
    {14, 0.05 },
    {15, 0.12 },
    {16, 0.06 }
};

// -----------------------------------------------------------------------------
// Calibration and map structures
// -----------------------------------------------------------------------------

struct PedestalHgEntry {
  double medianAdc = 0.0;
  double rmsAdc = 0.0;
};

using PedestalHgVector = std::vector<PedestalHgEntry>;
using AdcToGeVVector = std::vector<double>;
using FersKey = uint64_t;

struct SipmMapEntry {
  int calibIndex = -1;
  int boardID = -1;
  std::string type;
  int row = -1;
  int column = -1;
  double x = 0.0;
  double y = 0.0;
  int fersId = -1;
  int ch = -1;
  std::string module_name;
  double x_local = 0.0;
  double y_local = 0.0;
};

using SipmLookup = std::unordered_map<std::string, SipmMapEntry>;

// -----------------------------------------------------------------------------
// Utility helpers
// -----------------------------------------------------------------------------

std::string make_sipm_key(const std::string& tower,
                          const std::string& type,
                          int row,
                          int column)
{
  std::ostringstream os;
  os << tower << '|' << type << '|' << row << '|' << column;
  return os.str();
}

int grouped_column(int colID_original, unsigned int groupingValue)
{
  return static_cast<int>(colID_original / groupingValue);
}

FersKey make_fers_key(int boardID, int fersId)
{
  return (static_cast<FersKey>(static_cast<uint32_t>(boardID)) << 32) |
         static_cast<uint32_t>(fersId);
}

double GetFersThresholdGeV(int fersId)
{
  if (fersId < 0) {
    throw std::runtime_error("Invalid FERS id: " + std::to_string(fersId));
  }

  const auto it = fers_to_thr_map.find(static_cast<unsigned int>(fersId));
  if (it == fers_to_thr_map.end()) {
    throw std::runtime_error("Missing threshold for FERS id: " + std::to_string(fersId));
  }

  return it->second;
}

std::string simTower_to_tbTower(const std::string& simTower)
{
  return std::to_string(306 + std::stoi(simTower));
}

const SipmMapEntry* FindSipmInfo(const SipmLookup& lookup,
                                 const std::string& tower,
                                 const std::string& type,
                                 int row,
                                 int groupedCol)
{
  const auto key = make_sipm_key(tower, type, row, groupedCol);
  const auto it = lookup.find(key);
  return (it == lookup.end()) ? nullptr : &it->second;
}

// -----------------------------------------------------------------------------
// Calibration loaders
// -----------------------------------------------------------------------------

PedestalHgVector LoadHgPedestalsFromJson(const std::string& jsonPath)
{
  std::ifstream in(jsonPath);
  if (!in) {
    throw std::runtime_error("Cannot open pedestal JSON: " + jsonPath);
  }

  json j;
  in >> j;

  if (!j.is_object()) {
    throw std::runtime_error("Pedestal JSON must be an object: " + jsonPath);
  }

  std::size_t maxIndex = 0;
  for (auto it = j.begin(); it != j.end(); ++it) {
    const int idx = std::stoi(it.key());
    if (idx < 0) {
      throw std::runtime_error("Negative pedestal index in: " + jsonPath);
    }
    maxIndex = std::max(maxIndex, static_cast<std::size_t>(idx));
  }

  PedestalHgVector values(maxIndex + 1);
  std::vector<bool> found(maxIndex + 1, false);

  for (auto it = j.begin(); it != j.end(); ++it) {
    const int idx = std::stoi(it.key());
    const json& node = it.value();

    values[static_cast<std::size_t>(idx)].medianAdc = node.value("median_HG", 0.0);
    values[static_cast<std::size_t>(idx)].rmsAdc = node.value("iqr_eff_HG", 0.0);
    found[static_cast<std::size_t>(idx)] = true;
  }

  for (std::size_t i = 0; i < found.size(); ++i) {
    if (!found[i]) {
      std::ostringstream os;
      os << "Missing HG pedestal entry for channel " << i
         << " in " << jsonPath;
      throw std::runtime_error(os.str());
    }
  }

  return values;
}

AdcToGeVVector LoadAdcToGeVFromJson(const std::string& jsonPath)
{
  std::ifstream in(jsonPath);
  if (!in) {
    throw std::runtime_error("Cannot open ADC->GeV JSON: " + jsonPath);
  }

  json j;
  in >> j;

  if (!j.is_array()) {
    throw std::runtime_error("ADC->GeV JSON must be an array: " + jsonPath);
  }

  AdcToGeVVector values;
  values.reserve(j.size());

  for (const auto& item : j) {
    values.push_back(item.get<double>());
  }

  if (values.empty()) {
    throw std::runtime_error("ADC->GeV JSON is empty: " + jsonPath);
  }

  return values;
}

SipmLookup LoadSipmMap(const std::string& jsonPath)
{
  std::ifstream in(jsonPath);
  if (!in) {
    throw std::runtime_error("Cannot open SiPM map JSON: " + jsonPath);
  }

  json j;
  in >> j;

  SipmLookup lookup;
  for (auto it = j.begin(); it != j.end(); ++it) {
    const json& node = it.value();

    SipmMapEntry entry;
    entry.calibIndex  = std::stoi(it.key());
    entry.boardID     = node.value("boardID", -1);
    entry.type        = node.value("type", "");
    entry.row         = node.value("row", -1);
    entry.column      = node.value("column", -1);
    entry.x           = node.value("x", 0.0);
    entry.y           = node.value("y", 0.0);
    entry.fersId      = node.value("fersId", -1);
    entry.ch          = node.value("ch", -1);
    entry.module_name = node.value("module_name", "");
    entry.x_local     = node.value("x_local", 0.0);
    entry.y_local     = node.value("y_local", 0.0);

    lookup[make_sipm_key(entry.module_name, entry.type, entry.row, entry.column)] = entry;
  }

  return lookup;
}

void ValidateCalibrationSizes(const PedestalHgVector& pedestalsHg,
                              const AdcToGeVVector& adcToGeV)
{
  if (pedestalsHg.size() != TotalCalibChannels) {
    std::ostringstream os;
    os << "Expected " << TotalCalibChannels
       << " HG pedestal entries, got " << pedestalsHg.size();
    throw std::runtime_error(os.str());
  }

  if (adcToGeV.size() != TotalCalibChannels) {
    std::ostringstream os;
    os << "Expected " << TotalCalibChannels
       << " ADC->GeV entries, got " << adcToGeV.size();
    throw std::runtime_error(os.str());
  }
}

void ValidateMapEntry(const SipmMapEntry& entry)
{
  if (entry.calibIndex < 0 || entry.calibIndex >= static_cast<int>(TotalCalibChannels)) {
    std::ostringstream os;
    os << "Invalid calibIndex " << entry.calibIndex;
    throw std::runtime_error(os.str());
  }
}

// -----------------------------------------------------------------------------
// Conversion helpers
// -----------------------------------------------------------------------------

double ConvertPheToGeV(double rawPhe, const std::string& type)
{
  if (type == "S") {
    return rawPhe / sciPheGeV;
  }
  if (type == "C") {
    return rawPhe / cerPheGeV;
  }

  throw std::runtime_error("Unknown SiPM type for phe->GeV conversion: " + type);
}

double ConvertGeVToAdc(double energyGeV, double adcToGeV)
{
  if (adcToGeV <= 0.0) {
    throw std::runtime_error("Non-positive ADC->GeV calibration encountered");
  }
  return energyGeV / adcToGeV;
}

void AddMappedSignal(const SipmLookup& sipmLookup,
                     const std::string& tbTower,
                     const std::string& type,
                     int row,
                     int col,
                     double rawPhe,
                     const AdcToGeVVector& adcToGeV,
                     std::vector<double>& signalGeV,
                     std::vector<double>& signalAdc)
{
  const int groupedCol = grouped_column(col, grouping);
  const SipmMapEntry* info = FindSipmInfo(sipmLookup, tbTower, type, row, groupedCol);
  if (!info) {
    return;
  }

  ValidateMapEntry(*info);
  const std::size_t idx = static_cast<std::size_t>(info->calibIndex);

  const double energyGeV = ConvertPheToGeV(rawPhe, type);
  const double adcSignal = ConvertGeVToAdc(energyGeV, adcToGeV[idx]);

  signalGeV[idx] += energyGeV;
  signalAdc[idx] += adcSignal;
}

void BuildRawHgAdcWithActivation(const std::vector<double>& signalGeV,
                                 std::vector<double>& signalAdc,
                                 const PedestalHgVector& pedestalsHg,
                                 const AdcToGeVVector& adcToGeV,
                                 const std::vector<int>& boardID,
                                 const std::vector<int>& fersId,
                                 TRandom3& rng,
                                 std::vector<double>& hgAdc,
                                 std::vector<int>& hgAdcInt)
{
  (void)signalGeV; // intentionally preserved even when FERS is not activated

  std::vector<double> channelNoiseAdc(signalAdc.size(), 0.0);
  std::unordered_map<FersKey, unsigned int> channelsOverThreshold;

  for (std::size_t idx = 0; idx < signalAdc.size(); ++idx) {
    if (boardID[idx] < 0 || fersId[idx] <= 0) {
      continue;
    }

    const double noiseSigmaAdc = pedestalsHg[idx].rmsAdc;
    const double noiseAdc = rng.Gaus(0.0, noiseSigmaAdc);
    channelNoiseAdc[idx] = noiseAdc;

    const double thresholdGeV = GetFersThresholdGeV(fersId[idx]);
    const double thresholdAdc = ConvertGeVToAdc(thresholdGeV, adcToGeV[idx]);
    const double thresholdDecisionAdc = signalAdc[idx] + noiseAdc;

    if (thresholdDecisionAdc > thresholdAdc) {
      ++channelsOverThreshold[make_fers_key(boardID[idx], fersId[idx])];
    }
  }

  for (std::size_t idx = 0; idx < signalAdc.size(); ++idx) {
    if (boardID[idx] < 0 || fersId[idx] <= 0) {
      signalAdc[idx] = 0.0;
      hgAdc[idx] = 0.0;
      hgAdcInt[idx] = 0;
      continue;
    }

    const FersKey key = make_fers_key(boardID[idx], fersId[idx]);
    const auto it = channelsOverThreshold.find(key);
    const bool isActive = (it != channelsOverThreshold.end()) && (it->second >= FersMultiplicity);

    if (!isActive) {
      signalAdc[idx] = 0.0;
      hgAdc[idx] = 0.0;
      hgAdcInt[idx] = 0;
      continue;
    }

    const double pedestal = pedestalsHg[idx].medianAdc;
    const double rawHg = signalAdc[idx] + pedestal + channelNoiseAdc[idx];
    hgAdc[idx] = ClampNegativeAdcToZero ? std::max(0.0, rawHg) : rawHg;
    hgAdcInt[idx] = static_cast<int>(std::llround(hgAdc[idx]));
  }
}

// -----------------------------------------------------------------------------
// Main converter
// -----------------------------------------------------------------------------

void HidraTB25Converter(double energy, const std::string& input)
{
  const std::string sipmMapPath =
    "/home/apareti/HidraSim2025/TBDataPreparation/2025_SPS/MapAndCalibration/sipm_map.json";

  const std::string pedestalJsonPath =
    "/home/apareti/HidraSim2025/TBDataPreparation/2025_SPS/MapAndCalibration/SiPM_pedestals_v1.json";

  const std::string adcToGeVJsonPath =
    "/home/apareti/HidraSim2025/TBDataPreparation/2025_SPS/MapAndCalibration/SiPM_ADCtoGeV_v1.json";

  SipmLookup sipmLookup;
  PedestalHgVector pedestalsHg;
  AdcToGeVVector adcToGeV;

  try {
    sipmLookup = LoadSipmMap(sipmMapPath);
    pedestalsHg = LoadHgPedestalsFromJson(pedestalJsonPath);
    adcToGeV = LoadAdcToGeVFromJson(adcToGeVJsonPath);
    ValidateCalibrationSizes(pedestalsHg, adcToGeV);

    std::cout << "Loaded " << sipmLookup.size()
              << " SiPM map entries from " << sipmMapPath << std::endl;
    std::cout << "Loaded " << pedestalsHg.size()
              << " HG pedestal entries from " << pedestalJsonPath << std::endl;
    std::cout << "Loaded " << adcToGeV.size()
              << " ADC->GeV entries from " << adcToGeVJsonPath << std::endl;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return;
  }

  const std::string infile = "../build/" + input;
  std::cout << "Using file: " << infile << std::endl;

  TFile* simfile = TFile::Open(infile.c_str(), "READ");
  if (!simfile || simfile->IsZombie()) {
    std::cerr << "Cannot open input file " << infile << std::endl;
    return;
  }

  TTree* simtree = static_cast<TTree*>(simfile->Get("DREMTubesout"));
  if (!simtree) {
    std::cerr << "Cannot find TTree DREMTubesout in " << infile << std::endl;
    simfile->Close();
    delete simfile;
    return;
  }

  std::ostringstream enss;
  enss << energy;
  const std::string outfile = "HidraSimTB25_converted" + enss.str() + ".root";
  TFile fout(outfile.c_str(), "RECREATE");

  TTree* outtree = new TTree("Phys2025", "Simulation converted to TB-like HG ADC channels");
  TTree* metaTree = new TTree("ChannelMap", "Static channel metadata and calibrations");

  unsigned int event = 0;
  int pdg = 0;
  double primaryEnergy = 0.0;
  double beamX = 0.0;
  double beamY = 0.0;

  std::vector<double> signalGeV(TotalCalibChannels, 0.0);
  std::vector<double> signalAdc(TotalCalibChannels, 0.0);
  std::vector<double> hgAdc(TotalCalibChannels, 0.0);
  std::vector<int> hgAdcInt(TotalCalibChannels, 0);

  std::vector<int> boardID(TotalCalibChannels, -1);
  std::vector<int> fersId(TotalCalibChannels, -1);
  std::vector<int> chID(TotalCalibChannels, -1);
  std::vector<int> fiberIsCher(TotalCalibChannels, 0);
  std::vector<double> pedestalMedianHgAdc(TotalCalibChannels, 0.0);
  std::vector<double> pedestalRmsHgAdc(TotalCalibChannels, 0.0);

  for (const auto& [key, entry] : sipmLookup) {
    if (entry.calibIndex < 0 || entry.calibIndex >= static_cast<int>(TotalCalibChannels)) {
      continue;
    }
    const std::size_t idx = static_cast<std::size_t>(entry.calibIndex);
    boardID[idx] = entry.boardID;
    fersId[idx] = entry.fersId;
    chID[idx] = entry.ch;
    fiberIsCher[idx] = (entry.type == "C") ? 1 : 0;
  }

  outtree->Branch("event", &event);
  outtree->Branch("PrimaryPDGID", &pdg);
  outtree->Branch("PrimaryParticleEnergy", &primaryEnergy);
  outtree->Branch("PrimaryX", &beamX);
  outtree->Branch("PrimaryY", &beamY);

  outtree->Branch("signalGeV", &signalGeV);
  outtree->Branch("signalAdc", &signalAdc); // phe -> GeV -> ADC before pedestal, but zero-filled if FERS is off
  outtree->Branch("hgAdc", &hgAdc);         // final HG ADC after adding pedestal and smearing, zero-filled if FERS is off
  outtree->Branch("hgAdcInt", &hgAdcInt);   // final HG ADC as integer, zero-filled if FERS is off

  for (std::size_t idx = 0; idx < pedestalsHg.size(); ++idx) {
    pedestalMedianHgAdc[idx] = pedestalsHg[idx].medianAdc;
    pedestalRmsHgAdc[idx] = pedestalsHg[idx].rmsAdc;
  }

  metaTree->Branch("boardID", &boardID);
  metaTree->Branch("fersId", &fersId);
  metaTree->Branch("chID", &chID);
  metaTree->Branch("fiberIsCher", &fiberIsCher);
  metaTree->Branch("pedestalMedianHgAdc", &pedestalMedianHgAdc);
  metaTree->Branch("pedestalRmsHgAdc", &pedestalRmsHgAdc);
  metaTree->Branch("adcToGeV", &adcToGeV);
  metaTree->Fill();

  double lenergy = 0.0;
  double denergy = 0.0;
  double edep = 0.0;
  double Stot = 0.0;
  double Ctot = 0.0;
  double PSdep = 0.0;
  std::vector<double>* TowerE = nullptr;
  std::vector<double>* SPMT = nullptr;
  std::vector<double>* CPMT = nullptr;
  std::vector<double>* SSiPM = nullptr;
  std::vector<double>* CSiPM = nullptr;

  simtree->SetBranchAddress("PrimaryPDGID", &pdg);
  simtree->SetBranchAddress("PrimaryParticleEnergy", &primaryEnergy);
  simtree->SetBranchAddress("EscapedEnergyl", &lenergy);
  simtree->SetBranchAddress("EscapedEnergyd", &denergy);
  simtree->SetBranchAddress("EnergyTot", &edep);
  simtree->SetBranchAddress("NofPMTScinDet", &Stot);
  simtree->SetBranchAddress("NofPMTCherDet", &Ctot);
  simtree->SetBranchAddress("PSEnergy", &PSdep);
  simtree->SetBranchAddress("PrimaryX", &beamX);
  simtree->SetBranchAddress("PrimaryY", &beamY);
  simtree->SetBranchAddress("VecTowerE", &TowerE);
  simtree->SetBranchAddress("VecSPMT", &SPMT);
  simtree->SetBranchAddress("VecCPMT", &CPMT);
  simtree->SetBranchAddress("VectorSignals", &SSiPM);
  simtree->SetBranchAddress("VectorSignalsCher", &CSiPM);

  TRandom3 rng(NoiseRandomSeed);

  const Long64_t nentries = simtree->GetEntries();
  std::cout << "Entries " << nentries << std::endl;

  for (Long64_t i = 0; i < nentries; ++i) {
    simtree->GetEntry(i);
    event = static_cast<unsigned int>(i);

    std::fill(signalGeV.begin(), signalGeV.end(), 0.0);
    std::fill(signalAdc.begin(), signalAdc.end(), 0.0);
    std::fill(hgAdc.begin(), hgAdc.end(), 0.0);
    std::fill(hgAdcInt.begin(), hgAdcInt.end(), 0);

    for (unsigned int n = 0; n < SSiPM->size(); ++n) {
      const double rawPhe = SSiPM->at(n);
      if (rawPhe == 0.0) {
        continue;
      }

      const unsigned int towID  = static_cast<unsigned int>(n / (NofFiberscolumn * NofFibersrow / 2));
      const unsigned int SiPMID = n % (NofFiberscolumn * NofFibersrow / 2);
      const unsigned int colID  = NofFiberscolumn - static_cast<unsigned int>(SiPMID / (NofFibersrow / 2));
      const unsigned int rowID  = 2 * static_cast<unsigned int>(SiPMID % (NofFibersrow / 2));

      const int groupedCol = grouped_column(static_cast<int>(colID), grouping);
      const std::string tbTower = simTower_to_tbTower(std::to_string(towID));
      AddMappedSignal(
        sipmLookup,
        tbTower,
        "S",
        static_cast<int>(rowID),
        static_cast<int>(colID),
        rawPhe,
        adcToGeV,
        signalGeV,
        signalAdc
      );

        const SipmMapEntry* mapped =
        FindSipmInfo(sipmLookup, tbTower, "S", static_cast<int>(rowID), groupedCol);

        if(printMapCheck) {
            if (mapped) {
            std::cout
                << "Event " << event
                << " | SIM S:"
                << " tower=" << tbTower
                << " row=" << rowID
                << " colRaw=" << colID
                << " colGrouped=" << groupedCol
                << " rawPhe=" << rawPhe
                << "  -->  MAP:"
                << " calibIndex=" << mapped->calibIndex
                << " module_name=" << mapped->module_name
                << " type=" << mapped->type
                << " row=" << mapped->row
                << " column=" << mapped->column
                << " boardID=" << mapped->boardID
                << " fersId=" << mapped->fersId
                << " ch=" << mapped->ch
                << std::endl;
            } else {
            std::cout
                << "Event " << event
                << " | SIM S:"
                << " tower=" << tbTower
                << " row=" << rowID
                << " colRaw=" << colID
                << " colGrouped=" << groupedCol
                << " rawPhe=" << rawPhe
                << "  -->  MAP: NOT FOUND"
                << std::endl;
            }
        }


    }

    for (unsigned int n = 0; n < CSiPM->size(); ++n) {
      const double rawPhe = CSiPM->at(n);
      if (rawPhe == 0.0) {
        continue;
      }

      const unsigned int towID  = static_cast<unsigned int>(n / (NofFiberscolumn * NofFibersrow / 2));
      const unsigned int SiPMID = n % (NofFiberscolumn * NofFibersrow / 2);
      const unsigned int colID  = NofFiberscolumn - static_cast<unsigned int>(SiPMID / (NofFibersrow / 2));
      const unsigned int rowID  = 2 * static_cast<unsigned int>(SiPMID % (NofFibersrow / 2)) + 1;

      const std::string tbTower = simTower_to_tbTower(std::to_string(towID));
      AddMappedSignal(
        sipmLookup,
        tbTower,
        "C",
        static_cast<int>(rowID),
        static_cast<int>(colID),
        rawPhe,
        adcToGeV,
        signalGeV,
        signalAdc
      );
    }

    BuildRawHgAdcWithActivation(
      signalGeV,
      signalAdc,
      pedestalsHg,
      adcToGeV,
      boardID,
      fersId,
      rng,
      hgAdc,
      hgAdcInt
    );

    outtree->Fill();
    //break;
  }

  fout.Write();
  fout.Close();

  simfile->Close();
  delete simfile;

  std::cout << "Wrote converted TB-like file: " << outfile << std::endl;
}  