# HiDRa Analysis

This directory contains the shared geometry header and the standard ROOT macro
used to analyze HidraSim simulation ntuples.

## Digitisation

`HidraDigi_v0.C` reads the binned calorimeter hits from the `DREMTubesout`
tree and converts them to detected photoelectrons. Run it after the Geant4
simulation and before the higher-level analysis. The calorimeter sensitive
detector records raw timing and longitudinal-distance cells from all fibre
cores. The digitizer retains fibre-level channels for modules listed by
`SiPMMod` and groups all fibres in every other module into one scintillation
PMT channel and one Cherenkov PMT channel.

From the `analysis/` directory, the default command is:

```bash
root -l -b -q 'HidraDigi_v0.C("../build/DREMTubesout_Run0.root")'
```

This creates `HidraDigi_v0.root` in the current directory. Light attenuation
is disabled by default. The scintillation response is sampled with a Poisson
mean of 9.5 p.e./MeV of Birks-quenched visible energy, while the Cherenkov
response uses a mean of 0.153 p.e. per trapped photon. These are the values
used by the standard `DREMTubesSignalHelper` processing.

The complete call signature is:

```cpp
HidraDigi_v0(inputFile,
             outputFile,
             attenuateScintillation,
             scintillationAttenuationLength_mm,
             attenuateCherenkov,
             cherenkovAttenuationLength_mm,
             randomSeed,
             pmtThresholdPE)
```

For example, to enable 6.5 m scintillation attenuation and 9 m Cherenkov
attenuation with a reproducible random seed:

```bash
root -l -b -q 'HidraDigi_v0.C("../build/DREMTubesout_Run0.root","HidraDigi_attenuated.root",true,6500.,true,9000.,12345)'
```

The attenuation lengths are specified in millimetres and must be positive
when their corresponding switch is enabled. Each already-smeared
photoelectron survives with probability `exp(-distance/attenuationLength)`,
matching the ordering in the standard simulation processing. A seed of zero
uses ROOT's automatically generated seed; use a nonzero value for reproducible
digitisation. `pmtThresholdPE` is the cumulative detected-photoelectron
threshold used to extract the single PMT time of arrival. It defaults to one
photoelectron and must be positive.

The output tree is named `HidraDigi`. It retains all input branches and adds:

- `DigiScintillationPE`: digitised scintillation p.e. for each input hit cell.
- `DigiCherenkovPE`: digitised Cherenkov p.e. for each input hit cell.
- `DigiTotalScintillationPE`: event-wide scintillation p.e. sum after
  Poisson conversion and, when enabled, light attenuation.
- `DigiTotalCherenkovPE`: event-wide Cherenkov p.e. sum after Poisson
  conversion and, when enabled, light attenuation.
- `DigiTotalScintillationSiPMPE` and `DigiTotalCherenkovSiPMPE`: event-wide
  totals containing only SiPM-readout channels.
- `DigiTotalScintillationPMTPE` and `DigiTotalCherenkovPMTPE`: event-wide
  totals containing only PMT-readout channels.
- `DigiScintillationFiberTowerID` and `DigiScintillationFiberID`: tower and
  fibre identifiers for each scintillation-channel record.
- `DigiCherenkovFiberTowerID` and `DigiCherenkovFiberID`: tower and fibre
  identifiers for each Cherenkov-channel record.
- `DigiScintillationFiberX_mm`, `DigiScintillationFiberY_mm`,
  `DigiCherenkovFiberX_mm`, and `DigiCherenkovFiberY_mm`: channel positions at
  the sensor plane, decoded with the geometry in `HidraGeo.h`.
- `DigiScintillationFiberTotalPE` and `DigiCherenkovFiberTotalPE`: total
  converted and, when enabled, attenuated p.e. collected in each individual
  SiPM channel. `DigiScintillationSiPMTotalPE` and
  `DigiCherenkovSiPMTotalPE` are explicit aliases for these same
  channel-level vectors.
- `DigiScintillationFiberTimeOfArrival_ns` and
  `DigiCherenkovFiberTimeOfArrival_ns`: one arrival-time array per channel.
  Each element is the centre of the simulation arrival-time bin for one
  detected photoelectron. The values are time ordered, and each array length
  equals the corresponding channel entry in the matching `FiberTotalPE`
  vector.
- `DigiPMTTowerID`: module identifier for every PMT-readout module. This
  vector is parallel to all `Digi...PMT...` vectors.
- `DigiScintillationPMTTotalPE` and `DigiCherenkovPMTTotalPE`: one total
  detected signal per PMT module after all fibres have been grouped. When
  attenuation is enabled, only surviving photoelectrons enter these totals.
- `DigiScintillationPMTTimeOfArrival_ns` and
  `DigiCherenkovPMTTimeOfArrival_ns`: the centre of the first time bin whose
  cumulative grouped PMT signal reaches `pmtThresholdPE`. The value is `NaN`
  when that PMT does not reach threshold.

The scintillation `DigiScintillationFiber...` vectors are mutually parallel,
as are the Cherenkov `DigiCherenkovFiber...` vectors. Index `i` describes one
channel of the corresponding type. Channels that have simulated calorimeter
hits but fluctuate to zero detected photoelectrons are retained with an empty
arrival-time array and zero total signal. The arrival-time values inherit the
simulation's underflow, regular, and overflow time-bin convention.

For each signal type, the sum of the SiPM fibre totals and grouped PMT totals
equals the corresponding event-wide `DigiTotal...PE`. Attenuation is applied
to every raw time-distance cell before either fibre-level or PMT-level
aggregation, so it consistently affects the channel totals and timing.
The combined totals also satisfy
`DigiTotalScintillationPE = DigiTotalScintillationSiPMPE +
DigiTotalScintillationPMTPE`, with the analogous identity for Cherenkov.

The conversion constants, attenuation configuration, attenuation lengths, and
random seed are also written into the output ROOT file as metadata.
`PMTTimeThresholdPE` records the threshold used for PMT timing.

## Standard Simulation Analysis

`HidraAna.C` is the general analysis macro for HidraSim simulation ntuples.
It reads the simulation output tree, reconstructs the calorimeter response, and
produces ROOT histograms for the requested beam energy.

The workflow is:

1. Open the input ROOT file and read the `DREMTubes` simulation tree.
2. Load the detector geometry from `HidraGeo.h`.
3. Loop over events and accumulate PMT tower signals and SiPM fiber signals.
4. Convert scintillation and Cherenkov SiPM photoelectron yields to energy.
5. Reconstruct SiPM channel positions from the simulation vector index.
6. Fill calorimeter maps, SiPM coordinate maps, energy histograms, and optional
   event-display histograms.

Inputs:

- Beam energy, passed as the first macro argument.
- HidraSim ROOT ntuple, passed as the second macro argument.
- Geometry constants from `HidraGeo.h`.

Outputs:

- A ROOT histogram file named from the input energy, for example `hidra10.root`.
- Energy-response histograms for scintillation, Cherenkov, and combined signals.
- Calorimeter tower maps and SiPM coordinate maps.
- Optional event-display histograms written at the configured event interval.

Typical usage from this directory:

```bash
root -l -b -q 'HidraAna.C(energy, "input.root")'
```

The macro uses `HidraGeo.h` for the detector geometry constants. SiPM channels
are decoded from the simulation vector index using the same convention used when
the SiPM copy numbers are assigned in the Geant4 geometry. The resulting SiPM
coordinates are filled in the simulation coordinate system.

Generated ROOT files, logs, summaries, compiled ROOT dictionaries, and plots are
analysis outputs and should normally not be committed.

## TB25 Analysis

TB25-specific analysis files live in `tb25_analysis/`. See
`tb25_analysis/README.md` for the TB25 workflow, mapping inputs, and configurable
analysis switches.
