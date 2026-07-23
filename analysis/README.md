# HiDRa Analysis

This directory contains the shared geometry header and the standard ROOT macro
used to analyze HidraSim simulation ntuples.

## Digitisation

`HidraDigi_v0.C` reads the binned calorimeter hits from the `DREMTubesout`
tree and converts them to detected photoelectrons. Run it after the Geant4
simulation and before the higher-level analysis.

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
             randomSeed)
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
digitisation.

The output tree is named `HidraDigi`. It retains all input branches and adds:

- `DigiScintillationPE`: digitised scintillation p.e. for each input hit cell.
- `DigiCherenkovPE`: digitised Cherenkov p.e. for each input hit cell.
- `DigiTotalScintillationPE`: event-wide scintillation p.e. sum.
- `DigiTotalCherenkovPE`: event-wide Cherenkov p.e. sum.
- `DigiFiberTowerID` and `DigiFiberID`: tower and fibre identifiers for each
  fibre-level record.
- `DigiFiberIsCherenkov`: 0 for a scintillation fibre and 1 for a Cherenkov
  fibre.
- `DigiFiberX_mm` and `DigiFiberY_mm`: fibre position at the sensor plane,
  decoded with the geometry in `HidraGeo.h`.
- `DigiFiberTotalPE`: total converted and, when enabled, attenuated signal in
  that fibre.
- `DigiFiberTimeOfArrival_ns`: one arrival-time array per fibre. Each element
  is the centre of the simulation arrival-time bin for one detected
  photoelectron. The values are time ordered, and the array length equals
  `DigiFiberTotalPE`.

All `DigiFiber...` vectors are parallel: index `i` describes one fibre. Fibres
that have simulated calorimeter hits but fluctuate to zero detected
photoelectrons are retained with an empty arrival-time array and zero total
signal. The arrival-time values inherit the simulation's underflow, regular,
and overflow time-bin convention.

The conversion constants, attenuation configuration, attenuation lengths, and
random seed are also written into the output ROOT file as metadata.

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
