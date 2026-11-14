# NMF-RI Photometry

Blind spectral unmixing of two-component spectrally resolved fiber-photometry recordings
with NMF-RI (non-negative matrix factorization with reinitialization;
Jiménez-Sánchez et al., *Bioinformatics* 2020).

Given a mixed spectral time series recorded under simultaneous excitation and one initial
emission spectrum per fluorophore, the program estimates the emission spectrum of each
component (**A**) and its time course (**H**) such that **Y ≈ A·H**. The initial spectra only
start the factorization; the final spectra are estimated from the recording. No per-animal
single-excitation calibration recording is required.

<!-- Authors: fill in the MATLAB version(s) used for testing before release. -->

## Requirements

- MATLAB R2020b or later (`uifigure`, implicit expansion)
- Statistics and Machine Learning Toolbox (`ecdf`, called in `dataPreprocessing.m`)
- `NMF_RI.m` and `dataPreprocessing.m` from the original NMF-RI code (see below)

## NMF-RI code

`NMF_RI.m` and `dataPreprocessing.m` are not included in this repository. Download them from the
original NMF-RI repository and put them in this folder:

https://github.com/djimenezsanchez/NMF-RI

- Jiménez-Sánchez D, Ariz M, Morgado JM, Cortés-Domínguez I, Ortiz-de-Solórzano C.
  NMF-RI: blind spectral unmixing of highly mixed multispectral flow and image cytometry data.
  *Bioinformatics* 36(5):1590–1598 (2020). doi:10.1093/bioinformatics/btz751

## Files

| File | Description |
|---|---|
| `NMF_RI_Photometry.m` | Graphical interface: data loading, unmixing, plotting, export |
| `Initial spectrum/` | Emission spectra used as initial spectra (see the README in that folder) |
| `NMF_RI.m` | NMF-RI algorithm (download from the original repository) |
| `dataPreprocessing.m` | Dark-offset subtraction (download from the original repository) |

## Quick start

1. Put `NMF_RI_Photometry.m`, `NMF_RI.m` and `dataPreprocessing.m` in one folder on the MATLAB
   path.
2. Run `NMF_RI_Photometry`.
3. **Tab 1 (Data):** select the mixed recording (TXT) and one emission spectrum (CSV) for each
   component; set the wavelength range and component names; press
   *Load data and preview initial spectra*.
4. **Tab 2 (NMF-RI Settings):** keep the default parameters and press *Run NMF-RI*.
5. **Tab 3 (Results):** inspect the resolved spectra and time courses; press *Save All Results*.

Saved outputs: `<name>.csv` with the time course of each component (**H**), `<name>_spectra.csv`
with the resolved spectra (**A**), and one reconstructed TXT per component
(`<name>_<component>.txt`) in the original file format. Because each column of **A** sums to
one, **H** gives the intensity of each component summed over wavelengths. Each reconstructed
TXT contains A(:,k)·H(k,:) of the dark-offset-subtracted data; wavelengths outside the selected
range are written as zero.

## Input formats

**Mixed recording (TXT).** OceanView-style export. A header ending with the line
`>>>>>Begin Spectral Data<<<<<`, followed by one row of wavelengths, followed by one row per
time point: date, time, timestamp, then one intensity value per wavelength.

**Initial spectrum (CSV).** First column: wavelength (nm). The emission column is the first
column whose header contains `em` and does not contain `ex` or `2p`
(for example `EGFP em`, `jRGECO1a (calcium saturated) em`, `Emission`). Any source can be used
(published spectra such as CSV files exported from FPbase, https://www.fpbase.org, or in vitro
measurements). Values are baseline-subtracted, peak-normalized and interpolated onto the
recording's wavelength axis; wavelengths outside the CSV range are set to zero. Before unmixing,
values ≤ 0 are replaced by 1e-4 and each spectrum is normalized to unit sum. The spectra used in
this study are in `Initial spectrum/`.

## Parameters

| Parameter (GUI) | Default | Meaning |
|---|---|---|
| Wavelength range | 450–650 nm | Spectral window used for unmixing |
| theta1 | 0.01 | Intra-layer criterion: a reinitialization layer is closed when the change of **A** falls below theta1 or increases |
| theta2 | 0.01 | Inter-layer criterion: the algorithm stops when **A** changes by less than theta2 between layers |
| alphaA | 0.001 | Sparseness weight for **A** (relative to the mean of **Y**) |
| alphaX | 0.001 | Sparseness weight for **H** (relative to the mean of **Y**) |

The maximum number of iterations is 300. The number of iterations is printed in the MATLAB
command window.

## Method

Each recorded spectrum is modeled as a non-negative mixture **Y = A·H**. The minimum across
wavelengths is first subtracted from each spectrum to remove the dark offset
(`dataPreprocessing.m`). **A** is initialized with the initial spectra normalized to unit sum
(**A0**), and **H** with their least-squares solution. **A** and **H** are then updated with the
multiplicative updates and reinitialization layers of NMF-RI. The full time course is finally
computed from the resolved spectra by least squares.

## Citation

