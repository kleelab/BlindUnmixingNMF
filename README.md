# Blind spectral unmixing for spectrally resolved fiber photometry

MATLAB code for blind spectral unmixing of two fluorescent sensors recorded simultaneously with
spectrally resolved fiber photometry, and for its validation. The mixed emission spectrum
recorded at every time point is separated into the spectrum and time course of each sensor with
NMF-RI (non-negative matrix factorization with reinitialization; Jiménez-Sánchez et al.,
*Bioinformatics* 2020). The spectra are estimated from the recording, and the initial spectra
(for example published spectra from FPbase) only start the factorization.

<!-- Authors: fill in the MATLAB version(s) used for testing before release. -->

## Folders

| Folder | Purpose | Main files |
|---|---|---|
| [`NMF`](NMF) | Unmixing of mixed spectral recordings with a graphical interface | `NMF_RI_Photometry.m`, `NMF_RI.m`, `dataPreprocessing.m` |
| [`validation`](validation) | Validation against a known ground truth with in silico mixtures of recordings from mice expressing a single fluorophore | `NMF_Hybrid_Validation.m`, `NMF_Hybrid_Equivalence.m` |
| [`simulation`](simulation) | Accuracy of unmixing at different AUC ratios of the two sensors, with a power decoder built from spectra measured at several excitation powers | `NMF_Power_Decoder.m`, `NMF_Spike_Recovery.m` |

Each folder has its own README with the input formats, method, settings and outputs.

## Requirements

- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox (`ecdf`, called in `dataPreprocessing.m`)

## Installation

Download or clone this repository and add the three folders to the MATLAB path:

```matlab
addpath('<repository>/NMF', '<repository>/validation', '<repository>/simulation');
```

The scripts in `validation` and `simulation` call `NMF_RI.m` and `dataPreprocessing.m` from the
`NMF` folder, so the `NMF` folder must always be on the path.

## Workflow

1. **Unmixing (`NMF`).** Run `NMF_RI_Photometry`, select a mixed recording (TXT) and one initial
   emission spectrum (CSV) for each sensor, and run NMF-RI. The resolved spectra, the time
   course of each component and one reconstructed TXT per component are saved.
2. **Validation (`validation`).** Run `NMF_Hybrid_Validation` with recordings from mice that
   express a single sensor (or a sensor and a control fluorescent protein). The recordings are
   summed in silico, unmixed with the same pipeline, and compared with the original recordings
   before mixing. Run `NMF_Hybrid_Equivalence` on the saved `.mat` file for agreement (CCC),
   Bland–Altman and equivalence (TOST) statistics.
3. **Different AUC ratios (`simulation`).** Run `NMF_Power_Decoder` once for each sensor with
   spectra measured at several excitation powers. Then run `NMF_Spike_Recovery`, which increases
   the AUC of one sensor by known factors with the decoders, unmixes the mixtures, and analyzes
   how the difference from the ground truth depends on the AUC ratio of the two sensors.

## Data format

- **Spectral recordings (TXT).** OceanView export: a header ending with the line
  `>>>>>Begin Spectral Data<<<<<`, one row of wavelengths, then one row per time point (date,
  time, timestamp, then one intensity value per wavelength).
- **Initial emission spectra (CSV).** First column: wavelength (nm). The emission column is the
  first column whose header contains `em` and does not contain `ex` or `2p` (for example CSV
  files exported from FPbase, https://www.fpbase.org).
- **Folder layout and trial structure** for `validation` and `simulation`: see the README of
  each folder.

Spectral data and spectrum files are not included.

## Original NMF-RI code

`NMF_RI.m` and `dataPreprocessing.m` are from the original NMF-RI code
(https://github.com/djimenezsanchez/NMF-RI). `NMF_RI.m` is unmodified. `dataPreprocessing.m`
has an optional threshold argument added. See `NMF/README.md` for details.

## License

See `LICENSE`.

## Citation

If you use this code, please cite:

- Jiménez-Sánchez D, Ariz M, Morgado JM, Cortés-Domínguez I, Ortiz-de-Solórzano C.
  NMF-RI: blind spectral unmixing of highly mixed multispectral flow and image cytometry data.
  *Bioinformatics* 36(5):1590–1598 (2020). doi:10.1093/bioinformatics/btz751
- [This work: citation to be added]

## Contact

<!-- Authors: add the contact (e.g. the lead contact of the paper) before release. -->
