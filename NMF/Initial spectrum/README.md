# Initial emission spectra

Emission spectra used as the initial spectra (**A0**) of NMF-RI in `NMF_RI_Photometry.m`, and in
the validation and simulation scripts. Select one file per component in Tab 1 of
`NMF_RI_Photometry` (or pass it as `cfg.comp1Csv` / `cfg.comp2Csv` in `../../validation` and
`../../simulation`). The initial spectra only start the factorization; the final spectra are
estimated from the recording.

## Files

| File | Fluorophore | Source |
|---|---|---|
| `EGFP.csv` | EGFP | FPbase |
| `GCaMP6f.csv` | GCaMP6f | FPbase |
| `jRGECO1a.csv` | jRGECO1a | FPbase |
| `tdTomato.csv` | tdTomato | FPbase |

The spectrum of RdLight1 is not included (see below).

<!-- Authors: check that the file names above match the files in this folder. -->

## Sources

**FPbase.** The spectra of EGFP, GCaMP6f, jRGECO1a and tdTomato were downloaded from FPbase
(https://www.fpbase.org) on [date]. If you use these spectra, please cite FPbase and the original
publications of the fluorophores listed on FPbase.

- Lambert TJ. FPbase: a community-editable fluorescent protein database. *Nat Methods*
  16:277–278 (2019).

**RdLight1.** The spectrum of RdLight1 used in this study was provided by Lin Tian and is not
included in this repository. Please request it from Lin Tian.

- Patriarchi T, et al. An expanded palette of dopamine sensors for multiplex imaging in vivo.
  *Nat Methods* 17:1147–1155 (2020).

## Format

CSV with the wavelength (nm) in the first column. `NMF_RI_Photometry.m` uses the first column
whose header contains `em` and does not contain `ex` or `2p` (for example `EGFP em`) as the
emission spectrum. Each spectrum is baseline-subtracted, peak-normalized and interpolated onto
the wavelength axis of the recording; wavelengths outside the range of the file are set to zero.
Spectra from other sources can be used in the same format.
