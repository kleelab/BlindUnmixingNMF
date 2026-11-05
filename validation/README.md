# Validation of NMF-RI unmixing with in silico mixtures

Validation of blind NMF-RI unmixing against a known ground truth. Recordings from mice that
express a single fluorophore are summed in silico, the mixtures are unmixed with the same
NMF-RI pipeline as the experimental data (`../NMF`), and each unmixed component is compared
with the original recording before mixing (ground truth).

<!-- Authors: fill in the MATLAB version(s) used for testing before release. -->

## Requirements

- MATLAB R2020b or later (`uifigure`, `readtable` with `VariableNamingRule`, `writecell`)
- Statistics and Machine Learning Toolbox (`ecdf`, called in `dataPreprocessing.m`)
- `NMF_RI.m` and `dataPreprocessing.m` from the `NMF` folder on the MATLAB path

## Files

| File | Description |
|---|---|
| `NMF_Hybrid_Validation.m` | Builds the in silico mixtures, unmixes them with NMF-RI, compares unmixed and ground-truth spectra (R²) and AUCs, and saves figures, an Excel file and a `.mat` file |
| `NMF_Hybrid_Equivalence.m` | Agreement and equivalence statistics (CCC, Bland–Altman, proportional bias, TOST) and figures from the `.mat` file, without re-analysis |

## Quick start

1. Put the `NMF` and `validation` folders on the MATLAB path.
2. Run `NMF_Hybrid_Validation`. In the input window, select
   - the green slot (sensor folder or control dye file),
   - the red slot (sensor folder or control dye file),
   - one emission spectrum CSV per slot (initial spectra, as in `../NMF`),
   - the output folder,

   and press *Run*.
3. Run `NMF_Hybrid_Equivalence` and select the `.mat` file saved in step 2
   (`<green>+<red>_mixture_R2_v2.mat`).

Without the input window:

```matlab
cfg.greenPath = 'D:\data\GCaMP';            % sensor folder (or dye file)
cfg.redPath   = 'D:\data\jRGECO';           % sensor folder (or dye file)
cfg.comp1Csv  = 'D:\spectra\GCaMP6f.csv';   % initial spectrum, green slot
cfg.comp2Csv  = 'D:\spectra\jRGECO1a.csv';  % initial spectrum, red slot
cfg.outDir    = 'D:\analysis';
cfg.pairing   = 'all';                      % every green mouse with every red mouse
NMF_Hybrid_Validation(cfg);
NMF_Hybrid_Equivalence('D:\analysis\GCaMP+jRGECO_mixture_R2_v2.mat');
```

Optional fields: `cfg.greenName`, `cfg.redName` (default: folder or file name),
`cfg.wlRange` (default `[450 650]`), `cfg.theta1`, `cfg.theta2`, `cfg.alphaA`, `cfg.alphaX`
(defaults 0.01, 0.01, 0.001, 0.001).

## Input data

**Sensor (folder).** One subfolder per mouse. Each subfolder contains a spectral time series
(TXT, OceanView format as in `../NMF`) whose file name contains the slot tag: `G` for the green
slot (470 nm excitation) and `R` for the red slot (560 nm excitation). The tag must follow a
separator (`_`, `.`, `-`, space) or start the name and must not be followed by a letter, so `G`
does not match `GR`; upper and lower case are treated as the same. If the folder has no mouse
subfolders, TXT files with the tag directly in the folder are used, and the text before the tag
is taken as the mouse name.

**Control dye (file).** One TXT file of a fluorescent protein that does not report activity
(for example tdTomato or eGFP), either in the OceanView format or as two columns (wavelength,
intensity). All rows are averaged into one static spectrum.

**Trial structure.** Recordings are split into trials of 1,177 samples (100 Hz). Reward onset is
sample 400 (0 s). The reward period is samples 400 to 499 (0 to 1 s), and the baseline period is
samples 100 to 150 (−3 to −2.5 s). Incomplete trials at the end of a recording are not used.
These values are set in `defaultSettings` of `NMF_Hybrid_Validation.m`.

At least one slot must be a sensor folder. Modes: sensor + sensor, sensor + dye, dye + sensor.

## Method

**Ground truth.** For each recording, the offset of the file (the minimum of its time-averaged
spectrum, 450 to 650 nm) is subtracted. For a dye file, the minimum of the averaged spectrum is
subtracted. The offset-corrected recordings are the ground truth.

**Mixtures.**
- Sensor + sensor: trial k of the green recording is added to trial k of the red recording, over
  the smaller number of trials of the two. `pairing = 'one-to-one'` (default) pairs mice sorted
  by name once each; `pairing = 'all'` combines every green mouse with every red mouse.
- Sensor + dye: the static dye spectrum is added to every time point of each sensor recording
  (one mixture per sensor mouse).

**Unmixing.** Each mixture is unmixed with NMF-RI exactly as in `NMF_RI_Photometry.m`, with the
same initial spectra and parameters for every mixture. Component 1 is the green slot and
component 2 is the red slot.

**Spectra.** For each trial, the reward spectrum is the mean over the reward period, and the
delta spectrum is the reward spectrum minus the mean over the baseline period. The similarity of
the unmixed and ground-truth spectra is quantified for each trial as
R² = 1 − SS_res / SS_tot over 450 to 650 nm (ground truth = observed, unmixed = predicted), for
the reward and delta spectra, and is reported as mean ± SEM over trials.

**AUC.** AUC is measured on the trial-mean reward spectrum around the peak of each spectrum
(peak ± 10 nm, trapezoidal integration). Peaks are found after 2 nm smoothing, T1 within 500 to
540 nm (green sensor) and T2 within 580 to 620 nm (red sensor). Ground-truth and unmixed peaks
are found separately, so each spectrum is measured around its own peak. The T2 AUC of the green
sensor uses the T2 peak of the red slot. ΔAUC is the AUC of the delta spectrum in the same window.

**Agreement and equivalence (`NMF_Hybrid_Equivalence.m`).** Percentage difference =
(unmixed − ground truth) / ground truth × 100. Because each mouse contributes to several
mixtures, all statistics use mouse-level values (the mean over the mixtures that include the
mouse).
- Lin's concordance correlation coefficient (CCC) of the AUC, 95% CI by Fisher z (Lin, 1989)
- Bland–Altman bias and 95% limits of agreement, bias ± 1.96 SD (Bland & Altman, 1986)
- Proportional bias: least-squares regression of the difference on the ground-truth AUC, slope
  tested against 0 (t test, df = n − 2) (Bland & Altman, 1999)
- Equivalence: two one-sided tests (TOST) against ±10%, α = 0.05 for each test, equivalent when
  the 90% CI of the mean difference lies within the margin (Schuirmann, 1987); AUC and ΔAUC
- Trial-to-trial coefficient of variation of the ground-truth AUC of each mouse

`NMF_Hybrid_Equivalence.m` does not require the Statistics and Machine Learning Toolbox.

## Outputs

`NMF_Hybrid_Validation.m` (in the output folder):

| File | Content |
|---|---|
| `mixtures_v2/Mixture##_<label>_mixture_v2.pdf` | Ground truth (line) and unmixed components (filled) of each mixture, reward and delta spectra |
| `<green>+<red>_mixture_summary_v2.pdf` | Mean over mixtures (normalized), SEM band, R² |
| `<green>+<red>_mixture_AUC_peak_v2.pdf` | Ground truth vs unmixed AUC, ΔAUC and peak wavelength, one line per mixture |
| `<green>+<red>_mixture_R2_v2.xlsx` | Sheets README, R2, R2_trials, AUC_peak, Summary, Info |
| `<green>+<red>_mixture_R2_v2.mat` | Results (variable `R`) for `NMF_Hybrid_Equivalence.m` |

Unmixed spectra are not saved as files.

`NMF_Hybrid_Equivalence.m` (in the folder of the `.mat` file; figures as PDF and EPS):

| File | Content |
|---|---|
| `*_Agreement` | Ground truth vs unmixed AUC with the line of identity and CCC |
| `*_BlandAltman` | Percentage difference vs ground-truth AUC, bias, limits of agreement and proportional bias |
| `*_TOST` | Mean difference and 90% CI of AUC and ΔAUC with the equivalence margin |
| `*_Equivalence.xlsx` | Sheets README, Agreement, BlandAltman, ProportionalBias, TOST, MouseMeans, TrialCV |

Figure settings (axes, marker size, colors, margin) are in the SETTINGS block at the top of
`NMF_Hybrid_Equivalence.m`. Edit them and run again to redraw without re-analysis.

## Notes

- `dataPreprocessing.m` subtracts the minimum across wavelengths from each time point before
  NMF, as in the experimental pipeline. This offset is not part of the unmixed components. Its
  mean is reported in the Info sheet of the Excel file (Preprocessing offset removed (counts)).
- In sensor + dye mixtures, the ΔAUC of the dye is zero by construction. The ΔAUC assigned to
  the dye component is reported in the Excel file as leakage,
  100 × AUC(unmixed dye delta) / AUC(ground-truth sensor delta).

## References

- Bland JM, Altman DG. Statistical methods for assessing agreement between two methods of
  clinical measurement. *Lancet* 327:307–310 (1986).
- Bland JM, Altman DG. Measuring agreement in method comparison studies.
  *Stat Methods Med Res* 8:135–160 (1999).
- Lin LI. A concordance correlation coefficient to evaluate reproducibility.
  *Biometrics* 45:255–268 (1989).
- Schuirmann DJ. A comparison of the two one-sided tests procedure and the power approach for
  assessing the equivalence of average bioavailability.
  *J Pharmacokinet Biopharm* 15:657–680 (1987).

## Citation

- [This work: citation to be added]
