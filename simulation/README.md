# Unmixing at different AUC ratios (simulation)

Tests how the accuracy of blind NMF-RI unmixing depends on the AUC ratio of two sensors. The
total AUC of one sensor is increased by known factors with a power decoder built from spectra
measured at several excitation powers, the scaled recording is mixed in silico with the other
sensor at its measured AUC, and each mixture is unmixed with the same NMF-RI pipeline as the
experimental data (`../NMF`). Each unmixed component is compared with its ground truth.

<!-- Authors: fill in the MATLAB version(s) used for testing before release. -->

## Requirements

- MATLAB R2020b or later (`readtable` with `VariableNamingRule`, `writecell`)
- Statistics and Machine Learning Toolbox (`ecdf`, called in `dataPreprocessing.m`)
- `NMF_RI.m` and `dataPreprocessing.m` from the `NMF` folder on the MATLAB path

## Files

| File | Description |
|---|---|
| `NMF_Power_Decoder.m` | Builds one power decoder per mouse from spectra measured at several excitation powers, validates it by cross-validation, and saves it as a `.mat` file (one file per sensor) |
| `NMF_Spike_Recovery.m` | Scales the AUC of one sensor with the decoders, creates the in silico mixtures, unmixes them with NMF-RI, and analyzes the agreement and the dependence of the difference on the AUC ratio |

## Quick start

1. Put the `NMF` and `simulation` folders on the MATLAB path.
2. Run `NMF_Power_Decoder` once for each sensor (for example GCaMP and jRGECO). Select the sensor
   folder with the power subfolders and the output folder. The decoder is saved as
   `PowerDecoder_<sensor>.mat`.
3. Run `NMF_Spike_Recovery`. Select the green sensor folder, the red sensor folder, the two
   emission spectrum CSV files (initial spectra, as in `../NMF`), the output folder, the two names,
   and the two decoder files from step 2.

Without dialogs:

```matlab
NMF_Power_Decoder(struct('sensorDir', 'D:\power\GCaMP',  'outDir', 'D:\analysis'));
NMF_Power_Decoder(struct('sensorDir', 'D:\power\jRGECO', 'outDir', 'D:\analysis'));

cfg.greenPath    = 'D:\data\GCaMP';              % mouse subfolders with G files
cfg.redPath      = 'D:\data\jRGECO';             % mouse subfolders with R files
cfg.comp1Csv     = 'D:\spectra\GCaMP6f.csv';     % initial spectrum, green sensor
cfg.comp2Csv     = 'D:\spectra\jRGECO1a.csv';    % initial spectrum, red sensor
cfg.outDir       = 'D:\analysis';
cfg.greenDecoder = 'D:\analysis\PowerDecoder_GCaMP.mat';
cfg.redDecoder   = 'D:\analysis\PowerDecoder_jRGECO.mat';
NMF_Spike_Recovery(cfg);
```

Optional fields: `cfg.sensorName` (decoder; default: folder name), `cfg.greenName`,
`cfg.redName` (default: folder names), `cfg.wlRange` (default `[450 650]`), `cfg.theta1`,
`cfg.theta2`, `cfg.alphaA`, `cfg.alphaX` (defaults 0.01, 0.01, 0.001, 0.001).

Redraw figures and Excel files from saved results, without re-analysis:

```matlab
NMF_Power_Decoder('D:\analysis\PowerDecoder_GCaMP.mat');
NMF_Spike_Recovery('D:\analysis\GCaMP+jRGECO_SpikeRecovery.mat');
```

Figure and analysis settings are in the SETTINGS block at the top of each file.

## Input data

**Power decoder (one folder per sensor).**

```
<sensor folder>/power25/<mouse>_<...>.txt
<sensor folder>/power50/<mouse>_<...>.txt
<sensor folder>/power100/<mouse>_<...>.txt
```

The power is the number in the subfolder name (`power25`, `Power25` and `POWER25` give 25 µW).
The mouse name is the text before the first `_` of the file name; upper and lower case are
treated as the same. Each file is a static spectrum, either an OceanView time series (one row
per frame) or wavelength and value columns (several frames allowed). A folder for the reference
power (`P.refPower`, 25 µW, the power of the time-series recordings) is required. The mouse
names must match those of the sensor recordings used in `NMF_Spike_Recovery.m`.

**Sensor recordings.** The same layout as in `../validation`: one subfolder per mouse, each
with a spectral time series (TXT) whose file name contains `G` (green sensor, 470 nm excitation)
or `R` (red sensor, 560 nm excitation). Trials of 1,177 samples (100 Hz), reward onset at sample
400, reward period samples 400 to 499 (0 to 1 s), baseline samples 100 to 150 (−3 to −2.5 s).

## Method

**Power decoder (`NMF_Power_Decoder.m`).** For each mouse and power p, the frames of the file are
averaged, and the minimum of the mean spectrum is subtracted (offset). The intensity ratio k_p is
the total AUC (450 to 650 nm) at p divided by that at the reference power. At each wavelength λ,
the spectra at all powers are fitted by least squares with a straight line in k,

S(λ, k) = B(λ) + k E(λ),

where E(λ) changes in proportion to the total AUC and B(λ) does not change with it. To increase
the total AUC of a recording at the reference power by a factor k, the spectrum at every time
point is multiplied by the gain

G_k(λ) = [B(λ) + k E(λ)] / [B(λ) + E(λ)].

Where the model spectrum at the reference power, B(λ) + E(λ), is below 2% of its maximum, the
gain is k. The gain is smoothed with a 2 nm moving average. A mouse is excluded when its total
AUC does not rise by at least 1.1 times from one power to the next (quality check), or when it
has no spectrum at the reference power or only one power.

**Decoder validation.** For every mouse and every power, the spectrum at that power is predicted
without it (cross-validation) and compared with the measured spectrum over 480 to 630 nm.
- Decoder (animal): the line through the other powers of the same mouse
- Decoder (pooled): the fixed fraction B / (B + E) averaged over the other mice, applied to the
  nearest other power of the mouse
- Proportional: the nearest other power multiplied by the ratio of the total AUCs (baseline)
- Noise ceiling: half the difference between the means of the odd and even frames

Decoding accuracy = 1 − ‖predicted − measured‖ / ‖measured‖, and normalized RMSE (%) =
100 × ‖predicted − measured‖ / ‖measured‖. Mean ± SEM over mice. The methods are compared at
every power with a paired t test and an exact Wilcoxon signed rank test (n = mice).

**Mixtures (`NMF_Spike_Recovery.m`).** The ground truth of each recording is defined as in
`../validation` (offset of each file subtracted). Every green mouse is paired with every red mouse
(`P.pairing = 'all'`). For each pair, one sensor is scaled by the factors in `P.factors`
(default 1, 1.5, 2, 2.5 and 3 times the measured total AUC) while the other sensor stays at its
measured level (1×):
- A. green sensor varied, red sensor 1×
- B. red sensor varied, green sensor 1×

The 1× mixture is shared by A and B. With `P.scaling = 'decoder'` (default), every time point is
multiplied by the gain G_k(λ) of the mouse's decoder. With `P.scaling = 'proportional'`, the
whole recording is multiplied by k. The scaled ground truth is mixed trial by trial with the
other sensor, and each mixture is unmixed with NMF-RI as in `NMF_RI_Photometry.m`. Factors
outside the measured power range of a decoder are extrapolated and reported as a warning. Mice
without a decoder are skipped.

**Measures.** R² of the unmixed and ground-truth spectra and the AUC (T1 for the green sensor,
T2 for the red sensor, own peak ± 10 nm) are computed as in `../validation`. Percentage
difference = (unmixed − ground truth) / ground truth × 100.
- Agreement and equivalence for every condition, on mouse-level means: Lin's CCC (95% CI),
  Bland–Altman bias and 95% limits of agreement, proportional bias regression, TOST (±10%)
- Difference vs input ratio: x = total ground-truth AUC of the unchanged sensor (1×) / total
  ground-truth AUC of the varied sensor (n×), over all mixtures and factors. Least-squares fit of
  difference = c + b·x for the varied sensor and difference = c + b/x for the unchanged sensor.
  The shaded band is the standard error of the fitted mean (residual variance of the fit). The
  95% CI of b and c is estimated by bootstrap resampling of green and red mice independently
  (`P.nBoot`, default 2,000).
- Window scaling: ratio of the T1 or T2 AUC (n× / 1×) vs ratio of the total AUC (n× / 1×) of the
  time-averaged ground-truth spectra, with the slope and R² of the least-squares line over all
  points. A slope of 1 means that a fold change defined by the total AUC applies equally to the
  AUC within each window.

## Outputs

`NMF_Power_Decoder.m` (in the output folder):

| File | Content |
|---|---|
| `PowerDecoder_<sensor>.mat` | Decoder of every mouse (variable `D`), input for `NMF_Spike_Recovery.m` |
| `PowerDecoder_<sensor>_animals.pdf/.eps` | Relative gain (G_k / k) of every mouse at the measured powers |
| `PowerDecoder_<sensor>_summary.pdf/.eps` | Relative gain at the highest power, T2/T1 and peak shifts vs intensity ratio |
| `PowerDecoder_<sensor>_decoding_validation.pdf/.eps` | Decoding accuracy and normalized RMSE vs held-out power (mean ± SEM) |
| `PowerDecoder_<sensor>.xlsx` | Sheets README, Animals, CrossValidation, CVSummary, CVTests, RelGain |

`NMF_Spike_Recovery.m` (in the output folder):

| File | Content |
|---|---|
| `<green>+<red>_SpikeRecovery_spectra_<sensor>.pdf/.eps` | Ground-truth spectra of every mouse at every factor |
| `<green>+<red>_SpikeRecovery_AUC_<sensor>_varied.pdf/.eps` | Unmixed vs ground-truth AUC and difference vs ground-truth AUC at every factor, for the varied and the unchanged sensor |
| `<green>+<red>_SpikeRecovery_DiffRatio_<sensor>_varied.pdf/.eps` | Difference vs input ratio with the fitted curve |
| `<green>+<red>_SpikeRecovery_WindowScaling_<sensor>.pdf/.eps` | T1 and T2 AUC ratio vs total AUC ratio |
| `spectrum/` | Mixture, ground-truth and unmixed spectra of every mixture and condition, and `<green>+<red>_Spectrum_R2.xlsx` |
| `unmixing/` | Agreement and Bland–Altman figures of every condition, and `<green>+<red>_Unmixing_stats.xlsx` (README, Agreement, BlandAltman, ProportionalBias, TOST, MouseMeans) |
| `<green>+<red>_SpikeRecovery.xlsx` | Sheets README, Pairs, Mice, Levels, Stats, AUC, AUCStats, DiffRatio, DiffRatioFit, WindowScaling, WindowScalingFit, WindowScalingStats |
| `<green>+<red>_SpikeRecovery.mat` | All results (variable `R`), for redrawing |

Unmixed spectra are not saved as data files. With many mixtures, the `spectrum/` figures take
most of the run time; set `P.supSpectra = false` to skip them.

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
