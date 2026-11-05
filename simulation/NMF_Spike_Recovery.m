function varargout = NMF_Spike_Recovery(cfg)
% NMF_Spike_Recovery
% Spike-and-recovery / linearity-of-dilution test of blind NMF-RI unmixing, done
% separately for every animal (in silico mixtures of measured single-sensor data).
%
% Every green mouse is paired with every red mouse at their MEASURED brightness
% (no reference value). One sensor is then made k times brighter (brightness ratio k =
% total AUC / total AUC at the recording power) while the other stays at its measured level:
%     A. green varied: green at k = P.factors, red x 1
%     B. red varied:   red at k = P.factors,   green x 1
% (k = 1: the native mixture, shared by A and B and identical to the mixtures of
% NMF_Hybrid_Equivalence). How the recording is made k times brighter (P.scaling):
%   'decoder' (default): every wavelength (every time point) x the gain G_k(lambda) of that
%       animal's power decoder (NMF_Power_Decoder; measured at 25, 50, 100 uW), so that the
%       spectrum changes as measured when the brightness changes (not simply scaled)
%   'proportional': the whole recording x k (spectrum shape unchanged)
% Each mixture is unmixed blindly by NMF-RI exactly as NMF_RI_Photometry; AUC with v2
% windows (own peak +/- 10 nm).
%
% Per animal (the animal whose signal was varied; mean over its mixtures):
%   recovery line   unmixed AUC vs ground-truth AUC (absolute, counts nm) over the
%                   factors: slope (1 = proportional) and intercept (constant amount)
%   recovery (%)    (unmixed(f) - unmixed(x1)) / (GT(f) - GT(x1)) x 100: how much of the
%                   added / removed amount is found (100% = complete)
%   difference (%)  (unmixed - GT) / GT x 100 at every factor
%   crosstalk (%)   slope of the partner's unmixed AUC vs the varied sensor's GT AUC
%                   x 100: % of the added amount that appears in the other sensor
%                   (partner GT is constant)
%
% Run once per red sensor (e.g. GCaMP + jRGECO, then GCaMP + RdLight); compare the
% two .mat files with NMF_Spike_Recovery_Compare.
%
% Usage
%   NMF_Spike_Recovery              % dialogs: green folder, red folder, two CSVs, output folder,
%                                   % then the two decoder files (PowerDecoder_<green>.mat, _<red>.mat)
%   cfg.greenDecoder, cfg.redDecoder  % decoder files (P.scaling = 'decoder')
%   NMF_Spike_Recovery(cfg)         % cfg.greenPath, cfg.redPath, cfg.comp1Csv, cfg.comp2Csv,
%                                   % cfg.outDir, optional cfg.greenName, cfg.redName
%   NMF_Spike_Recovery('...\GCaMP+jRGECO_SpikeRecovery.mat')   % redraw figure and Excel
%                                   % from a saved result (no NMF-RI, takes seconds)
%
% Output (cfg.outDir), one set per sensor (green, red)
%   <green>+<red>_SpikeRecovery_spectra_<sensor>.pdf / .eps
%       A. every animal: time-averaged ground-truth spectra at every factor, overlaid
%          (light = x1 ... dark = largest factor), 450-650 nm
%   <green>+<red>_SpikeRecovery_AUC_<sensor>_varied.pdf / .eps   (one file per varied sensor)
%       title '<varied> x1 - x3, <other> unchanged'
%       row 1: the varied (scaled) sensor; row 2: the other (unscaled) sensor, from the same mixtures
%       B. unmixed vs ground-truth T1 (green) / T2 (red) AUC: line per factor through the animal means,
%          colour = factor of the varied sensor (light x1 ... dark); 1:1 dashed line; CCC [95% CI] per factor
%       C. difference (%) vs ground-truth T1 / T2 AUC: regression line per factor with R^2
%   <green>+<red>_SpikeRecovery_DiffRatio_<sensor>_varied.pdf / .eps   (one file per varied sensor)
%       difference (%) (T1 / T2 AUC) of the varied sensor (left) and of the unchanged sensor (right) vs the
%       input ratio x = total GT AUC of the unchanged sensor (1x) / of the varied sensor (nx) (reward spectra,
%       all wavelengths); light dots = every mixture and factor (not grouped); black: least squares fit over all
%       points (varied: c + b x; unchanged: c + b / x), grey band = +/- SEM of the fitted curve; b = difference
%       added at x = 1, c = difference without the partner (95% CI by animal bootstrap in DiffRatioFit)
%   <green>+<red>_SpikeRecovery_WindowScaling_<sensor>.pdf / .eps   (one file per sensor)
%       T1 and T2 AUC ratio (x n / x1) vs total AUC ratio of the ground-truth spectra (time average): light dots =
%       every animal and factor (not grouped), dark dots = mean per factor, dashed y = x; slope and R^2 of
%       the least squares line over all points
%   spectrum\<green mouse>_<red mouse>_<green>x<k>_<red>x<k>.pdf   (every mixture and condition)
%       left: mixture reward spectrum (counts), black = ground truth, fills = NMF-RI components (stacked),
%             NMF unmixing R^2; middle / right: ground truth (black) and unmixed (fill) of each sensor,
%             normalised to the ground-truth peak, R^2 (per trial, mean +/- SEM)
%   spectrum\<green>+<red>_Spectrum_R2.xlsx    R^2 of every mixture and condition, summary per condition
%   unmixing\<green>x<k>_<red>x<k>.pdf / .eps   (every condition; as NMF_Hybrid_Equivalence); axes the same
%       within each series (green varied / red varied) and as in _AUC_<sensor>_varied; x1 + x1 is drawn
%       for both series: <green>x1_<red>x1_<green>varied and ..._<red>varied
%       green agreement, green Bland-Altman, red agreement, red Bland-Altman (light dots = mixtures,
%       dark dots = mouse means; CCC [95% CI]; bias, 95% LoA, proportional-bias regression)
%   unmixing\<green>+<red>_Unmixing_stats.xlsx  README, Agreement, BlandAltman, ProportionalBias, TOST, MouseMeans
%   <green>+<red>_SpikeRecovery.xlsx   README, Pairs, Mice, Levels, Stats, AUC, AUCStats, DiffRatio,
%                                      DiffRatioFit, WindowScaling, WindowScalingFit, WindowScalingStats
%   <green>+<red>_SpikeRecovery.mat    all results (redraw: NMF_Spike_Recovery('...mat'))
%
% Required on the MATLAB path: NMF_RI.m, dataPreprocessing.m

%% ===================== SETTINGS =====================
P.factors        = [1 1.5 2 2.5 3];   % brightness ratios (must contain 1); with the decoder keep them within the
                                      % measured range of every animal (25-100 uW: about 1 to 3.1-4.1)
P.scaling        = 'decoder';     % 'decoder': every wavelength x the animal's power-decoder gain G_k(lambda)
                                  % (NMF_Power_Decoder); 'proportional': whole spectrum x k
P.pairing        = 'all';         % 'all' (every green with every red mouse) or 'one-to-one'
P.margin         = 10;            % +/- % reference lines and TOST margin of the recovery (100 +/- margin)
P.slopeMargin    = 0.10;          % TOST margin of the slope (1 +/- slopeMargin)
P.fontName       = 'Arial';
P.fontSize       = 8;
P.panelCm        = [4.2 4.2];     % panel width, height (cm)
P.dotSize        = 8;             % marker size (pt), no outline
P.lineWidth      = 0.8;           % pt, line of each animal
P.tintMax        = 0.6;           % animals drawn from the sensor colour (first) to 60% white (last)
P.refColor       = [0 0 0];       % identity / reference lines
% ---- axes of the AUC figures and of unmixing\ (the same within each series)
%   rows: varied sensor (1 = green varied, 2 = red varied); columns: sensor shown (1 = green T1, 2 = red T2)
%   P.axAUC  [lo hi tick], x 10^3 counts nm: agreement x and y, Bland-Altman x, panels B and C
%   P.axDiff [lo hi tick], %: difference axis, panel C and Bland-Altman y
%   [] = automatic (covers every mixture of the series). Values below: GCaMP + jRGECO (main figure)
P.axAUC  = {[0 150 50], [0 15 5];      % green varied: GCaMP T1, jRGECO T2
            [0 40 10],  [0 50 10]};    % red varied:   GCaMP T1, jRGECO T2
P.axDiff = {[0 1 0.2],  [0 12 2];
            [0 2.5 0.5], [0 2.5 0.5]};
% ---- difference vs input ratio (_DiffRatio_<sensor>_varied; every mixture and factor)
%   x = input ratio: total GT AUC of the unchanged sensor (1x) / total GT AUC of the varied sensor (nx);
%   y = difference (%) of the varied sensor (left panel, fit c + b x) and of the unchanged sensor (right, c + b / x)
%   [lo hi tick]; ticks at lo, hi and multiples of 'tick' (e.g. [-0.5 4 1]: -0.5, 0, 1, 2, 3, 4)
P.axRatioX    = {[], [0 6 1]};    % x per varied sensor (1 = green, 2 = red); [] = automatic
P.axRatioY    = {[], []; [-0.5 4 1], []}; % y (%); rows: varied sensor (1 = green, 2 = red),
                                  % columns: sensor shown (1 = green T1, 2 = red T2); [] = automatic
P.ratioDotSize = 5;               % marker size (pt) of the mixture dots
P.nBoot       = 2000;             % bootstrap resamples (over animals) for the 95% CI of b and c
P.axisLineWidth  = 0.6;           % pt
% ---- supplementary figures (subfolders of the output folder)
%   spectrum\  <green mouse>_<red mouse>_<green>x<k>_<red>x<k>   every mixture and condition
%   unmixing\  <green>x<k>_<red>x<k>                             every condition (all mixtures)
P.supSpectra       = true;          % spectrum\ figures (one per mixture and condition)
P.supSpectrumTypes = {'.pdf'};      % file types of the spectrum figures; {'.pdf', '.eps'} for both
P.supPanelCm       = [4.0 4.0];     % panel width, height (cm)
P.specXLim         = [450 650];     % wavelength axis (nm)
P.specXTick        = 50;
P.specLineWidth    = 1.0;           % ground-truth line (pt)
P.specYLim         = [-50 5000];    % mixture panel y axis (counts), the same in every file; lower limit < 0 leaves a
                                    % small gap above the x axis. [] = automatic, common to all files (covers every mixture)
P.specYTick        = 1000;          % mixture panel y tick step (counts)
P.pairTint         = 0.7;           % unmixing\: mixture dots = sensor colour mixed with white (dark dots = mouse means)
P.supRefColor      = [0.55 0.55 0.55];   % identity and zero lines
P.baMinSpan        = 3;             % Bland-Altman automatic y range: minimum span (percentage points)
P.baRegLine        = 'all';         % proportional-bias regression line: 'all', 'significant' or 'none'
P.alpha            = 0.05;          % TOST: alpha of each one-sided test (90% CI); margin = P.margin
%% ====================================================

    opt = defaultSettings();
    if nargin >= 1 && ischar(cfg) && numel(cfg) > 4 && strcmpi(cfg(end-3:end), '.mat')
        % redraw only
        S = load(cfg);
        R = S.R;
        [d, n] = fileparts(cfg);
        base = fullfile(d, n);
        fprintf('Redrawing %s\n', cfg);
        MS = spikeMetrics(R);
        ST = spikeStats(MS, P);
        spikeFigures(R, MS, P, base);
        writeSpikeExcel(R, MS, ST, P, [base '.xlsx']);
        printSummary(R, MS, ST);
        if nargout > 0, varargout{1} = R; end
        return;
    end
    if nargin < 1 || isempty(cfg)
        cfg = askInputs();
        if isempty(cfg), fprintf('Cancelled.\n'); return; end
    end
    cfg = fillCfg(cfg, opt);
    opt.wlRange = cfg.wlRange;
    if exist('NMF_RI', 'file') ~= 2 || exist('dataPreprocessing', 'file') ~= 2
        error('NMF_RI.m and dataPreprocessing.m must be on the MATLAB path.');
    end
    fac = P.factors(:)';
    if ~any(abs(fac - 1) < 1e-9), error('P.factors must contain 1.'); end
    useDec = strcmpi(P.scaling, 'decoder');
    DEC = {[], []};
    if useDec
        fld = {'greenDecoder', 'redDecoder'}; lbl = {'GREEN', 'RED'};
        for s = 1:2
            if ~isfield(cfg, fld{s}) || isempty(cfg.(fld{s}))
                [f, d] = uigetfile({'PowerDecoder_*.mat', 'Power decoder (*.mat)'}, ...
                    sprintf('Power decoder of the %s sensor (PowerDecoder_....mat)', lbl{s}), cfg.outDir);
                if ~ischar(f), fprintf('Cancelled.\n'); return; end
                cfg.(fld{s}) = fullfile(d, f);
            end
            S = load(cfg.(fld{s}));
            if ~isfield(S, 'D') || ~isfield(S.D, 'animals') || ~isfield(S.D.animals, 'fixed')
                error('%s is not a power decoder of NMF_Power_Decoder 2.0.', cfg.(fld{s}));
            end
            DEC{s} = S.D;
            fprintf('Decoder %s: %s (%d animals)\n', lbl{s}, cfg.(fld{s}), numel(S.D.animals));
        end
    end
    if exist(cfg.outDir, 'dir') ~= 7, mkdir(cfg.outDir); end
    names = {cfg.greenName, cfg.redName};
    cols = {sensorColor(names{1}, opt, opt.colorComp1), sensorColor(names{2}, opt, opt.colorComp2)};
    % conditions: [green factor, red factor, direction (1 = green varied, 2 = red varied), level index]
    nF = numel(fac);
    cond = [fac', ones(nF, 1), ones(nF, 1), (1:nF)'; ones(nF, 1), fac', 2 * ones(nF, 1), (1:nF)'];
    % unique conditions to compute (x 1 is shared by A and B)
    cr = round(cond(:, 1:2) * 1e6);
    uc = zeros(0, 2); cIdx = zeros(size(cr, 1), 1);
    for i = 1:size(cr, 1)
        j = find(round(uc(:, 1) * 1e6) == cr(i, 1) & round(uc(:, 2) * 1e6) == cr(i, 2), 1);
        if isempty(j), uc(end+1, :) = cond(i, 1:2); j = size(uc, 1); end %#ok<AGROW>
        cIdx(i) = j;
    end
    nU = size(uc, 1);
    tTotal = tic;
    warn = {};

    fprintf('\n==============================================================\n');
    fprintf(' NMF-RI spike-and-recovery: %s + %s\n', names{1}, names{2});
    fprintf(' factors of the measured AUC: %s\n', num2str(fac));
    fprintf(' A. %s varied (%s x 1)   B. %s varied (%s x 1)\n', names{1}, names{2}, names{2}, names{1});
    fprintf('==============================================================\n');

    % ----- mice
    tags = {opt.greenTag, opt.redTag};
    paths = {cfg.greenPath, cfg.redPath};
    src = cell(1, 2);
    for s = 1:2
        [src{s}, w] = listAnimals(paths{s}, tags{s}); warn = [warn, w]; %#ok<AGROW>
        fprintf('%s mice (%s file): %d\n', names{s}, tags{s}, numel(src{s}));
        if isempty(src{s}), error('No %s mice with a %s file found in %s.', names{s}, tags{s}, paths{s}); end
    end
    nG = numel(src{1}); nR = numel(src{2});
    % ----- decoder animal of every mouse (checked before the files are read)
    if useDec
        % decoders chosen in the wrong order (green / red)?
        sn = {'', ''};
        for s = 1:2, if isfield(DEC{s}, 'sensor'), sn{s} = lower(DEC{s}.sensor); end, end
        if ~isempty(sn{1}) && ~isempty(sn{2}) && ~isempty(strfind(sn{1}, lower(names{2}))) && ...
                ~isempty(strfind(sn{2}, lower(names{1}))) %#ok<STREMP>
            DEC = DEC([2 1]);
            fprintf('Decoders were given in the reverse order; swapped (%s: %s, %s: %s).\n', ...
                names{1}, DEC{1}.sensor, names{2}, DEC{2}.sensor);
        end
        for s = 1:2
            dn = {DEC{s}.animals.name};
            nMatch = 0;
            fprintf('%s decoder animals: %s\n', names{s}, strjoin(dn, ', '));
            for k = 1:numel(src{s})
                ia = matchAnimal(src{s}(k).name, dn);
                src{s}(k).dec = ia;
                if ia > 0
                    nMatch = nMatch + 1;
                    if ~strcmp(src{s}(k).name, dn{ia})
                        fprintf('   %s -> decoder %s\n', src{s}(k).name, dn{ia});
                    end
                else
                    fprintf('   %s: no decoder (skipped)\n', src{s}(k).name);
                end
            end
            if nMatch == 0
                error(['None of the %s mice (%s) has a decoder in the %s decoder file (%s).\n' ...
                    'Check that the decoder file belongs to %s and that the animal names match.'], names{s}, ...
                    strjoin({src{s}.name}, ', '), names{s}, strjoin(dn, ', '), names{s});
            end
        end
    end
    if strcmp(P.pairing, 'all')
        [jj, ii] = meshgrid(1:nR, 1:nG); ii = ii'; jj = jj'; ii = ii(:); jj = jj(:);
    else
        n = min(nG, nR); ii = (1:n)'; jj = (1:n)';
    end
    nPair = numel(ii);
    nRun = nPair * nU;
    fprintf('\nTOTAL: %d mixtures x %d conditions = %d NMF-RI runs\n', nPair, nU, nRun);

    hw = [];
    try
        hw = waitbar(0, 'Reading files', 'Name', sprintf('Spike recovery: %d runs', nRun));
    catch
    end
    barCleaner = onCleanup(@() closeBar(hw)); %#ok<NASGU>

    % ----- read files (offset-free ground truth, as v2)
    fprintf('\nReading %d files\n', nG + nR);
    lab = {'green', 'red'};
    for s = 1:2
        for k = 1:numel(src{s})
            progress(hw, 0, sprintf('Reading %s', src{s}(k).name));
            [src{s}(k), w] = loadAnimal(src{s}(k), opt, lab{s}); warn = [warn, w]; %#ok<AGROW>
        end
    end

    % ----- mixtures
    L = opt.trialLength; rw = opt.rewardSamples; bl = opt.baselineSamples;
    sm = @(wl, y) smoothNm(wl, y, opt.peakSmoothNm);
    pairs = struct('green', {}, 'red', {}, 'nTrials', {}, 'gAUC0', {}, 'rAUC0', {}, 'cond', {}, ...
        'raw', {}, 'nmf', {}, 'rawDelta', {}, 'nmfDelta', {}, 'iter', {}, 'extrap', {}, 'wl', {}, 'spec', {});
    run = 0; tRun = tic;
    SPEC = struct('dir', {}, 'animal', {}, 'wl', {}, 'S', {});
    for p = 1:nPair
        g = src{1}(ii(p)); r = src{2}(jj(p));
        ttl = sprintf('%s %s + %s %s', names{1}, g.name, names{2}, r.name);
        if g.nTrials < 1 || r.nTrials < 1
            warn{end+1} = sprintf('%s skipped (no complete trials).', ttl); %#ok<AGROW>
            run = run + nU;
            continue;
        end
        nTr = min(g.nTrials, r.nTrials); nS = nTr * L;
        wl = g.wl(:);
        Gm = sourceMatrix(g, wl, nS);
        [Rm, interp] = sourceMatrix(r, wl, nS);
        if interp, warn{end+1} = sprintf('%s: wavelength axis differs; interpolated.', ttl); end %#ok<AGROW>
        A0 = [loadEmissionSpectrum(cfg.comp1Csv, wl), loadEmissionSpectrum(cfg.comp2Csv, wl)];
        A0(A0 <= 0) = 1e-4;
        A0 = A0 * diag(1 ./ (sum(A0, 1) + eps));
        gTr = trialMeans(Gm, rw, nTr, L); rTr = trialMeans(Rm, rw, nTr, L);     % reward spectrum of every trial
        gRw = mean(gTr, 2); gDl = gRw - mean(trialMeans(Gm, bl, nTr, L), 2);
        rRw = mean(rTr, 2); rDl = rRw - mean(trialMeans(Rm, bl, nTr, L), 2);
        gAUC0 = winAuc(wl, gRw, bandPeak(wl, sm(wl, gRw), opt.bandT1), opt.aucHalfWidth);
        rAUC0 = winAuc(wl, rRw, bandPeak(wl, sm(wl, rRw), opt.bandT2), opt.aucHalfWidth);
        fprintf('\n==== Mixture %d/%d: %s  (%d trials; measured %s T1 AUC %.4g, %s T2 AUC %.4g)\n', ...
            p, nPair, ttl, nTr, names{1}, gAUC0, names{2}, rAUC0);
        % gains of the two animals at every factor (decoder) or the factor itself (proportional)
        Gain = {ones(numel(wl), nU), ones(numel(wl), nU)};
        extrap = false(nU, 2);
        an = {g.name, r.name}; iaS = [0 0];
        if useDec
            ok = true;
            for s = 1:2
                if s == 1, ia = g.dec; else, ia = r.dec; end
                iaS(s) = ia;
                if ia == 0
                    warn{end+1} = sprintf('%s skipped: no decoder for %s.', ttl, an{s}); %#ok<AGROW>
                    ok = false; break;
                end
                for u = 1:nU
                    [Gain{s}(:, u), extrap(u, s)] = decoderGainAt(DEC{s}, ia, uc(u, s), wl);
                end
            end
            if ~ok, run = run + nU; continue; end
        else
            for u = 1:nU, Gain{1}(:, u) = uc(u, 1); Gain{2}(:, u) = uc(u, 2); end
        end
        % time-averaged ground-truth spectra of the two animals at every factor (figure A), once per animal
        srcM = {Gm, Rm};
        for s = 1:2
            if any([SPEC.dir] == s & strcmpi({SPEC.animal}, an{s})), continue; end
            b0 = mean(srcM{s}, 2); Sx = zeros(numel(wl), numel(fac));
            for f = 1:numel(fac)
                if useDec, Gf = decoderGainAt(DEC{s}, iaS(s), fac(f), wl); else, Gf = fac(f); end
                Sx(:, f) = b0 .* Gf;
            end
            SPEC(end+1) = struct('dir', s, 'animal', an{s}, 'wl', wl, 'S', Sx); %#ok<AGROW>
        end
        clear srcM;
        raw = nan(nU, 2); nmf = nan(nU, 2); rawD = nan(nU, 2); nmfD = nan(nU, 2); iter = nan(nU, 1);
        spc = struct('truthG', {}, 'truthR', {}, 'compG', {}, 'compR', {}, 'r2Mix', {}, 'r2G', {}, 'r2R', {});
        for u = 1:nU
            run = run + 1;
            fg = uc(u, 1); fr = uc(u, 2);
            msg = sprintf('Mixture %d/%d, %s x%g / %s x%g', p, nPair, names{1}, fg, names{2}, fr);
            if run > 1
                el = toc(tRun); msg = sprintf('%s  (about %.1f min left)', msg, el / (run - 1) * (nRun - run + 1) / 60);
            end
            progress(hw, (run - 1) / nRun, msg);
            Y = bsxfun(@times, Gm, Gain{1}(:, u)) + bsxfun(@times, Rm, Gain{2}(:, u));
            [YP, ~, ~] = dataPreprocessing(A0, Y);
            clear Y;
            H0 = max(100 * eps, pinv(A0) * YP);
            [A, nIt] = callNmf(YP, A0, H0, cfg);
            H = max(1e6 * eps, pinv(A' * A) * A' * YP);
            clear YP H0;
            hTr = trialMeans(H, rw, nTr, L);
            hRw = mean(hTr, 2);
            hBl = mean(trialMeans(H, bl, nTr, L), 2);
            % reward spectra (mean over trials) and per-trial R^2 of the unmixing (supplementary spectrum figures)
            tG = bsxfun(@times, gTr, Gain{1}(:, u)); tR = bsxfun(@times, rTr, Gain{2}(:, u));
            cG = A(:, 1) * hTr(1, :); cR = A(:, 2) * hTr(2, :);
            spc(u) = struct('truthG', mean(tG, 2), 'truthR', mean(tR, 2), 'compG', mean(cG, 2), 'compR', mean(cR, 2), ...
                'r2Mix', r2Trials(tG + tR, cG + cR), 'r2G', r2Trials(tG, cG), 'r2R', r2Trials(tR, cR));
            clear tG tR cG cR hTr;
            T  = {gRw .* Gain{1}(:, u), rRw .* Gain{2}(:, u)};
            TD = {gDl .* Gain{1}(:, u), rDl .* Gain{2}(:, u)};
            N  = {A(:, 1) * hRw(1), A(:, 2) * hRw(2)};
            ND = {A(:, 1) * (hRw(1) - hBl(1)), A(:, 2) * (hRw(2) - hBl(2))};
            C = aucPeakCompare(wl, T, N, TD, ND, {'sensor', 'sensor'}, names, opt);
            keys = {'green_T1_AUC', 'red_T2_AUC'}; dkeys = {'green_T1_delta_AUC', 'red_T2_delta_AUC'};
            for k = 1:2
                m = C.m(strcmp({C.m.key}, keys{k}));
                if ~isempty(m), raw(u, k) = m(1).raw; nmf(u, k) = m(1).nmf; end
                m = C.m(strcmp({C.m.key}, dkeys{k}));
                if ~isempty(m), rawD(u, k) = m(1).raw; nmfD(u, k) = m(1).nmf; end
            end
            iter(u) = nIt;
            fprintf('   %s x%-5g %s x%-5g:  %s GT %8.1f unmixed %8.1f (%+6.2f %%)   %s GT %8.1f unmixed %8.1f (%+6.2f %%)  (%d it.)\n', ...
                names{1}, fg, names{2}, fr, names{1}, raw(u, 1), nmf(u, 1), pct(nmf(u, 1), raw(u, 1)), ...
                names{2}, raw(u, 2), nmf(u, 2), pct(nmf(u, 2), raw(u, 2)), nIt);
        end
        pairs(end+1) = struct('green', g.name, 'red', r.name, 'nTrials', nTr, 'gAUC0', gAUC0, 'rAUC0', rAUC0, ...
            'cond', uc, 'raw', raw, 'nmf', nmf, 'rawDelta', rawD, 'nmfDelta', nmfD, 'iter', iter, 'extrap', extrap, ...
            'wl', wl, 'spec', spc); %#ok<AGROW>
        if any(extrap(:))
            warn{end+1} = sprintf('%s: factor outside the measured power range of the decoder (extrapolated).', ttl); %#ok<AGROW>
        end
        clear Gm Rm gTr rTr;
    end
    if isempty(pairs)
        for k = 1:numel(warn), fprintf('   %s\n', warn{k}); end
        error('No mixture could be processed (see the messages above).');
    end
    fprintf('\nAll runs finished (%.1f min).\n', toc(tRun) / 60);

    % ----- results
    R.version = 'NMF_Spike_Recovery 2.0';
    R.windowRule = 'v2: own peak +/- 10 nm (ground truth: ground-truth peak; unmixed: unmixed peak)';
    R.names = names;
    R.colors = cols;
    R.factors = fac;
    R.scaling = P.scaling;
    R.decoderFiles = {'', ''};
    if useDec, R.decoderFiles = {cfg.greenDecoder, cfg.redDecoder}; end
    R.cond = cond;                 % [green factor, red factor, direction, level]
    R.condIdx = cIdx;              % row of cond -> computed condition (pairs(p).cond)
    R.pairs = pairs;
    R.spectra = SPEC;              % time-averaged ground-truth spectra of every animal at every factor
    R.inputs = cfg;
    R.settings = opt;
    R.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    base = fullfile(cfg.outDir, safeFileName(sprintf('%s+%s_SpikeRecovery', names{1}, names{2})));
    progress(hw, 1, 'Saving');
    save([base '.mat'], 'R');
    fprintf('   saved %s.mat\n', base);
    MS = spikeMetrics(R);
    ST = spikeStats(MS, P);
    spikeFigures(R, MS, P, base);
    writeSpikeExcel(R, MS, ST, P, [base '.xlsx']);
    printSummary(R, MS, ST);
    if ~isempty(warn)
        fprintf('\nWarnings (%d):\n', numel(warn)); fprintf('   %s\n', warn{:});
    end
    fprintf('Total time: %.1f min\nOutput folder: %s\n', toc(tTotal) / 60, cfg.outDir);
    if nargout > 0, varargout{1} = R; end
end


function ia = matchAnimal(name, list)
% Index of the animal in list (0 = none): exact (case-insensitive), then after removing
% separators and leading zeros of numbers (FGCa07 = FGCa7 = fgca-7), then a unique name that
% starts with the other and is followed by a non-digit (FGCa7 = FGCa7R25, but FGCa1 ~= FGCa10),
% then a unique name that ends with the other (FGCa7 = GCa7, FjR6 = jR6; never FGCa21 = a1).
    ia = find(strcmpi(list, name), 1);
    if ~isempty(ia), return; end
    nz = @(x) regexprep(regexprep(lower(x), '[^a-z0-9]', ''), '(?<![0-9])0+(?=[0-9])', '');
    a = nz(name); b = cellfun(nz, list, 'UniformOutput', false);
    ia = find(strcmp(b, a), 1);
    if ~isempty(ia), return; end
    hit = false(size(b));
    for k = 1:numel(b)
        [sh, lo] = deal(b{k}, a);
        if numel(sh) > numel(lo), [sh, lo] = deal(lo, sh); end
        if ~isempty(sh) && strncmp(lo, sh, numel(sh)) && ...
                ~(isstrprop(sh(end), 'digit') && isstrprop(lo(numel(sh) + 1), 'digit'))
            hit(k) = true;
        end
    end
    ia = find(hit);
    if numel(ia) == 1, return; end
    hit = false(size(b));
    for k = 1:numel(b)
        [sh, lo] = deal(b{k}, a);
        if numel(sh) > numel(lo), [sh, lo] = deal(lo, sh); end
        if numel(sh) >= 3 && any(isletter(sh)) && any(isstrprop(sh, 'digit')) && numel(lo) > numel(sh) && ...
                strcmp(lo(end - numel(sh) + 1:end), sh) && ...
                ~(isstrprop(sh(1), 'digit') && isstrprop(lo(end - numel(sh)), 'digit'))
            hit(k) = true;
        end
    end
    ia = find(hit);
    if numel(ia) ~= 1, ia = 0; end
end


function [G, extrap] = decoderGainAt(D, ia, k, wl)
% Gain of decoder animal ia at brightness ratio k, on the wavelength axis wl:
% G = (B + k E) / (B + E) (two-component model of NMF_Power_Decoder); k (proportional)
% where the reference spectrum is below minFrac of its maximum; smoothed as in the decoder.
    A = D.animals(ia);
    st = D.settings;
    extrap = k < 1 - 1e-9 || k > max(A.ratio) + 1e-9;
    if abs(k - 1) < 1e-12, G = ones(numel(wl), 1); return; end
    Sref = A.fixed + A.prop;
    ok = Sref >= st.minFrac * max(Sref);
    g = k * ones(size(Sref));
    g(ok) = (A.fixed(ok) + k * A.prop(ok)) ./ Sref(ok);
    g = smoothNm(D.wl, g, st.gainSmoothNm);
    G = interp1(D.wl(:), g(:), wl(:), 'linear', k);
end


function printSummary(R, MS, ST)
    fprintf('\nPer animal (mean over its mixtures):\n');
    for i = 1:numel(MS)
        M = MS(i);
        fprintf('   %s varied  %-10s measured AUC %8.1f  slope %.4f  intercept %+7.1f  recovery %6.1f %%  crosstalk into %s %+6.2f %%\n', ...
            R.names{M.dir}, M.mouse, M.nativeAUC, M.slope, M.intercept, M.recoveryMean, R.names{3 - M.dir}, M.crosstalk);
    end
    fprintf('\nAcross animals (mean [95%% CI], test against the ideal value):\n');
    for i = 1:numel(ST)
        S = ST(i);
        fprintf('   %s varied  %-38s n = %d  %10.4g [%10.4g, %10.4g]  vs %g: %s', R.names{S.dir}, S.measure, S.n, ...
            S.mean, S.ciLo, S.ciHi, S.ref, pStr(S.p));
        if isfinite(S.tostMargin), fprintf('   TOST +/-%g: %s', S.tostMargin, pStr(S.pTOST)); end
        fprintf('\n');
    end
end


%% =====================================================================
%  Inputs and settings
%  =====================================================================
function cfg = askInputs()
    cfg = [];
    d = uigetdir(pwd, 'Green sensor folder (mouse subfolders with G files), e.g. GCaMP');
    if ~ischar(d), return; end
    c.greenPath = d;
    d = uigetdir(fileparts(c.greenPath), 'Red sensor folder (mouse subfolders with R files), e.g. jRGECO or RdLight');
    if ~ischar(d), return; end
    c.redPath = d;
    [f, d] = uigetfile({'*.csv', 'CSV'}, 'Component 1 CSV (green sensor reference spectrum, A0)');
    if ~ischar(f), return; end
    c.comp1Csv = fullfile(d, f);
    [f, d] = uigetfile({'*.csv', 'CSV'}, 'Component 2 CSV (red sensor reference spectrum, A0)', d);
    if ~ischar(f), return; end
    c.comp2Csv = fullfile(d, f);
    d = uigetdir(fileparts(c.greenPath), 'Output folder');
    if ~ischar(d), return; end
    c.outDir = d;
    [~, gn] = fileparts(stripSep(c.greenPath));
    [~, rn] = fileparts(stripSep(c.redPath));
    a = inputdlg({'Green name (component 1)', 'Red name (component 2)'}, 'Names', 1, {gn, rn});
    if isempty(a), return; end
    c.greenName = strtrim(a{1});
    c.redName = strtrim(a{2});
    cfg = c;
end


function c = fillCfg(c, opt)
    def = {'wlRange', opt.wlRange; 'theta1', opt.theta1; 'theta2', opt.theta2; ...
           'alphaA', opt.alphaA; 'alphaX', opt.alphaX; 'greenName', ''; 'redName', ''};
    for k = 1:size(def, 1)
        if ~isfield(c, def{k, 1}) || isempty(c.(def{k, 1})), c.(def{k, 1}) = def{k, 2}; end
    end
    need = {'greenPath', 'redPath', 'comp1Csv', 'comp2Csv', 'outDir'};
    for k = 1:numel(need)
        if ~isfield(c, need{k}) || isempty(c.(need{k})), error('cfg.%s is missing.', need{k}); end
    end
    c.greenPath = stripSep(c.greenPath);
    c.redPath = stripSep(c.redPath);
    if isempty(c.greenName), [~, c.greenName] = fileparts(c.greenPath); end
    if isempty(c.redName), [~, c.redName] = fileparts(c.redPath); end
    if exist(c.greenPath, 'dir') ~= 7, error('Green folder not found: %s', c.greenPath); end
    if exist(c.redPath, 'dir') ~= 7, error('Red folder not found: %s', c.redPath); end
    if exist(c.comp1Csv, 'file') ~= 2, error('Component 1 CSV not found: %s', c.comp1Csv); end
    if exist(c.comp2Csv, 'file') ~= 2, error('Component 2 CSV not found: %s', c.comp2Csv); end
end


function s = defaultSettings()
% Same analysis settings as NMF_Hybrid_Validation_v2.
    s.wlRange         = [450 650];
    s.trialLength     = 1177;
    s.baselineSamples = 100:150;
    s.rewardSamples   = 400:499;
    s.greenTag        = 'G';
    s.redTag          = 'R';
    s.theta1 = 0.01; s.theta2 = 0.01; s.alphaA = 0.001; s.alphaX = 0.001;
    green = [26 141 60] / 255; red = [230 0 18] / 255; gray = [0.50 0.50 0.50];   % #1A8D3C, #E60012
    s.sensorColors = {'GCaMP', green; 'jRGECO', red; 'rdLight', red; 'tdTomato', gray; 'eGFP', gray};
    s.colorComp1 = green;
    s.colorComp2 = red;
    s.bandT1       = [500 540];
    s.bandT2       = [580 620];
    s.aucHalfWidth = 10;
    s.peakSmoothNm = 2;
end


function v = pct(a, b)
    v = 100 * (a - b) / b;
end


function v = ternary(c, a, b)
    if c, v = a; else, v = b; end
end


%% =====================================================================
%  Per-animal recovery metrics (shared with NMF_Spike_Recovery_Compare)
%  =====================================================================
function MS = spikeMetrics(R)
% One entry per (direction, animal whose signal was varied); values are means over
% the mixtures (pairs) that contain the animal.
%   dir 1: green varied, animal = green mouse; dir 2: red varied, animal = red mouse
    fac = R.factors(:);
    nF = numel(fac);
    i1 = find(abs(fac - 1) < 1e-9, 1);
    MS = struct('dir', {}, 'mouse', {}, 'nPairs', {}, 'factors', {}, 'nativeAUC', {}, 'gt', {}, 'unmixed', {}, ...
        'diffPct', {}, 'recovery', {}, 'added', {}, 'partnerChange', {}, 'partnerPct', {}, 'slope', {}, 'intercept', {}, ...
        'r2', {}, 'recoveryMean', {}, 'crosstalk', {}, 'invA', {}, 'invB', {}, 'invR2', {}, 'partnerNames', {});
    for d = 1:2
        o = 3 - d;
        rows = find(R.cond(:, 3) == d);
        [~, ord] = sort(R.cond(rows, 4));
        rows = rows(ord);
        u = R.condIdx(rows);
        if d == 1, who = {R.pairs.green}; else, who = {R.pairs.red}; end
        mice = unique(who);
        for m = 1:numel(mice)
            ps = find(strcmp(who, mice{m}));
            np = numel(ps);
            GT = nan(nF, np); UN = nan(nF, np); PC = nan(nF, np); PP = nan(nF, np); AD = nan(nF, np);
            sl = nan(1, np); ic = nan(1, np); r2 = nan(1, np); xt = nan(1, np);
            ia = nan(1, np); ib = nan(1, np); ir = nan(1, np);
            partners = cell(1, np);
            for q = 1:np
                Q = R.pairs(ps(q));
                t = Q.raw(u, d); n = Q.nmf(u, d); po = Q.nmf(u, o);
                GT(:, q) = t; UN(:, q) = n;
                AD(:, q) = t - t(i1); PC(:, q) = po - po(i1);
                PP(:, q) = 100 * (po - po(i1)) / Q.raw(u(i1), o);   % % of the other sensor's GT AUC
                [cf, rr] = linFit(t, n);
                sl(q) = cf(2); ic(q) = cf(1); r2(q) = rr;
                cx = linFit(t, po);
                xt(q) = 100 * cx(2);
                [ci, rr] = linFit(1 ./ t, 100 * (n - t) ./ t);   % difference (%) = b + a / AUC
                ia(q) = ci(2); ib(q) = ci(1); ir(q) = rr;
                if d == 1, partners{q} = Q.red; else, partners{q} = Q.green; end
            end
            gt = mean(GT, 2); un = mean(UN, 2);
            rec = 100 * (UN - repmat(UN(i1, :), nF, 1)) ./ (GT - repmat(GT(i1, :), nF, 1));
            rec(i1, :) = NaN;
            dp = 100 * (UN - GT) ./ GT;
            recM = mean(rec, 2);
            MS(end+1) = struct('dir', d, 'mouse', mice{m}, 'nPairs', np, 'factors', fac, 'nativeAUC', gt(i1), ...
                'gt', gt, 'unmixed', un, 'diffPct', mean(dp, 2), 'recovery', recM, 'added', mean(AD, 2), ...
                'partnerChange', mean(PC, 2), 'partnerPct', mean(PP, 2), 'slope', mean(sl), 'intercept', mean(ic), 'r2', mean(r2), ...
                'recoveryMean', mean(recM(isfinite(recM))), 'crosstalk', mean(xt), ...
                'invA', mean(ia), 'invB', mean(ib), 'invR2', mean(ir), ...
                'partnerNames', {partners}); %#ok<AGROW>
        end
    end
end


function [cf, r2] = linFit(x, y)
% y = cf(1) + cf(2) * x (least squares)
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    cf = [NaN; NaN]; r2 = NaN;
    if numel(x) < 2, return; end
    X = [ones(numel(x), 1), x];
    cf = X \ y;
    res = y - X * cf;
    r2 = 1 - sum(res .^ 2) / max(sum((y - mean(y)) .^ 2), eps);
end


function C = tintColors(col, n, tmax)
% n shades of col, from col (first) to tmax white (last)
    C = cell(1, n);
    for k = 1:n
        if n > 1, t = tmax * (k - 1) / (n - 1); else, t = 0; end
        C{k} = col * (1 - t) + t;
    end
end


function ST = spikeStats(MS, P)
% Tests across animals (n = animals whose signal was varied), per direction:
%   slope vs 1 (t-test) and equivalence to 1 +/- P.slopeMargin (TOST)
%   intercept vs 0; recovery vs 100 and TOST 100 +/- P.margin
%   a of difference = b + a / AUC vs 0 (trend of the difference with the AUC)
%   crosstalk vs 0
    ST = struct('dir', {}, 'key', {}, 'measure', {}, 'n', {}, 'mean', {}, 'sd', {}, 'ciLo', {}, 'ciHi', {}, ...
        'ref', {}, 't', {}, 'df', {}, 'p', {}, 'tostMargin', {}, 'pTOST', {}, 'equivalent', {}, 'r2', {});
    defs = {'slope', 'Slope (unmixed vs GT AUC)', 1, P.slopeMargin; ...
            'intercept', 'Intercept (counts nm)', 0, NaN; ...
            'recoveryMean', 'Recovery (%)', 100, P.margin; ...
            'invA', 'a of difference (%) = b + a / AUC', 0, NaN; ...
            'crosstalk', 'Crosstalk into the other sensor (%)', 0, NaN};
    for d = 1:2
        M = MS([MS.dir] == d);
        if isempty(M), continue; end
        for k = 1:size(defs, 1)
            v = [M.(defs{k, 1})];
            S = oneSample(v, defs{k, 3});
            pT = NaN; eq = false;
            if isfinite(defs{k, 4})
                T = tost(v(:) - defs{k, 3}, defs{k, 4}, 0.05);
                pT = T.p; eq = T.equivalent;
            end
            r2 = NaN;
            if strcmp(defs{k, 1}, 'invA'), r2 = mean([M.invR2]); end
            ST(end+1) = struct('dir', d, 'key', defs{k, 1}, 'measure', defs{k, 2}, 'n', S.n, 'mean', S.mean, ...
                'sd', S.sd, 'ciLo', S.ciLo, 'ciHi', S.ciHi, 'ref', defs{k, 3}, 't', S.t, 'df', S.df, 'p', S.p, ...
                'tostMargin', defs{k, 4}, 'pTOST', pT, 'equivalent', eq, 'r2', r2); %#ok<AGROW>
        end
    end
end


function S = oneSample(v, ref)
% mean, SD, 95% CI and two-sided one-sample t-test against ref
    v = v(isfinite(v));
    S.n = numel(v); S.mean = mean(v); S.sd = NaN; S.ciLo = NaN; S.ciHi = NaN;
    S.t = NaN; S.df = S.n - 1; S.p = NaN;
    if S.n < 2, return; end
    S.sd = std(v);
    se = S.sd / sqrt(S.n);
    tc = tinvL(0.975, S.df);
    S.ciLo = S.mean - tc * se; S.ciHi = S.mean + tc * se;
    if se > 0
        S.t = (S.mean - ref) / se;
        S.p = 2 * (1 - tcdfL(abs(S.t), S.df));
    end
end


function s = pStr(p)
    if ~isfinite(p), s = 'p = n/a';
    elseif p < 0.001, s = 'p < 0.001';
    else, s = sprintf('p = %.3g', p);
    end
end


%% =====================================================================
%  Figure (2 x 3): one line per animal (darker = larger measured AUC)
%  =====================================================================
function spikeFigure(R, MS, ST, P, base)
    pw = P.panelCm(1); ph = P.panelCm(2); gx = 1.75; gy = 2.9; legW = 2.6;
    x0 = 1.9; y0 = 1.9;
    W = x0 + 3 * pw + 2 * gx + 0.8 + legW;
    Hh = y0 + 2 * ph + gy + 1.2;
    f = newFigure(W, Hh);
    ov = axes(f, 'Units', 'centimeters', 'Position', [0 0 W Hh], 'XLim', [0 W], 'YLim', [0 Hh], 'Visible', 'off');
    hold(ov, 'on');
    tNames = {'T1', 'T2'};
    for d = 1:2
        o = 3 - d;
        M = MS([MS.dir] == d);
        nM = numel(M);
        if nM == 0, continue; end
        [~, o2] = sort([M.nativeAUC], 'descend');   % colour: darker = larger measured AUC
        M = M(o2);
        S = ST([ST.dir] == d);
        g = @(key) S(strcmp({S.key}, key));
        cs = tintColors(R.colors{d}, nM, P.tintMax);
        yb = y0 + (2 - d) * (ph + gy);                % row A on top
        nm = tx(R.names{d}); pn = tx(R.names{o});
        % ---- 1. recovery line
        ax = newAxes(f, [x0, yb, pw, ph], P);
        v = [vertcat(M.gt); vertcat(M.unmixed)] / 1e3;
        [lim, st] = niceAxis([0; v]); lim(1) = 0;
        plot(ax, lim, lim, '--', 'Color', P.refColor, 'LineWidth', 0.6);
        for k = 1:nM
            drawMouse(ax, M(k).gt / 1e3, M(k).unmixed / 1e3, cs{k}, P);
        end
        ylim1 = headroom(lim, st);
        set(ax, 'XLim', lim, 'YLim', ylim1, 'XTick', lim(1):st:lim(2), 'YTick', ylim1(1):st:ylim1(2));
        xlabel(ax, sprintf('%s %s AUC, ground truth (\\times10^{3})', nm, tNames{d}));
        ylabel(ax, sprintf('%s %s AUC, unmixed (\\times10^{3})', nm, tNames{d}));
        panelTitle(ax, 'Recovery line', P);
        a = g('slope'); b = g('intercept');
        panelNote(ax, {sprintf('slope %.3f [%.3f, %.3f]', a.mean, a.ciLo, a.ciHi), ...
            sprintf('  vs 1: %s;  TOST \\pm%g: %s', pStr(a.p), a.tostMargin, pStr(a.pTOST)), ...
            sprintf('intercept %.0f [%.0f, %.0f]', b.mean, b.ciLo, b.ciHi), ...
            sprintf('  vs 0: %s', pStr(b.p))}, P);
        % ---- 2. difference (%) vs ground-truth AUC
        ax = newAxes(f, [x0 + (pw + gx), yb, pw, ph], P);
        v = vertcat(M.diffPct);
        [yl, st] = niceAxis([v; -P.margin; P.margin]);
        yl = headroom(yl, st);
        plot(ax, lim, [0 0], '--', 'Color', P.refColor, 'LineWidth', 0.6);
        plot(ax, lim, -P.margin * [1 1], ':', 'Color', P.refColor, 'LineWidth', 0.6);
        plot(ax, lim, P.margin * [1 1], ':', 'Color', P.refColor, 'LineWidth', 0.6);
        for k = 1:nM
            drawMouse(ax, M(k).gt / 1e3, M(k).diffPct, cs{k}, P);
        end
        set(ax, 'XLim', lim, 'XTick', lim(1):niceStep(lim(2)):lim(2), 'YLim', yl, 'YTick', yl(1):st:yl(2));
        xlabel(ax, sprintf('%s %s AUC, ground truth (\\times10^{3})', nm, tNames{d}));
        ylabel(ax, 'difference (%)');
        panelTitle(ax, sprintf('%s difference', nm), P);
        a = g('invA');
        panelNote(ax, {'difference = b + a / AUC', sprintf('a = %.3g [%.3g, %.3g]', a.mean, a.ciLo, a.ciHi), ...
            sprintf('vs 0: %s;  R^2 %.2f', pStr(a.p), a.r2)}, P, 'right');
        % ---- 3. change of the other sensor
        ax = newAxes(f, [x0 + 2 * (pw + gx), yb, pw, ph], P);
        ya = vertcat(M.partnerPct);
        [yl, st] = niceAxis([ya; -P.margin; P.margin]);
        yl = headroom(yl, st);
        plot(ax, lim, [0 0], '--', 'Color', P.refColor, 'LineWidth', 0.6);
        plot(ax, lim, -P.margin * [1 1], ':', 'Color', P.refColor, 'LineWidth', 0.6);
        plot(ax, lim, P.margin * [1 1], ':', 'Color', P.refColor, 'LineWidth', 0.6);
        for k = 1:nM
            drawMouse(ax, M(k).gt / 1e3, M(k).partnerPct, cs{k}, P);
        end
        set(ax, 'XLim', lim, 'XTick', lim(1):niceStep(lim(2)):lim(2), 'YLim', yl, 'YTick', yl(1):st:yl(2));
        xlabel(ax, sprintf('%s %s AUC, ground truth (\\times10^{3})', nm, tNames{d}));
        ylabel(ax, sprintf('%s change (%% of its AUC)', pn));
        panelTitle(ax, sprintf('%s (other sensor)', pn), P);
        a = g('crosstalk');
        panelNote(ax, {sprintf('crosstalk %.2f [%.2f, %.2f] %%', a.mean, a.ciLo, a.ciHi), ...
            sprintf('  vs 0: %s', pStr(a.p))}, P);
        % ---- row label and animals
        text(ov, 0.3, yb + ph + 0.65, sprintf(['%s.  %s varied (%s at its measured level); one line per %s mouse ' ...
            '(legend: measured AUC x10^3; darker = larger)'], char('A' + d - 1), R.names{d}, R.names{o}, R.names{d}), 'FontName', P.fontName, 'FontSize', P.fontSize + 1, ...
            'FontWeight', 'bold', 'Interpreter', 'none', 'VerticalAlignment', 'bottom');
        lx = x0 + 3 * pw + 2 * gx + 0.8;
        for k = 1:nM
            yy = yb + ph - 0.2 - (k - 1) * 0.42;
            plot(ov, lx, yy, 'o', 'MarkerSize', P.dotSize, 'MarkerFaceColor', cs{k}, 'MarkerEdgeColor', 'none');
            text(ov, lx + 0.3, yy, sprintf('%s (%.1f)', M(k).mouse, M(k).nativeAUC / 1e3), 'FontName', P.fontName, 'FontSize', P.fontSize - 1, ...
                'Interpreter', 'none', 'VerticalAlignment', 'middle');
        end
    end
    text(ov, 0.3, 0.25, sprintf(['Values: mean [95%% CI] over animals (n = animals whose signal was varied); ' ...
        'one-sample t-tests against the ideal value; TOST: two one-sided tests (equivalence).']), ...
        'FontName', P.fontName, 'FontSize', P.fontSize - 1, 'Color', [0.35 0.35 0.35], 'Interpreter', 'none');
    saveFigure(f, base);
end


function yl = headroom(yl, st)
% extra space above the data for the statistics text
    yl(2) = yl(2) + ceil(0.4 * (yl(2) - yl(1)) / st) * st;
end


function drawMouse(ax, x, y, c, P)
    plot(ax, x, y, '-', 'Color', c, 'LineWidth', P.lineWidth);
    plot(ax, x, y, 'o', 'MarkerSize', P.dotSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'none', 'LineStyle', 'none');
end


function panelTitle(ax, s, P)
    text(ax, 0, 1.04, s, 'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontName', P.fontName, ...
        'FontSize', P.fontSize, 'FontWeight', 'bold', 'Interpreter', 'none');
end


function panelNote(ax, s, P, side)
    if nargin < 4, side = 'left'; end
    y = 0.97; va = 'top';
    if strcmp(side, 'bottomright'), side = 'right'; y = 0.03; va = 'bottom'; end
    if strcmp(side, 'right'), x = 0.97; else, x = 0.04; end
    text(ax, x, y, s, 'Units', 'normalized', 'VerticalAlignment', va, 'HorizontalAlignment', side, ...
        'FontName', P.fontName, 'FontSize', P.fontSize - 1, 'Color', [0.25 0.25 0.25]);
end


function s = facLabel(v)
    if abs(v - round(v)) < 1e-9, s = sprintf('%d', round(v)); else, s = sprintf('%g', v); end
end


function s = sdv(v)
    v = v(isfinite(v));
    if numel(v) > 1, s = std(v); else, s = 0; end
end


%% =====================================================================
%  Excel
%  =====================================================================
function writeSpikeExcel(R, MS, ST, P, xlsFile)
    if exist(xlsFile, 'file') == 2
        try, delete(xlsFile); catch, end %#ok<NOCOM>
    end
    scTxt = 'whole recording of the varied sensor x factor (spectrum and time course, shape unchanged)';
    if isfield(R, 'scaling') && strcmpi(R.scaling, 'decoder')
        scTxt = sprintf(['every wavelength of the varied recording x the power-decoder gain of that animal ' ...
            '(NMF_Power_Decoder): %s; %s'], R.decoderFiles{1}, R.decoderFiles{2});
    end
    try
        C = { ...
            'NMF-RI spike-and-recovery (in silico mixtures of measured single-sensor data)', ''; ...
            'Version', R.version; 'Date', R.date; ...
            'Green', sprintf('%s  %s', R.names{1}, R.inputs.greenPath); ...
            'Red', sprintf('%s  %s', R.names{2}, R.inputs.redPath); ...
            'Factors', sprintf('brightness ratios %s of the varied sensor (other sensor at its measured level)', num2str(R.factors)); ...
            'A / B', sprintf('A: %s varied (animal = %s mouse); B: %s varied (animal = %s mouse)', ...
                R.names{1}, R.names{1}, R.names{2}, R.names{2}); ...
            'Scaling', scTxt; ...
            'Ground truth', 'measured single-sensor data minus the file offset, scaled likewise'; ...
            'NMF-RI', 'as NMF_RI_Photometry (A0 from the CSVs, dataPreprocessing, NMF_RI, H by least squares)'; ...
            'AUC windows', R.windowRule; ...
            'Animal level', 'mean over the mixtures (partners) of the animal whose signal was varied'; ...
            'Slope, intercept', 'least squares unmixed AUC = intercept + slope x ground-truth AUC (absolute, counts nm) over the factors; slope 1 = proportional'; ...
            'Recovery (%)', '(unmixed(f) - unmixed(x1)) / (GT(f) - GT(x1)) x 100; mean over f ~= 1'; ...
            'Difference (%)', '(unmixed - GT) / GT x 100'; ...
            'Crosstalk (%)', 'slope of the other sensor''s unmixed AUC vs the varied sensor''s GT AUC x 100 (other sensor GT constant)'; ...
            'DiffRatio', ['every mixture and factor: x = input ratio = total AUC (all wavelengths) of the ground-truth ' ...
                'reward spectrum of the unchanged sensor (1x) / that of the varied sensor (nx), as the factors (total AUC); ' ...
                'y = difference (%) (T1 / T2 window AUC) of the varied sensor (least squares c + b x) and of the unchanged ' ...
                'sensor (c + b / x). b (%) = difference added at x = 1 (equal input); c (%) = difference without the partner. ' ...
                sprintf('95%% CI: bootstrap (%d resamples) over animals', P.nBoot)]; ...
            'WindowScaling', ['time-averaged ground-truth spectrum of every animal at every factor: total AUC (all ' ...
                'wavelengths), T1 and T2 AUC (+/- 10 nm around the peak of the x1 spectrum in the T1 / T2 band, the same ' ...
                'window at every factor), each as a ratio to x1 of the same animal. Slope of the T1 / T2 ratio vs the total ' ...
                'ratio through (1, 1) per animal (1 = the window changes as the total); mean [95% CI] and t-test vs 1']; ...
            'Stats', sprintf(['across animals: one-sample t-tests (slope vs 1, intercept vs 0, recovery vs 100, a vs 0, ' ...
                'crosstalk vs 0), 95%% CI; TOST: slope 1 +/- %g, recovery 100 +/- %g %%'], P.slopeMargin, P.margin)};
        writecell(clean(C), xlsFile, 'Sheet', 'README');

        C = {'Green mouse', 'Red mouse', 'Trials', sprintf('Measured %s T1 AUC', R.names{1}), ...
             sprintf('Measured %s T2 AUC', R.names{2}), 'Green factor', 'Red factor', 'NMF-RI iterations', ...
             sprintf('%s T1 AUC GT', R.names{1}), sprintf('%s T1 AUC unmixed', R.names{1}), sprintf('%s difference (%%)', R.names{1}), ...
             sprintf('%s T2 AUC GT', R.names{2}), sprintf('%s T2 AUC unmixed', R.names{2}), sprintf('%s difference (%%)', R.names{2}), ...
             sprintf('%s T1 delta AUC GT', R.names{1}), sprintf('%s T1 delta AUC unmixed', R.names{1}), ...
             sprintf('%s T2 delta AUC GT', R.names{2}), sprintf('%s T2 delta AUC unmixed', R.names{2})};
        for p = 1:numel(R.pairs)
            Q = R.pairs(p);
            for u = 1:size(Q.cond, 1)
                C(end+1, :) = {Q.green, Q.red, Q.nTrials, Q.gAUC0, Q.rAUC0, Q.cond(u, 1), Q.cond(u, 2), Q.iter(u), ...
                    Q.raw(u, 1), Q.nmf(u, 1), 100 * (Q.nmf(u, 1) - Q.raw(u, 1)) / Q.raw(u, 1), ...
                    Q.raw(u, 2), Q.nmf(u, 2), 100 * (Q.nmf(u, 2) - Q.raw(u, 2)) / Q.raw(u, 2), ...
                    Q.rawDelta(u, 1), Q.nmfDelta(u, 1), Q.rawDelta(u, 2), Q.nmfDelta(u, 2)}; %#ok<AGROW>
            end
        end
        writecell(clean(C), xlsFile, 'Sheet', 'Pairs');

        fac = R.factors(:)';
        hdr = {'Varied sensor', 'Animal', 'n mixtures', 'Measured AUC (x1)', 'Slope', 'Intercept (counts nm)', ...
               'R^2', 'Mean recovery (%)', 'Crosstalk into other sensor (%)'};
        for k = 1:numel(fac), hdr{end+1} = sprintf('Recovery x%g (%%)', fac(k)); end %#ok<AGROW>
        for k = 1:numel(fac), hdr{end+1} = sprintf('Difference x%g (%%)', fac(k)); end %#ok<AGROW>
        for k = 1:numel(fac), hdr{end+1} = sprintf('GT AUC x%g', fac(k)); end %#ok<AGROW>
        for k = 1:numel(fac), hdr{end+1} = sprintf('Unmixed AUC x%g', fac(k)); end %#ok<AGROW>
        for k = 1:numel(fac), hdr{end+1} = sprintf('Other sensor change x%g', fac(k)); end %#ok<AGROW>
        for k = 1:numel(fac), hdr{end+1} = sprintf('Other sensor change x%g (%% of its AUC)', fac(k)); end %#ok<AGROW>
        C = hdr;
        for i = 1:numel(MS)
            M = MS(i);
            C(end+1, :) = [{R.names{M.dir}, M.mouse, M.nPairs, M.nativeAUC, M.slope, M.intercept, M.r2, ...
                M.recoveryMean, M.crosstalk}, num2cell(M.recovery(:)'), num2cell(M.diffPct(:)'), ...
                num2cell(M.gt(:)'), num2cell(M.unmixed(:)'), num2cell(M.partnerChange(:)'), num2cell(M.partnerPct(:)')]; %#ok<AGROW>
        end
        writecell(clean(C), xlsFile, 'Sheet', 'Mice');

        C = {'Varied sensor', 'Factor', 'n animals', 'Recovery mean (%)', 'Recovery SD (%)', ...
             'Difference mean (%)', 'Difference SD (%)', 'Other sensor change mean (counts nm)', 'Other sensor change SD'};
        for d = 1:2
            M = MS([MS.dir] == d);
            if isempty(M), continue; end
            Rc = [M.recovery]; Dp = [M.diffPct]; Pc = [M.partnerChange];
            for k = 1:numel(fac)
                C(end+1, :) = {R.names{d}, fac(k), numel(M), mean(Rc(k, :)), sdv(Rc(k, :)), mean(Dp(k, :)), ...
                    sdv(Dp(k, :)), mean(Pc(k, :)), sdv(Pc(k, :))}; %#ok<AGROW>
            end
        end
        writecell(clean(C), xlsFile, 'Sheet', 'Levels');

        C = {'Varied sensor', 'Measure', 'n animals', 'Mean', 'SD', '95% CI lower', '95% CI upper', 'Ideal value', ...
             't (vs ideal)', 'df', 'p (vs ideal)', 'TOST margin (+/-)', 'p TOST', 'Equivalent', 'Mean R^2'};
        for i = 1:numel(ST)
            S = ST(i);
            eqs = ''; if isfinite(S.tostMargin), eqs = ternary(S.equivalent, 'equivalent', 'not equivalent'); end
            C(end+1, :) = {R.names{S.dir}, S.measure, S.n, S.mean, S.sd, S.ciLo, S.ciHi, S.ref, S.t, S.df, S.p, ...
                S.tostMargin, S.pTOST, eqs, S.r2}; %#ok<AGROW>
        end
        writecell(clean(C), xlsFile, 'Sheet', 'Stats');
        writeIncreaseSheets(R, MS, xlsFile);
        writeRatioSheets(R, P, xlsFile);
        writeWindowSheets(R, xlsFile);
        fprintf('   saved %s\n', xlsFile);
    catch ME
        fprintf('   Excel file not written: %s\n', ME.message);
    end
end


%% =====================================================================
%  Figures: A spectra; B T1 / T2 AUC (ground truth vs unmixed) and C difference, for the
%  varied sensor and for the other sensor
%  (one set per sensor; factor colours from light (x1) to dark (largest))
%  =====================================================================
function spikeFigures(R, MS, P, base)
    % sensor colours from the current settings (also when redrawing an older .mat);
    % darkest shade = this colour (green #1A8D3C, red #E60012)
    ds = defaultSettings();
    R.colors = {sensorColor(R.names{1}, ds, ds.colorComp1), sensorColor(R.names{2}, ds, ds.colorComp2)};
    for d = 1:2
        nm = safeFileName(R.names{d});
        if isfield(R, 'spectra') && ~isempty(R.spectra)
            spectraFigure(R, d, P, sprintf('%s_spectra_%s', base, nm));
        end
        aucFigure(R, d, P, sprintf('%s_AUC_%s_varied', base, nm));
        ratioFigure(R, d, P, sprintf('%s_DiffRatio_%s_varied', base, nm));
        windowFigure(R, d, P, sprintf('%s_WindowScaling_%s', base, nm));
    end
    supplementFigures(R, P, base);
end


function C = factorColors(col, n)
% n shades of col: lightest first (x1), darkest last
    C = cell(1, n);
    for k = 1:n
        t = 0.75 * (n - k) / max(n - 1, 1);
        C{k} = col * (1 - t) + t;
    end
end


function spectraFigure(R, d, P, fileBase)
% A. time-averaged ground-truth spectrum of every animal at every factor, overlaid
    Sp = R.spectra([R.spectra.dir] == d);
    if isempty(Sp), return; end
    fac = R.factors(:)'; nF = numel(fac);
    cs = factorColors(R.colors{d}, nF);
    nA = numel(Sp); nc = min(4, nA); nr = ceil(nA / nc);
    pw = 4.0; ph = 3.4; gx = 1.7; gy = 1.6; x0 = 1.8; y0 = 1.3;
    W = x0 + nc * pw + (nc - 1) * gx + 2.2; H = y0 + nr * ph + (nr - 1) * gy + 1.2;
    f = newFigure(W, H);
    for a = 1:nA
        r = ceil(a / nc); c = a - (r - 1) * nc;
        ax = newAxes(f, [x0 + (c - 1) * (pw + gx), y0 + (nr - r) * (ph + gy), pw, ph], P);
        wl = Sp(a).wl; S = Sp(a).S;
        for k = 1:nF
            plot(ax, wl, S(:, k), '-', 'Color', cs{k}, 'LineWidth', P.lineWidth);
        end
        [yl, st] = niceAxis([0; S(:)]); yl(1) = 0;
        set(ax, 'XLim', [450 650], 'XTick', 450:50:650, 'YLim', yl, 'YTick', yl(1):st:yl(2));
        text(ax, 0, 1.04, Sp(a).animal, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
            'FontName', P.fontName, 'FontSize', P.fontSize, 'FontWeight', 'bold', 'Interpreter', 'none');
        if c == 1, ylabel(ax, 'Intensity (counts)'); end
        if r == nr, xlabel(ax, 'Wavelength (nm)'); end
    end
    ov = axes(f, 'Units', 'centimeters', 'Position', [0 0 W H], 'XLim', [0 W], 'YLim', [0 H], 'Visible', 'off');
    hold(ov, 'on');
    lx = x0 + nc * pw + (nc - 1) * gx + 0.5;
    for k = 1:nF
        yy = y0 + nr * ph + (nr - 1) * gy - 0.2 - (k - 1) * 0.42;
        plot(ov, lx + [-0.3 0.3], [yy yy], '-', 'Color', cs{k}, 'LineWidth', 1.2);
        text(ov, lx + 0.45, yy, sprintf('%g%s', fac(k), char(215)), 'FontName', P.fontName, 'FontSize', P.fontSize - 1);
    end
    text(ov, 0.3, H - 0.3, sprintf('%s: ground-truth spectra (time average) at total AUC %s', R.names{d}, ...
        strjoin(arrayfun(@(v) sprintf('%g%s', v, char(215)), fac, 'UniformOutput', false), ', ')), ...
        'FontName', P.fontName, 'FontSize', P.fontSize - 1, 'VerticalAlignment', 'top', 'Interpreter', 'none');
    saveFigure(f, fileBase);
end


function Q = aucData(R, d, s)
% Animal-level AUC of sensor s (1: green T1, 2: red T2) when sensor d is varied, at every
% factor: ground truth, unmixed and difference (%); animal = mouse of sensor s, mean over the
% mixtures (partners) that contain it.
    fac = R.factors(:)'; nF = numel(fac);
    rows = find(R.cond(:, 3) == d);
    [~, ord] = sort(R.cond(rows, 4)); rows = rows(ord);
    u = R.condIdx(rows);
    if s == 1, who = {R.pairs.green}; else, who = {R.pairs.red}; end
    mice = unique(who);
    Q.fac = fac; Q.animals = mice; Q.d = d; Q.s = s;
    Q.gt = nan(numel(mice), nF); Q.un = Q.gt; Q.diff = Q.gt;
    for a = 1:numel(mice)
        ps = find(strcmp(who, mice{a}));
        G = nan(nF, numel(ps)); U = G;
        for q = 1:numel(ps)
            G(:, q) = R.pairs(ps(q)).raw(u, s); U(:, q) = R.pairs(ps(q)).nmf(u, s);
        end
        Q.gt(a, :) = mean(G, 2)'; Q.un(a, :) = mean(U, 2)';
        Q.diff(a, :) = mean(100 * (U - G) ./ G, 2)';
    end
    Q.ccc = lincc(Q.gt(:), Q.un(:));
    % per factor: Lin's CCC and the line unmixed = a + b * ground truth (panel B)
    Q.cccF = struct('k', {}, 'ccc', {}, 'lo', {}, 'hi', {}, 'n', {}, 'a', {}, 'b', {}, 'r2', {}, 'xMin', {}, 'xMax', {});
    for j = 1:nF
        xx = Q.gt(:, j); yy = Q.un(:, j);
        ok = isfinite(xx) & isfinite(yy);
        C = lincc(xx(ok), yy(ok));
        [cf, r2] = linFit(xx(ok), yy(ok));
        Q.cccF(end+1) = struct('k', fac(j), 'ccc', C.ccc, 'lo', C.lo, 'hi', C.hi, 'n', sum(ok), 'a', cf(1), ...
            'b', cf(2), 'r2', r2, 'xMin', min(xx(ok)), 'xMax', max(xx(ok)));
    end
    Q.fit = struct('k', {}, 'a', {}, 'b', {}, 'r2', {}, 'xMin', {}, 'xMax', {}, 'n', {});
    for j = 1:nF
        xx = Q.gt(:, j); yy = Q.diff(:, j);
        ok = isfinite(xx) & isfinite(yy);
        [cf, r2] = linFit(xx(ok), yy(ok));
        Q.fit(end+1) = struct('k', fac(j), 'a', cf(1), 'b', cf(2), 'r2', r2, 'xMin', min(xx(ok)), ...
            'xMax', max(xx(ok)), 'n', sum(ok));
    end
end


function aucFigure(R, d, P, fileBase)
% One file per varied sensor d: row 1 = the varied sensor, row 2 = the other sensor
%   B. unmixed vs ground-truth T1 / T2 AUC: line per factor through the animal means (colour =
%      factor of the varied sensor, light x1 ... dark), 1:1 dashed line, Lin's CCC [95% CI] per factor
%   C. difference (%) vs ground-truth T1 / T2 AUC: regression line per factor with R^2
    AX = axisGroups(R, P, false);
    tName = {'T1', 'T2'};
    pw = 4.6; ph = 4.6; gx = 2.2; gy = 2.6; x0 = 1.8; y0 = 1.4;
    W = x0 + 2 * pw + gx + 3.4; H = y0 + 2 * ph + gy + 2.8;
    f = newFigure(W, H);
    ov = axes(f, 'Units', 'centimeters', 'Position', [0 0 W H], 'XLim', [0 W], 'YLim', [0 H], 'Visible', 'off');
    hold(ov, 'on');
    order = [d, 3 - d];
    for row = 1:2
        s = order(row);
        Q = aucData(R, d, s);
        if isempty(Q.animals), continue; end
        nF = numel(Q.fac);
        cs = factorColors(R.colors{s}, nF);
        nm = tx(R.names{s});
        yb = y0 + (2 - row) * (ph + gy);
        % ---- B
        ax = newAxes(f, [x0, yb, pw, ph], P);
        lim = AX(d, s).lim; st = AX(d, s).st;          % same axes as the unmixing\ figures of this varied sensor
        plot(ax, lim, lim, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.6);
        labB = {};
        for j = 1:numel(Q.cccF)
            F = Q.cccF(j);
            if ~isfinite(F.b), continue; end
            xr = [F.xMin F.xMax];
            plot(ax, xr / 1e3, (F.a + F.b * xr) / 1e3, '-', 'Color', cs{j}, 'LineWidth', 1.5);
            labB{end+1} = sprintf('\\color[rgb]{%.2f %.2f %.2f}%g%s: %.4f [%.4f, %.4f]', cs{j}, F.k, char(215), ...
                F.ccc, F.lo, F.hi); %#ok<AGROW>
        end
        set(ax, 'XLim', lim, 'YLim', lim, 'XTick', lim(1):st:lim(2), 'YTick', lim(1):st:lim(2));
        xlabel(ax, sprintf('Ground-truth %s AUC (10^{3} counts \\cdot nm)', tName{s}));
        ylabel(ax, sprintf('Unmixed %s AUC (10^{3} counts \\cdot nm)', tName{s}));
        text(ax, 0, 1.04, sprintf('%s %s AUC', nm, tName{s}), 'Units', 'normalized', ...
            'VerticalAlignment', 'bottom', 'FontName', P.fontName, 'FontSize', P.fontSize, 'FontWeight', 'bold');
        text(ax, 0.98, 0.03, [{sprintf('CCC [95%% CI], n = %d mice', numel(Q.animals))}, labB], 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontName', P.fontName, 'FontSize', P.fontSize - 2);
        % ---- C
        ax = newAxes(f, [x0 + pw + gx, yb, pw, ph], P);
        yv = []; xv = 0; lab = {};
        for j = 1:numel(Q.fit)
            F = Q.fit(j);
            if ~isfinite(F.b), continue; end
            xr = [F.xMin F.xMax];
            if diff(xr) <= 0, xr = xr + [-0.01 0.01] * max(xr(1), 1); end
            plot(ax, xr / 1e3, F.a + F.b * xr, '-', 'Color', cs{j}, 'LineWidth', 1.5);
            yv = [yv, F.a + F.b * xr]; xv = [xv, xr / 1e3]; %#ok<AGROW>
            lab{end+1} = sprintf('\\color[rgb]{%.2f %.2f %.2f}%g%s: R^2 = %.2f', cs{j}, F.k, char(215), F.r2); %#ok<AGROW>
        end
        xl = AX(d, s).lim; sx = AX(d, s).st;
        yl = AX(d, s).baY; sy = AX(d, s).baSt;         % same as the Bland-Altman panels of unmixing\ (this series)
        if yl(1) < 0 && yl(2) > 0
            plot(ax, xl, [0 0], '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.6);
        end
        set(ax, 'XLim', xl, 'XTick', xl(1):sx:xl(2), 'YLim', yl, 'YTick', yl(1):sy:yl(2));
        xlabel(ax, sprintf('Ground-truth %s AUC (10^{3} counts \\cdot nm)', tName{s}));
        ylabel(ax, 'Difference (%)');
        text(ax, 0, 1.04, sprintf('%s %s difference', nm, tName{s}), 'Units', 'normalized', ...
            'VerticalAlignment', 'bottom', 'FontName', P.fontName, 'FontSize', P.fontSize, 'FontWeight', 'bold');
        text(ax, 0.97, 0.97, lab, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'FontName', P.fontName, 'FontSize', P.fontSize - 1);
        % row label and factor legend
        if row == 1, rl = 'scaled'; else, rl = 'unscaled'; end
        text(ov, 0.3, yb + ph + 0.75, sprintf('%s %s (%s)', nm, tName{s}, rl), 'FontName', P.fontName, ...
            'FontSize', P.fontSize, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
        lx = x0 + 2 * pw + gx + 0.6;
        for j = 1:nF
            yy = yb + ph - 0.2 - (j - 1) * 0.42;
            plot(ov, lx + [-0.25 0.15], [yy yy], '-', 'Color', cs{j}, 'LineWidth', 1.5);
            text(ov, lx + 0.3, yy, sprintf('%s %g%s', tx(R.names{d}), Q.fac(j), char(215)), 'FontName', P.fontName, ...
                'FontSize', P.fontSize - 1);
        end
    end
    text(ov, 0.3, H - 0.3, sprintf('%s %g%s - %g%s, %s unchanged', tx(R.names{d}), min(R.factors), char(215), ...
        max(R.factors), char(215), tx(R.names{3 - d})), 'FontName', P.fontName, 'FontSize', P.fontSize + 1, ...
        'FontWeight', 'bold', 'VerticalAlignment', 'top');
    text(ov, 0.3, H - 0.85, sprintf(['Total AUC of %s %s; animal means (mean over the mixtures of each animal); ' ...
        'colour = factor of %s'], R.names{d}, strjoin(arrayfun(@(v) sprintf('%g%s', v, char(215)), R.factors, 'UniformOutput', false), ', '), ...
        R.names{d}), 'FontName', P.fontName, 'FontSize', P.fontSize - 1, 'VerticalAlignment', 'top', 'Interpreter', 'none');
    saveFigure(f, fileBase);
end


function writeIncreaseSheets(R, MS, xlsFile) %#ok<INUSL>
% Excel sheets for redrawing B and C: AUC (every animal and factor) and AUCStats
    C = {'Sensor varied', 'Sensor shown', 'Animal', 'Factor (total AUC of the varied sensor)', ...
         'Ground-truth T1/T2 AUC', 'Unmixed T1/T2 AUC', 'Difference (%)'};
    S = {'Sensor varied', 'Sensor shown', 'Statistic', 'Factor', 'Value 1', 'Value 2', 'Value 3', 'n'};
    for d = 1:2
        for s = [d, 3 - d]
            Q = aucData(R, d, s);
            for a = 1:numel(Q.animals)
                for j = 1:numel(Q.fac)
                    C(end+1, :) = {R.names{d}, R.names{s}, Q.animals{a}, Q.fac(j), Q.gt(a, j), Q.un(a, j), Q.diff(a, j)}; %#ok<AGROW>
                end
            end
            S(end+1, :) = {R.names{d}, R.names{s}, 'CCC, 95% CI lower, upper', 'all', Q.ccc.ccc, Q.ccc.lo, Q.ccc.hi, Q.ccc.n}; %#ok<AGROW>
            for j = 1:numel(Q.cccF)
                F = Q.cccF(j);
                S(end+1, :) = {R.names{d}, R.names{s}, 'CCC, 95% CI lower, upper', F.k, F.ccc, F.lo, F.hi, F.n}; %#ok<AGROW>
                S(end+1, :) = {R.names{d}, R.names{s}, 'Unmixed vs ground truth: intercept, slope, R^2', F.k, F.a, F.b, ...
                    F.r2, F.n}; %#ok<AGROW>
            end
            for j = 1:numel(Q.fit)
                F = Q.fit(j);
                S(end+1, :) = {R.names{d}, R.names{s}, 'Difference vs AUC: intercept, slope (per counts nm), R^2', ...
                    F.k, F.a, F.b, F.r2, F.n}; %#ok<AGROW>
            end
        end
    end
    writecell(clean(C), xlsFile, 'Sheet', 'AUC');
    writecell(clean(S), xlsFile, 'Sheet', 'AUCStats');
end


%% =====================================================================
%  Difference vs input ratio (one figure per experiment = per varied sensor d, two panels)
%    x = input ratio: total AUC (all wavelengths) of the ground-truth reward spectrum of the unchanged
%        sensor (1x) / that of the varied sensor (nx); total AUC as the factors
%    y = difference (%) (T1 / T2 window AUC) of the varied sensor (left) and of the unchanged sensor (right);
%        every mixture and factor
%    fit: varied sensor   difference = c + b x     (partner / itself = x)
%         unchanged sensor difference = c + b / x   (partner / itself = 1 / x)
%       b (%) = difference added at x = 1 (equal input of the two sensors)
%       c (%) = difference without the partner (varied: x -> 0; unchanged: x -> infinity)
%    95% CI of b and c: bootstrap over animals (green and red mice resampled independently)
%  =====================================================================
function Z = ratioData(R, d, s, P)
% d = varied sensor; s = sensor whose difference is shown (s = d: the varied sensor itself;
% s = 3 - d: the unchanged sensor). x = input ratio, the same in both cases.
    o = 3 - d;                                     % unchanged sensor
    fac = R.factors(:)'; nF = numel(fac);
    rows = find(R.cond(:, 3) == d);
    [~, ord] = sort(R.cond(rows, 4)); rows = rows(ord);
    u = R.condIdx(rows);
    nP = numel(R.pairs); nPt = nP * nF;
    Z.d = d; Z.s = s; Z.o = o; Z.fac = fac;
    if s == d, Z.model = 'lin'; else, Z.model = 'recip'; end
    Z.green = cell(nPt, 1); Z.red = cell(nPt, 1);
    Z.pair = zeros(nPt, 1); Z.k = zeros(nPt, 1); Z.lev = zeros(nPt, 1);
    Z.gtSh = nan(nPt, 1); Z.unSh = nan(nPt, 1); Z.gtU = nan(nPt, 1); Z.gtV = nan(nPt, 1);
    Z.totU = nan(nPt, 1); Z.totV = nan(nPt, 1);
    i = 0;
    for p = 1:nP
        Q = R.pairs(p);
        for j = 1:nF
            i = i + 1;
            Z.green{i} = Q.green; Z.red{i} = Q.red; Z.pair(i) = p; Z.k(i) = fac(j); Z.lev(i) = j;
            Z.gtSh(i) = Q.raw(u(j), s); Z.unSh(i) = Q.nmf(u(j), s);
            Z.gtU(i) = Q.raw(u(j), o); Z.gtV(i) = Q.raw(u(j), d);
            if isfield(Q, 'spec') && numel(Q.spec) >= u(j) && ~isempty(Q.spec(u(j)).truthG)
                tot = [trapz(Q.wl(:), Q.spec(u(j)).truthG(:)), trapz(Q.wl(:), Q.spec(u(j)).truthR(:))];
                Z.totU(i) = tot(o); Z.totV(i) = tot(d);
            end
        end
    end
    % x = input ratio: total AUC (all wavelengths) of the ground-truth reward spectra, as the factors
    % (total AUC); older results without spectra: T1 / T2 window AUC
    Z.xKind = 'total';
    if all(isfinite(Z.totU) & isfinite(Z.totV))
        Z.x = Z.totU ./ Z.totV;
    else
        Z.x = Z.gtU ./ Z.gtV; Z.xKind = 'window';
        fprintf('   [note] no ground-truth spectra in this result: ratio from the T1 / T2 window AUC\n');
    end
    Z.y = 100 * (Z.unSh - Z.gtSh) ./ Z.gtSh;
    ok = isfinite(Z.x) & isfinite(Z.y) & Z.x > 0;
    Z.ok = ok; Z.n = sum(ok); Z.nPairs = numel(unique(Z.pair(ok)));
    [Z.c, Z.b, Z.r2] = ratioFit(Z.x(ok), Z.y(ok), Z.model);
    Z.fitted = Z.c + Z.b * ratioG(Z.x, Z.model);
    Z.cCI = [NaN NaN]; Z.bCI = [NaN NaN]; Z.nBoot = 0; Z.bootUnit = '';
    if ~isfield(P, 'nBoot') || P.nBoot < 1 || Z.nPairs < 3, return; end
    % bootstrap: every green mouse with every red mouse -> resample both sets of mice independently and keep
    % all their mixtures; otherwise resample mixtures
    gN = {R.pairs.green}; rN = {R.pairs.red};
    ug = unique(gN); ur = unique(rN);
    M = zeros(numel(ug), numel(ur));
    for p = 1:nP
        if any(Z.pair == p & ok), M(strcmp(ug, gN{p}), strcmp(ur, rN{p})) = p; end
    end
    full = all(M(:) > 0) && numel(ug) >= 2 && numel(ur) >= 2;
    pl = unique(Z.pair(ok));
    rowsOf = cell(nP, 1);
    for p = 1:nP, rowsOf{p} = find(Z.pair == p & ok); end
    st0 = [];
    try, st0 = rng; rng(1, 'twister'); catch, end %#ok<NOCOM>
    B = nan(P.nBoot, 2);
    for bb = 1:P.nBoot
        if full
            sel = M(randi(numel(ug), numel(ug), 1), randi(numel(ur), numel(ur), 1));
        else
            sel = pl(randi(numel(pl), numel(pl), 1));
        end
        idx = vertcat(rowsOf{sel(:)});
        [B(bb, 1), B(bb, 2)] = ratioFit(Z.x(idx), Z.y(idx), Z.model);
    end
    if ~isempty(st0), try, rng(st0); catch, end, end %#ok<NOCOM>
    Z.cCI = qCI(B(:, 1)); Z.bCI = qCI(B(:, 2)); Z.nBoot = P.nBoot;
    if full
        Z.bootUnit = sprintf('animals (%d green x %d red mice, resampled independently)', numel(ug), numel(ur));
    else
        Z.bootUnit = sprintf('mixtures (%d)', numel(pl));
    end
end


function g = ratioG(x, model)
% regressor of the fit: x (varied sensor: partner / itself = x) or 1 / x (unchanged sensor)
    if strcmp(model, 'lin'), g = x; else, g = 1 ./ x; end
end


function t = ratioModelText(model)
    if strcmp(model, 'lin'), t = 'difference = c + b \cdot x'; else, t = 'difference = c + b / x'; end
end


function [c, b, r2] = ratioFit(x, y, model)
% least squares y = c + b * g(x)
    c = NaN; b = NaN; r2 = NaN;
    x = x(:); y = y(:);
    if numel(x) < 3 || numel(unique(x)) < 2, return; end
    X = [ones(numel(x), 1), ratioG(x, model)];
    cf = X \ y; c = cf(1); b = cf(2);
    ss = sum((y - mean(y)) .^ 2);
    if ss > 0, r2 = 1 - sum((y - X * cf) .^ 2) / ss; end
end


function [yh, se] = ratioBand(x, y, xq, model)
% fitted y = c + b * g(x) at xq and its standard error (SEM of the fitted mean; n = numel(x))
    x = x(:); y = y(:); xq = xq(:);
    X = [ones(numel(x), 1), ratioG(x, model)];
    cf = X \ y;
    s2 = sum((y - X * cf) .^ 2) / max(numel(x) - 2, 1);
    V = s2 * inv(X' * X); %#ok<MINV>
    Xq = [ones(numel(xq), 1), ratioG(xq, model)];
    yh = Xq * cf;
    se = sqrt(max(sum((Xq * V) .* Xq, 2), 0));
end


function ci = qCI(v)
% 2.5th and 97.5th percentiles (no toolbox)
    v = sort(v(isfinite(v))); n = numel(v);
    if n < 20, ci = [NaN NaN]; return; end
    ci = [v(max(1, round(0.025 * n))), v(min(n, round(0.975 * n)))];
end


function [lim, st] = ratioAxis(P, fld, idx, v)
% axis from the settings (P.(fld){idx}) or automatic
    a = [];
    if isfield(P, fld)
        C = P.(fld);
        if numel(idx) == 1 && numel(C) >= idx, a = C{idx};
        elseif numel(idx) == 2 && size(C, 1) >= idx(1) && size(C, 2) >= idx(2), a = C{idx(1), idx(2)};
        end
    end
    if numel(a) == 3, lim = a(1:2); st = a(3); return; end
    v = v(isfinite(v));
    [lim, st] = niceAxis(v);
    if min(v) >= 0, lim(1) = 0; end              % no negative margin when every value is >= 0
end


function ratioFigure(R, d, P, fileBase)
% one file per experiment (varied sensor d); left: difference of the varied sensor, right: difference of
% the unchanged sensor, both vs the same input ratio x = total GT AUC unchanged (1x) / varied (nx).
% Light dots = every mixture and factor (no grouping); black = least squares fit over all points
% (varied: c + b x; unchanged: c + b / x); grey band = +/- SEM of the fitted curve
    X = char(215);
    tName = {'T1', 'T2'};
    o = 3 - d; nmO = tx(R.names{o}); nmD = tx(R.names{d});
    ms = P.dotSize; if isfield(P, 'ratioDotSize'), ms = P.ratioDotSize; end
    shown = [d, o];
    Zs = cell(1, 2);
    for t = 1:2, Zs{t} = ratioData(R, d, shown(t), P); end
    if Zs{1}.n < 1 && Zs{2}.n < 1, return; end
    xs = Zs{2}.x(Zs{2}.ok);
    [xl, sx] = ratioAxis(P, 'axRatioX', d, [0; xs]);
    pw = 4.6; ph = 4.6; gx = 2.2; x0 = 1.8; y0 = 1.7; top = 0.4 + 0.38 * 4;
    f = newFigure(x0 + 2 * pw + gx + 0.6, y0 + ph + top);
    for t = 1:2
        Z = Zs{t}; s = Z.s;
        if Z.n < 1, continue; end
        nmS = tx(R.names{s});
        col = R.colors{s}; light = col + (1 - col) * P.pairTint;
        ax = newAxes(f, [x0 + (t - 1) * (pw + gx), y0, pw, ph], P);
        xv = Z.x(Z.ok); yv = Z.y(Z.ok);
        [yl, sy] = ratioAxis(P, 'axRatioY', [d s], [0; yv]);
        if yl(1) < 0 && yl(2) > 0
            plot(ax, xl, [0 0], '--', 'Color', [0 0 0], 'LineWidth', 0.6);
        end
        plot(ax, xv, yv, 'o', 'MarkerSize', ms, 'MarkerFaceColor', light, 'MarkerEdgeColor', 'none', 'LineStyle', 'none');
        if isfinite(Z.b)
            xq = linspace(max(min(xv), xl(1) + eps), min(max(xv), xl(2)), 200)';
            [yh, se] = ratioBand(xv, yv, xq, Z.model);
            fill(ax, [xq; flipud(xq)], [yh + se; flipud(yh - se)], [0.72 0.72 0.72], 'EdgeColor', 'none');
            plot(ax, xq, yh, '-', 'Color', [0 0 0], 'LineWidth', 1.2);
        end
        set(ax, 'XLim', xl, 'XTick', limTicks(xl, sx), 'YLim', yl, 'YTick', limTicks(yl, sy));
        out = sum(xv < xl(1) | xv > xl(2) | yv < yl(1) | yv > yl(2));
        if out > 0
            fprintf('   [note] %s varied, %s difference vs ratio: %d of %d point(s) outside the axes (P.axRatioX / P.axRatioY)\n', ...
                R.names{d}, R.names{s}, out, Z.n);
        end
        if strcmp(Z.xKind, 'total')
            xlabel(ax, {'Input ratio (total GT AUC)', sprintf('%s 1%s / %s {\\itn}%s', nmO, X, nmD, X)});
        else
            xlabel(ax, {'GT AUC ratio', sprintf('%s %s 1%s / %s %s {\\itn}%s', nmO, tName{o}, X, nmD, tName{d}, X)});
        end
        ylabel(ax, sprintf('%s %s difference (%%)', nmS, tName{s}));
        hdr = {sprintf('\\color[rgb]{%.3f %.3f %.3f}%s %s', col, nmS, tName{s}), ...
            ['\color[rgb]{0 0 0}' ratioModelText(Z.model)], ...
            sprintf('\\color[rgb]{0 0 0}b = %.3g %%,  c = %.3g %%', Z.b, Z.c), ...
            sprintf('\\color[rgb]{0 0 0}R^2 = %.2f  (n = %d)', Z.r2, Z.n)};
        text(ax, 0, 1.03, hdr, 'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontName', P.fontName, ...
            'FontSize', P.fontSize);
        fprintf('   %s varied: %s %s %s   b = %.4g [%.4g, %.4g] %%   c = %.4g [%.4g, %.4g] %%   R^2 = %.3f   (n = %d)\n', ...
            R.names{d}, R.names{s}, tName{s}, strrep(ratioModelText(Z.model), '\cdot', '*'), Z.b, Z.bCI(1), Z.bCI(2), ...
            Z.c, Z.cCI(1), Z.cCI(2), Z.r2, Z.n);
    end
    saveFigure(f, fileBase);
end


function writeRatioSheets(R, P, xlsFile)
% DiffRatio: every mixture and factor; DiffRatioFit: fit per experiment and sensor shown
    tName = {'T1', 'T2'};
    C = {'Varied sensor', 'Sensor shown', 'Green mouse', 'Red mouse', 'Factor of the varied sensor', ...
         'Total GT AUC unchanged (1x)', 'Total GT AUC varied (nx)', 'Input ratio x = unchanged / varied (total AUC)', ...
         'GT T1/T2 AUC of the sensor shown', 'Unmixed T1/T2 AUC of the sensor shown', 'Difference (%)', 'Model', ...
         'Fitted (%)'};
    S = {'Varied sensor', 'Sensor shown', 'Difference of', 'Model', 'n points', 'n mixtures', 'b (%)', 'b 95% CI lower', ...
         'b 95% CI upper', 'c (%)', 'c 95% CI lower', 'c 95% CI upper', 'R^2', 'Bootstrap resamples', 'Bootstrap unit'};
    for d = 1:2
        for s = [d, 3 - d]
            Z = ratioData(R, d, s, P);
            mt = strrep(ratioModelText(Z.model), '\cdot', '*');
            for i = 1:numel(Z.x)
                C(end+1, :) = {R.names{d}, R.names{s}, Z.green{i}, Z.red{i}, Z.k(i), Z.totU(i), Z.totV(i), Z.x(i), ...
                    Z.gtSh(i), Z.unSh(i), Z.y(i), mt, Z.fitted(i)}; %#ok<AGROW>
            end
            S(end+1, :) = {R.names{d}, R.names{s}, [tName{s} ' AUC; x = ' Z.xKind ' AUC'], mt, Z.n, Z.nPairs, Z.b, ...
                Z.bCI(1), Z.bCI(2), Z.c, Z.cCI(1), Z.cCI(2), Z.r2, Z.nBoot, Z.bootUnit}; %#ok<AGROW>
        end
    end
    writecell(clean(C), xlsFile, 'Sheet', 'DiffRatio');
    writecell(clean(S), xlsFile, 'Sheet', 'DiffRatioFit');
end


%% =====================================================================
%  T1 / T2 window AUC vs total AUC (one figure per sensor): does the window AUC change by the
%  same ratio as the total AUC when the brightness is changed (power decoder)?
%    spectra: time-averaged ground-truth spectra of every animal at every factor (R.spectra)
%    x = total AUC (all wavelengths) at x n / at x1; y = T1 or T2 AUC at x n / at x1 (same animal)
%    T1 / T2 window: +/- aucHalfWidth around the peak of the x1 spectrum in bandT1 / bandT2,
%    the same window at every factor
%    slope through (1, 1) per animal; across animals: mean [95% CI], t-test vs 1
%  =====================================================================
function W = windowData(R, s)
    fac = R.factors(:)'; nF = numel(fac);
    i1 = find(abs(fac - 1) < 1e-9, 1); kk = setdiff(1:nF, i1);
    Sp = R.spectra([R.spectra.dir] == s);
    st = R.settings; bands = {st.bandT1, st.bandT2};
    nA = numel(Sp);
    W.s = s; W.fac = fac; W.i1 = i1; W.kk = kk; W.animals = {Sp.animal};
    W.tot = nan(nA, nF); W.win = nan(nA, nF, 2); W.centre = nan(nA, 2);
    for a = 1:nA
        wl = Sp(a).wl(:); S = Sp(a).S;
        y1 = smoothNm(wl, S(:, i1), st.peakSmoothNm);
        for t = 1:2, W.centre(a, t) = bandPeak(wl, y1, bands{t}); end
        for f = 1:nF
            W.tot(a, f) = trapz(wl, S(:, f));
            for t = 1:2, W.win(a, f, t) = winAuc(wl, S(:, f), W.centre(a, t), st.aucHalfWidth); end
        end
    end
    W.xr = bsxfun(@rdivide, W.tot, W.tot(:, i1));
    W.yr = bsxfun(@rdivide, W.win, W.win(:, i1, :));
    W.slope = nan(nA, 2);
    for t = 1:2
        for a = 1:nA
            dx = W.xr(a, kk) - 1; dy = W.yr(a, kk, t) - 1;
            ok = isfinite(dx) & isfinite(dy);
            if any(ok) && sum(dx(ok) .^ 2) > 0
                W.slope(a, t) = sum(dx(ok) .* dy(ok)) / sum(dx(ok) .^ 2);
            end
        end
        W.stat(t) = oneSample(W.slope(:, t), 1);
        % all animals and factors pooled: y = a + b x (least squares); mean per factor
        xx = W.xr; yy = W.yr(:, :, t);
        [cf, r2] = linFit(xx(:), yy(:));
        W.reg(t) = struct('a', cf(1), 'b', cf(2), 'r2', r2, 'n', sum(isfinite(xx(:)) & isfinite(yy(:))), ...
            'xMin', min(xx(:)), 'xMax', max(xx(:)));
        W.mx(:, t) = mean(xx, 1, 'omitnan')'; W.my(:, t) = mean(yy, 1, 'omitnan')';
    end
end


function windowFigure(R, s, P, fileBase)
% every animal and factor as light dots (no grouping); dark dots = mean per factor; dashed y = x;
% slope and R^2: least squares y = a + b x over all points
    if ~isfield(R, 'spectra') || isempty(R.spectra) || ~any([R.spectra.dir] == s), return; end
    W = windowData(R, s);
    col = R.colors{s}; light = col + (1 - col) * P.pairTint;
    tName = {'T1', 'T2'}; nm = tx(R.names{s});
    pw = 4.6; ph = 4.6; gx = 2.0; x0 = 1.8; y0 = 1.5; top = 0.4 + 0.38 * 3;
    f = newFigure(x0 + 2 * pw + gx + 0.6, y0 + ph + top);
    for t = 1:2
        ax = newAxes(f, [x0 + (t - 1) * (pw + gx), y0, pw, ph], P);
        yr = W.yr(:, :, t);
        [lim, st] = niceAxis([0; W.xr(:); yr(:)]);
        lim(1) = 0;
        plot(ax, lim, lim, '--', 'Color', [0 0 0], 'LineWidth', 0.8);
        plot(ax, W.xr(:), yr(:), 'o', 'MarkerSize', P.dotSize, 'MarkerFaceColor', light, 'MarkerEdgeColor', 'none', ...
            'LineStyle', 'none');
        plot(ax, W.mx(:, t), W.my(:, t), 'o', 'MarkerSize', P.dotSize, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', ...
            'LineStyle', 'none');
        set(ax, 'XLim', lim, 'YLim', lim, 'XTick', lim(1):st:lim(2), 'YTick', lim(1):st:lim(2));
        xlabel(ax, sprintf('Total AUC ratio ({\\itn}%s / 1%s)', char(215), char(215)));
        ylabel(ax, sprintf('%s AUC ratio ({\\itn}%s / 1%s)', tName{t}, char(215), char(215)));
        Rg = W.reg(t);
        hdr = {sprintf('\\color[rgb]{%.3f %.3f %.3f}%s %s', col, nm, tName{t}), ...
            sprintf('\\color[rgb]{0 0 0}slope = %.3f', Rg.b), sprintf('\\color[rgb]{0 0 0}R^2 = %.4f', Rg.r2)};
        text(ax, 0, 1.03, hdr, 'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontName', P.fontName, ...
            'FontSize', P.fontSize);
        fprintf('   %s %s AUC vs total AUC: y = %.4f + %.4f x, R^2 = %.4f (n = %d)\n', R.names{s}, tName{t}, ...
            Rg.a, Rg.b, Rg.r2, Rg.n);
    end
    saveFigure(f, fileBase);
end


function writeWindowSheets(R, xlsFile)
    if ~isfield(R, 'spectra') || isempty(R.spectra), return; end
    C = {'Sensor', 'Animal', 'Factor', 'Total AUC', 'T1 window centre (nm)', 'T1 AUC', 'T2 window centre (nm)', ...
         'T2 AUC', 'Total AUC ratio (x n / x1)', 'T1 AUC ratio', 'T2 AUC ratio', 'T1 / total', 'T2 / total'};
    S = {'Sensor', 'Window', 'Animal', 'Slope through (1,1)', 'SD', '95% CI lower', '95% CI upper', 't (vs 1)', 'df', 'p (vs 1)'};
    L = {'Sensor', 'Window', 'n points (animals x factors)', 'Slope', 'Intercept', 'R^2'};
    L2 = {'Sensor', 'Window', 'Factor', 'Mean total AUC ratio', 'Mean T1/T2 AUC ratio'};
    tName = {'T1', 'T2'};
    for s = 1:2
        if ~any([R.spectra.dir] == s), continue; end
        W = windowData(R, s);
        for a = 1:numel(W.animals)
            for j = 1:numel(W.fac)
                C(end+1, :) = {R.names{s}, W.animals{a}, W.fac(j), W.tot(a, j), W.centre(a, 1), W.win(a, j, 1), ...
                    W.centre(a, 2), W.win(a, j, 2), W.xr(a, j), W.yr(a, j, 1), W.yr(a, j, 2), ...
                    W.win(a, j, 1) / W.tot(a, j), W.win(a, j, 2) / W.tot(a, j)}; %#ok<AGROW>
            end
        end
        for t = 1:2
            for a = 1:numel(W.animals)
                S(end+1, :) = {R.names{s}, tName{t}, W.animals{a}, W.slope(a, t), '', '', '', '', '', ''}; %#ok<AGROW>
            end
            T = W.stat(t);
            S(end+1, :) = {R.names{s}, tName{t}, sprintf('mean (n = %d)', T.n), T.mean, T.sd, T.ciLo, T.ciHi, ...
                T.t, T.df, T.p}; %#ok<AGROW>
            Rg = W.reg(t);
            L(end+1, :) = {R.names{s}, tName{t}, Rg.n, Rg.b, Rg.a, Rg.r2}; %#ok<AGROW>
            for j = 1:numel(W.fac)
                L2(end+1, :) = {R.names{s}, tName{t}, W.fac(j), W.mx(j, t), W.my(j, t)}; %#ok<AGROW>
            end
        end
    end
    writecell(clean(C), xlsFile, 'Sheet', 'WindowScaling');
    writecell(clean([L; cell(1, size(L, 2)); [L2, cell(size(L2, 1), size(L, 2) - size(L2, 2))]]), xlsFile, ...
        'Sheet', 'WindowScalingFit');
    writecell(clean(S), xlsFile, 'Sheet', 'WindowScalingStats');
end


%% =====================================================================
%  Supplementary figures
%    spectrum\  every mixture and condition: mixture, green and red spectra (ground truth vs unmixed)
%    unmixing\  every condition: agreement and Bland-Altman of both sensors (as NMF_Hybrid_Equivalence)
%  =====================================================================
function supplementFigures(R, P, base)
    outDir = fileparts(base);
    pre = safeFileName(sprintf('%s+%s', R.names{1}, R.names{2}));
    if P.supSpectra
        supSpectra(R, P, fullfile(outDir, 'spectrum'), pre);
    end
    supUnmixing(R, P, fullfile(outDir, 'unmixing'), pre);
end


function s = condName(R, u)
% file-name label of computed condition u, e.g. GCaMPx2_jRGECOx1
    c = R.pairs(1).cond(u, :);
    s = safeFileName(sprintf('%sx%g_%sx%g', R.names{1}, c(1), R.names{2}, c(2)));
end


function t = condTitle(R, u)
% coloured title of condition u, e.g. GCaMP 2x + jRGECO 1x
    c = R.pairs(1).cond(u, :); cl = R.colors;
    x = char(215);                                  % multiplication sign (not TeX \times: font substitution in PDF)
    t = sprintf(['\\color[rgb]{%.3f %.3f %.3f}%s %g%s\\color[rgb]{0 0 0} + ' ...
        '\\color[rgb]{%.3f %.3f %.3f}%s %g%s\\color[rgb]{0 0 0}'], cl{1}, tx(R.names{1}), c(1), x, ...
        cl{2}, tx(R.names{2}), c(2), x);
end


function v = r2Trials(T, Pr)
% R^2 of every trial (columns) of the ground truth T vs the prediction Pr; [mean SEM n]
    n = size(T, 2);
    r2 = nan(n, 1);
    for k = 1:n
        y = T(:, k);
        ss = sum((y - mean(y)) .^ 2);
        if ss > 0, r2(k) = 1 - sum((y - Pr(:, k)) .^ 2) / ss; end
    end
    r2 = r2(isfinite(r2));
    if isempty(r2), v = [NaN NaN 0]; return; end
    se = 0; if numel(r2) > 1, se = std(r2) / sqrt(numel(r2)); end
    v = [mean(r2), se, numel(r2)];
end


% ----------------------------------------------------------------------- spectrum\
function supSpectra(R, P, folder, pre)
    if ~isfield(R.pairs, 'spec') || isempty(R.pairs(1).spec)
        fprintf(['\nSpectrum figures: this .mat file has no mixture spectra (saved by an older version); ' ...
            'run the analysis again to make them.\n']);
        return;
    end
    if exist(folder, 'dir') ~= 7, mkdir(folder); end
    nU = size(R.pairs(1).cond, 1); nP = numel(R.pairs);
    % highest mixture value of all files (450-650 nm): automatic common y axis, or check of P.specYLim
    top = 0;
    for p = 1:nP
        w = R.pairs(p).wl(:); in = w >= P.specXLim(1) & w <= P.specXLim(2);
        for u = 1:nU
            S = R.pairs(p).spec(u);
            top = max([top; S.truthG(in) + S.truthR(in); max(S.compG(in), 0) + max(S.compR(in), 0)]);
        end
    end
    if isempty(P.specYLim)
        st = niceStep(top); P.specYTick = st;
        P.specYLim = [-0.01 * ceil(top / st) * st, ceil(top / st) * st];
        fprintf('Spectrum y axis (common to all files): %g to %g counts\n', P.specYLim);
    elseif top > P.specYLim(2)
        fprintf(['   [note] the highest mixture value is %.0f counts, above the y axis limit %g (P.specYLim): ' ...
            'peaks are cut in some files. Set P.specYLim = [] (automatic) or a higher limit.\n'], top, P.specYLim(2));
    end
    fprintf('\nSpectrum figures: %d mixtures x %d conditions -> %s\n', nP, nU, folder);
    t0 = tic;
    for p = 1:nP
        for u = 1:nU
            fb = fullfile(folder, safeFileName(sprintf('%s_%s_%s', R.pairs(p).green, R.pairs(p).red, condName(R, u))));
            spectrumFigure(R, p, u, P, fb);
        end
        fprintf('   %s + %s: %d figures (%.0f s)\n', R.pairs(p).green, R.pairs(p).red, nU, toc(t0));
    end
    % R^2 table
    C = {'Green mouse', 'Red mouse', 'Condition', sprintf('%s factor', R.names{1}), sprintf('%s factor', R.names{2}), ...
         'R^2 mixture (mean)', 'R^2 mixture (SEM)', sprintf('R^2 %s (mean)', R.names{1}), sprintf('R^2 %s (SEM)', R.names{1}), ...
         sprintf('R^2 %s (mean)', R.names{2}), sprintf('R^2 %s (SEM)', R.names{2}), 'n trials'};
    for p = 1:nP
        for u = 1:nU
            S = R.pairs(p).spec(u); c = R.pairs(p).cond(u, :);
            C(end+1, :) = {R.pairs(p).green, R.pairs(p).red, condName(R, u), c(1), c(2), S.r2Mix(1), S.r2Mix(2), ...
                S.r2G(1), S.r2G(2), S.r2R(1), S.r2R(2), S.r2Mix(3)}; %#ok<AGROW>
        end
    end
    M = {'Condition', sprintf('%s factor', R.names{1}), sprintf('%s factor', R.names{2}), 'n mixtures', ...
         'R^2 mixture (mean over mixtures)', 'SEM', sprintf('R^2 %s', R.names{1}), 'SEM', sprintf('R^2 %s', R.names{2}), 'SEM'};
    for u = 1:nU
        v = zeros(nP, 3);
        for p = 1:nP
            S = R.pairs(p).spec(u); v(p, :) = [S.r2Mix(1), S.r2G(1), S.r2R(1)];
        end
        c = R.pairs(1).cond(u, :);
        row = {condName(R, u), c(1), c(2), nP};
        for q = 1:3
            w = v(isfinite(v(:, q)), q);
            se = NaN; if numel(w) > 1, se = std(w) / sqrt(numel(w)); end
            row = [row, {mean(w), se}]; %#ok<AGROW>
        end
        M(end+1, :) = row; %#ok<AGROW>
    end
    xls = fullfile(folder, [pre '_Spectrum_R2.xlsx']);
    if exist(xls, 'file') == 2, try, delete(xls); catch, end, end %#ok<NOCOM>
    try
        writecell(clean(M), xls, 'Sheet', 'Summary');
        writecell(clean(C), xls, 'Sheet', 'Mixtures');
        fprintf('   saved %s\n', xls);
    catch ME
        fprintf('   %s not saved: %s\n', xls, ME.message);
    end
end


function spectrumFigure(R, p, u, P, fileBase)
% left: mixture (counts), black = ground truth, fills = NMF-RI components stacked;
% middle / right: ground truth (black) and unmixed (fill) of each sensor, normalised to the
% ground-truth peak; R^2 per trial (mean +/- SEM)
    Q = R.pairs(p); S = Q.spec(u); wl = double(Q.wl(:)); cl = R.colors;
    pw = P.supPanelCm(1); ph = P.supPanelCm(2); gx = 1.6; x0 = 1.6; y0 = 1.2;
    W = x0 + 3 * pw + 2 * gx + 0.4; H = y0 + ph + 2.6;
    f = newFigure(W, H);
    ov = axes(f, 'Units', 'centimeters', 'Position', [0 0 W H], 'XLim', [0 W], 'YLim', [0 H], 'Visible', 'off');
    hold(ov, 'on');
    text(ov, 0.3, H - 0.25, sprintf('%s + %s:  %s', tx(Q.green), tx(Q.red), condTitle(R, u)), 'FontName', P.fontName, ...
        'FontSize', P.fontSize, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    z = zeros(size(wl));
    xt = P.specXLim(1):P.specXTick:P.specXLim(2);
    r2txt = @(v) sprintf('R^2 = %.4f %s %.4f', v(1), char(177), v(2));
    % ---- mixture
    ax = newAxes(f, [x0, y0, pw, ph], P);
    y1 = max(double(S.compG(:)), 0); y12 = y1 + max(double(S.compR(:)), 0);
    mix = double(S.truthG(:) + S.truthR(:));
    fill(ax, [wl; flipud(wl)], [y12; flipud(y1)], cl{2}, 'EdgeColor', 'none');
    fill(ax, [wl; flipud(wl)], [y1; z], cl{1}, 'EdgeColor', 'none');
    plot(ax, wl, mix, 'k-', 'LineWidth', P.specLineWidth);
    in = wl >= P.specXLim(1) & wl <= P.specXLim(2);
    yl = P.specYLim; st = P.specYTick;
    set(ax, 'XLim', P.specXLim, 'XTick', xt, 'YLim', yl, 'YTick', 0:st:yl(2));
    xlabel(ax, 'wavelength (nm)'); ylabel(ax, 'intensity (counts)');
    specLines(ax, {condTitle(R, u), 'NMF unmixing', r2txt(S.r2Mix)}, P);
    % ---- each sensor
    tr = {S.truthG, S.truthR}; cp = {S.compG, S.compR}; r2 = {S.r2G, S.r2R}; bnd = {'T1', 'T2'};
    for s = 1:2
        ax = newAxes(f, [x0 + s * (pw + gx), y0, pw, ph], P);
        t = double(tr{s}(:)); c = double(cp{s}(:));
        pk = max(t(in)); if ~(pk > 0), pk = 1; end
        fill(ax, [wl; flipud(wl)], [max(c, 0) / pk; z], cl{s}, 'EdgeColor', 'none');
        plot(ax, wl, t / pk, 'k-', 'LineWidth', P.specLineWidth);
        if s == 1, band = R.settings.bandT1; else, band = R.settings.bandT2; end
        mb = wl >= band(1) & wl <= band(2);
        if any(mb)
            [mx, im] = max(t(mb)); wb = wl(mb);
            text(ax, wb(im), mx / pk + 0.04, bnd{s}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontName', P.fontName, 'FontSize', P.fontSize);
        end
        set(ax, 'XLim', P.specXLim, 'XTick', xt, 'YLim', [0 1.2], 'YTick', 0:0.5:1);
        xlabel(ax, 'wavelength (nm)'); ylabel(ax, 'normalized intensity');
        nm = tx(R.names{s}); k = R.pairs(1).cond(u, s);
        specLines(ax, {sprintf('GT of %s %g%s', nm, k, char(215)), ...
            sprintf('\\color[rgb]{%.3f %.3f %.3f}unmixed %s', cl{s}, nm), r2txt(r2{s})}, P);
    end
    saveFigure2(f, fileBase, P.supSpectrumTypes, false);
end


function specLines(ax, lines, P)
% text lines above the axes (one text object per line)
    n = numel(lines);
    for j = 1:n
        text(ax, 0.02, 1.03 + 0.1 * (n - j), lines{j}, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
            'FontName', P.fontName, 'FontSize', P.fontSize, 'Clipping', 'off');
    end
end


function saveFigure2(f, fileBase, types, verbose)
% as saveFigure, chosen file types, optionally quiet
    for e = 1:numel(types)
        file = [fileBase types{e}];
        if exist(file, 'file') == 2
            try, delete(file); catch, end %#ok<NOCOM>
        end
        try
            if exist('exportgraphics') > 0 %#ok<EXIST>
                exportgraphics(f, file, 'ContentType', 'vector', 'BackgroundColor', 'white');
            elseif strcmpi(types{e}, '.pdf')
                print(f, file, '-dpdf', '-painters');
            else
                print(f, file, '-depsc', '-painters');
            end
            if verbose, fprintf('   saved %s\n', file); end
        catch ME
            fprintf('   %s not saved: %s\n', file, ME.message);
        end
    end
    close(f);
end


% ----------------------------------------------------------------------- unmixing\
function supUnmixing(R, P, folder, pre)
    if exist(folder, 'dir') ~= 7, mkdir(folder); end
    nU = size(R.pairs(1).cond, 1);
    band = {'T1', 'T2'};
    fprintf('\nUnmixing figures (agreement and Bland-Altman): %d conditions -> %s\n', nU, folder);
    U = struct('u', {}, 'name', {}, 'cond', {}, 'S', {});
    for u = 1:nU
        S = struct('slot', {}, 'name', {}, 'color', {}, 'band', {}, 'auc', {}, 'dauc', {}, 'ccc', {}, 'ba', {}, ...
            'tostAUC', {}, 'tostDAUC', {}, 'reg', {});
        for s = 1:2
            D = condData(R, u, s, false);
            DD = condData(R, u, s, true);
            S(s).slot = s; S(s).name = R.names{s}; S(s).color = R.colors{s}; S(s).band = band{s};
            S(s).auc = D; S(s).dauc = DD;
            S(s).ccc = lincc(D.animalRaw, D.animalNmf);
            S(s).ba = blandAltman(D.animalPct);
            S(s).tostAUC = tost(D.animalPct, P.margin, P.alpha);
            S(s).tostDAUC = tost(DD.animalPct, P.margin, P.alpha);
            S(s).reg = baRegression(D.animalRaw, D.animalPct);
        end
        U(u).u = u; U(u).name = condName(R, u); U(u).cond = R.pairs(1).cond(u, :); U(u).S = S;
    end
    % axes: the same within each varied sensor (green x1-x3 with red x1; red x1-x3 with green x1) and
    % the same as the _AUC_<sensor>_varied figure; x1 + x1 belongs to both series and is drawn for each
    AX = axisGroups(R, P, true);
    for d = 1:2
        us = groupConditions(R, d);
        for u = us(:)'
            fn = U(u).name;
            if all(abs(U(u).cond - 1) < 1e-9), fn = sprintf('%s_%svaried', fn, safeFileName(R.names{d})); end
            unmixingFigure(U(u), R, P, fullfile(folder, fn), AX(d, :));
            S = U(u).S;
            fprintf('   %-30s', fn);
            for s = 1:2
                fprintf('  %s: CCC %.4f [%.4f, %.4f], bias %+.2f %% (LoA %+.2f to %+.2f)', S(s).name, S(s).ccc.ccc, ...
                    S(s).ccc.lo, S(s).ccc.hi, S(s).ba.bias, S(s).ba.lo, S(s).ba.hi);
            end
            fprintf('\n');
        end
    end
    writeUnmixingStats(U, R, P, fullfile(folder, [pre '_Unmixing_stats.xlsx']));
end


function us = groupConditions(R, d)
% computed conditions of the series in which sensor d is varied (x1 ... largest factor), in order
    rows = find(R.cond(:, 3) == d);
    [~, o] = sort(R.cond(rows, 4));
    us = R.condIdx(rows(o));
    us = us(:)';
end


function AX = axisGroups(R, P, verbose)
% Axes shared by all figures of one varied sensor d (rows) for sensor s (columns):
%   lim, st   ground-truth / unmixed AUC axis (x 10^3 counts nm), agreement x and y, Bland-Altman x,
%             and panels B and C of <base>_AUC_<sensor>_varied: covers every mixture of the series
%   baY, baSt difference axis (%), Bland-Altman of unmixing\ and panel C of _AUC_<sensor>_varied: every
%             mixture, the limits of agreement and the regression lines of the series; 0-1 % (tick 0.2)
%             when all of them lie within 0 to 1 %
    AX = struct('lim', {}, 'st', {}, 'baY', {}, 'baSt', {});
    for d = 1:2
        us = groupConditions(R, d);
        for s = 1:2
            v = 0; b = 0;
            for u = us
                D = condData(R, u, s, false);
                B = blandAltman(D.animalPct);
                v = [v; D.pairRaw / 1e3; D.pairNmf / 1e3]; %#ok<AGROW>
                b = [b; D.pairPct; B.lo; B.hi]; %#ok<AGROW>
            end
            Q = aucData(R, d, s);                      % regression lines of panel C
            for k = 1:numel(Q.fit)
                F = Q.fit(k);
                if isfinite(F.b), b = [b; F.a + F.b * [F.xMin; F.xMax]]; end %#ok<AGROW>
            end
            [lim, st] = niceAxis(v); lim(1) = 0;
            bb = b(2:end); bb = bb(isfinite(bb));
            if ~isempty(bb) && min(bb) >= 0 && max(bb) <= 1
                yl = [0 1]; ys = 0.2;                  % small differences: 0-1 %
            else
                [yl, ys] = baLimits(b, P.baMinSpan);
            end
            % set axes (P.axAUC, P.axDiff); report values outside them (not shown)
            nm = sprintf('%s varied, %s', R.names{d}, R.names{s});
            if size(P.axAUC, 1) >= d && size(P.axAUC, 2) >= s && numel(P.axAUC{d, s}) == 3
                lim = P.axAUC{d, s}(1:2); st = P.axAUC{d, s}(3);
                n = sum(v(2:end) < lim(1) | v(2:end) > lim(2));
                if n > 0 && verbose
                    fprintf(['   [note] %s: %d AUC value(s) outside the axis %g-%g (P.axAUC; range of the data ' ...
                        '%.3g-%.3g x10^3)\n'], nm, n, lim, min(v(2:end)), max(v(2:end)));
                end
            end
            if size(P.axDiff, 1) >= d && size(P.axDiff, 2) >= s && numel(P.axDiff{d, s}) == 3
                yl = P.axDiff{d, s}(1:2); ys = P.axDiff{d, s}(3);
                n = sum(bb < yl(1) | bb > yl(2));
                if n > 0 && verbose
                    fprintf(['   [note] %s: %d difference value(s) (mixtures, LoA or lines) outside the axis %g to %g %%' ...
                        ' (P.axDiff; range of the data %.3g to %.3g %%)\n'], nm, n, yl, min(bb), max(bb));
                end
            end
            AX(d, s) = struct('lim', lim, 'st', st, 'baY', yl, 'baSt', ys);
        end
    end
end


function D = condData(R, u, s, delta)
% mixture values and mouse-level means of sensor s (1 green T1, 2 red T2) in condition u
    n = numel(R.pairs);
    D.pairRaw = nan(n, 1); D.pairNmf = nan(n, 1); D.pairAnimal = cell(n, 1);
    for i = 1:n
        Q = R.pairs(i);
        if delta
            D.pairRaw(i) = Q.rawDelta(u, s); D.pairNmf(i) = Q.nmfDelta(u, s);
        else
            D.pairRaw(i) = Q.raw(u, s); D.pairNmf(i) = Q.nmf(u, s);
        end
        if s == 1, D.pairAnimal{i} = Q.green; else, D.pairAnimal{i} = Q.red; end
    end
    ok = isfinite(D.pairRaw) & isfinite(D.pairNmf) & D.pairRaw ~= 0;
    D.pairRaw = D.pairRaw(ok); D.pairNmf = D.pairNmf(ok); D.pairAnimal = D.pairAnimal(ok);
    D.pairPct = 100 * (D.pairNmf - D.pairRaw) ./ D.pairRaw;
    D.animal = {};
    g = zeros(numel(D.pairAnimal), 1);
    for i = 1:numel(D.pairAnimal)
        j = find(strcmp(D.animal, D.pairAnimal{i}), 1);
        if isempty(j), D.animal{end+1, 1} = D.pairAnimal{i}; j = numel(D.animal); end
        g(i) = j;
    end
    na = numel(D.animal);
    D.animalRaw = zeros(na, 1); D.animalNmf = zeros(na, 1); D.animalPct = zeros(na, 1);
    D.animalNPairs = zeros(na, 1); D.animalSdNmf = zeros(na, 1);
    for a = 1:na
        w = g == a;
        D.animalRaw(a) = mean(D.pairRaw(w)); D.animalNmf(a) = mean(D.pairNmf(w));
        D.animalPct(a) = mean(D.pairPct(w)); D.animalNPairs(a) = sum(w);
        if sum(w) > 1, D.animalSdNmf(a) = std(D.pairNmf(w)); end
    end
end


function unmixingFigure(Uu, R, P, fileBase, AXd)
% one row: green agreement, green Bland-Altman, red agreement, red Bland-Altman
    pw = P.supPanelCm(1); ph = P.supPanelCm(2);
    gIn = 1.7; gOut = 2.2; x0 = 1.6; y0 = 2.5;
    W = x0 + 4 * pw + 2 * gIn + gOut + 0.4; H = y0 + ph + 2.7;
    f = newFigure(W, H);
    ov = axes(f, 'Units', 'centimeters', 'Position', [0 0 W H], 'XLim', [0 W], 'YLim', [0 H], 'Visible', 'off');
    hold(ov, 'on');
    text(ov, 0.3, H - 0.25, condTitle(R, Uu.u), 'FontName', P.fontName, 'FontSize', P.fontSize + 1, ...
        'FontWeight', 'bold', 'VerticalAlignment', 'top');
    for s = 1:2
        xa = x0 + (s - 1) * (2 * pw + gIn + gOut);
        drawAgreementPanel(newAxes(f, [xa, y0, pw, ph], P), Uu.S(s), P, AXd(s));
        drawBAPanel(newAxes(f, [xa + pw + gIn, y0, pw, ph], P), Uu.S(s), P, AXd(s));
    end
    saveFigure2(f, fileBase, {'.pdf', '.eps'}, false);
end


function drawAgreementPanel(ax, Sk, P, A)
    D = Sk.auc;
    sc = 1e3; unitTxt = '10^{3} counts \cdot nm';
    x = D.pairRaw / sc; y = D.pairNmf / sc;
    lim = A.lim; st = A.st; hi = lim(2);
    plot(ax, lim, lim, '--', 'Color', P.supRefColor, 'LineWidth', P.lineWidth);
    drawDots(ax, x, y, D.animalRaw / sc, D.animalNmf / sc, Sk.color, numel(D.pairRaw) > numel(D.animal), P);
    set(ax, 'XLim', lim, 'YLim', lim, 'XTick', 0:st:hi, 'YTick', 0:st:hi);
    xlabel(ax, ['Ground-truth AUC (' unitTxt ')']);
    ylabel(ax, ['Unmixed AUC (' unitTxt ')']);
    C = Sk.ccc;
    txt = {sprintf('CCC = %.4f', C.ccc), sprintf('95%% CI %.4f to %.4f', C.lo, C.hi), nTextSup(D)};
    for j = 1:numel(txt)
        text(ax, 0.04, 0.97 - 0.085 * (j - 1), txt{j}, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontName', P.fontName, 'FontSize', P.fontSize);
    end
    supPanelTitle(ax, Sk, 'AUC', P, 1.04);
end


function drawBAPanel(ax, Sk, P, A)
    D = Sk.auc; B = Sk.ba; G = Sk.reg;
    sc = 1e3; unitTxt = '10^{3} counts \cdot nm';
    x = D.pairRaw / sc;
    xl = A.lim; hi = xl(2); lo = 0; yl = A.baY; st = A.baSt;
    plot(ax, xl, [0 0], '-', 'Color', P.supRefColor, 'LineWidth', P.axisLineWidth);
    plot(ax, xl, [B.bias B.bias], '-', 'Color', Sk.color, 'LineWidth', P.lineWidth);
    plot(ax, xl, [B.lo B.lo], ':', 'Color', Sk.color, 'LineWidth', P.lineWidth);
    plot(ax, xl, [B.hi B.hi], ':', 'Color', Sk.color, 'LineWidth', P.lineWidth);
    drawDots(ax, x, D.pairPct, D.animalRaw / sc, D.animalPct, Sk.color, numel(D.pairRaw) > numel(D.animal), P);
    drawReg = strcmpi(P.baRegLine, 'all') || (strcmpi(P.baRegLine, 'significant') && G.p < 0.05);
    if drawReg && isfinite(G.slope)
        xr = [G.xMin G.xMax];
        plot(ax, xr / sc, G.intercept + G.slope * xr, '-', 'Color', [0 0 0], 'LineWidth', 1.0);
    end
    sx = A.st;
    set(ax, 'XLim', xl, 'YLim', yl, 'XTick', 0:sx:hi, 'YTick', tickList(yl, st));
    xlabel(ax, ['Ground-truth AUC (' unitTxt ')']);
    ylabel(ax, 'Difference (%)');
    txt = {sprintf('bias %+.2f%%', B.bias), sprintf('95%% LoA %+.2f to %+.2f%%', B.lo, B.hi)};
    for j = 1:numel(txt)
        text(ax, 0.5, 1.12 - 0.085 * (j - 1), txt{j}, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontName', P.fontName, 'FontSize', P.fontSize, 'Clipping', 'off');
    end
    supPanelTitle(ax, Sk, 'AUC', P, 1.24);
    if ~strcmpi(P.baRegLine, 'none') && isfinite(G.slope)
        e = round(log10(sc));
        dd = max(0, 1 - floor(log10(abs(sc * G.slope) + eps)));
        fm = sprintf('%%.%df', dd);
        rt = {sprintf(['slope ' fm ' [' fm ', ' fm '] %%/10^{%d}'], sc * G.slope, sc * G.ciLo, sc * G.ciHi, e), ...
              sprintf('R^2 = %.2f, %s', G.r2, pTextSup(G.p))};
        for j = 1:2
            text(ax, 0.5, -0.27 - 0.09 * (j - 1), rt{j}, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top', 'FontName', P.fontName, 'FontSize', P.fontSize, 'Clipping', 'off');
        end
    end
end


function drawDots(ax, xp, yp, xa, ya, col, showPairs, P)
% mixtures as light (tinted, opaque) dots underneath; mouse means as dark dots on top
    if showPairs
        light = col + (1 - col) * P.pairTint;
        plot(ax, xp, yp, 'o', 'MarkerSize', P.dotSize, 'MarkerFaceColor', light, 'MarkerEdgeColor', 'none', ...
            'LineStyle', 'none');
    end
    plot(ax, xa, ya, 'o', 'MarkerSize', P.dotSize, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', 'LineStyle', 'none');
end


function supPanelTitle(ax, Sk, what, P, y)
    c = Sk.color;
    text(ax, 0.5, y, sprintf('\\color[rgb]{%.3f %.3f %.3f}%s\\color[rgb]{0 0 0} %s %s', c(1), c(2), c(3), ...
        tx(Sk.name), Sk.band, what), 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontName', P.fontName, 'FontSize', P.fontSize + 1, 'Clipping', 'off');
end


function t = nTextSup(D)
    na = numel(D.animal); np = numel(D.pairRaw);
    if np > na, t = sprintf('n = %d mice (%d mixtures)', na, np); else, t = sprintf('n = %d mice', na); end
end


function s = pTextSup(p)
    if ~isfinite(p), s = 'p = n/a';
    elseif p < 0.001, s = 'p < 0.001';
    else, s = sprintf('p = %.3f', p);
    end
end


function G = baRegression(x, y)
% least-squares line y = intercept + slope * x on the mouse-level values (proportional bias,
% Bland & Altman 1999); slope tested against 0 (t-test, df = n - 2), 95% CI of the slope
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    G = struct('n', numel(x), 'slope', NaN, 'intercept', NaN, 'se', NaN, 'ciLo', NaN, 'ciHi', NaN, ...
        't', NaN, 'df', numel(x) - 2, 'p', NaN, 'r2', NaN, 'xMin', NaN, 'xMax', NaN);
    if G.n < 3, return; end
    X = [ones(G.n, 1), x];
    cf = X \ y;
    res = y - X * cf;
    sxx = sum((x - mean(x)) .^ 2);
    G.intercept = cf(1); G.slope = cf(2);
    G.se = sqrt(sum(res .^ 2) / G.df / sxx);
    G.r2 = 1 - sum(res .^ 2) / max(sum((y - mean(y)) .^ 2), eps);
    tc = tinvL(0.975, G.df);
    G.ciLo = G.slope - tc * G.se; G.ciHi = G.slope + tc * G.se;
    if G.se > 0
        G.t = G.slope / G.se;
        G.p = 2 * (1 - tcdfL(abs(G.t), G.df));
    end
    G.xMin = min(x); G.xMax = max(x);
end


function B = blandAltman(d)
% bias and 95% limits of agreement (bias +/- 1.96 SD)
    d = d(isfinite(d));
    B.n = numel(d); B.bias = mean(d); B.sd = std(d);
    B.lo = B.bias - 1.96 * B.sd; B.hi = B.bias + 1.96 * B.sd;
end


function [sc, txt] = auSc(mx)
% axis scale for AUC values: 10^3 steps
    e = 3 * floor(log10(max(mx, 1)) / 3);
    if e < 3, e = 0; end
    sc = 10 ^ e;
    if e == 0, txt = 'counts \cdot nm'; else, txt = sprintf('10^{%d} counts \\cdot nm', e); end
end


function [yl, st] = baLimits(v, minSpan)
% automatic y range on tick values (1-2-5 x 10^k, 4-6 ticks); at least minSpan percentage points
    v = v(isfinite(v));
    if isempty(v), yl = [-1 1]; st = 0.5; return; end
    if nargin < 2 || isempty(minSpan), minSpan = 0; end
    r = max(max(v) - min(v), 0.1);
    lo0 = min(v) - 0.015 * r; hi0 = max(v) + 0.015 * r;
    if hi0 - lo0 < minSpan
        c = (min(v) + max(v)) / 2;
        lo0 = c - 0.485 * minSpan; hi0 = c + 0.485 * minSpan;
    else
        if min(v) >= 0, lo0 = 0; end
        if max(v) <= 0, hi0 = 0; end
    end
    r = hi0 - lo0;
    best = inf; yl = [lo0 hi0]; st = niceStep(r);
    for p = 10 .^ (floor(log10(r)) - 2:floor(log10(r)) + 1)
        for m = [1 2 2.5 5]
            s = m * p;
            lo = floor(lo0 / s + 1e-9) * s; hi = ceil(hi0 / s - 1e-9) * s;
            nt = numel(tickList([lo hi], s));
            if nt < 4 || nt > 6, continue; end
            score = (hi - lo) / r + 0.15 * (m == 2.5);
            if score < best, best = score; yl = [lo hi]; st = s; end
        end
    end
end


function t = tickList(lim, st)
    t = (ceil(lim(1) / st - 1e-9):floor(lim(2) / st + 1e-9)) * st;
    t(abs(t) < 1e-9 * st) = 0;
end


function t = limTicks(lim, st)
% ticks at multiples of st plus both axis limits (e.g. [-0.5 4], st 1: -0.5, 0, 1, 2, 3, 4)
    t = unique([lim(1), tickList(lim, st), lim(2)]);
end


function v = niceUpper(x)
    s = niceStep(x);
    v = ceil(x / s) * s;
end


function writeUnmixingStats(U, R, P, xls)
    if exist(xls, 'file') == 2, try, delete(xls); catch, end, end %#ok<NOCOM>
    f1 = sprintf('%s factor', R.names{1}); f2 = sprintf('%s factor', R.names{2});
    try
        C = { ...
            'Agreement of unmixed vs ground-truth AUC in every condition of the spike-recovery test', ''; ...
            'Analysis', sprintf('%s; scaling: %s; AUC windows: %s', R.version, R.scaling, R.windowRule); ...
            'Condition', '<green>x<k>_<red>x<k>: the green and red recordings x their brightness factor, then mixed and unmixed by NMF-RI'; ...
            'Percentage difference', '(unmixed - ground truth) / ground truth x 100'; ...
            'Mouse level', 'one value per mouse: mean over the mixtures that include it; all statistics use mouse-level values'; ...
            'Agreement', 'Lin''s CCC of mouse-level AUC, 95% CI by Fisher z (Lin 1989)'; ...
            'Bland-Altman', 'bias = mean mouse-level difference; 95% limits of agreement = bias +/- 1.96 SD'; ...
            'Proportional bias', 'regression of the mouse-level difference (%) on the ground-truth AUC; slope t-test (df = n - 2)'; ...
            'TOST', sprintf('two one-sided t-tests against +/-%g%% (alpha = %g each, %g%% CI)', P.margin, P.alpha, 100 * (1 - 2 * P.alpha)); ...
            'Date', datestr(now, 'yyyy-mm-dd HH:MM:SS')};
        writecell(clean(C), xls, 'Sheet', 'README');
        A = {'Condition', f1, f2, 'Sensor', 'Measure', 'n mice', 'n mixtures', 'CCC', 'CCC 95% CI lower', 'CCC 95% CI upper'};
        B = {'Condition', f1, f2, 'Sensor', 'Measure', 'n mice', 'Bias (%)', 'SD (%)', '95% LoA lower (%)', '95% LoA upper (%)'};
        G = {'Condition', f1, f2, 'Sensor', 'Measure', 'n mice', 'Slope (% per 10^3 counts nm)', 'Slope 95% CI lower', ...
             'Slope 95% CI upper', 'Intercept (%)', 'R^2', 't', 'df', 'p (slope = 0)'};
        T = {'Condition', f1, f2, 'Sensor', 'Measure', 'n mice', 'Mean difference (%)', 'SD (%)', ...
             sprintf('%g%% CI lower (%%)', 100 * (1 - 2 * P.alpha)), sprintf('%g%% CI upper (%%)', 100 * (1 - 2 * P.alpha)), ...
             'Margin (+/- %)', 'p TOST', 'Equivalent'};
        M = {'Condition', f1, f2, 'Sensor', 'Measure', 'Mouse', 'n mixtures', 'Ground truth (mean)', 'Unmixed (mean)', ...
             'SD of unmixed across mixtures', 'Difference (%)'};
        for u = 1:numel(U)
            c = U(u).cond; nm = U(u).name;
            for s = 1:2
                Sk = U(u).S(s); ms = [Sk.band ' AUC'];
                A(end+1, :) = {nm, c(1), c(2), Sk.name, ms, numel(Sk.auc.animal), numel(Sk.auc.pairRaw), ...
                    Sk.ccc.ccc, Sk.ccc.lo, Sk.ccc.hi}; %#ok<AGROW>
                B(end+1, :) = {nm, c(1), c(2), Sk.name, ms, Sk.ba.n, Sk.ba.bias, Sk.ba.sd, Sk.ba.lo, Sk.ba.hi}; %#ok<AGROW>
                g = Sk.reg;
                G(end+1, :) = {nm, c(1), c(2), Sk.name, ms, g.n, 1e3 * g.slope, 1e3 * g.ciLo, 1e3 * g.ciHi, ...
                    g.intercept, g.r2, g.t, g.df, g.p}; %#ok<AGROW>
                for which = 1:2
                    if which == 1, t = Sk.tostAUC; mt = ms; Dm = Sk.auc; else, t = Sk.tostDAUC; mt = [Sk.band ' delta AUC']; Dm = Sk.dauc; end
                    T(end+1, :) = {nm, c(1), c(2), Sk.name, mt, t.n, t.mean, t.sd, t.ciLo, t.ciHi, t.margin, t.p, ...
                        ternary(t.equivalent, 'yes', 'no')}; %#ok<AGROW>
                    for a = 1:numel(Dm.animal)
                        M(end+1, :) = {nm, c(1), c(2), Sk.name, mt, Dm.animal{a}, Dm.animalNPairs(a), Dm.animalRaw(a), ...
                            Dm.animalNmf(a), Dm.animalSdNmf(a), Dm.animalPct(a)}; %#ok<AGROW>
                    end
                end
            end
        end
        writecell(clean(A), xls, 'Sheet', 'Agreement');
        writecell(clean(B), xls, 'Sheet', 'BlandAltman');
        writecell(clean(G), xls, 'Sheet', 'ProportionalBias');
        writecell(clean(T), xls, 'Sheet', 'TOST');
        writecell(clean(M), xls, 'Sheet', 'MouseMeans');
        fprintf('   saved %s\n', xls);
    catch ME
        fprintf('   %s not saved: %s\n', xls, ME.message);
    end
end


%% =====================================================================
%  Helpers (unchanged from NMF_Hybrid_Validation_v2 / NMF_AUC_Sweep)
%  =====================================================================
function s = tx(s)
% Escape TeX special characters in names.
    s = strrep(strrep(strrep(s, '\', '\\'), '_', '\_'), '^', '\^');
end


function [yl, st] = niceAxis(v)
% Limits on ticks of 1-2-5 x 10^k, 4-7 ticks, covering v.
    v = v(isfinite(v));
    lo = min(v); hi = max(v);
    if hi <= lo, hi = lo + 1; end
    best = [];
    for e = -2:4
        for m = [1 2 5]
            st = m * 10 ^ e;
            a = floor(lo / st - 1e-9) * st; b = ceil(hi / st + 1e-9) * st;
            nt = round((b - a) / st) + 1;
            if nt >= 4 && nt <= 7 && (isempty(best) || b - a < best(2) - best(1))
                best = [a b st];
            end
        end
    end
    if isempty(best), best = [lo hi (hi - lo) / 4]; end
    yl = best(1:2); st = best(3);
end


function s = niceStep(range)
% 1-2-2.5-5 x 10^k step giving about 4 intervals.
    raw = range / 4;
    if ~(raw > 0), s = 1; return; end
    p = 10 ^ floor(log10(raw));
    m = [1 2 2.5 5 10] * p;
    s = m(find(m >= raw, 1));
end


function A = emptySource()
    A = struct('name', '', 'file', '', 'type', 'sensor', 'wl', [], 'Y', [], 'spec', [], ...
        'nTrials', 0, 'nSamples', 0, 'nRows', 0, 'fileMin', NaN, 'dec', 0);
end


function [L, w] = listAnimals(folder, tag)
% Mouse subfolders containing a file with the tag; if there are none,
% TXT files with the tag directly in the folder.
    w = {};
    L = repmat(emptySource(), 0, 1);
    d = dir(folder);
    d = d([d.isdir]);
    subs = {d.name};
    subs = sort(subs(~cellfun(@(x) strncmp(x, '.', 1), subs)));
    for k = 1:numel(subs)
        [f, ww] = findFile(fullfile(folder, subs{k}), tag, subs{k});
        w = [w, ww]; %#ok<AGROW>
        if isempty(f), continue; end
        A = emptySource();
        A.name = subs{k};
        A.file = f;
        L(end+1, 1) = A; %#ok<AGROW>
    end
    if ~isempty(L)
        return;
    end
    t = dir(fullfile(folder, '*.txt'));
    fn = sort({t.name});
    for k = 1:numel(fn)
        s = tagStart(fn{k}, tag);
        if isempty(s), continue; end
        nm = regexprep(fn{k}(1:s - 1), '[_.\- ]+$', '');
        if isempty(nm)
            [~, nm] = fileparts(fn{k});
        end
        A = emptySource();
        A.name = nm;
        A.file = fullfile(folder, fn{k});
        L(end+1, 1) = A; %#ok<AGROW>
    end
    if isempty(L)
        w{end+1} = sprintf('No "%s" files found in %s', tag, folder);
    end
end


function [A, w] = loadAnimal(A, opt, label)
% Reads the recording, subtracts the file offset (min of the time-averaged spectrum) and
% keeps the complete trials (single precision to save memory).
    w = {};
    t0 = tic;
    [wl, Y, nDropped, fmt] = readDyeTxt(A.file, opt.wlRange);
    [~, fn, fe] = fileparts(A.file);
    if strcmp(fmt, 'column')
        w{end+1} = sprintf('%s (%s): %s is a single spectrum (wavelength / value columns), not a time series; animal skipped.', ...
            A.name, label, [fn fe]);
        fprintf(' %-6s %-14s %-28s SKIPPED: single spectrum, not a time series\n', label, A.name, [fn fe]);
        A.nTrials = 0;
        A.Y = [];
        return;
    end
    if isempty(wl) || isempty(Y)
        w{end+1} = sprintf('%s (%s): no values within %g-%g nm in %s; animal skipped.', ...
            A.name, label, opt.wlRange(1), opt.wlRange(2), [fn fe]);
        fprintf(' %-6s %-14s %-28s SKIPPED: no values within %g-%g nm\n', label, A.name, [fn fe], opt.wlRange(1), opt.wlRange(2));
        A.nTrials = 0;
        A.Y = [];
        return;
    end
    if nDropped > 0
        w{end+1} = sprintf('%s (%s): %d incomplete rows removed.', A.name, label, nDropped);
    end
    L = opt.trialLength;
    nT = size(Y, 2);
    nTr = floor(nT / L);
    A.wl = wl;
    A.nSamples = nT;
    A.nRows = nT;
    A.nTrials = nTr;
    if nTr < 1
        w{end+1} = sprintf('%s (%s): fewer time points (%d) than one trial.', A.name, label, nT);
        A.Y = [];
        return;
    end
    if mod(nT, L) ~= 0
        w{end+1} = sprintf('%s (%s): %d time points is not a multiple of %d; %d complete trials used.', ...
            A.name, label, nT, L, nTr);
    end
    % Offset removal, one value per file: the minimum of the time-averaged
    % spectrum of the whole recording (450-650 nm), i.e. the baseline offset
    % of that file. The offset-free data are the ground truth (starts at 0,
    % same reference as the unmixed components); the time course is unchanged.
    A.fileMin = min(mean(Y, 2));
    A.Y = single(Y(:, 1:nTr * L) - A.fileMin);
    fprintf(' %-6s %-14s %-28s %d trials, offset %.1f counts (min of the time-averaged spectrum) (%.1f s)\n', ...
        label, A.name, [fn fe], nTr, A.fileMin, toc(t0));
end


function [M, interp] = sourceMatrix(A, wl, nS)
% Sensor: first nS samples (double). Dye: static spectrum repeated nS times.
    if isempty(A.wl) || isempty(wl)
        error('%s: empty wavelength axis (no values within the wavelength range).', A.name);
    end
    interp = numel(A.wl) ~= numel(wl);
    if ~interp
        interp = max(abs(A.wl(:) - wl(:))) > 1e-6;
    end
    if strcmp(A.type, 'dye')
        sp = A.spec(:);
        if interp
            sp = interp1(A.wl(:), sp, wl, 'linear', 0);
        end
        M = repmat(sp, 1, nS);
    else
        M = double(A.Y(:, 1:nS));
        if interp
            M = interp1(A.wl(:), M, wl, 'linear', 0);
        end
    end
end


function [ANmf, nIter] = callNmf(YP, A0, H0, cfg)
% NMF_RI (unchanged); its console output is captured to keep the log short.
    ANmf = []; nIter = 0; Hunused = []; %#ok<NASGU>
    t1 = cfg.theta1; t2 = cfg.theta2; aA = cfg.alphaA; aX = cfg.alphaX; %#ok<NASGU>
    evalc('[ANmf, Hunused, nIter] = NMF_RI(YP, A0, H0, t1, t2, aA, aX);');
end


function M = trialMeans(X, samples, nTr, L)
    M = zeros(size(X, 1), nTr);
    for k = 1:nTr
        M(:, k) = mean(X(:, (k - 1) * L + samples), 2);
    end
end


function C = aucPeakCompare(wl, T, N, TD, ND, types, names, opt, TT)
% Raw (ground truth) vs NMF, trial means (counts). T, N: reward spectra;
% TD, ND: delta spectra (reward - baseline).
% v2: each spectrum is measured around its own peak (AUC and delta AUC):
%   raw window = raw peak +/- aucHalfWidth, unmixed window = unmixed peak +/- aucHalfWidth
%   T1 centre: green peak in bandT1; T2 centre: red peak in bandT2
%   (green T2 also uses the red T2 centre of the same side).
% Peaks: found after a short moving average.
    hw = opt.aucHalfWidth;
    sm = @(y) smoothNm(wl, y, opt.peakSmoothNm);
    C.winT1raw = NaN;  C.winT1nmf = NaN;
    if strcmp(types{1}, 'sensor')
        C.winT1raw = bandPeak(wl, sm(T{1}), opt.bandT1);
        C.winT1nmf = bandPeak(wl, sm(N{1}), opt.bandT1);
    end
    C.winT2raw = bandPeak(wl, sm(T{2}), opt.bandT2);
    C.winT2nmf = bandPeak(wl, sm(N{2}), opt.bandT2);
    wr1 = C.winT1raw; wn1 = C.winT1nmf;       % v2: separate windows
    wr2 = C.winT2raw; wn2 = C.winT2nmf;

    m = repmat(mkMetric('', '', '', 1, NaN, NaN, 1, 1), 0, 1);
    gS = strcmp(types{1}, 'sensor');
    rS = strcmp(types{2}, 'sensor');
    if gS, cR = 4; else, cR = 1; end          % column of the red sensor panels
    if gS
        r1 = winAuc(wl, T{1}, wr1, hw);  n1 = winAuc(wl, N{1}, wn1, hw);
        r2 = winAuc(wl, T{1}, wr2, hw);  n2 = winAuc(wl, N{1}, wn2, hw);
        m(end+1, 1) = mkMetric('green_T1_AUC', [names{1} ' T1 AUC'], 'auc', 1, r1, n1, 1, 1);
        m(end+1, 1) = mkMetric('green_T2_AUC', [names{1} ' T2 AUC'], 'auc', 1, r2, n2, 1, 2);
        m(end+1, 1) = mkMetric('green_T2_T1_ratio', [names{1} ' T2/T1'], 'ratio', 1, r2 / r1, n2 / n1, 1, 3);
    end
    if rS
        m(end+1, 1) = mkMetric('red_T2_AUC', [names{2} ' T2 AUC'], 'auc', 2, ...
            winAuc(wl, T{2}, wr2, hw), winAuc(wl, N{2}, wn2, hw), 1, cR);
    end
    % delta AUC: delta spectrum (reward - baseline), same windows as the AUC above
    if gS
        m(end+1, 1) = mkMetric('green_T1_delta_AUC', [names{1} ' T1 delta AUC'], 'dauc', 1, ...
            winAuc(wl, TD{1}, wr1, hw), winAuc(wl, ND{1}, wn1, hw), 2, 1);
        m(end+1, 1) = mkMetric('green_T2_delta_AUC', [names{1} ' T2 delta AUC'], 'dauc', 1, ...
            winAuc(wl, TD{1}, wr2, hw), winAuc(wl, ND{1}, wn2, hw), 2, 2);
    end
    if rS
        m(end+1, 1) = mkMetric('red_T2_delta_AUC', [names{2} ' T2 delta AUC'], 'dauc', 2, ...
            winAuc(wl, TD{2}, wr2, hw), winAuc(wl, ND{2}, wn2, hw), 2, cR);
    end
    if gS
        m(end+1, 1) = mkMetric('green_T1_peak', [names{1} ' T1 peak'], 'peak', 1, ...
            bandPeak(wl, sm(T{1}), opt.bandT1), bandPeak(wl, sm(N{1}), opt.bandT1), 3, 1);
    end
    if rS
        m(end+1, 1) = mkMetric('red_T2_peak', [names{2} ' T2 peak'], 'peak', 2, ...
            bandPeak(wl, sm(T{2}), opt.bandT2), bandPeak(wl, sm(N{2}), opt.bandT2), 3, cR);
    end
    % trial-level raw AUC (reward spectrum of each trial, same raw window),
    % used for the trial-to-trial coefficient of variation
    if nargin >= 9 && ~isempty(TT)
        for q = 1:numel(m)
            switch m(q).key
                case 'green_T1_AUC', Y = TT{1}; c = wr1;
                case 'green_T2_AUC', Y = TT{1}; c = wr2;
                case 'red_T2_AUC',   Y = TT{2}; c = wr2;
                otherwise, continue;
            end
            v = nan(size(Y, 2), 1);
            for k = 1:size(Y, 2)
                v(k) = winAuc(wl, Y(:, k), c, hw);
            end
            m(q).rawTrials = v;
        end
    end
    C.m = m;
end


function s = mkMetric(key, label, kind, slot, raw, nmf, row, col)
% row 1 AUC / ratio, row 2 delta AUC, row 3 peak; col = figure column
    s = struct('key', key, 'label', label, 'kind', kind, 'slot', slot, 'raw', raw, 'nmf', nmf, ...
        'row', row, 'col', col, 'rawTrials', []);
end


function y = smoothNm(wl, y, nm)
% Moving average over about nm nanometres (edges: average of the available points).
    y = y(:);
    d = median(diff(wl(:)));
    k = max(1, round(nm / d));
    if k <= 1
        return;
    end
    w = ones(k, 1);
    y = conv(y, w, 'same') ./ conv(ones(size(y)), w, 'same');
end


function p = bandPeak(wl, y, band)
% Wavelength of the highest local maximum within the band (band maximum if none).
    wl = wl(:); y = y(:);
    n = numel(y);
    p = NaN;
    lm = false(n, 1);
    lm(2:n-1) = y(2:n-1) > y(1:n-2) & y(2:n-1) >= y(3:n);
    inBand = wl >= band(1) & wl <= band(2);
    cand = find(inBand & lm);
    if isempty(cand)
        cand = find(inBand);
    end
    if isempty(cand)
        return;
    end
    [~, k] = max(y(cand));
    p = wl(cand(k));
end


function a = winAuc(wl, y, c, hw)
    a = NaN;
    if ~isfinite(c), return; end
    m = wl(:) >= c - hw & wl(:) <= c + hw;
    if sum(m) < 2, return; end
    y = y(:);
    a = trapz(wl(m), y(m));
end


function progress(hw, frac, msg)
    if isempty(hw) || ~ishandle(hw), return; end
    try
        waitbar(min(max(frac, 0), 1), hw, strrep(msg, '_', '\_'));
        drawnow;
    catch
    end
end


function closeBar(hw)
    if ~isempty(hw) && ishandle(hw)
        delete(hw);
    end
end


function c = sensorColor(name, opt, fallback)
    c = fallback;
    for k = 1:size(opt.sensorColors, 1)
        if ~isempty(strfind(lower(name), lower(opt.sensorColors{k, 1}))) %#ok<STREMP>
            c = opt.sensorColors{k, 2};
            return;
        end
    end
end


function em = loadEmissionSpectrum(csvFile, wlTarget)
% As NMF_RI_Photometry: emission column, min-subtracted, max-normalised,
% interpolated onto the target wavelengths (0 outside), negative -> 0.
    T = readtable(csvFile, 'VariableNamingRule', 'preserve');
    if width(T) < 2
        error('Spectrum CSV needs at least two columns: %s', csvFile);
    end
    wlCol = convertColumnToDouble(T{:, 1});
    colNames = string(T.Properties.VariableNames);
    emIdx = [];
    for i = 2:numel(colNames)
        name = lower(colNames(i));
        if contains(name, 'em') && ~contains(name, 'ex') && ~contains(name, '2p')
            emIdx = i;
            break;
        end
    end
    if isempty(emIdx)
        error('No emission column found in: %s', csvFile);
    end
    raw = convertColumnToDouble(T{:, emIdx});
    valid = isfinite(raw) & isfinite(wlCol);
    wlSrc = wlCol(valid);
    emSrc = raw(valid);
    if numel(wlSrc) < 2
        error('Not enough valid emission data in: %s', csvFile);
    end
    [wlSrc, sortIdx] = sort(wlSrc);
    emSrc = emSrc(sortIdx);
    [wlSrc, uniqueIdx] = unique(wlSrc, 'stable');
    emSrc = emSrc(uniqueIdx);
    emSrc = emSrc - min(emSrc);
    if max(emSrc) > 0
        emSrc = emSrc / max(emSrc);
    end
    em = interp1(wlSrc, emSrc, wlTarget, 'linear', 0);
    em = max(em, 0);
end


function x = convertColumnToDouble(raw)
    if isnumeric(raw)
        x = double(raw);
    else
        x = str2double(string(raw));
    end
    x = x(:);
end


function s = tagStart(fileName, tag)
% Start index of the tag in the file name: tag preceded by a separator
% (or, if none, at the start of the name) and not followed by a letter,
% so 'G' does not match 'GR'. Case-insensitive.
    e = regexptranslate('escape', tag);
    s = regexpi(fileName, ['(?<=[_.\- ])' e '(?![A-Za-z])'], 'once');
    if isempty(s)
        s = regexpi(fileName, ['^' e '(?![A-Za-z])'], 'once');
    end
end


function [f, w] = findFile(folder, tag, an)
    f = '';
    w = {};
    if isempty(folder) || exist(folder, 'dir') ~= 7
        return;
    end
    d = dir(fullfile(folder, '*.txt'));
    names = sort({d.name});
    hit = ~cellfun(@isempty, cellfun(@(n) tagStart(n, tag), names, 'UniformOutput', false));
    names = names(hit);
    if isempty(names)
        w{end+1} = sprintf('%s: no file with "%s" (*.txt) in %s', an, tag, folder);
        return;
    end
    if numel(names) > 1
        w{end+1} = sprintf('%s: %d files match "%s"; using %s', an, numel(names), tag, names{1});
    end
    f = fullfile(folder, names{1});
end


function [wl, Y, nDropped, fmt] = readDyeTxt(file, wlRange)
% Dye spectrum TXT in either format:
%   'column'  one spectrum, one "wavelength <tab> value" pair per line
%             (OceanView single-spectrum export)
%   'row'     same layout as the sensor recordings (wavelength line, then
%             one row per spectrum); read with readSpectralTxt
    marker = '>>>>>Begin Spectral Data<<<<<';
    % look at the first 50 lines after the marker only (sensor files are large)
    fid = fopen(file, 'r');
    if fid < 0
        error('Cannot open %s', file);
    end
    found = false;
    for k = 1:500
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        if ~isempty(strfind(ln, marker)) %#ok<STREMP>
            found = true;
            break;
        end
    end
    if ~found, frewind(fid); end
    two = 0; more = 0;
    for i = 1:50
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        s = strtrim(ln);
        if isempty(s) || s(1) == '#', continue; end
        if numel(sscanf(s, '%f')) == 2, two = two + 1; else, more = more + 1; end
    end
    fclose(fid);
    if two < 10 || more > 0
        fmt = 'row';
        [wl, Y, nDropped] = readSpectralTxt(file, wlRange);
        return;
    end
    txt = fileread(file);
    k = strfind(txt, marker);
    if ~isempty(k)
        body = txt(k(1) + numel(marker):end);
    else
        body = txt;
    end
    lines = regexp(body, '\r?\n', 'split');
    vals = nan(numel(lines), 2);
    n = 0;
    other = 0;
    for i = 1:numel(lines)
        s = strtrim(lines{i});
        if isempty(s) || s(1) == '#'
            continue;
        end
        v = sscanf(s, '%f');
        if numel(v) == 2
            n = n + 1;
            vals(n, :) = v(:)';
        else
            other = other + 1;
        end
    end
    vals = vals(1:n, :);
    if n >= 10 && other == 0 && all(diff(vals(:, 1)) > 0)
        fmt = 'column';
        keep = vals(:, 1) >= wlRange(1) & vals(:, 1) <= wlRange(2) & isfinite(vals(:, 2));
        wl = vals(keep, 1);
        Y = vals(keep, 2);
        nDropped = 0;
    else
        fmt = 'row';
        [wl, Y, nDropped] = readSpectralTxt(file, wlRange);
    end
end


function [wl, Y, nDropped] = readSpectralTxt(file, wlRange)
% Reads an OceanView-style spectral TXT (with or without the
% ">>>>>Begin Spectral Data<<<<<" marker). Lines starting with '#' after the
% marker are skipped. Each data row: text fields followed by one value per
% wavelength; the last nPix values of a row form the spectrum.
    fid = fopen(file, 'r');
    if fid < 0
        error('Cannot open %s', file);
    end
    cleaner = onCleanup(@() fclose(fid));

    marker = '>>>>>Begin Spectral Data<<<<<';
    found = false;
    for k = 1:500
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        if ~isempty(strfind(ln, marker)) %#ok<STREMP>
            found = true;
            break;
        end
    end
    if ~found
        frewind(fid);
    end

    wlLine = '';
    while true
        ln = fgetl(fid);
        if ~ischar(ln)
            error('No wavelength line found in %s', file);
        end
        s = strtrim(ln);
        if isempty(s) || s(1) == '#'
            continue;
        end
        wlLine = ln;
        break;
    end
    tok = strsplit(wlLine, sprintf('\t'));
    v = str2double(tok);
    wlAll = v(isfinite(v));
    wlAll = wlAll(:);
    nPix = numel(wlAll);
    if nPix < 2
        error('Could not read wavelengths in %s', file);
    end

    pos = ftell(fid);
    first = fgetl(fid);
    while ischar(first) && isempty(strtrim(first))
        pos = ftell(fid);
        first = fgetl(fid);
    end
    if ~ischar(first)
        error('No data rows in %s', file);
    end
    ftok = strsplit(first, sprintf('\t'), 'CollapseDelimiters', false);
    while ~isempty(ftok) && isempty(strtrim(ftok{end}))
        ftok(end) = [];
    end
    nPrefix = numel(ftok) - nPix;
    if nPrefix < 0
        error('Data rows have fewer values than wavelengths in %s', file);
    end
    fseek(fid, pos, 'bof');

    keep = wlAll >= wlRange(1) & wlAll <= wlRange(2);
    nKeep = sum(keep);
    block = 5000;
    Y = zeros(nKeep, block);
    nRow = 0;
    nDropped = 0;
    while true
        ln = fgetl(fid);
        if ~ischar(ln)
            break;
        end
        if isempty(strtrim(ln))
            continue;
        end
        if nPrefix > 0
            tabs = find(ln == sprintf('\t'), nPrefix);
            if numel(tabs) < nPrefix
                nDropped = nDropped + 1;
                continue;
            end
            ln = ln(tabs(end) + 1:end);
        end
        v = sscanf(ln, '%f');
        if numel(v) ~= nPix || any(~isfinite(v(keep)))
            nDropped = nDropped + 1;
            continue;
        end
        nRow = nRow + 1;
        if nRow > size(Y, 2)
            Y = [Y, zeros(nKeep, block)]; %#ok<AGROW>
        end
        Y(:, nRow) = v(keep);
    end
    Y = Y(:, 1:nRow);
    wl = wlAll(keep);
end


function p = stripSep(p)
    p = strtrim(p);
    while numel(p) > 1 && (p(end) == '\' || p(end) == '/')
        p(end) = [];
    end
end


function s = safeFileName(s)
    s = regexprep(s, '[<>:"/\\|?*]', '_');
end


function T = tost(d, margin, alpha)
% Two one-sided tests against +/- margin; (1 - 2 alpha) CI of the mean.
    d = d(isfinite(d));
    n = numel(d);
    T.n = n; T.margin = margin;
    T.mean = mean(d); T.sd = std(d);
    T.se = T.sd / sqrt(n); T.df = n - 1;
    T.ciLevel = 100 * (1 - 2 * alpha);
    T.tLower = NaN; T.pLower = NaN; T.tUpper = NaN; T.pUpper = NaN;
    T.p = NaN; T.ciLo = NaN; T.ciHi = NaN; T.equivalent = false;
    if n < 2 || T.se == 0, return; end
    tc = tinvL(1 - alpha, T.df);
    T.ciLo = T.mean - tc * T.se;
    T.ciHi = T.mean + tc * T.se;
    T.tLower = (T.mean + margin) / T.se;          % H0: mean <= -margin
    T.pLower = 1 - tcdfL(T.tLower, T.df);
    T.tUpper = (T.mean - margin) / T.se;          % H0: mean >= +margin
    T.pUpper = tcdfL(T.tUpper, T.df);
    T.p = max(T.pLower, T.pUpper);
    T.equivalent = T.p < alpha;
end


function p = tcdfL(t, df)
% Student t CDF without the Statistics Toolbox.
    x = df / (df + t ^ 2);
    tail = 0.5 * betainc(x, df / 2, 0.5);
    if t > 0, p = 1 - tail; else, p = tail; end
end


function t = tinvL(p, df)
    t = fzero(@(x) tcdfL(x, df) - p, [-1e3 1e3]);
end


function f = newFigure(w, h)
    f = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 w h], ...
        'Renderer', 'painters', 'InvertHardcopy', 'off');
    set(f, 'UserData', [w h], 'PaperUnits', 'centimeters', 'PaperSize', [w h], ...
        'PaperPosition', [0 0 w h]);
end


function ax = newAxes(f, posCm, P)
    ax = axes(f, 'Units', 'centimeters', 'Position', posCm);
    hold(ax, 'on');
    set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.02 0.02], 'LineWidth', P.axisLineWidth, ...
        'FontName', P.fontName, 'FontSize', P.fontSize, 'XColor', 'k', 'YColor', 'k', 'Layer', 'top', ...
        'Color', 'none');
    set(get(ax, 'XLabel'), 'FontSize', P.fontSize, 'FontName', P.fontName);
    set(get(ax, 'YLabel'), 'FontSize', P.fontSize, 'FontName', P.fontName);
    try
        set(ax, 'XTickLabelRotation', 0, 'YTickLabelRotation', 0);   % keep tick labels horizontal
    catch
    end
end


function saveFigure(f, fileBase)
    for ext = {'.pdf', '.eps'}
        file = [fileBase ext{1}];
        if exist(file, 'file') == 2
            try, delete(file); catch, end %#ok<NOCOM>
        end
        try
            if exist('exportgraphics') > 0 %#ok<EXIST>
                exportgraphics(f, file, 'ContentType', 'vector', 'BackgroundColor', 'white');
            elseif strcmp(ext{1}, '.pdf')
                print(f, file, '-dpdf', '-painters');
            else
                print(f, file, '-depsc', '-painters');
            end
            fprintf('   saved %s\n', file);
        catch ME
            fprintf('   %s not saved: %s\n', file, ME.message);
        end
    end
    close(f);
end


function C = clean(C)
    C(cellfun(@isempty, C)) = {''};
    for k = 1:numel(C)
        if islogical(C{k}), C{k} = double(C{k}); end
        if isnumeric(C{k}) && isscalar(C{k}) && ~isfinite(C{k}), C{k} = 'NaN'; end
    end
end


function C = lincc(x, y)
% Lin's concordance correlation coefficient with a 95% CI (Fisher z, Lin 1989).
    x = x(:); y = y(:);
    n = numel(x);
    C.n = n;
    C.ccc = NaN; C.lo = NaN; C.hi = NaN;
    if n < 3, return; end
    mx = mean(x); my = mean(y);
    sx = mean((x - mx) .^ 2); sy = mean((y - my) .^ 2);
    sxy = mean((x - mx) .* (y - my));
    rc = 2 * sxy / (sx + sy + (mx - my) ^ 2);
    r = sxy / sqrt(sx * sy);
    u = (mx - my) / (sx * sy) ^ 0.25;
    C.ccc = rc;
    C.r = r;
    if abs(rc) >= 1 || n < 4, return; end
    se = sqrt(((1 - r ^ 2) * rc ^ 2 / ((1 - rc ^ 2) * r ^ 2) ...
        + 2 * rc ^ 3 * (1 - rc) * u ^ 2 / (r * (1 - rc ^ 2) ^ 2) ...
        - rc ^ 4 * u ^ 4 / (2 * r ^ 2 * (1 - rc ^ 2) ^ 2)) / (n - 2));
    z = atanh(rc);
    C.lo = tanh(z - 1.96 * se);
    C.hi = tanh(z + 1.96 * se);
end
