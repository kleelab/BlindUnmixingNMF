function varargout = NMF_Hybrid_Validation_v2(cfg)
% NMF_Hybrid_Validation_v2
%
% *** VERSION 2: SEPARATE AUC WINDOWS FOR RAW AND UNMIXED ***
%   Only difference from NMF_Hybrid_Validation (v1): in the raw vs unmixed
%   AUC comparison (figure *_AUC_peak_v2.pdf, sheet AUC_peak) each spectrum
%   is measured around ITS OWN peak:
%     raw AUC      = raw spectrum,     window = raw peak     +/- 10 nm
%     unmixed AUC  = unmixed spectrum, window = unmixed peak +/- 10 nm
%   T1 window: green peak within 500-540 nm (raw green / unmixed green).
%   T2 window: red peak within 580-620 nm (raw red / unmixed red); the green
%   T2 AUC uses the same T2 centres (raw red peak for raw, unmixed red peak
%   for unmixed), because the green sensor has no own peak there.
%   v1 used one window (raw peak +/- 10 nm) for both raw and unmixed.
%   Everything else (mixtures, NMF-RI, R^2, figures) is identical to v1.
%   All output names end in _v2, so v1 results are not overwritten.
%
% Validation of blind NMF-RI unmixing with in silico mixtures of measured single-sensor data.
%
% Two inputs: a green slot and a red slot. Each slot is either
%   a FOLDER  = sensor time series (one subfolder per mouse; green slot: G/g
%               file, 470 nm excitation; red slot: R/r file, 560 nm excitation)
%   a FILE    = control dye (in vivo spectrum TXT, e.g. tdTomato, eGFP);
%               all rows are averaged into one static spectrum
%
% Modes
%   folder + folder  sensor + sensor   pairs of a green-only and a red-only
%                    animal (1:1 or all combinations); 6 panels
%   folder + file    sensor + dye      one mixture per sensor animal; 4 panels
%   file + folder    dye + sensor      (e.g. eGFP + jRGECO)
%
% Per mixture
%   1. Ground truth: measured single-sensor data minus the offset of each file
%      (one value per file: minimum of the time-averaged spectrum, 450-650 nm),
%      so the ground truth starts at 0, like the unmixed components.
%      Dye: mean of all rows of the dye file minus its minimum across wavelengths.
%   2. In silico mixture = green truth + red truth, trial k + trial k, first
%      min(n trials) trials. A dye spectrum is added unchanged to every time point.
%   3. Blind NMF-RI exactly as NMF_RI_Photometry: A0 from the two CSV spectra,
%      dataPreprocessing, NMF_RI (unchanged), H = max(1e6*eps, pinv(A'A) A' Y).
%      Component 1 = green slot, component 2 = red slot. No unmixed files saved.
%   4. Reward spectrum = mean over the reward samples of each trial;
%      delta spectrum = reward - baseline mean.
%      R^2 = 1 - SSres/SStot per trial over 450-650 nm, in counts
%      (truth = observed, NMF = predicted); mean +/- SEM over trials.
%
% Panels (row 1 reward, row 2 delta)
%   sensor + sensor: reward sum | green | red;  delta sum | green | red
%   sensor + dye:    reward sum | green | red;  delta: sensor only, below the
%                    sensor's reward panel (dye delta is 0 by construction).
%                    Dye delta leakage is written to Excel only:
%                    100 * AUC(NMF dye delta) / AUC(sensor truth delta),
%                    trial-mean spectra, signed AUC over 450-650 nm.
%
% Figures
%   Black line = ground truth; filled area = NMF component(s); sum panel:
%   component 1 from 0, component 2 stacked on top.
%   Mixture figures: reward row normalised (truth and NMF divided by the truth
%   peak; same y axis in every figure); delta row in counts (y axis per panel).
%   Summary figure: all panels normalised (delta: divided by max |truth|),
%   mean across mixtures; grey band = SEM across mixtures;
%   R^2 = mean +/- SEM of the per-mixture means.
%
% Usage
%   NMF_Hybrid_Validation_v2          % opens the input window
%   NMF_Hybrid_Validation_v2(cfg)     % runs without the window
%
%   cfg.greenPath   folder (green sensor) or file (green dye)
%   cfg.redPath     folder (red sensor) or file (red dye)
%   cfg.greenName, cfg.redName   optional; taken from the folder / file name
%   cfg.comp1Csv    emission spectrum CSV, green slot (initial A0, component 1)
%   cfg.comp2Csv    emission spectrum CSV, red slot   (initial A0, component 2)
%   cfg.outDir      output folder
%   cfg.wlRange     [450 650] (default)
%   cfg.theta1, cfg.theta2, cfg.alphaA, cfg.alphaX  (defaults 0.01 0.01 0.001 0.001)
%   cfg.pairing     sensor + sensor only: 'one-to-one' (default; animals sorted by
%                   name, each used once) or 'all' (every green with every red)
%
% Required on the MATLAB path: NMF_RI.m, dataPreprocessing.m
%
% Output (in cfg.outDir)
%   mixtures_v2\Mixture01_GCaMP_M1+jRGECO_R1_mixture_v2.pdf       (sensor + sensor)
%   mixtures_v2\Mixture01_GCaMP_M1+tdTomato_dye_mixture_v2.pdf  (sensor + dye)
%   <green>+<red>_mixture_summary_v2.pdf
%   <green>+<red>_mixture_AUC_peak_v2.pdf  raw vs NMF, one line per mixture:
%       reward spectra (trial means); v2: AUC window = own peak +/- 10 nm
%       (raw: raw peak, unmixed: unmixed peak; T1: green peak 500-540 nm,
%       T2: red peak 580-620 nm); peaks found separately (2 nm smoothing)
%       sensor + sensor : green T1 AUC, green T2 AUC, green T2/T1, red T2 AUC;
%                         delta AUC below each AUC (green T1, green T2, red T2);
%                         green T1 peak, red T2 peak
%       green sensor + red dye : green T1 / T2 AUC, T2/T1; green T1 / T2 delta AUC;
%                                green T1 peak
%       green dye + red sensor : red T2 AUC; red T2 delta AUC; red T2 peak
%       delta AUC = AUC of the delta spectrum (reward - baseline) in the AUC window
%   <green>+<red>_mixture_R2*.mat     analysis results for NMF_Hybrid_Equivalence (figures without re-analysis)
%   <green>+<red>_mixture_R2_v2.xlsx  sheets README (definitions, settings), R2,
%                                  R2_trials, AUC_peak, Summary, Info
%
% Note: dataPreprocessing subtracts the minimum of each time point before
%   NMF, as in the real-data pipeline. This offset is not part of the NMF
%   components; its mean is in the Excel file (Prep_offset_counts).

    if nargin < 1 || isempty(cfg)
        cfg = selectInputs();
        if isempty(cfg)
            fprintf('Cancelled.\n');
            if nargout > 0, varargout{1} = []; end
            return;
        end
    end
    S = runHybrid(cfg);
    if nargout > 0
        varargout{1} = S;
    end
end


%% =====================================================================
%  Settings
%  =====================================================================
function s = defaultSettings()
    s.wlRange         = [450 650];   % nm
    s.trialLength     = 1177;        % samples per trial (100 Hz)
    s.baselineSamples = 100:150;     % samples within each trial (1-based)
    s.rewardSamples   = 400:499;     % reward onset at sample 400, 1 s
    s.greenTag        = 'G';         % sensor file tag, green slot (470 nm excitation)
    s.redTag          = 'R';         % sensor file tag, red slot (560 nm excitation)

    % NMF-RI defaults (NMF_RI_Photometry, tab 2)
    s.theta1 = 0.01;
    s.theta2 = 0.01;
    s.alphaA = 0.001;
    s.alphaX = 0.001;

    green  = [0.13 0.55 0.25];
    red    = [0.80 0.15 0.15];
    gray   = [0.50 0.50 0.50];
    s.sensorColors = {'GCaMP', green; 'jRGECO', red; 'rdLight', red; ...
                      'tdTomato', gray; 'eGFP', gray};
    s.colorComp1 = green;     % fallback, green slot
    s.colorComp2 = red;       % fallback, red slot
    s.semAlpha   = 0.25;

    % normalised y axes
    % AUC / peak comparison (raw vs NMF, reward spectra, trial means)
    s.bandT1       = [500 540];   % nm, T1 peak search
    s.bandT2       = [580 620];   % nm, T2 peak search
    s.aucHalfWidth = 10;          % nm, v2: AUC window = own peak +/- aucHalfWidth
    s.peakSmoothNm = 2;           % nm, moving average before peak search

    s.yNormLim  = [-0.05 1];   % fixed y range of all normalised panels
    s.yTick     = 0.2;         % normalised y tick

    s.xTickSpectrum = 50;          % nm
    s.panelCm       = [4.2 3.6];   % panel width, height
    s.fontName      = 'Arial';
    s.fontSize      = 8;
end


%% =====================================================================
%  Input window
%  =====================================================================
function cfg = selectInputs()
    cfg = [];
    opt = defaultSettings();
    defRaw = 'D:\NMF data\Raw_data';
    defOut = 'D:\NMF data\analysis';

    fig = uifigure('Name', 'NMF-RI validation with in silico mixtures (v2: separate AUC windows)', ...
        'Position', [200 100 800 580], 'Resize', 'off');
    fig.UserData = struct('ok', false);
    fig.CloseRequestFcn = @(~,~) onCancel();

    uilabel(fig, 'Position', [20 545 760 24], 'FontWeight', 'bold', ...
        'Text', 'In silico mixtures of measured single-sensor data, unmixed blindly by NMF-RI');
    uilabel(fig, 'Position', [20 522 760 22], 'FontColor', [0.35 0.35 0.35], ...
        'Text', 'Folder = sensor time series (mouse subfolders; green: G file, red: R file).   File = control dye spectrum (TXT).');

    hG = pathRow(487, 'Green slot (folder or file)', 'green');
    hR = pathRow(452, 'Red slot (folder or file)',   'red');

    uilabel(fig, 'Position', [20 417 190 22], 'Text', 'Green name (component 1)');
    hGN = uieditfield(fig, 'text', 'Position', [220 417 130 22]);
    uilabel(fig, 'Position', [380 417 170 22], 'Text', 'Red name (component 2)');
    hRN = uieditfield(fig, 'text', 'Position', [540 417 130 22]);

    hC1  = pathRow(377, 'Component 1 CSV (green, A0)', 'csv');
    hC2  = pathRow(342, 'Component 2 CSV (red, A0)',   'csv');
    hOut = pathRow(302, 'Output folder',               'out');
    hOut.Value = defOut;

    uilabel(fig, 'Position', [20 262 190 22], 'Text', 'Wavelength range (nm)');
    uilabel(fig, 'Position', [220 262 40 22], 'Text', 'from');
    hW1 = uieditfield(fig, 'numeric', 'Position', [260 262 70 22], 'Value', opt.wlRange(1));
    uilabel(fig, 'Position', [345 262 25 22], 'Text', 'to');
    hW2 = uieditfield(fig, 'numeric', 'Position', [370 262 70 22], 'Value', opt.wlRange(2));

    uilabel(fig, 'Position', [20 222 300 22], 'Text', 'NMF-RI parameters', 'FontWeight', 'bold');
    uilabel(fig, 'Position', [20 192 220 22], 'Text', 'theta1 (intra-layer convergence)');
    hT1 = uieditfield(fig, 'numeric', 'Position', [260 192 90 22], 'Value', opt.theta1, 'ValueDisplayFormat', '%.6g');
    uilabel(fig, 'Position', [400 192 220 22], 'Text', 'theta2 (inter-layer convergence)');
    hT2 = uieditfield(fig, 'numeric', 'Position', [640 192 90 22], 'Value', opt.theta2, 'ValueDisplayFormat', '%.6g');
    uilabel(fig, 'Position', [20 160 230 22], 'Text', 'alphaA (spectral sparsity weight)');
    hAA = uieditfield(fig, 'numeric', 'Position', [260 160 90 22], 'Value', opt.alphaA, 'ValueDisplayFormat', '%.6g');
    uilabel(fig, 'Position', [400 160 235 22], 'Text', 'alphaX (concentration sparsity weight)');
    hAX = uieditfield(fig, 'numeric', 'Position', [640 160 90 22], 'Value', opt.alphaX, 'ValueDisplayFormat', '%.6g');

    uilabel(fig, 'Position', [20 115 200 22], 'Text', 'Pairing (sensor + sensor only)');
    hPair = uidropdown(fig, 'Position', [220 115 450 22], ...
        'Items', {'1:1 (animals sorted by name, each animal used once)', 'All combinations (every green with every red animal)'}, ...
        'ItemsData', {'one-to-one', 'all'}, 'Value', 'one-to-one');

    uibutton(fig, 'Position', [560 25 100 30], 'Text', 'Run', ...
        'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) onRun());
    uibutton(fig, 'Position', [680 25 100 30], 'Text', 'Cancel', ...
        'ButtonPushedFcn', @(~,~) onCancel());

    uiwait(fig);
    if isvalid(fig)
        if fig.UserData.ok
            cfg = fig.UserData.cfg;
        end
        delete(fig);
    end

    function h = pathRow(y, label, kind)
        uilabel(fig, 'Position', [20 y 200 22], 'Text', label);
        switch kind
            case {'green', 'red'}
                h = uieditfield(fig, 'text', 'Position', [220 y 360 22]);
                uibutton(fig, 'Position', [590 y 90 22], 'Text', 'Folder...', ...
                    'ButtonPushedFcn', @(~,~) onBrowse(h, 'dir', kind));
                uibutton(fig, 'Position', [690 y 90 22], 'Text', 'File (dye)...', ...
                    'ButtonPushedFcn', @(~,~) onBrowse(h, 'txt', kind));
            otherwise
                h = uieditfield(fig, 'text', 'Position', [220 y 460 22]);
                if strcmp(kind, 'csv'), k2 = 'csv'; else, k2 = 'dir'; end
                uibutton(fig, 'Position', [690 y 90 22], 'Text', 'Browse...', ...
                    'ButtonPushedFcn', @(~,~) onBrowse(h, k2, ''));
        end
    end

    function onBrowse(h, kind, slot)
        start = h.Value;
        if ~isempty(start) && exist(start, 'file') == 2
            start = fileparts(start);
        end
        if isempty(start) || exist(start, 'dir') ~= 7
            if strcmp(slot, '') && strcmp(kind, 'dir'), start = defOut; else, start = defRaw; end
        end
        if exist(start, 'dir') ~= 7, start = pwd; end
        p = '';
        switch kind
            case 'dir'
                d = uigetdir(start);
                if ischar(d), p = d; end
            case 'txt'
                [f, d] = uigetfile({'*.txt', 'Spectrum TXT (*.txt)'}, 'Select control dye spectrum TXT', start);
                if ischar(f), p = fullfile(d, f); end
            case 'csv'
                [f, d] = uigetfile({'*.csv', 'CSV files (*.csv)'}, 'Select emission spectrum CSV', start);
                if ischar(f), p = fullfile(d, f); end
        end
        if ~isempty(p)
            h.Value = p;
            nm = autoName(p, opt);
            if strcmp(slot, 'green') && isempty(strtrim(hGN.Value)), hGN.Value = nm; end
            if strcmp(slot, 'red')   && isempty(strtrim(hRN.Value)), hRN.Value = nm; end
        end
        figure(fig);
    end

    function onRun()
        c.greenPath = strtrim(hG.Value);
        c.redPath   = strtrim(hR.Value);
        c.greenName = strtrim(hGN.Value);
        c.redName   = strtrim(hRN.Value);
        c.comp1Csv  = strtrim(hC1.Value);
        c.comp2Csv  = strtrim(hC2.Value);
        c.outDir    = strtrim(hOut.Value);
        c.wlRange   = [hW1.Value hW2.Value];
        c.theta1    = hT1.Value;
        c.theta2    = hT2.Value;
        c.alphaA    = hAA.Value;
        c.alphaX    = hAX.Value;
        c.pairing   = hPair.Value;
        c = completeCfg(c, opt);
        msg = checkCfg(c);
        if ~isempty(msg)
            uialert(fig, msg, 'Invalid input');
            return;
        end
        fig.UserData = struct('ok', true, 'cfg', c);
        uiresume(fig);
    end

    function onCancel()
        fig.UserData = struct('ok', false);
        uiresume(fig);
    end
end


function nm = autoName(p, opt)
% Folder: folder name. File: a known sensor/dye name found in the file
% name (canonical spelling), otherwise the file name without extension.
    p = stripSep(p);
    [~, base] = fileparts(p);
    nm = base;
    if exist(p, 'file') == 2
        for k = 1:size(opt.sensorColors, 1)
            key = opt.sensorColors{k, 1};
            if ~isempty(regexpi(base, regexptranslate('escape', key), 'once'))
                nm = key;
                return;
            end
        end
    end
end


function c = completeCfg(c, opt)
    % accept the field names of the first version
    if ~isfield(c, 'greenPath') && isfield(c, 'greenDir'), c.greenPath = c.greenDir; end
    if ~isfield(c, 'redPath')   && isfield(c, 'redDir'),   c.redPath   = c.redDir;   end
    def = {'wlRange', opt.wlRange; 'theta1', opt.theta1; 'theta2', opt.theta2; ...
           'alphaA', opt.alphaA; 'alphaX', opt.alphaX; 'pairing', 'one-to-one'; ...
           'greenName', ''; 'redName', ''; 'greenPath', ''; 'redPath', ''; ...
           'comp1Csv', ''; 'comp2Csv', ''; 'outDir', ''};
    for k = 1:size(def, 1)
        if ~isfield(c, def{k, 1}) || isempty(c.(def{k, 1}))
            c.(def{k, 1}) = def{k, 2};
        end
    end
    c.greenPath = stripSep(c.greenPath);
    c.redPath   = stripSep(c.redPath);
    c.outDir    = stripSep(c.outDir);
    if isempty(c.greenName) && ~isempty(c.greenPath)
        c.greenName = autoName(c.greenPath, opt);
    end
    if isempty(c.redName) && ~isempty(c.redPath)
        c.redName = autoName(c.redPath, opt);
    end
end


function t = pathType(p)
    if isempty(p)
        t = '';
    elseif exist(p, 'dir') == 7
        t = 'sensor';
    elseif exist(p, 'file') == 2
        t = 'dye';
    else
        t = '';
    end
end


function msg = checkCfg(c)
    msg = '';
    p = [c.theta1 c.theta2 c.alphaA c.alphaX];
    gt = pathType(c.greenPath);
    rt = pathType(c.redPath);
    if isempty(gt)
        msg = ['Green slot: folder or file not found: ' c.greenPath];
    elseif isempty(rt)
        msg = ['Red slot: folder or file not found: ' c.redPath];
    elseif strcmp(gt, 'dye') && strcmp(rt, 'dye')
        msg = 'Both slots are files (dyes). At least one slot must be a sensor folder.';
    elseif isempty(c.comp1Csv) || exist(c.comp1Csv, 'file') ~= 2
        msg = ['Component 1 CSV not found: ' c.comp1Csv];
    elseif isempty(c.comp2Csv) || exist(c.comp2Csv, 'file') ~= 2
        msg = ['Component 2 CSV not found: ' c.comp2Csv];
    elseif isempty(c.outDir)
        msg = 'Please select an output folder.';
    elseif isempty(c.greenName) || isempty(c.redName) || strcmpi(c.greenName, c.redName)
        msg = 'Please enter two different names.';
    elseif numel(c.wlRange) ~= 2 || ~(c.wlRange(2) > c.wlRange(1))
        msg = 'Wavelength range: "to" must be larger than "from".';
    elseif any(~isfinite(p)) || any(p < 0)
        msg = 'NMF-RI parameters must be numbers >= 0.';
    elseif ~any(strcmp(c.pairing, {'one-to-one', 'all'}))
        msg = 'cfg.pairing must be ''one-to-one'' or ''all''.';
    elseif exist('NMF_RI', 'file') ~= 2
        msg = 'NMF_RI.m is not on the MATLAB path.';
    elseif exist('dataPreprocessing', 'file') ~= 2
        msg = 'dataPreprocessing.m is not on the MATLAB path.';
    end
end


%% =====================================================================
%  Main
%  =====================================================================
function S = runHybrid(cfg)
    opt = defaultSettings();
    cfg = completeCfg(cfg, opt);
    opt.wlRange = cfg.wlRange;
    msg = checkCfg(cfg);
    if ~isempty(msg)
        error(msg);
    end
    tTotal = tic;
    warn = {};

    types = {pathType(cfg.greenPath), pathType(cfg.redPath)};
    names = {cfg.greenName, cfg.redName};
    cols = {sensorColor(names{1}, opt, opt.colorComp1), sensorColor(names{2}, opt, opt.colorComp2)};
    isDyeMode = any(strcmp(types, 'dye'));
    slotLabel = names;
    for s = 1:2
        if strcmp(types{s}, 'dye'), slotLabel{s} = [names{s} '_dye']; end
    end
    if isDyeMode
        modeTxt = 'sensor + control dye';
        unit = 'Mixture';
    else
        modeTxt = 'sensor + sensor';
        unit = 'Mixture';
    end

    hybDir = fullfile(cfg.outDir, 'mixtures_v2');
    if exist(hybDir, 'dir') ~= 7
        mkdir(hybDir);
    end

    fprintf('\n==============================================================\n');
    fprintf(' NMF-RI validation with in silico mixtures  [v2: separate AUC windows for ground truth and unmixed]\n');
    fprintf(' Mode : %s\n', modeTxt);
    fprintf(' Green: %-10s (%s)  %s\n', names{1}, types{1}, cfg.greenPath);
    fprintf(' Red  : %-10s (%s)  %s\n', names{2}, types{2}, cfg.redPath);
    fprintf('==============================================================\n');

    % ----- sources
    tags = {opt.greenTag, opt.redTag};
    src = cell(1, 2);
    for s = 1:2
        if strcmp(types{s}, 'sensor')
            [src{s}, w] = listAnimals(ternary(s == 1, cfg.greenPath, cfg.redPath), tags{s});
            warn = [warn, w]; %#ok<AGROW>
            fprintf('%s mice (%s file): %d\n', names{s}, tags{s}, numel(src{s}));
            if isempty(src{s})
                error('No %s mice with a %s file found.', names{s}, tags{s});
            end
        end
    end

    % ----- mixtures (index into each slot's source list; 0 = dye)
    if ~isDyeMode
        nG = numel(src{1}); nR = numel(src{2});
        if strcmp(cfg.pairing, 'all')
            [jj, ii] = meshgrid(1:nR, 1:nG);
            ii = ii'; jj = jj';
            ii = ii(:); jj = jj(:);
        else
            n = min(nG, nR);
            ii = (1:n)'; jj = (1:n)';
            for k = n+1:nG
                warn{end+1} = sprintf('%s mouse %s not paired (1:1, fewer %s mice).', names{1}, src{1}(k).name, names{2}); %#ok<AGROW>
            end
            for k = n+1:nR
                warn{end+1} = sprintf('%s mouse %s not paired (1:1, fewer %s mice).', names{2}, src{2}(k).name, names{1}); %#ok<AGROW>
            end
        end
    elseif strcmp(types{1}, 'sensor')
        ii = (1:numel(src{1}))'; jj = zeros(size(ii));
    else
        jj = (1:numel(src{2}))'; ii = zeros(size(jj));
    end
    idx = {ii, jj};
    nHyb = numel(ii);

    fprintf('\n--------------------------------------------------------------\n');
    if isDyeMode
        fprintf(' TOTAL: %d mixtures  (one per sensor mouse)\n', nHyb);
    else
        fprintf(' TOTAL: %d mixtures  (pairing: %s)\n', nHyb, cfg.pairing);
    end
    fprintf('--------------------------------------------------------------\n');
    for p = 1:nHyb
        fprintf('   %s %2d: %s\n', unit, p, hybridTitle(src, idx, p, names, types));
    end

    % ----- progress window (optional; console output is always shown)
    used = cell(1, 2);
    nLoad = 0;
    for s = 1:2
        if strcmp(types{s}, 'sensor')
            used{s} = unique(idx{s});
        else
            used{s} = 1;
        end
        nLoad = nLoad + numel(used{s});
    end
    nSteps = nLoad + nHyb + 1;
    hw = [];
    try
        hw = waitbar(0, sprintf('Reading files (0/%d)', nLoad), ...
            'Name', sprintf('In silico mixtures: %d %ss', nHyb, lower(unit)));
    catch
    end
    barCleaner = onCleanup(@() closeBar(hw)); %#ok<NASGU>   % also closes it after an error

    % ----- read files
    fprintf('\nReading %d files\n', nLoad);
    step = 0;
    dye = cell(1, 2);
    for s = 1:2
        if strcmp(types{s}, 'dye')
            step = step + 1;
            progress(hw, step / nSteps, sprintf('Reading files (%d/%d): %s', step, nLoad, names{s}));
            fprintf('   [%d/%d]', step, nLoad);
            [dye{s}, w] = loadDye(cfg, s, names{s}, opt); warn = [warn, w]; %#ok<AGROW>
            continue;
        end
        lab = {'green', 'red'};
        for k = used{s}(:)'
            step = step + 1;
            progress(hw, step / nSteps, sprintf('Reading files (%d/%d): %s', step, nLoad, src{s}(k).name));
            fprintf('   [%d/%d]', step, nLoad);
            [src{s}(k), w] = loadAnimal(src{s}(k), opt, lab{s}); warn = [warn, w]; %#ok<AGROW>
        end
    end

    % ----- run every mixture
    hyb = struct([]);
    tRun = tic;
    for p = 1:nHyb
        ttl = hybridTitle(src, idx, p, names, types);
        msgP = sprintf('%s %d/%d: %s', unit, p, nHyb, ttl);
        fprintf('\n==== %s  (%.0f%% done', msgP, 100 * (p - 1) / nHyb);
        if p > 1
            el = toc(tRun);
            fprintf(', elapsed %.1f min, about %.1f min left)\n', el / 60, el / (p - 1) * (nHyb - p + 1) / 60);
        else
            fprintf(')\n');
        end
        progress(hw, (nLoad + p - 1) / nSteps, sprintf('%s - NMF-RI running', msgP));

        sl = cell(1, 2);
        for s = 1:2
            if strcmp(types{s}, 'dye')
                sl{s} = dye{s};
            else
                sl{s} = src{s}(idx{s}(p));
            end
        end
        if (strcmp(types{1}, 'sensor') && sl{1}.nTrials < 1) || (strcmp(types{2}, 'sensor') && sl{2}.nTrials < 1)
            warn{end+1} = sprintf('%s %d (%s) skipped (no complete trials).', unit, p, ttl); %#ok<AGROW>
            continue;
        end
        [H, w] = runOne(sl, types, names, cfg, opt);
        warn = [warn, w]; %#ok<AGROW>
        H.index = p;
        H.nTotal = nHyb;
        H.unit = unit;
        H.title = ttl;
        H.label = sprintf('%s_%s+%s_%s', names{1}, sl{1}.name, names{2}, sl{2}.name);
        if strcmp(types{1}, 'dye'), H.label = sprintf('%s_dye+%s_%s', names{1}, names{2}, sl{2}.name); end
        if strcmp(types{2}, 'dye'), H.label = sprintf('%s_%s+%s_dye', names{1}, sl{1}.name, names{2}); end
        if isempty(hyb)
            hyb = H;
        else
            hyb(end+1) = H; %#ok<AGROW>
        end
    end
    fprintf('\nAll %d %ss unmixed (%.1f min).\n', numel(hyb), lower(unit), toc(tRun) / 60);
    if isempty(hyb)
        error('No mixture could be processed.');
    end

    % ----- figures
    progress(hw, (nSteps - 1) / nSteps, 'Saving figures and Excel');
    ylReward = rewardLimits(hyb, opt);
    fprintf('\nSaving figures\n');
    for p = 1:numel(hyb)
        progress(hw, (nSteps - 1) / nSteps, sprintf('Saving figure %d/%d', p, numel(hyb)));
        file = fullfile(hybDir, sprintf('%s%02d_%s_mixture_v2.pdf', unit, hyb(p).index, safeFileName(hyb(p).label)));
        hyb(p).pdf = hybridFigure(hyb(p), names, types, cols, ylReward, modeTxt, opt, file);
    end

    Sm = summarise(hyb);
    base = safeFileName([slotLabel{1} '+' slotLabel{2}]);
    sumFile = summaryFigure(Sm, names, types, cols, cfg, modeTxt, unit, opt, ...
        fullfile(cfg.outDir, [base '_mixture_summary_v2.pdf']));
    cmpFile = aucPeakFigure(hyb, names, cols, modeTxt, unit, opt, ...
        fullfile(cfg.outDir, [base '_mixture_AUC_peak_v2.pdf']));

    xlsFile = writeExcel(hyb, Sm, names, types, cfg, opt, unit, fullfile(cfg.outDir, [base '_mixture_R2_v2.xlsx']));
    matFile = saveResultsMat(hyb, names, types, cols, cfg, opt, modeTxt, unit, fullfile(cfg.outDir, [base '_mixture_R2_v2.mat']));
    progress(hw, 1, 'Done');

    S.inputs = cfg;
    S.mode = modeTxt;
    S.settings = opt;
    S.hybrids = rmfield(hyb, 'wl');
    S.summary = Sm;
    S.files = struct('summaryPdf', sumFile, 'aucPeakPdf', cmpFile, 'excel', xlsFile, 'mat', matFile, 'hybridFolder', hybDir);
    S.warnings = warn(:);
    S.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    fprintf('\nR^2 across %ss (mean +/- SEM, n = %d)\n', lower(unit), numel(hyb));
    for q = 1:numel(Sm.panel)
        fprintf('   %-22s %.4f +/- %.4f\n', Sm.panel(q).metric, Sm.panel(q).r2Mean, Sm.panel(q).r2Sem);
    end
    if isDyeMode
        fprintf('   %-22s %.2f +/- %.2f %%\n', 'dye delta leakage', Sm.leakMean, Sm.leakSem);
    end
    if isempty(warn)
        fprintf('\nNo warnings.\n');
    else
        fprintf('\nWarnings (%d):\n', numel(warn));
        fprintf('   %s\n', warn{:});
    end
    fprintf('Total time: %.1f min\n', toc(tTotal) / 60);
    fprintf('Output folder: %s\n', cfg.outDir);
end


function t = hybridTitle(src, idx, p, names, types)
    part = cell(1, 2);
    for s = 1:2
        if strcmp(types{s}, 'dye')
            part{s} = [names{s} ' (dye)'];
        else
            part{s} = [names{s} ' ' src{s}(idx{s}(p)).name];
        end
    end
    t = [part{1} ' + ' part{2}];
end


%% =====================================================================
%  Sources
%  =====================================================================
function A = emptySource()
    A = struct('name', '', 'file', '', 'type', 'sensor', 'wl', [], 'Y', [], 'spec', [], ...
        'nTrials', 0, 'nSamples', 0, 'nRows', 0, 'fileMin', NaN);
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


function [A, w] = loadDye(cfg, s, name, opt)
% Control dye: all rows of the TXT averaged into one static spectrum;
% the minimum of the mean spectrum (450-650 nm) is subtracted (same rule as sensors).
    w = {};
    if s == 1, f = cfg.greenPath; else, f = cfg.redPath; end
    [wl, Y, nDropped, fmt] = readDyeTxt(f, opt.wlRange);
    if isempty(Y) || isempty(wl)
        error('No spectrum values within %g-%g nm in dye file %s', opt.wlRange(1), opt.wlRange(2), f);
    end
    if nDropped > 0
        w{end+1} = sprintf('%s (dye): %d incomplete rows removed.', name, nDropped);
    end
    A = emptySource();
    A.name = name;
    A.file = f;
    A.type = 'dye';
    A.wl = wl;
    A.nRows = size(Y, 2);
    A.nTrials = Inf;
    % static spectrum = mean of all rows; offset = its minimum across
    % wavelengths (same rule as the sensor data: min of the time-averaged spectrum)
    sp = mean(Y, 2);
    A.fileMin = min(sp);
    A.spec = sp - A.fileMin;
    [~, fn, fe] = fileparts(f);
    if strcmp(fmt, 'column')
        fprintf(' dye    %-14s %-28s single spectrum (wavelength / value columns), offset %.1f counts, peak %.1f counts\n', ...
            name, [fn fe], A.fileMin, max(A.spec));
    else
        fprintf(' dye    %-14s %-28s %d rows averaged, offset %.1f counts, peak %.1f counts\n', ...
            name, [fn fe], A.nRows, A.fileMin, max(A.spec));
    end
    if strcmp(fmt, 'row') && A.nRows >= opt.trialLength
        fprintf('        note: the dye file has %d rows (>= 1 trial); treated as a static spectrum (mean of all rows).\n', A.nRows);
    end
end


%% =====================================================================
%  One mixture: mixture, NMF-RI, comparison
%  =====================================================================
function [H, w] = runOne(sl, types, names, cfg, opt)
    w = {};
    t0 = tic;
    L = opt.trialLength;
    nTr = min(sl{1}.nTrials, sl{2}.nTrials);
    nS = nTr * L;
    if strcmp(types{1}, 'sensor'), wl = sl{1}.wl(:); else, wl = sl{2}.wl(:); end

    Yt = cell(1, 2);
    for s = 1:2
        [Yt{s}, interp] = sourceMatrix(sl{s}, wl, nS);
        if interp
            w{end+1} = sprintf('%s: wavelength axis of %s differs; interpolated.', sl{1}.name, names{s}); %#ok<AGROW>
        end
    end
    Y = Yt{1} + Yt{2};     % hybrid mixture

    % ----- NMF-RI, as NMF_RI_Photometry
    a1 = loadEmissionSpectrum(cfg.comp1Csv, wl);
    a2 = loadEmissionSpectrum(cfg.comp2Csv, wl);
    A0 = [a1(:), a2(:)];
    A0(A0 <= 0) = 1e-4;
    A0 = A0 * diag(1 ./ (sum(A0, 1) + eps));

    [YPrepared, ~, ~] = dataPreprocessing(A0, Y);
    prepOffset = mean(min(Y, [], 1));          % mean per-timepoint minimum removed
    H0 = max(100 * eps, pinv(A0) * YPrepared);
    [ANmf, nIter] = callNmf(YPrepared, A0, H0, cfg);
    HFull = max(1e6 * eps, pinv(ANmf' * ANmf) * ANmf' * YPrepared);
    clear YPrepared H0 Y;

    % ----- per-trial reward and delta spectra
    rw = opt.rewardSamples;
    bl = opt.baselineSamples;
    tRw = cell(1, 2); tDl = cell(1, 2); cRw = cell(1, 2); cDl = cell(1, 2);
    hRw = trialMeans(HFull, rw, nTr, L);
    hBl = trialMeans(HFull, bl, nTr, L);
    for s = 1:2
        tRw{s} = trialMeans(Yt{s}, rw, nTr, L);
        tDl{s} = tRw{s} - trialMeans(Yt{s}, bl, nTr, L);
        cRw{s} = ANmf(:, s) * hRw(s, :);
        cDl{s} = ANmf(:, s) * (hRw(s, :) - hBl(s, :));
    end

    sumT = sprintf('%s + %s', names{1}, names{2});
    pan = [ ...
        makePanel('reward_sum',            1, 1, tRw{1} + tRw{2}, cRw{1}, cRw{2}, false, {tRw{1}, tRw{2}}), ...
        makePanel(['reward_' names{1}],    1, 2, tRw{1},          cRw{1}, [],     false), ...
        makePanel(['reward_' names{2}],    1, 3, tRw{2},          [],     cRw{2}, false)];
    leak = NaN;
    if ~any(strcmp(types, 'dye'))
        pan = [pan, ...
            makePanel('delta_sum',         2, 1, tDl{1} + tDl{2}, cDl{1}, cDl{2}, true, {tDl{1}, tDl{2}}), ...
            makePanel(['delta_' names{1}], 2, 2, tDl{1},          cDl{1}, [],     true), ...
            makePanel(['delta_' names{2}], 2, 3, tDl{2},          [],     cDl{2}, true)];
    else
        sS = find(strcmp(types, 'sensor'));   % sensor slot
        sD = 3 - sS;                          % dye slot
        if sS == 1
            pan = [pan, makePanel(['delta_' names{1}], 2, 2, tDl{1}, cDl{1}, [], true)];
        else
            pan = [pan, makePanel(['delta_' names{2}], 2, 3, tDl{2}, [], cDl{2}, true)];
        end
        aT = trapz(wl, mean(tDl{sS}, 2));
        aD = trapz(wl, mean(cDl{sD}, 2));
        if aT ~= 0
            leak = 100 * aD / aT;
        end
    end

    H.sumTitle = sumT;
    H.green = sl{1}.name;
    H.red = sl{2}.name;
    H.greenFile = sl{1}.file;
    H.redFile = sl{2}.file;
    H.nTrialsGreen = sl{1}.nTrials;
    H.nTrialsRed = sl{2}.nTrials;
    H.nTrials = nTr;
    H.nIter = nIter;
    H.prepOffset = prepOffset;
    H.leakagePct = leak;
    H.cmp = aucPeakCompare(wl, {mean(tRw{1}, 2), mean(tRw{2}, 2)}, ...
        {mean(cRw{1}, 2), mean(cRw{2}, 2)}, {mean(tDl{1}, 2), mean(tDl{2}, 2)}, ...
        {mean(cDl{1}, 2), mean(cDl{2}, 2)}, types, names, opt, tRw);
    H.A = ANmf;
    H.wl = wl;
    H.panel = pan;
    H.pdf = '';
    H.seconds = toc(t0);

    if isfinite(sl{1}.nTrials) && isfinite(sl{2}.nTrials) && sl{1}.nTrials ~= sl{2}.nTrials
        fprintf('   trials: %s %d, %s %d -> first %d summed\n', names{1}, sl{1}.nTrials, names{2}, sl{2}.nTrials, nTr);
    end
    fprintf('   n = %d trials, NMF-RI %d iterations, %.1f s\n', nTr, nIter, H.seconds);
    for q = 1:numel(pan)
        fprintf('   R^2 %-22s %.4f +/- %.4f\n', pan(q).metric, pan(q).r2Mean, pan(q).r2Sem);
    end
    if isfinite(leak)
        fprintf('   dye delta leakage      %.2f %%\n', leak);
    end
    for q = 1:numel(H.cmp.m)
        mq = H.cmp.m(q);
        fprintf('   %-24s ground truth %10.4g | unmixed %10.4g\n', mq.label, mq.raw, mq.nmf);
    end
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


function Q = makePanel(metric, row, col, truthTr, c1Tr, c2Tr, isDelta, parts)
% Per-trial R^2 (counts) of truth vs NMF; trial means in counts and normalised.
% parts (sum panels only): {green, red} measured single-sensor spectra per trial;
% the figure shows the mixture as these two stacked (no NMF, no R^2 in the figure;
% the R^2 is still computed for the Excel file).
    pred = zeros(size(truthTr));
    if ~isempty(c1Tr), pred = pred + c1Tr; end
    if ~isempty(c2Tr), pred = pred + c2Tr; end
    n = size(truthTr, 2);
    r2 = nan(n, 1);
    for k = 1:n
        y = truthTr(:, k);
        ssTot = sum((y - mean(y)) .^ 2);
        if ssTot > 0
            r2(k) = 1 - sum((y - pred(:, k)) .^ 2) / ssTot;
        end
    end
    Q.metric = metric;
    Q.row = row;
    Q.col = col;
    Q.isDelta = isDelta;
    Q.r2Trials = r2;
    [Q.r2Mean, Q.r2Sem, Q.r2N] = meanSem(r2);
    Q.truth = mean(truthTr, 2);
    Q.c1 = [];
    Q.c2 = [];
    if ~isempty(c1Tr), Q.c1 = mean(c1Tr, 2); end
    if ~isempty(c2Tr), Q.c2 = mean(c2Tr, 2); end
    if isDelta
        p = max(abs(Q.truth));
    else
        p = max(Q.truth);
    end
    if ~(p > 0)
        p = NaN;
    end
    Q.peak = p;                % counts; normalisation value
    Q.truthN = Q.truth / p;
    Q.c1N = Q.c1 / p;
    Q.c2N = Q.c2 / p;
    Q.isSum = nargin >= 8 && ~isempty(parts);
    Q.p1 = []; Q.p2 = []; Q.p1N = []; Q.p2N = [];
    if Q.isSum
        Q.p1 = mean(parts{1}, 2);
        Q.p2 = mean(parts{2}, 2);
        Q.p1N = Q.p1 / p;
        Q.p2N = Q.p2 / p;
    end
end


function [m, s, n] = meanSem(v)
    v = v(isfinite(v));
    n = numel(v);
    if n == 0
        m = NaN; s = NaN;
    elseif n == 1
        m = v; s = 0;
    else
        m = mean(v); s = std(v) / sqrt(n);
    end
end


%% =====================================================================
%  Summary across mixtures
%  =====================================================================
function Sm = summarise(hyb)
    wl = hyb(1).wl;
    Sm.wl = wl;
    Sm.n = numel(hyb);
    for q = 1:numel(hyb(1).panel)
        T = []; C1 = []; C2 = []; P1 = []; P2 = []; r2 = [];
        for i = 1:numel(hyb)
            Q = hyb(i).panel(q);
            r2(end+1) = Q.r2Mean; %#ok<AGROW>
            if ~isfinite(Q.peak), continue; end
            T(:, end+1) = onAxis(hyb(i).wl, Q.truthN, wl); %#ok<AGROW>
            if ~isempty(Q.c1N), C1(:, end+1) = onAxis(hyb(i).wl, Q.c1N, wl); end %#ok<AGROW>
            if ~isempty(Q.c2N), C2(:, end+1) = onAxis(hyb(i).wl, Q.c2N, wl); end %#ok<AGROW>
            if isfield(Q, 'isSum') && Q.isSum
                P1(:, end+1) = onAxis(hyb(i).wl, Q.p1N, wl); %#ok<AGROW>
                P2(:, end+1) = onAxis(hyb(i).wl, Q.p2N, wl); %#ok<AGROW>
            end
        end
        Q1 = hyb(1).panel(q);
        D.metric = Q1.metric;
        D.row = Q1.row;
        D.col = Q1.col;
        D.isDelta = Q1.isDelta;
        if isempty(T)
            D.truth = nan(size(wl)); D.truthSem = zeros(size(wl));
        else
            D.truth = mean(T, 2);
            if size(T, 2) > 1
                D.truthSem = std(T, 0, 2) / sqrt(size(T, 2));
            else
                D.truthSem = zeros(size(wl));
            end
        end
        D.c1 = [];
        D.c2 = [];
        if ~isempty(C1), D.c1 = mean(C1, 2); end
        if ~isempty(C2), D.c2 = mean(C2, 2); end
        D.isSum = ~isempty(P1);
        D.p1 = []; D.p2 = [];
        if D.isSum, D.p1 = mean(P1, 2); D.p2 = mean(P2, 2); end
        [D.r2Mean, D.r2Sem, D.r2N] = meanSem(r2);
        if q == 1
            Sm.panel = D;
        else
            Sm.panel(q) = D;
        end
    end
    [Sm.leakMean, Sm.leakSem, Sm.leakN] = meanSem([hyb.leakagePct]);
end


%% =====================================================================
%  Figures
%  =====================================================================
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


function file = aucPeakFigure(hyb, names, cols, modeTxt, unit, opt, file)
% Raw vs NMF per mixture (paired). Row 1 AUC / ratio, row 2 delta AUC (below the
% matching AUC), row 3 peak (below the matching sensor).
    M = hyb(1).cmp.m;
    if isempty(M)
        file = '';
        return;
    end
    rowsUsed = unique([M.row]);
    nRow = numel(rowsUsed);
    nCol = max([M.col]);
    pw = 3.2; ph = 4.0; gx = 2.0; gy = 2.8;
    w = max(gx + nCol * (pw + gx), 20);          % wide enough for the header
    x0 = (w - nCol * (pw + gx) + gx) / 2;          % panels centred
    h = 1.2 + nRow * ph + (nRow - 1) * gy + 3.8;
    f = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 w h], ...
        'Name', 'AUC and peak');
    set(f, 'UserData', [w h]);
    for q = 1:numel(M)
        r = find(rowsUsed == M(q).row);
        y0 = 1.2 + (nRow - r) * (ph + gy);
        raw = arrayfun(@(H) H.cmp.m(q).raw, hyb);
        nmf = arrayfun(@(H) H.cmp.m(q).nmf, hyb);
        ax = axes(f, 'Units', 'centimeters', 'Position', [x0 + (M(q).col - 1) * (pw + gx) y0 pw ph]);
        set(ax, 'Units', 'normalized');
        drawPaired(ax, raw(:), nmf(:), cols{M(q).slot}, M(q), opt);
    end
    head = { ...
        sprintf('\\bfGround truth vs unmixed (v2, separate AUC windows):  %s + %s\\rm', tex(names{1}), tex(names{2})), ...
        sprintf('%s in silico mixtures, n = %d %ss; reward and \\Delta (reward - baseline) spectra, trial means; one line per %s', ...
            modeTxt, numel(hyb), lower(unit), lower(unit)), ...
        sprintf('v2: AUC and \\DeltaAUC window = own peak %s %g nm (ground truth: ground-truth peak, unmixed: unmixed peak; T1: green %d-%d nm, T2: red %d-%d nm); mean %s SEM', ...
            char(177), opt.aucHalfWidth, opt.bandT1(1), opt.bandT1(2), opt.bandT2(1), opt.bandT2(2), char(177))};
    figHeader(f, head, opt);
    file = savePdf(f, file);
end


function drawPaired(ax, raw, nmf, col, Mq, opt)
    hold(ax, 'on');
    for i = 1:numel(raw)
        plot(ax, [1 2], [raw(i) nmf(i)], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);
    end
    plot(ax, ones(size(raw)), raw, 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', [0.25 0.25 0.25], 'MarkerEdgeColor', 'none');
    plot(ax, 2 * ones(size(nmf)), nmf, 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none');
    [mr, sr] = meanSem(raw);
    [mn, sn] = meanSem(nmf);
    plot(ax, [0.72 0.88], [mr mr], 'k-', 'LineWidth', 2);
    plot(ax, [0.8 0.8], [mr - sr, mr + sr], 'k-', 'LineWidth', 1);
    plot(ax, [2.12 2.28], [mn mn], '-', 'Color', col, 'LineWidth', 2);
    plot(ax, [2.2 2.2], [mn - sn, mn + sn], '-', 'Color', col, 'LineWidth', 1);
    hold(ax, 'off');

    v = [raw; nmf];
    v = v(isfinite(v));
    switch Mq.kind
        case 'peak'
            ylab = 'peak wavelength (nm)';
            d = nmf - raw;
            [md, sd] = meanSem(d);
            dtxt = sprintf('unmixed - ground truth: %+.2f %s %.2f nm', md, char(177), sd);
            if isempty(v), v = [0 1]; end
            lo = min(v) - 1; hi = max(v) + 1;
            [~, st] = countLimits([0; hi - lo]);
            yl = [floor(lo / st) * st, ceil(hi / st) * st];
        otherwise
            switch Mq.kind
                case 'ratio', ylab = 'T2 / T1 AUC ratio';
                case 'dauc',  ylab = '\Delta AUC (counts x nm)';
                otherwise,    ylab = 'AUC (counts x nm)';
            end
            d = 100 * (nmf - raw) ./ raw;
            [md, sd] = meanSem(d);
            dtxt = sprintf('unmixed - ground truth: %+.1f %s %.1f %%', md, char(177), sd);
            [yl, st] = countLimits([v; 1.1 * max([v; 0])]);
    end
    set(ax, 'XLim', [0.5 2.5], 'XTick', [1 2], 'XTickLabel', {'ground truth', 'unmixed'}, ...
        'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.8, 'FontName', opt.fontName, ...
        'FontSize', opt.fontSize, 'XColor', 'k', 'YColor', 'k');
    if numel(yl) == 2 && all(isfinite(yl)) && yl(2) > yl(1)
        try
            set(ax, 'YLim', yl, 'YTick', tickVector(yl, st));
        catch
            set(ax, 'YLimMode', 'auto', 'YTickMode', 'auto');
        end
    end
    ylabel(ax, ylab, 'FontName', opt.fontName, 'FontSize', opt.fontSize);
    [tt, df, pp] = pairedT(raw, nmf);
    text(ax, 0.5, 1.26, colorName(Mq.label, col), 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontName', opt.fontName, 'FontSize', opt.fontSize, 'Clipping', 'off');
    text(ax, 0.5, 1.13, dtxt, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontName', opt.fontName, 'FontSize', opt.fontSize - 1, 'Clipping', 'off');
    text(ax, 0.5, 1.02, pairedText(tt, df, pp), 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontName', opt.fontName, 'FontSize', opt.fontSize - 1, 'Clipping', 'off');
end


function [t, df, p, n] = pairedT(a, b)
% Paired t-test of b vs a (two-sided), without the Statistics Toolbox:
% t = mean(d) / (SD(d) / sqrt(n)), d = b - a, df = n - 1,
% p = betainc(df / (df + t^2), df/2, 1/2).
    ok = isfinite(a) & isfinite(b);
    d = b(ok) - a(ok);
    n = numel(d);
    df = n - 1;
    t = NaN; p = NaN;
    if n < 2
        return;
    end
    sd = std(d);
    if sd == 0
        if mean(d) == 0
            t = 0; p = 1;               % identical in every pair
        end
        return;                          % constant non-zero difference: t undefined
    end
    t = mean(d) / (sd / sqrt(n));
    p = betainc(df / (df + t ^ 2), df / 2, 0.5);
end


function s = pairedText(t, df, p)
    if ~isfinite(p)
        s = 'paired t-test: n/a';
    elseif p < 0.001
        s = sprintf('paired t(%d) = %.2f, p < 0.001', df, t);
    else
        s = sprintf('paired t(%d) = %.2f, p = %.3f', df, t, p);
    end
end


function yl = rewardLimits(hyb, opt)
% Normalised reward row: one y range shared by all mixture figures.
    v = [];
    for i = 1:numel(hyb)
        for q = 1:numel(hyb(i).panel)
            Q = hyb(i).panel(q);
            if Q.row == 1
                v = [v; panelValues(Q.truthN, [], Q.c1N, Q.c2N)]; %#ok<AGROW>
            end
        end
    end
    yl = niceLimits(v, opt);
end


function v = panelValues(t, s, c1, c2)
    v = t(:);
    if ~isempty(s)
        v = [v; t(:) + s(:); t(:) - s(:)];
    end
    if ~isempty(c1) && ~isempty(c2)
        v = [v; c1(:); c2(:); max(c1(:), 0) + max(c2(:), 0)];
    elseif ~isempty(c1)
        v = [v; c1(:)];
    elseif ~isempty(c2)
        v = [v; c2(:)];
    end
    v = v(isfinite(v));
end


function yl = niceLimits(v, opt)
% Normalised axes: fixed range opt.yNormLim (default -0.05 to 1).
% Values outside the range are clipped by the axes; reported in the console.
    yl = opt.yNormLim;
    if isempty(v)
        return;
    end
    nOut = sum(v > yl(2) + 1e-9 | v < yl(1) - 1e-9);
    if nOut > 0
        fprintf('   [axis] %d normalised values outside %g..%g (range %.3f..%.3f) are clipped in the figure.\n', ...
            nOut, yl(1), yl(2), min(v), max(v));
    end
end


function [yl, tick] = countLimits(v)
% Count axes: 3-5 ticks of 1/2/2.5/5 x 10^k; 0 always included.
    if isempty(v)
        yl = [0 1]; tick = 0.5;
        return;
    end
    hi = max(0, max(v));
    lo = min(0, min(v));
    if lo > -0.03 * hi, lo = 0; end      % ignore tiny negative noise
    if hi <= lo, hi = lo + 1; end
    raw = (hi - lo) / 4;
    p10 = 10 ^ floor(log10(raw));
    m = [1 2 2.5 5 10] * p10;
    tick = m(find(m >= raw, 1));
    if isempty(tick) || ~isfinite(tick) || tick <= 0
        tick = m(end);
    end
    yl = [floor(lo / tick + 1e-9) * tick, ceil(hi / tick - 1e-9) * tick];
    if yl(2) <= yl(1), yl(2) = yl(1) + tick; end
end


function file = hybridFigure(H, names, types, cols, ylReward, modeTxt, opt, file)
    head = { ...
        sprintf('\\bf%s %d/%d:  %s\\rm', H.unit, H.index, H.nTotal, tex(H.title)), ...
        sprintf('%s in silico mixture, n = %d trials, NMF-RI %d iterations', modeTxt, H.nTrials, H.nIter), ...
        sprintf(['first panel: in silico mixture = measured single-sensor spectra, stacked; ' ...
                 'other panels: black line = ground truth, fill = NMF-RI component']), ...
        sprintf('reward: normalised to the mixture / ground-truth peak; \\Delta: counts; R^2 per trial, mean %s SEM', char(177))};
    f = gridFigure(opt, [H.label ' - mixture']);
    for q = 1:numel(H.panel)
        Q = H.panel(q);
        if Q.row == 1
            D = struct('truth', Q.truthN, 'truthSem', [], 'c1', Q.c1N, 'c2', Q.c2N, ...
                'r2Mean', Q.r2Mean, 'r2Sem', Q.r2Sem, 'isSum', Q.isSum, 'p1', Q.p1N, 'p2', Q.p2N);
            yl = ylReward; yt = opt.yTick;
            ylab = '';
            if Q.col == 1, ylab = 'reward spectrum (norm.)'; end
        else
            D = struct('truth', Q.truth, 'truthSem', [], 'c1', Q.c1, 'c2', Q.c2, ...
                'r2Mean', Q.r2Mean, 'r2Sem', Q.r2Sem, 'isSum', Q.isSum, 'p1', Q.p1, 'p2', Q.p2);
            [yl, yt] = countLimits(panelValues(Q.truth, [], Q.c1, Q.c2));
            ylab = '\Delta intensity (counts)';
        end
        ax = gridAxes(f, Q.row, Q.col, opt);
        drawPanel(ax, H.wl, D, yl, yt, panelTitle(Q, names, types, cols), ylab, cols, opt);
    end
    figHeader(f, head, opt);
    file = savePdf(f, file);
end


function file = summaryFigure(Sm, names, types, cols, cfg, modeTxt, unit, opt, file)
    if strcmp(modeTxt, 'sensor + sensor')
        if strcmp(cfg.pairing, 'all'), pr = ', all combinations'; else, pr = ', 1:1 pairing'; end
    else
        pr = '';
    end
    t = cell(1, 2);
    for s = 1:2
        if strcmp(types{s}, 'dye'), t{s} = [tex(names{s}) ' (dye)']; else, t{s} = tex(names{s}); end
    end
    head = { ...
        sprintf('\\bfAll %ss:  %s + %s\\rm', lower(unit), t{1}, t{2}), ...
        sprintf('%s in silico mixtures, n = %d %ss%s', modeTxt, Sm.n, lower(unit), pr), ...
        sprintf('first panel: in silico mixture = measured single-sensor spectra, stacked (mean across %ss)', lower(unit)), ...
        sprintf('other panels: black line = ground truth, mean %s SEM across %ss (grey); fill: NMF-RI, mean; R^2 mean %s SEM across %ss', ...
                 char(177), lower(unit), char(177), lower(unit)), ...
        'normalised to the ground truth (reward: truth peak = 1; \Delta: max |truth| = 1)'};
    f = gridFigure(opt, 'mixture summary');
    yl = zeros(2, 2);
    for row = 1:2
        v = [];
        for q = 1:numel(Sm.panel)
            D = Sm.panel(q);
            if D.row == row
                v = [v; panelValues(D.truth, D.truthSem, D.c1, D.c2)]; %#ok<AGROW>
            end
        end
        yl(row, :) = niceLimits(v, opt);
    end
    firstDelta = true;
    for q = 1:numel(Sm.panel)
        D = Sm.panel(q);
        ylab = '';
        if D.row == 1 && D.col == 1
            ylab = 'reward spectrum (norm.)';
        elseif D.row == 2 && firstDelta
            ylab = '\Delta spectrum (norm.)';
            firstDelta = false;
        end
        ax = gridAxes(f, D.row, D.col, opt);
        drawPanel(ax, Sm.wl, D, yl(D.row, :), opt.yTick, panelTitle(D, names, types, cols), ylab, cols, opt);
    end
    figHeader(f, head, opt);
    file = savePdf(f, file);
end


function t = panelTitle(Q, names, types, cols)
    if Q.col == 1
        t = [colorName(names{1}, cols{1}) ' + ' colorName(names{2}, cols{2})];
        return;
    end
    s = Q.col - 1;
    if strcmp(types{s}, 'dye')
        t = [colorName(names{s}, cols{s}) ' (dye)'];
    else
        t = [colorName(names{s}, cols{s}) '-only'];
    end
end


function t = colorName(name, col)
    t = sprintf('\\color[rgb]{%.3f %.3f %.3f}%s\\color[rgb]{0 0 0}', col(1), col(2), col(3), tex(name));
end


function f = gridFigure(opt, name)
    pw = opt.panelCm(1); ph = opt.panelCm(2);
    w = 1.6 + 3 * (pw + 1.6);
    h = 1.3 + 2 * ph + 2.6 + 3.6;
    f = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 w h], 'Name', name);
    set(f, 'UserData', [w h]);
end


function ax = gridAxes(f, row, col, opt)
    pw = opt.panelCm(1); ph = opt.panelCm(2);
    x0 = 1.6 + (col - 1) * (pw + 1.6);
    if row == 1
        y0 = 1.3 + ph + 2.6;
    else
        y0 = 1.3;
    end
    ax = axes(f, 'Units', 'centimeters', 'Position', [x0 y0 pw ph]);
    set(ax, 'Units', 'normalized');
end


function drawPanel(ax, x, D, yl, ytick, ttl, ylab, cols, opt)
    x = x(:);
    z = zeros(size(x));
    hold(ax, 'on');
    if yl(1) < 0
        plot(ax, [x(1) x(end)], [0 0], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    end
    sumMode = isfield(D, 'isSum') && D.isSum && ~isempty(D.p1);
    if sumMode
        % in silico mixture: measured green-sensor spectrum from 0, red sensor stacked on top
        y1 = D.p1(:);
        y12 = y1 + D.p2(:);
        fill(ax, [x; flipud(x)], [y12; flipud(y1)], cols{2}, 'EdgeColor', 'none');
        fill(ax, [x; flipud(x)], [y1; z], cols{1}, 'EdgeColor', 'none');
    elseif ~isempty(D.c1) && ~isempty(D.c2)
        % stacked: component 1 from 0, component 2 on top (positive parts)
        y1 = max(D.c1(:), 0);
        y12 = y1 + max(D.c2(:), 0);
        fill(ax, [x; flipud(x)], [y12; flipud(y1)], cols{2}, 'EdgeColor', 'none');
        fill(ax, [x; flipud(x)], [y1; z], cols{1}, 'EdgeColor', 'none');
    elseif ~isempty(D.c1)
        fill(ax, [x; flipud(x)], [D.c1(:); z], cols{1}, 'EdgeColor', 'none');
    elseif ~isempty(D.c2)
        fill(ax, [x; flipud(x)], [D.c2(:); z], cols{2}, 'EdgeColor', 'none');
    end
    if ~sumMode && ~isempty(D.truthSem) && any(D.truthSem(:) > 0)
        t = D.truth(:); s = D.truthSem(:);
        fill(ax, [x; flipud(x)], [t + s; flipud(t - s)], [0 0 0], ...
            'EdgeColor', 'none', 'FaceAlpha', opt.semAlpha);
    end
    if ~sumMode
        plot(ax, x, D.truth, 'k-', 'LineWidth', 2);
    end
    hold(ax, 'off');
    set(ax, 'XLim', double(opt.wlRange), ...
        'XTick', tickVector(opt.wlRange, opt.xTickSpectrum), ...
        'Box', 'off', 'TickDir', 'out', 'Layer', 'top', 'LineWidth', 0.8, ...
        'FontName', opt.fontName, 'FontSize', opt.fontSize, 'XColor', 'k', 'YColor', 'k');
    allv = [D.truth(:); D.c1(:); D.c2(:)];
    allv = allv(isfinite(allv));
    if isempty(allv) || max(abs(allv)) < 1e-6
        fprintf('   [note] panel "%s": all values are ~0 (max |value| = %g); check this recording.\n', ...
            regexprep(ttl, '\\color\[rgb\]\{[^}]*\}', ''), max([abs(allv); 0]));
    end
    yl = double(yl(:)');
    ytick = double(ytick);
    ok = numel(yl) == 2 && all(isfinite(yl)) && yl(2) > yl(1) && isscalar(ytick) && isfinite(ytick) && ytick > 0;
    why = 'invalid limits or tick';
    if ok
        try
            set(ax, 'YLim', yl);
            set(ax, 'YTick', tickVector(yl, ytick));
        catch ME
            ok = false;
            why = ME.message;
        end
    end
    if ~ok
        % diagnostic: report this panel's data and use automatic y axis
        vr = @(v) sprintf('%g..%g (%d non-finite)', min(v(isfinite(v))), max(v(isfinite(v))), sum(~isfinite(v(:))));
        fprintf('   [axis] panel "%s": y limits [%s], tick %s -> automatic y axis (%s)\n', ...
            regexprep(ttl, '\\color\[rgb\]\{[^}]*\}', ''), num2str(yl, '%g '), num2str(ytick, '%g'), why);
        fprintf('          truth %s', vr(D.truth(:)));
        if ~isempty(D.c1), fprintf(', comp1 %s', vr(D.c1(:))); end
        if ~isempty(D.c2), fprintf(', comp2 %s', vr(D.c2(:))); end
        fprintf('\n');
        set(ax, 'YLimMode', 'auto', 'YTickMode', 'auto');
    end
    xlabel(ax, 'wavelength (nm)', 'FontName', opt.fontName, 'FontSize', opt.fontSize);
    if ~isempty(ylab)
        ylabel(ax, ylab, 'FontName', opt.fontName, 'FontSize', opt.fontSize);
    end
    % panel name and R^2 above the axes (explicit text objects)
    if sumMode
        r2s = 'in silico mixture';
    else
        r2s = sprintf('R^2 = %.4f %s %.4f', D.r2Mean, char(177), D.r2Sem);
    end
    text(ax, 0.5, 1.20, ttl, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontName', opt.fontName, 'FontSize', opt.fontSize, 'Clipping', 'off');
    text(ax, 0.5, 1.04, r2s, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontName', opt.fontName, 'FontSize', opt.fontSize, 'Clipping', 'off');
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


function figHeader(f, lines, opt)
% Header lines at the top of the figure, on an invisible axes.
    sz = get(f, 'UserData');
    ax = axes(f, 'Units', 'centimeters', 'Position', [0 sz(2) - 2.2 sz(1) 2.1], 'Visible', 'off', ...
        'XLim', [0 1], 'YLim', [0 1]);
    n = numel(lines);
    for k = 1:n
        y = 1 - (k - 0.5) / n;
        text(ax, 0.5, y, lines{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontName', opt.fontName, 'FontSize', opt.fontSize, 'Interpreter', 'tex');
    end
end


function s = tex(s)
    s = strrep(s, '\', '\\');
    s = strrep(s, '_', '\_');
    s = strrep(s, '^', '\^');
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


function v = tickVector(lim, step)
% Ticks at multiples of step within lim; strictly increasing.
    n = floor(lim(2) / step + 1e-9) - ceil(lim(1) / step - 1e-9);
    if ~isfinite(n) || n < 0 || n > 40
        % step unsuitable for this range: use a 1/2/2.5/5 x 10^k step instead
        raw = diff(lim) / 4;
        p10 = 10 ^ floor(log10(raw));
        m = [1 2 2.5 5 10] * p10;
        step = m(find(m >= raw, 1));
        if isempty(step), step = m(end); end
    end
    v = (ceil(lim(1) / step - 1e-9):floor(lim(2) / step + 1e-9)) * step;
    v(abs(v) < 1e-6 * step) = 0;   % avoid a "-0" label
    v = unique(v);
end


function file = savePdf(f, file)
    sizeCm = get(f, 'UserData');
    if exist(file, 'file') == 2
        try
            delete(file);
        catch
        end
        if exist(file, 'file') == 2
            [p, n, e] = fileparts(file);
            file = fullfile(p, [n '_' datestr(now, 'yyyymmdd_HHMMSS') e]);
            fprintf('Existing PDF could not be replaced (open?). Writing %s\n', file);
        end
    end
    try
        if exist('exportgraphics') > 0 %#ok<EXIST>
            exportgraphics(f, file, 'ContentType', 'vector', 'BackgroundColor', 'white');
        else
            set(f, 'PaperUnits', 'centimeters', 'PaperSize', sizeCm, ...
                'PaperPosition', [0 0 sizeCm]);
            print(f, file, '-dpdf', '-painters');
        end
        fprintf('   saved %s\n', file);
    catch ME
        fprintf('PDF not saved: %s\n', ME.message);
        file = '';
    end
    close(f);   % figures are not displayed; closed after saving
end


%% =====================================================================
%  Excel
%  =====================================================================
function matFile = saveResultsMat(hyb, names, types, cols, cfg, opt, modeTxt, unit, matFile)
% Results for NMF_Hybrid_Equivalence (figures and statistics without re-analysis).
% R.mixtures(i).cmp.m(q): key, label, kind ('auc', 'dauc', 'ratio', 'peak'), slot
% (1 green, 2 red), raw, nmf, rawTrials (trial-level raw AUC, AUC metrics only).
    R.version = 'v2';
    R.windowRule = 'own peak +/- 10 nm (ground truth: ground-truth peak, unmixed: unmixed peak)';
    R.mode = modeTxt;
    R.unit = unit;
    R.names = names;
    R.types = types;
    R.colors = cols;
    R.pairing = cfg.pairing;
    R.inputs = cfg;
    R.settings = opt;
    R.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    keep = {'index', 'label', 'green', 'red', 'greenFile', 'redFile', 'nTrialsGreen', 'nTrialsRed', ...
            'nTrials', 'nIter', 'prepOffset', 'leakagePct', 'cmp'};
    for i = 1:numel(hyb)
        for k = 1:numel(keep)
            h.(keep{k}) = hyb(i).(keep{k});
        end
        h.r2 = struct('metric', {hyb(i).panel.metric}, 'r2Mean', {hyb(i).panel.r2Mean}, ...
                      'r2Sem', {hyb(i).panel.r2Sem});
        if i == 1
            R.hybrids = h;
        else
            R.hybrids(i) = h;
        end
    end
    try
        save(matFile, 'R', '-v7');
        fprintf('Results file (.mat): %s\n', matFile);
    catch ME
        fprintf('.mat file not written: %s\n', ME.message);
        matFile = '';
    end
end


function xlsFile = writeExcel(hyb, Sm, names, types, cfg, opt, unit, xlsFile)
% Sheets: README, R2, R2_trials, AUC_peak, Summary, Info (plain-English column names).
    if exist(xlsFile, 'file') == 2
        try
            delete(xlsFile);
        catch
        end
        if exist(xlsFile, 'file') == 2
            [p, n, e] = fileparts(xlsFile);
            xlsFile = fullfile(p, [n '_' datestr(now, 'yyyymmdd_HHMMSS') e]);
            fprintf('Existing Excel file could not be replaced (open?). Writing %s\n', xlsFile);
        end
    end
    isDye = any(strcmp(types, 'dye'));
    P = hyb(1).panel;
    nP = numel(P);
    r2Lab = cell(1, nP);
    for q = 1:nP
        r2Lab{q} = r2Label(P(q), names, types);
    end
    M = hyb(1).cmp.m;
    srcHdr = cell(1, 2);
    for s = 1:2
        if strcmp(types{s}, 'dye')
            srcHdr{s} = [names{s} ' (dye)'];
        else
            srcHdr{s} = [names{s} ' animal'];
        end
    end
    idHdr = {'No', unit, srcHdr{1}, srcHdr{2}};
    idRow = @(H) {H.index, H.label, srcName(H.green, types{1}), srcName(H.red, types{2})};
    nU = sprintf('n %ss', lower(unit));
    try
        % ----- README
        C = { ...
            'NMF-RI validation with in silico mixtures, VERSION 2', 'AUC windows: ground truth and unmixed each around their own peak (v1: one window, ground-truth peak, for both)'; ...
            'Mode', ternary(isDye, sprintf('sensor + control dye (%s + %s)', srcLabel(names{1}, types{1}), srcLabel(names{2}, types{2})), ...
                                   sprintf('sensor + sensor (%s + %s)', names{1}, names{2})); ...
            sprintf('Number of %ss', lower(unit)), numel(hyb); ...
            'Date', datestr(now, 'yyyy-mm-dd HH:MM:SS'); ...
            '', ''; ...
            'Sheet', 'Content'; ...
            'R2', sprintf('One row per %s. R2 of each figure panel: mean over trials, then the SEM over trials (SEM = SD / sqrt(n trials)).', lower(unit)); ...
            'R2_trials', 'One row per trial: R2 of each figure panel.'; ...
            'AUC_peak', sprintf('One row per %s. Ground truth vs unmixed (NMF): AUC, T2/T1 AUC ratio, delta AUC and peak wavelength.', lower(unit)); ...
            'Summary', sprintf('Mean and SEM across %ss (SEM = SD / sqrt(n %ss)) of every measure.', lower(unit), lower(unit)); ...
            'Info', sprintf('One row per %s: trial numbers, NMF-RI iterations, preprocessing offset, normalisation values, files.', lower(unit)); ...
            '', ''; ...
            'Definition', ''; ...
            'Ground truth', 'measured single-sensor data minus the offset of each file (one value per file: minimum of the time-averaged spectrum, 450-650 nm); dye: measured mean dye spectrum minus its minimum'; ...
            'Unmixed', 'NMF-RI component: A(:,k) * H(k,:); component 1 = green slot, component 2 = red slot'; ...
            'In silico mixture', 'green ground truth + red ground truth, trial k + trial k, first min(n trials) trials; a dye spectrum is added to every time point'; ...
            'Reward spectrum', sprintf('mean over reward samples %d-%d of each trial', opt.rewardSamples(1), opt.rewardSamples(end)); ...
            'Delta spectrum', sprintf('reward spectrum - mean over baseline samples %d-%d of the same trial', opt.baselineSamples(1), opt.baselineSamples(end)); ...
            'R2', '1 - SSres / SStot over 450-650 nm, per trial, in counts; ground truth = observed, unmixed = predicted'; ...
            'R2 of sum', 'ground-truth green + ground-truth red vs unmixed component 1 + component 2'; ...
            'AUC (v2)', sprintf('reward spectrum, mean over trials; ground truth: ground-truth peak +/- %g nm; unmixed: unmixed peak +/- %g nm (separate windows)', opt.aucHalfWidth, opt.aucHalfWidth); ...
            'T1 AUC window', sprintf('green peak within %d-%d nm: ground-truth green peak for ground truth, unmixed green peak for unmixed', opt.bandT1(1), opt.bandT1(2)); ...
            'T2 AUC window', sprintf('red peak within %d-%d nm: ground-truth red peak for ground truth, unmixed red peak for unmixed (also used for the green T2 AUC)', opt.bandT2(1), opt.bandT2(2)); ...
            'T2/T1 AUC ratio', 'T2 AUC / T1 AUC of the green sensor'; ...
            'Delta AUC', sprintf('AUC of the delta spectrum (reward - baseline, mean over trials) in the same windows as the AUC (ground truth: ground-truth peak +/- %g nm; unmixed: unmixed peak +/- %g nm)', opt.aucHalfWidth, opt.aucHalfWidth); ...
            'Peak', sprintf('highest local maximum after a %g nm moving average, found separately in ground truth and unmixed', opt.peakSmoothNm); ...
            'Difference', 'AUC and ratio: (unmixed - ground truth) / ground truth x 100 (%); peak: unmixed - ground truth (nm)'; ...
            'Paired t-test', sprintf(['unmixed vs ground truth, one pair of values per %s, two-sided; t = mean(d) / (SD(d) / sqrt(n)), ' ...
                'd = unmixed - ground truth (in the measure''s own unit), df = n - 1. n/a: fewer than 2 values or identical differences.'], lower(unit)); ...
            'Dye delta leakage', '100 x AUC of the unmixed dye delta spectrum / AUC of the ground-truth sensor delta spectrum (trial means, signed, 450-650 nm)'; ...
            'Preprocessing offset', 'mean per-time-point minimum removed from the mixture by dataPreprocessing (not part of the NMF components)'; ...
            'Normalisation value', 'value the figure spectra are divided by: reward = ground-truth peak; delta (summary figure) = max |ground truth|'; ...
            '', ''; ...
            'Setting', 'Value'; ...
            'Green slot (component 1)', sprintf('%s: %s', srcLabel(names{1}, types{1}), cfg.greenPath); ...
            'Red slot (component 2)',   sprintf('%s: %s', srcLabel(names{2}, types{2}), cfg.redPath); ...
            'Component 1 CSV (initial A0)', cfg.comp1Csv; ...
            'Component 2 CSV (initial A0)', cfg.comp2Csv; ...
            'Wavelength range (nm)', sprintf('%g-%g', cfg.wlRange(1), cfg.wlRange(2)); ...
            'theta1', cfg.theta1; 'theta2', cfg.theta2; 'alphaA', cfg.alphaA; 'alphaX', cfg.alphaX; ...
            'Pairing (sensor + sensor only)', cfg.pairing; ...
            'Trial length (samples)', opt.trialLength};
        writecell(cleanCell(C), xlsFile, 'Sheet', 'README');

        % ----- R2 (mean, then SEM over trials)
        C = [idHdr, {'Trials used'}, r2Lab, strcat({'SEM of '}, r2Lab)];
        for i = 1:numel(hyb)
            H = hyb(i);
            C(end+1, :) = [idRow(H), {H.nTrials}, num2cell([H.panel.r2Mean]), num2cell([H.panel.r2Sem])]; %#ok<AGROW>
        end
        writecell(cleanCell(C), xlsFile, 'Sheet', 'R2');

        % ----- R2 per trial
        C = [idHdr, {'Trial'}, r2Lab];
        for i = 1:numel(hyb)
            H = hyb(i);
            for k = 1:H.nTrials
                row = [idRow(H), {k}];
                for q = 1:nP
                    row{end+1} = H.panel(q).r2Trials(k); %#ok<AGROW>
                end
                C(end+1, :) = row; %#ok<AGROW>
            end
        end
        writecell(cleanCell(C), xlsFile, 'Sheet', 'R2_trials');

        % ----- AUC and peak: raw vs unmixed
        if ~isempty(M)
            hdr = [idHdr, {'Trials used'}];
            if strcmp(types{1}, 'sensor')
                hdr{end+1} = sprintf('T1 AUC window centre for ground truth (nm) = ground-truth %s T1 peak', names{1});
                hdr{end+1} = sprintf('T1 AUC window centre for unmixed (nm) = unmixed %s T1 peak', names{1});
            end
            hdr{end+1} = sprintf('T2 AUC window centre for ground truth (nm) = ground-truth %s T2 peak', srcLabel(names{2}, types{2}));
            hdr{end+1} = sprintf('T2 AUC window centre for unmixed (nm) = unmixed %s T2 peak', srcLabel(names{2}, types{2}));
            for q = 1:numel(M)
                cn = cmpNames(M(q));
                hdr = [hdr, cn(1:3)]; %#ok<AGROW>
            end
            C = hdr;
            for i = 1:numel(hyb)
                H = hyb(i);
                row = [idRow(H), {H.nTrials}];
                if strcmp(types{1}, 'sensor'), row = [row, {H.cmp.winT1raw, H.cmp.winT1nmf}]; end %#ok<AGROW>
                row = [row, {H.cmp.winT2raw, H.cmp.winT2nmf}]; %#ok<AGROW>
                for q = 1:numel(M)
                    mq = H.cmp.m(q);
                    row = [row, {mq.raw, mq.nmf, cmpDiff(mq)}]; %#ok<AGROW>
                end
                C(end+1, :) = row; %#ok<AGROW>
            end
            writecell(cleanCell(C), xlsFile, 'Sheet', 'AUC_peak');
        end

        % ----- Summary across mixtures (two tables)
        W = 13;
        C = cell(0, W);
        C(end+1, 1) = {sprintf('R2 (per %s: mean over trials; here: mean and SEM across %ss)', lower(unit), lower(unit))};
        C(end+1, 1:4) = {'Measure', nU, 'Mean', 'SEM'};
        for q = 1:nP
            D = Sm.panel(q);
            C(end+1, 1:4) = {r2Lab{q}, D.r2N, D.r2Mean, D.r2Sem}; %#ok<AGROW>
        end
        if isDye
            C(end+1, 1:4) = {'Dye delta leakage (%)', Sm.leakN, Sm.leakMean, Sm.leakSem};
        end
        if ~isempty(M)
            C(end+2, 1) = {sprintf('Ground truth vs unmixed (reward spectra and delta spectra, mean over trials; here: mean and SEM across %ss)', lower(unit))};
            C(end+1, :) = {'Measure', 'Unit', nU, 'Ground-truth mean', 'Ground-truth SEM', 'Unmixed mean', 'Unmixed SEM', ...
                           'Difference mean', 'Difference SEM', 'Difference unit', ...
                           'Paired t-test t', 'Paired t-test df', 'Paired t-test p (two-sided)'};
            for q = 1:numel(M)
                cn = cmpNames(M(q));
                raw = arrayfun(@(H) H.cmp.m(q).raw, hyb);
                nmf = arrayfun(@(H) H.cmp.m(q).nmf, hyb);
                d = arrayfun(@(H) cmpDiff(H.cmp.m(q)), hyb);
                [mr, sr, nr] = meanSem(raw);
                [mn, sn] = meanSem(nmf);
                [md, sd] = meanSem(d);
                [tt, df, pp] = pairedT(raw(:), nmf(:));
                C(end+1, :) = {cn{4}, cn{5}, nr, mr, sr, mn, sn, md, sd, cn{6}, tt, df, pp}; %#ok<AGROW>
            end
        end
        if ~isDye && strcmp(cfg.pairing, 'all')
            C(end+2, 1) = {'All combinations: mixtures share mice and are not independent.'};
        end
        writecell(cleanCell(C), xlsFile, 'Sheet', 'Summary');

        % ----- Info
        hdr = [idHdr, {sprintf('%s trials available', srcLabel(names{1}, types{1})), ...
                       sprintf('%s trials available', srcLabel(names{2}, types{2})), ...
                       'Trials used', 'NMF-RI iterations', 'Preprocessing offset removed (counts)'}];
        if isDye
            hdr{end+1} = 'Dye delta leakage (%)';
        end
        for q = 1:nP
            hdr{end+1} = ['Normalisation value (counts): ' strrep(r2Lab{q}, 'R2 of ', '')]; %#ok<AGROW>
        end
        hdr = [hdr, {[srcLabel(names{1}, types{1}) ' file'], [srcLabel(names{2}, types{2}) ' file']}];
        C = hdr;
        for i = 1:numel(hyb)
            H = hyb(i);
            row = [idRow(H), {trialCount(H.nTrialsGreen), trialCount(H.nTrialsRed), H.nTrials, H.nIter, H.prepOffset}];
            if isDye
                row{end+1} = H.leakagePct; %#ok<AGROW>
            end
            row = [row, num2cell([H.panel.peak]), {H.greenFile, H.redFile}]; %#ok<AGROW>
            C(end+1, :) = row; %#ok<AGROW>
        end
        writecell(cleanCell(C), xlsFile, 'Sheet', 'Info');
        fprintf('\nExcel file: %s\n', xlsFile);
    catch ME
        fprintf('Excel file not written: %s\n', ME.message);
        xlsFile = '';
    end
end


function t = srcLabel(name, type)
    if strcmp(type, 'dye'), t = [name ' (dye)']; else, t = name; end
end


function lab = r2Label(Q, names, types)
% e.g. 'R2 of GCaMP unmixing, reward spectrum'
    if Q.row == 1, per = 'reward spectrum'; else, per = 'delta spectrum'; end
    if Q.col == 1
        what = sprintf('sum (%s + %s)', srcLabel(names{1}, types{1}), srcLabel(names{2}, types{2}));
    else
        s = Q.col - 1;
        what = sprintf('%s unmixing', srcLabel(names{s}, types{s}));
    end
    lab = sprintf('R2 of %s, %s', what, per);
end


function c = cmpNames(mq)
% {raw column, unmixed column, difference column, summary name, unit, difference unit}
    b = mq.label;                     % e.g. 'GCaMP T1 AUC', 'GCaMP T2/T1', 'GCaMP T1 peak'
    switch mq.kind
        case {'auc', 'dauc'}
            c = {['ground-truth ' b ' (counts x nm)'], ['unmixed ' b ' (counts x nm)'], ...
                 [b ' difference, (unmixed - ground truth) / ground truth (%)'], b, 'counts x nm', '%'};
        case 'ratio'
            b = [b ' AUC ratio'];
            c = {['ground-truth ' b], ['unmixed ' b], [b ' difference, (unmixed - ground truth) / ground truth (%)'], b, 'ratio', '%'};
        otherwise
            c = {['ground-truth ' b ' (nm)'], ['unmixed ' b ' (nm)'], [b ' difference, unmixed - ground truth (nm)'], b, 'nm', 'nm'};
    end
end


function d = cmpDiff(mq)
    if strcmp(mq.kind, 'peak')
        d = mq.nmf - mq.raw;
    else
        d = 100 * (mq.nmf - mq.raw) / mq.raw;
    end
end


function v = srcName(name, type)
    if strcmp(type, 'dye'), v = [name ' (dye)']; else, v = name; end
end


function v = trialCount(n)
    if isfinite(n), v = n; else, v = 'dye (static spectrum)'; end
end


function v = ternary(c, a, b)
    if c, v = a; else, v = b; end
end


function C = cleanCell(C)
    C(cellfun(@isempty, C)) = {''};
    for k = 1:numel(C)
        if isnumeric(C{k}) && isscalar(C{k}) && ~isfinite(C{k})
            C{k} = 'NaN';
        end
    end
end


%% =====================================================================
%  Spectra and files
%  =====================================================================
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


function X = onAxis(wlFrom, X, wlTo)
    if numel(wlFrom) == numel(wlTo) && max(abs(wlFrom(:) - wlTo(:))) < 1e-6
        return;
    end
    X = interp1(wlFrom(:), X, wlTo(:), 'linear', 0);
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
