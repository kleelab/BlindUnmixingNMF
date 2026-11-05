function varargout = NMF_Power_Decoder(cfg)
% NMF_Power_Decoder
% Builds one spectral "power decoder" per animal from static spectra measured at
% several excitation powers (e.g. 25, 50, 100 uW), for one sensor, and saves all
% animals of that sensor in one .mat file.
%
% Why: when the brightness of a recording changes, the spectrum does not simply scale
% (T2 does not rise in proportion to T1, peaks can shift). The decoder stores how every
% wavelength changes with the overall brightness of that animal, measured.
%
% Decoder of an animal (two-component model)
%   Spectra S_p(lambda) at the powers p (mean over the frames of the file, offset removed
%   as in the time-series analysis). Brightness ratio k_p = total AUC (wlRange) at p /
%   total AUC at the reference power (P.refPower, the power of the time-series recordings).
%   Every wavelength is fitted over all powers by a straight line in k:
%       S(lambda, k) = B(lambda) + k E(lambda)
%   B: part that does not change with the brightness (background, residual offset);
%   E: part that is proportional to it (sensor emission). Fixed fraction B / (B + E).
%   Gain G_k(lambda) = (B + k E) / (B + E): the factor by which every wavelength of a
%   reference-power recording is multiplied to make it k times brighter (applied later to
%   the time series: Y_k(lambda, t) = Y(lambda, t) x G_k(lambda)). Where the reference
%   spectrum is below P.minFrac of its maximum the gain is k (proportional).
%   Relative gain G_k / k: 1 = the wavelength scales in proportion to the total AUC.
%
% Decoding validation (cross-validation; the held-out spectrum is never used to predict it)
%   For every animal and every power p, the spectrum at p is predicted by
%     decoder (animal): line through the animal's other powers (leave-one-power-out)
%     decoder (pooled): fixed fraction averaged over the OTHER animals (leave-one-animal-out),
%                       applied to the animal's nearest other power
%     proportional:     the animal's nearest other power x (k_p / k_base)  (baseline)
%   Decoding accuracy = 1 - ||predicted - measured|| / ||measured|| (P.evalRange; 1 = perfect);
%   prediction error = normalized RMSE (%) = 100 x (1 - accuracy); also T1 and T2 AUC errors.
%   Noise ceiling: half of the difference between the means of the odd and even frames.
%   Mean +/- SEM over animals at every power; paired t-test and exact Wilcoxon signed-rank
%   test (decoder vs proportional, decoder vs pooled), n = animals.
%   Also: peak wavelengths (T1, T2), T2/T1 and linearity (R^2 of the line) at every power.
%
% Folder layout (one sensor):
%   <sensor folder>\power25\<animal>_<...>     e.g. GCa4_G25_Subt2__0__17-29-06-052
%   <sensor folder>\power50\...   <sensor folder>\power100\...
%   Animal name = text before the first '_' of the file name (P.namePattern); upper and
%   lower case are treated as the same (GCa4 = gca4 = GCA4).
%   Power = number in the folder name (power25, Power25, POWER25 -> 25).
%
% Usage
%   NMF_Power_Decoder                 % dialogs: sensor folder (e.g. ...\power\GCaMP), output folder
%   NMF_Power_Decoder(cfg)            % cfg.sensorDir, cfg.outDir, optional cfg.sensorName
%   NMF_Power_Decoder('...\PowerDecoder_GCaMP.mat')   % redraw figures and Excel
%   Run once per sensor (GCaMP, jRGECO).
%
% Output (cfg.outDir)
%   PowerDecoder_<sensor>.mat      variable D (all animals of the sensor)
%   PowerDecoder_<sensor>_animals.pdf/.eps   relative gain of every animal
%   PowerDecoder_<sensor>_summary.pdf/.eps   relative gain, T2/T1, peaks
%   PowerDecoder_<sensor>_decoding_validation.pdf/.eps   decoding accuracy (0-1) and prediction
%       error (normalized RMSE, %) vs power, mean +/- SEM over animals
%   PowerDecoder_<sensor>.xlsx     README, Animals, CrossValidation, CVSummary, CVTests, RelGain

%% ===================== SETTINGS =====================
P.refPower       = 25;            % uW, power of the time-series recordings (gain = 1 there)
P.wlRange        = [450 650];     % nm, as the NMF analysis; total AUC over this range
P.namePattern    = '^([^_]+)_';   % animal name = first token of the file name
P.offset         = 'fileMin';     % 'fileMin': minus the minimum of the mean spectrum (as the
                                  % time-series analysis); 'none': as recorded
P.evalRange      = [480 630];     % nm, range of the decoding accuracy / RMSE
P.exclude        = {};            % animals left out, e.g. {'FjR3'} (upper / lower case ignored)
P.qcMinStep      = 1.1;           % quality check: the total AUC must rise by >= this factor from one
                                  % power to the next higher one; otherwise the animal is excluded
P.minFrac        = 0.02;          % gain only where the reference spectrum >= 2% of its maximum
P.gainSmoothNm   = 2;             % moving average of the gain (nm)
P.bandT1         = [500 540];     % nm, fixed T1 band (AUC) and T1 peak search
P.bandT2         = [580 620];     % nm, fixed T2 band (AUC) and T2 peak search
P.peakSmoothNm   = 2;             % smoothing before the peak search (nm)
P.fontName       = 'Arial';
P.fontSize       = 7;
P.panelCm        = [4.0 3.4];     % panel width, height (cm)
P.lineWidth      = 0.8;
P.dotSize        = 8;
P.axisLineWidth  = 0.5;
%% ====================================================

    if nargin >= 1 && ischar(cfg) && numel(cfg) > 4 && strcmpi(cfg(end-3:end), '.mat')
        S = load(cfg);
        D = S.D;
        [d, n] = fileparts(cfg);
        base = fullfile(d, n);
        fprintf('Redrawing %s\n', cfg);
        decoderFigures(D, P, base);
        writeDecoderExcel(D, P, [base '.xlsx']);
        if nargout > 0, varargout{1} = D; end
        return;
    end
    if nargin < 1 || isempty(cfg)
        cfg = struct();
        dd = uigetdir(pwd, 'Sensor folder with power25, power50, power100 (e.g. ...\power\GCaMP)');
        if ~ischar(dd), fprintf('Cancelled.\n'); return; end
        cfg.sensorDir = dd;
        dd = uigetdir(fileparts(cfg.sensorDir), 'Output folder');
        if ~ischar(dd), fprintf('Cancelled.\n'); return; end
        cfg.outDir = dd;
    end
    cfg.sensorDir = stripSep(cfg.sensorDir);
    if ~isfield(cfg, 'outDir') || isempty(cfg.outDir), cfg.outDir = fileparts(cfg.sensorDir); end
    if ~isfield(cfg, 'sensorName') || isempty(cfg.sensorName), [~, cfg.sensorName] = fileparts(cfg.sensorDir); end
    fprintf('\nInput : %s\nOutput: %s\n', cfg.sensorDir, cfg.outDir);
    if exist(cfg.sensorDir, 'dir') ~= 7, error('Folder not found: %s', cfg.sensorDir); end
    if exist(cfg.outDir, 'dir') ~= 7, mkdir(cfg.outDir); end
    warn = {};
    t0 = tic;

    % ----- power folders
    dd = dir(cfg.sensorDir);
    dd = dd([dd.isdir]);
    pw = []; pdir = {};
    for i = 1:numel(dd)
        tk = regexp(lower(dd(i).name), 'power\D*(\d+(\.\d+)?)', 'tokens', 'once');
        if isempty(tk), continue; end
        pw(end+1) = str2double(tk{1}); %#ok<AGROW>
        pdir{end+1} = fullfile(cfg.sensorDir, dd(i).name); %#ok<AGROW>
    end
    if isempty(pw), error('No power folders (power25, power50, ...) in %s', cfg.sensorDir); end
    [pw, o] = sort(pw); pdir = pdir(o);
    if ~any(abs(pw - P.refPower) < 1e-9)
        error('No folder for the reference power %g uW in %s', P.refPower, cfg.sensorDir);
    end
    fprintf('\n==============================================================\n');
    fprintf(' Power decoder: %s   powers %s uW (reference %g uW)\n', cfg.sensorName, num2str(pw), P.refPower);
    fprintf('==============================================================\n');

    % ----- read every file: mean spectrum
    rec = struct('animal', {}, 'power', {}, 'file', {}, 'wl', {}, 'spec', {}, 'hA', {}, 'hB', {}, 'nFrames', {}, 'offset', {});
    for ip = 1:numel(pw)
        ff = dir(pdir{ip});
        ff = ff(~[ff.isdir]);
        for i = 1:numel(ff)
            [~, nm, ext] = fileparts(ff(i).name);
            if any(strcmpi(ext, {'.mat', '.xlsx', '.xls', '.pdf', '.eps', '.png', '.jpg', '.fig', '.m', '.ini'})), continue; end
            tk = regexp(ff(i).name, P.namePattern, 'tokens', 'once');
            if isempty(tk), warn{end+1} = sprintf('%s: no animal name (P.namePattern); skipped.', ff(i).name); continue; end %#ok<AGROW>
            file = fullfile(pdir{ip}, ff(i).name);
            try
                [wl, Y, ~, fmt] = readDyeTxt(file, P.wlRange);
            catch ME
                warn{end+1} = sprintf('%s: not read (%s); skipped.', ff(i).name, ME.message); %#ok<AGROW>
                continue;
            end
            if isempty(wl) || isempty(Y)
                warn{end+1} = sprintf('%s: no values within %g-%g nm; skipped.', ff(i).name, P.wlRange); %#ok<AGROW>
                continue;
            end
            Y = double(Y);
            s = mean(Y, 2);
            hA = nan(size(s)); hB = nan(size(s));
            if size(Y, 2) >= 2, hA = mean(Y(:, 1:2:end), 2); hB = mean(Y(:, 2:2:end), 2); end
            off = 0;
            if strcmpi(P.offset, 'fileMin'), off = min(s); s = s - off; hA = hA - off; hB = hB - off; end
            rec(end+1) = struct('animal', tk{1}, 'power', pw(ip), 'file', [nm ext], 'wl', wl(:), 'spec', s(:), ...
                'hA', hA(:), 'hB', hB(:), 'nFrames', size(Y, 2), 'offset', off); %#ok<AGROW>
            fprintf(' %6g uW  %-10s %-45s %5d frame(s) (%s), offset %.1f\n', pw(ip), tk{1}, [nm ext], size(Y, 2), fmt, off);
        end
    end
    if isempty(rec), error('No spectra read.'); end

    % ----- common wavelength axis: the first reference-power file
    iRef = find([rec.power] == P.refPower, 1);
    wl = rec(iRef).wl;
    keys = lower({rec.animal});                      % animal names: upper / lower case treated as the same
    [ukeys, iFirst] = unique(keys);
    names = {rec(iFirst).animal};
    A = [];
    for a = 1:numel(names)
        rr = rec(strcmp(keys, ukeys{a}));
        ps = unique([rr.power]);
        if ~any(ps == P.refPower)
            warn{end+1} = sprintf('%s: no %g uW spectrum; animal skipped.', names{a}, P.refPower); %#ok<AGROW>
            continue;
        end
        if numel(ps) < 2
            warn{end+1} = sprintf('%s: only one power; animal skipped.', names{a}); %#ok<AGROW>
            continue;
        end
        Sx = nan(numel(wl), numel(ps)); HA = Sx; HB = Sx; files = cell(1, numel(ps)); nfr = zeros(1, numel(ps)); offs = zeros(1, numel(ps));
        for j = 1:numel(ps)
            q = rr([rr.power] == ps(j));
            if numel(q) > 1
                warn{end+1} = sprintf('%s %g uW: %d files, averaged.', names{a}, ps(j), numel(q)); %#ok<AGROW>
            end
            M = zeros(numel(wl), numel(q)); MA = M; MB = M;
            for i = 1:numel(q)
                s = q(i).spec; sa = q(i).hA; sb = q(i).hB;
                if numel(q(i).wl) ~= numel(wl) || max(abs(q(i).wl - wl)) > 1e-6
                    s = interp1(q(i).wl, s, wl, 'linear', 0);
                    sa = interp1(q(i).wl, sa, wl, 'linear', 0); sb = interp1(q(i).wl, sb, wl, 'linear', 0);
                end
                M(:, i) = s; MA(:, i) = sa; MB(:, i) = sb;
            end
            Sx(:, j) = mean(M, 2); HA(:, j) = mean(MA, 2); HB(:, j) = mean(MB, 2);
            files{j} = strjoin({q.file}, ' + ');
            nfr(j) = sum([q.nFrames]); offs(j) = mean([q.offset]);
        end
        if any(strcmpi(P.exclude, names{a}))
            warn{end+1} = sprintf('%s: excluded (P.exclude).', names{a}); %#ok<AGROW>
            continue;
        end
        tot = trapz(wl, Sx); rq = tot / tot(ps == P.refPower);
        if any(rq(2:end) ./ rq(1:end-1) < P.qcMinStep)
            warn{end+1} = sprintf(['%s: EXCLUDED by the quality check - total AUC does not rise with the power ' ...
                '(ratios %s at %s uW); check the files.'], names{a}, num2str(rq, '%.2f '), num2str(ps)); %#ok<AGROW>
            fprintf(' %s EXCLUDED: total AUC ratios %s at %s uW\n', names{a}, num2str(rq, '%.2f '), num2str(ps));
            continue;
        end
        A = [A, buildAnimal(names{a}, wl, ps, Sx, HA, HB, files, nfr, offs, P)]; %#ok<AGROW>
    end
    if isempty(A), error('No animal with the reference power and at least one other power.'); end
    A = crossValidate(A, wl, P);

    D.version = 'NMF_Power_Decoder 2.0';
    D.sensor = cfg.sensorName;
    D.sensorDir = cfg.sensorDir;
    D.wl = wl;
    D.powers = pw;
    D.refPower = P.refPower;
    D.settings = P;
    D.method = ['Two-component model of every wavelength over the powers: S = B + k E (k = total AUC / total ' ...
        'AUC at the reference power; B does not change with the brightness, E is proportional to it). ' ...
        'Spectrum at k = B + k E; gain = (B + k E) / (B + E).'];
    D.animals = A;
    D.cv = cvTable(A, pw);
    D.cvTests = cvTests(D.cv, pw);
    D.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    base = fullfile(cfg.outDir, safeFileName(sprintf('PowerDecoder_%s', cfg.sensorName)));
    save([base '.mat'], 'D');
    fprintf('   saved %s.mat\n', base);
    decoderFigures(D, P, base);
    writeDecoderExcel(D, P, [base '.xlsx']);
    printDecoderSummary(D);
    if ~isempty(warn)
        fprintf('\nWarnings (%d):\n', numel(warn)); fprintf('   %s\n', warn{:});
    end
    fprintf('Total time: %.1f s\n', toc(t0));
    if nargout > 0, varargout{1} = D; end
end


function printDecoderSummary(D)
    fprintf('\n%s (reference %g uW):\n', D.sensor, D.refPower);
    for a = 1:numel(D.animals)
        A = D.animals(a);
        fprintf('   %-10s powers %-14s ratio %-20s fixed part: T1 %5.1f %%, T2 %5.1f %%\n', A.name, num2str(A.powers), ...
            num2str(A.ratio, '%.2f '), 100 * A.fixedT1, 100 * A.fixedT2);
    end
    fprintf('\nCross-validation, decoding accuracy (mean +/- SEM over animals):\n');
    for p = D.powers(:)'
        fprintf('   %6g uW:', p);
        for m = {'animal', 'pooled', 'proportional', 'noise'}
            T = D.cv([D.cv.power] == p & strcmp({D.cv.method}, m{1}) & strcmp({D.cv.measure}, 'accuracy'));
            if ~isempty(T) && T.n > 0, fprintf('  %s %.4f +/- %.4f', m{1}, T.mean, T.sem); end
        end
        fprintf('\n');
    end
end


%% =====================================================================
%  Decoder of one animal
%  =====================================================================
function A = buildAnimal(name, wl, ps, S, HA, HB, files, nfr, offs, P)
% S: mean spectra (wl x powers); HA, HB: means of the odd and even frames (noise check)
    iR = find(ps == P.refPower, 1);
    tot = trapz(wl, S);
    ratio = tot / tot(iR);
    np = numel(ps);
    aT1 = zeros(1, np); aT2 = zeros(1, np); pk1 = zeros(1, np); pk2 = zeros(1, np);
    for j = 1:np
        aT1(j) = bandAuc(wl, S(:, j), P.bandT1);
        aT2(j) = bandAuc(wl, S(:, j), P.bandT2);
        sm = smoothNm(wl, S(:, j), P.peakSmoothNm);
        pk1(j) = bandPeak(wl, sm, P.bandT1);
        pk2(j) = bandPeak(wl, sm, P.bandT2);
    end
    A.name = name; A.files = files; A.nFrames = nfr; A.offsets = offs;
    A.powers = ps; A.spectra = S; A.halfA = HA; A.halfB = HB; A.ratio = ratio; A.totalAUC = tot;
    A.aucT1 = aT1; A.aucT2 = aT2; A.t2t1 = aT2 ./ aT1; A.peakT1 = pk1; A.peakT2 = pk2;
    % two-component model over all powers: S(lambda) = B(lambda) + ratio * E(lambda)
    %   B: part that does not change with the brightness, E: part proportional to it
    X = [ones(np, 1), ratio(:)];
    C = X \ S';                                  % 2 x wl
    A.fixed = C(1, :)';                          % B
    A.prop = C(2, :)';                           % E
    Sref = A.fixed + A.prop;                     % model at ratio 1
    A.fixedFrac = A.fixed ./ Sref;               % fraction of the reference that does not scale
    A.fixedFrac(Sref < P.minFrac * max(Sref)) = NaN;
    m1 = wl >= P.bandT1(1) & wl <= P.bandT1(2); m2 = wl >= P.bandT2(1) & wl <= P.bandT2(2);
    A.fixedT1 = bandAuc(wl, A.fixed, P.bandT1) / bandAuc(wl, Sref, P.bandT1);
    A.fixedT2 = bandAuc(wl, A.fixed, P.bandT2) / bandAuc(wl, Sref, P.bandT2);
    % linearity of every wavelength vs ratio (R^2 of the two-component line)
    A.linR2T1 = NaN; A.linR2T2 = NaN;
    if np >= 3
        res = S' - X * C;
        r2 = 1 - sum(res .^ 2, 1) ./ max(sum((S' - mean(S', 1)) .^ 2, 1), eps);
        A.linR2T1 = median(r2(m1)); A.linR2T2 = median(r2(m2));
    end
    % relative gain at the measured ratios (model): G / ratio, 1 = proportional
    A.relGain = nan(numel(wl), np);
    for j = 1:np
        g = decoderGain(A, ratio(j), wl, P);
        rg = g / ratio(j);
        rg(Sref < P.minFrac * max(Sref)) = NaN;
        A.relGain(:, j) = rg;
    end
    A.cv = [];                                   % filled by crossValidate (needs all animals)
end


function Sk = decoderSpectrum(A, k)
% Spectrum of animal A at brightness ratio k (two-component model; ratio 1 = reference power).
    Sk = max(0, A.fixed + k * A.prop);
end


function G = decoderGain(A, k, wl, P)
% Gain of animal A at brightness ratio k: multiply a reference-power spectrum (or every time
% point of a reference-power recording) by G to make it k times brighter. Where the reference
% is below P.minFrac of its maximum the gain is k (proportional).
    Sref = A.fixed + A.prop;
    ok = Sref >= P.minFrac * max(Sref);
    G = k * ones(size(Sref));
    G(ok) = (A.fixed(ok) + k * A.prop(ok)) ./ Sref(ok);
    G = smoothNm(wl, G, P.gainSmoothNm);
end


%% =====================================================================
%  Cross-validation
%  =====================================================================
function A = crossValidate(A, wl, P)
% For every animal and every power p (held out), the spectrum at p is predicted
%   decoder (animal):  two-component line through the animal's other powers
%   decoder (pooled):  fraction of the fixed part from the OTHER animals (leave-one-animal-out),
%                      applied to the animal's nearest other power
%   proportional:      the animal's nearest other power x (ratio p / ratio base)
% and compared with the measured spectrum. The brightness ratio of p is taken as measured.
% Noise ceiling: difference between the means of the odd and even frames of p.
% Accuracy = 1 - ||predicted - measured|| / ||measured|| over P.evalRange (1 = perfect).
    ev = wl >= P.evalRange(1) & wl <= P.evalRange(2);
    nA = numel(A);
    meth = {'animal', 'pooled', 'proportional'};
    for a = 1:nA
        ps = A(a).powers; S = A(a).spectra; r = A(a).ratio; np = numel(ps);
        % pooled fixed fraction from the other animals
        F = [];
        for b = setdiff(1:nA, a)
            if numel(A(b).powers) >= 2, F = [F, A(b).fixedFrac]; end %#ok<AGROW>
        end
        fPool = nan(numel(wl), 1);
        for i = 1:numel(wl)
            v = F(i, isfinite(F(i, :)));
            if ~isempty(v), fPool(i) = mean(v); end
        end
        fPool(~isfinite(fPool)) = 0;
        cv = struct('power', {}, 'ratio', {}, 'method', {}, 'accuracy', {}, 'nrmse', {}, 't1', {}, 't2', {}, 't2t1', {});
        for t = 1:np
            others = setdiff(1:np, t);
            if isempty(others), continue; end
            Sm = S(:, t); kt = r(t);
            [~, ib] = min(abs(log(r(others)) - log(kt))); base = others(ib);
            kb = r(base);
            pred = cell(1, 3);
            if numel(others) >= 2
                X = [ones(numel(others), 1), r(others)'];
                Cf = X \ S(:, others)';
                pred{1} = (Cf(1, :) + Cf(2, :) * kt)';
            else
                pred{1} = nan(size(Sm));
            end
            pred{2} = S(:, base) .* (fPool + (1 - fPool) * kt) ./ (fPool + (1 - fPool) * kb);
            pred{3} = S(:, base) * kt / kb;
            for m = 1:3
                cv(end+1) = cvEntry(wl, ev, pred{m}, Sm, ps(t), kt, meth{m}, P); %#ok<AGROW>
            end
            % noise ceiling: half-frame difference
            if all(isfinite(A(a).halfA(:, t))) && all(isfinite(A(a).halfB(:, t)))
                d = (A(a).halfA(:, t) - A(a).halfB(:, t)) / 2;
                cv(end+1) = cvEntry(wl, ev, Sm + d, Sm, ps(t), kt, 'noise', P); %#ok<AGROW>
            end
        end
        A(a).cv = cv;
    end
end


function e = cvEntry(wl, ev, Sp, Sm, p, k, meth, P)
    e.power = p; e.ratio = k; e.method = meth;
    if any(~isfinite(Sp))
        e.accuracy = NaN; e.nrmse = NaN; e.t1 = NaN; e.t2 = NaN; e.t2t1 = NaN; return;
    end
    nr = norm(Sp(ev) - Sm(ev)) / norm(Sm(ev));
    e.accuracy = max(0, 1 - nr);
    e.nrmse = 100 * nr;
    t1p = bandAuc(wl, Sp, P.bandT1); t1m = bandAuc(wl, Sm, P.bandT1);
    t2p = bandAuc(wl, Sp, P.bandT2); t2m = bandAuc(wl, Sm, P.bandT2);
    e.t1 = 100 * (t1p - t1m) / t1m;
    e.t2 = 100 * (t2p - t2m) / t2m;
    e.t2t1 = 100 * ((t2p / t1p) - (t2m / t1m)) / (t2m / t1m);
end


function T = cvTable(A, powers)
% mean, SEM and paired tests over animals for every power and method
    meth = {'animal', 'pooled', 'proportional', 'noise'};
    meas = {'accuracy', 'nrmse', 't1', 't2'};
    T = struct('power', {}, 'method', {}, 'measure', {}, 'n', {}, 'mean', {}, 'sem', {}, 'values', {}, 'animals', {});
    for p = powers(:)'
        for m = 1:numel(meth)
            for q = 1:numel(meas)
                v = []; nm = {};
                for a = 1:numel(A)
                    c = A(a).cv;
                    if isempty(c), continue; end
                    k = find([c.power] == p & strcmp({c.method}, meth{m}), 1);
                    if isempty(k), continue; end
                    x = c(k).(meas{q});
                    if q >= 3, x = abs(x); end              % AUC errors: absolute
                    if isfinite(x), v(end+1) = x; nm{end+1} = A(a).name; end %#ok<AGROW>
                end
                s = NaN; if numel(v) > 1, s = std(v) / sqrt(numel(v)); end
                mu = NaN; if ~isempty(v), mu = mean(v); end
                T(end+1) = struct('power', p, 'method', meth{m}, 'measure', meas{q}, 'n', numel(v), 'mean', mu, ...
                    'sem', s, 'values', v, 'animals', {nm}); %#ok<AGROW>
            end
        end
    end
end


function St = cvTests(T, powers)
% decoder (animal) vs proportional and vs pooled, paired over animals, at every power
    St = struct('power', {}, 'measure', {}, 'comparison', {}, 'n', {}, 'mean1', {}, 'mean2', {}, 't', {}, 'df', {}, ...
        'p', {}, 'pWilcoxon', {});
    meas = {'accuracy', 'nrmse', 't1', 't2'};
    cmp = {'animal', 'proportional'; 'animal', 'pooled'; 'pooled', 'proportional'};
    for p = powers(:)'
        for q = 1:numel(meas)
            for c = 1:size(cmp, 1)
                a = T([T.power] == p & strcmp({T.method}, cmp{c, 1}) & strcmp({T.measure}, meas{q}));
                b = T([T.power] == p & strcmp({T.method}, cmp{c, 2}) & strcmp({T.measure}, meas{q}));
                if isempty(a) || isempty(b), continue; end
                [nm, ia, ib] = intersect(a.animals, b.animals);
                x = a.values(ia); y = b.values(ib);
                [t, df, pv] = pairedT(y(:), x(:));
                pw = signRankExact(x(:) - y(:));
                St(end+1) = struct('power', p, 'measure', meas{q}, 'comparison', sprintf('%s vs %s', cmp{c, 1}, cmp{c, 2}), ...
                    'n', numel(nm), 'mean1', mean(x), 'mean2', mean(y), 't', t, 'df', df, 'p', pv, 'pWilcoxon', pw); %#ok<AGROW>
            end
        end
    end
end


function a = bandAuc(wl, y, band)
    m = wl >= band(1) & wl <= band(2);
    a = trapz(wl(m), y(m));
end


%% =====================================================================
%  Figures
%  =====================================================================
function decoderFigures(D, P, base)
    col = sensorCol(D.sensor);
    A = D.animals; nA = numel(A); wl = D.wl;
    vis = wl >= 470 & wl <= 630;                     % range used for automatic y limits
    pAll = D.powers(D.powers ~= D.refPower);
    pc = cell(1, numel(pAll));
    for j = 1:numel(pAll)
        t = 0.55 * (numel(pAll) - j) / max(numel(pAll) - 1, 1);   % highest power darkest
        pc{j} = col * (1 - t) + t;
    end
    % ---------- 1. every animal: relative gain at the measured powers
    nc = min(4, nA); nr = ceil(nA / nc);
    pw = P.panelCm(1); ph = P.panelCm(2); gx = 1.5; gy = 1.9; x0 = 1.6; y0 = 1.3;
    W = x0 + nc * pw + (nc - 1) * gx + 0.5; H = y0 + nr * ph + (nr - 1) * gy + 1.3;
    f = newFigure(W, H);
    v = [];
    for a = 1:nA
        m = A(a).powers ~= D.refPower;
        v = [v; reshape(A(a).relGain(vis, m), [], 1)]; %#ok<AGROW>
    end
    [yl, st] = niceAxis([v(isfinite(v)); 1]);
    for a = 1:nA
        r = ceil(a / nc); c = a - (r - 1) * nc;
        ax = newAxes(f, [x0 + (c - 1) * (pw + gx), y0 + (nr - r) * (ph + gy), pw, ph], P);
        bandLines(ax, yl, P);
        plot(ax, wl([1 end]), [1 1], '--', 'Color', [0 0 0], 'LineWidth', 0.5);
        lab = {};
        for j = 1:numel(pAll)
            i = find(A(a).powers == pAll(j), 1);
            if isempty(i), continue; end
            plot(ax, wl, A(a).relGain(:, i), '-', 'Color', pc{j}, 'LineWidth', P.lineWidth);
            lab{end+1} = sprintf('\\color[rgb]{%.2f %.2f %.2f}%g \\muW (\\times%.2f)', pc{j}, pAll(j), A(a).ratio(i)); %#ok<AGROW>
        end
        set(ax, 'XLim', P.wlRange, 'XTick', P.wlRange(1):50:P.wlRange(2), 'YLim', yl, 'YTick', yl(1):st:yl(2));
        text(ax, 0, 1.04, A(a).name, 'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontName', P.fontName, ...
            'FontSize', P.fontSize, 'FontWeight', 'bold', 'Interpreter', 'none');
        text(ax, 0.04, 0.97, lab, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontName', P.fontName, 'FontSize', P.fontSize - 1);
        if c == 1, ylabel(ax, 'relative gain (G / ratio)'); end
        if r == nr, xlabel(ax, 'wavelength (nm)'); end
    end
    ov = axes(f, 'Units', 'centimeters', 'Position', [0 0 W H], 'XLim', [0 W], 'YLim', [0 H], 'Visible', 'off');
    text(ov, 0.3, H - 0.3, sprintf(['%s: change of every wavelength relative to the total AUC (1 = proportional), ' ...
        'vs %g \\muW; dotted: T1 and T2 bands'], D.sensor, D.refPower), 'FontName', P.fontName, 'FontSize', P.fontSize, ...
        'VerticalAlignment', 'top', 'Interpreter', 'tex');
    saveFigure(f, [base '_animals']);

    % ---------- 2. summary
    cs = cell(1, nA);
    for a = 1:nA, t = 0.6 * (a - 1) / max(nA - 1, 1); cs{a} = col * (1 - t) + t; end
    pw2 = 4.0; ph2 = 3.6; gx2 = 1.6; x0 = 1.6; y0 = 1.3;
    W = x0 + 4 * pw2 + 3 * gx2 + 2.4; H = y0 + ph2 + 1.4;
    f = newFigure(W, H);
    % (a) relative gain at the highest power of every animal
    ax = newAxes(f, [x0, y0, pw2, ph2], P);
    G = nan(numel(wl), nA);
    for a = 1:nA
        [~, i] = max(A(a).ratio);
        G(:, a) = A(a).relGain(:, i);
    end
    gv = reshape(G(vis, :), [], 1);
    [yl, st] = niceAxis([gv(isfinite(gv)); 1]);
    bandLines(ax, yl, P);
    plot(ax, wl([1 end]), [1 1], '--', 'Color', [0 0 0], 'LineWidth', 0.5);
    for a = 1:nA, plot(ax, wl, G(:, a), '-', 'Color', cs{a}, 'LineWidth', 0.5); end
    Gm = nan(numel(wl), 1);
    for i = 1:numel(wl), v = G(i, isfinite(G(i, :))); if ~isempty(v), Gm(i) = mean(v); end, end
    plot(ax, wl, Gm, '-', 'Color', [0 0 0], 'LineWidth', 1.2);
    set(ax, 'XLim', P.wlRange, 'XTick', P.wlRange(1):50:P.wlRange(2), 'YLim', yl, 'YTick', yl(1):st:yl(2));
    xlabel(ax, 'wavelength (nm)'); ylabel(ax, 'relative gain (G / ratio)');
    ttl(ax, 'Relative gain, highest power', P);
    % (b) T2/T1 relative to the reference power
    ax = newAxes(f, [x0 + (pw2 + gx2), y0, pw2, ph2], P);
    xs = []; ys = [];
    for a = 1:nA
        [r, o] = sort(A(a).ratio); y = A(a).t2t1(o) / A(a).t2t1(A(a).powers(o) == D.refPower);
        drawLineDots(ax, r, y, cs{a}, P); xs = [xs, r]; ys = [ys, y]; %#ok<AGROW>
    end
    [xl, sx] = niceAxis([xs(:); 1]); xl(1) = min(xl(1), 0); [yl, st] = niceAxis([ys(:); 1]);
    plot(ax, xl, [1 1], '--', 'Color', [0 0 0], 'LineWidth', 0.5);
    set(ax, 'XLim', xl, 'XTick', xl(1):sx:xl(2), 'YLim', yl, 'YTick', yl(1):st:yl(2));
    xlabel(ax, 'brightness ratio (total AUC)'); ylabel(ax, 'T2/T1 (rel. to reference)');
    ttl(ax, 'T2/T1 vs brightness', P);
    % (c), (d) peak shifts
    for b = 1:2
        ax = newAxes(f, [x0 + (1 + b) * (pw2 + gx2), y0, pw2, ph2], P);
        ys = [];
        for a = 1:nA
            [r, o] = sort(A(a).ratio);
            if b == 1, pk = A(a).peakT1; else, pk = A(a).peakT2; end
            pk = pk(o); y = pk - pk(A(a).powers(o) == D.refPower);
            drawLineDots(ax, r, y, cs{a}, P); ys = [ys, y]; %#ok<AGROW>
        end
        [yl, st] = niceAxis([ys(:); -1; 1]);
        plot(ax, xl, [0 0], '--', 'Color', [0 0 0], 'LineWidth', 0.5);
        set(ax, 'XLim', xl, 'XTick', xl(1):sx:xl(2), 'YLim', yl, 'YTick', yl(1):st:yl(2));
        xlabel(ax, 'brightness ratio (total AUC)'); ylabel(ax, 'peak shift (nm)');
        if b == 1, ttl(ax, 'T1 peak shift', P); else, ttl(ax, 'T2 peak shift', P); end
    end
    % animals
    ov = axes(f, 'Units', 'centimeters', 'Position', [0 0 W H], 'XLim', [0 W], 'YLim', [0 H], 'Visible', 'off');
    hold(ov, 'on');
    lx = x0 + 4 * pw2 + 3 * gx2 + 0.5;
    for a = 1:nA
        yy = y0 + ph2 - 0.15 - (a - 1) * 0.38;
        plot(ov, lx, yy, 'o', 'MarkerSize', P.dotSize * 0.8, 'MarkerFaceColor', cs{a}, 'MarkerEdgeColor', 'none');
        text(ov, lx + 0.25, yy, A(a).name, 'FontName', P.fontName, 'FontSize', P.fontSize - 1, 'Interpreter', 'none');
    end
    text(ov, 0.3, H - 0.3, sprintf('%s power decoder (reference %g \\muW)', D.sensor, D.refPower), ...
        'FontName', P.fontName, 'FontSize', P.fontSize, 'VerticalAlignment', 'top');
    saveFigure(f, [base '_summary']);

    % ---------- 3. decoding validation
    validationFigure(D, P, [base '_decoding_validation']);
end


function validationFigure(D, P, fileBase)
% (a) decoding accuracy (0-1) and (b) prediction error (normalized RMSE, %) vs the held-out
% power, mean +/- SEM over animals: decoder (animal) in the sensor colour. Circles P.dotSize,
% no outline. (Other methods: Excel only.)
    col = sensorCol(D.sensor);
    meth = {'animal'};
    lab = {'decoder (animal)'};
    cm = {col};
    pws = D.powers(:)'; xp = 1:numel(pws);
    pw = 4.6; ph = 4.0; gx = 2.0; x0 = 1.7; y0 = 1.3;
    W = x0 + 2 * pw + gx + 0.6; H = y0 + ph + 1.8;
    f = newFigure(W, H);
    meas = {'accuracy', 'nrmse'};
    yl = {'decoding accuracy', 'prediction error (normalized RMSE, %)'};
    tt = {'Decoding accuracy', 'Prediction error'};
    cap = 0.08;
    for q = 1:2
        ax = newAxes(f, [x0 + (q - 1) * (pw + gx), y0, pw, ph], P);
        allv = 0;
        for m = numel(meth):-1:1                  % control first, decoder on top
            mu = nan(size(pws)); se = mu;
            for i = 1:numel(pws)
                T = D.cv([D.cv.power] == pws(i) & strcmp({D.cv.method}, meth{m}) & strcmp({D.cv.measure}, meas{q}));
                if ~isempty(T) && T.n > 0, mu(i) = T.mean; se(i) = T.sem; end
            end
            ok = isfinite(mu);
            if ~any(ok), continue; end
            plot(ax, xp(ok), mu(ok), '-', 'Color', cm{m}, 'LineWidth', 1.0);
            for i = find(ok & isfinite(se))
                plot(ax, xp(i) * [1 1], mu(i) + se(i) * [-1 1], '-', 'Color', cm{m}, 'LineWidth', 0.8);
                plot(ax, xp(i) + cap * [-1 1], (mu(i) - se(i)) * [1 1], '-', 'Color', cm{m}, 'LineWidth', 0.8);
                plot(ax, xp(i) + cap * [-1 1], (mu(i) + se(i)) * [1 1], '-', 'Color', cm{m}, 'LineWidth', 0.8);
            end
            plot(ax, xp(ok), mu(ok), 'o', 'MarkerSize', P.dotSize, 'MarkerFaceColor', cm{m}, ...
                'MarkerEdgeColor', 'none', 'LineStyle', 'none');
            allv = [allv, mu(ok) + max(se(ok), 0)]; %#ok<AGROW>
        end
        if q == 1
            set(ax, 'YLim', [0 1], 'YTick', 0:0.2:1);
        else
            [ylm, st] = niceAxis(allv(isfinite(allv))');
            ylm(1) = 0;
            set(ax, 'YLim', ylm, 'YTick', 0:st:ylm(2));
        end
        set(ax, 'XLim', [0.5 numel(pws) + 0.5], 'XTick', xp, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), pws, ...
            'UniformOutput', false));
        xlabel(ax, 'held-out power (\muW)'); ylabel(ax, yl{q});
        ttl(ax, tt{q}, P);
    end
    ov = axes(f, 'Units', 'centimeters', 'Position', [0 0 W H], 'XLim', [0 W], 'YLim', [0 H], 'Visible', 'off');
    hold(ov, 'on');
    nA = numel(D.animals);
    text(ov, 0.3, H - 0.3, {sprintf('%s decoder: spectrum at each power predicted without it (cross-validation)', D.sensor), ...
        sprintf('mean \\pm SEM, n = %d animals; accuracy = 1 - ||pred. - meas.|| / ||meas.|| (%g-%g nm)', nA, P.evalRange)}, 'FontName', P.fontName, 'FontSize', P.fontSize - 1, 'VerticalAlignment', 'top');
    saveFigure(f, fileBase);
end


function bandLines(ax, yl, P)
    for b = [P.bandT1, P.bandT2]
        plot(ax, [b b], yl, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    end
end


function drawLineDots(ax, x, y, c, P)
    plot(ax, x, y, '-', 'Color', c, 'LineWidth', P.lineWidth);
    plot(ax, x, y, 'o', 'MarkerSize', P.dotSize * 0.75, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'none');
end


function ttl(ax, s, P)
    text(ax, 0, 1.04, s, 'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontName', P.fontName, ...
        'FontSize', P.fontSize, 'FontWeight', 'bold', 'Interpreter', 'none');
end


function c = sensorCol(name)
    n = lower(name);
    if ~isempty(strfind(n, 'gcamp')) || ~isempty(strfind(n, 'gfp')) || ~isempty(strfind(n, 'green')) %#ok<STREMP>
        c = [0.13 0.55 0.25];
    else
        c = [0.80 0.15 0.15];
    end
end


%% =====================================================================
%  Excel
%  =====================================================================
function writeDecoderExcel(D, P, xlsFile)
    if exist(xlsFile, 'file') == 2
        try, delete(xlsFile); catch, end %#ok<NOCOM>
    end
    try
        C = {'Power decoder (static spectra at several excitation powers)', ''; ...
             'Version', D.version; 'Date', D.date; 'Sensor', D.sensor; 'Folder', D.sensorDir; ...
             'Powers (uW)', num2str(D.powers); 'Reference power (uW)', D.refPower; ...
             'Spectrum', sprintf('mean over the frames of each file; offset: %s', P.offset); ...
             'Brightness ratio', sprintf('total AUC (%g-%g nm) / total AUC at the reference power', P.wlRange); ...
             'Method', D.method; ...
             'Gain', sprintf('spectrum(k) / reference spectrum; k (proportional) where the reference < %g of its maximum; smoothed %g nm', ...
                P.minFrac, P.gainSmoothNm); ...
             'Relative gain', 'gain / ratio: 1 = the wavelength changes in proportion to the total AUC'; ...
             'T1 / T2', sprintf('fixed bands %g-%g and %g-%g nm (AUC and peak search)', P.bandT1, P.bandT2); ...
             'Linearity', 'median R^2 of a straight line of every wavelength vs the ratio (>= 3 powers)'; ...
             'Fixed part', 'two-component model S = B + k E; fixed fraction = B / (B + E) in the T1 and T2 bands'; ...
             'Cross-validation', sprintf(['spectrum at each power predicted without it: decoder (animal) = line through the ' ...
                'animal''s other powers; decoder (pooled) = fixed fraction from the other animals applied to the nearest other ' ...
                'power; proportional = nearest other power x ratio; noise ceiling = half the odd/even frame difference. ' ...
                'Decoding accuracy = 1 - ||pred - meas|| / ||meas|| (%g-%g nm); normalized RMSE = 100 x (1 - accuracy). ' ...
                'AUC errors in CVSummary are absolute values. Tests: paired t-test and exact Wilcoxon signed-rank, n = animals.'], P.evalRange)};
        writecell(clean(C), xlsFile, 'Sheet', 'README');
        np = numel(D.powers);
        hdr = {'Animal', 'Fixed part T1 (fraction)', 'Fixed part T2 (fraction)', 'Linearity R^2 T1', 'Linearity R^2 T2'};
        fl = {'Files', 'Frames', 'Offset', 'Ratio (total AUC)', 'Total AUC', 'T1 AUC', 'T2 AUC', 'T2/T1', 'T1 peak (nm)', 'T2 peak (nm)'};
        for i = 1:numel(fl)
            for j = 1:np, hdr{end+1} = sprintf('%s %g uW', fl{i}, D.powers(j)); end %#ok<AGROW>
        end
        C = hdr;
        for a = 1:numel(D.animals)
            A = D.animals(a);
            row = {A.name, A.fixedT1, A.fixedT2, A.linR2T1, A.linR2T2};
            vals = {A.files, A.nFrames, A.offsets, A.ratio, A.totalAUC, A.aucT1, A.aucT2, A.t2t1, A.peakT1, A.peakT2};
            for i = 1:numel(vals)
                for j = 1:np
                    k = find(A.powers == D.powers(j), 1);
                    if isempty(k), row{end+1} = ''; %#ok<AGROW>
                    elseif iscell(vals{i}), row{end+1} = vals{i}{k}; %#ok<AGROW>
                    else, row{end+1} = vals{i}(k); %#ok<AGROW>
                    end
                end
            end
            C(end+1, :) = row; %#ok<AGROW>
        end
        writecell(clean(C), xlsFile, 'Sheet', 'Animals');
        C = {'Animal', 'Held-out power (uW)', 'Ratio', 'Method', 'Decoding accuracy', 'Normalized RMSE (%)', ...
             'T1 AUC error (%)', 'T2 AUC error (%)', 'T2/T1 error (%)'};
        for a = 1:numel(D.animals)
            c = D.animals(a).cv;
            for i = 1:numel(c)
                C(end+1, :) = {D.animals(a).name, c(i).power, c(i).ratio, c(i).method, c(i).accuracy, c(i).nrmse, ...
                    c(i).t1, c(i).t2, c(i).t2t1}; %#ok<AGROW>
            end
        end
        writecell(clean(C), xlsFile, 'Sheet', 'CrossValidation');
        C = {'Held-out power (uW)', 'Method', 'Measure', 'n animals', 'Mean', 'SEM'};
        for i = 1:numel(D.cv)
            T = D.cv(i);
            C(end+1, :) = {T.power, T.method, T.measure, T.n, T.mean, T.sem}; %#ok<AGROW>
        end
        writecell(clean(C), xlsFile, 'Sheet', 'CVSummary');
        C = {'Held-out power (uW)', 'Measure', 'Comparison', 'n animals', 'Mean (first)', 'Mean (second)', ...
             't (paired)', 'df', 'p (paired t-test)', 'p (Wilcoxon signed-rank, exact)'};
        for i = 1:numel(D.cvTests)
            T = D.cvTests(i);
            C(end+1, :) = {T.power, T.measure, T.comparison, T.n, T.mean1, T.mean2, T.t, T.df, T.p, T.pWilcoxon}; %#ok<AGROW>
        end
        writecell(clean(C), xlsFile, 'Sheet', 'CVTests');
        C = {'Wavelength (nm)'}; M = D.wl(:);
        for a = 1:numel(D.animals)
            A = D.animals(a);
            for j = 1:numel(A.powers)
                if A.powers(j) == D.refPower, continue; end
                C{end+1} = sprintf('%s %g uW', A.name, A.powers(j)); %#ok<AGROW>
                M = [M, A.relGain(:, j)]; %#ok<AGROW>
            end
        end
        writecell(clean([C; num2cell(M)]), xlsFile, 'Sheet', 'RelGain');
        fprintf('   saved %s\n', xlsFile);
    catch ME
        fprintf('   Excel file not written: %s\n', ME.message);
    end
end


%% =====================================================================
%  Spectrum file reader (decoder version: also reads multi-frame two-column exports)
%  =====================================================================
function [wl, Y, nDropped, fmt] = readDyeTxt(file, wlRange)
% Reads a spectrum file in either layout:
%   'column'  "wavelength <tab> value" lines; text lines (e.g. "Time (s): 0.00, Frame: 1",
%             "Wavelength (nm)  Fluorescence Intensity") are ignored. Several frames in one
%             file are separated where the wavelength starts again from the beginning.
%   'row'     OceanView time series (wavelength line, then one row per spectrum);
%             read with readSpectralTxt
    nDropped = 0;
    fid = fopen(file, 'r');
    if fid < 0, error('Cannot open %s', file); end
    two = 0; more = 0;
    for i = 1:300
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        v = sscanf(strrep(strtrim(ln), ',', ' '), '%f');
        if isempty(v) || ~isempty(regexp(strtrim(ln), '[A-DF-Za-df-z]', 'once')), continue; end   % text line
        if numel(v) == 2, two = two + 1; elseif numel(v) > 2, more = more + 1; end
    end
    fclose(fid);
    if two < 10 || more > two
        fmt = 'row';
        [wl, Y, nDropped] = readSpectralTxt(file, wlRange);
        return;
    end
    fmt = 'column';
    txt = fileread(file);
    tk = regexp(txt, '(?m)^\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)[ \t,;]+([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*$', 'tokens');
    if isempty(tk), error('No wavelength / value pairs in %s', file); end
    tk = [tk{:}];
    v = str2double(reshape(tk, 2, []))';         % n x 2
    v = v(all(isfinite(v), 2), :);
    % frames: a new frame starts where the wavelength does not increase
    st = [1; find(diff(v(:, 1)) <= 0) + 1];
    en = [st(2:end) - 1; size(v, 1)];
    len = en - st + 1;
    L = mode(len);
    keepF = find(len == L);
    nDropped = numel(len) - numel(keepF);
    wAll = v(st(keepF(1)):en(keepF(1)), 1);
    keep = wAll >= wlRange(1) & wAll <= wlRange(2);
    wl = wAll(keep);
    Y = zeros(sum(keep), numel(keepF));
    for k = 1:numel(keepF)
        blk = v(st(keepF(k)):en(keepF(k)), :);
        if max(abs(blk(:, 1) - wAll)) > 1e-6
            Y(:, k) = interp1(blk(:, 1), blk(:, 2), wl, 'linear', NaN);
        else
            Y(:, k) = blk(keep, 2);
        end
    end
    Y = Y(:, all(isfinite(Y), 1));
end


%% =====================================================================
%  Helpers (unchanged from NMF_Hybrid_Validation_v2 / NMF_Hybrid_Equivalence)
%  =====================================================================
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


function p = stripSep(p)
    p = strtrim(p);
    while numel(p) > 1 && (p(end) == '\' || p(end) == '/')
        p(end) = [];
    end
end


function s = safeFileName(s)
    s = regexprep(s, '[<>:"/\\|?*]', '_');
end


function [yl, st] = niceAxis(v)
% Limits on ticks of 1-2-5 x 10^k, 4-7 ticks, covering v.
    v = v(isfinite(v));
    lo = min(v); hi = max(v);
    if hi <= lo, hi = lo + 1; end
    best = [];
    for e = -6:6
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


function p = signRankExact(d)
% Exact two-sided Wilcoxon signed-rank test (zero differences dropped).
    d = d(isfinite(d) & d ~= 0);
    n = numel(d);
    p = NaN;
    if n < 1, return; end
    [~, o] = sort(abs(d)); r = zeros(n, 1); r(o) = 1:n;
    % average ranks for ties
    ad = abs(d);
    for i = 1:n, r(i) = mean(find(sort(ad) == ad(i))); end
    W = sum(r(d > 0));
    signs = dec2bin(0:2 ^ n - 1) == '1';
    Wall = signs * r;
    mu = sum(r) / 2;
    p = min(1, mean(abs(Wall - mu) >= abs(W - mu) - 1e-12));
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
