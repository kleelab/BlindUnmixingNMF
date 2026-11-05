function NMF_Hybrid_Equivalence(matFile)
% NMF_Hybrid_Equivalence
% Agreement and equivalence of unmixed vs ground-truth values, from the .mat
% file saved by NMF_Hybrid_Validation (v1) or NMF_Hybrid_Validation_v2.
% No re-analysis: edit the SETTINGS block below and run again to redraw.
%
% Usage
%   NMF_Hybrid_Equivalence                 % choose the .mat file
%   NMF_Hybrid_Equivalence('D:\...\GCaMP+jRGECO_mixture_R2_v2.mat')
%
% Figures (each saved as PDF and EPS, vector, for Illustrator)
%   *_Agreement     ground truth vs unmixed AUC, one panel per sensor; light dots =
%                   in silico mixtures, dark dots = mouse means; dashed = identity (y = x);
%                   Lin's CCC (95% CI) on mouse-level means
%   *_BlandAltman   (unmixed - ground truth) / ground truth (%) vs ground-truth AUC,
%                   one panel per sensor; light dots = pairs, dark dots = mouse means; solid = bias,
%                   dotted = 95% limits of agreement (bias +/- 1.96 SD), mouse level
%   *_TOST          equivalence (two one-sided tests), AUC and delta AUC of each sensor:
%                   mean percentage difference and 90% CI (mouse level);
%                   grey band = equivalence margin
% Excel  *_Equivalence.xlsx   README, Agreement, BlandAltman, TOST, MouseMeans, TrialCV
%
% Sensors: green-sensor metrics use the T1 window, red-sensor metrics the T2 window.
% Sensor + sensor: two sensors (two columns). Sensor + control dye: one sensor.
% Mouse level: every mouse contributes one value, the mean over the mixture
% pairs that include it (pairs are not independent when mice are combined
% with several partners).

%% ===================== SETTINGS (edit, then run again) =====================
% Agreement and Bland-Altman: light dots = individual in silico mixtures,
% dark dots = mouse-level means (the values used for CCC, bias and LoA).
% Light dots are drawn only when there are more pairs than mice (sensor + sensor);
% sensor + control dye (one mixture per mouse): dark dots only.
P.dotSizeAgreementPair     = 8;   % marker size (pt), pairs, agreement plots
P.dotSizeAgreementAnimal   = 8;   % marker size (pt), animal means, agreement plots
P.dotSizeBlandAltmanPair   = 8;   % marker size (pt), pairs, Bland-Altman plots
P.dotSizeBlandAltmanAnimal = 8;   % marker size (pt), animal means, Bland-Altman plots
P.pairTint           = 0.7;    % pair dots: sensor colour mixed with white (0 = sensor colour, 1 = white);
                               % opaque, looks like 50% transparency, keeps PDF/EPS fully vector
P.dotSizeTOST        = 8;      % marker size (pt), mean in the TOST plot
P.dotEdgeColor       = 'none'; % marker edge colour ('none' = no edge, colour only); e.g. [1 1 1] = white edge
P.dotEdgeWidth       = 0.3;    % marker edge width (pt)
P.lineWidth          = 0.75;   % identity, bias, limits of agreement (pt)
P.ciLineWidth        = 1.0;    % TOST confidence-interval bar (pt)
P.axisLineWidth      = 0.5;    % axes (pt)
P.capHeight          = 0.30;   % TOST cap height (fraction of the row spacing)
P.fontName           = 'Arial';
P.fontSize           = 7;      % pt
P.panelCm            = [4.0 4.0];  % agreement and Bland-Altman panel: width, height (cm)
P.agreeTick          = {[], 2.5};  % agreement plot tick step (x and y, axis unit x10^3): {green sensor, red sensor}
                                   % [] = automatic. Red sensor: 0, 2.5, 5, ...
% Bland-Altman axes: {green sensor, red sensor} (applied by sensor, also when only one sensor is present); {} = automatic
% (each panel fitted to its own data). Example: P.baYLim = {[-4 1], [-6 14]};
P.baYLim             = {[-1 2], []};   % y limits (%); [] = automatic (green -1 to 2, red automatic)
P.baMinSpan          = 3;      % automatic axes: minimum y range (percentage points); avoids
                               % exaggerating very small differences (0 = fit the data tightly)
P.baYTick            = {0.5, 1};      % y tick step (%), e.g. {1, 5}; [] = automatic
P.baXFromZero        = true;   % Bland-Altman x axis (ground-truth AUC) starts at 0
P.baXLim             = {};     % x limits per panel in the axis unit (x10^3), e.g. {[0 40], [0 15]}; {} = automatic
P.baXTick            = {[], 2.5};     % x tick step per panel (x10^3); [] = automatic (red: 0, 2.5, 5, ...)
P.baRegLine          = 'all';  % Bland-Altman regression of the mouse-level difference on the ground-truth AUC
                               % (proportional bias): 'all' panels, 'significant' (p < 0.05) only, or 'none'
P.baRegColor         = [0 0 0];    % regression line colour
P.baRegLineWidth     = 1.0;        % regression line width (pt)
P.tostWidthCm        = 4.5;    % TOST plot width (cm)
P.tostRowCm          = 0.6;    % TOST height per row (cm)
P.margin             = 10;     % equivalence margin (+/- %)
P.alpha              = 0.05;   % per one-sided test -> 90% confidence interval
P.marginColor        = [0.90 0.90 0.90];   % equivalence band (solid; EPS has no transparency)
P.refColor           = [0.55 0.55 0.55];   % identity and zero lines
P.colorGreen         = [];     % [] = sensor colour stored in the .mat file
P.colorRed           = [];
P.outDir             = '';     % '' = folder of the .mat file
%% ============================================================================

    if nargin < 1 || isempty(matFile)
        [f, p] = uigetfile({'*.mat', 'Mixture validation results (*.mat)'}, 'Select the .mat file');
        if ~ischar(f), fprintf('Cancelled.\n'); return; end
        matFile = fullfile(p, f);
    end
    L = load(matFile);
    if ~isfield(L, 'R')
        error('%s does not contain the results variable R.', matFile);
    end
    R = L.R;
    [mDir, mBase] = fileparts(matFile);
    if isempty(P.outDir), P.outDir = mDir; end
    if exist(P.outDir, 'dir') ~= 7, mkdir(P.outDir); end
    base = fullfile(P.outDir, mBase);

    cols = R.colors;
    if ~isempty(P.colorGreen), cols{1} = P.colorGreen; end
    if ~isempty(P.colorRed),   cols{2} = P.colorRed;   end

    fprintf('\n==============================================================\n');
    fprintf(' Agreement and equivalence: %s\n', mBase);
    fprintf(' %s, %s, AUC windows: %s\n', R.version, R.mode, R.windowRule);
    fprintf('==============================================================\n');

    % ----- sensors present and their metrics
    keyAUC  = {'green_T1_AUC', 'red_T2_AUC'};
    keyDAUC = {'green_T1_delta_AUC', 'red_T2_delta_AUC'};
    band    = {'T1', 'T2'};
    S = struct('slot', {}, 'name', {}, 'color', {}, 'auc', {}, 'dauc', {}, 'cv', {});
    for s = 1:2
        if ~strcmp(R.types{s}, 'sensor'), continue; end
        k = numel(S) + 1;
        S(k).slot = s;
        S(k).name = R.names{s};
        S(k).color = cols{s};
        S(k).band = band{s};
        S(k).auc = metricData(R, keyAUC{s}, s);
        S(k).dauc = metricData(R, keyDAUC{s}, s);
        S(k).cv = trialCV(R, keyAUC{s}, s);
        fprintf('%s: %d %ss, %d mice; trial-to-trial CV of ground-truth AUC %.1f +/- %.1f %% (mean +/- SEM over mice)\n', ...
            S(k).name, numel(S(k).auc.pairRaw), lower(R.unit), numel(S(k).auc.animal), S(k).cv.meanCV, S(k).cv.semCV);
    end
    if isempty(S)
        error('No sensor in %s.', matFile);
    end

    % ----- statistics (mouse level)
    for k = 1:numel(S)
        S(k).ccc = lincc(S(k).auc.animalRaw, S(k).auc.animalNmf);
        S(k).ba = blandAltman(S(k).auc.animalPct);
        S(k).tostAUC = tost(S(k).auc.animalPct, P.margin, P.alpha);
        S(k).tostDAUC = tost(S(k).dauc.animalPct, P.margin, P.alpha);
        S(k).reg = baRegression(S(k).auc.animalRaw, S(k).auc.animalPct);
        fprintf('%s %s AUC : CCC %.4f (95%% CI %.4f to %.4f); bias %+.2f %% (LoA %+.2f to %+.2f %%)\n', ...
            S(k).name, S(k).band, S(k).ccc.ccc, S(k).ccc.lo, S(k).ccc.hi, S(k).ba.bias, S(k).ba.lo, S(k).ba.hi);
        fprintf('   TOST AUC       : %+.2f %% (90%% CI %+.2f to %+.2f), p = %.4g -> %s\n', S(k).tostAUC.mean, ...
            S(k).tostAUC.ciLo, S(k).tostAUC.ciHi, S(k).tostAUC.p, yesno(S(k).tostAUC.equivalent));
        fprintf('   TOST delta AUC : %+.2f %% (90%% CI %+.2f to %+.2f), p = %.4g -> %s\n', S(k).tostDAUC.mean, ...
            S(k).tostDAUC.ciLo, S(k).tostDAUC.ciHi, S(k).tostDAUC.p, yesno(S(k).tostDAUC.equivalent));
        G = S(k).reg;
        fprintf('   Proportional bias (difference vs ground-truth AUC, mouse level): slope %.3g %%/10^3 (95%% CI %.3g to %.3g), R^2 %.2f, p = %.4g (n = %d)\n', ...
            1e3 * G.slope, 1e3 * G.ciLo, 1e3 * G.ciHi, G.r2, G.p, G.n);
    end

    % ----- figures
    figAgreement(S, R, P, [base '_Agreement']);
    figBlandAltman(S, R, P, [base '_BlandAltman']);
    figTOST(S, P, [base '_TOST']);

    % ----- Excel
    writeStats(S, R, P, matFile, [base '_Equivalence.xlsx']);
    fprintf('Done. Output folder: %s\n', P.outDir);
end


%% =====================================================================
%  Data
%  =====================================================================
function D = metricData(R, key, slot)
% Pair values and mouse-level means of one metric.
    H = R.hybrids;
    n = numel(H);
    D.pairRaw = nan(n, 1);
    D.pairNmf = nan(n, 1);
    D.pairAnimal = cell(n, 1);
    for i = 1:n
        m = H(i).cmp.m;
        q = find(strcmp({m.key}, key), 1);
        if isempty(q), continue; end
        D.pairRaw(i) = m(q).raw;
        D.pairNmf(i) = m(q).nmf;
        if slot == 1, D.pairAnimal{i} = H(i).green; else, D.pairAnimal{i} = H(i).red; end
    end
    ok = isfinite(D.pairRaw) & isfinite(D.pairNmf) & D.pairRaw ~= 0;
    D.pairRaw = D.pairRaw(ok);
    D.pairNmf = D.pairNmf(ok);
    D.pairAnimal = D.pairAnimal(ok);
    D.pairPct = 100 * (D.pairNmf - D.pairRaw) ./ D.pairRaw;
    % mice in order of first appearance
    D.animal = {};
    g = zeros(numel(D.pairAnimal), 1);
    for i = 1:numel(D.pairAnimal)
        j = find(strcmp(D.animal, D.pairAnimal{i}), 1);
        if isempty(j)
            D.animal{end+1, 1} = D.pairAnimal{i};
            j = numel(D.animal);
        end
        g(i) = j;
    end
    na = numel(D.animal);
    D.animalRaw = zeros(na, 1);
    D.animalNmf = zeros(na, 1);
    D.animalPct = zeros(na, 1);
    D.animalNPairs = zeros(na, 1);
    D.animalSdNmf = zeros(na, 1);
    for a = 1:na
        w = g == a;
        D.animalRaw(a) = mean(D.pairRaw(w));
        D.animalNmf(a) = mean(D.pairNmf(w));
        D.animalPct(a) = mean(D.pairPct(w));
        D.animalNPairs(a) = sum(w);
        if sum(w) > 1, D.animalSdNmf(a) = std(D.pairNmf(w)); else, D.animalSdNmf(a) = 0; end
    end
end


function C = trialCV(R, key, slot)
% Trial-to-trial coefficient of variation of the ground-truth AUC, per mouse
% (pair with the most trials of that mouse), then mean +/- SEM across mice.
    H = R.hybrids;
    names = {};
    cv = [];
    nTr = [];
    for i = 1:numel(H)
        m = H(i).cmp.m;
        q = find(strcmp({m.key}, key), 1);
        if isempty(q) || ~isfield(m, 'rawTrials') || isempty(m(q).rawTrials), continue; end
        v = m(q).rawTrials(isfinite(m(q).rawTrials));
        if numel(v) < 2, continue; end
        if slot == 1, an = H(i).green; else, an = H(i).red; end
        j = find(strcmp(names, an), 1);
        if isempty(j)
            names{end+1} = an; %#ok<AGROW>
            cv(end+1) = 100 * std(v) / mean(v); %#ok<AGROW>
            nTr(end+1) = numel(v); %#ok<AGROW>
        elseif numel(v) > nTr(j)
            cv(j) = 100 * std(v) / mean(v);
            nTr(j) = numel(v);
        end
    end
    C.animal = names(:);
    C.cv = cv(:);
    C.nTrials = nTr(:);
    [C.meanCV, C.semCV] = meanSem(cv);
end


%% =====================================================================
%  Statistics
%  =====================================================================
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


function G = baRegression(x, y)
% Least-squares line y = intercept + slope * x on the mouse-level values (proportional bias,
% Bland & Altman 1999); slope tested against 0 (t-test, df = n - 2), 95% CI of the slope.
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
% Bias and 95% limits of agreement (bias +/- 1.96 SD).
    d = d(isfinite(d));
    B.n = numel(d);
    B.bias = mean(d);
    B.sd = std(d);
    B.lo = B.bias - 1.96 * B.sd;
    B.hi = B.bias + 1.96 * B.sd;
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


function [m, s] = meanSem(v)
    v = v(isfinite(v));
    if isempty(v), m = NaN; s = NaN;
    elseif numel(v) == 1, m = v; s = 0;
    else, m = mean(v); s = std(v) / sqrt(numel(v));
    end
end


%% =====================================================================
%  Figures
%  =====================================================================
function figAgreement(S, R, P, fileBase)
    nS = numel(S);
    gx = 1.6; x0 = 1.5; y0 = 1.2; top = 1.4;
    w = x0 + nS * P.panelCm(1) + (nS - 1) * gx + 0.5;
    h = y0 + P.panelCm(2) + top;
    f = newFigure(w, h);
    for k = 1:nS
        D = S(k).auc;
        q = S(k).slot;                       % 1 = green sensor, 2 = red sensor (axis settings)
        [sc, unitTxt] = auSc(max([D.pairRaw; D.pairNmf]));
        x = D.pairRaw / sc; y = D.pairNmf / sc;
        hi = niceUpper(1.05 * max([x; y]));
        lim = [0 hi];
        st = niceStep(hi);
        if numel(P.agreeTick) >= q && ~isempty(P.agreeTick{q})
            st = P.agreeTick{q};
            hi = ceil(hi / st - 1e-9) * st;
            lim = [0 hi];
        end
        ax = newAxes(f, [x0 + (k - 1) * (P.panelCm(1) + gx), y0, P.panelCm], P);
        plot(ax, lim, lim, '--', 'Color', P.refColor, 'LineWidth', P.lineWidth);
        drawDots(ax, x, y, D.animalRaw / sc, D.animalNmf / sc, S(k).color, ...
            P.dotSizeAgreementPair, P.dotSizeAgreementAnimal, numel(D.pairRaw) > numel(D.animal), P);
        set(ax, 'XLim', lim, 'YLim', lim, 'XTick', 0:st:hi, 'YTick', 0:st:hi);
        xlabel(ax, ['Ground-truth AUC (' unitTxt ')']);
        ylabel(ax, ['Unmixed AUC (' unitTxt ')']);
        C = S(k).ccc;
        txt = {sprintf('CCC = %.4f', C.ccc), sprintf('95%% CI %.4f to %.4f', C.lo, C.hi), nText(D, R)};
        for j = 1:numel(txt)          % one text object per line (clean in Illustrator)
            text(ax, 0.04, 0.97 - 0.085 * (j - 1), txt{j}, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'FontName', P.fontName, 'FontSize', P.fontSize);
        end
        panelTitle(ax, S(k), 'AUC', P);
    end
    saveFigure(f, fileBase);
end


function figBlandAltman(S, R, P, fileBase)
    nS = numel(S);
    gx = 1.6; x0 = 1.5; y0 = 2.4; top = 1.9;
    w = x0 + nS * P.panelCm(1) + (nS - 1) * gx + 0.5;
    h = y0 + P.panelCm(2) + top;
    f = newFigure(w, h);
    for k = 1:nS
        D = S(k).auc; B = S(k).ba;
        q = S(k).slot;                       % 1 = green sensor, 2 = red sensor (axis settings)
        % y range of this panel: automatic (data, limits of agreement and 0)
        % unless set in P.baYLim / P.baYTick
        [yl, st] = baLimits([D.pairPct; B.lo; B.hi; 0], P.baMinSpan);
        if numel(P.baYLim) >= q && numel(P.baYLim{q}) == 2
            yl = P.baYLim{q};
            st = niceStep(diff(yl));
        end
        if numel(P.baYTick) >= q && ~isempty(P.baYTick{q})
            st = P.baYTick{q};
        end
        [sc, unitTxt] = auSc(max(D.pairRaw));
        x = D.pairRaw / sc;
        hi = niceUpper(1.05 * max(x));
        if P.baXFromZero
            lo = 0;
        else
            lo = max(0, niceLower(0.95 * min(x), niceStep(hi)));
        end
        if numel(P.baXLim) >= q && numel(P.baXLim{q}) == 2
            lo = P.baXLim{q}(1); hi = P.baXLim{q}(2);
        end
        xl = [lo hi];
        ax = newAxes(f, [x0 + (k - 1) * (P.panelCm(1) + gx), y0, P.panelCm], P);
        plot(ax, xl, [0 0], '-', 'Color', P.refColor, 'LineWidth', P.axisLineWidth);
        plot(ax, xl, [B.bias B.bias], '-', 'Color', S(k).color, 'LineWidth', P.lineWidth);
        plot(ax, xl, [B.lo B.lo], ':', 'Color', S(k).color, 'LineWidth', P.lineWidth);
        plot(ax, xl, [B.hi B.hi], ':', 'Color', S(k).color, 'LineWidth', P.lineWidth);
        drawDots(ax, x, D.pairPct, D.animalRaw / sc, D.animalPct, S(k).color, ...
            P.dotSizeBlandAltmanPair, P.dotSizeBlandAltmanAnimal, numel(D.pairRaw) > numel(D.animal), P);
        % regression of the mouse-level difference on the ground-truth AUC (within the data range)
        G = S(k).reg;
        drawReg = strcmpi(P.baRegLine, 'all') || (strcmpi(P.baRegLine, 'significant') && G.p < 0.05);
        if drawReg && isfinite(G.slope)
            xr = [G.xMin G.xMax];
            plot(ax, xr / sc, G.intercept + G.slope * xr, '-', 'Color', P.baRegColor, 'LineWidth', P.baRegLineWidth);
        end
        sx = niceStep(hi - lo);
        if numel(P.baXTick) >= q && ~isempty(P.baXTick{q}), sx = P.baXTick{q}; end
        set(ax, 'XLim', xl, 'YLim', yl, 'XTick', ceil(lo / sx) * sx:sx:hi, 'YTick', tickList(yl, st));
        xlabel(ax, ['Ground-truth AUC (' unitTxt ')']);
        ylabel(ax, 'Difference (%)');
        txt = {sprintf('bias %+.2f%%', B.bias), sprintf('95%% LoA %+.2f to %+.2f%%', B.lo, B.hi)};
        for j = 1:numel(txt)
            text(ax, 0.5, 1.12 - 0.085 * (j - 1), txt{j}, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontName', P.fontName, 'FontSize', P.fontSize, 'Clipping', 'off');
        end
        panelTitle(ax, S(k), 'AUC', P, 1.24);
        % regression statistics and n below the x-axis label (keeps the plot area free)
        if ~strcmpi(P.baRegLine, 'none') && isfinite(G.slope)
            e = round(log10(sc));
            dd = max(0, 1 - floor(log10(abs(sc * G.slope) + eps)));      % 2 significant digits of the slope
            fm = sprintf('%%.%df', dd);
            rt = {sprintf(['slope ' fm ' [' fm ', ' fm '] %%/10^{%d}'], sc * G.slope, sc * G.ciLo, sc * G.ciHi, e), ...
                  sprintf('R^2 = %.2f, %s', G.r2, pText(G.p))};
            for j = 1:2
                text(ax, 0.5, -0.27 - 0.09 * (j - 1), rt{j}, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'top', 'FontName', P.fontName, 'FontSize', P.fontSize, 'Clipping', 'off');
            end
        end
        text(ax, 1, -0.49, nText(D, R), 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'FontName', P.fontName, 'FontSize', P.fontSize - 1, ...
            'Color', P.refColor, 'Clipping', 'off');
    end
    saveFigure(f, fileBase);
end


function figTOST(S, P, fileBase)
% Rows (top to bottom): sensor 1 AUC, sensor 1 delta AUC, sensor 2 AUC, sensor 2 delta AUC.
    rows = struct('label', {}, 'T', {}, 'color', {});
    for k = 1:numel(S)
        rows(end+1) = struct('label', sprintf('%s %s AUC', S(k).name, S(k).band), ...
            'T', S(k).tostAUC, 'color', S(k).color); %#ok<AGROW>
        rows(end+1) = struct('label', sprintf('%s %s \\DeltaAUC', S(k).name, S(k).band), ...
            'T', S(k).tostDAUC, 'color', S(k).color); %#ok<AGROW>
    end
    nR = numel(rows);
    xLab = 3.4; xRight = 2.2; y0 = 1.3; top = 0.6;
    w = xLab + P.tostWidthCm + xRight;
    hAx = nR * P.tostRowCm;
    h = y0 + hAx + top;
    f = newFigure(w, h);
    ax = newAxes(f, [xLab, y0, P.tostWidthCm, hAx], P);
    m = P.margin;
    ciAll = [arrayfun(@(r) r.T.ciLo, rows), arrayfun(@(r) r.T.ciHi, rows)];
    xm = max([1.5 * m, abs(ciAll(isfinite(ciAll))) * 1.1]);
    st = niceStep(2 * xm);
    xm = ceil(xm / st) * st;
    yl = [0.5, nR + 0.5];
    patch(ax, [-m m m -m], [yl(1) yl(1) yl(2) yl(2)], P.marginColor, 'EdgeColor', 'none');
    plot(ax, [0 0], yl, '-', 'Color', P.refColor, 'LineWidth', P.axisLineWidth);
    cap = P.capHeight / 2;
    for i = 1:nR
        yy = nR + 1 - i;
        T = rows(i).T; c = rows(i).color;
        if isfinite(T.ciLo)
            plot(ax, [T.ciLo T.ciHi], [yy yy], '-', 'Color', c, 'LineWidth', P.ciLineWidth);
            plot(ax, [T.ciLo T.ciLo], yy + [-cap cap], '-', 'Color', c, 'LineWidth', P.ciLineWidth);
            plot(ax, [T.ciHi T.ciHi], yy + [-cap cap], '-', 'Color', c, 'LineWidth', P.ciLineWidth);
        end
        plot(ax, T.mean, yy, 'o', 'MarkerSize', P.dotSizeTOST, 'MarkerFaceColor', c, ...
            'MarkerEdgeColor', 'none', 'LineStyle', 'none');
        if T.p < 0.001, ptxt = 'p < 0.001'; else, ptxt = sprintf('p = %.3f', T.p); end
        text(ax, xm * 1.08, yy, sprintf('%s (n = %d)', ptxt, T.n), 'FontName', P.fontName, ...
            'FontSize', P.fontSize, 'VerticalAlignment', 'middle', 'Clipping', 'off');
    end
    % row labels (sensor colour) left of the plot; no y axis
    for i = 1:nR
        c = rows(i).color;
        text(ax, -xm - 0.06 * 2 * xm, nR + 1 - i, rows(i).label, 'Color', c, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'FontName', P.fontName, 'FontSize', P.fontSize, 'Clipping', 'off');
    end
    set(ax, 'XLim', [-xm xm], 'YLim', yl, 'XTick', -xm:st:xm, 'YTick', [], 'YColor', 'none');
    xlabel(ax, 'Difference (%), mean and 90% CI');
    text(ax, 0, yl(2), sprintf('equivalence margin \\pm%g%%', m), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontName', P.fontName, 'FontSize', P.fontSize - 1, ...
        'Color', P.refColor, 'Clipping', 'off');
    saveFigure(f, fileBase);
end


function drawDots(ax, xp, yp, xa, ya, col, sizePair, sizeAnimal, showPairs, P)
% Pairs as light (tinted, opaque) dots underneath; mouse means as dark dots on top.
    if showPairs
        light = col + (1 - col) * P.pairTint;
        plot(ax, xp, yp, 'o', 'MarkerSize', sizePair, 'MarkerFaceColor', light, ...
            'MarkerEdgeColor', P.dotEdgeColor, 'LineWidth', P.dotEdgeWidth, 'LineStyle', 'none');
    end
    plot(ax, xa, ya, 'o', 'MarkerSize', sizeAnimal, 'MarkerFaceColor', col, ...
        'MarkerEdgeColor', P.dotEdgeColor, 'LineWidth', P.dotEdgeWidth, 'LineStyle', 'none');
end


function s = num2s(v)
% two significant digits, e.g. -0.28, -0.0041, -3
    if ~isfinite(v), s = 'NaN'; return; end
    s = sprintf('%.2g', v);
    if any(s == 'e'), s = sprintf('%.2f', v); end
end


function s = pText(p)
    if ~isfinite(p), s = 'p = n/a';
    elseif p < 0.001, s = 'p < 0.001';
    else, s = sprintf('p = %.3f', p);
    end
end


function t = nText(D, R)
    na = numel(D.animal); np = numel(D.pairRaw);
    if np > na
        t = sprintf('n = %d mice (%d %ss)', na, np, lower(R.unit));
    else
        t = sprintf('n = %d mice', na);
    end
end


function panelTitle(ax, Sk, what, P, y)
    if nargin < 5, y = 1.04; end
    c = Sk.color;
    text(ax, 0.5, y, sprintf('\\color[rgb]{%.3f %.3f %.3f}%s\\color[rgb]{0 0 0} %s %s', ...
        c(1), c(2), c(3), Sk.name, Sk.band, what), 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontName', P.fontName, ...
        'FontSize', P.fontSize + 1, 'Clipping', 'off');
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


function [sc, txt] = auSc(mx)
% Axis scale for AUC values: 10^3 steps.
    e = 3 * floor(log10(max(mx, 1)) / 3);
    if e < 3, e = 0; end
    sc = 10 ^ e;
    if e == 0
        txt = 'counts \cdot nm';
    else
        txt = sprintf('10^{%d} counts \\cdot nm', e);
    end
end


function [yl, st] = baLimits(v, minSpan)
% Automatic y range. Limits always fall on tick values; if all values are
% >= 0 (or <= 0) the axis starts (ends) at 0. Tick steps 1-2-5 x 10^k
% (2.5 x 10^k only if nothing else fits), 4-6 ticks, range closest to the data.
% Very small differences are shown on at least minSpan percentage points.
    v = v(isfinite(v));
    if isempty(v), yl = [-1 1]; st = 0.5; return; end
    if nargin < 2 || isempty(minSpan), minSpan = 0; end
    r = max(max(v) - min(v), 0.1);
    lo0 = min(v) - 0.015 * r;
    hi0 = max(v) + 0.015 * r;
    if hi0 - lo0 < minSpan
        c = (min(v) + max(v)) / 2;
        lo0 = c - 0.485 * minSpan;
        hi0 = c + 0.485 * minSpan;
    else
        if min(v) >= 0, lo0 = 0; end       % no negative padding below 0
        if max(v) <= 0, hi0 = 0; end
    end
    r = hi0 - lo0;
    best = inf; yl = [lo0 hi0]; st = niceStep(r);
    e = floor(log10(r)) - 2:floor(log10(r)) + 1;
    for p = 10 .^ e
        for m = [1 2 2.5 5]
            s = m * p;
            lo = floor(lo0 / s + 1e-9) * s;
            hi = ceil(hi0 / s - 1e-9) * s;
            nt = numel(tickList([lo hi], s));
            if nt < 4 || nt > 6, continue; end
            score = (hi - lo) / r + 0.15 * (m == 2.5);
            if score < best
                best = score; yl = [lo hi]; st = s;
            end
        end
    end
end


function t = tickList(lim, st)
    t = (ceil(lim(1) / st - 1e-9):floor(lim(2) / st + 1e-9)) * st;
    t(abs(t) < 1e-9 * st) = 0;
end


function s = niceStep(range)
% 1-2-2.5-5 x 10^k step giving about 4 intervals.
    raw = range / 4;
    if ~(raw > 0), s = 1; return; end
    p = 10 ^ floor(log10(raw));
    m = [1 2 2.5 5 10] * p;
    s = m(find(m >= raw, 1));
end


function v = niceUpper(x)
    s = niceStep(x);
    v = ceil(x / s) * s;
end


function v = niceLower(x, s)
    v = floor(x / s) * s;
end


function r = rangeOf(v)
    r = max(v) - min(v);
end


function t = yesno(b)
    if b, t = 'equivalent'; else, t = 'not equivalent'; end
end


%% =====================================================================
%  Excel
%  =====================================================================
function writeStats(S, R, P, matFile, xlsFile)
    if exist(xlsFile, 'file') == 2
        try, delete(xlsFile); catch, end %#ok<NOCOM>
    end
    try
        C = { ...
            'Agreement and equivalence of unmixed vs ground-truth values (in silico mixtures of measured single-sensor data)', ''; ...
            'Source (.mat)', matFile; ...
            'Analysis version', sprintf('%s (AUC windows: %s)', R.version, R.windowRule); ...
            'Mode', R.mode; ...
            'Date', datestr(now, 'yyyy-mm-dd HH:MM:SS'); ...
            '', ''; ...
            'Percentage difference', '(unmixed - ground truth) / ground truth x 100'; ...
            'Mouse level', 'one value per mouse: mean over the in silico mixtures that include the mouse; all statistics use mouse-level values'; ...
            'Agreement', 'Lin''s concordance correlation coefficient (CCC) of mouse-level AUC, 95% CI by Fisher z (Lin 1989)'; ...
            'Bland-Altman', 'bias = mean mouse-level percentage difference; 95% limits of agreement = bias +/- 1.96 SD (Bland & Altman 1986)'; ...
            'TOST', sprintf('two one-sided t-tests against +/-%g%% (alpha = %g each); equivalent when the %g%% CI lies within the margin (Schuirmann 1987)', ...
                P.margin, P.alpha, 100 * (1 - 2 * P.alpha)); ...
            'Trial-to-trial CV', 'SD / mean x 100 of the trial-level ground-truth AUC of each mouse (mixture with the most trials), then mean and SEM across mice'; ...
            'Proportional bias', 'least-squares regression of the mouse-level difference (%) on the mouse-level ground-truth AUC; slope tested against 0 (t-test, df = n - 2), 95% CI of the slope (Bland & Altman 1999)'; ...
            'Sheets', 'Agreement, BlandAltman, ProportionalBias, TOST, MouseMeans, TrialCV'};
        writecell(clean(C), xlsFile, 'Sheet', 'README');

        C = {'Sensor', 'Measure', 'n mice', 'n mixtures', 'CCC', 'CCC 95% CI lower', 'CCC 95% CI upper', 'Pearson r'};
        for k = 1:numel(S)
            C(end+1, :) = {S(k).name, [S(k).band ' AUC'], numel(S(k).auc.animal), numel(S(k).auc.pairRaw), ...
                S(k).ccc.ccc, S(k).ccc.lo, S(k).ccc.hi, getf(S(k).ccc, 'r')}; %#ok<AGROW>
        end
        writecell(clean(C), xlsFile, 'Sheet', 'Agreement');

        C = {'Sensor', 'Measure', 'n mice', 'Bias (%)', 'SD (%)', '95% LoA lower (%)', '95% LoA upper (%)'};
        for k = 1:numel(S)
            B = S(k).ba;
            C(end+1, :) = {S(k).name, [S(k).band ' AUC'], B.n, B.bias, B.sd, B.lo, B.hi}; %#ok<AGROW>
        end
        writecell(clean(C), xlsFile, 'Sheet', 'BlandAltman');

        C = {'Sensor', 'Measure', 'n mice', 'Slope (% per 10^3 counts nm)', 'Slope 95% CI lower', 'Slope 95% CI upper', ...
             'Intercept (%)', 'R^2', 't', 'df', 'p (slope = 0)', 'Significant (p < 0.05)', 'AUC min', 'AUC max'};
        for k = 1:numel(S)
            G = S(k).reg;
            C(end+1, :) = {S(k).name, [S(k).band ' AUC'], G.n, 1e3 * G.slope, 1e3 * G.ciLo, 1e3 * G.ciHi, G.intercept, ...
                G.r2, G.t, G.df, G.p, ternaryS(G.p < 0.05, 'yes', 'no'), G.xMin, G.xMax}; %#ok<AGROW>
        end
        writecell(clean(C), xlsFile, 'Sheet', 'ProportionalBias');

        C = {'Sensor', 'Measure', 'n mice', 'Mean difference (%)', 'SD (%)', 'SE (%)', 'df', ...
             sprintf('%g%% CI lower (%%)', 100 * (1 - 2 * P.alpha)), sprintf('%g%% CI upper (%%)', 100 * (1 - 2 * P.alpha)), ...
             'Margin (+/- %)', 't (H0: mean <= -margin)', 'p lower', 't (H0: mean >= +margin)', 'p upper', ...
             'p TOST', 'Equivalent'};
        for k = 1:numel(S)
            for which = 1:2
                if which == 1, T = S(k).tostAUC; ms = [S(k).band ' AUC']; else, T = S(k).tostDAUC; ms = [S(k).band ' delta AUC']; end
                C(end+1, :) = {S(k).name, ms, T.n, T.mean, T.sd, T.se, T.df, T.ciLo, T.ciHi, T.margin, ...
                    T.tLower, T.pLower, T.tUpper, T.pUpper, T.p, yesno(T.equivalent)}; %#ok<AGROW>
            end
        end
        writecell(clean(C), xlsFile, 'Sheet', 'TOST');

        C = {'Sensor', 'Measure', 'Mouse', 'n mixtures', 'Ground-truth mean', 'Unmixed mean', ...
             'SD of unmixed across mixtures', 'Difference (%), mean over mixtures'};
        for k = 1:numel(S)
            for which = 1:2
                if which == 1, D = S(k).auc; ms = [S(k).band ' AUC']; else, D = S(k).dauc; ms = [S(k).band ' delta AUC']; end
                for a = 1:numel(D.animal)
                    C(end+1, :) = {S(k).name, ms, D.animal{a}, D.animalNPairs(a), D.animalRaw(a), ...
                        D.animalNmf(a), D.animalSdNmf(a), D.animalPct(a)}; %#ok<AGROW>
                end
            end
        end
        writecell(clean(C), xlsFile, 'Sheet', 'MouseMeans');

        C = {'Sensor', 'Measure', 'Mouse', 'n trials', 'CV of ground-truth AUC (%)'};
        for k = 1:numel(S)
            V = S(k).cv;
            for a = 1:numel(V.animal)
                C(end+1, :) = {S(k).name, [S(k).band ' AUC'], V.animal{a}, V.nTrials(a), V.cv(a)}; %#ok<AGROW>
            end
            C(end+1, :) = {S(k).name, [S(k).band ' AUC'], 'mean across mice', numel(V.cv), V.meanCV}; %#ok<AGROW>
            C(end+1, :) = {S(k).name, [S(k).band ' AUC'], 'SEM across mice', numel(V.cv), V.semCV}; %#ok<AGROW>
        end
        writecell(clean(C), xlsFile, 'Sheet', 'TrialCV');
        fprintf('   saved %s\n', xlsFile);
    catch ME
        fprintf('   Excel file not written: %s\n', ME.message);
    end
end


function v = ternaryS(c, a, b)
    if c, v = a; else, v = b; end
end


function v = getf(s, f)
    if isfield(s, f), v = s.(f); else, v = NaN; end
end


function C = clean(C)
    C(cellfun(@isempty, C)) = {''};
    for k = 1:numel(C)
        if islogical(C{k}), C{k} = double(C{k}); end
        if isnumeric(C{k}) && isscalar(C{k}) && ~isfinite(C{k}), C{k} = 'NaN'; end
    end
end
