function NMF_RI_Photometry
    % =====================================================================
    % NMF-RI Photometry v2
    % Blind spectral unmixing for multi-component fiber photometry
    % =====================================================================
    % Two-component blind spectral unmixing of mixed emission spectra.
    %   - Initial spectra A0 are supplied by the user as CSV files
    %     (any source: theoretical spectra, in vitro measurements,
    %      single-expression reference recordings, etc.)
    %   - NMF-RI estimates the final spectra A and time-series H
    %
    % Required files
    %   NMF_RI.m
    %   dataPreprocessing.m
    %   Two emission-spectrum CSV files (one per component)
    % =====================================================================

    nComp = 2;

    fig = uifigure( ...
        'Name', 'NMF-RI v2 Blind Spectral Unmixing - Fiber Photometry', ...
        'Position', [70, 50, 1320, 850]);

    tabGroup = uitabgroup(fig, 'Position', [20, 20, 1280, 810]);
    tab1 = uitab(tabGroup, 'Title', '1. Data');
    tab2 = uitab(tabGroup, 'Title', '2. NMF-RI Settings');
    tab3 = uitab(tabGroup, 'Title', '3. Results');

    %% App state
    appState = struct();
    appState.mixedFile       = '';
    appState.comp1File       = '';
    appState.comp2File       = '';
    appState.Y               = [];   % [W x T]
    appState.wavelengths     = [];   % cropped wavelength
    appState.allWavelengths  = [];   % original wavelength
    appState.wlMask          = [];   % crop mask in original wavelength
    appState.rawTimeInfo     = [];   % Date / Time / Timestamp
    appState.A               = [];   % [W x 2]
    appState.H               = [];   % [2 x T]
    appState.A0              = [];   % [W x 2]
    appState.nIter           = 0;
    appState.YPrepared       = [];

    %% ===================== Tab 1: Data =====================
    uilabel(tab1, ...
        'Position', [30, 755, 600, 26], ...
        'Text', 'Data input and initial spectra', ...
        'FontWeight', 'bold', ...
        'FontSize', 15);

    % Mixed TXT
    uilabel(tab1, ...
        'Position', [30, 715, 500, 24], ...
        'Text', '1. Mixed spectral time series (raw TXT)', ...
        'FontWeight', 'bold');

    uibutton(tab1, ...
        'Position', [30, 684, 220, 28], ...
        'Text', 'Browse mixed TXT...', ...
        'ButtonPushedFcn', @(b,e) browse_mixed());

    lblMixed = uilabel(tab1, ...
        'Position', [265, 684, 570, 28], ...
        'Text', 'No file selected', ...
        'FontColor', [0.5 0.5 0.5]);

    % Component 1 CSV
    uilabel(tab1, ...
        'Position', [30, 645, 650, 24], ...
        'Text', '2. Component 1 CSV (emission spectrum, initial A0)', ...
        'FontWeight', 'bold');

    uibutton(tab1, ...
        'Position', [30, 614, 220, 28], ...
        'Text', 'Browse Component 1 CSV...', ...
        'ButtonPushedFcn', @(b,e) browse_comp1());

    lblComp1 = uilabel(tab1, ...
        'Position', [265, 614, 570, 28], ...
        'Text', 'No file selected', ...
        'FontColor', [0.5 0.5 0.5]);

    % Component 2 CSV
    uilabel(tab1, ...
        'Position', [30, 575, 650, 24], ...
        'Text', '3. Component 2 CSV (emission spectrum, initial A0)', ...
        'FontWeight', 'bold');

    uibutton(tab1, ...
        'Position', [30, 544, 220, 28], ...
        'Text', 'Browse Component 2 CSV...', ...
        'ButtonPushedFcn', @(b,e) browse_comp2());

    lblComp2 = uilabel(tab1, ...
        'Position', [265, 544, 570, 28], ...
        'Text', 'No file selected', ...
        'FontColor', [0.5 0.5 0.5]);

    % Wavelength range
    uilabel(tab1, ...
        'Position', [30, 500, 200, 24], ...
        'Text', 'Wavelength range (nm)', ...
        'FontWeight', 'bold');

    uilabel(tab1, 'Position', [30, 470, 45, 24], 'Text', 'Min');
    txtWLMin = uieditfield(tab1, 'text', ...
        'Position', [80, 470, 95, 24], ...
        'Value', '450');

    uilabel(tab1, 'Position', [195, 470, 45, 24], 'Text', 'Max');
    txtWLMax = uieditfield(tab1, 'text', ...
        'Position', [245, 470, 95, 24], ...
        'Value', '650');

    % Component names
    uilabel(tab1, ...
        'Position', [30, 430, 200, 24], ...
        'Text', 'Component names', ...
        'FontWeight', 'bold');

    uilabel(tab1, 'Position', [30, 400, 100, 24], 'Text', 'Component 1');
    txtName1 = uieditfield(tab1, 'text', ...
        'Position', [135, 400, 160, 24], ...
        'Value', 'Component 1');

    uilabel(tab1, 'Position', [30, 370, 100, 24], 'Text', 'Component 2');
    txtName2 = uieditfield(tab1, 'text', ...
        'Position', [135, 370, 160, 24], ...
        'Value', 'Component 2');

    % Preview axes
    axInit = uiaxes(tab1, 'Position', [455, 330, 790, 420]);
    title(axInit, 'Initial spectra A0');
    xlabel(axInit, 'Wavelength (nm)');
    ylabel(axInit, 'Sum-normalized intensity');

    uibutton(tab1, ...
        'Position', [30, 315, 380, 38], ...
        'Text', 'Load data and preview initial spectra', ...
        'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(b,e) load_data());

    txtDataStatus = uilabel(tab1, ...
        'Position', [30, 255, 390, 50], ...
        'Text', '', ...
        'WordWrap', 'on', ...
        'FontColor', [0.25 0.25 0.25]);

    uibutton(tab1, ...
        'Position', [30, 210, 230, 36], ...
        'Text', 'Continue to Settings >>', ...
        'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(b,e) set(tabGroup, 'SelectedTab', tab2));

    %% ===================== Tab 2: NMF-RI Settings =====================
    uilabel(tab2, ...
        'Position', [30, 755, 520, 26], ...
        'Text', 'NMF-RI parameters (Jimenez-Sanchez et al. 2020)', ...
        'FontWeight', 'bold', ...
        'FontSize', 15);

    uilabel(tab2, ...
        'Position', [30, 690, 470, 56], ...
        'WordWrap', 'on', ...
        'FontColor', [0.35 0.35 0.35], ...
        'Text', ['A is estimated from the recording, then the final H ' ...
                 'is computed over the entire time series.']);

    labelX = 30;
    editX  = 335;
    editW  = 115;

    uilabel(tab2, 'Position', [labelX, 640, 295, 24], ...
        'Text', 'theta1 (intra-layer convergence)');
    txtTheta1 = uieditfield(tab2, 'text', ...
        'Position', [editX, 640, editW, 24], ...
        'Value', '0.01');

    uilabel(tab2, 'Position', [labelX, 602, 295, 24], ...
        'Text', 'theta2 (inter-layer convergence)');
    txtTheta2 = uieditfield(tab2, 'text', ...
        'Position', [editX, 602, editW, 24], ...
        'Value', '0.01');

    uilabel(tab2, 'Position', [labelX, 564, 295, 24], ...
        'Text', 'alphaA (spectral sparsity weight)');
    txtAlphaA = uieditfield(tab2, 'text', ...
        'Position', [editX, 564, editW, 24], ...
        'Value', '0.001');

    uilabel(tab2, 'Position', [labelX, 526, 295, 24], ...
        'Text', 'alphaX (concentration sparsity weight)');
    txtAlphaX = uieditfield(tab2, 'text', ...
        'Position', [editX, 526, editW, 24], ...
        'Value', '0.001');

    uibutton(tab2, ...
        'Position', [30, 470, 420, 42], ...
        'Text', 'Run NMF-RI', ...
        'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(b,e) run_nmf_ri());

    txtRunStatus = uilabel(tab2, ...
        'Position', [30, 400, 430, 58], ...
        'Text', '', ...
        'WordWrap', 'on', ...
        'FontColor', [0.25 0.25 0.25]);

    uibutton(tab2, ...
        'Position', [30, 352, 230, 36], ...
        'Text', 'Continue to Results >>', ...
        'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(b,e) go_to_results());

    axAcompare = uiaxes(tab2, 'Position', [510, 260, 735, 470]);
    title(axAcompare, 'A0 (dashed) vs estimated A (solid)');
    xlabel(axAcompare, 'Wavelength (nm)');
    ylabel(axAcompare, 'Sum-normalized intensity');

    %% ===================== Tab 3: Results =====================
    axSpec = uiaxes(tab3, 'Position', [30, 445, 590, 320]);
    title(axSpec, 'Resolved Spectra (A)');
    xlabel(axSpec, 'Wavelength (nm)');
    ylabel(axSpec, 'Sum-normalized intensity');

    axConc = uiaxes(tab3, 'Position', [650, 445, 590, 320]);
    title(axConc, 'Concentration Profile (H)');
    xlabel(axConc, 'Timepoint');
    ylabel(axConc, 'Intensity (a.u.)');

    axRecon = uiaxes(tab3, 'Position', [30, 105, 1210, 300]);
    title(axRecon, 'Observed vs Reconstructed (Timepoint 1)');
    xlabel(axRecon, 'Wavelength (nm)');
    ylabel(axRecon, 'Intensity');

    uilabel(tab3, ...
        'Position', [30, 65, 160, 24], ...
        'Text', 'Output file name');

    txtOutName = uieditfield(tab3, 'text', ...
        'Position', [180, 65, 280, 24], ...
        'Value', 'NMF_RI_output');

    uibutton(tab3, ...
        'Position', [480, 59, 220, 34], ...
        'Text', 'Save All Results', ...
        'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(b,e) download_results());

    txtOutStatus = uilabel(tab3, ...
        'Position', [720, 55, 510, 42], ...
        'Text', '', ...
        'FontColor', [0.25 0.25 0.25], ...
        'WordWrap', 'on');

    %% ===================== Nested callbacks =====================

    function browse_mixed()
        [file, path] = uigetfile({'*.txt', 'TXT files (*.txt)'}, ...
            'Select mixed raw TXT');
        if isequal(file, 0), return; end
        appState.mixedFile = fullfile(path, file);
        lblMixed.Text = file;
    end

    function browse_comp1()
        [file, path] = uigetfile({'*.csv', 'CSV files (*.csv)'}, ...
            'Select Component 1 emission spectrum CSV');
        if isequal(file, 0), return; end
        appState.comp1File = fullfile(path, file);
        lblComp1.Text = file;
    end

    function browse_comp2()
        [file, path] = uigetfile({'*.csv', 'CSV files (*.csv)'}, ...
            'Select Component 2 emission spectrum CSV');
        if isequal(file, 0), return; end
        appState.comp2File = fullfile(path, file);
        lblComp2.Text = file;
    end

    function load_data()
        fprintf('\n');
        fprintf('============================================================\n');
        fprintf('NMF-RI v2 : data loading started\n');
        fprintf('Start time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf('============================================================\n');

        if isempty(appState.mixedFile) || ...
                isempty(appState.comp1File) || ...
                isempty(appState.comp2File)
            fprintf('[ERROR] Mixed TXT and both component CSV files are required.\n\n');
            uialert(fig, 'Please select all three files.', 'Missing files');
            return;
        end

        wlMin = str2double(txtWLMin.Value);
        wlMax = str2double(txtWLMax.Value);

        if ~isfinite(wlMin) || ~isfinite(wlMax) || wlMin >= wlMax
            fprintf('[ERROR] Invalid wavelength range: Min=%s, Max=%s\n\n', ...
                txtWLMin.Value, txtWLMax.Value);
            uialert(fig, 'Enter a valid wavelength Min/Max.', 'Wavelength range error');
            return;
        end

        txtDataStatus.Text = 'Reading mixed TXT...';
        drawnow;
        tTotal = tic;

        try
            fprintf('[1/5] Input files\n');
            fprintf('      Mixed TXT      : %s\n', appState.mixedFile);
            fprintf('      Component 1 CSV: %s\n', appState.comp1File);
            fprintf('      Component 2 CSV: %s\n', appState.comp2File);
            fprintf('      Requested range: %.3f - %.3f nm\n\n', wlMin, wlMax);

            fprintf('[2/5] Parsing mixed TXT\n');
            tStep = tic;
            [wlAll, spectra, timeInfo] = readMixedTxtLogged(appState.mixedFile);
            fprintf('      Done: spectra [%d x %d]\n', size(spectra,1), size(spectra,2));
            fprintf('      Elapsed: %.2f s\n\n', toc(tStep));

            fprintf('[3/5] Cropping wavelengths and building Y\n');
            tStep = tic;
            wlMask = wlAll >= wlMin & wlAll <= wlMax;
            if ~any(wlMask)
                error('No wavelengths fall within the selected range.');
            end

            wlCrop = wlAll(wlMask);
            Y = spectra(:, wlMask)';
            [W, T] = size(Y);

            if W < 3
                error('Only %d wavelengths selected. Widen the range.', W);
            end

            fprintf('      Actual range: %.3f - %.3f nm\n', wlCrop(1), wlCrop(end));
            fprintf('      Y size      : %d x %d [W x T]\n', W, T);
            fprintf('      Elapsed     : %.2f s\n\n', toc(tStep));

            clear spectra;

            fprintf('[4/5] Reading component emission spectra\n');
            tStep = tic;
            aComp1 = loadEmissionSpectrum(appState.comp1File, wlCrop);
            aComp2 = loadEmissionSpectrum(appState.comp2File, wlCrop);
            fprintf('      Component 1 peak: %.3f nm\n', wlCrop(find(aComp1 == max(aComp1), 1, 'first')));
            fprintf('      Component 2 peak: %.3f nm\n', wlCrop(find(aComp2 == max(aComp2), 1, 'first')));
            fprintf('      Elapsed: %.2f s\n\n', toc(tStep));

            fprintf('[5/5] Building initial spectra A0\n');
            tStep = tic;
            A0 = [aComp1(:), aComp2(:)];
            A0(A0 <= 0) = 1e-4;
            A0 = A0 * diag(1 ./ (sum(A0, 1) + eps));

            appState.Y               = Y;
            appState.wavelengths     = wlCrop;
            appState.allWavelengths  = wlAll;
            appState.wlMask          = wlMask;
            appState.rawTimeInfo     = timeInfo;
            appState.A0              = A0;
            appState.A               = [];
            appState.H               = [];
            appState.nIter           = 0;
            appState.YPrepared       = [];

            fprintf('      A0 size: %d x %d\n', size(A0,1), size(A0,2));
            fprintf('      Elapsed: %.2f s\n\n', toc(tStep));

            names = {txtName1.Value, txtName2.Value};
            cla(axInit);
            hold(axInit, 'on');
            plot(axInit, wlCrop, A0(:,1), 'g-', 'LineWidth', 2, ...
                'DisplayName', [names{1}, ' (A0)']);
            plot(axInit, wlCrop, A0(:,2), 'r-', 'LineWidth', 2, ...
                'DisplayName', [names{2}, ' (A0)']);
            hold(axInit, 'off');
            legend(axInit, 'show', 'Location', 'best');
            title(axInit, 'Initial spectra A0');

            elapsed = toc(tTotal);
            txtDataStatus.Text = sprintf( ...
                'Done: Y [%d x %d], %.3f-%.3f nm, %.2f s', ...
                W, T, wlCrop(1), wlCrop(end), elapsed);

            fprintf('============================================================\n');
            fprintf('Data loading complete: %.2f s\n', elapsed);
            fprintf('============================================================\n\n');

        catch ME
            fprintf('\n[ERROR] Data loading failed\n');
            fprintf('Message: %s\n', ME.message);
            printErrorStack(ME);
            fprintf('\n');
            txtDataStatus.Text = ['Load error: ', ME.message];
            uialert(fig, sprintf('Data loading error:\n%s', ME.message), 'Load failed');
        end
    end

    function run_nmf_ri()
        fprintf('\n');
        fprintf('============================================================\n');
        fprintf('NMF-RI v2 : run started\n');
        fprintf('Start time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf('============================================================\n');

        if isempty(appState.Y) || isempty(appState.A0)
            fprintf('[ERROR] Load data in Tab 1 first.\n\n');
            uialert(fig, 'Please load data in Tab 1 first.', 'No data');
            return;
        end

        if exist('NMF_RI', 'file') ~= 2
            fprintf('[ERROR] NMF_RI.m not found.\n\n');
            uialert(fig, 'NMF_RI.m is not on the MATLAB path.', 'NMF_RI.m missing');
            return;
        end

        if exist('dataPreprocessing', 'file') ~= 2
            fprintf('[ERROR] dataPreprocessing.m not found.\n\n');
            uialert(fig, 'dataPreprocessing.m is not on the MATLAB path.', ...
                'dataPreprocessing.m missing');
            return;
        end

        theta1 = str2double(txtTheta1.Value);
        theta2 = str2double(txtTheta2.Value);
        alphaA = str2double(txtAlphaA.Value);
        alphaX = str2double(txtAlphaX.Value);

        params = [theta1, theta2, alphaA, alphaX];
        if any(~isfinite(params))
            fprintf('[ERROR] Non-numeric parameter value.\n\n');
            uialert(fig, 'All parameters must be numeric.', 'Parameter error');
            return;
        end

        Y  = appState.Y;
        A0 = appState.A0;

        fprintf('[Settings]\n');
        fprintf('      Y       : %d x %d [W x T]\n', size(Y,1), size(Y,2));
        fprintf('      A0      : %d x %d\n', size(A0,1), size(A0,2));
        fprintf('      theta1  : %.8f\n', theta1);
        fprintf('      theta2  : %.8f\n', theta2);
        fprintf('      alphaA  : %.8f\n', alphaA);
        fprintf('      alphaX  : %.8f\n\n', alphaX);

        d = uiprogressdlg(fig, ...
            'Title', 'Running NMF-RI', ...
            'Message', '[1/3] Preparing data...', ...
            'Indeterminate', 'on');

        tTotal = tic;

        try
            % ---------- [1/3] Data preparation ----------
            fprintf('[1/3] Data preparation\n');
            tStep = tic;

            [YPrepared, ~, ~] = dataPreprocessing(A0, Y);

            tPre = toc(tStep);
            totalCount = size(YPrepared, 2);

            fprintf('      Timepoints: %d\n', totalCount);
            fprintf('      Elapsed   : %.2f s\n\n', tPre);

            % ---------- [2/3] Spectral estimation ----------
            d.Message = '[2/3] Estimating spectra...';
            drawnow;

            fprintf('[2/3] NMF-RI optimization\n');
            fprintf('      No further output until this step completes.\n');
            fprintf('      Start time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
            tStep = tic;

            H0 = max(100 * eps, pinv(A0) * YPrepared);

            [ANmf, ~, nIter] = NMF_RI( ...
                YPrepared, A0, H0, theta1, theta2, alphaA, alphaX);

            tNmf = toc(tStep);
            fprintf('      End time  : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
            fprintf('      Iterations: %d\n', nIter);
            fprintf('      Elapsed   : %.2f s (%.2f min)\n\n', tNmf, tNmf/60);

            % ---------- [3/3] Full time-series H ----------
            d.Message = '[3/3] Computing full time-series H...';
            drawnow;

            fprintf('[3/3] Full time-series H\n');
            tStep = tic;
            HFull = max(1e6 * eps, ...
                pinv(ANmf' * ANmf) * ANmf' * YPrepared);
            tH = toc(tStep);

            fprintf('      H size : %d x %d\n', size(HFull,1), size(HFull,2));
            fprintf('      Minimum: %.6e\n', min(HFull(:)));
            fprintf('      Maximum: %.6f\n', max(HFull(:)));
            fprintf('      Elapsed: %.2f s\n\n', tH);

            % ---------- Spectral change from A0 ----------
            angDeg = zeros(1, nComp);
            for k = 1:nComp
                v0 = A0(:,k);  v1 = ANmf(:,k);
                denom = norm(v0) * norm(v1);
                if denom > 0
                    angDeg(k) = acosd(min(1, max(-1, dot(v0, v1) / denom)));
                end
            end
            fprintf('[Spectral change from A0]\n');
            for k = 1:nComp
                fprintf('      Component %d: %.3f deg\n', k, angDeg(k));
            end
            fprintf('\n');

            appState.A         = ANmf;
            appState.H         = HFull;
            appState.nIter     = nIter;
            appState.YPrepared = YPrepared;

            if isvalid(d), close(d); end

            totalElapsed = toc(tTotal);
            txtRunStatus.Text = sprintf( ...
                'Done: %d iterations, NMF %.2f min, total %.2f min', ...
                nIter, tNmf/60, totalElapsed/60);

            names = {txtName1.Value, txtName2.Value};
            cla(axAcompare);
            hold(axAcompare, 'on');
            plot(axAcompare, appState.wavelengths, A0(:,1), 'g--', ...
                'LineWidth', 1.2, 'DisplayName', [names{1}, ' A0']);
            plot(axAcompare, appState.wavelengths, ANmf(:,1), 'g-', ...
                'LineWidth', 2.0, 'DisplayName', [names{1}, ' result']);
            plot(axAcompare, appState.wavelengths, A0(:,2), 'r--', ...
                'LineWidth', 1.2, 'DisplayName', [names{2}, ' A0']);
            plot(axAcompare, appState.wavelengths, ANmf(:,2), 'r-', ...
                'LineWidth', 2.0, 'DisplayName', [names{2}, ' result']);
            hold(axAcompare, 'off');
            legend(axAcompare, 'show', 'Location', 'best');
            title(axAcompare, sprintf( ...
                ['A0 (dashed) vs estimated A (solid) - %d iterations, ' ...
                 'shift %.1f/%.1f deg'], nIter, angDeg(1), angDeg(2)));

            fprintf('============================================================\n');
            fprintf('NMF-RI run complete\n');
            fprintf('Preparation  : %.2f s\n', tPre);
            fprintf('Optimization : %.2f s\n', tNmf);
            fprintf('Full H       : %.2f s\n', tH);
            fprintf('Total        : %.2f s (%.2f min)\n', totalElapsed, totalElapsed/60);
            fprintf('============================================================\n\n');

        catch ME
            if exist('d', 'var') && isvalid(d), close(d); end
            fprintf('\n[ERROR] NMF-RI run failed\n');
            fprintf('Message: %s\n', ME.message);
            printErrorStack(ME);
            fprintf('\n');
            txtRunStatus.Text = ['NMF-RI error: ', ME.message];
            uialert(fig, sprintf('NMF-RI run error:\n%s', ME.message), ...
                'NMF-RI run error');
        end
    end

    function go_to_results()
        if isempty(appState.A) || isempty(appState.H)
            uialert(fig, 'Please run NMF-RI in Tab 2 first.', 'No results');
            return;
        end

        tabGroup.SelectedTab = tab3;
        try
            plot_results();
        catch ME
            uialert(fig, sprintf('Plotting error:\n%s', ME.message), 'Plotting error');
        end
    end

    function plot_results()
        A  = appState.A;
        H  = appState.H;
        wl = appState.wavelengths;

        if isempty(appState.YPrepared)
            YPlot = appState.Y - min(appState.Y(:));
        else
            YPlot = appState.YPrepared;
        end

        names = {txtName1.Value, txtName2.Value};
        T = size(H, 2);

        cla(axSpec);
        hold(axSpec, 'on');
        plot(axSpec, wl, A(:,1), 'g-', 'LineWidth', 2, 'DisplayName', names{1});
        plot(axSpec, wl, A(:,2), 'r-', 'LineWidth', 2, 'DisplayName', names{2});
        hold(axSpec, 'off');
        legend(axSpec, 'show', 'Location', 'best');
        title(axSpec, 'Resolved Spectra (A)');

        cla(axConc);
        hold(axConc, 'on');
        plot(axConc, 1:T, H(1,:), 'g-', 'LineWidth', 1, 'DisplayName', names{1});
        plot(axConc, 1:T, H(2,:), 'r-', 'LineWidth', 1, 'DisplayName', names{2});
        hold(axConc, 'off');
        legend(axConc, 'show', 'Location', 'best');
        title(axConc, 'Concentration Profile (H) - full time series');

        YRecon = A * H;
        cla(axRecon);
        hold(axRecon, 'on');
        plot(axRecon, wl, YPlot(:,1), 'k-', 'LineWidth', 1.5, ...
            'DisplayName', 'Observed / prepared (t=1)');
        plot(axRecon, wl, YRecon(:,1), 'b--', 'LineWidth', 1.5, ...
            'DisplayName', 'Reconstructed (t=1)');
        hold(axRecon, 'off');
        legend(axRecon, 'show', 'Location', 'best');
        title(axRecon, 'Observed vs Reconstructed (Timepoint 1)');
    end

    function download_results()
        if isempty(appState.A) || isempty(appState.H)
            uialert(fig, 'No results to save.', 'Cannot save');
            return;
        end

        defaultName = strtrim(txtOutName.Value);
        if isempty(defaultName)
            defaultName = 'NMF_RI_output';
        end

        [file, path] = uiputfile('*.csv', 'Save results', [defaultName, '.csv']);
        if isequal(file, 0)
            fprintf('[SAVE] Save cancelled by user.\n');
            return;
        end

        [~, baseName, ~] = fileparts(file);
        mainCsv = fullfile(path, [baseName, '.csv']);
        spectraCsv = fullfile(path, [baseName, '_spectra.csv']);

        names = {txtName1.Value, txtName2.Value};
        validNames = cell(1, nComp);
        for k = 1:nComp
            validNames{k} = matlab.lang.makeValidName(names{k});
        end

        fprintf('\n');
        fprintf('============================================================\n');
        fprintf('NMF-RI v2 : saving results\n');
        fprintf('Folder   : %s\n', path);
        fprintf('Base name: %s\n', baseName);
        fprintf('============================================================\n');

        tTotal = tic;

        try
            T = size(appState.H, 2);

            fprintf('[1/3] Saving time-series CSV\n');
            tStep = tic;
            timePoint = (0:T-1)';

            if ~isempty(appState.rawTimeInfo) && size(appState.rawTimeInfo,1) == T
                dateCol = appState.rawTimeInfo(:,1);
                timeCol = appState.rawTimeInfo(:,2);
                timestampCol = appState.rawTimeInfo(:,3);

                Tbl = table(timePoint, dateCol, timeCol, timestampCol, ...
                    appState.H(1,:)', appState.H(2,:)', ...
                    'VariableNames', ...
                    {'TimePoint', 'Date', 'Time', 'Timestamp', ...
                     validNames{1}, validNames{2}});
            else
                Tbl = table(timePoint, appState.H(1,:)', appState.H(2,:)', ...
                    'VariableNames', {'TimePoint', validNames{1}, validNames{2}});
            end

            writetable(Tbl, mainCsv);
            fprintf('      File : %s\n', mainCsv);
            fprintf('      Size : %d rows x %d columns\n', height(Tbl), width(Tbl));
            fprintf('      Elapsed: %.2f s\n\n', toc(tStep));

            fprintf('[2/3] Saving estimated spectra CSV\n');
            tStep = tic;
            SpecTbl = table(appState.wavelengths(:), ...
                appState.A(:,1), appState.A(:,2), ...
                'VariableNames', {'Wavelength_nm', validNames{1}, validNames{2}});
            writetable(SpecTbl, spectraCsv);
            fprintf('      File : %s\n', spectraCsv);
            fprintf('      Size : %d rows x %d columns\n', height(SpecTbl), width(SpecTbl));
            fprintf('      Elapsed: %.2f s\n\n', toc(tStep));

            fprintf('[3/3] Saving reconstructed component TXT files\n');
            compPaths = strings(nComp,1);
            for k = 1:nComp
                tStep = tic;
                compPaths(k) = string(fullfile(path, ...
                    sprintf('%s_%s.txt', baseName, validNames{k})));

                fprintf('      [%d/%d] Writing %s\n', k, nComp, compPaths(k));
                writeComponentTxtLogged(char(compPaths(k)), ...
                    appState.mixedFile, ...
                    appState.A(:,k), ...
                    appState.H(k,:), ...
                    appState.wlMask, ...
                    appState.rawTimeInfo, ...
                    validNames{k});
                fprintf('      [%d/%d] Done: %.2f s\n', k, nComp, toc(tStep));
            end

            elapsed = toc(tTotal);
            txtOutStatus.Text = sprintf( ...
                'Saved: 2 CSV files + %d component TXT files, %.2f s', ...
                nComp, elapsed);

            fprintf('============================================================\n');
            fprintf('Results saved\n');
            fprintf('1) %s\n', mainCsv);
            fprintf('2) %s\n', spectraCsv);
            for k = 1:nComp
                fprintf('%d) %s\n', k + 2, compPaths(k));
            end
            fprintf('Total save time: %.2f s\n', elapsed);
            fprintf('============================================================\n\n');

        catch ME
            fprintf('\n[ERROR] Saving results failed\n');
            fprintf('Message: %s\n', ME.message);
            printErrorStack(ME);
            fprintf('\n');
            txtOutStatus.Text = ['Save error: ', ME.message];
            uialert(fig, sprintf('Save error:\n%s', ME.message), 'Save error');
        end
    end
end

%% ===================== Local helper functions =====================

function em = loadEmissionSpectrum(csvFile, wlTarget)
    % Locate the emission column in a spectrum CSV and interpolate it onto
    % the target wavelength grid. The CSV may come from any source.
    % VariableNamingRule='preserve' avoids automatic header renaming warnings.

    T = readtable(csvFile, 'VariableNamingRule', 'preserve');

    if width(T) < 2
        error('Spectrum CSV needs at least two columns: %s', csvFile);
    end

    wlCol = convertColumnToDouble(T{:,1});
    colNames = string(T.Properties.VariableNames);

    emIdx = [];
    for i = 2:numel(colNames)
        name = lower(colNames(i));
        if contains(name, 'em') && ...
                ~contains(name, 'ex') && ...
                ~contains(name, '2p')
            emIdx = i;
            break;
        end
    end

    if isempty(emIdx)
        error('No emission column found in: %s', csvFile);
    end

    raw = convertColumnToDouble(T{:,emIdx});
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
    elseif iscell(raw)
        x = str2double(string(raw));
    else
        x = str2double(string(raw));
    end
    x = x(:);
end

function [wavelengths, spectra, timeInfo] = readMixedTxtLogged(filename)
    % ================================================================
    % Parsing strategy (unchanged from the original implementation):
    %   1) read the file line by line with readlines
    %   2) the first valid line after the marker holds the wavelengths
    %   3) split each data row on whitespace
    %   4) the last nPix values of each row form the spectrum
    % ================================================================

    totalTimer = tic;

    fprintf('      [READ 1/6] Reading file\n');
    tStep = tic;

    lines = readlines(filename);

    fprintf('                 Total lines: %d\n', numel(lines));
    fprintf('                 Elapsed    : %.2f s\n', toc(tStep));

    fprintf('      [READ 2/6] Locating Begin Spectral Data marker\n');
    tStep = tic;

    startIdx = find(contains(lines, '>>>>>Begin Spectral Data<<<<<'), ...
        1, 'first');

    if isempty(startIdx)
        error('Cannot find Begin Spectral Data marker.');
    end

    fprintf('                 Marker line: %d\n', startIdx);
    fprintf('                 Elapsed    : %.2f s\n', toc(tStep));

    fprintf('      [READ 3/6] Collecting valid data lines\n');
    tStep = tic;

    dataLines = lines(startIdx + 1:end);
    dataLines = dataLines(strlength(strtrim(dataLines)) > 0);

    if numel(dataLines) < 2
        error('Not enough spectral data after the marker.');
    end

    fprintf('                 Valid data lines: %d\n', numel(dataLines));
    fprintf('                 Elapsed         : %.2f s\n', toc(tStep));

    fprintf('      [READ 4/6] Parsing wavelength header\n');
    tStep = tic;

    wavelengthTokens = split(strtrim(dataLines(1)));
    wavelengths = str2double(wavelengthTokens);
    wavelengths = wavelengths(~isnan(wavelengths));
    wavelengths = wavelengths(:)';

    nPix = numel(wavelengths);

    if nPix == 0
        error('Could not parse the wavelength header.');
    end

    fprintf('                 Wavelengths: %d\n', nPix);
    fprintf('                 Range      : %.3f - %.3f nm\n', ...
        wavelengths(1), wavelengths(end));
    fprintf('                 Elapsed    : %.2f s\n', toc(tStep));

    fprintf('      [READ 5/6] Parsing spectral rows\n');
    tStep = tic;

    maxRows = numel(dataLines) - 1;

    % Preallocate to the maximum possible size
    spectra = nan(maxRows, nPix);
    dateList = strings(maxRows, 1);
    timeList = strings(maxRows, 1);
    timestampList = strings(maxRows, 1);

    validCount = 0;
    skippedCount = 0;
    reportStep = max(1, floor(maxRows / 20));  % report about every 5%

    for i = 2:numel(dataLines)
        row = split(strtrim(dataLines(i)));

        if numel(row) >= nPix + 3
            % The last nPix values of the row form the spectrum
            numericPart = str2double(row(end-nPix+1:end));

            if all(~isnan(numericPart))
                validCount = validCount + 1;

                spectra(validCount, :) = numericPart(:)';
                dateList(validCount) = string(row(1));
                timeList(validCount) = string(row(2));
                timestampList(validCount) = string(row(3));
            else
                skippedCount = skippedCount + 1;
            end
        else
            skippedCount = skippedCount + 1;
        end

        processed = i - 1;

        if mod(processed, reportStep) == 0 || processed == maxRows
            fprintf(['                 [%6.2f%%] %d / %d rows', ...
                ' | valid=%d | skipped=%d | %.2f s\n'], ...
                100 * processed / maxRows, ...
                processed, maxRows, validCount, skippedCount, toc(tStep));
            drawnow;
        end
    end

    if validCount == 0
        error('No valid spectral data rows were found.');
    end

    fprintf('      [READ 6/6] Trimming preallocated arrays\n');
    shrinkTimer = tic;

    spectra = spectra(1:validCount, :);
    dateList = dateList(1:validCount);
    timeList = timeList(1:validCount);
    timestampList = timestampList(1:validCount);

    timeInfo = [dateList, timeList, timestampList];

    fprintf('                 Spectra size: %d x %d\n', ...
        size(spectra, 1), size(spectra, 2));
    fprintf('                 Valid rows  : %d\n', validCount);
    fprintf('                 Skipped rows: %d\n', skippedCount);
    fprintf('                 Elapsed     : %.2f s\n', toc(shrinkTimer));
    fprintf('      [READ] Total parsing time: %.2f s\n', toc(totalTimer));
end

function writeComponentTxtLogged(outPath, sourceTxtPath, ACol, HRow, ...
        wlMask, timeInfo, componentName)
    % Write a reconstructed component TXT that keeps the original
    % wavelength axis; wavelengths outside the crop are written as zero.

    totalTimer = tic;

    fprintf('              [TXT 1/4] Reading source header and wavelengths\n');
    srcLines = readlines(sourceTxtPath);
    startIdx = find(contains(srcLines, '>>>>>Begin Spectral Data<<<<<'), 1, 'first');
    if isempty(startIdx)
        error('Begin Spectral Data marker not found in the source TXT.');
    end

    if startIdx + 1 > numel(srcLines)
        error('No wavelength header in the source TXT.');
    end

    wlHeaderLine = srcLines(startIdx + 1);
    wlTokens = split(strtrim(wlHeaderLine));
    wlAll = str2double(wlTokens);
    nPix = sum(isfinite(wlAll));

    if nPix ~= numel(wlMask)
        error('Source wavelength count (%d) differs from wlMask length (%d).', nPix, numel(wlMask));
    end

    fprintf('              [TXT 2/4] Computing component reconstruction\n');
    reconCrop = ACol(:) * HRow(:)';
    nTime = size(reconCrop, 2);
    reconFull = zeros(nPix, nTime);
    reconFull(wlMask, :) = reconCrop;

    fprintf('              [TXT 3/4] Opening output file\n');
    fid = fopen(outPath, 'w');
    if fid < 0
        error('Cannot write file: %s', outPath);
    end
    cleanupObj = onCleanup(@() fcloseIfOpen(fid));

    for i = 1:startIdx
        fprintf(fid, '%s\n', srcLines(i));
    end
    fprintf(fid, '# Reconstructed component: %s (A(:,k) * H(k,:))\n', componentName);
    fprintf(fid, '%s\n', wlHeaderLine);

    fprintf('              [TXT 4/4] Writing %d timepoints\n', nTime);
    hasTimeInfo = ~isempty(timeInfo) && size(timeInfo,1) == nTime;
    reportStep = max(1, floor(nTime/10));

    % Write row by row; the numeric spectrum is formatted in one sprintf call.
    for t = 1:nTime
        if hasTimeInfo
            prefix = sprintf('%s\t%s\t%s', ...
                char(timeInfo(t,1)), char(timeInfo(t,2)), char(timeInfo(t,3)));
        else
            prefix = sprintf('0\t0\t%d', t);
        end

        spectrumText = sprintf('\t%.6f', reconFull(:,t));
        fprintf(fid, '%s%s\n', prefix, spectrumText);

        if mod(t, reportStep) == 0 || t == nTime
            fprintf('                      %6.2f%% (%d/%d) | %.2f s\n', ...
                100*t/nTime, t, nTime, toc(totalTimer));
            drawnow;
        end
    end

    fclose(fid);
    fid = -1;
    clear cleanupObj;

    fprintf('              TXT done: %.2f s\n', toc(totalTimer));
end

function fcloseIfOpen(fid)
    if isnumeric(fid) && isscalar(fid) && fid > 0
        try
            fclose(fid);
        catch
        end
    end
end

function printErrorStack(ME)
    for i = 1:min(5, numel(ME.stack))
        fprintf('        %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end
