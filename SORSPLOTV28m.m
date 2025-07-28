%% === Enhanced SORS Analysis Script - FINAL CORRECTED VERSION ===
% This script correctly analyzes and visualizes SORS datasets where spectra from
% multiple offsets are concatenated into a single matrix. It includes robust
% multi-file handling, user-driven dataset selection, and detailed plotting
% with clear annotations and labeling for mixture analysis.

clear; close all; clc;

%% --- USER SETTINGS ---
selectedOffset = 1;          % Which detector to analyze (1-6 corresponds to 0-5 mm)
showAllOffsetsSeparate = true; % Set to true to plot all offsets in subplots
zoom_peak = 'mtbe';             % Options: 'anisole', 'mtbe', 'pcm', or 'none'
zoom_range = 80;                % Zoom range in cm⁻¹ around the peak

annotatePeaks = true;        % Show peak annotations on plots
showZoomBox = true;          % Show a box around the zoom region on the main plot
showZoomedPlot = true;       % Show a separate, zoomed-in plot

applySNV = true;               % Apply Standard Normal Variate correction
savePlots = true;               % Save plots to files
saveCSV = true;                 % Save peak analysis results to a CSV file
outputFolder = 'SORS_Analysis_Results_filtrate&mixture'; % Name for the output folder

% Peak definitions (cm⁻¹) for annotation and analysis
anisole_peak = 1002;
mtbe_peak = 730;
pcm_peaks = [860, 1600];
trackPeaks = [anisole_peak, mtbe_peak, pcm_peaks];
peakSearchWindow = 80;

%% --- 1. SETUP AND LOAD DATA ---
% Create the output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Let the user select multiple .mat files
[fileList, dataPath] = uigetfile('*.mat', 'Select your .mat files...', 'MultiSelect', 'on');
if ischar(fileList), fileList = {fileList}; end
if isempty(fileList) || (iscell(fileList) && isempty(fileList{1}))
    error('No files selected. Analysis cancelled.');
end
nFiles = numel(fileList);

% Initialize storage for all loaded datasets
SpectrumSORS_all = cell(1, nFiles);
Ramanshift_all = cell(1, nFiles);
indexStartEnd_all = cell(1, nFiles);
dataLabels = cell(1, nFiles);
dataTags = cell(1, nFiles);

fprintf('=== LOADING DATASETS ===\n');
for i = 1:nFiles
    filePath = fullfile(dataPath, fileList{i});
    fprintf('Loading file %d/%d: %s\n', i, nFiles, fileList{i});
    data = load(filePath);
    
    % Validate that the necessary variables exist in the file
    if ~isfield(data, 'SpectrumSORS') || ~isfield(data, 'Ramanshift') || ~isfield(data, 'indexStartEnd')
        error('File %s is missing one or more required variables (SpectrumSORS, Ramanshift, indexStartEnd).', fileList{i});
    end
    
    SpectrumSORS_all{i} = data.SpectrumSORS;
    Ramanshift_all{i} = data.Ramanshift;
    indexStartEnd_all{i} = data.indexStartEnd;
    dataLabels{i} = erase(fileList{i}, '.mat');
    
    % Ask user to classify the dataset type (e.g., Mixture, Wet Cake)
    choice = questdlg(sprintf('What type of data is in "%s"?', fileList{i}), ...
        'Dataset Classification', 'Mixture', 'Wet Cake', 'Filtrate', 'Mixture');
    if isempty(choice), choice = 'Unknown'; end
    dataTags{i} = strrep(choice, ' ', '_');
end

% If more than one file was loaded, ask the user to select the primary one for analysis
if nFiles > 1
    selectionList = cellfun(@(label, tag) sprintf('%s (%s)', label, strrep(tag, '_', ' ')), ...
        dataLabels, dataTags, 'UniformOutput', false);
    [selectedDatasetIdx, isOK] = listdlg('PromptString', 'Select primary dataset for detailed analysis:', ...
        'SelectionMode', 'single', 'ListString', selectionList, 'ListSize', [400, 300]);
    if ~isOK, error('No dataset selected.'); end
else
    selectedDatasetIdx = 1;
end

% Set the primary dataset for the main analysis
SpectrumSORS = SpectrumSORS_all{selectedDatasetIdx};
Ramanshift = Ramanshift_all{selectedDatasetIdx};
indexStartEnd = indexStartEnd_all{selectedDatasetIdx};
currentDataLabel = dataLabels{selectedDatasetIdx};

fprintf('\nPrimary dataset selected: %s\n', currentDataLabel);

%% --- 2. ANALYZE DATA STRUCTURE AND APPLY CORRECTIONS ---
numSamples = size(SpectrumSORS, 1);
numOffsets = size(indexStartEnd, 1);

colorMap = distinguishable_colors(numSamples); % Define colors for each sample

function colors = distinguishable_colors(n_colors)
    % Generate n_colors visually distinct colors using Lab color space
    colors = zeros(n_colors, 3);
    cmap = lines(max(7,n_colors));  % Start with MATLAB's lines
    for i = 1:n_colors
        colors(i,:) = cmap(mod(i-1, size(cmap,1)) + 1, :);
    end
end
% Define the labels for the 5 mixture samples (Anisole:MTBE ratio)
sampleTitles = {'100% Anisole: 0% MTBE', '75% Anisole: 25% MTBE', '50% Anisole: 50% MTBE', '25 Anisole: 75% MTBE', '0% Anisole: 100 MTBE'};
if numSamples ~= 5
    % If not 5 samples, create generic labels
    sampleTitles = arrayfun(@(x) sprintf('Sample %d', x), 1:numSamples, 'UniformOutput', false);
    fprintf('Warning: Found %d samples, expected 5. Using generic labels.\n', numSamples);
end

offsetTitles = arrayfun(@(x) sprintf('%d mm', x - 1), 1:numOffsets, 'UniformOutput', false);

% Apply SNV Correction if enabled
if applySNV
    fprintf('\nApplying SNV Correction...\n');
    for s = 1:numSamples
        for o = 1:numOffsets
            idx_start = indexStartEnd(o, 1);
            idx_end = indexStartEnd(o, 2);
            spectrum_seg = SpectrumSORS(s, idx_start:idx_end);
            if std(spectrum_seg) > 0
                SpectrumSORS(s, idx_start:idx_end) = (spectrum_seg - mean(spectrum_seg)) / std(spectrum_seg);
            end
        end
    end
end

%% --- 3. MAIN PLOTTING (CORRECTED INDEXING) ---
fprintf('\nCreating Main Plots...\n');
mainFig = figure('Name', sprintf('SORS Spectra - %s', currentDataLabel), 'Position', [100, 100, 1200, 800]);

% Determine zoom range
zoom_center = [];
if strcmpi(zoom_peak, 'anisole'), zoom_center = anisole_peak; end
if strcmpi(zoom_peak, 'mtbe'), zoom_center = mtbe_peak; end
if strcmpi(zoom_peak, 'pcm'), zoom_center = pcm_peaks(1); end
if ~isempty(zoom_center)
    zoom_x_min = zoom_center - zoom_range;
    zoom_x_max = zoom_center + zoom_range;
end

if showAllOffsetsSeparate
    % Plot all detectors in separate subplots
    for o = 1:numOffsets
        subplot(numOffsets, 1, o);
        hold on;
        idx_start = indexStartEnd(o, 1);
        idx_end = indexStartEnd(o, 2);
        
        for s = 1:numSamples
            % CORRECTED PLOTTING: Plot segment against the FULL Ramanshift axis
            plot(Ramanshift, SpectrumSORS(s, idx_start:idx_end), 'LineWidth', 1.2, 'Color', colorMap(s,:));
        end
        
        title(sprintf('Detector %s - %s', offsetTitles{o}, currentDataLabel));
        common_plot_settings(Ramanshift, trackPeaks, annotatePeaks, zoom_center, zoom_x_min, zoom_x_max, showZoomBox);
        if o == 1, legend(sampleTitles, 'Location', 'best'); end
    end
else
    % Plot only the selected detector
    hold on;
    idx_start = indexStartEnd(selectedOffset, 1);
    idx_end = indexStartEnd(selectedOffset, 2);
    
    for s = 1:numSamples
        % CORRECTED PLOTTING
        plot(Ramanshift, SpectrumSORS(s, idx_start:idx_end), 'LineWidth', 1.5, 'Color', colorMap(s,:), 'DisplayName', sampleTitles{s});
    end
    
    title(sprintf('All Samples at Detector %s - %s', offsetTitles{selectedOffset}, currentDataLabel));
    common_plot_settings(Ramanshift, trackPeaks, annotatePeaks, zoom_center, zoom_x_min, zoom_x_max, showZoomBox);
    legend('Location', 'best');
end

if savePlots, saveas(mainFig, fullfile(outputFolder, sprintf('%s_Main_Spectra.png', currentDataLabel))); end


%% --- 4. ZOOMED PLOT (CORRECTED INDEXING) ---
if showZoomedPlot && ~isempty(zoom_center)
    fprintf('Creating Zoomed Plot...\n');
    zoomFig = figure('Name', sprintf('Zoomed View: %s', zoom_peak), 'Position', [200, 200, 800, 600]);
    hold on;
    
    idx_start = indexStartEnd(selectedOffset, 1);
    idx_end = indexStartEnd(selectedOffset, 2);
    
    for s = 1:numSamples
        % CORRECTED PLOTTING
        plot(Ramanshift, SpectrumSORS(s, idx_start:idx_end), 'LineWidth', 2, 'Color', colorMap(s,:), 'DisplayName', sampleTitles{s});
    end
    
    title(sprintf('Zoomed View: %.0f cm^{-1} - Detector %s', zoom_center, offsetTitles{selectedOffset}));
    common_plot_settings(Ramanshift, [], false, [], [], [], false); % No extra annotations on zoom
    xlim([zoom_x_min zoom_x_max]); % Apply zoom
    ylims = ylim;
    line([zoom_center zoom_center], ylims, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
    legend('Location', 'best');
    if savePlots, saveas(zoomFig, fullfile(outputFolder, sprintf('%s_Zoomed_Plot_%s.png', currentDataLabel, zoom_peak))); end
end

%% --- 5. PEAK ANALYSIS AND CSV EXPORT (CORRECTED INDEXING) ---
if saveCSV
    fprintf('\nPerforming Peak Analysis for CSV Export...\n');
    results = {};
    for p = 1:length(trackPeaks)
        peakVal = trackPeaks(p);
        for o = 1:numOffsets
            idx_start = indexStartEnd(o, 1);
            idx_end = indexStartEnd(o, 2);
            for s = 1:numSamples
                % CORRECTED ANALYSIS: Use full Ramanshift vector with the spectral segment
                [pos, inten, fwhm, area] = find_peak_metrics(Ramanshift, SpectrumSORS(s, idx_start:idx_end), peakVal, peakSearchWindow);
                results(end+1,:) = {sampleTitles{s}, offsetTitles{o}, peakVal, pos, inten, fwhm, area};
            end
        end
    end
    
    T = cell2table(results, 'VariableNames', {'Sample_Ratio','Detector','Expected_Peak','Found_Peak','Intensity','FWHM','Area'});
    csvPath = fullfile(outputFolder, sprintf('%s_peak_analysis.csv', currentDataLabel));
    writetable(T, csvPath);
    fprintf('Peak analysis results exported to: %s\n', csvPath);
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');


%% === HELPER FUNCTIONS ===

function common_plot_settings(x_axis, peaks_to_annotate, show_annotations, zoom_c, zoom_min, zoom_max, show_zoom_box)
    % Standard settings for all plots
    xlabel('Raman Shift (cm^{-1})');
    ylabel('Intensity (a.u.)');
    xlim([min(x_axis) max(x_axis)]);
    grid on; box on;
    
    % Add peak annotations if requested
    if show_annotations && ~isempty(peaks_to_annotate)
        ylims = ylim;
        peak_names = {'Anisole', 'MTBE', 'PCM-860', 'PCM-1600'};
        for i = 1:length(peaks_to_annotate)
            peak = peaks_to_annotate(i);
            line([peak peak], ylims, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
            text(peak, ylims(2)*0.9, peak_names{i}, 'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
        end
    end
    
    % Add zoom box if requested
    if show_zoom_box && ~isempty(zoom_c)
        ylims = ylim;
        patch([zoom_min zoom_max zoom_max zoom_min], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
              'red', 'EdgeColor', 'red', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.5);
    end
end

function [peak_pos, peak_int, width, area] = find_peak_metrics(x, y, center, window)
    % Simple function to find peak metrics within a search window
    idx = (x >= center - window/2) & (x <= center + window/2);
    if ~any(idx)
        [peak_pos, peak_int, width, area] = deal(NaN);
        return;
    end
    
    x_sub = x(idx);
    y_sub = y(idx);
    
    [peak_int, rel_idx] = max(y_sub);
    peak_pos = x_sub(rel_idx);
    
    half_max = peak_int / 2;
    left_idx = find(y_sub(1:rel_idx) <= half_max, 1, 'last');
    if isempty(left_idx), left_idx = 1; end
    
    right_idx_temp = find(y_sub(rel_idx:end) <= half_max, 1, 'first');
    if isempty(right_idx_temp), right_idx_temp = length(y_sub) - rel_idx + 1; end
    right_idx = rel_idx + right_idx_temp - 1;
    
    width = x_sub(right_idx) - x_sub(left_idx);
    area = trapz(x_sub(left_idx:right_idx), y_sub(left_idx:right_idx));
end