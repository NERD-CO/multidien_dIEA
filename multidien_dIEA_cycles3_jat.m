clear
close all
%%

%% INDEX FOOOF DATA
cd('C:\Users\Admin\Documents\Github\multidien_dIEA\fooofFILES')
filenames_fooof = dir('*_fooof.mat');

%% INDEX HOURLY HISTOGRAMS
addpath('Histograms') % add folder with histogram files to path
filenames_histograms = dir('C:\Users\Admin\Documents\Github\multidien_dIEA\Histograms/*_Histogram_Hourly.csv');

for iHistogram = 1:length(filenames_histograms) % label histogram file list with PDMS IDs
    filenames_histograms_ID(iHistogram,1) = str2double(extract(filenames_histograms(iHistogram).name,digitsPattern));
end
clear iHistogram

load('aperiodic_periods.mat');

%% PER-PARTICIPANT ANALYSIS

for iParticipant = 1:length(filenames_foof)
    % iParticipant = 1; %%%% temp
    data_fooof = load(filenames_fooof(iParticipant).name);
    % data_fooof = load(fooof_fileN);
    histogram_row = find(filenames_histograms_ID==data_fooof.ID);
    data_histogram = readtable(strcat('Histograms/',filenames_histograms(histogram_row).name));
    % data_histogram = readtable(histo_fileN);
    % clear histogram_row
    histo_days_since_implant = days(data_histogram.RegionStartTime - data_fooof.metadata.Implant);
    % don't worry that some values are negative - these are pre-implant
    histo_detections_preZ = data_histogram.EpisodeStarts;

    plot(histo_days_since_implant,histo_detections_preZ) %%%
    xlabel('days since implant')
    ylabel('detections')

    epoch_edges = data_fooof.epoch_days;
    epoch_edges = vertcat(epoch_edges,min(histo_days_since_implant),max(histo_days_since_implant)+1);
    epoch_edges = sort(epoch_edges);


    histo_detections_zScored = [];
    for iEpoch = 1:(length(epoch_edges)-1)
        idx = histo_days_since_implant >= epoch_edges(iEpoch) & ...
            histo_days_since_implant < epoch_edges(iEpoch+1);

        preZ = histo_detections_preZ(idx);
        postZ = (preZ - mean(preZ,"omitnan")) / std(preZ,"omitnan");
        % can't use "zscore" function because of NaN values
        histo_detections_zScored = vertcat(histo_detections_zScored,postZ);

        clear idx preZ postZ
    end
    % clear iEpoch
    % clear epoch_edges histo_detections_preZ

    iTime = histo_days_since_implant;
    iValu = histo_detections_zScored;

    % get rid of data before implant
    idx = histo_days_since_implant > 0;
    iTime = iTime(idx);
    iValu = iValu(idx);
    clear idx

    % get rid of duplicate times
    idx(1) = true;
    for iBin = 2:length(iTime)
        idx(iBin) = iTime(iBin)>iTime(iBin-1);
    end
    clear iBin
    iTime = iTime(idx);
    iValu = iValu(idx);
    clear idx

    peakIndices = periodogram_peaks.periods_days{periodogram_peaks.ID==data_fooof.ID};

    inverseWaveForm_PLV_fun(iValu , peakIndices)

end

