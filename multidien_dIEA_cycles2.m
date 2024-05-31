clear
close all

%% INDEX FOOOF DATA
filenames_fooof = dir('*_fooof.mat');

%% INDEX HOURLY HISTOGRAMS
addpath('Histograms') % add folder with histogram files to path
filenames_histograms = dir('Histograms/*_Histogram_Hourly.csv');

for iHistogram = 1:length(filenames_histograms) % label histogram file list with PDMS IDs
    filenames_histograms_ID(iHistogram,1) = str2double(extract(filenames_histograms(iHistogram).name,digitsPattern));
end
clear iHistogram

%% PER-PARTICIPANT ANALYSIS

% for iParticipant = 1:length(filenames_foof)
iParticipant = 1; %%%% temp
    data_fooof = load(filenames_fooof(iParticipant).name);
    histogram_row = find(filenames_histograms_ID==data_fooof.ID);
    data_histogram = readtable(strcat('Histograms/',filenames_histograms(histogram_row).name));
    clear histogram_row
    histo_days_since_implant = days(data_histogram.RegionStartTime - data_fooof.metadata.Implant);
        % don't worry that some values are negative - these are pre-implant
    histo_detections_preZ = data_histogram.EpisodeStarts;

    plot(histo_days_since_implant,histo_detections_preZ) %%%
        xlabel('days since implant')
        ylabel('detections')

    epoch_edges = data_fooof.epoch_days;
    epoch_edges = vertcat(epoch_edges,min(histo_days_since_implant),max(histo_days_since_implant)+1);
    epoch_edges = sort(epoch_edges);

    %% z-score histogram counts within each epoch

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
    clear iEpoch

    clear epoch_edges histo_detections_preZ

    %%

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

    period_vec = 3:60;
    fvec = 1./period_vec;

    [pxx,f] = plomb(iValu,iTime,fvec);
    figure
    plot(period_vec,pxx)
    xlabel('Period (days)')
    ylabel('PSD')
    title('Periodogram with prespecified periods')


    clear pxx f

    Pfa = [50 10 1 0.01]/100;
    Pd = 1-Pfa;

    [pxx,f,pth] = plomb(iValu,iTime,'normalized','Pd',Pd);
    figure
    hold on
    plot(1./f,pxx)
    plot([3 60],[pth(1) pth(1)])
    plot([3 60],[pth(2) pth(2)])
    plot([3 60],[pth(3) pth(3)])
    plot([3 60],[pth(4) pth(4)])
    xlim([0 60])
    xlabel('Period (days)')
    ylabel('PSD')
    legend({'periodogram' 'P(FA)=0.5' 'P(FA)=0.1' 'P(FA)=0.01' 'P(FA)=0.0.0001'})
    title({'Periodogram with pdvec';'with False Alarm Probabilities'})

    % 
    % plot(f,pxx,f,pth*ones(size(f')))
    % xlabel('f')
    % text(0.3*[1 1 1 1],pth-.5,[repmat('P_{fa} = ',[4 1]) num2str(Pfa')])



    % %% say for the sake of argument, period is 27 d
    % % limit to continuous period
    % idx = iTime>250 & iTime < 370;
    % xData = iTime(idx);
    % yData = iValu(idx);
    % 
    % yData_HP = highpass(yData,1/28,24);
    % yData_HP_LP = lowpass(yData_HP,1/26,24);
    % 
    % figure
    % hold on
    % plot(xData,yData)
    % plot(xData,yData_HP_LP)


% % % for ii=1:(length(histo_days_since_implant)-1)
% % % gap(ii)=histo_days_since_implant(ii+1)-histo_days_since_implant(ii);
% % % end
% % % idx = gap>0.05;
% % % sum(idx)


    %%

    % clear data_fooof data_histo
% end
clear iParticipant

clear filenames_fooof filenames_histograms filenames_histograms_ID