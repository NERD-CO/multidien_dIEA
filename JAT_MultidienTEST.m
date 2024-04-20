

pcNAME = getenv('COMPUTERNAME');

switch pcNAME
    case 'DESKTOP-EGFQKAI'
        mainLOC = 'C:\Users\johna\OneDrive\Documents\Github\multidien_dIEA';

    otherwise

end

cd(mainLOC)

%%
fooof_fileN = '10162654_fooof.mat';
% load(fooof_fileN)

%%
histo_fileN = 'UColorado_AAB_10162654_Histogram_Hourly.csv';

%%
data_fooof = load(fooof_fileN);
data_histogram = readtable(histo_fileN);

%%

histo_days_since_implant = days(data_histogram.RegionStartTime - data_fooof.metadata.Implant);
% don't worry that some values are negative - these are pre-implant
histo_detections_preZ = data_histogram.EpisodeStarts;

plot(histo_days_since_implant,histo_detections_preZ) %%%
xlabel('days since implant')
ylabel('detections')

%%

epoch_edges = data_fooof.epoch_days;
epoch_edges = vertcat(epoch_edges,min(histo_days_since_implant),max(histo_days_since_implant)+1);
epoch_edges = sort(epoch_edges);

%% z-score histogram counts within each epoch

histo_detections_zScored = [];
for iEpoch = 1:(length(epoch_edges)-1)
    idx = histo_days_since_implant >= epoch_edges(iEpoch) & ...
        histo_days_since_implant < epoch_edges(iEpoch+1);
    preZ = histo_detections_preZ(idx);

    if all(isnan(preZ))
        continue
    elseif any(isnan(preZ))
        preZ = preZ(~isnan(preZ));
    end

    postZ = (preZ - mean(preZ,"omitnan")) / std(preZ,"omitnan");
    
    % can't use "zscore" function because of NaN values
    histo_detections_zScored = vertcat(histo_detections_zScored,postZ);
end

%%

hdzS = smoothdata(histo_detections_zScored,'gaussian',30);

close all
% disp(any(isnan(histo_detections_zScored)))
plot(hdzS)

%%
close all
[pxx,w] = periodogram(histo_detections_zScored);
pow2plot = pow2db(pxx);
pow2plotS = smoothdata(pow2plot,'sgolay',120);

% remove 0.05 
ind2keep = w >= 0.1;
powDATA = pow2plotS(ind2keep);
wDATA = w(ind2keep);


plot(w,pow2plotS)
xlim([0.1 3])

meanPOW = mean(pow2plotS);
stdPOW = std(pow2plotS);
thrFac = 1.25;
thrPOW = meanPOW + (stdPOW * thrFac);

yline(thrPOW,'-k','threshold')

%%
close all
coefficients = polyfit(wDATA, powDATA, 2);
y_fit = polyval(coefficients, wDATA);
% x_fit = linspace(min(wDATA), max(wDATA), 100); % Generate points for a smooth plot
% y_fit = polyval(coefficients, x_fit); % Evaluate polynomial

residuals = powDATA - y_fit;

meanPOW = mean(residuals);
stdPOW = std(residuals);
thrFac = 1.25;
thrPOW = meanPOW + (stdPOW * thrFac);




plot(wDATA,residuals)

yline(thrPOW,'-k','threshold')
%%
close all
detS = detrend(powDATA)
plot(detS)




%%
close all
periodogram(histo_detections_zScored);
%1) hourly histogram event counts
%2) z-score between programming changes
%3) periodogram
%4) find significant periodogram peaks
%5) inverse wavelet transform to reconstruct each peak (effectively bandpassing at significant frequencies)
%6) find period with highest phase locking value
%7) use hilbert transform to calculate phase angle of each scheduled recording
%8) classify scheduled recordings as peak/trough/rising/falling phase
% The paper from Penn has some specifics but is missing other details and has things I would do differently.