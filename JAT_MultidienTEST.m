

pcNAME = getenv('COMPUTERNAME');

switch pcNAME
    case 'DESKTOP-EGFQKAI'
        mainLOC = 'C:\Users\johna\OneDrive\Documents\Github\multidien_dIEA';
    case 'DESKTOP-FAGRV5G'
        mainLOC = 'C:\Users\Admin\Documents\Github\multidien_dIEA';

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
skipDATA = floor(diff(histo_days_since_implant));
startIND = find(skipDATA) + 1;
stopIND = startIND + skipDATA(find(skipDATA));

part1 = histo_detections_preZ(1:startIND(1)-1);
part1a = nan(stopIND(1)-startIND(1),1)
part2 = histo_detections_preZ(stopIND(1)+1:startIND(2)-1)
part2a = nan(stopIND(2)-startIND(2),1)
part3 = histo_detections_preZ(stopIND(2)+1:startIND(3)-1)
part3a = nan(stopIND(3)-startIND(3),1)
part4 = histo_detections_preZ(stopIND(3)+1:end)

allParts = [part1 ; part1a ; part2 ; part2a ; part3 ; part3a; part4]


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

hdzS = smoothdata(histo_detections_zScored,'gaussian',50);
close
[u1,d] = envelope(hdzS,50,'rms');

[u2,d] = envelope(hdzS,10,'peak');

figure;
% disp(any(isnan(histo_detections_zScored)))
plot(hdzS)
hold on
plot(u1)
plot(u2)

yline(mean(u2)-(std(u2)*0.3))

peakENV1 = u2;
peakENV1(u2 <  mean(u2)-(std(u2)*0.3)) = 0;
plot(peakENV1)
legend('Raw','RMS','Peak')
%%
plot(peakENV1)

N = numel(peakENV1);   % Number of data points
fftResult = fft(peakENV1);   % Compute FFT

% Calculate the two-sided spectrum and then the single-sided spectrum
P2 = abs(fftResult/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

Fs = 0.24;  % Sampling frequency (times per day)
f = Fs*(0:(N/2))/N;  % Frequency vector

% [peakValues, peakLocations] = findpeaks(P1, 'MinPeakHeight', threshold);
[peakValues, peakLocations] = findpeaks(P1, 'MinPeakHeight', 0.3);
peakFrequencies = f(peakLocations);  % Frequencies corresponding to the found peaks


figure;
plot(f, P1);
title('Single-Sided Amplitude Spectrum of Data');
xlabel('Frequency (cycles per day)');
ylabel('|P1(f)|');
hold on;
plot(peakFrequencies, peakValues, 'r*', 'MarkerSize', 10); % Mark peaks
hold off;


%%
[pxx,f] = periodogram(peakENV1,[],[],1);
pow2plot = pow2db(pxx);
% plot(f,10*log10(pxx))
plot(f(5:end),pow2plot(5:end))
pow2plotS = smoothdata(pow2plot,'sgolay',120)

xlabel('Cycles/Day')
ylabel('dB / (Cycles/Year)')

%%
histo_days_since_implant2 = histo_days_since_implant(8:end)
allParts2 = allParts(5:end)
histo_days_since_implant2 = days(histo_days_since_implant2)

oneday = seconds(days(1));

[pxx,f] = plomb(allParts2,histo_days_since_implant2,0.0417,10,'normalized');

f = f*oneday;

[pmax,lmax] = max(pxx);
f0 = f(lmax);

plot(f,pxx,f0,pmax,'o')
xlabel('Frequency (day^{-1})')

%%
close all
[pxx,w] = periodogram(hdzS);
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