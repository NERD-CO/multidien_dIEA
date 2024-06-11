% Sample signal (replace this with your actual signal)
% signal = randn(1, 1000); 
signal = iValu;
% fs = 1000; % 
fs = 250;

% Perform Continuous Wavelet Transform (CWT)
signal = signal(~isnan(signal));
[cfs, frequencies] = cwt(signal, 'amor', fs);

% cfs = pxx;
% frequencies = f;

% Assume peakIndices is a vector containing the indices of the peaks in the original signal
% Convert these indices to the corresponding frequency indices in the CWT
% peakIndices = [50, 100, 150]; % Replace with your actual peak indices
peakIndices = [20, 35];

% Find the corresponding scales/frequencies for the peak indices
% Note: this is a simplified way; in practice, you may need to identify the correct scales more accurately
% Here we assume the peakIndices correspond directly to frequencies
numFrequencies = size(cfs, 1); % Number of frequency bins in the CWT
% freqIndices = round(peakIndices * numFrequencies / length(signal));
freqIndices = peakIndices;

% Initialize a matrix to store the isolated peak components
isolatedPeaks = zeros(size(cfs));

% Reconstruct the signal from each isolated peak component using ICWT
reconstructedPeaks = zeros(length(freqIndices), length(signal));

for i = 1:length(freqIndices)
    freqIndex = freqIndices(i);
    if freqIndex > 0 && freqIndex <= numFrequencies
        % Create a CWT matrix with only the current peak component
        peakCWT = zeros(size(cfs));
        peakCWT(freqIndex, :) = cfs(freqIndex, :);
        
        % Perform Inverse Continuous Wavelet Transform (ICWT)
        reconstructedSignal = icwt(peakCWT, 'amor');
        
        % Store the reconstructed signal
        reconstructedPeaks(i, :) = reconstructedSignal;
    end
end

% Plot the original signal and the reconstructed signals
%%
close all

figure;
subplot(3, 1, 1);
plot(signal);
title('Original Signal');

for i = 1:length(freqIndices)
    subplot(3, 1, i + 1);
    % plot(reconstructedPeaks(i, :));

    analyticSignal = envelope(reconstructedPeaks(1, :));
    sm_analyticSignal = smoothdata(analyticSignal,'sgolay',250);
    % instantaneousPhase = angle(analyticSignal);
    plot(sm_analyticSignal)

    title(['Reconstructed Signal from Peak ', num2str(peakIndices(i))]);
end


%%

% Select one of the reconstructed signals for PLV analysis
% reconstructedSignal = reconstructedPeaks(1, :); % Choose the desired signal
reconstructedSignal = sm_analyticSignal;

% Compute the phase of the reconstructed signal using the Hilbert transform
analyticSignal = hilbert(reconstructedSignal);
instantaneousPhase = angle(analyticSignal);

% Define the period length (e.g., 100 samples, adjust as needed)
% periodLength = 100;
periodLength = 31;

% Compute the PLV for each period
numPeriods = floor(length(instantaneousPhase) / periodLength);
plvValues = zeros(1, numPeriods);

for k = 1:numPeriods
    % Extract the phase values for the current period
    phaseValues = instantaneousPhase((k-1)*periodLength + 1 : k*periodLength);
    
    % Calculate the phase locking value (PLV)
    plvValues(k) = abs(mean(exp(1i * phaseValues)));
end

% Identify the period with the highest PLV
[maxPLV, maxPLVIndex] = max(plvValues);

% Output the period with the highest PLV
fprintf('The period with the highest PLV is from sample %d to %d with PLV = %.4f\n', ...
    (maxPLVIndex-1)*periodLength + 1, maxPLVIndex*periodLength, maxPLV);

figure;
subplot(5, 1, 1);
plot(signal);
title('Original Signal');
subplot(5,1,2)
plot(sm_analyticSignal)
title('Reconstructed inverse signl - Envelope for peak at 22 days');
subplot(5,1,3)
plot(instantaneousPhase)
title('Instantaneous Phase')
subplot(5,1,4)
plot(plvValues)
title('Phase Locking values')


%%
% Initialize arrays to store the classification
numPeriods = floor(length(instantaneousPhase) / periodLength);
phaseClassification = strings(1, length(instantaneousPhase));

% Classify each period based on phase
for k = 1:numPeriods
    startIdx = (k-1)*periodLength + 1;
    endIdx = k*periodLength;
    phaseValues = instantaneousPhase(startIdx:endIdx);
    
    % Classify phase region
    for j = 1:length(phaseValues)
        phase = phaseValues(j);
        if abs(phase) < pi/4
            phaseClassification(startIdx + j - 1) = "Peak";
        elseif abs(phase - pi) < pi/4
            phaseClassification(startIdx + j - 1) = "Trough";
        elseif phase > 0 && phase < pi
            phaseClassification(startIdx + j - 1) = "Rising";
        else
            phaseClassification(startIdx + j - 1) = "Falling";
        end
    end
end

subplot(5,1,5)
% plot(signal);
hold on;
for i = 1:length(instantaneousPhase)
    if phaseClassification(i) == "Peak"
        plot(i, signal(i), 'r.'); % Red dot for peak
    elseif phaseClassification(i) == "Trough"
        plot(i, signal(i), 'b.'); % Blue dot for trough
    elseif phaseClassification(i) == "Rising"
        plot(i, signal(i), 'g.'); % Green dot for rising
    elseif phaseClassification(i) == "Falling"
        plot(i, signal(i), 'm.'); % Magenta dot for falling
    end
end
title('Signal with Classified Phases');
% legend('Signal', 'Peak', '','','Trough', 'Rising', 'Falling');

% ylim([0 11])
