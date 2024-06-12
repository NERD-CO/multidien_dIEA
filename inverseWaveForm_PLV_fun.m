function [] = inverseWaveForm_PLV_fun(iValu , peakIndices)
signal = iValu;
fs = 250;

signal = signal(~isnan(signal));
[cfs, ~] = cwt(signal, 'amor', fs);

% peakIndices needs to be a row vector of Frequency Identities
% peakIndices = [20, 35]; 

% Find the corresponding scales/frequencies for the peak indices
% i'm assuming the peakIndices correspond directly to frequencies
numFrequencies = size(cfs, 1); % Number of frequency bins in the CWT
freqIndices = peakIndices;

% Reconstruct the signal from each isolated peak component using ICWT
reconstructedPeaks = zeros(length(freqIndices), length(signal));

for i = 1:length(freqIndices)
    freqIndex = freqIndices(i);
    if freqIndex > 0 && freqIndex <= numFrequencies
        % Create a CWT matrix with only the current peak component
        peakCWT = zeros(size(cfs));
        peakCWT(freqIndex, :) = cfs(freqIndex, :);
        
        reconstructedSignal = icwt(peakCWT, 'amor');
        
        % Store the reconstructed signal
        reconstructedPeaks(i, :) = reconstructedSignal;
    end
end

periodLength = 31; % might need to make an input argument - talk to David********

analyticSignals = zeros(size(reconstructedPeaks));
sm_analyticSignals = zeros(size(reconstructedPeaks)); 
hilbTrans = zeros(size(reconstructedPeaks));
instantaneousPhaseAll = zeros(size(reconstructedPeaks));
plvValuesAll = cell(1,height(reconstructedPeaks));
maxPLVall = zeros(1,height(reconstructedPeaks));
maxPLVIndexAll = zeros(1,height(reconstructedPeaks));
highestPLVval = zeros(1,height(reconstructedPeaks));
highestPLVloc  = zeros(1,height(reconstructedPeaks));
for i = 1:length(freqIndices)

    analyticSignals(i,:) = envelope(reconstructedPeaks(i,:));
    sm_analyticSignals(i,:) = smoothdata(analyticSignals(i,:),'sgolay',250);
    hilbTrans(i,:) = hilbert(sm_analyticSignals(i,:));
    instantaneousPhaseAll(i,:) = angle(hilbTrans(i,:));

    numPeriods = floor(length(instantaneousPhaseAll(i,:)) / periodLength);
    plvValues = zeros(1, numPeriods);

    % Compute the PLV for each period
    for k = 1:numPeriods
        phaseValues = instantaneousPhaseAll(i,(k-1)*periodLength + 1 : k*periodLength);

        plvValues(k) = abs(mean(exp(1i * phaseValues)));
    end

    plvValuesAll{i} = plvValues;

    [maxPLV, maxPLVIndex] = max(plvValues);

    maxPLVall(i) = maxPLV;
    maxPLVIndexAll(i) = maxPLVIndex;

    highestPLVval(i) = (maxPLVIndex-1)*periodLength + 1;
    highestPLVloc(i) =  maxPLVIndex*periodLength;

    figure;
    subplot(5, 1, 1);
    plot(signal);
    title('Original Signal');
    subplot(5,1,2)
    plot(sm_analyticSignals(1,:))
    title('Reconstructed inverse signl - Envelope for peak at 22 days');
    subplot(5,1,3)
    plot(instantaneousPhaseAll(1,:))
    title('Instantaneous Phase')
    subplot(5,1,4)
    plot(plvValues)
    title('Phase Locking values')

    phaseClass = strings(1, length(instantaneousPhaseAll(i,:)));
    for k = 1:numPeriods
        startIdx = (k-1)*periodLength + 1;
        endIdx = k*periodLength;
        phaseValues = instantaneousPhaseAll(i,startIdx:endIdx);

        % Classify phase region
        for j = 1:length(phaseValues)
            phase = phaseValues(j);
            if abs(phase) < pi/4
                phaseClass(startIdx + j - 1) = "Peak";
            elseif abs(phase - pi) < pi/4
                phaseClass(startIdx + j - 1) = "Trough";
            elseif phase > 0 && phase < pi
                phaseClass(startIdx + j - 1) = "Rising";
            else
                phaseClass(startIdx + j - 1) = "Falling";
            end
        end
    end

    subplot(5,1,5)
    % plot(signal);
    hold on;
    for ii = 1:length(instantaneousPhaseAll(i,:))
        if phaseClass(ii) == "Peak"
            plot(ii, signal(ii), 'r.'); % Red dot for peak
        elseif phaseClass(i) == "Trough"
            plot(ii, signal(ii), 'b.'); % Blue dot for trough
        elseif phaseClass(i) == "Rising"
            plot(ii, signal(ii), 'g.'); % Green dot for rising
        elseif phaseClass(i) == "Falling"
            plot(ii, signal(ii), 'm.'); % Magenta dot for falling
        end
    end
    title('Signal with Classified Phases');



end





end