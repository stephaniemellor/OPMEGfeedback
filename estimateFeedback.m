function [fedback] = estimateFeedback(feedback_bin,digital_chans_tsv, D, interpMethod)
% Read field that was fedback to the OPMs, reorder the channels and 
% interpolate to match the analogue output 

% If no method specified, use linear interpolation
if nargin < 4
    interpMethod = 'linear';
end

% Read digital channels
digital_channels = spm_load(digital_chans_tsv);
digital_channels = fieldnames(digital_channels);
digital_channels = cell2mat(digital_channels);
digital_channels = string(digital_channels(:, end-1:end));

% Read feedback binary file
dat = fopen(feedback_bin);
data = fread(dat,Inf,'double',0,'b');
fclose(dat);
data = reshape(data, size(digital_channels, 1), 3, []);

% Reorder and interpolate to match analogue channels
analogue_chans = D.chanlabels;
fedback = zeros(size(D,1), size(D,2));
for chan = 1:length(analogue_chans)
    if endsWith(analogue_chans(chan), 'RAD') || endsWith(analogue_chans(chan), '-Z')
        fedback(chan,:) = interp1(linspace(1, size(D,2), size(data,3)),...
            squeeze(data(contains(digital_channels, analogue_chans{chan}(4:5)), 3, :)), 1:size(D,2), interpMethod);
    elseif endsWith(analogue_chans(chan), 'TAN') || endsWith(analogue_chans(chan), '-Y')
        fedback(chan,:) = interp1(linspace(1, size(D,2), size(data,3)),...
            squeeze(data(contains(digital_channels, analogue_chans{chan}(4:5)), 2, :)), 1:size(D,2), interpMethod);
    end
end

end

