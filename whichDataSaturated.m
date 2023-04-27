function [saturatedFrames] = whichDataSaturated(data, varargin)
% function [saturatedFrames] = whichDataSaturated(data)
% Find where the OPMs have railed/saturated during a recording
% Will return all the frames where any OPM has saturated
%
% It's quite fiddly so you may have to change the options to appropriately
% threshold the data. 
% 
% The function creates a histogram of the data for each
% channel, with a step of "edge_step" between the edges of the histogram
% bins. If there are more than double the number of datapoints in the last
% "n_sat_bins" bins than the previous "n_good_bins", any datapoints in the
% last "n_sat_bins" are marked as bad. 
% 
% You can also introduce a threshold so that if value above which data is 
% rejected is less than "min_sat_thresh", the data will not be marked as 
% saturated. This is less sophisticated than changing the histogram 
% parameters but a straightforward fix to avoid removing good data and 
% useful for recordings with very few saturated channels.
%
% INPUT:
%   - data (nchan x nt double array). In nT!! The data that may or may not be
%       saturated. Dimension order should be channels then time.
%   - options:
%       - plotbool (boolean)        Default: true
%       - edge_step (double)        Default: 0.001 (nT)
%       - n_sat_bins (double)       Default: 3
%       - n_good_bins (double)      Default: 3
%       - min_sat_thresh (double)   Default: 1 (nT)
%
% OUTPUT:
%   - saturatedFrames (cell array). Every datapoint at which each OPM
%       saturated.
%
% FUNCTION DEPENDENCIES:
%   - None
%
% TODO:
%   - add option to output saturated frames for each channel, rather than
%   all frames where any OPM has railed
%
% AUTHOR INFO
% Stephanie Mellor, stephanie.mellor.17@ucl.ac.uk
% 21/06/21

% Parse inputs
defaults = struct('plotbool', true, 'edge_step', 0.001, 'n_sat_bins', 3, ...
    'n_good_bins', 3, 'min_sat_thresh', 1);
params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end
clear defaults

n_sat_bins = params.n_sat_bins;
n_good_bins = params.n_good_bins;

% Choose edges of histogram bins based on max of data
datmax = max(abs(data(:)));

% Form intensity histogram of measured fields
edges = -datmax-params.edge_step:params.edge_step:datmax+params.edge_step;
hist = discretize(data, edges);
hist = hist';
    
% Find saturated data
saturatedFrames = cell(1,size(data,1));
for j = 1:size(data,1)
    saturatedFrames{j} = [];
    % Find first and last non-empty bins
    lb = max(hist(:,j));
    fb = min(hist(:,j));
    % Check if fb < -min_sat_thresh or lb > min_sat_thresh
    if edges(fb) < -params.min_sat_thresh
        % Check whether the number of values in the end n_sat_bins bins is more than
        % double the number in the next n_good_bins
        if sum(hist(:,j) <= fb+n_sat_bins-1) > 2*sum((hist(:,j) <= fb+n_good_bins+n_sat_bins-1).*(hist(:,j) > fb+n_sat_bins-1))
            % Find saturated points
            saturatedFrames{j} = cat(1, saturatedFrames{j}, find(hist(:,j) <= fb+n_sat_bins-1));
        end
    end
    if edges(lb) > params.min_sat_thresh
        if sum(hist(:,j) >= lb-n_sat_bins+1) > 2*sum((hist(:,j) >= lb-n_sat_bins-n_good_bins+1).*(hist(:,j) < lb-n_sat_bins+1))
            % Find saturated points
            saturatedFrames{j} = cat(1, saturatedFrames{j}, find(hist(:,j) >= lb-n_sat_bins+1));
        end
    end
end

% Plot
if params.plotbool
    all_sat_frames = cell2mat(saturatedFrames');
    all_sat_frames = unique(all_sat_frames);
    figure; hold on;
    plot(data');
    plot(all_sat_frames, datmax*ones(1,length(all_sat_frames)), 'rx');
    plot(all_sat_frames, -datmax*ones(1,length(all_sat_frames)), 'rx');
end

end