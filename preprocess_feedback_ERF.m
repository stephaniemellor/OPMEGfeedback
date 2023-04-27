function [D, movement_data, feedback_applied, saturated_frames, beep_start_stop] = ...
    preprocess_feedback_ERF(raw_data_filepath_root, analysed_data_path, results_save_path, MRI_file_name, subject, session, task, run, varargin)
% Preprocess auditory ERF data for feedback experiment. 
% Steps: 
%   load data,
%   resample to 1 kHz, 
%   estimate feedback applied, 
%   load and sync optitrack,
%   regress position and rotation (off by default), 
%   HFC, 
%   filter, 
%   epoch, 
%   set trials containing saturated data to bad
%
% [D, movement_data, feedback_applied, saturated_frames] = 
%    preprocess_feedback_ERF(raw_data_filepath_root, subject, session, task, run, vargarin)
%
% Input:
%   - raw_data_filepath_root: file path to raw data folder (before
%   specifying subject) (string)
%   - MRI_file_name: file name and path of participant MRI (string)
%   - subject: subject name (string)
%   - session: session name (string)
%   - task: task name (string)
%   - run: run name (string)
% Input Options, string value pair:
%   - 'position_regress' (bool), default false
%   - 'HFC' (bool), default true
%   - 'filter' (bool), default true
%   - 'optitrack_R' (double matrix), rotation to apply to optitrack to put
%       it into the MSR coordinate frame. Default identity
%   - 'optitrack_T' (double vector), translation to apply to optitrack to put
%       it into the MSR coordinate frame. Default [0;0;0]
%   - 'optitrack_trigger' (string), default 'NI-TRIG-7'
%   - 'stim_trigger' (string), default 'NI-TRIG-1'
%   - 'calibration_scaling' (double), scaling on optitrack calibration,
%   default 1
%
% Function Dependencies:
%   - SPM12
%   - OPM toolbox       
%   - linspecer
%
% N.B. if an SPM file with the correct name already exists, it will be
% loaded rather than recreated. If you wish to recreate the file, delete it
% before running function. 

%% Parse inputs
defaults = struct('position_regress', false, ...
    'HFC', true, 'filter', true, 'optitrack_R', eye(3), 'optitrack_T', [0;0;0], ...
    'optitrack_trigger', 'NI-TRIG-7', 'stim_trigger', 'NI-TRIG-1', 'calibration_scaling', 1);  % define default values
params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end
clear defaults

%% Load OPM data

cd(analysed_data_path);

fpath = [raw_data_filepath_root, '\sub-', subject, '\ses-', session, '\meg\']; 

fnameRoot = ['sub-', subject, '_ses-', session, '_task-', task, '_run-', run];

if isfile([fnameRoot, '_meg.dat'])
    D = spm_eeg_load([fnameRoot, '_meg.dat']);
else
    S = [];
    S.data = [fpath, fnameRoot, '_meg.bin'];
    S.channels = [fpath, fnameRoot, '_channels.tsv'];
    S.meg = [fpath, fnameRoot, '_meg.json'];
    S.positions = [raw_data_filepath_root, '\sub-', subject, '\ses-', session, '\positions.tsv'];
    S.sMRI = MRI_file_name;
    S.path = pwd;
    D = spm_opm_create(S);
end

%% Resample to 1000 Hz

if isfile(['d', D.fname])
    D = spm_eeg_load(['d', D.fname]);
else
    S = [];
    S.D = D;
    S.fsample_new = 1000;
    D = spm_eeg_downsample(S);
end

%% Estimate feedback applied

fedback = estimateFeedback([fpath, fnameRoot, '_feedback.bin'],...
        [fpath, fnameRoot, '_digitalChannels.tsv'], D);

% Convert from pT to fT to match D
fedback = 1e3*fedback;

% Plot PSD - taken from spm_opm_psd
chans = selectchannels(D, D.sensors('MEG').label);
labs = chanlabels(D,chans);
fs = D.fsample();

N = round(3e3);
begs = 1:N:size(fedback,2);
ends= begs+N-1;

if(ends(end)>size(fedback,2))
    ends(end)=[];
    begs(end)=[];
end
nepochs=length(begs);

wind  = window(@hanning,N);
Nf= ceil((N+1)/2);
coFac= max(wind)/mean(wind);
wind = repmat(wind,1,length(chans));

%- create PSD
%--------------------------------------------------------------------------
freq = 0:fs/N:fs/2;
odd=mod(N,2)==1;
psdx= zeros(length(freq),length(chans));

for j = 1:nepochs
    inds = begs(j):ends(j);
    Btemp=fedback(chans,inds,1)';

    % window data and correct for amplitude loss;
    Btemp = Btemp.*wind*coFac;
    fzf=Btemp;

    % fourier transform data and get RMS
    xdft = fft(fzf);
    xdft = xdft(1:floor(N/2+1),:);
    tmppsd = abs(xdft)./sqrt(N*fs);

    if(odd)
        tmppsd(2:end) = tmppsd(2:end);
    else
        tmppsd(2:end-1) = tmppsd(2:end-1);
    end

    % accumulate avearge PSD to limit memory usage
    psdx=psdx + tmppsd/nepochs;
end

%- plot
po = psdx;
f= figure();
hold on
for i = 1:size(po,2)
    tag = [labs{i}, ', Index: ' num2str(indchannel(D,labs{i}))];
    plot(freq,po(:,i)','LineWidth',2,'tag',tag);
end
set(gca,'yscale','log')
set(gca,'xscale','log')

% xp2 =0:round(freq(end));
% yp2=ones(1,round(freq(end))+1)*15;
% p2 =plot(xp2,yp2,'--k');
% p2.LineWidth=2;
p3=semilogy(freq,median(po,2),'LineWidth',2,'tag','Median');
p3.Color='k';
xlabel('Frequency (Hz)')
labY = ['$$PSD (fT\sqrt[-1]{Hz}$$)'];
ylabel(labY,'interpreter','latex')
grid on
ax = gca; % current axes
ax.FontSize = 16;
ax.TickLength = [0.02 0.02];
fig= gcf;
fig.Color=[1,1,1];
xlim([0,100]);
datacursormode on
dcm = datacursormode(gcf);
set(dcm,'UpdateFcn',@getLabel)
xlim([0 150]);
xticks(gca, 10.^(0:1:2))
xticklabels(string(10.^(0:1:2)))

export_fig(fullfile(results_save_path, sprintf('fedback_psd_participant_%s_session_%s_run_%s', subject, session, run)), '-png', '-painters')

%% Load optitrack

movement_data_orig = csv2mat_sm(fullfile(fpath, [fnameRoot, '_optitrack.csv']));

% Synchronise with opm data
[movement_data, D, D_trim_samples] = syncOptitrackAndOPMdata(movement_data_orig, ...
    D, 'TriggerChannelName', params.optitrack_trigger);

% Trim feedback
feedback_applied = fedback(:,D_trim_samples(1):D_trim_samples(2));

% Apply any input rotations and translations to the optitrack
movement_data = applyTransformationToOptiData(movement_data, params.optitrack_R, params.optitrack_T);

% Apply correction for calibration scaling
S = params.calibration_scaling*eye(3);
Tt = [0;0;0];
movement_data = applyTransformationToOptiData(movement_data, S, Tt);

% Adjust for ground plane being on the floor, so centre of room is 1.1m
% above it
R = eye(3);
T = [0; -1.1; 0];
movement_data = applyTransformationToOptiData(movement_data, R, T);

%% Identify saturated periods of OPM data

saturated_frames = whichDataSaturated(1e-6*D(indchantype(D, 'MEGMAG', 'GOOD'),:,1), 'n_sat_bins', 5, 'n_good_bins', 5, 'edge_step', 0.001);
title(['Run ', run])

all_sat_frames = cell2mat(saturated_frames');
all_sat_frames = unique(all_sat_frames);

%% Regress position and rotation

if params.position_regress

    if isfile(['mr_', D.fname])
        D = spm_eeg_load(['mr_', D.fname]);
    else
        % Get unsaturated frames
        unsatFrames = setdiff(1:size(D,2), all_sat_frames);
    
        % Create regressors
        ref                 = cat(2, quat2eul(movement_data.rigidbodies.data(unsatFrames,[4,1:3])), ...
            movement_data.rigidbodies.data(unsatFrames,5:7));
        % LP-filter optitrack data
        [ref]          = ft_preproc_lowpassfilter(ref', D.fsample, 2, 5);
        ref            = ref';
    
        D = spm_regress_motive_OPMdata(D, ref, 10, unsatFrames);
    
        % Shielding factor
        figure;
        S = [];
        S.dB = 1;
        S.plot = 1;
        Dfname = D.fname;
        S.D1 = spm_eeg_load(fullfile(D.path, Dfname(4:end)));
        S.D2 = D;
        S.triallength = 50*1e3;
        S.channels = intersect(D.sensors('MEG').label, ...
            D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD')));
        spm_opm_rpsd(S);
        grid on;
        xlim([0 5]);
        ylim([-20 40]);
        legend('off')
        title(['Shielding Factor from movement regression, block ', run])
        export_fig(fullfile(results_save_path, sprintf('shielding_factor_movreg_%s', Dfname)), '-png', '-painters')
    end
end

%% HFC

if params.HFC
    if isfile(['h', D.fname])
        D = spm_eeg_load(['h', D.fname]);
    else
        S = [];
        S.D = D;
        D = spm_opm_hfc(S);
    
        % Shielding factor
        figure;
        S = [];
        S.dB = 1;
        S.plot = 1;
        Dfname = D.fname;
        S.D1 = spm_eeg_load(fullfile(D.path, Dfname(2:end)));
        S.D2 = D;
        S.triallength = 5*1e3;
        S.channels = intersect(D.sensors('MEG').label, ...
            D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD')));
        spm_opm_rpsd(S);
        grid on;
        xlim([0 150]);
        ylim([-20 40]);
        legend('off')
        title(['Shielding Factor from HFC, block ', run])
        export_fig(fullfile(results_save_path, sprintf('shielding_factor_HFC_%s', Dfname)), '-png', '-painters')
    end
end

%% Filter

if params.filter
    if isfile(['fffff', D.fname])
        D = spm_eeg_load(['fffff', D.fname]);
    else
        S = [];
        S.band = 'stop';
        S.freq = [48 52];
        S.D = D;
        D = spm_eeg_filter(S);
    
        S = [];
        S.band = 'stop';
        S.freq = [118 122];
        S.D = D;
        D = spm_eeg_filter(S);
    
        S = [];
        S.band = 'stop';
        S.freq = [81 85];
        S.D = D;
        D = spm_eeg_filter(S);
    
        S = [];
        S.band = 'high';
        S.freq = 2;
        S.D = D;
        D = spm_eeg_filter(S);
    
        S = [];
        S.band = 'low';
        S.freq = 40;
        S.D = D;
        S.order = 6;
        D = spm_eeg_filter(S);
    
    end
end

%% Epoch

S = [];
S.D = D;
S.timewin = [-200 500];
S.condlabels = {'tone'};
S.triggerChannels = {params.stim_trigger};
[D, trl] = spm_opm_epoch_trigger(S);

beep_start_stop = [trl(1,1), trl(end,2)];

%% Plot time series and optitrack

figure; 
t = tiledlayout(7,1,'TileSpacing','Compact');

% OPM
C = linspecer(length(indchantype(D, 'MEGMAG', 'GOOD')));

Dfname = D.fname;
start_of_base_name = strfind(Dfname, 't_dsub');

% Load data
DD = spm_eeg_load(Dfname(start_of_base_name:end));
    
% Plot
ax1 = nexttile(t, [3,1]); grid on; hold on; box on;

% Plot saturated datapoints
if ~isempty(all_sat_frames)
    times = DD.time(all_sat_frames)-DD.time(beep_start_stop(1));
    p = plot(times, 1.8*ones(size(times)), 'k.');
end

% OPM data
tt = linspace(0, DD.time(beep_start_stop(2)) - DD.time(beep_start_stop(1)), ...
    length(beep_start_stop(1):beep_start_stop(2)));
plot(tt, 1e-6*DD(indchantype(DD, 'MEGMAG', 'GOOD'), ...
    beep_start_stop(1):beep_start_stop(2), 1), 'LineWidth', 1.5);
ylim([-2 2])
set(ax1, 'ColorOrder', C)
set(ax1, 'FontSize', 13)
set(ax1, 'xticklabel', []);
ylabel(ax1, 'B (nT)')
if ~isempty(all_sat_frames)
    legend(p, 'Saturated Datapoints', 'location', 'southeast');
end

% Optitrack
good_opt_data = movement_data.rigidbodies.data(:,8) > 0;

good_opt_data(1:beep_start_stop(1)) = 0;
good_opt_data(beep_start_stop(2)+1:end) = 0;
first_good_ind = find(good_opt_data, 1);

% Find gaps
[~, step_out] = findpeaks(diff(~good_opt_data));
[~, step_in] = findpeaks(-diff(~good_opt_data));
step_in = step_in + 1;

% Take off 8 datapoints for interpolation
step_in = step_in + 8;

tt = movement_data.time - movement_data.time(beep_start_stop(1));

% Translation
ax2 = nexttile(t, [2, 1]); grid on; hold on; box on;
c = get(ax2,'ColorOrder');
for gap = 1:length(step_in)
    for ll = 1:3
        plot(tt(step_in(gap):step_out(gap)), movement_data.rigidbodies.data(step_in(gap):step_out(gap),ll+4) - ...
            movement_data.rigidbodies.data(first_good_ind,ll+4), 'color', c(ll,:), 'LineWidth', 2);
    end
end
ylim([-1.6 1.6])
set(ax2, 'FontSize', 13)
set(ax2, 'xticklabel', []);
ylabel('Transalation (m)')
% legend({'X, left-right', 'Y, up-down', 'Z, forward-back'}, 'location', 'best');

% Rotation
ax3 = nexttile(t, [2, 1]); grid on; hold on; box on;
eulZYX = quat2eul(movement_data.rigidbodies.data(:,[4, 1:3]));
eulZYX = eulZYX - eulZYX(first_good_ind, :);
eulZYX = eulZYX*180/pi;
for gap = 1:length(step_in)
    counter = 3;
    for ll = [1,3,2]
        counter = counter + 1;
        plot(tt(step_in(gap):step_out(gap)), eulZYX(step_in(gap):step_out(gap),ll), 'color', c(counter,:), 'LineWidth', 2);
    end
end
ylim([-75 75])
linkaxes([ax1, ax2, ax3], 'x');
xlim([0 movement_data.time(beep_start_stop(2)) - movement_data.time(beep_start_stop(1))]);
set(gca, 'FontSize', 13)
xlabel('Time (s)')
label = ylabel(['Rotation (',char(176),')']);
label.Position(1) = label.Position(1)-6.0;
set(gcf, 'color', 'w');
set(gcf, 'Position', [680   673   711   425]);
% legend({'Yaw', 'Pitch', 'Roll'}, 'location', 'bestoutside');
export_fig(fullfile(results_save_path, sprintf('raw_time_series_participant_%s_session_%s_run_%s', subject, session, run)), '-png', '-painters');


%% Set Trials containing saturated data to bad

badtriallist = [];
for tt = 1:size(trl,1)
    if any((all_sat_frames < trl(tt,2)).*...
            (all_sat_frames > trl(tt,1)))
        badtriallist = cat(1, badtriallist, tt);
    end
end
D = badtrials(D, badtriallist, 1);
save(D);

%% Get number of badtrials

n_sat_trials = length(badtriallist);
n_total_trials = size(trl,1);

fprintf('\n Run %s, no trials = %.f, no sat trial = %.f\n\n', run, n_total_trials, n_sat_trials);

% Estimate time saturated
saturatedFramesClipped = all_sat_frames(find(all_sat_frames > trl(1,1), 1 ):...
    find(all_sat_frames < trl(end,2), 1, 'last' ));

time_sat = length(saturatedFramesClipped);
time_sat = time_sat/D.fsample;
time_total = diff(beep_start_stop);
time_total = time_total/D.fsample;

fprintf('\n Run %s, total time = %.f, saturated time = %.f\n\n ', run, time_total, time_sat);


end