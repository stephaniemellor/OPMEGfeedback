%% HouseKeeping
clear; close all
addpath('D:\Documents\Software\spm12\PublicVersion');
spm('defaults', 'eeg')
addpath('D:\Documents\Software\FIL-OPMEG\OPM');
addpath('D:\Documents\Software\linspecer');
addpath('D:\Documents\Software\export_fig');

%% Config

resultsSavePath = '.\EnvironmentalNoise\results';
rawDataPath = '.\EnvironmentalNoise\sub-001\ses-001\meg';

%% Read 

% ------------ Feeding back to selected channels

% No filter / technically 100 Hz moving average filter
cd(rawDataPath);

S = [];
S.data = 'sub-001_ses-001_task-noise_run-002_meg.bin';
S.positions = 'positions.tsv';
S.path = '.\EnvironmentalNoise\analysedData';
D = spm_opm_create(S);
cd(S.path)

% Set DQ to bad
D = badchannels(D, selectchannels(D, 'regexp_G2-DQ.*'), 1);
save(D);

% Downsample
S = [];
S.D = D;
S.fsample_new = 1000;
D = spm_eeg_downsample(S);

% Save
data_selected_channels.("noFilter").data = fullfile(D.path, D.fname);
data_selected_channels.("noFilter").feedback_channels = ...
    intersect(D.chanlabels(selectchannels(D, 'regexp_G2-[A|D].*')), ...
    D.sensors('MEG').label);

% 1 Hz filter
cd(rawDataPath);
S = [];
S.data = 'sub-001_ses-001_task-noise_run-001_meg.bin';
S.positions = 'positions.tsv';
S.path = '.\EnvironmentalNoise\analysedData';
D = spm_opm_create(S);
cd(S.path)

% Set DQ to bad
D = badchannels(D, selectchannels(D, 'regexp_G2-DQ.*'), 1);
save(D);

% Downsample
S = [];
S.D = D;
S.fsample_new = 1000;
D = spm_eeg_downsample(S);

% Save
data_selected_channels.("filter").data = fullfile(D.path, D.fname);
data_selected_channels.("filter").feedback_channels = ...
    data_selected_channels.("noFilter").feedback_channels;
data_selected_channels.("filter").filter_cutoff = 1;

%% Time series and power spectrum for feedback with filter

% Get data
D = spm_eeg_load(data_selected_channels.filter.data);

meginds = intersect(indchannel(D, D.sensors('MEG').label), ...
    indchantype(D, 'MEGMAG', 'GOOD'));
feedbackinds = intersect(meginds, ...
    indchannel(D, data_selected_channels.filter.feedback_channels));
meginds = setdiff(meginds, feedbackinds);

% Plot
figure;
C = linspecer(size(D,1)*3, 'sequential');

% Time series
subplot(1,3,1); 
grid on; box on; hold on;
for ind = 1:length(meginds)
    if ind == floor(length(meginds)/2)
        h = plot(D.time, 1e-6*D(meginds(ind),:,1), '-', 'LineWidth', 0.5, 'color', C(ind,:));
    else
        plot(D.time, 1e-6*D(meginds(ind),:,1), '-', 'LineWidth', 0.5, 'color', C(ind,:));
    end
end
for ind = 1:length(feedbackinds)
    if ind == floor(length(feedbackinds)/2)
        g = plot(D.time, 1e-6*D(feedbackinds(ind),:,1), '-', 'LineWidth', 0.5, 'color', C(end-ind+1,:));
    else
        plot(D.time, 1e-6*D(feedbackinds(ind),:,1), '-', 'LineWidth', 0.5, 'color', C(end-ind+1,:));
    end
end
yl = [-2.5 2.5];
plot([300 300], yl, 'k--');
plot([max(D.time)-300, max(D.time)-300], yl, 'k--')
ylim(yl);
xlim([min(D.time), max(D.time)]);
set(gca, 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 18);
ylabel('B (nT)', 'FontSize', 18);
legend([h, g], {'No feedback', 'With feedback'}, 'location', 'northeast', 'FontSize', 16);
title('Time Series');
set(gcf, 'Position', [680   585   808   411]); 

% Cut to just central 30 minutes
S = [];
S.D = D;
S.timewin = [5*60*D.fsample, 35*60*D.fsample];
D = spm_eeg_crop(S);

% PSD
S=[];
S.D=D;
S.triallength=20*1e3;
[p,f] = spm_opm_psd(S);

subplot(1,3,2); hold on; box on; grid on;
for ind = 1:length(meginds)
    if ind == floor(length(meginds)/2)
        h = semilogy(f,p(:,meginds(ind)),'-','LineWidth',0.5, 'color', C(ind,:));
    else
        semilogy(f,p(:,meginds(ind)),'-','LineWidth',0.5, 'color', C(ind,:));
    end
end
semilogy(f,median(p(:,meginds),2),'k-','LineWidth',1.5);

for ind = 1:length(feedbackinds)
    if ind == floor(length(feedbackinds)/2)
        g = semilogy(f,p(:,feedbackinds(ind)),'-','LineWidth',0.5, 'color', C(end-ind+1,:));
    else
        semilogy(f,p(:,feedbackinds(ind)),'-','LineWidth',0.5, 'color', C(end-ind+1,:));
    end
end
semilogy(f,median(p(:,feedbackinds),2),'k-','LineWidth',1.5);
plot(f, 15*ones(size(f)), 'k--')
xlabel('Frequency (Hz)')
labY = ['$$PSD (fT\sqrt[-1]{Hz}$$)'];
ylabel(labY,'interpreter','latex')
xlim([0,250])
ylim([1,1e6])
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend([h, g], {'No feedback', 'With feedback'}, 'location', 'southeast')
title('Power Spectral Density (PSD)');
xticks(gca, 10.^(-1:1:2))
xticklabels(string(10.^(-1:1:2)))

% Estimate rpsd (isn't true as don't have the same channels
mpf = median(p(:,feedbackinds),2);
sepf = 1.2533*std(p(:,feedbackinds),[],2)./sqrt(length(feedbackinds));
mpr = median(p(:,meginds),2);
sepr = 1.2533*std(p(:,meginds),[],2)./sqrt(length(meginds));

rpsd = 20*log10(mpr./mpf);
se_rpsd = 20*sqrt((sepf./mpf).^2 + (sepr./mpr).^2)./log(10);

subplot(1,3,3); hold on; box on; grid on;
fill([f(2:end)';flipud(f(2:end)')], [rpsd(2:end)-se_rpsd(2:end); flipud(rpsd(2:end)+se_rpsd(2:end))], [0, 0, 0],...
    'linestyle', 'none', 'FaceAlpha', 0.4)
plot(f, rpsd, 'color', [0, 0, 0], 'LineWidth', 1.5);
xlim([0 250]);
yl = [-40 40];
ylim([-max(abs(yl)), max(abs(yl))])
grid on;
xl = xlim;
xl(1) = f(2);
plot(xl, [0 0], 'k--')
xlim(xl);
xlabel('Frequency (Hz)');
ylabel('Relative median psd (dB)')
title('Relative median PSD')
text(f(2),rpsd(2),['\leftarrow ', sprintf('%.f',rpsd(2)), ' \pm ', sprintf('%.f',se_rpsd(2)), ' dB'], ...
    'FontSize', 20);
set(gca, 'XScale', 'log');
xticks(gca, 10.^(-1:1:2))
xticklabels(string(10.^(-1:1:2)))

for plt = 1:3
    axes(subplot(1,3,plt))
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [2 380 1785 720])
set(gcf, 'color', 'w');
export_fig(fullfile(resultsSavePath, 'filter_demonstration_log'), '-png', '-painters')

%% Time series and power spectrum for feedback without filter

% Get data
D = spm_eeg_load(data_selected_channels.noFilter.data);

meginds = intersect(indchannel(D, D.sensors('MEG').label), ...
    indchantype(D, 'MEGMAG', 'GOOD'));
feedbackinds = intersect(meginds, ...
    indchannel(D, data_selected_channels.noFilter.feedback_channels));
meginds = setdiff(meginds, feedbackinds);

% Plot
figure;
C = linspecer(size(D,1)*3, 'sequential');

% Time series
subplot(1,3,1); 
grid on; box on; hold on;
for ind = 1:length(meginds)
    if ind == floor(length(meginds)/2)
        h = plot(D.time, 1e-6*D(meginds(ind),:,1), '-', 'LineWidth', 0.5, 'color', C(ind,:));
    else
        plot(D.time, 1e-6*D(meginds(ind),:,1), '-', 'LineWidth', 0.5, 'color', C(ind,:));
    end
end
for ind = 1:length(feedbackinds)
    if ind == floor(length(feedbackinds)/2)
        g = plot(D.time, 1e-6*D(feedbackinds(ind),:,1), '-', 'LineWidth', 0.5, 'color', C(end-ind+1,:));
    else
        plot(D.time, 1e-6*D(feedbackinds(ind),:,1), '-', 'LineWidth', 0.5, 'color', C(end-ind+1,:));
    end
end
yl = [-2.5 2.5];
plot([300 300], yl, 'k--');
plot([max(D.time)-300, max(D.time)-300], yl, 'k--')
ylim(yl);
xlim([min(D.time), max(D.time)]);
set(gca, 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 18);
ylabel('B (nT)', 'FontSize', 18);
legend([h, g], {'No feedback', 'With feedback'}, 'location', 'northeast', 'FontSize', 16);
title('Time Series');
set(gcf, 'Position', [680   585   808   411]); 

% Cut to just central 30 minutes
S = [];
S.D = D;
S.timewin = [5*60*D.fsample, 35*60*D.fsample];
D = spm_eeg_crop(S);

% PSD
S=[];
S.D=D;
S.triallength=20*1e3;
[p,f] = spm_opm_psd(S);

subplot(1,3,2); hold on; box on; grid on;
for ind = 1:length(meginds)
    if ind == floor(length(meginds)/2)
        h = semilogy(f,p(:,meginds(ind)),'-','LineWidth',0.5, 'color', C(ind,:));
    else
        semilogy(f,p(:,meginds(ind)),'-','LineWidth',0.5, 'color', C(ind,:));
    end
end
semilogy(f,median(p(:,meginds),2),'k-','LineWidth',1.5);

for ind = 1:length(feedbackinds)
    if ind == floor(length(feedbackinds)/2)
        g = semilogy(f,p(:,feedbackinds(ind)),'-','LineWidth',0.5, 'color', C(end-ind+1,:));
    else
        semilogy(f,p(:,feedbackinds(ind)),'-','LineWidth',0.5, 'color', C(end-ind+1,:));
    end
end
semilogy(f,median(p(:,feedbackinds),2),'k-','LineWidth',1.5);
plot(f, 15*ones(size(f)), 'k--')
xlabel('Frequency (Hz)')
labY = ['$$PSD (fT\sqrt[-1]{Hz}$$)'];
ylabel(labY,'interpreter','latex')
xlim([0,250])
ylim([1,1e6])
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend([h, g], {'No feedback', 'With feedback'}, 'location', 'southeast')
title('Power Spectral Density (PSD)');
xticks(gca, 10.^(-1:1:2))
xticklabels(string(10.^(-1:1:2)))

% Estimate rpsd (isn't true as don't have the same channels
mpf = median(p(:,feedbackinds),2);
sepf = 1.2533*std(p(:,feedbackinds),[],2)./sqrt(length(feedbackinds));
mpr = median(p(:,meginds),2);
sepr = 1.2533*std(p(:,meginds),[],2)./sqrt(length(meginds));

rpsd = 20*log10(mpr./mpf);
se_rpsd = 20*sqrt((sepf./mpf).^2 + (sepr./mpr).^2)./log(10);

subplot(1,3,3); hold on; box on; grid on;
fill([f(2:end)';flipud(f(2:end)')], [rpsd(2:end)-se_rpsd(2:end); flipud(rpsd(2:end)+se_rpsd(2:end))], [0, 0, 0],...
    'linestyle', 'none', 'FaceAlpha', 0.4)
plot(f, rpsd, 'color', [0, 0, 0], 'LineWidth', 1.5);
xlim([0 250]);
yl = [-40 40];
ylim([-max(abs(yl)), max(abs(yl))])
grid on;
xl = xlim;
xl(1) = f(2);
plot(xl, [0 0], 'k--')
xlim(xl);
xlabel('Frequency (Hz)');
ylabel('Relative median psd (dB)')
title('Relative median PSD')
text(f(2),rpsd(2),['\leftarrow ', sprintf('%.f',rpsd(2)), ' \pm ', sprintf('%.f',se_rpsd(2)), ' dB'], ...
    'FontSize', 20);
set(gca, 'XScale', 'log');
xticks(gca, 10.^(-1:1:2))
xticklabels(string(10.^(-1:1:2)))

for plt = 1:3
    axes(subplot(1,3,plt))
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [2 380 1785 720])
set(gcf, 'color', 'w');
export_fig(fullfile(resultsSavePath, 'no_filter_demonstration_log'), '-png', '-painters')

%% rPSD with selected sensors feeding back

psd_trial_length = 10; % s

% Compare different filters
C = linspecer(4);
h = gobjects(1+length(data_selected_channels.filter),1);

% Just plot rpsd
figure; hold on; grid on; box on;

% Plot standard errors first
% No filter
% Load data
D = spm_eeg_load(data_selected_channels.("noFilter").data);

% Get meg and feedback channel inds
meginds = intersect(indchannel(D, D.sensors('MEG').label), ...
    indchantype(D, 'MEGMAG', 'GOOD'));
feedbackinds = intersect(meginds, ...
    indchannel(D, data_selected_channels.("noFilter").feedback_channels));
meginds = setdiff(meginds, feedbackinds);

% Cut to just central 30 minutes
S = [];
S.D = D;
S.timewin = [5*60*D.fsample, 35*60*D.fsample];
D = spm_eeg_crop(S);

% PSD
S=[];
S.D=D;
S.triallength=psd_trial_length*1e3;
[p,f] = spm_opm_psd(S);

% estimate rPSD
mpf = median(p(:,feedbackinds),2);
sepf = std(p(:,feedbackinds),[],2);
mpr = median(p(:,meginds),2);
sepr = std(p(:,meginds),[],2);

rpsd = 20*log10(mpr./mpf);
se_rpsd = 20*sqrt((sepf./mpf).^2 + (sepr./mpr).^2)./log(10);

% Plot
fill([f';flipud(f')], [rpsd-se_rpsd; flipud(rpsd+se_rpsd)], [0,0,0],...
    'linestyle', 'none', 'FaceAlpha', 0.2)
h(1) = plot(f, rpsd, 'color', [0, 0, 0], 'LineWidth', 1.5, 'DisplayName', 'No filter');
hCounter = 2;

% filters
for fc = 1:length(data_selected_channels.filter)

    % Load data
    D = spm_eeg_load(data_selected_channels.filter(fc).data);

    % Cut to just central 30 minutes
    S = [];
    S.D = D;
    S.timewin = [5*60*D.fsample, 35*60*D.fsample];
    D = spm_eeg_crop(S);

    % Get meg and feedback channel inds
    meginds = intersect(indchannel(D, D.sensors('MEG').label), ...
        indchantype(D, 'MEGMAG', 'GOOD'));
    feedbackinds = intersect(meginds, ...
        indchannel(D, data_selected_channels.filter(fc).feedback_channels));
    meginds = setdiff(meginds, feedbackinds);

    % PSD
    S=[];
    S.D=D;
    S.triallength=psd_trial_length*1e3;
    [p,f] = spm_opm_psd(S);

    % estimate rPSD
    mpf = median(p(:,feedbackinds),2);
    sepf = std(p(:,feedbackinds),[],2);
    mpr = median(p(:,meginds),2);
    sepr = std(p(:,meginds),[],2);
    
    rpsd = 20*log10(mpr./mpf);
    se_rpsd = 20*sqrt((sepf./mpf).^2 + (sepr./mpr).^2)./log(10);

    % Plot
    fill([f';flipud(f')], [rpsd-se_rpsd; flipud(rpsd+se_rpsd)], C(fc,:),...
        'linestyle', 'none', 'FaceAlpha', 0.2)
    h(hCounter) = plot(f, rpsd, 'color', C(fc,:), 'LineWidth', 1.5, 'DisplayName',...
        ['Point-by-point filter cut-off = ' num2str(data_selected_channels.filter(fc).filter_cutoff), ' Hz']);
    hCounter = hCounter + 1;
end

% Move standard errors to the back
chH = get(gca,'Children');
set(gca,'Children',[chH(1:2:end);chH(2:2:end)])

% Other appearance options
xlim([0 300])
yl = [-40 40];
ylim([-max(abs(yl)), max(abs(yl))])
grid on;
plot([0 300], [0 0], 'k--')
xlabel('Frequency');
ylabel('Relative median psd (dB)')
title('Relative median PSD')
legend(h);

% Save
set(gca, 'FontSize', 16)
set(gcf, 'color', 'w');
set(gcf, 'Position', [680   585   642   513]);
export_fig(fullfile(resultsSavePath, 'filter_comparison'), '-png', '-nocrop', '-painters')

xlim([0 15]); ylim([-20 30])
export_fig(fullfile(resultsSavePath, 'filter_comparison_zoomed_in'), '-png', '-nocrop', '-painters')
