% External coil frequency sweep

clear; close all
addpath('D:\Documents\Software\spm12\PublicVersion');
spm('defaults', 'eeg')
addpath(genpath('D:\Documents\Software\FIL-OPMEG\OPM'));
addpath('D:\Documents\Software\linspecer');
addpath('D:\Documents\Software\export_fig');


%% Config

resultsSavePath = '.\ExternalCoils\results';
analysedDataPath = '.\ExternalCoils\analysedData';

%% Load data

frequency = {'(1)', '(5)', '1', '2', '4', '6', '8', '10'};
run = {'1', '2', '3'};

cd  '.\ExternalCoils\sub-001\ses-001\meg'

fedback = cell(length(frequency), length(run));

for f = 1:length(frequency)
    for r = 1:length(run)
        S = [];
        S.data= sprintf('sub-001_ses-001_task-%sHz_run-00%s_meg.bin', frequency{f}, run{r});
        S.positions='..\positions.tsv';
        S.path = analysedDataPath;
        D = spm_opm_create(S);

        fedback{f,r} = estimateFeedback(sprintf('sub-001_ses-001_task-%sHz_run-00%s_feedback.bin', frequency{f}, run{r}),...
            sprintf('sub-001_ses-001_task-%sHz_run-00%s_digitalChannels.tsv', frequency{f}, run{r}), D);
    end
end

cd(analysedDataPath);

%% Time series without feedback, with feedback, with feedback with filter

chaninds = intersect(indchantype(D, 'MEGMAG', 'GOOD'), indchannel(D, D.sensors('MEG').label));
chanind_subtracted = randi(length(chaninds)); % 48 used for thesis
chanind_D = chaninds(chanind_subtracted);
chanind = chanind_D;

figure;
frequencies = frequency;%{'(5)', '4', '10'};
run_plot = run(1:2);

h = gobjects(length(frequencies), length(run_plot));
t = tiledlayout(length(frequencies),length(run_plot));

for f = 1:length(frequencies)
    for r = 1:length(run_plot)
        D = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-00%s_meg.mat', frequencies{f}, run_plot{r}));

        h(f,r) = nexttile; 
        grid(h(f,r), 'on'); hold(h(f,r), 'on'); box(h(f,r), 'on');
        plot(h(f,r), D.time, 1e-6*(D(chaninds,:,1) - mean(D(chaninds,1:2*D.fsample,1),2)));
        
        % Get frequency of signal
        if contains(frequencies{f}, '(')
            freq = str2double(['0.', frequencies{f}(2)]);
        else
            freq = str2double(frequencies{f});
        end
        xlim([0 10/freq])
        set(h(f,r), 'FontSize', 14);
    end
    ylabel(h(f,1), sprintf('%.1f Hz', freq), 'Rotation', 0, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
end
linkaxes(h,'y');
yl = ylim(h(1));
ylim(h(1), [-max(abs(yl)), max(abs(yl))])

% Add shared title and axis labels
xlabel(t,'Time (s)', 'FontSize', 18)
ylabel(t,'B (nT)', 'FontSize', 18)

% ylabel(h(1,1), '0.5 Hz')
% ylabel(h(2,1), '4 Hz')
% ylabel(h(3,1), '10 Hz')

title(h(1,1), 'Feedback Off');
title(h(1,2), 'Feedback On')
% title(h(1,3), 'Feedback, 5 Hz filter')

% Move plots closer together
% xticklabels(h(1:2,:),{})
yticklabels(h(:,2),{})
t.TileSpacing = 'compact';

% Figure size
set(gcf, 'color', 'w');
set(gcf, 'Position', [680, 76, 840, 1021]);

export_fig(fullfile(resultsSavePath, 'coil_frequencies_time_series'), '-png', '-nocrop', '-painters')
%print(fullfile(resultsSavePath, 'coil_frequencies_time_series'), '-dpng', '-r300')

%% rPSD - feedback no filter vs no feedback

figure; hold on;
frequencies = frequency;%{'(5)', '4', '10'};
run_plot = run(1:2);
C = linspecer(length(frequencies));
h = gobjects(length(frequencies), 1);

for ff = 1:length(frequencies)
    D1 = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-001_meg.mat', frequencies{ff}));
    D2 = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-002_meg.mat', frequencies{ff}));

    S = [];
    S.D1 = D1;
    S.D2 = D2;
    S.triallength = 4*1e3;
    [shield, f] = spm_opm_rpsd(S);

    rpsd = median(shield,2);
    se_rpsd = std(shield,[],2);
    fill([f';flipud(f')], [rpsd-se_rpsd; flipud(rpsd+se_rpsd)], [0,0,0],...
        'linestyle', 'none', 'FaceAlpha', 0.1)

    if contains(frequencies{ff}, '(')
        freq = str2double(['0.', frequencies{ff}(2)]);
    else
        freq = str2double(frequencies{ff});
    end
    h(ff) = plot(f, rpsd, 'color', C(ff,:), 'DisplayName', sprintf('%.2f Hz', freq), 'LineWidth', 2);
end
xl = xlim;
plot(xl, [0, 0], 'k--')
set(gca, 'FontSize', 16)
lgd = legend(h);
lgd.Position = [0.7060, 0.5403, 0.1994, 0.3835];
xlim([0 12]);
ylim([-70 70]);
grid on; box on;

xlabel('Time (s)', 'FontSize', 18)
ylabel('Shielding Factor (dB)', 'FontSize', 18)

% Figure size
set(gcf, 'color', 'w');
set(gcf, 'Position', [598   437   642   541]);

export_fig(fullfile(resultsSavePath, 'coil_frequencies_rpsd'), '-png', '-nocrop', '-painters')

%% PSD without feedback, with feedback, with feedback with filter

C = linspecer(length(chaninds));

figure;
frequencies = {'(5)', '4', '10'};

h = gobjects(length(frequencies), length(run));
t = tiledlayout(length(frequencies),length(run));

for ff = 1:length(frequencies)
    for r = 1:length(run)
        D = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-00%s_meg.mat', frequencies{ff}, run{r}));

        S = [];
        S.D = D;
        S.triallength = 3000;
        [p,f] = spm_opm_psd(S);

        h(ff,r) = nexttile; 
        loglog(f,p(:,chaninds),'-','LineWidth',0.5);
        grid(h(ff,r), 'on'); hold(h(ff,r), 'on'); box(h(ff,r), 'on');
        plot(h(ff,r), f, 15*ones(size(f)), 'k--');
        if ff == 1
            plot(h(ff,r), [0.5, 0.5], [1 5e5], 'k--');
        elseif ff == 2
            plot(h(ff,r), [4, 4], [1 5e5], 'k--');
        elseif ff == 3
            plot(h(ff,r), [10, 10], [1 5e5], 'k--');
        end

        colororder(h(ff,r), C)
        
        xlim([0 20])
        set(h(ff,r), 'FontSize', 14);
    end
end
linkaxes(h,'xy');

% Add shared title and axis labels
xlabel(t,'Time (s)', 'FontSize', 18)
labY = ['\textsf{PSD }$\mathsf{(fT\sqrt[-1]{Hz})}$'];
ylabel(t, labY, 'interpreter', 'latex', 'FontSize', 18)

ylabel(h(1,1), '0.5 Hz')
ylabel(h(2,1), '4 Hz')
ylabel(h(3,1), '10 Hz')

title(h(1,1), 'No feedback');
title(h(1,2), 'Feedback, no filter')
title(h(1,3), 'Feedback, 5 Hz filter')

% Move plots closer together
xticklabels(h(1:2,:),{})
yticklabels(h(:,2:3),{})
t.TileSpacing = 'compact';

% Figure size
set(gcf, 'color', 'w');
set(gcf, 'Position', [680 430 840 630]);

%export_fig(fullfile(resultsSavePath, 'coil_frequencies_psd'), '-png', '-nocrop', '-painters')
print(fullfile(resultsSavePath, 'coil_frequencies_psd'), '-dpng', '-r300')


%% Subtract what was intended to be removed

figure;
frequencies = {'(5)', '4', '10'};

h = gobjects(length(frequencies), 3);
t = tiledlayout(length(frequencies), 3);

for f = 1:length(frequencies)
    D = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-001_meg.mat', frequencies{f}));

    f_ind = find(contains(frequency, frequencies{f}));
    subtracted = D(:,:,1)-mean(D(:,1:2*D.fsample,1),2) - 1e3*(fedback{f_ind,1}-mean(fedback{f_ind,1}(:,1:2*D.fsample),2));

    % Time Series
    % Plot original
    h(f,1) = nexttile((f-1)*3+1);
    hold(h(f,1), 'on'); 
    plot(h(f,1), D.time, 1e-6*(D(chanind,:,1)-mean(D(chanind,1:2*D.fsample,1),2)));
    plot(h(f,1), D.time, 1e-6*(1e3*fedback{f_ind,1}(chanind,:) - mean(1e3*fedback{f_ind,1}(chanind,1:2*D.fsample),2)));
    
%     yl = ylim;
%     k = kurtosis(1e-6*(D(chanind,1:2*D.fsample,1)-mean(D(chanind,1:2*D.fsample,1),2)));
%     text(h(f,1), min(D.time)+0.1, 0.22, sprintf('Kurtosis = %.2f', k), 'FontSize', 16);

    % Plot subtracted
    h(f,2) = nexttile((f-1)*3+2);
    plot(h(f,2), D.time, 1e-6*subtracted(chanind,:));

%     k = kurtosis(1e-6*subtracted(chanind,1:2*D.fsample));
%     text(h(f,2), min(D.time)+0.1, 0.22, sprintf('Kurtosis = %.2f', k), 'FontSize', 16);

    % Plot fedback
    D = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-002_meg.mat', frequencies{f}));
    h(f,3) = nexttile((f-1)*3+3);
    plot(h(f,3), D.time, 1e-6*(D(chanind,:)-mean(D(chanind,1:2*D.fsample),2)));

%     k = kurtosis(1e-6*(D(chanind,1:2*D.fsample)-mean(D(chanind,1:2*D.fsample),2)));
%     text(h(f,3), min(D.time)+0.1, 0.22, sprintf('Kurtosis = %.2f', k), 'FontSize', 16);

    for r = 1:3
        grid(h(f,r), 'on'); hold(h(f,r), 'on'); box(h(f,r), 'on');
        xlim(h(f,r), [0 2])
        set(h(f,r), 'FontSize', 14);
    end
end
linkaxes(h,'xy');

xlabel(t,'Time (s)', 'FontSize', 18)
ylabel(t, 'B (nT)', 'FontSize', 18)

ylabel(h(1,1), '0.5 Hz')
ylabel(h(2,1), '4 Hz')
ylabel(h(3,1), '10 Hz')

title(h(1,1), 'No Feedback');
title(h(1,2), 'Subtract Model')
title(h(1,3), 'Rec. With Feedback')

% Legend
leg = legend(h(3,1), {'Recorded', 'Model'});
leg_pos = leg.Position;
set(leg, 'Position', [0.175,0,leg_pos(3:4)]);

% Move plots closer together
xticklabels(h(1:2,:),{})
yticklabels(h(:,2:3),{})
t.TileSpacing = 'compact';

% Figure size
set(gcf, 'color', 'w');
set(gcf, 'Position', [680 430 840 630]);

export_fig(fullfile(resultsSavePath, 'subtracting_model_time_series'), '-png', '-nocrop', '-painters')
%print(fullfile(resultsSavePath, 'subtracting_model_time_series'), '-dpng', '-r300')

%% Try to add in time delay until recreate what was recorded (as well as possible)

% Find phase shift
D1 = spm_eeg_load('sub-001_ses-001_task-1Hz_run-001_meg.mat');
phase_shift = zeros(length(frequency), length(indchantype(D1, 'MEGMAG', 'GOOD')));
phi2 = zeros(size(phase_shift));
time_to_right = zeros(size(phase_shift));
frequencies = frequency;
data1_A = zeros(size(phase_shift));
data2_A = zeros(size(phase_shift));

for f = 3:length(frequencies)

    %f_ind = find(contains(frequency, frequencies{f}));
    f_ind = f;

    D1 = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-001_meg.mat', frequencies{f}));
    D2 = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-002_meg.mat', frequencies{f}));

    data1 = detrend(D1(indchantype(D1, 'MEGMAG', 'GOOD'),:,1)',1)';
    data2 = detrend(D2(indchantype(D2, 'MEGMAG', 'GOOD'),:,1)',1)';

    % Select model
    model = detrend(1e3*fedback{f_ind,1}(indchantype(D1, 'MEGMAG','GOOD'),:)', 1)';

    % Get frequency
    if contains(frequencies{f}, '(')
        freq = str2double(['0.', frequencies{f}(2)]);
    else
        freq = str2double(frequencies{f});
    end

    % Fit sine waves to data1, data2 and model
    data1_phi = zeros(1,size(data1,1));
    model_A = zeros(1, size(model, 1));
    for chan = 1:size(data1,1)
        data1_fit = fit(transpose(D1.time), data1(chan,:)', 'sin1', 'Lower', [0,floor(2*pi*freq),-pi], 'Upper', [Inf,ceil(2*pi*freq),pi]);
        data1_A(f,chan) = data1_fit.a1;
        data1_phi(chan) = data1_fit.c1;
    
        model_fit = fit(transpose(D1.time), model(chan,:)', 'sin1', 'Lower', [0,floor(2*pi*freq),-pi], 'Upper', [Inf,ceil(2*pi*freq),pi]);
        model_A(chan) = model_fit.a1;

        data2_fit = fit(transpose(D2.time), data2(chan,:)', 'sin1', 'Lower', [0,floor(2*pi*freq),-pi], 'Upper', [Inf,ceil(2*pi*freq),pi]);
        data2_A(f,chan) = data2_fit.a1;
    end

    % Find phase needed on model to get amplitude of data2
    for chan = 1:size(data1,1)
        phase_shift(f, chan) = -acos((data2_A(f,chan)^2 - model_A(chan)^2 - data1_A(f,chan)^2)/(2*data1_A(f,chan)*model_A(chan)));
    end
    phi2(f,:) = data1_phi - phase_shift(f,:);
    time_to_right(f,:) = phi2(f,:)./(2*pi*freq);
end

f = fit(cellfun(@str2double, frequencies(3:end))', median(data2_A(3:end,:)./data1_A(3:end,:), 2), 'poly1');
figure; hold on; box on;
plot(f, cellfun(@str2double, frequencies(3:end)), median(data2_A(3:end,:)./data1_A(3:end,:), 2))
errorbar(cellfun(@str2double, frequencies(3:end)), median(data2_A(3:end,:)./data1_A(3:end,:), 2), ...
    std(data2_A(3:end,:)./data1_A(3:end,:),[],2)/sqrt(size(data1_A,2)), 'b.');
legend(gca, 'off');
xlabel('External Coil Frequency (Hz)');
ylabel('A_3/A_1');
grid on; set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');
set(gcf, 'Position', [537 357 489 420]);


frequency_at_which_double = (2-f.p2)/f.p1;
time_to_right = 1/(2*frequency_at_which_double);
xlim([0 ceil(frequency_at_which_double)+0.6]);
ylim([0 2.1]);
plot([str2double(frequencies{end}), frequency_at_which_double], f([str2double(frequencies{end}), frequency_at_which_double]), 'r--')
plot(frequency_at_which_double, f(frequency_at_which_double), 'rx')

export_fig(fullfile(resultsSavePath, 'amplitude_vs_freq_best_fit'), '-png', '-nocrop', '-painters')


% % Insist that model lags data, not other way around
% if phase_shift > 0
%     phase_shift = -2*pi+phase_shift;
% end
% 
% model_phi = data1_phi - phase_shift;

% time_to_right = phase_shift./(2*pi*10);

samples_to_right = -round(time_to_right*D1.fsample);


%% Plot
figure;
frequencies = {'(5)', '4', '10'};

h = gobjects(length(frequencies), 4);
t = tiledlayout(length(frequencies), 4);

for f = 1:length(frequencies)

    D1 = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-001_meg.mat', frequencies{f}));
    data1 = detrend(D1(indchantype(D1, 'MEGMAG', 'GOOD'),:,1)',1)';
    f_ind = find(contains(frequency, frequencies{f}));
    model = detrend(1e3*fedback{f_ind,1}(indchantype(D1, 'MEGMAG','GOOD'),:)', 1)';

    % Test by shifting and subtracting
    data1_shifted = data1(:,-samples_to_right+1:end);
    model_shifted = model(:,1:end+samples_to_right);
    subtracted = data1_shifted-model_shifted;

    % Plot original
    h(f, 1) = nexttile((f-1)*4+1);
    hold(h(f,1), 'on');
    plot(h(f,1), D1.time(1:end-abs(samples_to_right)), 1e-6*(data1_shifted(chanind_subtracted,:)-...
        mean(data1_shifted(chanind_subtracted,1:2*D.fsample),2)));
    plot(h(f,1), D1.time(1:end-abs(samples_to_right)), 1e-6*(model_shifted(chanind_subtracted,:)-...
        mean(model_shifted(chanind_subtracted,1:2*D.fsample),2)));

    % Plot subtracted - no delay
    h(f,2) = nexttile((f-1)*4+2);
    plot(h(f,2), D1.time, 1e-6*(data1(chanind_subtracted,:)-model(chanind_subtracted,:)));

    % Plot subtracted - delay
    h(f,3) = nexttile((f-1)*4+3);
    plot(h(f,3), D1.time(1:end-abs(samples_to_right)), 1e-6*subtracted(chanind_subtracted,:));

    % Plot recorded with feedback
    D = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-002_meg.mat', frequencies{f}));
    h(f,4) = nexttile((f-1)*4+4);
    plot(h(f,4), D.time, 1e-6*(D(chanind,:)-mean(D(chanind,1:2*D.fsample),2)));

    for r = 1:4
        grid(h(f,r), 'on'); hold(h(f,r), 'on'); box(h(f,r), 'on');
        xlim(h(f,r), [0 2])
        set(h(f,r), 'FontSize', 14);
        ylim(h(f,r), [-0.2 0.2])
    end
end
linkaxes(h,'xy');

xlabel(t,'Time (s)', 'FontSize', 18)
ylabel(t, 'B (nT)', 'FontSize', 18)

ylabel(h(1,1), '0.5 Hz')
ylabel(h(2,1), '4 Hz')
ylabel(h(3,1), '10 Hz')

title(h(1,1), 'No Feedback');
title(h(1,2), 'Sub. Model - No delay')
title(h(1,3), 'Sub. Model - 42 ms delay')
title(h(1,4), 'Rec. With Feedback')

% Move plots closer together
yticklabels(h(:,2:end),{})
xticklabels(h(1:end-1,:),{})
t.TileSpacing = 'compact';

% Legend
leg = legend(h(3,1), {'Recorded', 'Model'});
leg_pos = leg.Position;
set(leg, 'Position', [0.1347,0.0016,0.1115,0.0737]);

% Figure size
set(gcf, 'color', 'w');
set(gcf, 'Position', [522 366 1275 630]);

export_fig(fullfile(resultsSavePath, 'lagged_model_time_series'), '-png', '-nocrop', '-painters')
%print(fullfile(resultsSavePath, 'lagged_delayed_model_time_series'), '-dpng', '-r300')

%% Looking just at 10 Hz case, repeat without shifting

D1 = spm_eeg_load('sub-001_ses-001_task-10Hz_run-001_meg.mat');
D2 = spm_eeg_load('sub-001_ses-001_task-10Hz_run-002_meg.mat');

data1 = detrend(D1(indchantype(D1, 'MEGMAG', 'GOOD'),:,1)',0)';
data2 = detrend(D2(indchantype(D2, 'MEGMAG', 'GOOD'),:,1)',0)';

% Select model
f_ind = find(contains(frequency, '10'));
model = detrend(1e3*fedback{f_ind,1}(indchantype(D1, 'MEGMAG','GOOD'),:)', 0)';

subtracted = data1-model;

% Then add three different "on times"
on_length = [5, 10, 50]*1e-3*D.fsample; % Time before feedback updates, in samples
subtracted = cat(3, subtracted, zeros([size(subtracted), length(on_length)]));
model_delayed = cat(3, model, zeros([size(model), length(on_length)]));
for ol = 1:length(on_length)
    model_tmp = model(:, 1:on_length(ol):end);
    F = griddedInterpolant(D1.time(1:on_length(ol):end), model_tmp', 'nearest');
    model_delayed(:,:,ol+1) = transpose(F(D1.time));
    subtracted(:,:,ol+1) = data1 - transpose(F(D1.time));
end
on_length = [1, on_length];

% Plot

figure;
t = tiledlayout(length(on_length), 2);
h = gobjects(length(on_length),2);

for ol = 1:length(on_length)
    % Plot original
    h(ol,1) = nexttile((ol-1)*2+1);
    hold(h(ol,1), 'on');
    plot(h(ol,1), D1.time, 1e-6*(data1(chanind_subtracted,:)-...
        mean(data1(chanind_subtracted,1:2*D.fsample),2)), 'LineWidth', 1.5);
    plot(h(ol,1), D1.time, 1e-6*(model_delayed(chanind_subtracted,:,ol)-...
        mean(model_delayed(chanind_subtracted,1:2*D.fsample,ol),2)), 'LineWidth', 1.5);


    % Plot subtracted
    h(ol,2) = nexttile((ol-1)*2+2);
    plot(h(ol,2), D1.time, 1e-6*subtracted(chanind_subtracted,:,ol), 'LineWidth', 1.5);

%     k = kurtosis(1e-6*subtracted(chanind_subtracted,:,ol));
%     text(h(ol,2), min(D.time)+0.1, 0.15, sprintf('Kurtosis = %.2f', k), 'FontSize', 16);

%     % Plot fedback
%     h(ol,3) = nexttile((ol-1)*3+3);
%     plot(h(ol,3), D2.time, 1e-6*(D2(chanind_D,:)-mean(D2(chanind_D,1:2*D.fsample),2)));
% 
    for r = 1:2
        grid(h(ol,r), 'on'); hold(h(ol,r), 'on'); box(h(ol,r), 'on');
        xlim(h(ol,r), [0 1])
        set(h(ol,r), 'FontSize', 14);
        ylim(h(ol,r), [-0.1 0.1])
    end
    if ol == 1
        ylabel(h(ol,1), 'Every sample', 'Rotation', 0, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
    else
        ylabel(h(ol,1), sprintf('%.f ms', on_length(ol)*1e3/D.fsample), 'Rotation', 0, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
    end
end
linkaxes(h,'xy');

xlabel(t,'Time (s)', 'FontSize', 18)
ylabel(t, 'B (nT)', 'FontSize', 18)

title(h(1,1), 'No Feedback');
title(h(1,2), 'Subtract Model')
% title(h(1,3), 'Rec. With Feedback')

% Move plots closer together
yticklabels(h(:,2:end),{})
xticklabels(h(1:end-1,:),{})
t.TileSpacing = 'compact';

% Legend
leg = legend(h(3,1), {'Recorded', 'Model'});
leg_pos = leg.Position;
set(leg, 'Position', [0.3336,0.0135,leg_pos(3:4)]);

% Figure size
set(gcf, 'color', 'w');
set(gcf, 'Position', [680, 304, 840, 793]);

export_fig(fullfile(resultsSavePath, 'delayed_model_time_series'), '-png', '-nocrop', '-painters')
%print(fullfile(resultsSavePath, 'delayed_model_time_series'), '-dpng', '-r300')

%% Plot

% chaninds = indchantype(D, 'MEGMAG', 'GOOD');
% chanind = chaninds(randi(length(chaninds)));
% 
% for f = 1:length(frequency)
%     for r = 1:length(run)
%         D = spm_eeg_load(sprintf('sub-001_ses-001_task-%sHz_run-00%s_meg.mat', frequency{f}, run{r}));
% 
%         figure; hold on; grid on; box on;
%         plot(D.time, D(chanind,:,1) - mean(D(chanind,:,1),2), ...
%             D.time, 1e3*fedback{f,r}(chanind,:) - mean(1e3*fedback{f,r}(chanind,:),2));
%         xlim([0 max(D.time)])
%         %ylim([-2e5 2e5])
% 
%         xlabel('Time (s)');
%         ylabel('B (fT)');
%         if r == 1
%             legend('Recorded', 'What would be fed back', 'location', 'best');
%         else
%             legend('Recorded', 'What was fed back', 'location', 'best');
%         end
% 
%     end
% end