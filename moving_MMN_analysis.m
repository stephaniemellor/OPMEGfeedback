%% Housekeeping
clear; close all;

addpath('D:\Documents\Software\spm12\FILversion');
spm('defaults', 'eeg');

addpath('D:\Documents\Scripts\opmeg_fieldmapping_code\HelperFunctions\icp')
addpath('D:\Documents\Software\FIL-OPMEG\optitrack');
addpath('D:\Documents\Software\FIL-OPMEG\OPM');
addpath('D:\Documents\Software\export_fig')

addpath('D:\Documents\Software\fieldtrip-lite-20180603\fieldtrip-20180603\external\brewermap')
colormap123 = colormap(flipud(brewermap(64,'RdBu')));

%% Config

% File save paths
resultsSavePath = '.\AEF\results';
analysedDataPath = '.\AEF\analysedData';
rawDataPath = '.\AEF';
cd(analysedDataPath);

% Block numbers for each condition. Row = participant/session
seated_no_feedback = [NaN, NaN; 2, 3; 1, 3];
seated_feedback = [NaN, NaN; 1, 4; 2, 4];
walking_no_feedback = [1, 4; 7, 8; 7, 8];
walking_feedback = [2, 3; 5, 6; 5, 6];

% Order of participants and sessions to analyse
process_order(1) = struct('sub', '001', 'ses', '001', 'stimTrig', 'NI-TRIG-4');
process_order(2) = struct('sub', '001', 'ses', '002', 'stimTrig', 'NI-TRIG-1');
process_order(3) = struct('sub', '002', 'ses', '001', 'stimTrig', 'NI-TRIG-1');

% Intialise cell arrays of meeg objects, optitrack data, feedback applied, 
% saturated frames and beep stop start datapoints for different
% participants
D = cell(8,length(process_order));
opti_data = cell(8,length(process_order));
feedback_applied = cell(8, length(process_order));
saturated_frames = cell(8, length(process_order));
beep_start_stop = cell(8, length(process_order));

% D with different preprocessing steps
D_pos_reg_no_HFC = cell(8,length(process_order));
D_pos_reg = cell(8,length(process_order));
D_no_HFC = cell(8,length(process_order));

% Merged D objects
cD = cell(4,length(process_order)); % rows go seated no feedback, seated feedback, walking no feedback, walking feedback
cD_pos_reg_no_HFC = cell(4,length(process_order));
cD_pos_reg = cell(4,length(process_order));
cD_no_HFC = cell(4,length(process_order));

%% Load OPM data

task = 'tone';
run_cell = {'001', '002', '003', '004', '005', '006', '007', '008'};

for ss = 1:length(process_order)

    sub = process_order(ss).sub;
    ses = process_order(ss).ses;
    sMRI = [rawDataPath, '\sub-', sub, '\sub-', sub, '.nii'];

    % Handle differences between participant 1 session 1 and other
    % recordings
    if ss == 1
        % Correct optitrack coordinate system. Rather than calibrating in 
        % Motive, a recording was made with the calibration square on the floor. 
        groundPlane = csv2mat_sm(fullfile([rawDataPath, '\sub-', sub, '\ses-', ses], 'GroundPlane.csv'));
        
        GroundMarker1_Camera = nanmean(groundPlane.markers.rigidbodymarkers(1).data(:,1:3), 1);
        GroundMarker2_Camera = nanmean(groundPlane.markers.rigidbodymarkers(2).data(:,1:3), 1);
        GroundMarker3_Camera = nanmean(groundPlane.markers.rigidbodymarkers(3).data(:,1:3), 1);
        groundPlane_Camera = [GroundMarker1_Camera', GroundMarker2_Camera', GroundMarker3_Camera'];
        
        % Get rotation and translation matrices to convert to room coords from
        % Camera
        [T, R] = calibrateGroundPlaneOnFloor(groundPlane_Camera);
        optiTrig = 'NI-TRIG-1'; 
        calibScale = 0.5;
    else
        T = [0;0;0];
        R = eye(3);
        optiTrig = 'NI-TRIG-7'; 
        calibScale = 1;
    end
    
    % Preprocess data
    run_list = [seated_no_feedback(ss,:), seated_feedback(ss,:), walking_no_feedback(ss,:), walking_feedback(ss,:)];
    for run_counter = 1:length(run_list)
        if ~isnan(run_list(run_counter))
            % Using default steps
            [D{run_counter, ss}, opti_data{run_counter, ss}, feedback_applied{run_counter, ss}, ...
                saturated_frames{run_counter, ss}, beep_start_stop{run_counter, ss}] =...
                    preprocess_feedback_ERF(rawDataPath, analysedDataPath, resultsSavePath, sMRI, sub, ses, task, run_cell{run_list(run_counter)}, ...
                    'optitrack_R', R, 'optitrack_T', T, 'optitrack_trigger', optiTrig, 'stim_trigger', process_order(ss).stimTrig, 'calibration_scaling', calibScale);

            % Preprocess data with no HFC
            D_no_HFC{run_counter, ss} =...
                preprocess_feedback_ERF(rawDataPath, analysedDataPath, resultsSavePath, sMRI, sub, ses, task, run_cell{run_list(run_counter)}, ...
                'optitrack_R', R, 'optitrack_T', T, 'optitrack_trigger', optiTrig, 'stim_trigger', process_order(ss).stimTrig, 'calibration_scaling', calibScale,...
                'HFC', false);

            % Preprocess data with movement regression
            D_pos_reg{run_counter, ss} =...
                preprocess_feedback_ERF(rawDataPath, analysedDataPath, resultsSavePath, sMRI, sub, ses, task,run_cell{run_list(run_counter)}, ...
                'optitrack_R', R, 'optitrack_T', T, 'optitrack_trigger', optiTrig, 'stim_trigger', process_order(ss).stimTrig, 'calibration_scaling', calibScale,...
                'position_regress', true);

            % Preprocess data with movement regression and no HFC
            D_pos_reg_no_HFC{run_counter, ss} =...
                preprocess_feedback_ERF(rawDataPath, analysedDataPath, resultsSavePath, sMRI, sub, ses, task, run_cell{run_list(run_counter)}, ...
                'optitrack_R', R, 'optitrack_T', T, 'optitrack_trigger', optiTrig, 'stim_trigger', process_order(ss).stimTrig, 'calibration_scaling', calibScale,...
                'position_regress', true, 'HFC', false);
        end
    end

    %% Merge feedback on/off runs

    for condition = 1:4
        if ~isempty(D{2*condition-1,ss}) && ~isempty(D{2*condition, ss})
            % Default preprocessing
            cD{condition, ss} = merge_feedback_Dobjects(D{2*condition-1, ss}, ...
                D{2*condition, ss});

            % No HFC
            cD_no_HFC{condition, ss} = merge_feedback_Dobjects(D_no_HFC{2*condition-1, ss}, ...
                D_no_HFC{2*condition, ss});

            % Position regression
            cD_pos_reg{condition, ss} = merge_feedback_Dobjects(D_pos_reg{2*condition-1, ss}, ...
                D_pos_reg{2*condition, ss});

            % Position regression, no HFC
            cD_pos_reg_no_HFC{condition, ss} = merge_feedback_Dobjects(D_pos_reg_no_HFC{2*condition-1, ss}, ...
                D_pos_reg_no_HFC{2*condition, ss});
        end
    end

    close all;
end

%% Plot median psd at each pre-processing step for each run

for ss = 1:length(process_order)
    for run_counter = 1:length(run_cell)
        if ~isempty(D{run_counter,ss})
            S = [];
            S.channels = D{run_counter,ss}.chanlabels(indchantype(D{run_counter,ss}, 'MEGMAG', 'GOOD'));
            S.plot = 0;
            S.triallength = 10e3;

            % Get start of base name
            Dfname = D{run_counter,ss}.fname;
            start_of_base_name = strfind(Dfname, 't_dsub');
    
            % Base data
            Dtmp = spm_eeg_load(Dfname(start_of_base_name:end));
            S.D = Dtmp;
            [psd{1}, freq] = spm_opm_psd(S);

            % HFC
            Dtmp = spm_eeg_load(['h', Dfname(start_of_base_name:end)]);
            S.D = Dtmp;
            [psd{2}, ~] = spm_opm_psd(S);
    
            % Filtering
            Dtmp = spm_eeg_load(['fffffh', Dfname(start_of_base_name:end)]);
            S.D = Dtmp;
            [psd{3}, ~] = spm_opm_psd(S);

            % Plot
            figure; grid on; hold on; box on;
            c = get(gca, 'colororder');
            for pp = 1:length(psd)
                err = 1.2533*squeeze(std(psd{pp},[],2)./sqrt(size(psd{pp},2)));
                po = squeeze(median(psd{pp},2));
                to_plot = [po-err; flipud(po+err)];
                to_plot(to_plot < 0) = min(to_plot(to_plot > 0));
                fill([freq';flipud(freq')], to_plot, c(pp,:),...
                    'LineStyle', 'none', 'FaceAlpha', 0.5);
                p(pp) = plot(freq,median(psd{pp},2),'LineWidth',1.5, 'color', c(pp,:));
            end
            plot([min(freq), max(freq)], [15 15], 'k--', 'LineWidth', 2);
            xlim([0 150])
            ylim([0.1, 1e6]);
            set(gca, 'YScale', 'log');
            legend(p, {'Original Data', 'HFC', 'Filtered'});
            xlabel('Frequency (Hz)')
            labY = ['$$PSD (fT\sqrt[-1]{Hz}$$)'];
            ylabel(labY,'interpreter','latex')
            set(gca, 'FontSize', 16);
            title(sprintf('Participant %s, Block %s', process_order(ss).sub, run_cell{run_counter}));
            set(gcf, 'color', 'w');
            
            export_fig(fullfile(resultsSavePath, sprintf('mpsd_processing_steps_participant_%s_session_%s_run_%s',...
                process_order(ss).sub, process_order(ss).ses, run_cell{run_counter})), '-png', '-painters')
        end
    end
end

%% Repeat with movement regression

for ss = 1:length(process_order)
    for run_counter = 1:length(run_cell)
        if ~isempty(D{run_counter,ss})
            S = [];
            S.channels = D_pos_reg{run_counter, ss}.chanlabels(indchantype(D_pos_reg{run_counter, ss}, 'MEGMAG', 'GOOD'));
            S.plot = 0;
            S.triallength = 10e3;

            % Get start of base name
            Dfname = D_pos_reg{run_counter, ss}.fname;
            start_of_base_name = strfind(Dfname, 't_dsub');

            % Base data
            Dtmp = spm_eeg_load(Dfname(start_of_base_name:end));
            S.D = Dtmp;
            [psd{1}, freq] = spm_opm_psd(S);

            % Movement regression
            Dtmp = spm_eeg_load(['mr_', Dfname(start_of_base_name:end)]);
            S.D = Dtmp;
            [psd{2}, ~] = spm_opm_psd(S);

            % HFC
            Dtmp = spm_eeg_load(['hmr_', Dfname(start_of_base_name:end)]);
            S.D = Dtmp;
            [psd{3}, ~] = spm_opm_psd(S);
    
            % Filtering
            Dtmp = spm_eeg_load(['fffffhmr_', Dfname(start_of_base_name:end)]);
            S.D = Dtmp;
            [psd{4}, ~] = spm_opm_psd(S);

            % Plot
            figure; grid on; hold on; box on;
            c = get(gca, 'colororder');
            for pp = 1:length(psd)
                err = 1.2533*squeeze(std(psd{pp},[],2)./sqrt(size(psd{pp},2)));
                po = squeeze(median(psd{pp},2));
                to_plot = [po-err; flipud(po+err)];
                to_plot(to_plot < 0) = min(to_plot(to_plot > 0));
                fill([freq';flipud(freq')], to_plot, c(pp,:),...
                    'LineStyle', 'none', 'FaceAlpha', 0.5);
                p(pp) = plot(freq,median(psd{pp},2),'LineWidth',1.5, 'color', c(pp,:));
            end
            plot([min(freq), max(freq)], [15 15], 'k--', 'LineWidth', 2);
            xlim([0 10])
            ylim([0.1, 1e6]);
            set(gca, 'YScale', 'log');
            legend(p, {'Original Data', 'Movement Regression', 'HFC', 'Filtered'}, 'location', 'best');
            xlabel('Frequency (Hz)')
            labY = ['$$PSD (fT\sqrt[-1]{Hz}$$)'];
            ylabel(labY,'interpreter','latex')
            set(gca, 'FontSize', 16);
            title(sprintf('Participant %s, Block %s', process_order(ss).sub, run_cell{run_counter}));
            set(gcf, 'color', 'w');
        
            export_fig(fullfile(resultsSavePath, sprintf('mpsd_processing_steps_mov_reg_participant_%s_session_%s_run_%s',...
                process_order(ss).sub, process_order(ss).ses, run_cell{run_counter})), '-png', '-painters')
        end
    end
end

%% Sensor level AEF

% Include bad trials?
incl_bad_trials = 0;

% Colour lines by left/right and distance from auditory cortex?
split_by_lr = 0;
colour_by_distance = 1;

for ss = 1:length(process_order)

    table_of_info = fullfile(rawDataPath, sprintf('sub-%s', process_order(ss).sub), sprintf('ses-%s', process_order(ss).ses), 'table_of_info.csv');
    slot_to_sens = fullfile(rawDataPath, sprintf('sub-%s', process_order(ss).sub), sprintf('ses-%s', process_order(ss).ses), 'slot2sens.csv');

    % Load layout
    load(fullfile(rawDataPath, sprintf('sub-%s', process_order(ss).sub), sprintf('ses-%s', process_order(ss).ses), 'layout.mat'));

    % Run for combined and separate runs and different preprocessing steps
    for DD = {cD, cD_no_HFC, cD_pos_reg, cD_pos_reg_no_HFC}
        DD = DD{1};

        % Average - each run individually
        for run_counter = 1:size(DD,1)

            Dtmp = DD{run_counter, ss};
            if ~isempty(Dtmp)
                data = ftraw(Dtmp);
            
                % Select equal number of trials per condition
                rng('default');
                rng(95);
                if size(Dtmp, 3) < 1000 % If runs have not been combined
                    trials = randperm(size(Dtmp,3), 540);
                else
                    trials = randperm(size(Dtmp,3), 1120);
                end
                
                if ~incl_bad_trials
                    % Just select good trials
                    cfg = [];
                    cfg.trials = intersect(indtrial(Dtmp, 'Cond1', 'GOOD'), trials);
                    data = ft_selectdata(cfg, data);
                else
                    cfg = [];
                    cfg.trials = trials;
                    data = ft_selectdata(cfg, data);
                end
            
                % Average
                cfg = [];
                cfg.channel = data.grad.label;
                avdata = ft_timelockanalysis(cfg, data);
            
                % T-stat
                epoched_dataset = [];
                for i = 1:length(data.trial)
                    epoched_dataset(:,:,i) = data.trial{1,i}(indchannel(Dtmp, avdata.label),:);
                end
                SE              = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
                avdata.t_value  = avdata.avg./SE;
    
                % Determine which sensors are in each hemisphere, and which are
                % closest to the auditory cortices
                clear chans
                if split_by_lr || colour_by_distance
        
                    % Get sensor positions in MNI coordinates
                    chanpos = Dtmp.sensors('MEG').chanpos;
                    chanpos_MNI = Dtmp.inv{1}.datareg.toMNI*cat(1, chanpos', ones(1, size(chanpos,1)));
                    chanpos_MNI = chanpos_MNI(1:3,:)';
        
                    % Separate into left and right sensors based on x
                    % position value
                    lchans = chanpos_MNI(:,1) < 0;
        
                    % Get distance between each channel and MNI auditory
                    % cortex coordinates
                    aud = [-48 -22 4; 48 -22 4];
                    ldist = sqrt(sum((chanpos_MNI - aud(1,:)).^2, 2));
                    rdist = sqrt(sum((chanpos_MNI - aud(2,:)).^2, 2));
        
                    [~, linds] = sort(ldist, 'descend');
                    linds = intersect(linds, find(lchans), 'stable');
                    [~, rinds] = sort(rdist, 'descend');
                    rinds = intersect(rinds, find(~lchans), 'stable');
        
                    % Create lists of channel names, split by left vs right
                    % and ordered by distance from auditory cortex
                    leftchans = Dtmp.sensors('MEG').label(linds);
                    rightchans = Dtmp.sensors('MEG').label(rinds);
                    
                    % Set up channels
                    if split_by_lr
                        plots = 2;
                    else
                        plots = 1;
                    end
                    chans{1} = zeros(length(leftchans),1);
                    chans{2} = zeros(length(rightchans),1);
                    for cc = 1:length(leftchans)
                        chans{1}(cc) = find(contains(avdata.label, leftchans{cc}));
                    end
                    for cc = 1:length(rightchans)
                        chans{2}(cc) = find(contains(avdata.label, rightchans{cc}));
                    end
        
                    % Assuming leftchans and rightchans have been ordred by
                    % increasing distance from primary auditory cortex, 
                    % (i.e. first is furthest), should colour them in 
                    % darkening shades as nearer to the auditory cortex
                    color{1} = interp1(linspace(1/min(ldist(linds)).^3, 1/max(ldist(linds)).^3, 10), linspace(0, 0.7, 10), ...
                        1./ldist(linds).^3);
                    color{2} = interp1(linspace(1/min(rdist(rinds)).^3, 1/max(rdist(rinds)).^3, 10), linspace(0, 0.7, 10), ...
                        1./rdist(rinds).^3);
                    for pp = 1:2
                        color{pp} = repmat(color{pp}, 1, 3);
                    end
                else
                    plots = 1;
                    chans{1} = avdata.label;
                    color{1} = repmat([1,1,1], length(chans{1}), 1);
                end
    
                if ss == 1
                    figure; hold on
                    if split_by_lr
                        t = tiledlayout(4,7);
                        set(gcf, 'Position', [442   549   981   377]);
                    else
                        t = tiledlayout(2,7);
                        set(gcf, 'Position', [442   549   981   377]);
                    end
    
                    % Topoplot - 90 ms
                    cfg = [];
                    cfg.layout    = layout;
                    cfg.zlim      = [-4 4];
                    cfg.colormap  = colormap123;
                    cfg.xlim = [0.085, 0.095];
                    cfg.comment = 'no';
                    cfg.parameter = 't_value';
                    cfg.interplimits = 'electrodes';
        
                    if split_by_lr
                        % Plot t-value topoplot below
                        nexttile(18, [2,2]);
                        cfg.parameter = 't_value';
                        cfg.colorbartext = 't-value';
                        cfg.figure = gca;
                        ft_topoplotER(cfg, avdata)
                        set(gca, 'FontSize', 14)
                
                        % Plot average topoplot at top
                        nexttile(4, [2,2]);
                        cfg.parameter = 'avg';
                        cfg.colorbartext = 'B (fT)';
                        cfg.figure = gca;
                    else
                        % If only 2 rows of figure, just plot t-value topoplot
                        nexttile(4, [2,2]);
                        cfg.figure = gca;
                    end
                    ft_topoplotER(cfg, avdata)
                    cb = colorbar;
                    set(cb, 'Position', [0.669314157866628,0.324668435687833,0.017471436642701,0.410079574922246]);
                    title(cb, 't-value');
                    set(gca, 'FontSize', 14)
                    title('90 ms')
    
                    % Topoplot - 110 ms
                    cfg.xlim = [0.105, 0.115];
                    if split_by_lr
                        % Plot t-value topoplot below
                        nexttile(20, [2,2]);
                        cfg.parameter = 't_value';
                        cfg.colorbartext = 't-value';
                        cfg.figure = gca;
                        ft_topoplotER(cfg, avdata)
                        set(gca, 'FontSize', 14)
                
                        % Plot average topoplot at top
                        nexttile(6, [2,2]);
                        cfg.parameter = 'avg';
                        cfg.colorbartext = 'B (fT)';
                        cfg.figure = gca;
                    else
                        % If only 2 rows of figure, just plot t-value topoplot
                        nexttile(6, [2,2]);
                        cfg.figure = gca;
                    end
                    ft_topoplotER(cfg, avdata)
                    cb = colorbar;
                    set(cb, 'Position', [0.905056595977079,0.324668435687833,0.017471436642701,0.410079574922246]);
                    title(cb, 't-value');
                    set(gca, 'FontSize', 14)
                    title('110 ms')
    
                    % Time series
                    for pp = 1:plots
                        % Plot average first
                        nexttile((pp-1)*5+1, [1,3]);
                        hold on; box on; grid on;
    
                        for cc = 1:length(chans{pp})
                            plot(avdata.time, avdata.avg(chans{pp}(cc),:), '-', 'Color', color{pp}(cc,:));
                        end
                        ylim([-300, 300])
                        xlim([-0.1, 0.4]);
                        xlabel('Time (s)', 'FontSize', 16);
                        ylabel('B (fT)', 'FontSize', 16);
                        set(gca, 'FontSize', 14);
            
                        % Then plot t-value below
                        nexttile((pp-1+plots)*7+1, [1,3]);
                        hold on; box on; grid on;
            
                        for cc = 1:length(chans{pp})
                            plot(avdata.time, avdata.t_value(chans{pp}(cc),:), '-', 'Color', color{pp}(cc,:));
                        end
                        ylim([-7 7])
                        xlim([-0.1, 0.4]);
                        xlabel('Time (s)', 'FontSize', 16);
                        ylabel('t-value', 'FontSize', 16);
                        set(gca, 'FontSize', 14)
                    end
                else
                    figure; hold on
                    if split_by_lr
                        t = tiledlayout(4,5);
                        set(gcf, 'Position', [442   549   710   277]);
                    else
                        t = tiledlayout(2,5);
                        set(gcf, 'Position', [442   549   710   277]);
                    end
            
                    % Topoplot - 100 ms
                    cfg = [];
                    cfg.layout    = layout;
                    cfg.colorbar  = 'EastOutside';
                    cfg.zlim      = 'maxabs';
                    cfg.colormap  = colormap123;
                    cfg.xlim = [0.095, 0.105];
                    cfg.comment = 'no';
                    cfg.parameter = 't_value';
                    cfg.colorbartext = 't-value';
                    cfg.interplimits = 'electrodes';
                    if run_counter < 3      % use higher colour scale limit if seated
                        cfg.zlim = [-13 13];
                    else
                        cfg.zlim = [-8 8];
                    end
        
                    if split_by_lr
                        % Plot t-value topoplot below
                        nexttile(14, [2,2]);
                        cfg.parameter = 't_value';
                        cfg.colorbartext = 't-value';
                        cfg.figure = gca;
                        ft_topoplotER(cfg, avdata)
                        set(gca, 'FontSize', 14)
            
                        % Plot average topoplot at top
                        nexttile(4, [2,2]);
                        cfg.parameter = 'avg';
                        cfg.colorbartext = 'B (fT)';
                        cfg.zlim = 'maxabs';
                        cfg.figure = gca;
                    else
                        % If only 2 rows of figure, just plot t-value topoplot
                        nexttile(4, [2,2]);
                        cfg.figure = gca;
                    end
                    ft_topoplotER(cfg, avdata)
                    set(gca, 'FontSize', 14)
        
                    % Time series
                    for pp = 1:plots
            
                        % Plot average first
                        nexttile((pp-1)*5+1, [1,3]);
                        hold on; box on; grid on;
            
                        for cc = 1:length(chans{pp})
                            plot(avdata.time, avdata.avg(chans{pp}(cc),:), '-', 'Color', color{pp}(cc,:));
                        end
                        yl = ylim;
                        ylim([-max(abs(yl)), max(abs(yl))])
                        xlim([-0.1, 0.4]);
                        xlabel('Time (s)', 'FontSize', 16);
                        ylabel('B (fT)', 'FontSize', 16);
                        set(gca, 'FontSize', 14);
            
                        % Then plot t-value below
                        nexttile((pp-1+plots)*5+1, [1,3]);
                        hold on; box on; grid on;
            
                        for cc = 1:length(chans{pp})
                            plot(avdata.time, avdata.t_value(chans{pp}(cc),:), '-', 'Color', color{pp}(cc,:));
                        end
                        if run_counter < 3
                            ylim([-13 13]);
                        else
                            ylim([-8 8]);
                        end
                        xlim([-0.1, 0.4]);
                        xlabel('Time (s)', 'FontSize', 16);
                        ylabel('t-value', 'FontSize', 16);
                        set(gca, 'FontSize', 14)
                    end
                end
        
                % Compact subplot spacing
                t.TileSpacing = 'compact';
        
                % Save image
                if incl_bad_trials
                    bad_str = 'including_sat_trials';
                else
                    bad_str = 'excluding_sat_trials';
                end
        
                set(gcf, 'color', 'w');
                export_fig([resultsSavePath, '\', sprintf('AEF_%s_filename_%s', bad_str, Dtmp.fname)], '-png', '-painters')
            end
        end
    end
end

%% Get number of saturated trials for each participant and walking condition (supplementary table)

for ss = 1:length(process_order)

    % Load layout
    load(fullfile(rawDataPath, sprintf('sub-%s', process_order(ss).sub), sprintf('ses-%s', process_order(ss).ses), 'layout.mat'));

    % Get trl matrix

    % No feedback
    for run_counter = 5:6
        Dfname = D{run_counter, ss}.fname;
        D{run_counter, ss} = spm_eeg_load(fullfile(D{run_counter, ss}.path, Dfname(3:end)));
        S = [];
        S.D = D{run_counter, ss};
        S.timewin = [-200 500];
        S.condlabels = {'tone'};
        S.triggerChannels = {process_order(ss).stimTrig};
        [D{run_counter, ss}, trl{run_counter}] = spm_opm_epoch_trigger(S);
    end

    % With feedback
    for run_counter = 7:8
        Dfname = D{run_counter, ss}.fname;
        D{run_counter, ss} = spm_eeg_load(fullfile(D{run_counter, ss}.path, Dfname(3:end)));
        S = [];
        S.D = D{run_counter, ss};
        S.timewin = [-200 500];
        S.condlabels = {'tone'};
        S.triggerChannels = {process_order(ss).stimTrig};
        [D{run_counter,ss}, trl{run_counter}] = spm_opm_epoch_trigger(S);
    end

    % Find which trial each saturated frame is in
    saturated_trials = cell(4, length(saturated_frames{8, ss}));
    for chan = 1:size(saturated_trials, 2)
        for run_counter = 5:8
            if ~isempty(saturated_frames{run_counter,ss}{chan})
                [~, col] = find((saturated_frames{run_counter,ss}{chan} >= trl{run_counter}(:,1)').*...
                    (saturated_frames{run_counter,ss}{chan} <= trl{run_counter}(:,2)'));
                saturated_trials{run_counter-4, chan} = unique(col);
            end
        end
    end

    % Combine runs
    % No feedback
    saturated_trials(2,:) = cellfun(@(x)x+size(trl{5},1), saturated_trials(2,:), 'UniformOutput',false);
    saturated_trials(4,:) = cellfun(@(x)x+size(trl{7},1), saturated_trials(4,:), 'UniformOutput',false);
    for chan = 1:size(saturated_trials,2)
        saturated_trials{1, chan} = cat(1, saturated_trials{1,chan}, saturated_trials{2,chan});
        saturated_trials{3, chan} = cat(1, saturated_trials{3,chan}, saturated_trials{4,chan});
    end
    saturated_trials_walking_no_feedback = saturated_trials(1,:);
    saturated_trials_walking_feedback = saturated_trials(3,:);
    clear saturated_trials

    % Only keep trials within the 1120 randomly selected trials
    rng('default');
    rng(95);
    trials = randperm(size(cD{3,ss},3), 1120);
    saturated_trials_walking_no_feedback = cellfun(@(x)intersect(x, trials), saturated_trials_walking_no_feedback, 'UniformOutput',false);
    rng('default');
    rng(95);
    trials = randperm(size(cD{4,ss},3), 1120);
    saturated_trials_walking_feedback = cellfun(@(x)intersect(x, trials), saturated_trials_walking_feedback, 'UniformOutput',false);

    % Save number of saturated trials for each channel
    nsat_trials_no_feedback{ss} = cell2mat(cellfun(@length, saturated_trials_walking_no_feedback, 'UniformOutput',false));
    nsat_trials_feedback{ss} = cell2mat(cellfun(@length, saturated_trials_walking_feedback, 'UniformOutput',false));

    % Plot spatial map of which channels saturate the most
    satchans = cD{3,ss}.chanlabels(indchantype(cD{3,ss}, 'MEGMAG', 'GOOD'));
    layout_chans = layout.label;

    % Y channels
    [seldat, sellay] = match_str(satchans, layout_chans);

    % No feedback
    figure; hold on;
    scatter(layout.pos(sellay,1), layout.pos(sellay,2), 300, 100*nsat_trials_no_feedback{ss}(seldat)/1120, 'filled', 'square');
    ft_plot_layout(layout, 'chanindx', [], 'box', 'no', 'label', 'no', 'pointsymbol', '.', 'pointcolor', 'w', 'pointsize', 8);
    daspect([1 1 1])
    colormap(colormap123);
    caxis([-35 35]);
    colorbar;
    set(gcf, 'color','w');
    set(gca, 'Visible', 'off');
    export_fig(fullfile(resultsSavePath, sprintf('sub-%s_ses-%s_percentage_saturated_trials_no_feedback_Ychans', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters');

    % Feedback on
    figure; hold on;
    scatter(layout.pos(sellay,1), layout.pos(sellay,2), 300, 100*nsat_trials_feedback{ss}(seldat)/1120, 'filled', 'square');
    ft_plot_layout(layout, 'chanindx', [], 'box', 'no', 'label', 'no', 'pointsymbol', '.', 'pointcolor', 'w', 'pointsize', 8);
    daspect([1 1 1])
    colormap(colormap123);
    caxis([-35 35]);
    colorbar;
    set(gcf, 'color','w');
    set(gca, 'Visible', 'off');
    export_fig(fullfile(resultsSavePath, sprintf('sub-%s_ses-%s_percentage_saturated_trials_feedback_Ychans', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters');

    % Z channels
    layout_chans(endsWith(layout_chans, '-Y')) = cellfun(@(x)strrep(x, '-Y', '-Z'), layout_chans(endsWith(layout_chans, '-Y')), 'UniformOutput',false);
    [seldat, sellay] = match_str(satchans, layout_chans);

    % No feedback
    figure; hold on;
    scatter(layout.pos(sellay,1), layout.pos(sellay,2), 300, 100*nsat_trials_no_feedback{ss}(seldat)/1120, 'filled', 'square');
    ft_plot_layout(layout, 'chanindx', [], 'box', 'no', 'label', 'no', 'pointsymbol', '.', 'pointcolor', 'w', 'pointsize', 8);
    daspect([1 1 1])
    colormap(colormap123);
    caxis([-35 35]);
    colorbar;
    set(gcf, 'color','w');
    set(gca, 'Visible', 'off');
    export_fig(fullfile(resultsSavePath, sprintf('sub-%s_ses-%s_percentage_saturated_trials_no_feedback_Zchans', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters');

    % Feedback on
    figure; hold on;
    scatter(layout.pos(sellay,1), layout.pos(sellay,2), 300, 100*nsat_trials_feedback{ss}(seldat)/1120, 'filled', 'square');
    ft_plot_layout(layout, 'chanindx', [], 'box', 'no', 'label', 'no', 'pointsymbol', '.', 'pointcolor', 'w', 'pointsize', 8);
    daspect([1 1 1])
    colormap(colormap123);
    caxis([-35 35]);
    colorbar;
    set(gcf, 'color','w');
    set(gca, 'Visible', 'off');
    export_fig(fullfile(resultsSavePath, sprintf('sub-%s_ses-%s_percentage_saturated_trials_feedback_Zchans', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters');
end

%% Look at decrease in SNR over time - check that time since calibration isn't having too big of an impact

clear percent_increase_walking percent_decrease_seated

colormap_rb = colormap(flipud(brewermap(4,'RdBu')));
colormap_pg = colormap(flipud(brewermap(4,'PRGn')));

for ss = 2:3

    C = [colormap_rb([2,1,3,4],:); colormap_pg([2,1,3,4],:)];

    figure; 
    h1 = subplot(1,1,1); hold on; grid on; box on;

    % Find sensor with maximum t-value in feedback off seated trials at 100 ms
    megchans = indchannel(cD{1,ss}, intersect(cD{1,ss}.sensors('MEG').label, ...
        cD{1,ss}.chanlabels(indchantype(cD{1,ss}, 'MEGMAG', 'GOOD'))));
    SE = std(cD{1,ss}(megchans,:,:),[],3)/sqrt(size(cD{1,ss},3));
    t_value  = mean(cD{1,ss}(megchans,:,:), 3)./SE;
    toi = logical((cD{1,ss}.time < 105*1e-3).*(cD{1,ss}.time > 95*1e-3));
    [~, peak_chan] = max(max(abs(t_value(:, toi)), [], 2));

    % Look at reduction in SNR between feedback on seated and feedback off
    % seated
    SNR = zeros(1,4);
    for run_counter = 1:4
        % Just select 1120 trials
        full_trial_list = indtrial(cD{run_counter,ss}, cD{run_counter,ss}.condlist{1}, 'GOOD');
        rng('default');
        rng(95);
        trials = randperm(size(cD{run_counter,ss},3), 1120);
        trials = intersect(full_trial_list, trials);

        SE = std(cD{run_counter,ss}(megchans,:,trials),[],3)/sqrt(length(trials));
        t_value  = mean(cD{run_counter,ss}(megchans,:,trials), 3)./SE;
        toi = logical((cD{run_counter,ss}.time < 105*1e-3).*(cD{run_counter,ss}.time > 95*1e-3));
        SNR(run_counter) = max(max(abs(t_value(:, toi)), [], 2));
        SNR(run_counter)  = SNR(run_counter).^2;
    end
    percent_decrease_seated(ss-1) = 100*(SNR(1) - SNR(2))/SNR(1);

    fprintf('Seated, Participant %s, SNR no feedback %.2f, SNR feedback %.2f, Percent Decrease from feedback = %.2f\n', ...
        process_order(ss).sub, SNR(1), SNR(2), percent_decrease_seated(ss-1));

    % Look at increases in SNR between feedback on walking and feedback off
    % walking
    percent_increase_walking(ss-1) = 100*(SNR(4) - SNR(3))/SNR(3);

    fprintf('Walking, Participant %s, SNR no feedback %.2f, SNR feedback %.2f, Percent Increase from feedback = %.2f\n', ...
        process_order(ss).sub, SNR(3), SNR(4), percent_increase_walking(ss-1));

    % AEF, peak chan, all runs
    clear trials
    for run_counter = 1:length(run_cell)
        % Just select 540 trials
        full_trial_list = indtrial(D{run_counter,ss}, D{run_counter,ss}.condlist{1}, 'GOOD');
        rng('default');
        rng(95);
        trials{run_counter} = randperm(size(D{run_counter,ss},3), 540);
        trials{run_counter} = intersect(full_trial_list, trials{run_counter});

        SE = std(D{run_counter,ss}(peak_chan,:,trials{run_counter}),[],3)/sqrt(length(trials{run_counter}));
        t_value  = mean(D{run_counter,ss}(peak_chan,:,trials{run_counter}), 3)./SE;
        plot(h1, D{run_counter,ss}.time, t_value, 'color', C(run_counter,:));
    end
    xlabel(h1, 'Time (s)');
    ylabel(h1, 't-value');
    xlim(h1, [-0.1 0.4]);
    legend(sprintf('Block %.f', seated_no_feedback(ss,1)), ...
        sprintf('Block %.f', seated_no_feedback(ss,2)), ...
        sprintf('Block %.f', seated_feedback(ss,1)), ...
        sprintf('Block %.f', seated_feedback(ss,2)), ...
        sprintf('Block %.f', walking_no_feedback(ss,1)), ...
        sprintf('Block %.f', walking_no_feedback(ss,2)), ...
        sprintf('Block %.f', walking_feedback(ss,1)), ...
        sprintf('Block %.f', walking_feedback(ss,2)), 'location', 'southeast');
    set(gcf,'color', 'w');
    set(gca, 'FontSize', 14);
    yl = ylim;
    ylim([-max(abs(yl)), max(abs(yl))]);
    export_fig(fullfile(resultsSavePath, sprintf('sub-%s_ses-%s_sensor_level_t_value_over_runs', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')

    % Only looking at seated runs - change in SNR over time
    figure; hold on; grid on; box on;
    for run_counter = 1:4
        % get mean t value for each 180 trial section of each run
        res = D{run_counter, ss}(megchans,:,trials{run_counter});
        res = reshape(res, size(res,1), size(res,2), 180, 3);
        tres = squeeze(mean(res, 3)./(std(res,[],3)/sqrt(size(res,3))));
        peakt = abs(mean(tres(peak_chan, toi, :), 2));
        peakt_err = std(tres(peak_chan, toi,:), [], 2)./sqrt(sum(toi));
        errorbar(squeeze(peakt), squeeze(peakt_err), '.-')
    end
    xlabel('180 Trial Set')
    ylabel('Mean t-value')
    legend(sprintf('Block %.f, no feedback', seated_no_feedback(ss,1)), ...
        sprintf('Block %.f, no feedback', seated_no_feedback(ss,2)), ...
        sprintf('Block %.f, feedback', seated_feedback(ss,1)), ...
        sprintf('Block %.f, feedback', seated_feedback(ss,2)), 'location', 'best');

    set(gcf, 'color', 'w');
    set(gca, 'FontSize', 14);
    export_fig(fullfile(resultsSavePath, sprintf('sub-%s_ses-%s_sensor_level_t_value_over_time', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')
    clear trials
end

fprintf('Average percent decrease in SNR = %.2f\n', mean(percent_decrease_seated))
fprintf('Average percent increase in SNR = %.2f\n', mean(percent_increase_walking))

clear tres res SE t_value peak_chan

%% Source localisation - Region of interest

% Get lead fields
M = gifti(fullfile(spm('dir'), 'canonical', 'cortex_5124.surf.gii'));
pos = [-48 -22 4; 48 -22 4]; % Auditory cortices
[~, leftind] = min(sqrt(sum((M.vertices - pos(1,:)).^2,2)));
[~, rightind] = min(sqrt(sum((M.vertices - pos(2,:)).^2,2)));
Dlist = {cD, cD_no_HFC, cD_pos_reg, cD_pos_reg_no_HFC};
Dlistnames = {'cD', 'cD_no_HFC', 'cD_pos_reg', 'cD_pos_reg_no_HFC'};

for ss = 1:length(process_order)
    Dcounter = 0;
    for DD = Dlist
        Dcounter = Dcounter + 1;
        Dtmp = DD{1};

        L{ss} = full(spm_eeg_lgainmat(Dtmp{4,ss},[leftind, rightind]));
        figure; 
        for ii = 1:2
            for jj = 1:2
                h(ii, jj) = subplot(2,2,(ii-1)*2+jj); hold on; grid on; box on;
            end
        end

        sigt = zeros(1, 4);
        for run_counter = 1:4

            if ~isempty(Dtmp{run_counter,ss})
            
                % Just select 1120 trials
                full_trial_list = indtrial(Dtmp{run_counter,ss}, 'Cond1', 'GOOD');
                rng('default');
                rng(95);
                trials = randperm(size(Dtmp{run_counter,ss},3), 1120);
                trials = intersect(full_trial_list, trials);
                fprintf('Participant %s, condition %.f, %.f trials\n', process_order(ss).sub, run_counter, length(trials));
        
                % Generate data
                for tt = 1:length(trials)
                    X{ss, run_counter}(:,:,tt) = pinv(L{ss})*Dtmp{run_counter,ss}(indchantype(Dtmp{run_counter,ss},'MEGMAG','GOOD'),:,...
                        trials(tt));
                end
                
                for ii = 1:2
                    plot(h(ii,1), -200:500, mean(X{ss, run_counter}(ii,:,:), 3), 'LineWidth', 2);
                    plot(h(ii,2), -200:500, mean(X{ss, run_counter}(ii,:,:), 3)./(std(X{ss, run_counter}(ii,:,:), [], 3)/sqrt(length(trials))), 'LineWidth', 2);
                end
        
                % Get significant t value
                df = length(trials)-1;
                alpha = 0.025/size(X{ss, run_counter}, 2);
                sigt(run_counter) = tinv(1-alpha, df);
            end
        end
        for ii = 1:2
            plot(h(ii,2), [-200 500], [max(sigt), max(sigt)], 'k--', [-200 500], [-max(sigt), -max(sigt)], 'k--', 'LineWidth', 2);
            for jj = 1:2
                xlabel(h(ii,jj), 'Time (ms)');
                xlim(h(ii,jj), [-200 500]);
                yl = ylim(h(ii,jj));
                ylim(h(ii,jj), [-max(abs(yl)), max(abs(yl))]);
            end
            ylabel(h(ii,1), 'Source Amplitude (nAm)');
            ylabel(h(ii,2), 't-value');
        end
        set(gcf, 'Position', [401, 334, 967, 701])
        set(gcf, 'color', 'w');
        export_fig(fullfile(resultsSavePath, sprintf('source_level_AEF_participant_%s_session_%s_%s', ...
            process_order(ss).sub, process_order(ss).ses, Dlistnames{Dcounter}(3:end))), '-png', '-painters')
    end
    % Save one with legend
    if ss == 1
        legend(h(2,2), {'Walking, no feedback', 'Walking, feedback', 'Bonferroni threshold'}, 'location', 'bestoutside');
    else
        legend(h(2,2), {'Seated, no feedback', 'Seated, feedback', 'Walking, no feedback', 'Walking, feedback', 'Bonferroni threshold'}, 'location', 'bestoutside');
    end
    export_fig(fullfile(resultsSavePath, sprintf('source_level_AEF_participant_%s_session_%s_legend', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')
end

clear L Dlist Dlistnames X M

%% Source level change in SNR with feedback

for ss = 1:length(process_order)
    
    L{ss} = full(spm_eeg_lgainmat(cD{4,ss},[leftind, rightind]));
    SNR = zeros(4, 2);

    for run_counter = 1:4
            
        if ~isempty(cD{run_counter,ss})
            % Just select 1120 trials
            full_trial_list = indtrial(cD{run_counter,ss}, cD{run_counter, ss}.condlist{1}, 'GOOD');
            rng('default');
            rng(95);
            trials = randperm(size(cD{run_counter,ss},3), 1120);
            trials = intersect(full_trial_list, trials);
        
            % Generate data
            X{ss, run_counter} = zeros(2, size(cD{run_counter,ss},2), length(trials));
            for tt = 1:length(trials)
                X{ss, run_counter}(:,:,tt) = pinv(L{ss})*cD{run_counter,ss}(indchantype(cD{4,ss},'MEGMAG','GOOD'),:,...
                    trials(tt));
            end
    
            % Calculate t-stat
            SE = std(X{ss, run_counter},[],3)/sqrt(size(X{ss, run_counter},3));
            t_value  = mean(X{ss, run_counter}, 3)./SE;
            if ss == 1
                toi = logical((cD{run_counter,ss}.time < 95*1e-3).*(cD{run_counter,ss}.time > 85*1e-3));
            else
                toi = logical((cD{run_counter,ss}.time < 105*1e-3).*(cD{run_counter,ss}.time > 95*1e-3));
            end
            SNR(run_counter,:) = max(abs(t_value(:, toi)), [], 2)';
            SNR(run_counter,:) = SNR(run_counter,:).^2;
        else
            SNR(run_counter,:) = NaN;
        end
    end
    percent_decrease_seated(ss,:) = 100*(SNR(1,:) - SNR(2,:))./SNR(1,:);

    for ac = 1:2
        fprintf('Seated, Participant %s, Session %s, Auditory Cortex %.f, SNR no feedback %.2f, SNR feedback %.2f, Percent Decrease = %.2f\n', ...
            process_order(ss).sub, process_order(ss).ses, ac, SNR(1,ac), SNR(2,ac), percent_decrease_seated(ss,ac));
    end

    percent_increase_walking(ss,:) = 100*(SNR(4,:) - SNR(3,:))./SNR(3,:);

    for ac = 1:2
        fprintf('Walking, Participant %s, Session %s, Auditory Cortex %.f, SNR no feedback %.2f, SNR feedback %.2f, Percent Increase from feedback = %.2f\n', ...
            process_order(ss).sub, process_order(ss).ses, ac, SNR(3,ac), SNR(4,ac), percent_increase_walking(ss,ac));
    end
end

fprintf('Average percent decrease in SNR = %.2f\n', mean(percent_decrease_seated, 'omitnan'))
fprintf('Average percent increase in SNR = %.2f\n', mean(percent_increase_walking(2:3,:), 'omitnan'))

clear L X leftind rightind

%% Source localisation - Minimum Norm

incl_bad_trials = 0;

% Load and interpolate atlas to spm canonical mesh
atlas = ft_read_atlas('D:\Documents\Software\fieldtrip-master\template\atlas\aal\ROI_MNI_V4.nii');
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
mesh_tissue = ft_read_headshape(fullfile(spm('dir'), 'canonical', 'cortex_5124.surf.gii'));
mesh_tissue = ft_sourceinterpolate(cfg, atlas, mesh_tissue);

% Inflate mesh for displaying source space
M_original = gifti(fullfile(spm('dir'), 'canonical', 'cortex_5124.surf.gii'));
M = spm_mesh_inflate(spm_mesh_inflate(M_original));
mesh = ft_read_headshape(fullfile(spm('dir'), 'canonical', 'cortex_5124.surf.gii'));
mesh.pos = M.vertices;
clear M

% Get smoothing kernel
[~,Di] = spm_mesh_neighbours(M_original,1);
muNeighbour = mean(mean(Di));
n = round((8/muNeighbour)^2);

% Create arrays to save t-stats, significance thresholds and tissues of peak t values
% Dlist = {D, D_no_HFC, D_pos_reg, D_pos_reg_no_HFC, cD, cD_no_HFC, cD_pos_reg, cD_pos_reg_no_HFC};
Dlist = {cD, cD_no_HFC, cD_pos_reg, cD_pos_reg_no_HFC};
t = cell(length(Dlist), length(run_cell), length(process_order));
sigt = cell(length(Dlist), length(run_cell), length(process_order));
peak_tissues = cell(length(Dlist), length(run_cell), length(process_order), 2);

% Run for all preprocessing steps and combined and not combined runs
Dcounter = 0;
for DD = Dlist
    Dcounter = Dcounter+1;
    DD = DD{1};

    for ss = 1:length(process_order)
        clear tinds
        for run_counter = 1:size(DD,1)
            if ~isempty(DD{run_counter,ss})
                Dtmp = DD{run_counter, ss};
    
                % If including bad trials, set all trials to good
                if incl_bad_trials
                    badtriallist = badtrials(Dtmp);
                    Dtmp = badtrials(Dtmp, 1:size(Dtmp,3), 0);
                    save(Dtmp);
                    Dbad_trial_copy = Dtmp;
                end
    
                % Just select subset of trials for analysis
                full_trial_list = indtrial(Dtmp, 'Cond1', 'GOOD');
                rng('default');
                rng(95);
                if size(Dtmp, 3) < 1000 % If runs have not been combined
                    trials = randperm(size(Dtmp,3), 540);
                else
                    trials = randperm(size(Dtmp,3), 1120);
                end
                trials_to_set_to_bad = setdiff(full_trial_list, trials);
                Dtmp = badtrials(Dtmp, trials_to_set_to_bad, 1);
                S = [];
                S.D = Dtmp;
                Dtmp = spm_eeg_remove_bad_trials(S);
                Dfname = Dtmp.fname;
    
                % Only run minimum norm if not already run
                if ~isfield(Dtmp.inv{1}, 'inverse')
        
                    matlabbatch = [];
                    matlabbatch{1}.spm.meeg.source.invert.D = {Dfname};
                    matlabbatch{1}.spm.meeg.source.invert.val = 1;
                    matlabbatch{1}.spm.meeg.source.invert.whatconditions.all = 1;
                    matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.invtype = 'IID';
                    matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.woi = [-Inf Inf];
                    matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.foi = [2 40];
                    matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
                    matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
                    matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
                    matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
                    matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
                    matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
                    matlabbatch{1}.spm.meeg.source.invert.modality = {'All'};
                    matlabbatch{2}.spm.meeg.source.results.D(1) = cfg_dep('Source inversion: M/EEG dataset(s) after imaging source reconstruction', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
                    matlabbatch{2}.spm.meeg.source.results.val = 1;
                    matlabbatch{2}.spm.meeg.source.results.woi = [0 200];
                    matlabbatch{2}.spm.meeg.source.results.foi = [2 40];
                    matlabbatch{2}.spm.meeg.source.results.ctype = 'evoked';
                    matlabbatch{2}.spm.meeg.source.results.space = 0;
                    matlabbatch{2}.spm.meeg.source.results.format = 'mesh';
                    matlabbatch{2}.spm.meeg.source.results.smoothing = 8;
            
                    a = spm_jobman('run',matlabbatch);
        
                    % Display
                    spm_mesh_render('Disp', a{2}.files{:});
                    spm_mesh_render('clim', gca, [0,8.5]);
            
                    H = spm_mesh_render('View',gca, [-90,10]);
                    spm_mesh_inflate(H.patch,Inf,1);
                    set(gcf, 'color', 'w');
                    if incl_bad_trials
                        export_fig(fullfile(resultsSavePath, sprintf('min_norm_incl_sat_trials_filename_%s_left', Dfname)), '-png', '-painters')
                    else
                        export_fig(fullfile(resultsSavePath, sprintf('min_norm_excl_sat_trials_filename_%s_left', Dfname)), '-png', '-painters')
                    end
                    spm_mesh_render('View',gca, [90,10]);
                    if incl_bad_trials
                        export_fig(fullfile(resultsSavePath, sprintf('min_norm_incl_sat_trials_filename_%s_right', Dfname)), '-png', '-painters')
                    else
                        export_fig(fullfile(resultsSavePath, sprintf('min_norm_excl_sat_trials_filename_%s_right', Dfname)), '-png', '-painters')
                    end
                    if incl_bad_trials
                        Dbad_trial_copy = badtrials(Dbad_trial_copy, badtriallist, 1);
                        save(Dbad_trial_copy);
                    end
                    Dtmp = spm_eeg_load(fullfile(Dtmp));
                end
    
                % Create source time courses from Minimum Norm results
            
                % Get weights
                U = Dtmp.inv{1}.inverse.U{1};
                weights = Dtmp.inv{1}.inverse.M;
        
                % Create each time course
                chaninds = selectchannels(Dtmp, Dtmp.inv{1}.forward.channels);
                if ss == 1
                    [~, tinds(1)] = min(abs(Dtmp.time - 92*1e-3));
                    [~, tinds(2)] = min(abs(Dtmp.time - 110*1e-3));
                else
                    [~, tinds] = min(abs(Dtmp.time - 100*1e-3));
                end
                dat = Dtmp(chaninds,tinds,:);
                nt = Dtmp.ntrials;
                vedata = cell(1, nt);
        
                for trial=1:nt
                    data = dat(:,:,trial);
                    vedata{trial} = transpose((U*data)'*weights');
                end
        
                % Rearrange vedata
                vedata_tp = cell(1, size(dat,2)); % One for each time point
                for tp = 1:length(tinds)
                    vedata_tp{tp} = zeros(size(weights,1), size(Dtmp,3));
                    for trial = 1:nt
                        vedata_tp{tp}(:,trial) = vedata{trial}(:,tp);
                    end
                end
        
                % T-stat over source time courses
                t{Dcounter,run_counter,ss} = ...
                    cell2mat(cellfun(@(x) mean(x,2)./(std(x,[],2)/sqrt(size(x,2))), vedata_tp, 'UniformOutput', false));
    
                % Smooth image
                t{Dcounter,run_counter,ss} = spm_mesh_smooth(M_original, t{Dcounter,run_counter,ss}, n);
            
                % Find significant values
                df = nt-1;
                alpha = 0.025/(size(Dtmp.inv{1}.inverse.M,1)*length(tinds));
                sigt{Dcounter,run_counter,ss} = tinv([alpha 1-alpha], df);
        
                % Display spatial maps
                cfg = [];
                cfg.facecolor = [0.4 0.4 0.4];
                cfg.vertexcolor = 'none';
    
                % Plot for each time point of interest
                for toi = 1:length(tinds)
                    figure;
                    mesh.pow = abs(t{Dcounter,run_counter,ss}(:,toi));
                    mesh.mask = mesh.pow > sigt{Dcounter,run_counter,ss}(2);
    
                    if max(mesh.pow) > sigt{Dcounter,run_counter,ss}
                        cfg.method         = 'surface';
                        if ss == 1
                            cfg.funcolorlim    = [4.4433 9.8450];
                        elseif ss == 2
                            cfg.funcolorlim    = [4.4433 15];
                        else
                            cfg.funcolorlim = [4.4433 20];
                        end
                        cfg.funparameter   = 'pow';
                        cfg.maskparameter  = 'mask';
                        cfg.funcolormap    = 'hot';
                        cfg.colorbartext = 't-stat';
                        ft_sourceplot(cfg, mesh);
                    else
                        surf.pos = mesh.pos;
                        surf.tri = mesh.tri;
                        ft_plot_mesh(surf,'edgecolor', 'none', 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
                        lighting gouraud
                        camlight
                    end
                    view ([-90 0])             % rotate the object in the view
                    camlight('headlight')
                    set(gcf, 'color', 'w');
                    set(gca, 'FontSize', 26);
                    material dull
            
                    if incl_bad_trials
                        export_fig(fullfile(resultsSavePath, sprintf('source_level_t_incl_sat_trials_%.fms_%s_left', 1e3*Dtmp.time(tinds(toi)), Dfname)), '-png', '-painters');
                    else
                        export_fig(fullfile(resultsSavePath, sprintf('source_level_t_excl_sat_trials_%.fms_%s_left', 1e3*Dtmp.time(tinds(toi)), Dfname)), '-png', '-painters');
                    end
                    view ([90 0])   
                    camlight('headlight')
                    if incl_bad_trials
                        export_fig(fullfile(resultsSavePath, sprintf('source_level_t_incl_sat_trials_%.fms_%s_right', 1e3*Dtmp.time(tinds(toi)), Dfname)), '-png', '-painters');
                    else
                        export_fig(fullfile(resultsSavePath, sprintf('source_level_t_excl_sat_trials_%.fms_%s_right', 1e3*Dtmp.time(tinds(toi)), Dfname)), '-png', '-painters');
                    end
    
                    % Compare localisation with auditory cortex in atlas
                
                    % Find peak t value
                    [~, maxind] = max(abs(t{Dcounter,run_counter,ss}(:,toi)));
                    if mesh_tissue.tissue(maxind) ~= 0
                        peak_tissues{Dcounter,run_counter,ss,toi} = mesh_tissue.tissuelabel{mesh_tissue.tissue(maxind)};
                    else
                        peak_tissues{Dcounter,run_counter,ss,toi} = maxind;
                    end
                end
            end
        end
    end
end

clear vedata vedata_tp weights U atlas dat data M_original mesh mesh_tissue

%% Get head position and orientation in room space

for ss = 1:length(process_order)
    if ss == 1
        marker_slots = [40, 44, 25, 12, 46];
    elseif ss == 2
        marker_slots = [47, 4, 39, 57, 8, 43];
    else
        marker_slots = [55, 25, 63, 19, 51, 22];
    end
    table_of_info = fullfile(rawDataPath, sprintf('sub-%s', process_order(ss).sub), sprintf('ses-%s', process_order(ss).ses), 'table_of_info.csv');
    clear R T
    
    % Get marker positions in MRI coordinates
    opti_data_mm{8,ss} = optitrack_convert_units(opti_data{8,ss}, 'mm');
    MarkerPosMRI = getMarkerPosInMRIcoords(table_of_info,...
        marker_slots, 55.27*ones(1,length(marker_slots)), opti_data_mm{8,ss});
    
    % Get transformation matrices to go from MRI to room coordinates
    for run_counter = 1:length(run_cell)
        if ~isempty(D{run_counter,ss})
            opti_data_mm{run_counter,ss} = optitrack_convert_units(opti_data{run_counter,ss}, 'mm');
            [~, ~, R{run_counter}, T{run_counter}] = getMagPosOriOverTime(opti_data_mm{run_counter,ss}, MarkerPosMRI, D{run_counter,ss});
        end
    end
    
    % Find head center over time
    ctx = gifti(D{8,ss}.inv{1}.mesh.tess_ctx);
    cortex_center_MRI = transpose(mean(ctx.vertices,1));
    clear ctx
    cortex_center_Room{ss} = cell(1, length(run_cell));
    for run_counter = 1:length(run_cell)
        if ~isempty(T{run_counter})
            cortex_center_Room{ss}{run_counter} = zeros(3, size(T{run_counter},2));
            for tt = 1:size(T{run_counter},2)
                cortex_center_Room{ss}{run_counter}(:,tt) = R{run_counter}(:,:,tt)*cortex_center_MRI + T{run_counter}(:,tt);
            end
        end
    end
end

clear T R opti_data_mm

%% Interpolate missing position data for cortex centre and head orientation

for ss = 1:length(process_order)
    good_opt_data{ss} = cell(1,length(run_cell));
    interp_cortex_center{ss} = cell(1,length(run_cell));
    
    for run_number = 1:length(run_cell)

        if ~isempty(cortex_center_Room{ss}{run_number})
    
            % Interpolate optitrack data as best as possible
            cortex_center = transpose(cortex_center_Room{ss}{run_number}*1e-3);
            
            good_opt_data{ss}{run_number} = opti_data{run_number,ss}.rigidbodies.data(:,8) > 0;
            interp_cortex_center{ss}{run_number} = cortex_center;
    
            % First, assume any missing data at the beginning of the recording by
            % repeating first non-zero datapoint
            first_good_val = find(good_opt_data{ss}{run_number}, 1);
            interp_cortex_center{ss}{run_number}(1:first_good_val-1,:) = repmat(cortex_center(first_good_val,:), first_good_val-1, 1);
        
            % Then interpolate later missing data
            % - choose mising data points to interpolate over
            [~, step_in] = findpeaks(diff(~good_opt_data{ss}{run_number}));
            [~, step_out] = findpeaks(-diff(~good_opt_data{ss}{run_number}));
            step_in = step_in + 1;
            if first_good_val ~= 1
                step_out = step_out(2:end);
            end
        
            % Linearly interpolate gaps which are smaller than 200 samples
            gapsize = step_out - step_in;
            for gap = 1:length(step_in)
                if gapsize(gap) < 200
                    interp_cortex_center{ss}{run_number}(step_in(gap):step_out(gap),:) = interp1([step_in(gap)-1, step_out(gap)+1], ...
                        cortex_center([step_in(gap)-1, step_out(gap)+1],:), step_in(gap):step_out(gap), 'linear');
        
                    % Mark that this data is okay
                    good_opt_data{ss}{run_number}(step_in(gap):step_out(gap)) = 1;
        
                    % reset the missing data boundaries (to effectively mark this
                    % data as interpolated)
                    step_in(gap) = NaN;
                    step_out(gap) = NaN;
                end
            end
            step_in = step_in(~isnan(step_in));
            step_out = step_out(~isnan(step_out));
        
            % Estimate which data is within 0.5 m radius cylinder in room center
            % If optitrack data is marked as missing/bad, assuming the participant
            % was outside of the central 0.5 m radius cylinder
            data_within_cylinder{ss}{run_number} = find(sqrt(sum(interp_cortex_center{ss}{run_number}(:,[1,3]).^2,2))<0.5 & good_opt_data{ss}{run_number});
            
            % Shift step_out back by 8 points (consistently too small, I think due
            % to the previous interpolation of the optitrack data from sampling at
            % 120 Hz to 1000 Hz, which is roughly 8 samples)
            for gap = 1:length(step_out)
                good_opt_data{ss}{run_number}(step_out(gap):step_out(gap)+8) = 0;
            end
            step_out = step_out + 8;

            good_opt_data{ss}{run_number} = int8(good_opt_data{ss}{run_number});
            if ss == 1
                % Add in some datapoints estimated from video to make sure front of the room is included
                % Mark by setting good_opt_data to 2 at these timepoints
                if run_number == 8 % block 3 in recording so 2nd walking feedback on block
                    interp_cortex_center{ss}{run_number}(57000,:) = [0, interp_cortex_center{ss}{run_number}(57000,2),2-0.8];
                    interp_cortex_center{ss}{run_number}(1e3*(2*60+15),:) = [-0.1, interp_cortex_center{ss}{run_number}(1e3*(2*60+15),2),2-0.8];
                    interp_cortex_center{ss}{run_number}(1e3*(3*60+20),:) = [0, interp_cortex_center{ss}{run_number}(1e3*(3*60+20),2),2-0.8];
                    interp_cortex_center{ss}{run_number}(1e3*(4*60+27),:) = [0, interp_cortex_center{ss}{run_number}(1e3*(4*60+27),2),2-0.7];
                    good_opt_data{ss}{run_number}(57000) = 2; 
                    good_opt_data{ss}{run_number}(1e3*(2*60+15)) = 2; 
                    good_opt_data{ss}{run_number}(1e3*(3*60+20)) = 2; 
                    good_opt_data{ss}{run_number}(1e3*(4*60+27)) = 2; 
                elseif run_number == 7 % block 2 in recording so 1st walking feedback on block
                    interp_cortex_center{ss}{run_number}(49000,:) = [0, interp_cortex_center{ss}{run_number}(49000,2),2-0.7];
                    interp_cortex_center{ss}{run_number}(1e3*(1*60+15),:) = [-1, interp_cortex_center{ss}{run_number}(1e3*(1*60+15),2),0.2];
                    interp_cortex_center{ss}{run_number}(1e3*(2*60+3),:) = [0, interp_cortex_center{ss}{run_number}(1e3*(2*60+3),2),2-0.7];
                    interp_cortex_center{ss}{run_number}(1e3*(3*60+9),:) = [-0.1, interp_cortex_center{ss}{run_number}(1e3*(3*60+9),2),2-0.7];
                    interp_cortex_center{ss}{run_number}(1e3*(4*60+20),:) = [-0.2, interp_cortex_center{ss}{run_number}(1e3*(4*60+20),2),2-0.7];
                    good_opt_data{ss}{run_number}(49000) = 2; 
                    good_opt_data{ss}{run_number}(1e3*(1*60+15)) = 2; 
                    good_opt_data{ss}{run_number}(1e3*(2*60+3)) = 2; 
                    good_opt_data{ss}{run_number}(1e3*(3*60+9)) = 2; 
                    good_opt_data{ss}{run_number}(1e3*(4*60+20)) = 2; 
                elseif run_number == 6 % block 4 in recording but 2nd walking feedback off block
                    interp_cortex_center{ss}{run_number}(54000,:) = [-0.1, interp_cortex_center{ss}{run_number}(54000,2),2-0.7];
                    interp_cortex_center{ss}{run_number}(1e3*(2*60+3),:) = [-0.1, interp_cortex_center{ss}{run_number}(1e3*(2*60+3),2),2-0.7];
                    interp_cortex_center{ss}{run_number}(1e3*(3*60+14),:) = [-0.1, interp_cortex_center{ss}{run_number}(1e3*(3*60+14),2),2-0.6];
                    interp_cortex_center{ss}{run_number}(1e3*(4*60+8),:) = [-0.1, interp_cortex_center{ss}{run_number}(1e3*(4*60+8),2),2-0.75];
                    interp_cortex_center{ss}{run_number}(1e3*(5*60+8),:) = [-0.1, interp_cortex_center{ss}{run_number}(1e3*(4*60+8),2),2-0.7];
                    good_opt_data{ss}{run_number}(54000) = 2;  
                    good_opt_data{ss}{run_number}(1e3*(2*60+3)) = 2; 
                    good_opt_data{ss}{run_number}(1e3*(3*60+14)) = 2; 
                    good_opt_data{ss}{run_number}(1e3*(4*60+8)) = 2;  
                    good_opt_data{ss}{run_number}(1e3*(5*60+8)) = 2; 
                end
            end
        
            % Fit gap
            for kk = 1:length(step_in)
                k = step_in(kk):step_out(kk);
        
                if step_in(kk)-6e3 > 0 && step_out(kk)+6e3 < size(interp_cortex_center{ss}{run_number},1)
                    x = step_in(kk)-6e3:step_out(kk)+6e3;
                    y =interp_cortex_center{ss}{run_number}(step_in(kk)-6e3:step_out(kk)+6e3,:);
                elseif step_in(kk)-6e3 < 0
                    x = 1:step_out(kk)+6e3;
                    y = interp_cortex_center{ss}{run_number}(1:step_out(kk)+6e3,:);
                elseif step_out(kk)+6e3 > size(interp_cortex_center{ss}{run_number},1)
                    x = step_in(kk)-6e3:size(interp_data{id},1);
                    y = interp_cortex_center{ss}{run_number}(step_in(kk)-6e3:end,:);
                end
                
                % Just select good data
                good_data_in_section = good_opt_data{ss}{run_number}(x) >= 1;
                x = x(good_data_in_section);
                y = y(good_data_in_section,:);
            
                for cc = 1:3
                    fobj = fit(x', y(:,cc), 'pchip');
                    interp_cortex_center{ss}{run_number}(k, cc) = fobj(k);
                end
            end
        end
    end
end

clear cortex_center_Room cortex_center_MRI cortex_center kk k step_in step_out good_data_in_section

%% Estimate distance between furthest points walked

for ss = 1:length(process_order)
    for run_number = 1:length(run_cell)
        if ~isempty(interp_cortex_center{ss}{run_number})
            % Find most extreme points
            mpos = mean(interp_cortex_center{ss}{run_number},1);
            dist_from_mpos = sqrt(sum((interp_cortex_center{ss}{run_number} - mpos).^2, 2));
            clear mpos
            % Take 10% of data furthest from mean position
            [~, ind] = sort(dist_from_mpos);
            clear dist_from_mpos
            opti_data_subset = interp_cortex_center{ss}{run_number}(ind(floor(0.9*length(ind)):end),:);
            dists = pdist(opti_data_subset);
            maxdist = max(dists);
            
            fprintf('Participant %s, Session %s, Run %s: max translation: %.2f cm\n', process_order(ss).sub, process_order(ss).ses, run_cell{run_number}, maxdist*100);
        end
    end
end

clear ind dists maxdist opti_data_subset

%% Trajectory of participant through MSR (figure 5)

clear H

for ss = 1:length(process_order)

    % Set up figures
    
    % modelled field
    figure; hold on; grid on; box on;
    H(1) = subplot(1,1,1);
    
    % moved trajectory/where optitrack is poorly sampled
    figure; hold on; grid on; box on;
    H(2) = subplot(1,1,1);
    
    % sat vs unsat, feedback off
    figure; hold on; grid on; box on;
    H(3) = subplot(1,1,1);
    
    % sat vs unsat, feedback on
    figure; hold on; grid on; box on;
    H(4) = subplot(1,1,1);

    % sat vs unsat, Y channels, feedback off
    figure; hold on; grid on; box on;
    H(5) = subplot(1,1,1);

    % sat vs unsat, Z channels, feedback off
    figure; hold on; grid on; box on;
    H(6) = subplot(1,1,1);

    % sat vs unsat, Y channels, feedback on
    figure; hold on; grid on; box on;
    H(7) = subplot(1,1,1);

    % sat vs unsat, Z channels, feedback on
    figure; hold on; grid on; box on;
    H(8) = subplot(1,1,1);

    for run_number = 5:8

        % Additional figure to save position just for this run
        figure; hold on; grid on; box on;
        G = subplot(1,1,1);
        
        % plot (modelled) field at each position
        
        % Estimate field at each position
        if run_number > 6
            Dfname = D{run_number, ss}.fname;
            Dfname = Dfname(strfind(Dfname,'t_dsub-'):end);
            DD = spm_eeg_load(Dfname);
            clear Dfname
            X = DD.sensors('MEG').chanori;
            Y = feedback_applied{run_number, ss}(selectchannels(DD, DD.sensors('MEG').label),:);
            B = Y'*pinv(X');
            Bmag = sqrt(sum(B.^2, 2));
            scatter3(H(1), interp_cortex_center{ss}{run_number}(:,1), ...
                interp_cortex_center{ss}{run_number}(:,2), ...
                interp_cortex_center{ss}{run_number}(:,3), 5, 1e-6*Bmag, ...
                'filled');
        end
    
        % plot where positions were interpolated
        p = plot3(H(2), interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==1,1), ...
            interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==1,2), ...
            interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==1,3), '.', 'color', [0.2666, 0.4471, 0.7686]);
        if any(good_opt_data{ss}{run_number} == 0)
            q = plot3(H(2), interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==0,1), ...
                interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==0,2), ...
                interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==0,3), ...
                '.', 'color', [0.5647, 0.6706, 0.8627], 'MarkerSize', 0.5);
        end
        if any(good_opt_data{ss}{run_number} == 2)
            u = plot3(H(2), interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==2,1), ...
                interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==2,2), ...
                interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==2,3), 'rx');
        end

        % Repeat just for this run
        plot3(G, interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==1,1), ...
            interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==1,2), ...
            interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==1,3), '.', 'color', [0.2666, 0.4471, 0.7686]);
        if any(good_opt_data{ss}{run_number} == 0)
            plot3(G, interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==0,1), ...
                interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==0,2), ...
                interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==0,3), ...
                '.', 'color', [0.5647, 0.6706, 0.8627], 'MarkerSize', 0.5);
        end
        if any(good_opt_data{ss}{run_number} == 2)
            plot3(G, interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==2,1), ...
                interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==2,2), ...
                interp_cortex_center{ss}{run_number}(good_opt_data{ss}{run_number}==2,3), 'rx');
        end
    
        % plot where sensors saturated
        all_sat_frames = cell2mat(saturated_frames{run_number, ss}');
        all_sat_frames = unique(all_sat_frames);
        unSatFrames = setdiff(1:size(interp_cortex_center{ss}{run_number},1), all_sat_frames);
        
        if run_number < 7
            plot_no = 3;
        else
            plot_no = 4;
        end
        
        r = plot3(H(plot_no), interp_cortex_center{ss}{run_number}(unSatFrames,1), ...
            interp_cortex_center{ss}{run_number}(unSatFrames,2), ...
            interp_cortex_center{ss}{run_number}(unSatFrames,3), '.', 'color', [194, 192, 148]./255);
        s = plot3(H(plot_no), interp_cortex_center{ss}{run_number}(all_sat_frames,1), ...
            interp_cortex_center{ss}{run_number}(all_sat_frames,2), ...
            interp_cortex_center{ss}{run_number}(all_sat_frames,3), '.', 'color', [172, 56, 56]./255);
        legend(H(plot_no), [r, s], {'Unsaturated', 'Saturated'}, 'location', 'south')

        % Plot where sensors saturated, Y channels
        chans = D{run_number, ss}.chanlabels(indchantype(D{run_number,ss}, 'MEGMAG', 'GOOD'));
        Ychans = endsWith(chans, '-Y');
        Y_sat_frames = cell2mat(saturated_frames{run_number, ss}(Ychans)');
        Y_sat_frames = unique(Y_sat_frames);
        unSatFrames = setdiff(1:size(interp_cortex_center{ss}{run_number},1), Y_sat_frames);

        if run_number < 7
            plot_no = 5;
        else
            plot_no = 7;
        end

        r = plot3(H(plot_no), interp_cortex_center{ss}{run_number}(unSatFrames,1), ...
            interp_cortex_center{ss}{run_number}(unSatFrames,2), ...
            interp_cortex_center{ss}{run_number}(unSatFrames,3), '.', 'color', [194, 192, 148]./255);
        s = plot3(H(plot_no), interp_cortex_center{ss}{run_number}(Y_sat_frames,1), ...
            interp_cortex_center{ss}{run_number}(Y_sat_frames,2), ...
            interp_cortex_center{ss}{run_number}(Y_sat_frames,3), '.', 'color', [172, 56, 56]./255);
        legend(H(plot_no), [r, s], {'Unsaturated', 'Saturated'}, 'location', 'south')

        % Plot where sensors saturated, Z channels
        Zchans = endsWith(chans, '-Z');
        Z_sat_frames = cell2mat(saturated_frames{run_number, ss}(Zchans)');
        Z_sat_frames = unique(Z_sat_frames);
        unSatFrames = setdiff(1:size(interp_cortex_center{ss}{run_number},1), Z_sat_frames);

        if run_number < 7
            plot_no = 6;
        else
            plot_no = 8;
        end

        r = plot3(H(plot_no), interp_cortex_center{ss}{run_number}(unSatFrames,1), ...
            interp_cortex_center{ss}{run_number}(unSatFrames,2), ...
            interp_cortex_center{ss}{run_number}(unSatFrames,3), '.', 'color', [194, 192, 148]./255);
        s = plot3(H(plot_no), interp_cortex_center{ss}{run_number}(Z_sat_frames,1), ...
            interp_cortex_center{ss}{run_number}(Z_sat_frames,2), ...
            interp_cortex_center{ss}{run_number}(Z_sat_frames,3), '.', 'color', [172, 56, 56]./255);
        legend(H(plot_no), [r, s], {'Unsaturated', 'Saturated'}, 'location', 'south')

        good_opt_data_in_beeps = good_opt_data{ss}{run_number};
        good_opt_data_in_beeps(1:beep_start_stop{run_number,ss}(1)-1) = 0;
        good_opt_data_in_beeps(beep_start_stop{run_number,ss}(2)+1:end) = 0;
    
        fprintf('Run %s:\n', run_cell{run_number})
        fprintf('Range of movement (cortex centre) = (%.2f, %.2f, %.2f) m\n', ...
            range(interp_cortex_center{ss}{run_number}(:,1)), range(interp_cortex_center{ss}{run_number}(:,2)), ...
            range(interp_cortex_center{ss}{run_number}(:,3)));
        fprintf('Range of movement (rigid body centre) = (%.2f, %.2f, %.2f) m\n', ...
            range(opti_data{run_number,ss}.rigidbodies.data(good_opt_data_in_beeps==1,5)), range(opti_data{run_number,ss}.rigidbodies.data(good_opt_data_in_beeps==1,6)), ...
            range(opti_data{run_number,ss}.rigidbodies.data(good_opt_data_in_beeps==1,7)));
        eulZYX = quat2eul(opti_data{run_number,ss}.rigidbodies.data(good_opt_data_in_beeps==1, [4, 1:3]));
        eulZYX = eulZYX*180/pi;
        fprintf('Yaw, pitch, roll = (%.2f, %.2f, %.2f) degrees\n', range(eulZYX(:,2)), range(eulZYX(:,3)), range(eulZYX(:,1)));        
        rbspeed = sqrt(sum(diff(opti_data{run_number,ss}.rigidbodies.data(:,5:7)).^2, 2))./diff(opti_data{run_number,ss}.time)';
        fprintf('Average speed (rigid body centre) = %.2f m/s\n', mean(rbspeed));
        rbspeed = sqrt(sum(diff(interp_cortex_center{ss}{run_number}).^2, 2))./diff(opti_data{run_number,ss}.time)';
        fprintf('Average speed (cortex centre) = %.2f m/s\n', mean(rbspeed));

        % Save trajectory for this run
        zlim(G, [-2 2]);
        ylim(G, [0, 1.1]);
        xlim(G, [-1.5, 1.5]);
        view(G, 180,0);
        daspect(G, [1 1 1]);
        set(G.Parent, 'Position', [680   670   497   308]);
        set(G, 'FontSize', 14);
        xlabel(G, 'X, left-right (m)')
        ylabel(G, 'Y, up-down (m)')
        zlabel(G, 'Z, forward-back (m)');
        set(G.Parent, 'color', 'w');  
        figure(G.Parent);
        export_fig(fullfile(resultsSavePath, sprintf('trajectory_participant_%s_session_%s_run_%s', process_order(ss).sub, process_order(ss).ses, run_cell{run_number})), '-png', '-painters')   

        % Repeat in YZ plane
        view(G, -90, 0);
        set(G, 'CameraUpVector', [0 1 0]);
        set(G, 'Position', [0.13,0.11,0.775,0.815]);
        set(G, 'CameraViewAngle', 6.311869756640776);
        export_fig(fullfile(resultsSavePath, sprintf('trajectory_participant_%s_session_%s_run_%s_YZ', process_order(ss).sub, process_order(ss).ses, run_cell{run_number})), '-png', '-painters') 
    end

    % Adjust colour scale for estimated field strength
    colormap(H(1), brewermap(64,'OrRd'));
    caxis(H(1), [0 4]);
    cb = colorbar(H(1));
    set(cb, 'FontSize', 14);
    title(cb, 'B (nT)', 'FontSize', 14);
    set(cb, 'Position', [0.7201    0.1871    0.0384    0.7188]);

    % Add a 50 cm circle (for inside/outside)
    for ii = 3:8
        theta = linspace(0, 2*pi, 100);
        z = sin(theta)*0.5;
        x = cos(theta)*0.5;
        plot3(H(ii), x, ones(1,length(theta)), z, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end

    % Add legends
    if ~isempty(q)
        if exist('u', 'var')
            legend(H(2), [p, q, u], {'Recorded', 'Interpolated', 'Manually Added'}, 'location', 'south')
        else
            legend(H(2), [p, q], {'Recorded', 'Interpolated'}, 'location', 'south')
        end
    end
    
    % Other figure settings
    for ii = 1:length(H)
        zlim(H(ii), [-2 2]);
        xlim(H(ii), [-1.5, 1.5]);
        ylim(H(ii), [0, 1.1]);
        view(H(ii), 180,0);
        daspect(H(ii), [1 1 1]);
        set(H(ii).Parent, 'Position', [680   670   497   308]);
        set(H(ii), 'FontSize', 14);
        xlabel(H(ii), 'X, left-right (m)')
        ylabel(H(ii), 'Y, up-down (m)')
        zlabel(H(ii), 'Z, forward-back (m)');
        set(H(ii).Parent, 'color', 'w');
    end
    
    % Save
    figure(H(2).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('combined_trajectories_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')   
    figure(H(1).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('estimated_field_strength_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')   
    figure(H(3).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_off_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')    
    figure(H(4).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_on_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')   

    figure(H(5).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_off_Ychans_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')    
    figure(H(6).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_off_Zchans_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')    
    figure(H(7).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_on_Ychans_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')    
    figure(H(8).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_on_Zchans_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters') 

    % Change view to side on and re-save
    for ii = 1:length(H)
        view(H(ii), -90, 0);
        set(H(ii), 'CameraUpVector', [0 1 0]);
        set(H(ii), 'Position', [0.13,0.11,0.775,0.815]);
        set(H(ii), 'CameraViewAngle', 6.311869756640776);
    end
    set(cb, 'Position', [0.897162374245484,0.266233766233766,0.034427162977856,0.568237662337668]);
    figure(H(2).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('combined_trajectories_YZ_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')   
    figure(H(1).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('estimated_field_strength_YZ_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')   
    figure(H(3).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_off_YZ_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')    
    figure(H(4).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_on_YZ_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters') 

    figure(H(5).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_off_Ychans_YZ_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')    
    figure(H(6).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_off_Zchans_YZ_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')    
    figure(H(7).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_on_Ychans_YZ_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters')    
    figure(H(8).Parent);
    export_fig(fullfile(resultsSavePath, sprintf('saturated_frames_feedback_on_Zchans_YZ_participant_%s_session_%s', process_order(ss).sub, process_order(ss).ses)), '-png', '-painters') 

    clear u
end

clear B Bmag chans Ychans Y_sat_frames Z_sat_frames

%% Number of trials within central 1m vs without

for ss = 1:length(process_order)
    within_cylinder_trial_list = cell(1, length(run_cell));
    for run_counter = 1:length(run_cell)

        if ~isempty(D{run_counter,ss})
            % Get trl matrix
            Dfname = D{run_counter, ss}.fname;
            D{run_counter,ss} = spm_eeg_load(fullfile(D{run_counter,ss}.path, Dfname(3:end)));
            S = [];
            S.D = D{run_counter,ss};
            S.timewin = [-200 500];
            S.condlabels = {'tone'};
            S.triggerChannels = {process_order(ss).stimTrig};
            [D{run_counter,ss}, trl{run_counter}] = spm_opm_epoch_trigger(S);
    
            % Check if trial is within cylinder
            for tt = 1:size(trl{run_counter},1)
                if length(intersect(trl{run_counter}(tt,1):trl{run_counter}(tt,2), data_within_cylinder{ss}{run_counter})) ==...
                        length(trl{run_counter}(tt,1):trl{run_counter}(tt,2))
                    within_cylinder_trial_list{run_counter} = cat(1, within_cylinder_trial_list{run_counter}, tt);
                end
            end
        end
    end
    within_cylinder_trial_list_merged{1} = [within_cylinder_trial_list{1}; ...
        size(trl{1},1)+within_cylinder_trial_list{2}];
    within_cylinder_trial_list_merged{2} = [within_cylinder_trial_list{3}; ...
        size(trl{3},1)+within_cylinder_trial_list{4}];
    within_cylinder_trial_list_merged{3} = [within_cylinder_trial_list{5}; ...
        size(trl{5},1)+within_cylinder_trial_list{6}];
    within_cylinder_trial_list_merged{4} = [within_cylinder_trial_list{7}; ...
        size(trl{7},1)+within_cylinder_trial_list{8}];

    for run_counter = 1:4
        if ~isempty(cD{run_counter,ss})
            % Randomly select 1120 trials
            rng('default');
            rng(95);
            trials = randperm(size(cD{run_counter, ss},3), 1120);
        
            within_50_tot = intersect(trials, within_cylinder_trial_list_merged{run_counter});
            within_50_sat = intersect(within_50_tot, indtrial(cD{run_counter, ss}, 'Cond1', 'BAD'));
            fprintf('Run %.f, Within 50 cm: saturated %.f, unsaturated %.f, total %.f\n', run_counter, ...
                length(within_50_sat), length(within_50_tot)-length(within_50_sat), length(within_50_tot));
        
            out_50_tot = setdiff(trials, within_cylinder_trial_list_merged{run_counter});
            out_50_sat = intersect(out_50_tot, indtrial(cD{run_counter, ss}, 'Cond1', 'BAD'));
            fprintf('Run %.f, Outside 50 cm: saturated %.f, unsaturated %.f, total %.f\n', run_counter, ...
                length(out_50_sat), length(out_50_tot)-length(out_50_sat), length(out_50_tot));
        end
    end
end
