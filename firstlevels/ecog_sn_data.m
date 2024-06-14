classdef ecog_sn_data < ecog_data
% ECOG_SN_DATA 
% TODO - description
    
properties
    %% ---- LANG ELECS ----
    s_vs_n_sig              
    s_vs_n_p_ratio
    s_vs_n_ops

    %% ---- DATA ----
    langloc_save_path
    langloc_crunched_file_name
    langloc_crunched_file_path
    preproc_class_file_name
    preproc_class_file_path
    
end
    
methods
    %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sn_obj = ecog_sn_data(...
            langloc_save_path,...
            langloc_crunched_file_name,...
            langloc_crunched_file_path,...
            preproc_class_file_name,...
            preproc_class_file_path)
    
        % load in preprocessed signal
        load(langloc_crunched_file_name) % var name : obj

        sn_obj@ecog_data(obj.for_preproc,...
                         obj.subject,...
                         obj.experiment,...
                         obj.crunched_file_name,... % obj.save_name,...
                         obj.crunched_file_path,... % obj.save_path,...
                         obj.raw_file_name,... % obj.file_name
                         obj.raw_file_path,... % obj.file_path
                         obj.elec_ch_label,...
                         obj.elec_ch,...
                         obj.elec_ch_prelim_deselect,...
                         obj.elec_ch_type...
        );

        % define ecog_sn_data specific properties
        sn_obj.langloc_save_path          = langloc_save_path;
        sn_obj.langloc_crunched_file_name = langloc_crunched_file_name;
        sn_obj.langloc_crunched_file_path = langloc_crunched_file_path;
        sn_obj.preproc_class_file_name    = preproc_class_file_name;
        sn_obj.preproc_class_file_path    = preproc_class_file_path;

        % add ecog_data properties not in constructor
        sn_obj.elec_data     = obj.elec_data;
        sn_obj.bip_elec_data = obj.bip_elec_data;
        sn_obj.stitch_index  = obj.stitch_index;
        sn_obj.sample_freq   = obj.sample_freq;
        sn_obj.trial_data    = obj.trial_data;

        sn_obj.trial_timing = obj.trial_timing;
        sn_obj.events_table = obj.events_table;
        sn_obj.condition    = obj.condition;
        sn_obj.session      = obj.session;
        
        sn_obj.elec_ch_with_IED      = obj.elec_ch_with_IED;
        sn_obj.elec_ch_with_noise    = obj.elec_ch_with_noise;
        sn_obj.elec_ch_user_deselect = obj.elec_ch_user_deselect;
        sn_obj.elec_ch_clean         = obj.elec_ch_clean;
        sn_obj.elec_ch_valid         = obj.elec_ch_valid;
        sn_obj.bip_ch                = obj.bip_ch;
        sn_obj.bip_ch_label          = obj.bip_ch_label;
        sn_obj.bip_ch_valid          = obj.bip_ch_valid;
        sn_obj.bip_ch_grp            = obj.bip_ch_grp;
        sn_obj.bip_ch_label_grp      = obj.bip_ch_label_grp;

    end 


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIND LANGUAGE RESPONSIVE ELECTRODES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function test_s_vs_n(obj,varargin)
        % TODO - description

        p = inputParser();
        addParameter(p,'words',1:12);
        addParameter(p,'S_condition_flag','S')
        addParameter(p,'N_condition_flag','N')
        addParameter(p,'n_rep',1000);
        addParameter(p,'corr_type','Spearman');
        addParameter(p,'threshold',0.01);
        addParameter(p,'side','both');
        addParameter(p,'sessions',[]);
        addParameter(p,'do_plot',false);
        parse(p, varargin{:});
        ops = p.Results;
        
        obj.s_vs_n_ops = ops;
            
        fprintf(1,'\n> Finding electrodes that respond significantly more to %s than %s ...\n',ops.S_condition_flag,ops.N_condition_flag);

        % only keep certain sessions for MITLanloc (e.g., remove speed tests)
        if ~isempty(ops.sessions)
            keep_trials = any(obj.session==ops.sessions,2);
        else
            keep_trials = [];
        end

        % S condition
        [S_ave_tbl,S_table] = obj.get_ave_cond_trial('words',ops.words,'condition',ops.S_condition_flag,'keep_trials',keep_trials);
        B = S_ave_tbl;
        S_ave = table2cell(B(:,~ismember(B.Properties.VariableNames,{'key','string'})));

        % N condition
        [N_ave_tbl,N_table] = obj.get_ave_cond_trial('words',ops.words,'condition',ops.N_condition_flag,'keep_trials',keep_trials);
        B = N_ave_tbl;
        keys = B(:,ismember(B.Properties.VariableNames,'key'));
        N_ave = table2cell(B(:,~ismember(B.Properties.VariableNames,{'key','string'})));

        % --- UNCOMMENT IF DESIRED ---
        % take ODD runs ONLY! (even will be used for plotting)
        % S_ave = cellfun(@(x) x(:,1:2:size(x,2)),S_ave,'uni',false);
        % N_ave = cellfun(@(x) x(:,1:2:size(x,2)),N_ave,'uni',false);

        % --- CORRELATION ---
        S_N_comb = cellfun(@(x,y) horzcat(x,y),S_ave,N_ave,'uni',false);
        S_N_flag = cellfun(@(x,y) horzcat(x*0+1,y*0-1),S_ave,N_ave,'uni',false);
        S_N_corr = cellfun(@(x,y) diag(corr(x',y','Type',ops.corr_type)),S_N_comb,S_N_flag,'uni',false);
        names_no_string = ~ismember(B.Properties.VariableNames,{'string'});
        S_N_corr_table = cell2table(horzcat('s_vs_n_corr',S_N_corr),'VariableNames',B.Properties.VariableNames(names_no_string));

        % --- PERMUTATION ---
        fprintf(1,'\n> Running permutation analysis ...\n');
        fprintf(1,'[');

        n_rep = ops.n_rep;
        S_N_corr_rnd_all = cell(size(S_N_corr));

        for k=1:n_rep
            if ~mod(k,100) 
                fprintf(1,'.');
            end
                
            S_N_flag_rnd_idx = cellfun(@(x) randperm(size(x,2)),S_N_flag,'uni',false);
            S_N_flag_rnd = cellfun(@(x,y) x(:,y),S_N_flag,S_N_flag_rnd_idx,'uni',false);
            S_N_corr_rnd = cellfun(@(x,y) diag(corr(x',y','Type','Spearman')),S_N_comb,S_N_flag_rnd,'uni',false);
            S_N_corr_rnd_all = cellfun(@(x,y) horzcat(x,y),S_N_corr_rnd_all,S_N_corr_rnd,'UniformOutput',false);
            
        end

        fprintf(1,'] done\n');
            
        S_N_corr_rnd_table = cell2table(horzcat('s_vs_n_corr_rnd',S_N_corr_rnd_all),'VariableNames',B.Properties.VariableNames(names_no_string));

            
        % --- SIGNFICANCE ---
        switch ops.side
            case 'left'
                S_N_p_ratio = cellfun(@(x,y) sum(x<y,2)/size(y,2),S_N_corr,S_N_corr_rnd_all,'uni',false);
                S_N_p_is_sig = cellfun(@(x) (x<(ops.threshold)),S_N_p_ratio,'uni',false);
            case 'right'
                S_N_p_ratio = cellfun(@(x,y) sum(x>y,2)/size(y,2),S_N_corr,S_N_corr_rnd_all,'uni',false);
                S_N_p_is_sig = cellfun(@(x) (x>(1-ops.threshold)),S_N_p_ratio,'uni',false);
            case 'both'
                S_N_p_ratio = cellfun(@(x,y) sum(x>y,2)/size(y,2),S_N_corr,S_N_corr_rnd_all,'uni',false);
                S_N_p_is_sig = cellfun(@(x)  ((x<(ops.threshold)) | (x>(1-ops.threshold) )),S_N_p_ratio,'uni',false);
        end

        S_N_p_ratio_tbl = cell2table(horzcat('s_vs_n_p_ratio',S_N_p_ratio),'VariableNames',B.Properties.VariableNames(names_no_string));
            
        S_N_p_is_sig{1} = obj.elec_ch_valid & S_N_p_is_sig{1}; % make sure only clean channels marked sig
        s_vs_n_sig = cell2table(horzcat('s_vs_n_sig',S_N_p_is_sig),'VariableNames',B.Properties.VariableNames(names_no_string));
        obj.s_vs_n_sig = s_vs_n_sig;
        obj.s_vs_n_p_ratio = S_N_p_ratio_tbl;


        % --- PLOTTING ---
        if ops.do_plot 
            obj.plot_s_vs_n(ops,...
                            'elec_data',...
                            obj.elec_ch_label,...
                            'S_table',S_table,...
                            'N_table',N_table,...
                            'S_N_corr_table',S_N_corr_table,...
                            'S_N_corr_rnd_table',S_N_corr_rnd_table,...
                            's_vs_n_sig',s_vs_n_sig,...
                            'S_N_p_ratio_tbl',S_N_p_ratio_tbl...
            );

            if ~isempty(obj.bip_elec_data)
                obj.plot_s_vs_n(ops,...
                                'bip_elec_data',...
                                obj.bip_ch_label,...
                                'S_table',S_table,...
                                'N_table',N_table,...
                                'S_N_corr_table',S_N_corr_table,...
                                'S_N_corr_rnd_table',S_N_corr_rnd_table,...
                                's_vs_n_sig',s_vs_n_sig,...
                                'S_N_p_ratio_tbl',S_N_p_ratio_tbl...
                );
            end
        end

    end

        
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AVERAGE S VS. N FOR LANGUAGE ELECTRODES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function lang_resp_plots(obj,varargin)
        % Plots average of significant electrodes. 

        p = inputParser();
        addParameter(p,'words',1:12);
        addParameter(p,'S_condition_flag','S');
        addParameter(p,'N_condition_flag','N');
        addParameter(p,'W_condition_flag',[]);
        addParameter(p,'J_condition_flag',[]);
        addParameter(p,'sessions',[]);
        addParameter(p,'subAverage',false);
        addParameter(p,'bipolarByShank',false);
        parse(p, varargin{:});
        ops = p.Results;

        if ops.subAverage
            obj.zscore_signal();
            obj.make_trials();
        end

        % -- TIMECOURSE ---
        % unipolar
        [S_values_uni,S_values_uni_sem] = obj.get_timecourses('words',ops.words,...
                                                              'condition',ops.S_condition_flag,...
                                                              'sessions',ops.sessions,...
                                                              'allElecs',ops.bipolarByShank);
        [N_values_uni,N_values_uni_sem] = obj.get_timecourses('words',ops.words,...
                                                              'condition',ops.N_condition_flag,...
                                                              'sessions',ops.sessions);
        if ops.W_condition_flag % SWJN
            [W_values_uni,W_values_uni_sem] = obj.get_timecourses('words',ops.words,...
                                                                  'condition',ops.W_condition_flag,...
                                                                  'sessions',ops.sessions);
            [J_values_uni,J_values_uni_sem] = obj.get_timecourses('words',ops.words,...
                                                                  'condition',ops.J_condition_flag,...
                                                                  'sessions',ops.sessions);
        end

        % bipolar
        if ~isempty(obj.bip_elec_data)
            [S_values_bip,S_values_bip_sem] = obj.get_timecourses('words',ops.words,...
                                                                  'condition',ops.S_condition_flag,...
                                                                  'sessions',ops.sessions,...
                                                                  'signalType','bipolar',...
                                                                  'allElecs',ops.bipolarByShank);
            [N_values_bip,N_values_bip_sem] = obj.get_timecourses('words',ops.words,...
                                                                  'condition',ops.N_condition_flag,...
                                                                  'sessions',ops.sessions,...
                                                                  'signalType','bipolar',...
                                                                  'allElecs',ops.bipolarByShank);
            if ops.W_condition_flag % SWJN
                [W_values_bip,W_values_bip_sem] = obj.get_timecourses('words',ops.words,...
                                                                      'condition',ops.W_condition_flag,...
                                                                      'sessions',ops.sessions,...
                                                                      'signalType','bipolar');
                [J_values_bip,J_values_bip_sem] = obj.get_timecourses('words',ops.words,...
                                                                      'condition',ops.J_condition_flag,...
                                                                      'sessions',ops.sessions,...
                                                                      'signalType','bipolar');
            end
        end

        % format into cell
        if ~isempty(obj.bip_elec_data)
            S_values = [{S_values_uni},{S_values_bip}];
            S_values_sem = [{S_values_uni_sem},{S_values_bip_sem}];
            N_values = [{N_values_uni},{N_values_bip}];
            N_values_sem = [{N_values_uni_sem},{N_values_bip_sem}];
            W_values = []; W_values_sem = [];
            J_values = []; J_values_sem = [];
            if ops.W_condition_flag % SWJN
                W_values = [{W_values_uni},{W_values_bip}];
                W_values_sem = [{W_values_uni_sem},{W_values_bip_sem}];
                J_values = [{J_values_uni},{J_values_bip}];
                J_values_sem = [{J_values_uni_sem},{J_values_bip_sem}];
            end

        else
            S_values = {S_values_uni};
            S_values_sem = {S_values_uni_sem};
            N_values = {N_values_uni};
            N_values_sem = {N_values_uni_sem};
            W_values = []; W_values_sem = [];
            J_values = []; J_values_sem = [];
            if ops.W_condition_flag % SWJN
                W_values = {W_values_uni};
                W_values_sem = {W_values_uni_sem};
                J_values = {J_values_uni};
                J_values_sem = {J_values_uni_sem};
            end
        end

        % clear variables
        clearvars S_values_uni S_values_bip S_values_uni_sem S_values_bip_sem
        clearvars N_values_uni N_values_bip N_values_uni_sem N_values_bip_sem
        clearvars W_values_uni W_values_bip W_values_uni_sem W_values_bip_sem
        clearvars J_values_uni J_values_bip J_values_uni_sem J_values_bip_sem
        
        % -- AVERAGED BY WORD ---
        [S_values_ave,S_values_ave_sem,nTrialsS] = obj.get_word_averages('words',ops.words,...
                                                            'condition',ops.S_condition_flag,...
                                                            'sessions',ops.sessions,...
                                                            'allElecs',ops.bipolarByShank);
        [N_values_ave,N_values_ave_sem,nTrialsN] = obj.get_word_averages('words',ops.words,...
                                                            'condition',ops.N_condition_flag,...
                                                            'sessions',ops.sessions,...
                                                            'allElecs',ops.bipolarByShank);
        W_values_ave = []; W_values_ave_sem = []; nTrialsW = [];
        J_values_ave = []; J_values_ave_sem = []; nTrialsJ = [];
        if ops.W_condition_flag % SWJN
            [W_values_ave,W_values_ave_sem,nTrialsW] = obj.get_word_averages('words',ops.words,...
                                                                'condition',ops.W_condition_flag,...
                                                                'sessions',ops.sessions);
            [J_values_ave,J_values_ave_sem,nTrialsJ] = obj.get_word_averages('words',ops.words,...
                                                                'condition',ops.J_condition_flag,...
                                                                'sessions',ops.sessions);
        end

        
        % create structure with all data 
        data = struct;
        data.S_values         = S_values;
        data.W_values         = W_values;
        data.J_values         = J_values;
        data.N_values         = N_values;
        data.S_values_sem     = S_values_sem;
        data.W_values_sem     = W_values_sem;
        data.J_values_sem     = J_values_sem;
        data.N_values_sem     = N_values_sem;
        data.S_values_ave     = S_values_ave;
        data.W_values_ave     = W_values_ave;
        data.J_values_ave     = J_values_ave;
        data.N_values_ave     = N_values_ave;
        data.S_values_ave_sem = S_values_ave_sem;
        data.W_values_ave_sem = W_values_ave_sem;
        data.J_values_ave_sem = J_values_ave_sem;
        data.N_values_ave_sem = N_values_ave_sem;
        data.nTrialsS         = nTrialsS;
        data.nTrialsW         = nTrialsW;
        data.nTrialsJ         = nTrialsJ;
        data.nTrialsN         = nTrialsN


        % -- PLOT --
        if ops.bipolarByShank
            obj.plot_timecourse_bip_by_shank(data,'words',ops.words);
        else
            obj.plot_timecourse(data,'subAverage',ops.subAverage,'words',ops.words);
            if ops.subAverage % only plot barplots with z-scored signal
                obj.plot_barplot(data,'subAverage',true,'words',ops.words);
                obj.plot_barplot(data,'subAverage',false,'words',ops.words);
            end
        end

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GET TIMECOURSES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [values,values_sem]=get_timecourses(obj,varargin)
        % TODO - description

        p = inputParser();
        addParameter(p,'words',[]);
        addParameter(p,'condition',[]);
        addParameter(p,'signalType','unipolar'); % 'bipolar'
        addParameter(p,'sessions',[]);
        addParameter(p,'useLangElecs',true);
        addParameter(p,'allElecs',false); % separate from useLangElecs
        addParameter(p,'split',[]); % empty 'odd' or 'even'
        addParameter(p,'average',true); % average over trials
        parse(p, varargin{:});
        ops = p.Results;

        if strcmp(ops.signalType,'unipolar')
            data = obj.elec_data;
            sig_elecs = logical(obj.s_vs_n_sig.elec_data{1,1});
            valid = obj.elec_ch_valid;
        elseif strcmp(ops.signalType,'bipolar')
            data = obj.bip_elec_data;
            sig_elecs = logical(obj.s_vs_n_sig.bip_elec_data{1,1});
            valid = obj.bip_ch_valid;
        end

        starts_n_stops = cellfun(@(x) [x.start(1), x.end(length(ops.words))],obj.trial_timing,'UniformOutput',false);

        % values is cell array (nTrials,1), within each cell (nElecs,time)
        if ops.useLangElecs && ~ops.allElecs
            values = cellfun(@(x) data(sig_elecs,x(1):x(2)),starts_n_stops,'UniformOutput',false); 
        elseif ops.allElecs % will only run for bipolar so don't need to check if clean elecs
            values = cellfun(@(x) data(:,x(1):x(2)),starts_n_stops,'UniformOutput',false);
        else
            values = cellfun(@(x) data((~sig_elecs & valid),x(1):x(2)),starts_n_stops,'UniformOutput',false);
        end

        if ~isempty(ops.sessions)
            values = values(strcmp(obj.condition,ops.condition) & any(obj.session==ops.sessions,2));
        else
            values = values(strcmp(obj.condition,ops.condition));
        end

        % split the data if specified
        if strcmp(ops.split,'odd')
            values = values(1:2:length(values)); 
        elseif strcmp(ops.split,'even')
            values = values(2:2:length(values));
        end

        % average over trials
        if ops.average
            values_done = zeros(size(values{1,1},1),size(values{1,1},2)); % nElecs x time
            values_sem = zeros(size(values{1,1},1),size(values{1,1},2));
            for i=1:size(values{1,1},1)
                tmp = zeros(length(values),size(values{1,1},2));
                for j=1:length(values)
                    tmp(j,:) = values{j,1}(i,:);
                end
                values_done(i,:) = mean(tmp,1);
                values_sem(i,:) = std(tmp,1)/sqrt(size(tmp,1));
            end
        % return all trials
        else
            values_done = zeros(length(values)*size(values{1,1},1),size(values{1,1},2)); % nElecs*nTrials x time
            values_sem = []; % no averaging so no SEM
            idx_curr = 1;
            for i=1:size(values{1,1},1)
                tmp = zeros(length(values),size(values{1,1},2));
                for j=1:length(values)
                    tmp(j,:) = values{j,1}(i,:);
                end
                values_done(idx_curr:idx_curr+length(values)-1,:) = tmp(:,:);
                idx_curr = idx_curr + length(values);
            end
        end
        values = values_done;

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GET WORD AVERAGES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [values_ave,values_ave_sem,nTrials]=get_word_averages(obj,varargin)
        % TODO - description

        p = inputParser();
        addParameter(p,'words',[]);
        addParameter(p,'condition',[]);
        addParameter(p,'sessions',[]);
        addParameter(p,'useLangElecs',true);
        addParameter(p,'allElecs',false); % separate from useLangElecs
        parse(p, varargin{:});
        ops = p.Results;

        if ~isempty(ops.sessions)
            keep_trials = any(obj.session==ops.sessions,2);
        else
            keep_trials = [];
        end
        
        [~,cond_table] = obj.get_ave_cond_trial('words',ops.words,'condition',ops.condition,'keep_trials',keep_trials);
        values_ave = table2cell(cond_table(:,~ismember(cond_table.Properties.VariableNames,{'key','string'})));
            
        % --- UNCOMMENT IF DESIRED ---
        % EVEN trials ONLY! (odd trials were used for chan identification)
        % values_ave = cellfun(@(x) x(:,2:2:size(x,2),:),values_ave,'UniformOutput',false);

        % average over trials
        nTrials = size(values_ave{1},2);
        values_ave_sem = cellfun(@(x) squeeze(std(x,[],2))/sqrt(nTrials),values_ave,'UniformOutput',false);
        values_ave = cellfun(@(x) squeeze(mean(x,2)),values_ave,'UniformOutput',false);

        % subset electrodes
        sig_elecs = table2cell(obj.s_vs_n_sig(:,~ismember(obj.s_vs_n_sig.Properties.VariableNames,{'key'})));
        if ops.useLangElecs && ~ops.allElecs
            elecs_use = sig_elecs;
        elseif ops.allElecs  
            elecs_use = cellfun(@(x) logical(ones(size(x))),sig_elecs,'UniformOutput',false);
        else
            elecs_use = cellfun(@(x) ~x,sig_elecs,'UniformOutput',false);
        end
        values_ave = cellfun(@(x,y) x(y,:),values_ave,elecs_use,'UniformOutput',false);
        values_ave_sem = cellfun(@(x,y) x(y,:),values_ave_sem,elecs_use,'UniformOutput',false);

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUMMARY STATISTICS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function summary=get_summary_statistics(obj,varargin)
        % Returns table of summary statistics.

        p = inputParser();
        addParameter(p,'sessions',[]); % sessions that are analyzed
        parse(p, varargin{:});
        ops = p.Results;

        if ~isempty(ops.sessions)
            keep_trials = find(any(obj.session==ops.sessions,2));
        else
            keep_trials = 1:length(obj.session);
        end

        columns = {'subject',...
                   'num_sig',...
                   'num_sig_left',...
                   'num_sig_right',...
                   'num_clean',...
                   'num_clean_left',...
                   'num_clean_right',...
                   'num_prelim_deselect',...
                   'num_w_IED',...
                   'num_w_noise',...
                   'num_user_deselect',...
                   'num_total_ecog_seeg',...
                   'num_total_ecog_seeg_left',...
                   'num_total_ecog_seeg_right',...
                   'num_total',...
                   'percent_sig',...
                   'chans_sig',...
                   'names_sig',...
                   'p_values_sig',...
                   'has_bipolar',...
                   'num_sig_bipolar',...
                   'num_total_bipolar',...
                   'percent_sig_bipolar',...
                   'chans_sig_bipolar',...
                   'names_sig_bipolar',...
                   'p_values_sig_bipolar',...
                   'native_sample_freq',...
                   'elecs_per_amp',...
                   'runs_analyzed',...
                   'num_runs_total',...
                   'num_words',...
                   'presentation_rate',...
                   'num_trials_S',...
                   'num_trials_N',...
                   'num_trials_W',...
                   'num_trials_J'...
        };

        % removed channel 193 from clean channels in BJH006
        % need to repreprocess because skull EEG channels weren't prelim removed
        % but for now need to manually remove this channel (11/3/2022)
        if strcmp(obj.subject,'BJH006')
            obj.elec_ch_clean = obj.elec_ch_clean(obj.elec_ch_clean~=193);
        end

        % --- unipolar ---
        mapped_hemi = cellfun(@(x) obj.anatomy.hemisphere(x),... % TODO - make mapping cleaner
                               obj.anatomy.mapping,'UniformOutput',false);
        num_sig = {sum(obj.s_vs_n_sig.elec_data{1})};
        num_sig_left = {sum(cell2mat(arrayfun(@(x) strcmp(mapped_hemi{x}{1},'left'),...
                        find(obj.s_vs_n_sig.elec_data{1}),'UniformOutput',false)))};
        num_sig_right = {sum(cell2mat(arrayfun(@(x) strcmp(mapped_hemi{x}{1},'right'),...
                         find(obj.s_vs_n_sig.elec_data{1}),'UniformOutput',false)))};
        num_clean = {length(obj.elec_ch_clean)};
        num_clean_left = {sum(cell2mat(arrayfun(@(x) strcmp(mapped_hemi{x}{1},'left'),...
                          obj.elec_ch_clean,'UniformOutput',false)))};
        num_clean_right = {sum(cell2mat(arrayfun(@(x) strcmp(mapped_hemi{x}{1},'right'),...
                           obj.elec_ch_clean,'UniformOutput',false)))};
        num_prelim_deselect = {length(obj.elec_ch_prelim_deselect)};
        num_w_IED = {length(obj.elec_ch_with_IED)};
        num_w_noise = {length(obj.elec_ch_with_noise)};
        num_user_deselect = {length(obj.elec_ch_user_deselect)};
        num_total_ecog_seeg = {sum(cell2mat(cellfun(@(x) (contains(x,'ecog') | contains(x,'seeg')),...
                                                  obj.elec_ch_type,'UniformOutput',false)))};
        num_total_ecog_seeg_left = {sum(cell2mat(cellfun(@(x,y) ((contains(x,'ecog') | contains(x,'seeg')) & ...
                                   strcmp(y,'left')),obj.elec_ch_type,mapped_hemi,'UniformOutput',false)))};
        num_total_ecog_seeg_right = {sum(cell2mat(cellfun(@(x,y) ((contains(x,'ecog') | contains(x,'seeg')) & ...
                                    strcmp(y,'right')),obj.elec_ch_type,mapped_hemi,'UniformOutput',false)))};
        num_total = {size(obj.elec_data,1)};
        percent_sig = {round((num_sig{1}/num_clean{1})*100,2)};
        chans_sig = {find(obj.s_vs_n_sig.elec_data{1})};
        names_sig = {obj.elec_ch_label(obj.s_vs_n_sig.elec_data{1})};
        ratio_tmp = obj.s_vs_n_p_ratio.elec_data{1};
        p_values_sig = {ratio_tmp(obj.s_vs_n_sig.elec_data{1})};

        % --- bipolar ---
        has_bipolar = ~isempty(obj.bip_elec_data) && sum(obj.s_vs_n_sig.bip_elec_data{1})>0;
        if has_bipolar
            num_sig_bipolar = {sum(obj.s_vs_n_sig.bip_elec_data{1})};
            num_total_bipolar = {length(obj.bip_ch_label)};
            percent_sig_bipolar = {round((num_sig_bipolar{1}/num_total_bipolar{1})*100,2)};
            chans_sig_bipolar = {find(obj.s_vs_n_sig.bip_elec_data{1})};
            names_sig_bipolar = {obj.bip_ch_label(obj.s_vs_n_sig.bip_elec_data{1})};
            ratio_tmp = obj.s_vs_n_p_ratio.bip_elec_data{1};
            p_values_sig_bipolar = {ratio_tmp(obj.s_vs_n_sig.bip_elec_data{1})};
        else
            num_sig_bipolar = {};
            num_total_bipolar = {};
            percent_sig_bipolar = {};
            chans_sig_bipolar = {};
            names_sig_bipolar = {};
            p_values_sig_bipolar = {};
        end

        % --- other ---
        native_sample_freq = {obj.for_preproc.sample_freq_raw};
        elecs_per_amp = {obj.for_preproc.elecs_per_amp};
        if ops.sessions
            runs_analyzed = {ops.sessions};
        else
            runs_analyzed = {length(obj.stitch_index)};
        end
        num_runs_total =  {length(obj.stitch_index)};
        num_words = sum(cell2mat(cellfun(@(x) contains(x,'word'),obj.trial_timing{keep_trials(1),1}.key,'UniformOutput',false)));
        presentation_rate = {(obj.trial_timing{keep_trials(1),1}.end(1) - ...
                            obj.trial_timing{keep_trials(1),1}.start(1)+1) / obj.sample_freq};
        if strcmp(obj.experiment,'MITSWJNTask')
            num_trials_S = {sum(cell2mat(cellfun(@(x) strcmp(x,'SENTENCES'),obj.condition(keep_trials),'UniformOutput',false)))};
            num_trials_N = {sum(cell2mat(cellfun(@(x) strcmp(x,'NONWORDS'),obj.condition(keep_trials),'UniformOutput',false)))};
            num_trials_J = {sum(cell2mat(cellfun(@(x) strcmp(x,'JABBERWOCKY'),obj.condition(keep_trials),'UniformOutput',false)))};
            num_trials_W = {sum(cell2mat(cellfun(@(x) strcmp(x,'WORDS'),obj.condition(keep_trials),'UniformOutput',false)))};
        else
            num_trials_S = {sum(cell2mat(cellfun(@(x) strcmp(x,'Sentences'),obj.condition(keep_trials),'UniformOutput',false)))};
            num_trials_N = {sum(cell2mat(cellfun(@(x) strcmp(x,'Jabberwocky'),obj.condition(keep_trials),'UniformOutput',false)))};
            num_trials_J = {};
            num_trials_W = {};
        end

        % return table with one row
        new_row = {obj.subject,...
                   num_sig,...
                   num_sig_left,...
                   num_sig_right,...
                   num_clean,...
                   num_clean_left,...
                   num_clean_right,...
                   num_prelim_deselect,...
                   num_w_IED,...
                   num_w_noise,...
                   num_user_deselect,...
                   num_total_ecog_seeg,...
                   num_total_ecog_seeg_left,...
                   num_total_ecog_seeg_right,...
                   num_total,...
                   percent_sig,...
                   chans_sig,...
                   names_sig,...
                   p_values_sig,...
                   has_bipolar,...
                   num_sig_bipolar,...
                   num_total_bipolar,...
                   percent_sig_bipolar,...
                   chans_sig_bipolar,...
                   names_sig_bipolar,...
                   p_values_sig_bipolar,...
                   native_sample_freq,...
                   elecs_per_amp,...
                   runs_analyzed,...
                   num_runs_total,...
                   num_words,...
                   presentation_rate,...
                   num_trials_S,...
                   num_trials_N,...
                   num_trials_W,...
                   num_trials_J...
        };
        summary = cell2table(new_row,'VariableNames',columns);

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT TIMECOURSE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_timecourse(obj,varargin)

        p = inputParser();
        addRequired(p,'data');
        addParameter(p,'subAverage',false);
        addParameter(p,'elec_flag','unipolar');
        addParameter(p,'doBipolar',true);
        addParameter(p,'words',1:12);
        addParameter(p,'noDisplay',true);
        parse(p, varargin{:});
        ops = p.Results;


        % variable definition
        if ~isempty(ops.data.W_values) % SWJN
            conditions = {'S','W','J','N'};
            labels = {'Sentences','Words','Jabberwocky','Nonwords'};
            colors = {'r',[0.4660 0.6740 0.1880],[1 0.725 0],'b'};
        else
            conditions = {'S','N'};
            labels = {'Sentences','Nonwords'};
            colors = {'r','b'};
        end

        % extract channel names
        if strcmp(ops.elec_flag,'unipolar')
            chan_names = obj.elec_ch_label(logical(obj.s_vs_n_sig.elec_data{1,1}));
            chan_nums = obj.elec_ch(logical(obj.s_vs_n_sig.elec_data{1,1}));
            idx = 1;
        elseif strcmp(ops.elec_flag,'bipolar')
            chan_names = obj.bip_ch_label(logical(obj.s_vs_n_sig.bip_elec_data{1,1}));
            chan_nums = obj.bip_ch(logical(obj.s_vs_n_sig.bip_elec_data{1,1}));
            idx = 2;
        else
            error('Signal flag not recognized');
        end
        chan_names = cellfun(@(x) strrep(x,'_',''),chan_names,'UniformOutput',false);
        
        % number of plots to produce
        if ops.subAverage
            nPlots = 1;
        else % all lang channels
            nPlots = size(ops.data.S_values{idx},1);
        end


        for i=1:nPlots

            % initialize plotting
            close all
            if ops.noDisplay
                set(0, 'DefaultFigureVisible', 'off')
            end
            f = figure; 
            hold on;
            set(f,'position',[1123 29 3000 1275])

            % vertical dashed lines at word onset
            nSamples = length(ops.data.S_values{idx});
            x = ((1:nSamples/length(ops.words):nSamples))/obj.sample_freq;
            for j=1:length(x)
                xline(x(j),'--');
            end
            xlim([0 nSamples/obj.sample_freq]);
            % ylim([-.2 .5]); 
            

            % --- AVERAGED BY WORD ---
            % xtick at halfway point between words
            x = ((1:nSamples/length(ops.words):nSamples) + nSamples/length(ops.words)/2) / obj.sample_freq;
            xticks(x);
            xticklabels(ops.words);
            % loop through conditions
            for j=1:length(conditions)
                if ops.subAverage
                    y = eval(strcat('ops.data.',conditions{j},'_values_ave{idx};'));
                    y_sem = std(y,[],1) ./ sqrt(size(y,1)); % over lang channels
                    y = mean(y,1);
                else
                    y = eval(strcat('ops.data.',conditions{j},'_values_ave{idx}(i,:);'));
                    y_sem = eval(strcat('ops.data.',conditions{j},'_values_ave_sem{idx}(i,:);')); % over trials
                end

                b1(j) = obj.plot_condition(x,y,y_sem,colors{j},labels{j},'isAveraged',true); 

            end


            % --- TIMECOURSE ---
            x = (1:nSamples) / obj.sample_freq;
            % loop through conditions
            for j=1:length(conditions)
                if ops.subAverage
                    y = eval(strcat('ops.data.',conditions{j},'_values{idx};'));
                    y_sem = std(y,[],1) ./ sqrt(size(y,1)); % over lang channels
                    y = mean(y,1);
                else
                    y = eval(strcat('ops.data.',conditions{j},'_values{idx}(i,:);'));
                    y_sem = eval(strcat('ops.data.',conditions{j},'_values_sem{idx}(i,:);')); % over trials
                end

                b2 = obj.plot_condition(x,y,y_sem,colors{j},labels{j},'isAveraged',false); % b2 will go unused

            end
        

            % other plotting parameters
            set(gca,'FontSize',18,'LineWidth',1.5,'Box','off');
            legend(b1,'Location','northwest','NumColumns',2,'FontSize',24)
            xlabel('Word Position');
            ylabel('High Gamma Envelope (a.u.)');
            if ops.subAverage
                subtxt = sprintf('Average of language-responsive channels (n = %d, %s)',size(ops.data.S_values{idx},1),ops.elec_flag);
                analysis_path = strcat(obj.langloc_save_path,'lang_electrodes',filesep);
                fname = sprintf('%s_average_s_v_n_%s_timecourse.png',obj.subject,ops.elec_flag);
            else
                subtxt = sprintf('%s (chan%d, nTrials = %d, %s)',chan_names{i},chan_nums(i),ops.data.nTrialsS,ops.elec_flag);
                analysis_path = strcat(obj.langloc_save_path,'lang_electrodes',filesep,obj.subject,filesep,'timecourse');
                fname = sprintf('%s_%s_%s.png',obj.subject,ops.elec_flag,chan_names{i});
            end
            title(obj.subject);
            subtitle(subtxt)
        

            % --- SAVE ---  
            if ~exist(analysis_path)
                mkdir(analysis_path);
            end   
            saveas(gcf,strcat(analysis_path,'/',fname));
            set(0, 'CurrentFigure', f);
            close(f);

        end

        % bipolar
        if ~isempty(obj.bip_elec_data) & ops.doBipolar
            obj.plot_timecourse(ops.data,'subAverage',ops.subAverage,'words',ops.words,'elec_flag','bipolar','doBipolar',false);
        end

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT TIMECOURSE BIPOLAR ELECTRODES BY SHANK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_timecourse_bip_by_shank(obj,varargin)
        % Plots grid of S vs N plots for bipolar electrodes by shank

        p = inputParser();
        addRequired(p,'data');
        addParameter(p,'words',1:12);
        addParameter(p,'noDisplay',true);
        parse(p, varargin{:});
        ops = p.Results;


        % variable definition
        conditions = {'S','N'};
        labels = {'Sentences','Nonwords'};
        colors = {'r','b'};

        chan_names = obj.bip_ch_label;
        chan_nums = obj.bip_ch;
        isLang = logical(obj.s_vs_n_sig.bip_elec_data{1,1});
        idx = 2; % bipolar
        chan_names = cellfun(@(x) strrep(x,'_',''),chan_names,'UniformOutput',false);

        % number of plots to produce
        nPlots = size(ops.data.S_values{idx},1);

        % placement of plots in grid
        prev_shank_name = 'PREV';
        plot_placement = zeros(1,nPlots); % value corresponds to row in grid, new shank starts with 1
        for i=1:nPlots
            individ_elecs = split(chan_names{i},'-');
            assert(length(individ_elecs)==2);
            shank_name = extract(individ_elecs,lettersPattern);
            assert(strcmp(shank_name{1},shank_name{2}));
            if ~strcmp(shank_name{1},prev_shank_name)
                row = 1;
            end   
            plot_placement(1,i) = row;
            row = row + 1;
            prev_shank_name = shank_name{1};
        end

        % size of grid
        nSubplots_x = sum(plot_placement==1);
        nSubplots_y = max(plot_placement);

        % initialize plotting
        close all
        if ops.noDisplay
            set(0, 'DefaultFigureVisible', 'off')
        end
        f = figure; 
        hold on;
        set(f,'position',[10 10 5000 2800]);

        fprintf(1,'\n[');
        shank_num = 0;
        assert(length(isLang)==nPlots);
        for i=1:nPlots

            % subplot
            if plot_placement(1,i)==1
                shank_num = shank_num + 1; 
            end
            ax = subplot(nSubplots_y,nSubplots_x,((plot_placement(1,i)-1)*nSubplots_x)+shank_num);
            hold on;

            % vertical dashed lines at word onset
            nSamples = length(ops.data.S_values{idx});
            x = ((1:nSamples/length(ops.words):nSamples))/obj.sample_freq;
            % for j=1:length(x)
            %     xline(x(j),'--');
            % end
            xlim([0 nSamples/obj.sample_freq]);
            % ylim([-.2 .5]); 
            

            % --- AVERAGED BY WORD ---
            % xtick at halfway point between words
            x = ((1:nSamples/length(ops.words):nSamples) + nSamples/length(ops.words)/2) / obj.sample_freq;
            xticks(x);
            xticklabels(ops.words);
            % loop through conditions
            for j=1:length(conditions)
                y = eval(strcat('ops.data.',conditions{j},'_values_ave{idx}(i,:);'));
                y_sem = eval(strcat('ops.data.',conditions{j},'_values_ave_sem{idx}(i,:);')); % over trials

                b1(j) = obj.plot_condition(x,y,y_sem,colors{j},labels{j},'isAveraged',true,...
                                    'linewidthAveraged',2,'markersizeAveraged',7,'linewidthErrorAveraged',1); 

            end


            % --- TIMECOURSE ---
            x = (1:nSamples) / obj.sample_freq;
            % loop through conditions
            % for j=1:length(conditions)
            %     y = eval(strcat('ops.data.',conditions{j},'_values{idx}(i,:);'));
            %     y_sem = eval(strcat('ops.data.',conditions{j},'_values_sem{idx}(i,:);')); % over trials

            %     b2 = obj.plot_condition(x,y,y_sem,colors{j},labels{j},'isAveraged',false); % b2 will go unused

            % end

            % make lang elec axes bold
            if isLang(i)
                set(gca,'FontSize',18,'LineWidth',4,'Box','on');
            else
                set(gca,'FontSize',18,'LineWidth',1.5,'Box','off');
            end

            % other plotting parameters
            set(gca,'xticklabels',[]);
            % legend(b1,'Location','northwest','NumColumns',2,'FontSize',24)
            % xlabel('Word Position');
            % ylabel('High Gamma Envelope (a.u.)');

            fprintf(1,'.');

        end

        analysis_path = strcat(obj.langloc_save_path,'bipolar_by_shank',filesep);
        fname = sprintf('%s_%s.png',obj.subject,'bipolar_electrodes_by_shank');
        sgtitle(obj.subject,'fontsize',24);
        fprintf(1,'] done\n');

        % --- SAVE ---  
        if ~exist(analysis_path)
            mkdir(analysis_path);
        end   
        saveas(gcf,strcat(analysis_path,'/',fname));
        set(0, 'CurrentFigure', f);
        close(f);

    end

    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT BARPLOT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_barplot(obj,varargin)

        p = inputParser();
        addRequired(p,'data');
        addParameter(p,'subAverage',false);
        addParameter(p,'elec_flag','unipolar');
        addParameter(p,'doBipolar',true);
        addParameter(p,'words',1:12);
        addParameter(p,'noDisplay',true);
        parse(p, varargin{:});
        ops = p.Results;


        % variable definition
        if ~isempty(ops.data.W_values) % SWJN
            conditions = {'S','W','J','N'};
            labels = {'Sentences','Words','Jabberwocky','Nonwords'};
            colors = {'r',[0.4660 0.6740 0.1880],[1 0.725 0],'b'};
        else
            conditions = {'S','N'};
            labels = {'Sentences','Nonwords'};
            colors = {'r','b'};
        end

        % extract channel names
        if strcmp(ops.elec_flag,'unipolar')
            chan_names = obj.elec_ch_label(logical(obj.s_vs_n_sig.elec_data{1,1}));
            chan_nums = obj.elec_ch(logical(obj.s_vs_n_sig.elec_data{1,1}));
            idx = 1;
        elseif strcmp(ops.elec_flag,'bipolar')
            chan_names = obj.bip_ch_label(logical(obj.s_vs_n_sig.bip_elec_data{1,1}));
            chan_nums = obj.bip_ch(logical(obj.s_vs_n_sig.bip_elec_data{1,1}));
            idx = 2;
        else
            error('Signal flag not recognized');
        end
        chan_names = cellfun(@(x) strrep(x,'_',''),chan_names,'UniformOutput',false);
        
        % number of plots to produce
        if ops.subAverage
            nPlots = 1;
        else % all lang channels
            nPlots = size(ops.data.S_values{idx},1);
        end


        for i=1:nPlots

            % initialize plotting
            close all
            if ops.noDisplay
                set(0, 'DefaultFigureVisible', 'off')
            end
            f = figure; 
            hold on;
            set(f,'position',[1123 29 500 600])

            % loop through conditions
            for j=1:length(conditions)
                if ops.subAverage
                    y = eval(strcat('mean(ops.data.',conditions{j},'_values_ave{idx},2);'));
                    y_sem = std(y,[],1) ./ sqrt(size(y,1)); % over lang channels
                    y = mean(y,1);
                else
                    y = eval(strcat('mean(ops.data.',conditions{j},'_values_ave{idx}(i,:),2);'));
                    y_sem = eval(strcat('mean(ops.data.',conditions{j},'_values_ave_sem{idx}(i,:),2);')); % over trials
                end

                % bar plot
                h(j) = bar(j,y,'displayname',labels{j});
                set(h(j), 'FaceColor', colors{j})
                errorbar(j,y,y_sem,'k','LineWidth',1.5,'Capsize',0);
            end
        
            % other plotting parameters
            ax1 = gca;                
            ax1.XAxis.Visible = 'off'; 
            ylim([-0.5,1])    
            set(gca,'FontSize',10,'LineWidth',1.5,'Box','off');
            legend(h,'FontSize',10,'Location','southwest');
            legend boxoff
            ylabel('Z-Scored High Gamma Envelope (a.u.)','FontSize',10)
            if ops.subAverage
                subtxt = sprintf('Average of language-responsive channels (n = %d, %s)',size(ops.data.S_values{idx},1),ops.elec_flag);
                analysis_path = strcat(obj.langloc_save_path,'lang_electrodes',filesep);
                fname = sprintf('%s_average_s_v_n_%s_barplot.png',obj.subject,ops.elec_flag);
            else
                subtxt = sprintf('%s (chan%d, nTrials = %d, %s)',chan_names{i},chan_nums(i),ops.data.nTrialsS,ops.elec_flag);
                analysis_path = strcat(obj.langloc_save_path,'lang_electrodes',filesep,obj.subject,filesep,'barplot');
                fname = sprintf('%s_%s_%s.png',obj.subject,ops.elec_flag,chan_names{i});
            end
            title(obj.subject);
            subtitle(subtxt,'FontSize',8)
        

            % --- SAVE ---  
            if ~exist(analysis_path)
                mkdir(analysis_path);
            end   
            saveas(gcf,strcat(analysis_path,'/',fname));
            set(0, 'CurrentFigure', f);
            close(f);

        end

        % bipolar
        if ~isempty(obj.bip_elec_data) & ops.doBipolar
            obj.plot_barplot(ops.data,'subAverage',ops.subAverage,'words',ops.words,'elec_flag','bipolar','doBipolar',false);
        end

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT CONDITION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function b=plot_condition(obj,varargin)

        p = inputParser();
        addRequired(p,'x');
        addRequired(p,'y');
        addRequired(p,'y_sem');
        addRequired(p,'c'); % color
        addRequired(p,'label');
        addParameter(p,'isAveraged',false); 
        addParameter(p,'linewidthAveraged',9);
        addParameter(p,'markersizeAveraged',20);
        addParameter(p,'linewidthErrorAveraged',3);
        parse(p, varargin{:});
        ops = p.Results;


        % --- PARAMETERS ---
        if ops.isAveraged
            % plot params
            color           = ops.c;
            linewidth       = ops.linewidthAveraged;
            marker          = 's';
            markersize      = ops.markersizeAveraged;
            markerfacecolor = ops.c;
            markeredgecolor = ops.c;
            label           = ops.label;
            % error params
            linestyle_err = 'none';
            linewidth_err = ops.linewidthErrorAveraged;
            capsize_err   = 0;
        else
            % plot params
            color           = ops.c;
            linewidth       = 2;
            marker          = 'default';
            markersize      = 'default';
            markerfacecolor = 'default';
            markeredgecolor = 'default';
            label           = '';
            % error params
            linestyle_err = 'none';
            facealpha_err = 0.1;
        end


        % --- SIGNAL ---
        b = plot(ops.x,ops.y,...
                        'Color',color,...
                        'LineWidth',linewidth,...
                        'Marker',marker,...
                        'MarkerSize',markersize,...
                        'MarkerFaceColor',markerfacecolor,...
                        'MarkerEdgeColor',markerfacecolor,...
                        'DisplayName',label...
        );


        % --- ERROR ---
        if ops.isAveraged
            errorbar(ops.x,ops.y,ops.y_sem,...
                     'Color',color,...
                     'LineStyle',linestyle_err,...
                     'LineWidth',linewidth_err,...
                     'CapSize',capsize_err...
            );
        else
            patch([ops.x fliplr(ops.x)], [ops.y+ops.y_sem fliplr(ops.y-ops.y_sem)],...
                  color,...
                  'LineStyle',linestyle_err,...
                  'FaceAlpha',facealpha_err...
            );
        end


    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOTTING ALL ELECTRODES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_s_vs_n(obj,prev_ops,varargin)
        % TODO - Description

        p = inputParser();
        addRequired(p,'elec_flag');
        addRequired(p,'channel_labels')
        addParameter(p,'S_table',[])
        addParameter(p,'N_table',[]);
        addParameter(p,'S_N_corr_table',[]);
        addParameter(p,'S_N_corr_rnd_table',[]);
        addParameter(p,'s_vs_n_sig',[]);
        addParameter(p,'S_N_p_ratio_tbl',[]);
        addParameter(p,'noDisplay',true)
        parse(p, varargin{:});
        ops = p.Results;

        analysis_path = strcat(obj.langloc_save_path,'all_electrodes/');

        % extract data of interest
        elec_flag = ops.elec_flag;
        S_dat = ops.S_table.(elec_flag){1};
        N_dat = ops.N_table.(elec_flag){1};
        SN_corr_dat = ops.S_N_corr_table.(elec_flag){1};
        SN_cor_rnd_dat = ops.S_N_corr_rnd_table.(elec_flag){1};
        SN_sig_dat = ops.s_vs_n_sig.(elec_flag){1};
        SN_p_ratio_dat = ops.S_N_p_ratio_tbl.(elec_flag){1};

        % plotting parameters
        num_rows = 3;
        num_columns = 3;
        nbins = 50;
        total_plots = num_rows*num_columns;
        pp = 0;

        fprintf(1,'\n> Plotting channels from %s ...\n',elec_flag);
        fprintf(1,'[');

        if strcmp(ops.elec_flag,'elec_data')
            elec_flag = 'unipolar';
        elseif strcmp(ops.elec_flag,'bip_elec_data')
            elec_flag = 'bipolar';
        else
            error('Signal flag not recognized');
        end

        % plot
        close all
        if ops.noDisplay
            set(0, 'DefaultFigureVisible', 'off')
        end
        f = figure;
        set(f,'position',[1123 29 1266 1275])
            
        % go through all electrodes
        for i=1:size(S_dat,1)
            s_electrode_resp = squeeze(S_dat(i,:,:))';
            n_electrode_resp = squeeze(N_dat(i,:,:))';
            s_n_rho = SN_corr_dat(i);
            s_n_rho_rnd = SN_cor_rnd_dat(i,:);
            is_sig = SN_sig_dat(i);
            p_ratio = SN_p_ratio_dat(i);
            word_pos = repmat(1:size(s_electrode_resp,1),size(s_electrode_resp,2),1)';
              
            % skip electrodes with empty (0-ed) channels
            if isnan(s_n_rho)
                continue
            end

            sup_title = (strcat(ops.channel_labels{i,1}));

            % plot the distbirtubiotn
            ax = subplot(num_rows,num_columns,3*(i-num_rows*fix((i-1)/num_rows)-1)+1);
            h1 = histogram(s_n_rho_rnd,nbins);
            hold on
            xline(s_n_rho,'linewidth',3);
            ax.Box = 'off';
            h1.EdgeColor = 'w';
            if ~ops.noDisplay; shg; end
            ax.XAxis.LineWidth = 2;
            ax.YAxis.LineWidth = 2;
            ax.Title.String=sprintf('corr=%f,\n p_{ratio}=%0.4f sig=%d',s_n_rho,p_ratio,is_sig);

            % plot the average signal by word for both conditions
            ax = subplot(num_rows,num_columns,3*(i-num_rows*fix((i-1)/num_rows)-1)+2);
            b1 = plot(mean(word_pos,2)+.1,mean(s_electrode_resp,2),'color',[1,.5,.5],'linewidth',2,'marker','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'displayname','S');
            hold on
            b1 = plot(mean(word_pos,2)-.1,mean(n_electrode_resp,2),'color',[.5,.5,1],'linewidth',2,'marker','o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'displayname','N');
                
            y = s_electrode_resp;
            x = word_pos;
            bl = errorbar(mean(x,2)+.1,mean(y,2),std(y,[],2)./sqrt(size(y,2)));
            bl.LineStyle = 'none';
            bl.Color = [1,.5,.5];
            bl.LineWidth = 2;
            bl.CapSize = 2;
            hAnnotation = get(bl,'Annotation');
            hLegendEntry = get(hAnnotation,'LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off');
                
            y = n_electrode_resp;
            x = word_pos;
            bl = errorbar(mean(x,2)-.1,mean(y,2),std(y,[],2)./sqrt(size(y,2)));
            bl.LineStyle = 'none';
            bl.Color = [.5,.5,1];
            bl.LineWidth = 2;
            bl.CapSize = 2;
            hAnnotation = get(bl,'Annotation');
            hLegendEntry = get(hAnnotation,'LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off');
                    
            ax.XAxis.Visible = 'on';
            ax.XTick = 1:max(word_pos(:));
            ax.XLim = [0,max(word_pos(:))+1];
            all_points = [s_electrode_resp(:);n_electrode_resp(:)];
            y_quantile = quantile(all_points,10);
            h = get(ax,'children');
            ax.FontSize = 12;
            set(ax,'ydir', 'normal','box','off','ylim',[y_quantile(1),y_quantile(end)]);
            ah = get(ax,'children');
            arrayfun(@(x) set(ah(x),'DisplayName',''),[1:2]);
            ah(3).DisplayName = 'N';
            ah(4).DisplayName = 'S';
            set(ax,'children',ah);
            legend(ah(3:4),'Location','northwest','NumColumns',2)
            xlabel('word position');
            ax.YLabel.String='High Gamma (a.u.)';
            ax.XAxis.LineWidth = 2;
            ax.YAxis.LineWidth = 2;
            ax.Title.String = erase(sup_title,'_');

            % save file if last row in figure
            if ~mod(i,num_rows) | i==size(S_dat,1)
                pp = pp+1;

                if ~exist(strcat(analysis_path,obj.subject))
                    mkdir(strcat(analysis_path,obj.subject));
                end
                
                set(gcf,'PaperPosition',[.25 .25 8 6])
                set(gcf,'PaperOrientation','landscape');
                fname = sprintf('%s_s_v_n_words-%d-%d_p_%0.2f_%s_%s_%s.pdf',...
                                obj.subject,...
                                min(prev_ops.words),...
                                max(prev_ops.words),...
                                prev_ops.threshold,...
                                prev_ops.side,...
                                elec_flag,...
                                num2str(num2str(pp))...
                );
                        
                print(f, '-bestfit','-dpdf','-opengl', strcat(analysis_path,obj.subject,'/',fname));
                set(0, 'CurrentFigure', f);
                clf reset;
                    
            end    
                
            fprintf(1,'.');

        end
                
        close(f);
        
        fprintf(1,'] done\n')
               
    end

        %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FORMAT FOR CDR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [event,resp]=format_for_CDR(obj,varargin)
        % Formats stimuli / signal for CDR and saves to csv
        

        % --- EVENTS ---
        % loop over all trials
        event = {};
        subject = obj.subject;
        for t=1:length(obj.trial_timing)

            curr_trial = obj.trial_timing{t};
            trial_num = obj.events_table.trial(t);
            condition = obj.condition{t}

            % loop over all events in trial
            for e=1:size(curr_trial,1)

                stim_type = curr_trial.key{e};
                if contains(stim_type,'word')
                    word_num = extract(stim_type,digitsPattern);
                    word_num = str2num(word_num{1}); 
                    duration = length(curr_trial.start(e):curr_trial.end(e))/obj.sample_freq;

                    % add new row to event table
                    row = {{time},...
                           {subject},...
                           {}
                    };
                else
                    continue
                end
            end
        end

    end
           


end

end

