function formatData(varargin)
    % Formats data for clustering
    % NOT FUNCTIONAL - ONLY PROVIDED FOR TRANSPARENCY

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'suffix','crunched');
    addParameter(p,'doZScore',true);
    addParameter(p,'useLangElecs',true); % if false will save clean, non-lang channels
    addParameter(p,'useWandJ',true); % for MITSWJNTask only
    addParameter(p,'resampleRate',60);
    addParameter(p,'split',[]); % empty 'odd' or 'even'
    addParameter(p,'averageTrials',true); % average condition responses over trials
    addParameter(p,'alignPresentationRate',true);
    addParameter(p,'signalType','unipolar'); 
    parse(p, varargin{:});
    ops = p.Results;


    % -------------------
    %% --- INITIALIZE ---
    % -------------------

    % paths
    [CLUSTER_PATH,SAVE_PATH] = initialize(ops.saveName);
    CRUNCHED_PATH = '/PATH/TO/CRUNCHED/DATAFILES/';

    % list of crunched files
    d = dir([CLUSTER_PATH '*.mat']); % empty structure (i.e., shouldn't find any files)
    for i=1:length(CRUNCHED_PATH)
        d_curr = dir([CRUNCHED_PATH{i} '*' ops.suffix '.mat']);
        d = [d; d_curr];
    end
    
    % format data files into folder/name
    d_files = transpose(arrayfun(@(x) {strcat(d(x).folder,filesep,d(x).name)},1:length(d)));

    % load subject session map
    if strcmp(ops.experiment,'both') | strcmp(ops.experiment,'MITLangloc')
        sub_sess_map = subjectSessionMap();
    end

    % experiments
    if strcmp(ops.experiment,'both')
        ops.experiment = {'MITSWJNTask','MITLangloc'};
        expt_suffix = 'MITSWJNTask_MITLangloc';
    else
        expt_suffix = ops.experiment;
    end

    % file naming conventions
    if ops.split
        split_string = ['_' ops.split];
    else
        split_string = '';
    end
    if ops.useLangElecs
        elecType = 'langElecs';
    else
        elecType = 'nonLangElecs';
    end
    if ~ops.averageTrials
        trial_string = '_all_trials';
    else
        trial_string = '';
    end
    saveDir = [SAVE_PATH 'data' filesep];
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end

    % load max/min values for normalizing individual trials
    if ~ops.averageTrials
        saveDir = [SAVE_PATH 'data' filesep];
        try
            all_X_table_averaged = readtable([saveDir expt_suffix '_' elecType '_labels_for_clustering' split_string '.csv'],'Delimiter',',');
        catch
            error('Before extracting trial-level responses, you must extract averaged responses so that the trial responses can be normalized by electrode');
        end
    end

    % -------------------
    %% --- FORMATTING ---
    % -------------------

    all_X = {}; % matrix to save 
    all_X_table = {}; % mapping of electrodes to subject
    if ops.averageTrials
        columns = {'subject','channel_number','channel_name'};
    else
        columns = {'subject','channel_number','channel_name','stimulus'};
    end

    % go through subjects
    for i=1:length(d_files)

        % extract subject name from file name
        subj = split(d_files{i},'/'); 
        subj = split(subj{end},'_'); 
        subj = subj{1};
        fprintf(1,'\nFormatting %s, subject %d of %d ...\n',subj,i,length(d_files));

        % load object
        load(d_files{i}); % loads as 'obj' or 'sn_data'
        
        % initial variables
        signalType = ops.signalType;
        srate = obj.sample_freq;
        cond = obj.condition;
        islang = cell2mat(obj.s_vs_n_sig.elec_data); % logical array
        isnotlang = (~islang & obj.elec_ch_valid);
        chan_labels = obj.elec_ch_label;

        % experiment specific variables
        if strcmp(ops.experiment,'MITSWJNTask')
            words = [1:8];
            nConditions = 4;
            S_condition_flag = 'SENTENCES';
            N_condition_flag = 'NONWORDS';
            W_condition_flag = 'WORDS';
            J_condition_flag = 'JABBERWOCKY';
            sessions = [];
        elseif strcmp(ops.experiment,'MITLangloc')
            words = [1:8]; % actually 12 words but only want 8 to match SWJN
            nConditions = 2;
            S_condition_flag = 'Sentences';
            N_condition_flag = 'Jabberwocky'; % actually N but labeled as J in MITLangloc
            W_condition_flag = [];
            J_condition_flag = [];
            sessions = sub_sess_map(subj);
        end
        all_cond_flags = [{S_condition_flag},{W_condition_flag},{J_condition_flag},{N_condition_flag}];

        % --- FOMRAT DATA FOR CLUSTERING ---

        % z-score
        if ops.doZScore
            obj.zscore_signal();
        end

        % extract average timecourse per condition for lang electrodes
        thing = true; % TODO - remove
        if thing
            [S_data,S_sem] = obj.get_timecourses('words',words,...
                                                 'condition',S_condition_flag,...
                                                 'useLangElecs',ops.useLangElecs,...
                                                 'sessions',sessions,...
                                                 'split',ops.split,...
                                                 'average',ops.averageTrials,...
                                                 'signalType',signalType);
            [N_data,N_sem] = obj.get_timecourses('words',words,...
                                                 'condition',N_condition_flag,...
                                                 'useLangElecs',ops.useLangElecs,...
                                                 'sessions',sessions,...
                                                 'split',ops.split,...
                                                 'average',ops.averageTrials,...
                                                 'signalType',signalType);
            % S_data = movmean(S_data,obj.sample_freq/10,2);
            % N_data = movmean(N_data,obj.sample_freq/10,2);
            S_data_p = permute(S_data,[3,2,1]);
            N_data_p = permute(N_data,[3,2,1]);
            % remove 80th S-W-J trials from AMC026 whose 80th N trial wasn't presented
            if ~ops.averageTrials & (size(N_data_p,3)<size(S_data_p,3))
                S_data_p = S_data_p(:,:,1:size(N_data_p,3));
            end
            if W_condition_flag & ops.useWandJ % SWJN
                [W_data,W_sem] = obj.get_timecourses('words',words,...
                                                     'condition',W_condition_flag,...
                                                     'useLangElecs',ops.useLangElecs,...
                                                     'sessions',sessions,...
                                                     'split',ops.split,...
                                                     'average',ops.averageTrials,...
                                                     'signalType',signalType);
                [J_data,J_sem] = obj.get_timecourses('words',words,...
                                                     'condition',J_condition_flag,...
                                                     'useLangElecs',ops.useLangElecs,...
                                                     'sessions',sessions,...
                                                     'split',ops.split,...
                                                     'average',ops.averageTrials,...
                                                     'signalType',signalType);
                % W_data = movmean(W_data,obj.sample_freq/10,2);
                % J_data = movmean(J_data,obj.sample_freq/10,2);
                W_data_p = permute(W_data,[3,2,1]);
                J_data_p = permute(J_data,[3,2,1]);
                % remove 80th S-W-J trials from AMC026 whose 80th N trial wasn't presented
                if ~ops.averageTrials & (size(N_data_p,3)<size(W_data_p,3))
                    W_data_p = W_data_p(:,:,1:size(N_data_p,3));
                end
                if ~ops.averageTrials & (size(N_data_p,3)<size(J_data_p,3))
                    J_data_p = J_data_p(:,:,1:size(N_data_p,3));
                end
                data_all = [S_data_p; W_data_p; J_data_p; N_data_p]; %cond,sample,elec
            else
                data_all = [S_data_p; N_data_p]; %cond,sample,elec
            end
        end
        choose = data_all;
        choose = permute(choose,[3,1,2]); %elec,cond,sample : have to permute in case 1 electrode

        t=0:1/srate:(size(choose,3)-1)/srate;
        length_words_samples = size(choose,3);
        
        % resample
        srate_new = ops.resampleRate; % was 100 in Tamar's original script
        tmp = resample(squeeze(choose(:,1,:))',srate_new,srate);
        if size(tmp,1)==1, tmp = tmp'; end % when only one electrode
        data_choose = nan(size(choose,1),size(choose,2),size(tmp,1));
        for ii=1:size(data_choose,2)
            tmp_resamp = resample(squeeze(choose(:,ii,:))',srate_new,srate)';
            if size(tmp_resamp)==1, tmp_resamp = tmp_resamp'; end % when only one electrode
            data_choose(:,ii,:) = tmp_resamp;
        end

        % compress time for longer presentation rates
        % TODO - need to look into this for MITLangloc
        if ops.alignPresentationRate
            word_len_seconds = length_words_samples/length(words)/srate;
            if word_len_seconds ~= 0.45

                sr_factor = word_len_seconds/0.45;
                t = 0:1/srate_new:(size(data_choose,3)-1)/srate_new;
                srate_interp = srate_new/sr_factor;
                t_interp = 0:1/srate_interp:t(end);

                tmp = interp1(t,squeeze(data_choose(:,1,:))',t_interp,'linear');
                if size(tmp,1)==1, tmp = tmp'; end % when only one electrode
                data_interp = nan(size(data_choose,1),size(data_choose,2),size(tmp,1));
                for j=1:size(data_choose,2)
                    tmp_interp = interp1(t,squeeze(data_choose(:,j,:))',t_interp,'linear')';
                    if size(tmp_interp,1)==1, tmp_interp = tmp_interp'; end % when only one electrode
                    data_interp(:,j,:) = tmp_interp;
                end

                data_choose = data_interp;
            
            end
        end

        % collapse across conditions (data_for_clustering is nConds*nSample x nElecs)
        data_for_clustering = nan(size(data_choose,3)*size(data_choose,2),size(data_choose,1));
        for ii=1:size(data_choose,2)
            data_for_clustering(((ii-1)*size(data_choose,3)+1):((ii-1)*size(data_choose,3)+size(data_choose,3)),:) = squeeze(data_choose(:,ii,:))';
        end
        t_new = 0:1/srate_new:(size(data_for_clustering,1)-1)/srate_new;

        % normalization
        saveDir = [SAVE_PATH 'data' filesep];
        if ~ops.zscoreByCondition & ops.averageTrials
            dc = data_for_clustering;
            MIN = min(dc); MAX = max(dc); 
            X = (dc-MIN)./(MAX-MIN); 
            X = X'; % elecs x samples*conds
        elseif ~ops.averageTrials
            % want to normalize by *electrode* not trial
            dc = data_for_clustering;
            subj_idxs = strcmp(all_X_table_averaged.subject,subj);
            MIN = all_X_table_averaged.min(subj_idxs)';
            MAX = all_X_table_averaged.max(subj_idxs)';
            if ops.useLangElecs
                MIN = repelem(MIN,size(dc,2)/sum(islang));
                MAX = repelem(MAX,size(dc,2)/sum(islang));
            else
                MIN = repelem(MIN,size(dc,2)/sum(isnotlang));
                MAX = repelem(MAX,size(dc,2)/sum(isnotlang));   
            end
            X = (dc-MIN)./(MAX-MIN); 
            X = X'; % elecs x samples*conds
        else
            X = data_for_clustering';
        end

        % extract stimuli (only added to csv file when trials not averaged)
        if ~ops.averageTrials % include trial by trial stimulus
            all_stimuli = cell(1,length(all_cond_flags));
            for c=1:length(all_cond_flags)
                condition = all_cond_flags(c);
                if ~isempty(condition)
                    if ~isempty(sessions)
                        stimuli = obj.events_table.stimulus_string(strcmp(obj.condition,condition{1}) & any(obj.session==sessions,2));
                    else
                        stimuli = obj.events_table.stimulus_string(strcmp(obj.condition,condition{1}));
                    end
                else
                    continue
                end
                all_stimuli{c} = stimuli;
            end
            % remove W and J if no stimuli
            all_stimuli(cell2mat(cellfun(@isempty,all_stimuli,'UniformOutput',false))) = []; 
            lengths = cell2mat(cellfun(@length,all_stimuli,'UniformOutput',false)); 
            %  remove 80th S-W-J trials from AMC026 whose 80th N trial wasn't presented
            if range(lengths)>0 
                for c=1:length(all_stimuli)
                    curr_stim = all_stimuli{c};
                    all_stimuli{c} = curr_stim(1:min(lengths),:);
                end
            end
            % combine stimuli into same cell - in form of {Strial1 / Wtrial1 / Jtrial1 / Ntrial1}
            all_stimuli_combined = cell(min(lengths),1);
            for tri=1:min(lengths)
                curr_stim = [];
                for con=1:length(all_stimuli)
                    curr_con = all_stimuli{con};
                    curr_con_stim = curr_con{tri};
                    if ~isempty(curr_stim)
                        curr_stim = strcat(curr_stim," / ",curr_con_stim);
                    else 
                        curr_stim = curr_con_stim;
                    end
                end
                all_stimuli_combined{tri} = curr_stim;
            end 
        end

        % make table with subject --> electrode mapping 
        subject_cell = repmat({subj},size(X,1),1);
        if ops.useLangElecs
            chan_num_cell = repelem(num2cell(find(islang)),size(X,1)/sum(islang));
            chan_label_cell = repelem(chan_labels(islang),size(X,1)/sum(islang));
        else
            chan_num_cell = repelem(num2cell(find(isnotlang)),size(X,1)/sum(isnotlang));
            chan_label_cell = repelem(chan_labels(isnotlang),size(X,1)/sum(isnotlang));   
        end
        if ~ops.averageTrials % include trial by trial stimulus
            stimuli_cell = repmat(all_stimuli_combined,size(X,1)/length(all_stimuli_combined),1);
            curr_X_table = cell2table([subject_cell,chan_num_cell,chan_label_cell,stimuli_cell],'VariableNames',columns);
        else
            curr_X_table = cell2table([subject_cell,chan_num_cell,chan_label_cell],'VariableNames',columns);
            curr_X_table.min = MIN';
            curr_X_table.max = MAX';
        end
        all_X_table{i} = curr_X_table;

        % append subject data to other subjects' data
        all_X{i} = X;
    end

    % group subjects with same presentation rate
    all_X_orig = all_X;
    all_X_table_orig = all_X_table;
    lengths = cell2mat(cellfun(@(x) size(X,2),all_X_orig,'UniformOutput',false));
    unique_lengths = unique(lengths);
    for i=1:length(unique_lengths)
        if unique_lengths(i)==0, continue; end
        mask = lengths==unique_lengths(i);
        all_X = cell2mat(all_X_orig(mask)');
        all_X_table = all_X_table_orig(mask)';
        all_X_table = vertcat(all_X_table{:});
        
        if ~ops.alignPresentationRate
            rate_string = ['_' num2str(round((size(all_X,2)/(ops.resampleRate*nConditions*length(words)))*1000,3)) 'ms'];
        else
            rate_string = '';
        end
           
        filename_table = [saveDir expt_suffix '_' elecType '_labels_for_clustering' split_string trial_string rate_string '.csv'];
        filename_matrix = [saveDir expt_suffix '_' elecType '_data_for_clustering' split_string trial_string rate_string '.mat'];
        writetable(all_X_table,filename_table);
        save(filename_matrix,'all_X','-v7.3');
    end

end