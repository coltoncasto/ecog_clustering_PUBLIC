function crunch(varargin)

    p = inputParser();
    addRequired(p,'subject');
    addRequired(p,'experiment');
    addParameter(p,'fromScratch',false);
    addParameter(p,'preEnvelopeExtraction',false);
    addParameter(p,'addAnatomy',false); % for when adding anatomy after
    addParameter(p,'decimation_factor',4);
    addParameter(p,'deselectElecsPrelim',true);
    addParameter(p,'doneVisualInspection',false);
    addParameter(p,'order','defaultECOG');
    addParameter(p,'isPlotVisible',true);
    parse(p, varargin{:});
    ops = p.Results;

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETUP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MINDHIVE_PATH = '/mindhive/evlab/u/Shared/ECoG';
    addpath(genpath([MINDHIVE_PATH filesep 'merged_ecog_pipeline' filesep 'utils']));
    DATA_PATH = [MINDHIVE_PATH filesep 'DATA' filesep ops.experiment];
    RAW_PATH = [DATA_PATH filesep ops.subject filesep ops.experiment filesep 'ECOG001'];
    SAVE_PATH = [MINDHIVE_PATH filesep 'crunched' filesep ops.experiment filesep 'crunched_for_merged_pipeline' filesep];
    OP_INFO_PATH = [MINDHIVE_PATH filesep 'merged_ecog_pipeline' filesep 'info' filesep];
    LOG_PATH = [MINDHIVE_PATH filesep 'merged_ecog_pipeline' filesep 'logs' filesep];
    
    save_filename = [SAVE_PATH ops.subject '_' ops.experiment '_crunched.mat'];
    % save_filename = [SAVE_PATH ops.subject '_' ops.experiment '_crunched_' ops.order '.mat'];

    formatOut = 'yyyymmdd_HHMM';
    date_string = datestr(now,formatOut);
    log_filename = [LOG_PATH ops.subject '_' ops.experiment '_' ops.order '_' date_string '.txt'];
    eval(strcat("diary ",log_filename));

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START PREPROCESSING FROM SCRATCH
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ops.fromScratch

        % --------------------
        % LOAD RAW DATA FILES
        % --------------------
        d = dir([RAW_PATH filesep 'ECOG*.dat']);
        d_files = transpose(arrayfun(@(x) {strcat(d(x).folder,filesep,d(x).name)},1:length(d)));
        
        % throw error if there are no files in the raw directory
        if isempty(d_files) 
            error(['Error: experiment folder *' ops.experiment, '* or subject folder *' ops.subject '* not found']);
        end

        % clean up list of files
        % (i.e., excludes bad/aborted/crashed files)
        temp_d_files = [];
        for i=1:length(d_files)
            current_d_file = d_files(i);
            % specific to AMC048 in MIT_SWJN_audio_visual, run 5 and AMC099 in MITNLengthSentences, run 4. AMC082 MITLangloc run1.
            if ~contains(current_d_file, "aborted")
                if ~contains(current_d_file,"64_blocks_switched")
                    % specific to run 7 AMC058 MIT_SWJN_audio_visual_task
                    if ~contains(current_d_file, "crashed")
                        temp_d_files = [temp_d_files;current_d_file];
                    end
                end
            end   
        end
        d_files = temp_d_files;

        % construct signal data array
        signal = [];
        curr_stitch = 1;
        stitch_index = []; % first sample in each data file
        for file_idx = 1:length(d_files)
            [signal_,~,parameters_] = load_bcidat(d_files{file_idx});
            signal = [signal, signal_'];

            stitch_index = [stitch_index; curr_stitch];
            curr_stitch = curr_stitch + size(signal_,1);
        end
        
        % load subject operation info
        op_info_file = [OP_INFO_PATH ops.subject '_op_info.mat'];
        load(op_info_file); % e.g., AMC001_op_info
        op_info_var_name = [ops.subject '_op_info'];
        % only need electrode labels and channel types
        ch_label = eval([op_info_var_name '.channel_label;']);
        ch_type  = eval([op_info_var_name '.channel_type;']);
        ch_nums  = 1:length(ch_label); 

        % remove empty rows of signal (no label) 
        % TODO - ask if they are really empty
        signal = signal(1:length(ch_label),:);

        % define electrodes to remove from the start
        if ops.deselectElecsPrelim
            prelim_deselect = [eval([op_info_var_name '.GND']),... 
                               eval([op_info_var_name '.REF']),... 
                               eval([op_info_var_name '.bad_channels']),...
                               eval([op_info_var_name '.skull_eeg_channels']),...
                               eval([op_info_var_name '.microphone_channels']),...
                               eval([op_info_var_name '.EMG_channels'])...
            ];

            % some fields can be strings (e.g., 'DigitalInput4')
            if isa(eval([op_info_var_name '.visual_trigger']),'double')
                prelim_deselect = [prelim_deselect, eval([op_info_var_name '.visual_trigger'])];
            end
            if isa(eval([op_info_var_name '.button_trigger']),'double')
                prelim_deselect = [prelim_deselect, eval([op_info_var_name '.button_trigger'])];
            end
            if isa(eval([op_info_var_name '.buzzer_trigger']),'double')
                prelim_deselect = [prelim_deselect, eval([op_info_var_name '.buzzer_trigger'])];
            end
            if isa(eval([op_info_var_name '.audio_trigger']),'double')
                prelim_deselect = [prelim_deselect, eval([op_info_var_name '.audio_trigger'])];
            end

            prelim_deselect = intersect([1:length(ch_label)], prelim_deselect);
        else
            prelim_deselect = [];
        end
        fprintf(1,'Electrodes removed upon preliminary review (i.e., ground/reference electrodes, non-cortical channels, other user-defined noisy channels): ')
        fprintf(1,'%d ', prelim_deselect(:)); fprintf('\n');

        % sample frequency
        sample_freq = parameters_.SamplingRate.NumericValue;

        % define number of electrodes per amp
        amplifiers = {'gUSBampADC',...
                      'gHIampADC',...
                      'NatusADC',...
                      'NihonKohdenADC',...
        };
        electrodes = {16, 64, 64, 64};
        amp_to_numElecs = containers.Map(amplifiers,electrodes);
        elecs_per_amp = amp_to_numElecs(parameters_.SignalSourceFilterChain.Value{2,1});

        % construct structure to be used for preprocessing
        for_preproc = struct;
        for_preproc.elec_data_raw     = signal;
        for_preproc.log_file_name     = log_filename;
        for_preproc.stitch_index_raw  = stitch_index;
        for_preproc.stitch_index_dec  = ceil(stitch_index / ops.decimation_factor);
        for_preproc.sample_freq_raw   = sample_freq;
        for_preproc.decimation_freq   = sample_freq / ops.decimation_factor;
        for_preproc.decimation_factor = ops.decimation_factor;
        for_preproc.elecs_per_amp     = elecs_per_amp;


        % ---------------------------------
        % BUILD OBJECT & PREPROCESS SIGNAL
        % ---------------------------------
        obj = ecog_data(for_preproc,...
                        ops.subject,...
                        ops.experiment,...         
                        save_filename,...
                        SAVE_PATH,...
                        d_files,...
                        RAW_PATH,...   
                        ch_label,...
                        ch_nums',...
                        prelim_deselect',...
                        ch_type...
        );

        % mark channels to remove during visual inspection if already done
        if ops.doneVisualInspection
            filename = 'visual_inspection_working.csv'; 
            visual_inspection = readtable(filename,'Delimiter',',','NumHeaderLines',0);
            sub_idx = find(strcmp(visual_inspection.subject,obj.subject));
            user_deselect = visual_inspection{sub_idx,2:end};
            user_deselect = user_deselect(~isnan(user_deselect)) % trim array
            obj.subject
            if ~isempty(user_deselect) % some subjects won't have any channels removed
                obj.elec_ch_user_deselect = user_deselect';
                obj.define_clean_channels();
            end
        end
                    
        obj.preprocess_signal('order',ops.order,'isPlotVisible',ops.isPlotVisible,'doneVisualInspection',ops.doneVisualInspection);

        % add anatomy
        ANATOMY_PATH = [MINDHIVE_PATH filesep 'DATA' filesep '_VERA_ANT' filesep];
        if ~strcmp(obj.subject,'BJH014') & ~strcmp(obj.subject,'BJH015')
            obj.add_anatomy(ANATOMY_PATH);
        end


        % --------------------------
        % FILL IN TRIAL INFORMATION
        % --------------------------

        if strcmp(ops.experiment,'MITNaturalisticStoriesTask') % special instructions for filling in
            % data to fill in
            trial_timing_raw = cell(length(d_files),1); % StimType, caption (stimulus), onset, offset
            trial_timing_dec  = cell(length(d_files),1);
            condition        = cell(length(d_files),1); % condition name (and sometimes number)
            session          = cell(length(d_files),1); % file/session number/name
            events_table     = cell(length(d_files),1); % consolidated table including probe and RT info

            trial_timing_keys = {'key',...
                                'string',...
                                'start',...
                                'end'...
            };
            story_names = {'Boar',...
                        'Aqua',...
                        'MatchstickSeller',...
                        'KingOfBirds',...
                        'Elvis',...
                        'MrSticky',...
                        'HighSchool',...
                        'Roswell',...
                        'Tulips',...
                        'Tourettes',...
            };

            % load stimulus timing file 
            STIM_PATH = [MINDHIVE_PATH filesep 'merged_ecog_pipeline' filesep 'expts' filesep ops.experiment filesep];
            stim_filename = [STIM_PATH 'naturalstories.t.itemmeasures.csv'];
            
            stim_table = readtable(stim_filename);
            stim_table.onsetsampleraw    = round(stim_table.onsettime*obj.for_preproc.sample_freq_raw);
            stim_table.offsetsampleraw   = round(stim_table.offsettime*obj.for_preproc.sample_freq_raw);
            stim_table.midpointsampleraw = round(stim_table.midpointtime*obj.for_preproc.sample_freq_raw);
            stim_table.pausedurationraw  = round(stim_table.pausetime*obj.for_preproc.sample_freq_raw);
            stim_table.onsetsampledec    = round(stim_table.onsettime*obj.for_preproc.decimation_freq);
            stim_table.offsetsampledec   = round(stim_table.offsettime*obj.for_preproc.decimation_freq);
            stim_table.midpointsampledec = round(stim_table.midpointtime*obj.for_preproc.decimation_freq);
            stim_table.pausedurationdec  = round(stim_table.pausetime*obj.for_preproc.decimation_freq);

            % go through each file again
            for file_idx = 1:length(d_files)
                fprintf(1,'\n> Processing trial info from file %d of %d ... \n',file_idx,length(d_files));
                fprintf(1,'[')


                % --- LOAD PARAMETERS AND STATES ---
                [~,states_,parameters_] = load_bcidat(d_files{file_idx});


                % --- EXTRACT TRIAL INFO FROM STIMULI TABLE --
                % TODO - description
                %
                % Returns:
                % -------- 
                %   trial_timing_raw_ :
                %   trial_timing_dec_ : 
                %   condition_        :
                %   session_          :
                %   events_table_     :

                % outputs
                trial_timing_raw_ = cell(1,1); 
                trial_timing_dec_ = cell(1,1);
                condition_        = cell(1,1); 
                session_          = zeros(1,1);
                events_table_     = []; % will add rows incrementally
                trial_num         = 1;

                % --- stimuli ----
                % stimulus string / story presented (e.g., Elvis)
                story_num = states_.StimulusCode(find(diff(states_.StimulusCode)>0)+1);
                story_name = story_names(story_num);
                    
                % stimulus type (e.g., story_1, story_2)
                stimulus_type = {['story_' num2str(story_num)]};


                % --- timing ---
                % initialize stimuli range arrays and downsample StimulusCode
                stimuli_range_raw = [];
                stimuli_range_dec = [];
                decimation_factor = obj.for_preproc.decimation_factor;
                stimCode_raw = states_.StimulusCode;
                stimCode_dec = downsample(double(states_.StimulusCode),decimation_factor);

                % raw sample freq
                stimulus_index_raw = find(stimCode_raw); % SAMPLES RAW
                assert(isempty(find(diff(stimulus_index_raw)>1)),"Story stopped and restarted in stimulus code");
                stimuli_range_raw = [min(stimulus_index_raw), max(stimulus_index_raw)]; 
                        
                % downsampled freq
                stimulus_index_dec = find(stimCode_dec); % SAMPLES DS
                stimuli_range_dec = [stimuli_range_dec; [min(stimulus_index_dec), max(stimulus_index_dec)]]; 

                stimuli_range_raw = num2cell(stimuli_range_raw);
                stimuli_range_dec = num2cell(stimuli_range_dec);

                % construct trial timing tables 
                trial_timing_cell_raw = [stimulus_type,story_name,stimuli_range_raw]; % raw
                trial_timing_raw_(trial_num) = {cell2table(trial_timing_cell_raw,'VariableNames',trial_timing_keys)};
                trial_timing_cell_dec = [stimulus_type,story_name,stimuli_range_dec]; % ds
                trial_timing_dec_(trial_num) = {cell2table(trial_timing_cell_dec,'VariableNames',trial_timing_keys)};


                % --- condition ---
                % story name
                condition_(trial_num) = story_name;
                    

                % --- session ---
                % session number 
                session_(trial_num) = file_idx;


                % --- events table ---
                % NOTE : different than events table for other experiments
                % add rows from story to events table
                events_table_ = stim_table(strcmp(stim_table.docid,story_name),:);


                fprintf(1,'.');

                trial_timing_raw(file_idx) = {trial_timing_raw_};
                trial_timing_dec(file_idx) = {trial_timing_dec_};
                condition(file_idx)        = {condition_};
                session(file_idx)          = {session_};
                events_table(file_idx)     = {events_table_};

                fprintf(1,'] done\n');

            end


        else % all other experiments (e.g., MITLangloc, MITConstituentBounds) 
            % data to fill in
            trial_timing_raw = cell(length(d_files),1); % StimType, caption (stimulus), onset, offset
            trial_timing_dec  = cell(length(d_files),1);
            condition        = cell(length(d_files),1); % condition name (and sometimes number)
            session          = cell(length(d_files),1); % file/session number/name
            events_table     = cell(length(d_files),1); % consolidated table including probe and RT info

            trial_timing_keys = {'key',...
                                'string',...
                                'start',...
                                'end'...
            };
            events_table_keys = {'session',...
                                'trial',...
                                'trial_onset',...
                                'condition',...
                                'stimulus_string',...
                                'stimulus_offset',...
                                'probe',...
                                'probe_result',...
                                'response',...
                                'RT',...
                                'accuracy'...
            };

            % go through each file again
            for file_idx = 1:length(d_files)
                fprintf(1,'\n> Processing trial info from file %d of %d ... \n',file_idx,length(d_files));
                fprintf(1,'[')


                % --- LOAD PARAMETERS AND STATES ---
                [~,states_,parameters_] = load_bcidat(d_files{file_idx});
            
                % extract stimulus and sequence info from file 
                stimuli_sequence = parameters_.Sequence.NumericValue;
                trials_value     = parameters_.Stimuli.NumericValue;
                stimuli_value    = parameters_.Stimuli.Value;

                % pull out row in Stimuli with name 'TrialNumber'
                TrialNumber_indx = cell2mat(cellfun(@(x) strcmp(x,'TrialNumber'),parameters_.Stimuli.RowLabels,'UniformOutput',false)); % logical array
                trial_for_stimuli_seq = trials_value(TrialNumber_indx,:);

                % list of trial numbers without NAs
                trials = rmmissing(unique(trial_for_stimuli_seq));
                

                % --- EXTRACT TRIAL SEQUENCE VALUES ---
                % Sequence is the order of the events/trials, NOT onsets/offsets in samples)
                % 
                % Returns:
                % --------
                % trial_seq_cell: nTrials x 2 cell array 
                %   Sequence values that correspond to trial (e.g., 3:17), trial number)

                trial_seq_cell = {};
                for trial_num = 1:length(trials)

                    % indices of trial in Stimuli cell array (and Sequence)
                    trial_stimuli_sequence = find(trial_for_stimuli_seq == trials(trial_num));

                    % starting index of trial in sequence
                    trial_instance_in_sequence = strfind(stimuli_sequence',trial_stimuli_sequence); 

                    % only one instance of trial in sequence (expected)
                    if length(trial_instance_in_sequence) == 1 
                        trial_seq_cell{trial_num,1} = stimuli_sequence(trial_instance_in_sequence+[0:length(trial_stimuli_sequence)-1]); % full sequence values of trial
                        trial_seq_cell{trial_num,2} = trials(trial_num); % trial num
                    % no instance of trial in sequence
                    elseif isempty(trial_instance_in_sequence)
                        fprintf('the stimuli for trial %d was not found in the parameter.Sequence \n',trials(trial_num))
                    % more than one instance of trial in sequence
                    else 
                        fprintf('more than one instance of trial found\n');
                    end
                end

                % remove fixations from trial_seq_cell for AMC026
                if strcmp(obj.subject,'AMC026')
                    trial_seq_cell_new = [];
                    for trial_num = 1:length(trial_seq_cell)
                        StimType_indx=cell2mat(cellfun(@(x) strcmp(x,'StimType'),parameters_.Stimuli.RowLabels,'UniformOutput',false));
                        stimulus_type = stimuli_value(find(StimType_indx),trial_for_stimuli_seq == trial_seq_cell{trial_num,2});
                        if ~strcmp(stimulus_type,'fixation')
                            trial_seq_cell_new = [trial_seq_cell_new; trial_seq_cell(trial_num,:)];
                        end
                    end
                    trial_seq_cell = trial_seq_cell_new;
                end


                % --- EXTRACT TRIAL INFO FROM STIMULI TABLE --
                % TODO - description
                %
                % Returns:
                % -------- 
                %   trial_timing_raw_ :
                %   trial_timing_dec_ : 
                %   condition_        :
                %   session_          :
                %   events_table_     :

                % outputs
                trial_timing_raw_ = cell(length(trial_seq_cell),1); 
                trial_timing_dec_ = cell(length(trial_seq_cell),1);
                condition_        = cell(length(trial_seq_cell),1); 
                session_          = zeros(length(trial_seq_cell),1);
                events_table_     = []; % will add rows incrementally

                % rows of interest (logical arrays)
                caption_indx=cell2mat(cellfun(@(x) strcmp(x,'caption'),parameters_.Stimuli.RowLabels,'UniformOutput',false)); 
                StimType_indx=cell2mat(cellfun(@(x) strcmp(x,'StimType'),parameters_.Stimuli.RowLabels,'UniformOutput',false)); 
                ConditionName_indx=cell2mat(cellfun(@(x) strcmp(x,'ConditionName'),parameters_.Stimuli.RowLabels,'UniformOutput',false));
                if strcmp(obj.subject,'AMC026')
                    ConditionName_indx=cell2mat(cellfun(@(x) strcmp(x,'WordType'),parameters_.Stimuli.RowLabels,'UniformOutput',false));
                end
                IsRight_indx=cell2mat(cellfun(@(x) strcmp(x,'IsRight'),parameters_.Stimuli.RowLabels,'UniformOutput',false)) | ...
                    cell2mat(cellfun(@(x) strcmp(x,'IsProbeCorrect'),parameters_.Stimuli.RowLabels,'UniformOutput',false)) | ...
                    cell2mat(cellfun(@(x) strcmp(x,'ComprehensionAnswer'),parameters_.Stimuli.RowLabels,'UniformOutput',false)); % row name depends on expt

                % go through each trial to extract info from parameters.Stimuli
                for trial_num = 1:length(trial_seq_cell)

                    % indices in sequence of trial
                    trial_indx = trial_seq_cell{trial_num}; 
                    

                    % --- stimuli ----
                    % stimuli presented (i.e., caption)
                    caption = stimuli_value(find(caption_indx),trial_for_stimuli_seq == trial_seq_cell{trial_num,2});
                    
                    % stimulus type (e.g., word, preprobe, question)
                    stimulus_type = stimuli_value(find(StimType_indx),trial_for_stimuli_seq == trial_seq_cell{trial_num,2});
                    word_mask = cell2mat(cellfun(@(x) strcmp(x,'word'),stimulus_type,'UniformOutput',false));
                    stimulus_type(word_mask) = cell(arrayfun(@(x) [stimulus_type{x} '_' num2str(x)],find(word_mask),'UniformOutput',false)); % append word number 

                    % stimulus string
                    stimulus_string = {strjoin(caption(word_mask))};


                    % --- timing ---
                    % initialize stimuli range arrays and downsample StimulusCode
                    stimuli_range_raw = [];
                    stimuli_range_dec = [];
                    decimation_factor = obj.for_preproc.sample_freq_raw / obj.for_preproc.decimation_freq;
                    stimCode_raw = states_.StimulusCode;
                    stimCode_dec = downsample(double(states_.StimulusCode),decimation_factor);

                    % go through each event in trial
                    for event_num = 1:length(trial_indx)

                        % raw sample freq
                        stimulus_index_raw = find(stimCode_raw == trial_indx(event_num)); % SAMPLES RAW
                        stimuli_range_raw = [stimuli_range_raw; [min(stimulus_index_raw), max(stimulus_index_raw)]]; 
                        
                        % downsampled freq
                        stimulus_index_dec = find(stimCode_dec == trial_indx(event_num)); % SAMPLES DS
                        stimuli_range_dec = [stimuli_range_dec; [min(stimulus_index_dec), max(stimulus_index_dec)]]; 

                    end
                    stimuli_range_raw = num2cell(stimuli_range_raw);
                    stimuli_range_dec = num2cell(stimuli_range_dec);

                    % construct trial timing tables 
                    trial_timing_cell_raw = [stimulus_type',caption',stimuli_range_raw]; % raw
                    trial_timing_raw_(trial_num) = {cell2table(trial_timing_cell_raw,'VariableNames',trial_timing_keys)};
                    trial_timing_cell_dec = [stimulus_type',caption',stimuli_range_dec]; % ds
                    trial_timing_dec_(trial_num) = {cell2table(trial_timing_cell_dec,'VariableNames',trial_timing_keys)};


                    % --- condition ---
                    % one condition name for each row in Stimuli with condition info (e.g., Condition, ConditionName)
                    condition_name = stimuli_value(find(ConditionName_indx),trial_for_stimuli_seq == trial_seq_cell{trial_num,2});
                    condition_name = unique(condition_name); % (e.g., 'threebound_47end'), 
                    if length(condition_name)>1
                        condition_name = condition_name(cellfun(@(x) ~isempty(x),condition_name));
                    end
                    assert(length(condition_name)==1,'more than one, or no, condition label for trial %d \n',trial_num);
                    condition_(trial_num) = condition_name;


                    % --- session ---
                    % session number (not from Stimuli table)
                    session_(trial_num) = file_idx;


                    % --- probe/question ---
                    % the probe/question itself
                    probe_mask = cell2mat(cellfun(@(x) strcmp(x,'probe'),stimulus_type,'UniformOutput',false)) | ...
                                cell2mat(cellfun(@(x) strcmp(x,'question'),stimulus_type,'UniformOutput',false));
                    assert(sum(probe_mask)==1,'more than one probe/question');
                    probe = caption(probe_mask); 

                    % whether the probe/question is correct 
                    % (not whether the subject got it right or wrong)
                    assert(sum(IsRight_indx) == 1,'more than one, or no, sequence number with IsRight info for trial %d \n',trial_num);
                    probe_result_full_trial = stimuli_value(find(IsRight_indx),trial_for_stimuli_seq == trial_seq_cell{trial_num,2});
                    probe_result_mask = cell2mat(cellfun(@(x) ~isempty(x),probe_result_full_trial,'UniformOutput',false)); % logical array
                    assert(sum(probe_result_mask)==1,'more than one probe or question result for trial %d \n',trial_num);
                    probe_result = probe_result_full_trial(probe_result_mask);

                    % change probe result to Y/N if 1/0 (e.g., MITLangloc)
                    % OR to Y/N if RIGHT/WRONG (e.g., MITSWJNTask AMC026)
                    if strcmp(probe_result{1},'0') || strcmp(probe_result{1},'WRONG')
                        probe_result = {'N'};
                    elseif strcmp(probe_result{1},'1') || strcmp(probe_result{1},'RIGHT')
                        probe_result = {'Y'};
                    end


                    % --- response ---
                    % subject response search window
                    index_isright_start = trial_timing_cell_raw{probe_result_mask,3}; % onset
                    index_isright_stop = trial_timing_cell_raw{probe_result_mask,4}; % offset 
                    buffer_before = obj.for_preproc.sample_freq_raw * 1; % 1 sec 
                    buffer_after  = obj.for_preproc.sample_freq_raw * 2; % 2 sec
                    start_index = max(index_isright_start-buffer_before,1);
                    end_index = min(index_isright_stop+buffer_after,length(states_.KeyDown));

                    % look for subject response
                    KeyDown = unique(states_.KeyDown(start_index:end_index));
                    KeyDown = intersect(KeyDown,[67,77]); % ,99,109
                    if length(KeyDown) ~= 1                % too many key's pressed or incorrect response
                        responseTime = nan;
                        TrialResponse = 'TOO_MANY_KEYS';
                    elseif KeyDown == 67 % || KeyDown == 99  % response is yes (1)
                        KeyDown_idxs = find(states_.KeyDown(start_index:end_index)==KeyDown);
                        responseTime = (KeyDown_idxs(1)-(obj.for_preproc.sample_freq_raw+1)) / obj.for_preproc.sample_freq_raw;
                        TrialResponse = 'Y'; 
                    elseif KeyDown == 77 % || KeyDown == 109 % response is no  (2) 
                        KeyDown_idxs = find(states_.KeyDown(start_index:end_index)==KeyDown);
                        responseTime = (KeyDown_idxs(1)-(obj.for_preproc.sample_freq_raw+1)) / obj.for_preproc.sample_freq_raw;
                        TrialResponse = 'N'; 
                    else                                   % incorrect response 
                        responseTime = nan;
                        TrialResponse = 'INCORRECT_KEY';
                    end

                    % deterimine if subject was right or wrong
                    if strcmp(probe_result{1},TrialResponse) 
                        accuracy = 1; % correct
                    elseif (strcmp(probe_result{1},'N') && (strcmp(TrialResponse,'Y') || strcmp(TrialResponse,'INCORRECT_KEY') || strcmp(TrialResponse,'TOO_MANY_KEYS'))) || ...
                        (strcmp(probe_result{1},'Y') && (strcmp(TrialResponse,'N') || strcmp(TrialResponse,'INCORRECT_KEY') || strcmp(TrialResponse,'TOO_MANY_KEYS')))
                        accuracy = 0; % incorrect 
                    else
                        error('probe result is %s, not Y or N\n', probe_result{1});
                    end

                    
                    % --- events table ---
                    % onsets/offsets for events table
                    trial_onset = trial_timing_cell_raw{1,3} / obj.for_preproc.sample_freq_raw; % start of trial
                    word_trial_timing = trial_timing_cell_raw(word_mask,:);
                    word_offset = word_trial_timing{end,end} / obj.for_preproc.sample_freq_raw; % end of last word before probe

                    % add row to events table
                    new_row = [{file_idx},...
                            {trial_num},...
                            {trial_onset},...
                            condition_name,...
                            stimulus_string,...
                            {word_offset},...
                            probe,...
                            probe_result,...
                            {TrialResponse},...
                            {responseTime},...
                            {accuracy}...
                    ]; % cell array

                    events_table_ = [events_table_; cell2table(new_row,'VariableNames',events_table_keys)];

                    fprintf(1,'.');

                end

                trial_timing_raw(file_idx) = {trial_timing_raw_};
                trial_timing_dec(file_idx)  = {trial_timing_dec_};
                condition(file_idx)        = {condition_};
                session(file_idx)          = {session_};
                events_table(file_idx)     = {events_table_};

                fprintf(1,'] done\n');

            end
        end

        fprintf(1,'\nDONE PROCESSING TRIAL INFO \n');


        % --- ADD TRIAL INFO TO OBJECT ---
        % TODO - description

        if length(d_files) == 1 % only one file with signal

            % no need for 1x1 cell arrays
            obj.condition    = condition_;
            obj.session      = session_;
            obj.events_table = events_table_;

            % add both versions of trial timing (raw & ds) 
            obj.for_preproc.trial_timing_raw = trial_timing_raw_;
            obj.for_preproc.trial_timing_dec = trial_timing_dec_;

            % update obj.trial_timing (trial timing for current version of signal)
            if obj.sample_freq == obj.for_preproc.decimation_freq % downsampling in preproc
                obj.trial_timing = trial_timing_dec_;
            else % no downsampling in preproc
                obj.trial_timing = trial_timing_raw_;
            end

        else % multiple files with signal

            % NOTE - cells will need to be combined
            obj.condition    = condition;
            obj.session      = session;
            obj.events_table = events_table;
            
            % add both versions of trial timing (raw & ds) 
            obj.for_preproc.trial_timing_raw = trial_timing_raw;
            obj.for_preproc.trial_timing_dec = trial_timing_dec;

            % update obj.trial_timing (trial timing for current version of signal)
            if obj.sample_freq == obj.for_preproc.decimation_freq % downsampling in preproc
                obj.trial_timing = trial_timing_dec;
            else % no downsampling in preproc
                obj.trial_timing = trial_timing_raw;
            end

            % combine separate data files
            obj.combine_data_files();
            
        end
        

        % --- SAVE OUTPUT ---
        fprintf(1,'\n> Saving object ... \n');
        diary off

        % make save directory
        if ~exist(SAVE_PATH, 'dir')
            mkdir(SAVE_PATH);
        end
    
        save(save_filename,'obj','-v7.3');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPROCESS FROM EXISTING OBJECT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        
        % add ecog_sn_data object path to matlab path
        % ONLY NECESSARY IF LANGLOC HAS ALREADY BEEN APPLIED
        ANALYSIS_PATH = [MINDHIVE_PATH filesep 'ANALYSIS' filesep 'MITLangloc'];
        addpath(genpath(ANALYSIS_PATH)); 

        % load object
        if isfile(save_filename)
            fprintf(1,'\n> Loading existing object file ...\n');
            load(save_filename); % obj
        else
            error('Object file does not exist: \n%s',save_filename);
        end


        % -------------------------------------
        % CREATE SIGNAL PRE ENVELOPE EXTRACTION
        % -------------------------------------
        if ops.preEnvelopeExtraction

            obj.preprocess_signal('order',ops.order);

            % build separate structure to save
            pre_extract = struct;
            pre_extract.elec_data     = obj.elec_data;
            pre_extract.bip_elec_data = obj.bip_elec_data;
            pre_extract.stitch_index  = obj.stitch_index; 
            pre_extract.sample_freq   = obj.sample_freq;
            pre_extract.for_preproc   = obj.for_preproc;
            pre_extract.subject       = obj.subject;
            pre_extract.experiment    = obj.experiment;

            % clear raw signal
            pre_extract.for_preproc.elec_data_raw = [];

            % overwrite obj so save names are the same 
            % NOTE - don't do this if you want to save the entire object
            obj = pre_extract;

            % save
            PRE_EXTRACT_SAVE_PATH = [SAVE_PATH 'structures' filesep 'pre_envelope_extraction' filesep];
            pre_extract_save_filename = [PRE_EXTRACT_SAVE_PATH ops.subject '_' ops.experiment '_crunched_pre_envelope_extraction.mat'];
            save(pre_extract_save_filename,'obj','-v7.3');

        end


        % -----------
        % ADD ANATOMY
        % -----------
        if ops.addAnatomy

            ANATOMY_PATH = [MINDHIVE_PATH filesep 'DATA' filesep '_VERA_ANT' filesep];
            obj.add_anatomy(ANATOMY_PATH);

            fprintf(1,'\n> Saving object ... \n');
            save(save_filename,'obj','-v7.3');

        end

    end

end