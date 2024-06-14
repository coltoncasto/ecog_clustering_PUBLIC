function analyze_MITLangloc(varargin)

    p = inputParser();
    % type of analysis
    addParameter(p,'fromScratch',false);
    addParameter(p,'plotLangResponsive',false);
    addParameter(p,'plotLangResponsiveBipolarByShank',false);
    addParameter(p,'summaryStatistics',false);
    % subjects on which analysis is performed
    addParameter(p,'doOneSub',[]);
    addParameter(p,'doOneSubGroup',[]);
    % output
    addParameter(p,'toStruct',false);
    addParameter(p,'formatCDR',false);
    addParameter(p,'renameAMC026Conditions',false);
    % other parameters
    addParameter(p,'experiment','MITLangloc');
    addParameter(p,'suffix','crunched*');
    parse(p, varargin{:});
    ops = p.Results;


    % --- INITIALIZE ---
    MINDHIVE_PATH = '/mindhive/evlab/u/Shared/ECoG';
    PREPROC_CLASS_PATH = [MINDHIVE_PATH filesep 'merged_ecog_pipeline' filesep];
    CRUNCHED_PATH = [MINDHIVE_PATH filesep 'crunched' filesep ops.experiment filesep 'crunched_for_merged_pipeline' filesep];
    ANALYSIS_PATH = [MINDHIVE_PATH filesep 'ANALYSIS' filesep 'MITLangloc'];
    SAVE_PATH = [ANALYSIS_PATH filesep 'output' filesep];

    addpath(genpath(PREPROC_CLASS_PATH)); 
    addpath(genpath([ANALYSIS_PATH filesep 'utils'])); 

    % crunched files
    if ops.doOneSub
        d = dir([CRUNCHED_PATH ops.doOneSub '*' ops.suffix '.mat']);
        assert(length(d)==1,'More than one crunched file for subject %s, experiment %s',ops.doOneSub,ops.experiment)
    elseif ops.doOneSubGroup
        d = dir([CRUNCHED_PATH ops.doOneSubGroup '*' ops.suffix '.mat']);
    else
        d = dir([CRUNCHED_PATH '*' ops.suffix '.mat']);
    end
    d_files = transpose(arrayfun(@(x) {strcat(d(x).folder,filesep,d(x).name)},1:length(d)));

    % summary statistics table
    if ops.summaryStatistics
        all_summary_stats = [];
    end

    % define parameters by experiment
    if strcmp(ops.experiment,'MITSWJNTask')
        ops.words = [1:8];
        ops.SConditionName = 'SENTENCES';
        ops.NConditionName = 'NONWORDS';
        ops.WConditionName = 'WORDS';
        ops.JConditionName = 'JABBERWOCKY';
    elseif strcmp(ops.experiment,'MITLangloc')
        ops.words = [1:12];
        ops.SConditionName = 'Sentences';
        ops.NConditionName = 'Jabberwocky'; % actually N but labeled as J in MITLangloc
    end

    % runs to keep by subject
    subjects = {'AMC026',...
                'AMC029',...
                'AMC031',...
                'AMC037',...
                'AMC038',...
                'AMC044',...
                'AMC082',...
                'AMC086',...
                'AMC088',...
                'AMC091',...
                'AMC092',...
                'AMC096',...
                'AMC097',...
                'AMC099',...
                'MCJ011',...
                'MCJ014',...
                'BJH006',...
                'BJH007',...
                'BJH008',...
                'BJH011',...
                'BJH012',...
                'BJH016',...
                'SLCH002'...
    };
    sessions = {[];...
                [];...
                [];...
                [];...
                [];...
                [];...
                [2:4];...
                [3:5];...
                [1:2];...
                [2:3];...
                [2:3];...
                [2:3];...
                [2:3];...
                [2:3];...
                [1];...
                [1:2];...
                [2:3];...
                [1:2];...
                [1:2];...
                [1:2];...
                [1:2];...
                [1:2];...
                [1:2]...
    };
    sub_sess_map = containers.Map(subjects,sessions);


    % --- ANALYSIS ---
    % go through subjects
    for i=1:length(d_files)
        subject = split(d_files{i},'/'); subject = subject{end}(1:6); % 6 = length('AMC099')
        fprintf(1,'\nAnalyzing %s, subject %d of %d ...\n',subject,i,length(d_files));
        

        % --- S VS N ---
        if ops.fromScratch

            % construct object 
            obj = ecog_sn_data(SAVE_PATH,...
                               d_files{i},...
                               CRUNCHED_PATH,...
                               [PREPROC_CLASS_PATH 'ecog_data.m'],...
                               PREPROC_CLASS_PATH...
            )
        
            % cut into trials if it hasn't already been done
            if isempty(obj.trial_data)
                obj.make_trials();
            end

            % test S vs. N
            obj.test_s_vs_n('words',ops.words,...
                            'S_condition_flag',ops.SConditionName,...
                            'N_condition_flag',ops.NConditionName,...
                            'n_rep',10000,...
                            'threshold',0.05,...
                            'side','right',...
                            'sessions',sub_sess_map(obj.subject),...
                            'do_plot',false...
            );

            % clear trial data (saves space)
            if ~isempty(obj.trial_data)
                obj.trial_data = [];
            end

            % save object again
            fprintf(1,'\n> Saving object wth language responsive electrodes ... \n');
            save(d_files{i},'obj','-v7.3');

        end

        % --- PLOT LANG RESPONSIVE CHANNELS ---
        if ops.plotLangResponsive

            % load object if not already in workspace
            if ~ops.fromScratch
                fprintf(1,'\n> Loading already constructed SN object ...\n');
                load(d_files{i})
            end

            % plot lang resp channels
            fprintf(1,'\n> Plotting language-responsive channels ... \n');
            true_false = {false,true}; % BE CAREFUL!! 'true' will z-score the signal in the object
            if strcmp(ops.experiment,'MITLangloc')
                for j=1:length(true_false)
                    obj.lang_resp_plots('words',ops.words,...
                                        'S_condition_flag',ops.SConditionName,...
                                        'N_condition_flag',ops.NConditionName,...
                                        'sessions',sub_sess_map(obj.subject),...
                                        'subAverage',true_false{j}...
                    );
                end

            elseif strcmp(ops.experiment,'MITSWJNTask')
                for j=1:length(true_false)
                    obj.lang_resp_plots('words',ops.words,...
                                        'S_condition_flag',ops.SConditionName,...
                                        'W_condition_flag',ops.WConditionName,...
                                        'N_condition_flag',ops.NConditionName,...
                                        'J_condition_flag',ops.JConditionName,...
                                        'subAverage',true_false{j}...
                    );
                end
            end
            
        end

        % --- PLOT LANG RESPONSIVE BIPOLAR CHANNELS BY SHANK ---
        if ops.plotLangResponsiveBipolarByShank

            % load object if not already in workspace
            if ~ops.fromScratch
                fprintf(1,'\n> Loading already constructed SN object ...\n');
                load(d_files{i})
            end

            % plot lang resp channels
            if ~isempty(obj.bip_elec_data)
                fprintf(1,'\n> Plotting bipolar language-responsive channels by shank ... \n');
                obj.lang_resp_plots('words',ops.words,...
                                    'S_condition_flag',ops.SConditionName,...
                                    'N_condition_flag',ops.NConditionName,...
                                    'sessions',sub_sess_map(obj.subject),...
                                    'bipolarByShank',true...
                );
            end
      
        end


        % --- SUMMARY STATISTICS ---
        if ops.summaryStatistics
        
            % load object if not already in workspace
            if ~ops.fromScratch && ~ops.plotLangResponsive
                fprintf(1,'\n> Loading already constructed SN object ...\n');
                load(d_files{i})
            end

            % get summary statistics for subject
            fprintf(1,'\n> Summarizing language responsive electrodes ... \n')
            sub_stats = obj.get_summary_statistics('sessions',sub_sess_map(obj.subject));

            all_summary_stats = [all_summary_stats; sub_stats];

        end  


        % --- FORMAT FOR CDR ---
        if ops.formatCDR
        
            % load object if not already in workspace
            if ~ops.fromScratch && ~ops.plotLangResponsive
                fprintf(1,'\n> Loading already constructed SN object ...\n');
                load(d_files{i})
            end

            % get summary statistics for subject
            fprintf(1,'\n> Formatting data for CDR ... \n')
            [event, resp] = obj.format_for_CDR();

        end

        
        % --- OBJECTS TO STRUCTS ---
        if ops.toStruct

            % load object
            fprintf(1,'\n> Loading obj ...\n');
            load(d_files{i}); % obj

            % make save directory
            CRUNCHED_PATH_STRUCT = [CRUNCHED_PATH 'structures' filesep];
            if ~exist(CRUNCHED_PATH_STRUCT, 'dir')
                mkdir(CRUNCHED_PATH_STRUCT);
            end

            % make save directory
            obj = struct(obj);
            fprintf(1,'\n> Saving struct...');
            save([CRUNCHED_PATH_STRUCT d(i).name],'obj','-v7.3');
            
        end

        % --- RENAME AMC026 CONDITIONS ---
        if ops.renameAMC026Conditions
        
            % load object if not already in workspace
            if ~ops.fromScratch && ~ops.plotLangResponsive
                fprintf(1,'\n> Loading already constructed SN object ...\n');
                load(d_files{i})
            end

            % get summary statistics for subject
            fprintf(1,'\n> Renaming conditions for subject AMC026 ... \n');
            orig_conds = {'1','2','3','4'};
            new_conds = {'SENTENCES','WORDS','NONWORDS','JABBERWOCKY'}
            condMap = containers.Map(orig_conds,new_conds)
            
            % update condition property of object
            obj.condition = cellfun(@(x) condMap(x),obj.condition,'UniformOutput',false)

            % update event table
            obj.events_table.condition = obj.condition;

            % save object again
            fprintf(1,'\n> Saving object wth language responsive electrodes ... \n');
            save(d_files{i},'obj','-v7.3');

        end 

        fprintf(1,'\nDONE WITH SUBJECT\n');

    end

    
    % --- SUMMARY STATISTICS --
    filename = [SAVE_PATH 'summary_statistics_' ops.experiment '.mat'];
    if ops.summaryStatistics
        save(filename,'all_summary_stats');
    end


    % --- AVERAGE OVER ALL SUBJECTS ----


end