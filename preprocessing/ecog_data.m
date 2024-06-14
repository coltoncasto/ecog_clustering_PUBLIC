classdef ecog_data < dynamicprops
% ECOG_DATA 
% TODO - SUMMARY

properties
    %% ---- DATA ----
    elec_data
    bip_elec_data
    stitch_index
    sample_freq
    for_preproc             % preproc
    trial_data              % trial

    %% ---- INFO ----
    subject
    experiment
    trial_timing
    events_table
    condition
    session
    crunched_file_name      % output
    crunched_file_path
    raw_file_name           % input
    raw_file_path

    %% ---- LABELS ----
    elec_ch                 % unipolar
    elec_ch_label 
    elec_ch_prelim_deselect
    elec_ch_with_IED
    elec_ch_with_noise
    elec_ch_user_deselect
    elec_ch_clean
    elec_ch_valid
    elec_ch_type
    bip_ch                  % bipolar
    bip_ch_label
    bip_ch_valid
    bip_ch_grp             
    bip_ch_label_grp
    

    %% ---- ANATOMY ----
    anatomy
    
end


methods
    %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = ecog_data(...
            for_preproc,...             % DATA
            subject,...                 % INFO
            experiment,...
            crunched_file_name,...
            crunched_file_path,...
            raw_file_name,...
            raw_file_path,...
            elec_ch_label,...           % LABELS
            elec_ch,...
            elec_ch_prelim_deselect,...
            elec_ch_type)
        
        % Construct an instance of this class

        %% ---- DATA ----
        obj.for_preproc=for_preproc;

        %% ---- INFO ----
        obj.subject=subject;
        obj.experiment=experiment;
        obj.crunched_file_name=crunched_file_name;
        obj.crunched_file_path=crunched_file_path;
        obj.raw_file_name=raw_file_name;
        obj.raw_file_path=raw_file_path;

        %% ---- LABELS ----
        obj.elec_ch=elec_ch;
        obj.elec_ch_label=elec_ch_label;
        obj.elec_ch_prelim_deselect=elec_ch_prelim_deselect;
        obj.elec_ch_type=elec_ch_type;
        
    end


    %% FUNCTION TO RUN PIPELINE


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPROCESS SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function preprocess_signal(obj,varargin)
        % TODO - description

        p = inputParser();
        addParameter(p,'order','defaultBOTH');
        addParameter(p,'isPlotVisible',true);
        addParameter(p,'doneVisualInspection',false);
        parse(p, varargin{:});
        ops = p.Results;


        % ---------------------------
        % DEFINE PREPROCESSING ORDER
        % ---------------------------
        % this pipeline was constructed for the 'default' ordering
        % all new preprocessing orders should be checked for problems
        if strcmp(ops.order,'defaultECOG')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'CAR',... 
                     'GaussianFilterExtraction',...
                     'removeOutliers',...
                     'downsample'...
            };

        elseif strcmp(ops.order,'defaultSEEGorBOTH')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'CAR',... 
                     'BipolarReferencing'...
                     'GaussianFilterExtraction',...
                     'removeOutliers',...
                     'downsample'...       
            };

        elseif strcmp(ops.order,'SEEGorBOTHbyShank')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'CAR',...
                     'ShankCSR',... 
                     'BipolarReferencing'...
                     'GaussianFilterExtraction',...
                     'removeOutliers',...
                     'downsample'...       
            };

        elseif strcmp(ops.order,'defaultMCJandBJH') 
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'GlobalMeanRemoval',... % could z-score here
                     'IEDRemoval',...
                     'visualInspection',...
                     'BipolarReferencing'...
                     'GaussianFilterExtraction',...
                     'removeOutliers',...
                     'downsample'...
            };

        elseif strcmp(ops.order,'preEnvelopeExtractionECOG')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'CAR',...
                     'downsample'...
            };

        elseif strcmp(ops.order,'preEnvelopeExtractionSEEGorBOTH')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'GlobalMeanRemoval',...
                     'BipolarReferencing'...
                     'downsample'...
            };
    
        elseif strcmp(ops.order,'preEnvelopeExtractionMCJandBJH') 
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'GlobalMeanRemoval',... % could z-score here
                     'IEDRemoval',...
                     'visualInspection',...
                     'BipolarReferencing'...      
                     'downsample'...
            };

        elseif strcmp(ops.order,'visuallyInspectFromRaw')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'GlobalMeanRemoval',...
                     'visualInspection',...
            };

        elseif strcmp(ops.order,'oldPipeline') % as close as we can get
            order = {'visualInspection',...
                     'highpassFilter',... 
                     'CAR',...
                     'notchFilter',...
                     'GaussianFilterExtraction',...
                     'downsample'...
            };

        elseif strcmp(ops.order,'test')
            order = {'downsample'};

        else
            error('Preprocessing order not specified');
        end

        obj.for_preproc.order = order;
        obj.for_preproc.isPlotVisible = ops.isPlotVisible; 

        % mappping from preprocessing name to its function in ecog_data
        names = {'highpassFilter',...
                 'notchFilter',...
                 'IEDRemoval',...
                 'visualInspection',...
                 'GlobalMeanRemoval',...
                 'CAR',...
                 'ShankCSR',...
                 'BipolarReferencing',...
                 'GaussianFilterExtraction',...
                 'BandpassExtraction',...
                 'zscore',...
                 'downsample',...
                 'removeOutliers'...
        };
        functions = {'highpass_filter',...
                     'notch_filter',...
                     'remove_IED',...
                     'visual_inspection',...
                     'reference_signal',...
                     'reference_signal',...
                     'reference_signal',...
                     'reference_signal',...
                     'extract_high_gamma',...
                     'extract_high_gamma',...
                     'zscore_signal',...
                     'downsample_signal',...
                     'remove_outliers'...
        };

        name_to_function = containers.Map(names,functions);
        function_order = cellfun(@(x) name_to_function(x),order,'UniformOutput',false);

        
        % -----------------------------
        % CALL PREPROCESSING FUNCTIONS
        % -----------------------------
        fprintf(1,'\nSTARTING TO PREPROCESS SIGNAL\n');

        obj.first_step('doneVisualInspection',ops.doneVisualInspection);

        prev_step = '';
        for i=1:length(function_order)
            step = function_order{i};

            % skip current step if same as last step 
            if strcmp(step,prev_step)
                continue
            end

            % add flags to reference_signal function
            if strcmp(step,'reference_signal')
                j = i; flags = ""; 
                while  1
                    flags = strcat(flags,"'do",order{j},"',true,");
                    j = j+1;
                    try
                        next_step = function_order{j};
                    catch
                        break;
                    end
                    if ~strcmp(step,next_step)
                        break;
                    end
                end
                flags = char(flags); flags = flags(1:end-1); % remove comma at end
                eval(strcat("obj.",step,"(",flags,");"));

            % add flags to extract_high_gamma function
            elseif strcmp(step,'extract_high_gamma')
                flags = strcat("'do",order{i},"',true");
                eval(strcat("obj.",step,"(",flags,");"));

            % visual inspection
            elseif strcmp(step,'visual_inspection')
                eval(strcat("obj.",step,"('doneVisualInspection',",num2str(ops.doneVisualInspection),");"));
                
            % functions without input flags/arguments
            else
                eval(strcat("obj.",step,"();"));
            end

            prev_step = step;

        end

        fprintf(1,'\nDONE PREPROCESSING SIGNAL \n');

    end


    %% MAIN SIGNAL PREPROCESSING METHODS


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HIGHPASS FILTER SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function highpass_filter(obj)
        % Highpass filters signal
        %
        %

        signal = obj.elec_data';

        highpass = obj.for_preproc.highpass;

        for k=1:length(obj.stitch_index) % number of separate data files with signal
            fprintf(1, '\n> Highpass filtering signal from file %d of %d ... \n',k,length(obj.stitch_index));
            fprintf(1,'[');

            if k == length(obj.stitch_index) % signal for file stops at end of matrix
                stop = size(signal,1);
            else % signal for file stops before stitch index of next file
                stop = obj.stitch_index(k+1)-1;
            end

            signal_ = signal(obj.stitch_index(k):stop,:);

            % highpass filter signal
            parfor idx_channel=1:size(signal_,2)
                warning('off', 'signal:filtfilt:ParseSOS');
                signal_(:,idx_channel) = filtfilt(highpass.sos,highpass.g,double(signal_(:,idx_channel)));
                fprintf(1,'.');
            end

            signal(obj.stitch_index(k):stop,:) = signal_;

            fprintf(1,'] done\n');
        end

        obj.elec_data = signal';

    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTCH FILTER SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function notch_filter(obj)
        % Notch filters signal and marks channels with significant noise

        signal = obj.elec_data';

        notch = obj.for_preproc.notch;

        % measure line noise prior to notch filtering
        signal_noise_before = obj.measure_line_noise(signal);

        for k=1:length(obj.stitch_index) % number of separate data files with signal
            fprintf(1, '\n> Notch filtering signal from file %d of %d ... \n',k,length(obj.stitch_index));
            fprintf(1,'[');

            if k == length(obj.stitch_index) % signal for file stops at end of matrix
                stop = size(signal,1);
            else % signal for file stops before stitch index of next file
                stop = obj.stitch_index(k+1)-1;
            end

            signal_ = signal(obj.stitch_index(k):stop,:);

            % notch filter
            parfor idx_channel=1:size(signal_,2)
                signal_preliminary = double(signal_(:,idx_channel));
                % remove all harmonics of line-noise
                for idx = 1:length(obj.for_preproc.filter_params.notch.fcenter) %#ok<PFBNS>
                    signal_preliminary = filtfilt(notch{idx}.b,notch{idx}.a,signal_preliminary); %#ok<PFBNS>
                end 
                signal_(:,idx_channel) = signal_preliminary;
                fprintf(1,'.');
            end

            signal(obj.stitch_index(k):stop,:) = signal_;

            fprintf(1,'] done\n');
        end

        obj.elec_data = signal';

        % measure line noise after notch filtering
        signal_noise_after = obj.measure_line_noise(signal);

        % mark channels with sig line noise after notch filtering
        obj.elec_ch_with_noise = obj.elec_ch(signal_noise_after(:,2) > (mean(signal_noise_after(:,2))+5*std(signal_noise_after(:,2))));
        obj.elec_ch_with_noise = intersect(obj.elec_ch_clean,obj.elec_ch_with_noise); % don't mark already noisy electrodes
        obj.define_clean_channels();

        obj.for_preproc.notchFilter_results.signal_noise_before_notch = signal_noise_before;
        obj.for_preproc.notchFilter_results.mean_signal_noise_before_notch = mean(signal_noise_before(obj.elec_ch_clean,2));
        obj.for_preproc.notchFilter_results.signal_noise_after_notch = signal_noise_after;
        obj.for_preproc.notchFilter_results.mean_signal_noise_after_notch = mean(signal_noise_after(obj.elec_ch_clean,2));

        fprintf(1,'\nReduced 60 Hz noise from %.2f to %.2f uV\n',mean(signal_noise_before(obj.elec_ch_clean,2)),mean(signal_noise_after(obj.elec_ch_clean,2)));
        fprintf(1,'Electrodes with significant line noise: ');
        fprintf(1,'%d ', obj.elec_ch_with_noise(:)); fprintf('\n');

        % plot line noise to select additional channels to remove
        if obj.for_preproc.isPlotVisible
            f = obj.plot_line_noise(signal_noise_before,signal_noise_after);
            prompt1 = '\nUSER INPUT REQUIRED: \nAdditional channels to remove due to significant line noise? (format: [1,2]) - ';

            more_line_noise = input(prompt1)';
            obj.elec_ch_with_noise = union(obj.elec_ch_with_noise,more_line_noise);
            obj.define_clean_channels();
            close(f);
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIND EPILECTIC ELECTRODES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function remove_IED(obj)
        % This function marks electrodes with significant Interictal 
        % Epileptiform Discharges (IEDs) using a protocol adopted from 
        % Radek Janca (see directory below for sample papers).
        % Once marked, these electrodes will be removed from subsequent 
        % analyses.
        %
        % This is the only function in the pipeline where the parameters
        % are hardcoded into the function itself instead of defined in the 
        % define_parameters function. That is because this script was 
        % designed by another lab and should only be modified with extreme
        % caution. Similarly, this is the only part of the pipeline that 
        % calls on functions not contained within this ecog_data.m file. 
        % This is, again, because the procedure was designed by another 
        % lab and it was determined that their script should be self-
        % contained. 
        %
        % NOTE - the following directory must be added to your MATLAB path
        % /mindhive/evlab/u/Shared/merged_ecog_pipeline/utils/JancaCodePapers 

        signal = double(obj.elec_data');

        fprintf(1, '\n> Finding electrodes with significant Interictal Epileptiform Discharges (IEDs) ... \n');
        
        detectionIEDs          = []; % output from Janca et al. script - 1 cell array per segment
        detectionIEDs.settings = '-k1 3.65 -h 60 -dec 200 -dt 0.005 -pt 0.12 -ti 1'; % if you change "-dec 200" here, do not forget to change in selectChannels_Using ... below
        detectionIEDs.segments = [];
            
        % NOT DEFINED IN THIS CLASS
        detectionIEDs = automaticSpikeDetection_UsingJancaMethod(signal,obj.sample_freq,detectionIEDs.settings);
        
        if obj.for_preproc.isPlotVisible
            obj.plot_channels(detectionIEDs.envelope,...
                            obj.elec_ch_label,...
                            obj.elec_ch_clean,...
                            obj.elec_ch_valid,...
                            't_len',100,...
                            'sample_freq',200,...
                            'plotIEDs',true,...
                            'chanIEDs',detectionIEDs.out.chan,...
                            'posIEDs',detectionIEDs.out.pos...
            ); 
        end
            
        % From this automatic assessment, selection of the final pool of channels
        detectionIEDs.tableChanSelection = [];  % info regarding chan selection
        detectionIEDs.threshold = 6.5; % channels with IEDs higher than threshold are removed  
            
        currIEDs.fs = 200; % default downsampling during automatic detection (-dec 200)
        currIEDs.discharges.MV = [];
        currIEDs.numSamples = 0;
            
        currIEDs.discharges.MV = detectionIEDs.discharges.MV;
        currIEDs.numSamples= currIEDs.numSamples + size(detectionIEDs.d_decim, 1);
            
        % Compute number of detected spike per channel: [c x 1] where c channels
        numSpikes     = []; numSpikes     = sum(currIEDs.discharges.MV==1, 1);
        totalDuration = []; totalDuration = (currIEDs.numSamples / currIEDs.fs) / 60; % in minutes
        numSpikes_min = []; numSpikes_min = numSpikes / totalDuration;
        numSpikes     = transpose(numSpikes);
        numSpikes_min = transpose(numSpikes_min);
            
        % Select channels with IEDs / minute below threshold - [c x 1] where c channels
        indChanSelected = [];
        indChanSelected = find(numSpikes_min < detectionIEDs.threshold);
        tableChanSelection.numSpikesAll           = numSpikes_min;
        tableChanSelection.indChansSelected       = indChanSelected;
        tableChanSelection.indChansDeselected     = setdiff(obj.elec_ch,indChanSelected);
        tableChanSelection.nameChansSelected      = transpose(obj.elec_ch_label(indChanSelected));
        tableChanSelection.numSpikesChansSelected = numSpikes_min(indChanSelected);
            
            
        obj.for_preproc.IEDRemoval_results=tableChanSelection;
        obj.for_preproc.IEDRemoval_results.threshold=detectionIEDs.threshold;
  
        % check if too many electrodes were removed
        if length(tableChanSelection.indChansDeselected) > ceil(size(obj.elec_ch,1)/3)
            fprintf(1,'Too many electrodes with significant IEDs, SKIPPING STEP\n')

            new_order_mask = cell2mat(cellfun(@(x) strcmp(x,'IEDRemoval'),obj.for_preproc.order,'UniformOutput',false));
            obj.for_preproc.order = obj.for_preproc.order(~new_order_mask);

        else 
            obj.elec_ch_with_IED = tableChanSelection.indChansDeselected;
            obj.elec_ch_with_IED = intersect(obj.elec_ch_clean,obj.elec_ch_with_IED); % don't mark already noisy electrodes
            obj.define_clean_channels();
            
            fprintf(1,'Electrodes with significant IEDs: ');
            fprintf(1,'%d ', obj.elec_ch_with_IED(:)); fprintf('\n');

        end
             
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VISUALLY INSPECT SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function visual_inspection(obj,varargin)
        % Visually inspect electrodes without sig line noise
        % User will be asked to provide additional electrodes to remove
        % Should be run after highpass and notch filtering signal

        p = inputParser();
        addParameter(p,'doneVisualInspection',false);
        parse(p, varargin{:});
        ops = p.Results;

        signal = obj.elec_data';

        fprintf(1,'\n> Visually inspecting signal ...\n');

        if obj.for_preproc.isPlotVisible
            obj.plot_channels(signal,...
                            obj.elec_ch_label,...
                            obj.elec_ch_clean,...
                            obj.elec_ch_valid,...
                            'stitch_index',obj.stitch_index,...
                            'sample_freq',obj.sample_freq,...
                            'downsample',true,...
                            'decimation_freq',obj.for_preproc.decimation_freq...
            );
        end

        % only do visual inspection if hasn't been done already
        if ~ops.doneVisualInspection
            prompt1 = '\nUSER INPUT REQUIRED: \nAdditional channels to remove from visual inspection? (format: [1,2]) - ';
            prompt2 = 'Your name please :) - ';

            obj.elec_ch_user_deselect = input(prompt1)';
        
            % output structure for object
            vi_ops.inspected = 1;
            vi_ops.inspected_by = input(prompt2,'s');
            vi_ops.inspection_date = datestr(now, 'yyyy/mm/dd-HH:MM');

            obj.for_preproc.visualInspection_results = vi_ops;
        end
        
        % update clean/valid electrodes
        obj.define_clean_channels();

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REMOVE COMMMON NOISE AND REFERENCE SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function reference_signal(obj,varargin)
        % Should be run after visually_inspect
        % TODO - decide what to do with the EEG channels
        
        p = inputParser();
        addParameter(p,'doGlobalMeanRemoval',false)
        addParameter(p,'doCAR',false);
        addParameter(p,'doShankCSR',false);
        addParameter(p,'doBipolarReferencing',false);
        parse(p, varargin{:});
        ops = p.Results;

        signal = obj.elec_data';
        signal_for_bip = obj.elec_data';

        % split ECoG (grids+srips) and SEEG channels
        ecog_chans = find(strcmp(obj.elec_ch_type, 'ecog_grid') | strcmp(obj.elec_ch_type, 'ecog_strip'));
        seeg_chans = find(strcmp(obj.elec_ch_type, 'seeg'));

        % initialize bipolar signal
        signal_bipolar = zeros(size(signal)); % might not be used

        for k=1:length(obj.stitch_index) % number of separate data files with signal
            fprintf(1, '\n> Referencing signal from file %d of %d ... \n',k,length(obj.stitch_index));
            
            if k == length(obj.stitch_index) % signal for file stops at end of matrix
                stop = size(signal,1);
            else % signal for file stops before stitch index of next file
                stop = obj.stitch_index(k+1)-1;
            end


            % --------------------
            % GLOBAL MEAN REMOVAL
            % --------------------
            if ops.doGlobalMeanRemoval
                fprintf(1, '\n>> Removing global mean of signal \n');

                signal_ = signal(obj.stitch_index(k):stop,:);

                overall_mean = mean(signal_(:,obj.elec_ch_clean),2);
                signal_ = signal_ - repmat(overall_mean,1,size(signal_,2));

                signal(obj.stitch_index(k):stop,:) = signal_;

            end


            % ---------------------------------------
            % COMMON AVERAGE REFERENCING (ECoG only) <-- NOT ANYMORE
            % ---------------------------------------
            if ops.doCAR && ~isempty(ecog_chans)
                fprintf(1, '\n>> Common average filtering signal \n');
                fprintf(1,'[');

                signal_ = signal(obj.stitch_index(k):stop,:);

                % determine the number of amps
                num_amps = ceil(size(signal_,2) / obj.for_preproc.elecs_per_amp);

                % exclude the channels that were noisy or not ECoG
                eligible_channels_ecog = intersect(ecog_chans,obj.elec_ch_clean);
                % eligible_channels = obj.elec_ch_clean;

                % for each of these amps
                for idx_amp = 1:num_amps
                    idx_low  = (idx_amp-1)*obj.for_preproc.elecs_per_amp+1;
                    idx_high = min((idx_amp-0)*obj.for_preproc.elecs_per_amp+0, max(eligible_channels_ecog));
        
                    % exclude the channels that were not on this amp
                    list_channels = intersect(eligible_channels_ecog,idx_low:idx_high);
        
                    % check if any channels are left
                    if ~isempty(list_channels) && length(list_channels)>1
                        % calculate the common average reference signal
                        signal_mean = mean(signal_(:,list_channels),2);
                        % subtract the common average signal from each channel of amp
                        for idx_ch = list_channels
                            signal_(:,idx_ch) = signal_(:,idx_ch) - signal_mean;
                            fprintf(1,'.');
                        end
                    else
                        % don't change the signal at all
                        for idx_ch = list_channels
                            signal_(:,idx_ch) = signal_(:,idx_ch);
                            fprintf(1,'.');
                        end
                    end
                end

                fprintf(1,'] done\n');

                signal(obj.stitch_index(k):stop,:) = signal_;

            % elseif ops.doCAR
            %     error('No ECoG channels to perform CAR on')
            end


            % ----------------------------------------
            % SHANK COMMON SOURCE REMOVAL (SEEG only)
            % ----------------------------------------
            if ops.doShankCSR && ~isempty(seeg_chans)
                fprintf(1, '\n>> Shank common source removal \n');
                fprintf(1,'[');

                signal_ = signal(obj.stitch_index(k):stop,:);

                % exclude the channels that were noisy or not SEEG
                eligible_channels_seeg = intersect(seeg_chans,obj.elec_ch_clean);

                if ops.doCAR & ~isempty(ecog_chans)
                    assert(isempty(intersect(eligible_channels_ecog,eligible_channels_seeg)),"Some channels are labelled as both ECoG and SEEG.")
                end

                % split labels into shanks
                [shank_locs,~,~] = obj.extract_shanks();
                
                % for each shank
                for idx_shk = 1:size(shank_locs,1);

                    % exlude channels that were not on this shank
                    same_shank = shank_locs{idx_shk};
                    same_shank = intersect(eligible_channels_seeg,same_shank);

                    % check if any channels are left
                    if ~isempty(same_shank) && length(same_shank)>1
                        % calculate the common average reference signal
                        signal_mean = mean(signal_(:,same_shank),1);
                        % subtract the common average signal from each channel of shank
                        for idx_ch = same_shank
                            signal_(:,idx_ch) = signal_(:,idx_ch) - signal_mean;
                            fprintf(1,'.');
                        end
                    else
                        % don't change the signal at all 
                        % (either none clean, only one clean, or none SEEG)
                        for idx_ch = same_shank
                            signal_(:,idx_ch) = signal_(:,idx_ch);
                            fprintf(1,'.');
                        end
                    end
                end

                signal(obj.stitch_index(k):stop,:) = signal_;

                fprintf(1,'] done\n');

            elseif ops.doShankCSR 
                error('No SEEG channels to perform shank CSR on')
            end


            % --------------------------------
            % BIPOLAR REFERENCING (SEEG only)
            % --------------------------------
            if ops.doBipolarReferencing && ~isempty(seeg_chans)
                fprintf(1, '\n>> Bipolar referencing signal \n');
                fprintf(1,'[');
                
                signal_ = signal_for_bip(obj.stitch_index(k):stop,:);

                % TODO - take the unneccesary bipolar stuff out of the loop

                % exclude the channels that were noisy or not SEEG
                eligible_channels_seeg = intersect(seeg_chans,obj.elec_ch_clean);
                
                if ops.doCAR & ~isempty(ecog_chans)
                    assert(isempty(intersect(eligible_channels_ecog,eligible_channels_seeg)),"Some channels are labelled as both ECoG and SEEG.")
                end

                % split labels into shanks 
                [shank_locs,~,shank_nums] = obj.extract_shanks();

                % find bipolar pairs of clean, SEEG electrodes
                % (i.e., electrodes on the same SEEG shank that are one apart and clean)
                chan_idx_for_bip = cellfun(@(x) intersect(x,eligible_channels_seeg,'stable'),shank_locs,'uni',false);
                chan_num_for_bip = cellfun(@(x) shank_nums(x),chan_idx_for_bip,'uni',false);
                % as cells 
                bipolar_diffs_idx_grp = cellfun(@(x,y) [y(find(diff(x)==1 & diff(y)==1))+1,y(find(diff(x)==1 & diff(y)==1))],chan_num_for_bip,chan_idx_for_bip,'uni',false)';
                bipolar_diffs_name_grp = cellfun(@(x) obj.elec_ch_label(x),bipolar_diffs_idx_grp,'uni',false);
                % as mats
                bipolar_diffs_idx = cell2mat(bipolar_diffs_idx_grp');
                bipolar_diffs_name = obj.elec_ch_label(bipolar_diffs_idx);
                
                bipolar_idxs = 1:size(bipolar_diffs_idx,1);
                bipolar_valid = ones(size(bipolar_diffs_idx,1),1);

                % create a biopolar version of the signal
                signal_bipolar_= double([]); 
                for bipolar_id = 1:size(bipolar_diffs_idx,1)
                    bipol_ch_1 = bipolar_diffs_idx(bipolar_id,1);
                    bipol_ch_2 = bipolar_diffs_idx(bipolar_id,2);
                    bip_ch_name1 = bipolar_diffs_name{bipolar_id,1};
                    bip_ch_name2 = bipolar_diffs_name{bipolar_id,2};
              
                    % basic assertions
                    B = cell2mat(extract(bip_ch_name2,lettersPattern));
                    A =cell2mat(extract(bip_ch_name1,lettersPattern));
                    assert(all(A==B),"Some channels are not on the same shank");
                    B =str2num(cell2mat(extract(bip_ch_name2,digitsPattern)));
                    A =str2num(cell2mat(extract(bip_ch_name1,digitsPattern)));
                    assert(A-B==1, "Some channels are more than one apart");
                    
                    % subtract channel 2 from channel 1 in pair
                    signal_bipolar_(:,bipolar_id) = signal_(:,bipol_ch_1)-signal_(:,bipol_ch_2);
                    fprintf(1,'.');
                end

                bipolar_diffs_name = arrayfun(@(x) {[bipolar_diffs_name{x,1} '-' bipolar_diffs_name{x,2}]},[1:size(bipolar_diffs_name,1)])';

                signal_bipolar(obj.stitch_index(k):stop,1:size(bipolar_diffs_idx,1)) = signal_bipolar_;

                fprintf(1,'] done\n');

            elseif ops.doBipolarReferencing
                error('No SEEG channels to perform bipolar referencing on')
            end
        end

        
        % unipolar 
        if ops.doGlobalMeanRemoval || ops.doCAR || ops.doShankCSR
            obj.elec_data = signal';
        end

        % bipolar
        if ops.doBipolarReferencing 

            % remove unused rows in bipolar signal
            signal_bipolar = signal_bipolar(:,1:size(bipolar_diffs_idx,1));

            obj.bip_elec_data    = signal_bipolar';
            obj.bip_ch           = bipolar_idxs';
            obj.bip_ch_label     = bipolar_diffs_name; 
            obj.bip_ch_valid     = bipolar_valid;
            obj.bip_ch_grp       = bipolar_diffs_idx_grp';
            obj.bip_ch_label_grp = bipolar_diffs_name_grp';
              
            if obj.for_preproc.isPlotVisible
                % plot to make sure signal looks right
                obj.plot_channels(signal_bipolar,...
                                obj.bip_ch_label,...
                                obj.bip_ch,...
                                obj.bip_ch_valid,...
                                'stitch_index',obj.stitch_index,...
                                'sample_freq',obj.sample_freq,...
                                'downsample',true,...
                                'decimation_freq',obj.for_preproc.decimation_freq...
                );
            end

        end

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT HIGH GAMMA SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function extract_high_gamma(obj,varargin)
        % Extracts high gamma signal using one of two methods: 
        %   (1) Chang Lab method (gaussian filters)
        %           hilbert_transform.m is a special function 
        %           See 'HELPER SIGNAL PREPROCESSING FUNCTIONS'
        %   (2) Standard bandpass filter
        %           hilbert.m is a standard MATLAB function
        % 
        % TODO - descriptions of above methods

        p = inputParser();
        addParameter(p,'doGaussianFilterExtraction',false);
        addParameter(p,'doBandpassExtraction',false);
        parse(p, varargin{:});
        ops = p.Results;

        % ensure that at least one, and only one, extraction method is specified
        if ops.doGaussianFilterExtraction == ops.doBandpassExtraction
            error('Must specify one - and only one - method for extracting the high gamme signal')
        end

        signal = obj.elec_data';
        if ~isempty(obj.bip_elec_data)
            signal_bipolar = obj.bip_elec_data';
        end

        % define filter params if it hasn't already been done
        if ~isfield(obj.for_preproc,'gaussian') || ~isfield(obj.for_preproc,'bandpass')
            obj.define_parameters()
        end


        for k=1:length(obj.stitch_index) % number of separate data files with signal
            fprintf(1, '\n> Extracting high gamma signal from file %d of %d ... \n',k,length(obj.stitch_index));

            if k == length(obj.stitch_index) % signal for file stops at end of matrix
                stop = size(signal,1);
            else % signal for file stops before stitch index of next file
                stop = obj.stitch_index(k+1)-1;
            end

            signal_ = signal(obj.stitch_index(k):stop,:);


            % ------------------------------------
            % EXTRACTION USING GUASSIAN FILTERING
            % ------------------------------------
            if ops.doGaussianFilterExtraction

                cfs = obj.for_preproc.gaussian.cfs;
                sds = obj.for_preproc.gaussian.sds;

                % define gaussian filters
                % filter_bank = cell(numel(cfs),1);
                filter_bank={};
                for s=1:length(cfs)
                    filter_bank{s} = obj.gaussian_filter(transpose(signal_(:,1)),obj.sample_freq,cfs(s),sds(s));
                end
                obj.for_preproc.gaussian.filter_banks = filter_bank;
                % size(filter_bank)

                % --- UNIPOLAR ---
                fprintf(1, '\n>> Extracting unipolar high gamma envelope based on gaussian filtering \n');
                fprintf(1,'[');

                signal_hilbert = nan*signal_;

                for kk=1:size(signal_,2)
                    signal_hilbert_all = cell2mat(cellfun(@abs,obj.hilbert_transform(double(transpose(signal_(:,kk))),obj.sample_freq,filter_bank),'UniformOutput',false));
                    % size(signal_hilbert_all)
                    signal_hilbert(:,kk) = transpose(mean(signal_hilbert_all,1));

                    fprintf(1,'.');

                end 

                fprintf(1,'] done\n');

                signal(obj.stitch_index(k):stop,:) = signal_hilbert;


                % --- BIPOLAR ---
                if ~isempty(obj.bip_elec_data)
                    fprintf(1, '\n>> Extracting bipolar high gamma envelope based on gaussian filtering \n');
                    fprintf(1,'[');

                    signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                    signal_hilbert_bipolar = nan*signal_bipolar_;

                    for kk=1:size(signal_bipolar_,2)
                        signal_hilbert_bipolar_all = cell2mat(cellfun(@abs,obj.hilbert_transform(double(transpose(signal_bipolar_(:,kk))),obj.sample_freq,filter_bank),'UniformOutput',false));
                        
                        signal_hilbert_bipolar(:,kk) = transpose(mean(signal_hilbert_bipolar_all,1));

                        fprintf(1,'.');

                    end

                    fprintf(1,'] done\n')

                    signal_bipolar(obj.stitch_index(k):stop,:) = signal_hilbert_bipolar;
                    
                end


            % ---------------------------------------------
            % EXTRACTION USING STANDARD BANDPASS FILTERING
            % ---------------------------------------------
            elseif ops.doBandpassExtraction

                B = obj.for_preproc.bandpass.B;
                A = obj.for_preproc.bandpass.A;

                % --- UNIPOLAR ---
                fprintf(1, '\n>> Extracting unipolar high gamma envelope based on bandpass filtering \n');

                % apply filter
                signal_hilbert = filtfilt(B,A,double(signal_));

                % measure envelope
                signal_hilbert = abs(hilbert(signal_hilbert));

                % truncate
                signal_hilbert(signal_hilbert < 0) = 0;


                % --- BIPOLAR ---
                if ~isempty(obj.bip_elec_data)
                    fprintf(1, '\n>> Extracting bipolar high gamma envelope based on bandpass filtering \n');

                    signal_hilbert_bipolar = filtfilt(B,A,double(signal_bipolar_));
                    signal_hilbert = abs(hilbert(signal_hilbert));
                    signal_hilbert(signal_hilbert < 0) = 0;

                    signal_bipolar(obj.stitch_index(k):stop,:) = signal_hilbert_bipolar;
                end
            
            end

        end

        obj.elec_data = signal';

        if ~isempty(obj.bip_elec_data)
            obj.bip_elec_data = signal_bipolar';
        end

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Z-SCORE SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function zscore_signal(obj)
        % Zscores signal using default MATLAB zscore function.

        signal = obj.elec_data';
        if ~isempty(obj.bip_elec_data)
            signal_bipolar = obj.bip_elec_data';
        end

        for k=1:length(obj.stitch_index) % number of separate data files with signal
            fprintf(1, '\n> Computing z-score of signal from file %d of %d ... \n',k,length(obj.stitch_index));


            if k == length(obj.stitch_index) % signal for file stops at end of matrix
                stop = size(signal,1);
            else % signal for file stops before stitch index of next file
                stop = obj.stitch_index(k+1)-1;
            end

            % --- UNIPOLAR ---
            signal_ = signal(obj.stitch_index(k):stop,:);
            signal_ = zscore(signal_);
            signal(obj.stitch_index(k):stop,:) = signal_;
        
            % --- BIPOLAR ---
            if ~isempty(obj.bip_elec_data)
                signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                signal_bipolar_ = zscore(signal_bipolar_);
                signal_bipolar(obj.stitch_index(k):stop,:) = signal_bipolar_;
            end

        end
        
        obj.elec_data = signal';

        if ~isempty(obj.bip_elec_data)
            obj.bip_elec_data = signal_bipolar';
        end


    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DOWNSAMPLE SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function downsample_signal(obj,varargin)
        % Downsamples signal and all other relevant parts of object
        % 

        p = inputParser();
        addParameter(p,'decimationFreq',obj.for_preproc.decimation_freq)
        parse(p, varargin{:});
        ops = p.Results;
        
        signal = obj.elec_data';
        if ~isempty(obj.bip_elec_data)
            signal_bipolar = obj.bip_elec_data';
        end

        % assert(obj.sample_freq > obj.for_preproc.decimation_freq,'signal has already been downsampled');

        signal_dec = [];
        signal_bipolar_dec = [];

        curr_stitch = 1;
        stitch_index = []; % first sample in each data file

        for k=1:length(obj.stitch_index) % number of separate data files with signal
            fprintf(1, '\n> Resampling signal from file %d of %d ... \n',k,length(obj.stitch_index));

            if k == length(obj.stitch_index) % signal for file stops at end of matrix
                stop = size(signal,1);
            else % signal for file stops before stitch index of next file
                stop = obj.stitch_index(k+1)-1;
            end

            % --- UNIPOLAR ---
            signal_ = signal(obj.stitch_index(k):stop,:);
            signal_ = resample(double(signal_),ops.decimationFreq,obj.sample_freq);
            signal_dec = [signal_dec; signal_];
        
            % --- BIPOLAR ---
            if ~isempty(obj.bip_elec_data)
                signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                signal_bipolar_ = resample(double(signal_bipolar_),ops.decimationFreq,obj.sample_freq);
                signal_bipolar_dec = [signal_bipolar_dec; signal_bipolar_];
            end

            % update stitch index
            stitch_index = [stitch_index; curr_stitch];
            curr_stitch = curr_stitch + size(signal_,1);

        end

        obj.elec_data = signal_dec';

        if ~isempty(obj.bip_elec_data)
            obj.bip_elec_data = signal_bipolar_dec';
        end

        % construct new trial timing table if downsample rate is not the preset downsample rate
        % (NOT ADVISED)
        if ~isempty(obj.trial_timing) && (ops.decimationFreq ~= obj.for_preproc.decimation_freq)
            decimation_factor = obj.sample_freq / ops.decimationFreq;
            trial_timing = cell(size(obj.trial_timing));
            for i=1:size(trial_timing,1)
                tmp_table = obj.trial_timing{i,1};
                tmp_table.start = round(tmp_table.start / decimation_factor);
                tmp_table.end = round(tmp_table.end / decimation_factor);
                trial_timing{i,1} = tmp_table;
            end
            obj.trial_timing = trial_timing;
        end

        % set trial timing table to the downsampled version if exists
        if ~isempty(obj.trial_timing) && (ops.decimationFreq == obj.for_preproc.decimation_freq)
            obj.trial_timing = obj.for_preproc.trial_timing_dec;
        end

        obj.stitch_index = stitch_index;
        obj.sample_freq = ops.decimationFreq;

    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REMOVE OUTLIERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function remove_outliers(obj)
        % Removes outliers from envelope.
        %
        % CANNOT be done as the first preprocessing step. MUST be done 
        % after high gamma extraction (need envelopes).

        ops.trimmed      = obj.for_preproc.outlier.trimmed;
        ops.threshold    = obj.for_preproc.outlier.threshold;
        ops.percentile   = obj.for_preproc.outlier.percentile;
        ops.buffer       = obj.for_preproc.outlier.buffer;
        ops.interpMethod = obj.for_preproc.outlier.interpMethod;

        samples_to_remove = ops.trimmed * obj.sample_freq;

        signal = obj.elec_data';
        signal = obj.zero_out_signal(signal);

        outlierRemoval_results.idxs    = [];
        outlierRemoval_results.prcnts  = [];
        outlierRemoval_results.ignored = [];

        if ~isempty(obj.bip_elec_data)
            signal_bipolar = obj.bip_elec_data';
            signal_bipolar = obj.zero_out_signal(signal_bipolar);

            outlierRemoval_results.idxs_bipolar    = [];
            outlierRemoval_results.prcnts_bipolar  = [];
            outlierRemoval_results.ignored_bipolar = [];
        end
        

        for k=1:length(obj.stitch_index) % number of separate data files with signal
            fprintf(1, '\n> Removing outliers from file %d of %d ... \n',k,length(obj.stitch_index));
            fprintf(1,'[');

            if k == length(obj.stitch_index) % signal for file stops at end of matrix
                stop = size(signal,1);
            else % signal for file stops before stitch index of next file
                stop = obj.stitch_index(k+1)-1;
            end

            % --- UNIPOLAR ---
            signal_ = signal(obj.stitch_index(k):stop,:);
            [signal_,idxs,prcnts,ignored] = obj.envelope_outliers(signal_,ops);
            signal(obj.stitch_index(k):stop,:) = signal_;
            fprintf(1,'] done\n');

            outlierRemoval_results.idxs = [outlierRemoval_results.idxs, idxs];
            outlierRemoval_results.prcnts = [outlierRemoval_results.prcnts, prcnts];
            outlierRemoval_results.ignored = [outlierRemoval_results.ignored; ismember(obj.elec_ch,ignored)'];
            
            ignored_not_noisy = intersect(obj.elec_ch_clean,ignored);
            fprintf(1,'Unequal rise and fall of outliers, ignored unipolar channels: ');
            fprintf(1,'%d ',ignored_not_noisy(:)); fprintf('\n');

            % --- BIPOLAR ---
            if ~isempty(obj.bip_elec_data)
                fprintf(1,'[');
                signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                [signal_bipolar_,idxs_bipolar,prcnts_bipolar,ignored_bipolar] = obj.envelope_outliers(signal_bipolar_,ops);
                signal_bipolar(obj.stitch_index(k):stop,:) = signal_bipolar_;
                fprintf(1,'] done\n');

                outlierRemoval_results.idxs_bipolar = [outlierRemoval_results.idxs_bipolar, idxs_bipolar];
                outlierRemoval_results.prcnts_bipolar = [outlierRemoval_results.prcnts_bipolar, prcnts_bipolar];
                outlierRemoval_results.ignored_bipolar = [outlierRemoval_results.ignored_bipolar, ismember(obj.bip_ch,ignored_bipolar)'];

                ignored_bipolar_not_noisy = intersect(obj.bip_ch,ignored_bipolar);
                fprintf(1,'Unequal rise and fall of outliers, ignored bipolar channels: ');
                fprintf(1,'%d ',ignored_bipolar_not_noisy(:)); fprintf('\n');
            end

        end

        obj.elec_data = signal';

        if ~isempty(obj.bip_elec_data)
            obj.bip_elec_data = signal_bipolar';
        end

        obj.for_preproc.outlierRemoval_results = outlierRemoval_results;

        if obj.for_preproc.isPlotVisible
            obj.plot_channels(signal,...
                            obj.elec_ch_label,...
                            obj.elec_ch_clean,...
                            obj.elec_ch_valid,...
                            'stitch_index',obj.stitch_index,...
                            't_len',60,...
                            'sample_freq',obj.sample_freq...
            );
        end

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADD ANATOMY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function add_anatomy(obj,anatomy_path,varargin)
        % Adds anatomy files to object, including a mapping from the channels
        % in the anatomy file to the channels in the object.
        %
        % Assumes particular naming convention for files 
        % TODO: change to be more general

        p = inputParser();
        addParameter(p,'veraFolder',true); % if anatomy_path refers to directory with vera folders or directly to anatomical files
        addParameter(p,'subdirName','VERA_'); % subject name is appended to end of this prefix
        addParameter(p,'templateName','cvs_avg35_inMNI152'); % same folder and file name (with .mat appended)
        parse(p, varargin{:});
        ops = p.Results;

        % subdirectory name (if applicable)
        if ops.veraFolder
            folder = [ops.subdirName obj.subject filesep];
            template_folder = [ops.templateName filesep];
        else
            folder = ''; % path refers to folder with all anatomical files
            template_folder = '';
        end

        % subject-specific space
        filename = [anatomy_path folder obj.subject '_brain.mat'];
        subject_space = load(filename);
        obj.anatomy.subject_space = subject_space;

        % MNI space
        filename = [anatomy_path folder obj.subject '_MNI_brain.mat'];
        mni_space = load(filename);
        obj.anatomy.mni_space = mni_space;

        % template brain
        filename = [anatomy_path template_folder ops.templateName '.mat'];
        template_brain = load(filename);
        obj.anatomy.template_brain = template_brain;

        % add mapping from anatomical files to object
        [mapping,labels] = obj.channel_mapping_anatomical(subject_space);
        obj.anatomy.mapping = mapping;
        obj.anatomy.labels = labels;
        
        % add hemisphere labels
        hemisphere = cell(size(subject_space.tala.electrodes,1),1);
        right_idxs = find(subject_space.tala.electrodes(:,1) > 0);
        hemisphere(right_idxs,1) = {'right'};
        left_idxs = find(subject_space.tala.electrodes(:,1) < 0);
        hemisphere(left_idxs,1) = {'left'};
        obj.anatomy.hemisphere = hemisphere;

    end


    %% HELPER SIGNAL PREPROCESSING METHODS


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIRST STEP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = first_step(obj,varargin)
        % Reset preprocessed parameters to raw values and re-define 
        % parameters (e.g., filter params).
        % 
        % Returns : obj

        p = inputParser();
        addParameter(p,'doneVisualInspection',false);
        parse(p, varargin{:});
        ops = p.Results;

        obj.stitch_index = obj.for_preproc.stitch_index_raw;
        obj.sample_freq  = obj.for_preproc.sample_freq_raw;

        % re-define parameters
        obj.define_parameters();

        % % check if visual inspection needs to be done again
        % if isfield(obj.for_preproc,'visual_inspection') % if it has been done before
        %     if any(cellfun(@(x) strcmp(x,'visualInspection'),obj.for_preproc.order,'UniformOutput',false))
        %         obj.for_preproc.visual_inspection = []; % only clear if visual inspection will be done again
        %         obj.elec_ch_user_deselect = [];
        %     else % use old visual inspection 
        %         fprintf(1,'\nChannels to be removed because of previously completed visual inspection: ');
        %         fprintf(1,'%d ',obj.elec_ch_user_deselect(:)); fprintf(1,'\n');
        %     end
        % end

        % reset noisy electrodes
        obj.elec_ch_with_IED = [];
        obj.elec_ch_with_noise = [];
        if ~ops.doneVisualInspection
            obj.elec_ch_user_deselect = [];
        end

        % define clean electrodes
        obj.define_clean_channels()

        % clear results from preprocessing steps
        obj.for_preproc.notchFilter_results = [];
        obj.for_preproc.IEDRemoval_results = [];
        obj.for_preproc.visualInspection_results = [];
        obj.for_preproc.outlierRemoval_results = [];

        % reset all bipolar info
        obj.bip_elec_data    = [];
        obj.bip_ch           = [];
        obj.bip_ch_label     = [];
        obj.bip_ch_valid     = [];
        obj.bip_ch_grp       = [];
        obj.bip_ch_label_grp = [];

        % update trial timing table to correspond to the original signal 
        if isfield(obj.for_preproc,'trial_timing_raw') % might not exist if first time preproc
            obj.trial_timing = obj.for_preproc.trial_timing_raw;
        end

        obj.elec_data = obj.for_preproc.elec_data_raw;

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINE PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function define_parameters(obj)
        % Define parameters to be used in preprocessing and store them 
        % in obj.for_preproc. This is the ONLY function where parameters 
        % should be modified by the user.
        % 
        % Returns : obj


        % SET FILTER PARAMETERS

        param = struct;
            
        % --- highpass filter ---
        param.highpass.Wp = 0.50; % Hz
        param.highpass.Ws = 0.05; % Hz
        param.highpass.Rp = 3;    % dB
        param.highpass.Rs = 30;   % dB
            
        % --- IIR peak filter ---
        param.peak.fcenter = [55,60,65];
        param.peak.bw      = ones(1,length(param.peak.fcenter)).*0.001;
            
        % --- notch filter ---
        param.notch.fcenter = [60,120,180,240];
        param.notch.bw = ones(1,length(param.notch.fcenter)).*0.001;

        % --- gaussian filter ---
        param.gaussian.f_gamma_low = 70;
        param.gaussian.f_gamma_high = 150;

        % --- bandpass filter ---
        param.bandpass.f_gamma_low = 70;
        param.bandpass.f_gamma_high = 150;
        param.bandpass.filter_order = 6;


        % CONSTRUCT FILTERS (shouldn't need to edit)

        % --- highpass filter ---
        highpass.Wp = param.highpass.Wp/(obj.sample_freq/2); 
        highpass.Ws = param.highpass.Ws/(obj.sample_freq/2);
        highpass.Rp = param.highpass.Rp; 
        highpass.Rs = param.highpass.Rs;

        [highpass.n,highpass.Wn] = buttord(highpass.Wp,highpass.Ws,highpass.Rp,highpass.Rs);
        highpass.n = highpass.n + rem(highpass.n,2);

        % caclulate the filter coefficients in Zero-Pole-Gain design
        [highpass.z,highpass.p,highpass.k] = butter(highpass.n,highpass.Wn,'high');
        [highpass.sos,highpass.g] = zp2sos(highpass.z,highpass.p,highpass.k);
        highpass.h = dfilt.df2sos(highpass.sos,highpass.g);

        % --- IIR peak filter ---
        % calculate the IIR-peak filter coefficients in a,b format
        for idx = 1:length(param.peak.fcenter)
            peak{idx}.wo = param.peak.fcenter(idx)/(obj.sample_freq/2);
            peak{idx}.bw = param.peak.bw(idx);  
            [peak{idx}.b,peak{idx}.a] = iirpeak(peak{idx}.wo,peak{idx}.bw);
        end

        % --- notch filter ---
        % calculate the IIR-peak filter coefficients in a,b format 
        for idx = 1:length(param.notch.fcenter)
            notch{idx}.wo = param.notch.fcenter(idx)/(obj.sample_freq/2);  
            notch{idx}.bw = param.notch.bw(idx);
            [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);  
        end 

        % --- gaussian filter ---
        %
        % See:
        %   Dichter, Benjamin K., Jonathan D. Breshears, Matthew K. Leonard, and Edward F. Chang. 2018.
        %   The Control of Vocal Pitch in Human Laryngeal Motor Cortex. Cell 174 (1).  
        [gaussian.cfs,gaussian.sds] = obj.get_filter_param_chang_lab(param.gaussian.f_gamma_low,param.gaussian.f_gamma_high);
        gaussian.cfs(1) = 73.0; % correction to stay away from 60 hz
        guassian.f_gamma_low = param.gaussian.f_gamma_low;
        gaussian.f_gamma_high = param.gaussian.f_gamma_high;

        % --- bandpass filter ---
        % Construct an FDESIGN object and call its BUTTER method.
        bandpass.h = fdesign.bandpass('N,F3dB1,F3dB2',param.bandpass.filter_order,param.bandpass.f_gamma_low,param.bandpass.f_gamma_high,obj.sample_freq);
        bandpass.Hd = design(bandpass.h,'butter');
        [bandpass.B, bandpass.A] = sos2tf(bandpass.Hd.sosMatrix,bandpass.Hd.scaleValues);
        bandpass.filter_order = param.bandpass.filter_order;
        bandpass.f_gamma_low = param.bandpass.f_gamma_low;
        bandpass.f_gamma_high = param.bandpass.f_gamma_high;


        % OTHER PARAMETERS

        % zero out signal
        zero = 1; % seconds

        % outlier removal
        outlier.trimmed = 1; % second
        outlier.threshold = 5;
        outlier.percentile = 0.9;
        outlier.buffer = 20;
        outlier.interpMethod = 'linear';


         % add fields to obj.info
        % WILL OVERWRITE OLD FILTER PARAMETERS
        obj.for_preproc.filter_params   = param;
        obj.for_preproc.highpass        = highpass;
        obj.for_preproc.peak            = peak;
        obj.for_preproc.notch           = notch;
        obj.for_preproc.gaussian        = gaussian;
        obj.for_preproc.bandpass        = bandpass;
        obj.for_preproc.zero            = zero;
        obj.for_preproc.outlier         = outlier;

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ZERO OUT SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function signal=zero_out_signal(obj,signal)
        % Zeros out first couple seconds of each file. 
        % 
        % Returns : signal

        samples_to_remove = obj.sample_freq*obj.for_preproc.zero;

        for k=1:length(obj.stitch_index) % number of separate data files with signal

            if k == length(obj.stitch_index) % signal for file stops at end of matrix
                stop = size(signal,1);
            else % signal for file stops before stitch index of next file
                stop = obj.stitch_index(k+1)-1;
            end

            signal_ = signal(obj.stitch_index(k):stop,:);
            mean_signal = mean(signal_,1);
            signal_(1:samples_to_remove,:) = zeros(samples_to_remove,size(signal_,2));
            signal_((end-samples_to_remove+1):end,:) = zeros(samples_to_remove,size(signal_,2));
            signal(obj.stitch_index(k):stop,:) = signal_;

        end

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEASURE LINE-NOISE POWER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function signal_noise=measure_line_noise(obj,signal)
        % 
        % 
        % Returns : measurement of signal noise

        fprintf(1, '\n> Measuring 60Hz noise power ...\n');
        fprintf(1,'[');

        peak = obj.for_preproc.peak;
    
        signal_noise = zeros(size(signal,2),length(peak));

        % calculate average root-mean-square of the line-noise
        for idx_channel=1:size(signal,2)
            % signal_noise(idx_channel,1) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
            for idx_filter=1:length(peak)
                signal_noise(idx_channel,idx_filter) = mean(abs(filter(peak{idx_filter}.b,peak{idx_filter}.a,signal(:,idx_channel))));
            end
            fprintf(1,'.');
        end

        fprintf(1,'] done\n');

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT LINE-NOISE POWER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function f=plot_line_noise(obj,noise_before,noise_after)
        % 

        fprintf(1, '\n> Plotting 60Hz noise power ...\n');

        close all
        f = figure;
        set(gcf,'position',[30,30,2300,900]);
        c= [0.4660 0.6740 0.1880];

        x = find(~obj.elec_ch_valid);
        idxs = ~obj.elec_ch_valid;

        % --- ALL NOISE BEFORE ---
        currsub = subplot(2,2,1);       
        stem(noise_before,'filled'); 
        axis tight; hold on;
        stem(x,noise_before(idxs,:),'filled','Color','k')
        legend({'55Hz noise','60Hz noise','65Hz noise','MARKED NOISY'},'Location','best','FontSize',16,'Box','off');
        ylabel('Noise (uV)','FontSize',18);
        title('BEFORE NOTCH FILTERING','FontSize',22);
        obj.update_position(currsub);

        % --- NOISE RATIO BEFORE ---
        currsub = subplot(2,2,3);
        stem(noise_before(:,2)./mean(noise_before(:,[1,3]),2),'filled','Color',c); 
        axis tight; hold on;
        stem(x,noise_before(idxs,2)./mean(noise_before(idxs,[1,3]),2),'filled','Color','k')
        xlabel('Channel #','FontSize',18);
        ylabel('60Hz noise / mean 55Hz+65Hz noise','FontSize',18)
        obj.update_position(currsub);

        % --- ALL NOISE AFTER ---
        currsub = subplot(2,2,2);
        stem(noise_after,'filled'); 
        axis tight; hold on;
        stem(x,noise_after(idxs,:),'filled','Color','k')
        title('AFTER NOTCH FILTERING','FontSize',22)
        obj.update_position(currsub);

        % --- NOISE RATIO AFTER ---
        currsub = subplot(2,2,4);
        stem(noise_after(:,2)./mean(noise_after(:,[1,3]),2),'filled','Color',c); 
        axis tight; hold on;
        stem(x,noise_after(idxs,2)./mean(noise_after(idxs,[1,3]),2),'filled','Color','k')
        xlabel('Channel #','FontSize',18)
        obj.update_position(currsub);

        % --- SAVE ---
        PATH = 'plots/line_noise/';
        filename = split(obj.for_preproc.log_file_name,'/');
        filename = split(filename{end},'.');
        filename = [filename{1} '_line_noise.png'];
        saveas(gcf,strcat(PATH,filename));
        set(0, 'CurrentFigure', f);

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UPDATE SUPLOT POSITION FOR PLOTTING LINE NOISE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function update_position(obj,currsub)
        
        pos = get(currsub, 'Position');
        new_pos = pos + [-0.05 -0.05 0.07 0.07];
        set(currsub, 'Position', new_pos)

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINE CLEAN CHANNELS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function define_clean_channels(obj)
        % Defines the clean channels to be used for analysis.
        %
        % Returns : obj
        
        bad_elecs = union(obj.elec_ch_prelim_deselect,obj.elec_ch_with_IED);
        bad_elecs = union(bad_elecs,obj.elec_ch_with_noise);
        bad_elecs = union(bad_elecs,obj.elec_ch_user_deselect);

        obj.elec_ch_clean = setdiff(obj.elec_ch,bad_elecs,'stable');
        obj.elec_ch_valid = ismember(obj.elec_ch,obj.elec_ch_clean);
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT SHANK INFORMATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [locs,ch_tags,ch_num]=extract_shanks(obj)
        
        ch_num = cellfun(@str2num,extract(obj.elec_ch_label,digitsPattern)); % nums
        ch_tags = extract(obj.elec_ch_label,lettersPattern); % letters from label
        % ch_tags = cellfun(@(x) [x '_'],ch_tags,'uni',false); % append underscore 
        [tags,~,~] = unique(ch_tags,'stable');
        locs = cellfun(@(x) find(ismember(ch_tags,x))',tags,'uni',false);
    end
    

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINE CHANG LAB FILTER PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [cfs,sds]=get_filter_param_chang_lab(obj,f_low,f_high)
        % TODO - description

        % Standard signal bands for neuroscience
        %   bands = ['theta', 'alpha', 'beta', 'gamma', 'high gamma']
        %   min_freqs = [4., 8., 15., 30., 70.]
        %   max_freqs = [7., 14., 29., 59., 150.]
        %   HG_freq = 200.

        cfs_round_factor = 1;
        cfs_round_val = 10^(cfs_round_factor);

        sds_round_factor = 2;
        sds_round_val = 10^(sds_round_factor);
        
        fq_min = 4.0749286538265;
        fq_max = 200.;
        scale = 7.;

        cfs = 2.^((log2(fq_min) * scale:1: log2(fq_max) * scale) / scale);
        sds = 10.^( log10(.39) + .5 * (log10(cfs)));
        sds = (sds) * sqrt(2.);
    
        sds=round(sds*sds_round_val)/sds_round_val;
        cfs=round(cfs*cfs_round_val)/cfs_round_val;
    
        if nargin<2
            sds=sds;
            cfs=cfs;
        else 
            index=(cfs<f_high) & (cfs>f_low);
            cfs=cfs(index);
            sds=sds(index);
        end 

    end
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINE GAUSSIAN FILTERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function k=gaussian_filter(obj,X,rate,center,sd)
        % TODO - description

        N = size(X,2);
        d = 1./rate;
        freq = [0:ceil(N/2-1), ceil(-(N)/2):1:-1]/(d*N);
        k = exp((-(abs(freq) - center).^2)/(2 * (sd^2)));
        k = k/norm(k);
    end 


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HILBERT TRANSFORM (HIGH GAMMA EXTRACTION USING GAUSSIAN FILTERS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Xh_out,X_fft_h_out]=hilbert_transform(obj,X,rate,varargin)
        % https://github.com/BouchardLab/ecogVIS/blob/master/ecogvis/signal_processing/hilbert_transform.py
        % Based on : Computing the Discrete-Time Analytic Signal via FFTS.
        % Lawrence  Marple,  Jr. IEEE 1999
        %
        % Arguments
        % ---------
        % X : ndarray (n_channels, n_time)
        %   Input data, dimensions
        % rate : float
        %   Number of samples per second.
        % filters : filter or cell of filters (optional)
        %   One or more bandpass filters
        % 
        % Returns
        % -------
        % Xh : ndarray, complex
        %   Bandpassed analytic signal
        %

        switch nargin 
        case 4
            filters=varargin{1};
            phase='None';
            X_fft_h='None';
        case 5
            filters=varargin{1};
            phase=varargin{2};
            X_fft_h='None';
        case 6
            filters=varargin{1};
            phase=varargin{2};
            X_fft_h=varargin{3};
        otherwise 
            filters='None';
            phase='None';
            X_fft_h='None';        
        end 
            
        if ~iscell(filters)
            filters = {filters};
        end 

        time = size(X,2);
        d=1./rate;
        freq=[0:ceil(time/2-1), ceil(-(time)/2):1:-1]/(d*time);
        Xh=cell(length(filters),1);
        %Xh = zeros((len(filters),) + X.shape, dtype=np.complex)

        % compute the one sided analytic signal in frequency domain 
        if strcmp(X_fft_h, 'None')
            % Heavyside filter
            h = zeros(1,length(freq));
            h(freq > 0) = 2.;
            h(1) = 1.;
            % X_fft_h is the single sided spectrum of X
            X_fft_h = fft(X) .* h;

            if ~strcmp(phase,'None')
                X_fft_h=X_fft.*phase;
            end
        end
            
        for ii=1:size(filters,2)
            f=filters{ii};
            if strcmp(f,'None')
                Xh{ii} = ifft(X_fft_h);
            else
                % filter the signal in different frequency bands 
                f = f / max(f);
                    
                Xh{ii} = ifft(X_fft_h .* f);
            end 
        end

        if size(Xh,1) == 1
            Xh_out= Xh{1};
            X_fft_h_out=X_fft_h;
        else 
            Xh_out= Xh;
            X_fft_h_out=X_fft_h;
        end 
        
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REMOVE OUTLIERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [envelopes_rjct,outliers_idx,outlier_prcnt,ignored_channels]=envelope_outliers(obj,envelopes,ops,varargin)
        % Detects outliers in signal envelopes. The difference between the median
        % and the Qth percentile (default: 90%) of the envelope distribution for
        % each electrode is measured. Envelopes are considered outliers if they
        % fall a certain number of units above this value.
        % 
        % 2016-08-15 - Created, Sam NH
        % 
        % 2019-01-21 - Cosmetic changes, Sam NH
        %
        % 2021-11-29 - update to fit with Evlab ecog pipeline 

        n_channels = size(envelopes,2);

        % make threshold parameters electrode specific if not already
        if isvector(ops.threshold)
            ops.threshold = repmat(ops.threshold(:), 1, n_channels);
        else
            assert(size(ops.threshold) == n_channels);
        end

        % % trim first and last bit of envelope
        % samples_to_remove = obj.sample_freq * 1;
        % envelopes = envelopes(samples_to_remove:size(envelopes,1)-samples_to_remove,:);

        % normalize by 90% of the distribution
        s = quantile(envelopes, ops.percentile);
        Z = envelopes ./ repmat(s, size(envelopes,1),1);
        clear s;

        outliers = Z > ops.threshold;
        envelopes_rjct = nan*envelopes;
        outlier_prcnt = nan*zeros(size(envelopes,2),1);
        outliers_idx = {};
        
        ignored_channels = [];

        % loop through channels
        for q=1:n_channels
            outlie = outliers(:,q);
            outlier_prcnt(q) = sum(outlie)./numel(outlie)*1e2;
            
            envl = envelopes(:,q);
            envl_rejct = envl;
            
            rise = find(diff(outlie)==1).';
            fall = find(diff(outlie)==-1).';

            try 
                assert(numel(rise)==numel(fall));
                process = true;
            catch err 
                ignored_channels = [ignored_channels, q];
                process = false;
            end 

            outlie_idx ={ };
            if numel(rise)>0 & process
                for r=rise % go through each outlier
                    
                    f = fall(find(fall>r,1,'first'));
                    
                    % find window to look at for interpolation
                    interp_win = r:f;
                    interp_val = envl(interp_win).';
                    lookup_window = max([r-ops.buffer,1]):min([(f+ops.buffer),numel(outlie)]);
                    main_sig = envl(lookup_window);
                    [C,ia,ib] = intersect(lookup_window,interp_win);

                    aux = main_sig;
                    fixed_sig = main_sig;
                    
                    % interpolate
                    aux(ia) = [];
                    intrp_sig = interp1(setdiff(lookup_window,interp_win,'stable'),aux,interp_win,ops.interpMethod);
                    fixed_sig(ia) = intrp_sig;

                    % add interpolated values to signal
                    envl_rejct(lookup_window) = fixed_sig;
                    outlie_idx = [outlie_idx,[interp_win;interp_val]];
                end 
            end

            envelopes_rjct(:,q) = envl_rejct;
            outliers_idx = [outliers_idx;{outlie_idx}];
            
            fprintf(1,'.');

        end 

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMBINE DATA FILES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function combine_data_files(obj)
        % Combines trial info from multiple data files
        %
        % This should be the last step in preprocessing after trial info is 
        % added if there is more than one data file.

        if isa(obj.trial_timing{1},'table') % should be a cell if multiple files
            fprintf(1,'No need to combine data files!\n')
            return
        end

        % output
        trial_timing_raw = [];
        trial_timing_dec = [];
        condition        = [];
        session          = [];
        events_table     = [];

        fprintf(1, '\n> Combining trial info from data files \n');
        fprintf(1,'[');

        for i=1:length(obj.stitch_index) % number of separate data files with signal
            
            samples_to_add_raw = obj.for_preproc.stitch_index_raw(i)-1;
            samples_to_add_dec = obj.for_preproc.stitch_index_dec(i)-1;
            trial_timing_raw_  = obj.for_preproc.trial_timing_raw{i}; % all trial timing tables for one file
            trial_timing_dec_  = obj.for_preproc.trial_timing_dec{i}; % all trial timing tables for one file

            % go through all trials in file
            for j=1:size(trial_timing_raw_,1)
                trial_timing_raw__ = trial_timing_raw_{j}; % one trial timing table
                trial_timing_dec__ = trial_timing_dec_{j}; % one trial timing table

                trial_timing_raw__(:,'start') = table(trial_timing_raw__.start + samples_to_add_raw); 
                trial_timing_raw__(:,'end')   = table(trial_timing_raw__.end + samples_to_add_raw); 
                trial_timing_dec__(:,'start') = table(trial_timing_dec__.start + samples_to_add_dec); 
                trial_timing_dec__(:,'end')   = table(trial_timing_dec__.end + samples_to_add_dec); 

                trial_timing_raw_{j} = trial_timing_raw__;
                trial_timing_dec_{j} = trial_timing_dec__;

                fprintf(1,'.');

            end

            trial_timing_raw = [trial_timing_raw; trial_timing_raw_];
            trial_timing_dec = [trial_timing_dec; trial_timing_dec_];

            % no need to update timing information
            condition    = [condition; obj.condition{i}];
            session      = [session; obj.session{i}];
            events_table = [events_table; obj.events_table{i}]; % based on within trial timing

        end

        fprintf(1,'] done\n');
            
        % add both versions of trial timing (raw & ds) 
        obj.for_preproc.trial_timing_raw = trial_timing_raw;
        obj.for_preproc.trial_timing_dec = trial_timing_dec;

        % update obj.trial_timing (trial timing for current version of signal)
        if obj.sample_freq == obj.for_preproc.decimation_freq % downsampling in preproc
            obj.trial_timing = trial_timing_dec;
        else % no downsampling in preproc
            obj.trial_timing = trial_timing_raw;
        end

        obj.condition    = condition;
        obj.session      = session;
        obj.events_table = events_table;

    end 


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHANNEL MAPPING FROM ANATOMICAL FILES TO OBJECT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [chan_mapping,chan_labels]=channel_mapping_anatomical(obj,vera_mat)
        % helper function for mapping channels in datafiles to channels in anatomical files
        % TODO - modify to be less manual
    
        fprintf(1, '\n> Adding anatomical files to object\n');

        % load channel labels
        chan_labels = obj.elec_ch_label;
        chan_types = obj.elec_ch_type;

        if strcmp(obj.subject,'BJH006')
            % map localized channels to channels in object
            reduced_localized_names = cellfun(@(x) split(x,'-'),vera_mat.electrodeNames,'UniformOutput',false);
            localized_names = cellfun(@(x) strcat(x{1},'_',extract(x{end},digitsPattern)),reduced_localized_names,'UniformOutput',false); 
            localized_names = cellfun(@(x) x{1}, localized_names,'UniformOutput',false);
            chan_mapping = cellfun(@(x) find(strcmp(x,localized_names)),chan_labels,'UniformOutput',false);
            chan_labels = cellfun(@(x) localized_names(x),chan_mapping,'UniformOutput',false);

        elseif strcmp(obj.subject,'SLCH002')
            % map localized channels to channels in object
            reduced_localized_names = cellfun(@(x) split(x,'^'),vera_mat.electrodeNames,'UniformOutput',false);
            localized_names = cellfun(@(x) strcat(x{1},'_',extract(x{2},digitsPattern)),reduced_localized_names,'UniformOutput',false); 
            localized_names = cellfun(@(x) x{1}, localized_names,'UniformOutput',false);
            chan_mapping = cellfun(@(x) find(strcmp(x,localized_names)),chan_labels,'UniformOutput',false);
            chan_labels = cellfun(@(x) localized_names(x),chan_mapping,'UniformOutput',false);

        elseif contains(obj.subject,'BJH') 
            % map localized channels to channels in object
            first_thing = cellfun(@(x) x(1:3),vera_mat.electrodeNames,'UniformOutput',false);
            second_thing = cellfun(@(x) extract(x,digitsPattern),vera_mat.electrodeNames,'UniformOutput',false); 
            second_thing = cellfun(@(x) x{end},second_thing,'UniformOutput',false);
            localized_names = cellfun(@(x,y) strcat(x([1,3]),'_',y),first_thing,second_thing,'UniformOutput',false);
            % localized_names = cellfun(@(x) x{1}, localized_names,'UniformOutput',false);
            if strcmp(obj.subject,'BJH008')
                for i=1:length(localized_names)
                    if contains(localized_names{i},'ER')
                        localized_names(i) = {['O' localized_names{i}(3:end)]};
                    elseif contains(localized_names{i},'FR')
                        localized_names(i) = {['P' localized_names{i}(3:end)]};
                    elseif contains(localized_names{i},'FR')
                        localized_names(i) = {['P' localized_names{i}(3:end)]};
                    else
                        localized_names(i) = {localized_names{i}([1,3:end])};
                    end
                end
            end
            chan_mapping = cellfun(@(x) find(strcmp(x,localized_names)),chan_labels,'UniformOutput',false);
            chan_labels = cellfun(@(x) localized_names(x),chan_mapping,'UniformOutput',false);

        else % AMC and MCJ subjects
            num_chans = sum(contains(chan_types,'ecog') | strcmp('seeg',chan_types));
            % assert(num_chans==length(vera_mat.electrodeNames),'Mapping may not be correct!');
            chan_mapping = num2cell(1:num_chans)';
            chan_labels = chan_labels(1:num_chans);
            % add empy cells to end of mapping/labels
            empty_to_add = cell(size(obj.elec_data,1)-num_chans,1);
            chan_mapping = [chan_mapping; empty_to_add];
            chan_labels = [chan_labels; empty_to_add];
        end
        
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT DATA STRUCTURES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function output_data_structures(obj)
        % Outputs a MATLAB structure and a python Xarray of object.
        %
        %

    end 


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT CHANNELS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_channels(obj,varargin)
        % TODO - description

        p = inputParser();
        addRequired(p,'signal');
        addRequired(p,'channel_labels');
        addRequired(p,'clean_channels');
        addRequired(p,'valid_channels');
        addParameter(p,'stitch_index',[1]);
        addParameter(p,'t_len',60); % size of viewing window in seconds
        addParameter(p,'sample_freq',1200);
        addParameter(p,'downsample',false);
        addParameter(p,'decimation_freq',300)
        addParameter(p,'plotIEDs',false);
        addParameter(p,'chanIEDs',[]);
        addParameter(p,'posIEDs',[]);
        addParameter(p,'save',false); 
        parse(p, varargin{:});
        ops = p.Results;
        
        curr_sample_freq = ops.sample_freq;
        D = ops.signal;
        

        % ------------------------------
        % FORMAT & NORMALIZE SIGNAL
        % ------------------------------
        x_norm_cell = [];
        for k=1:length(ops.stitch_index) % number of separate data files with signal

            if k == length(ops.stitch_index) % signal for file stops at end of matrix
                stop = size(D,1);
            else % signal for file stops before stitch index of next file
                stop = ops.stitch_index(k+1)-1;
            end 

            D_ = D(ops.stitch_index(k):stop,:);

            % formatting
            x_cell = mat2cell(D_',ones(1,size(D_,2))); % make signal from each electrode a cell

            % downsampling
            if ops.downsample
                decimation_factor = ops.sample_freq/ops.decimation_freq;
                x_cell = cellfun(@(x) downsample(x,decimation_factor),x_cell,'uni',false);
                curr_sample_freq = ops.decimation_freq;
            end

            % normalizing
            min_max = mean(cell2mat(cellfun(@(y) prctile(y,[3 97]),x_cell,'UniformOutput',false))); % mean 5th and 95th percentiles of ALL electrodes
            x_norm_cell_ = cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false); % normalize
            x_norm_cell_ = arrayfun(@(x) x_norm_cell_{x}*ops.valid_channels(x),1:size(x_norm_cell_,1),'uni',false)'; % set signal of noisy channels to 0
            
            if length(ops.stitch_index) > 1
                x_norm_cell = [x_norm_cell, x_norm_cell_];
            else 
                x_norm_cell = x_norm_cell_;
            end

        end

        % combine separate x_norm_cell columns
        if length(ops.stitch_index) > 1
            for k=1:length(obj.stitch_index)-1
                x_norm_cell(:,k+1) = arrayfun(@(x) {[x_norm_cell{x,k}, x_norm_cell{x,k+1}]},[1:size(x_norm_cell,1)])';
            end
            x_norm_cell = x_norm_cell(:,length(ops.stitch_index));
        end
        assert(size(x_norm_cell,2)==1,'x_norm_cell not in the correct format');


        % ------------
        % PLOT SIGNAL
        % ------------
        t_length = ops.t_len*curr_sample_freq;

        col_inf=inferno(floor(.8*size(x_norm_cell,1)));
        col_vir=viridis(floor(.8*size(x_norm_cell,1)));
        colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];

        close all
        figure(1);
        clf;
        set(gcf,'position',[31,1,1700,900]); % 1713, 1010
        ax = axes('position',[.05,.1,.93,.88]);
        hold on
        time_stamps = [1:size(x_norm_cell{1},2)]/curr_sample_freq;
        hold on
        H = arrayfun(@(x) plot(time_stamps,x_norm_cell{x}+x,'color',colors(x,:),'tag',sprintf('ch %d, tag %s',x,ops.channel_labels{x})),[1:size(x_norm_cell,1)]);
           

        if ops.plotIEDs
            spike_chan = ops.chanIEDs;
            [spike_chan_sort,sort_idx] = sort(spike_chan);
            spike_times_sort = ops.posIEDs(sort_idx);
            
            yval = cell2mat(arrayfun(@(x) x_norm_cell{spike_chan_sort(x)}(floor(spike_times_sort(x)*200))+spike_chan_sort(x),1:length(spike_chan_sort),'uni',false));
            H1 = scatter(spike_times_sort,yval,50,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 .7 .7],'LineWidth',2,'MarkerFaceAlpha',.5);
        end
            
        set(gcf,'doublebuffer','on');
        set(ax,'ytick',[1:size(x_norm_cell,1)]);
        set(ax,'yticklabel','');
        arrayfun(@(x) text(0,x,num2str(x),'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle'),ops.clean_channels);
           
        set(ax,'ylim',[0,size(x_norm_cell,1)+4]);
        ax.XAxis.TickLength = [0.005,0.01];
        ax.YAxis.TickLength = [0.005,0.01];
        set(ax,'xlim',[0 ops.t_len]);
        pos = get(ax,'position');
        Newpos = [pos(1) pos(2)-0.1 pos(3) 0.05];
        xmax=max(time_stamps);
        S = ['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(ops.t_len) '])'];
            
        h = uicontrol('style','slider','units','normalized','position',Newpos,'callback',S,'min',0,'max',xmax-ops.t_len);
        datacursormode on
            
        waitfor(findobj('type','figure','number',1));

    end


    %% TRIAL DATA METHODS


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAKE TRIAL DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function make_trials(obj)
        % TODO - description

        if ~isempty(obj.bip_elec_data)
            trial_keys = {'key','string','elec_data','bip_elec_data'};
        else % no bipolar data
            trial_keys = {'key','string','elec_data'};
        end

        fprintf(1, '\n> Cutting signal into trial data ... \n');
        fprintf(1,'[');

        % go through each trial timing table
        for k = 1:size(obj.trial_timing)
            trial_time_tbl = obj.trial_timing{k};

            trial_elec_data = arrayfun(@(x) obj.elec_data(:,trial_time_tbl(x,:).start:trial_time_tbl(x,:).end),1:size(trial_time_tbl,1),'uni',false)';
        
            if ~isempty(obj.bip_elec_data)
                trial_bip_elec_data = arrayfun(@(x) obj.bip_elec_data(:,trial_time_tbl(x,:).start:trial_time_tbl(x,:).end),1:size(trial_time_tbl,1),'uni',false)';
                obj.trial_data{k,1} = table(trial_time_tbl.key,trial_time_tbl.string,trial_elec_data,trial_bip_elec_data,'VariableNames',trial_keys);
            else % no bipolar data 
                obj.trial_data{k,1} = table(trial_time_tbl.key,trial_time_tbl.string,trial_elec_data,'VariableNames',trial_keys);
            end

            fprintf(1,'.');

        end

        fprintf(1,'] done\n');

    end

    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT CONDITION TRIAL DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cond_data=get_cond_resp(obj,condition,varargin)
        % Extracts trial data of given condition. 

        p = inputParser();
        addParameter(p,'keep_trials',[]);
        parse(p, varargin{:});
        ops = p.Results;

        if isempty(obj.trial_data)
            obj.make_trials();
        end
        
        if ~isempty(ops.keep_trials)
            assert(length(obj.condition)==length(ops.keep_trials),'Logical array of trials to keep has incorrect dimensions');
            cond_id = find(cell2mat(arrayfun(@(x) (strcmp(obj.condition{x},condition) && ops.keep_trials(x)),1:length(obj.condition),'UniformOutput',false)));
        else
            cond_id = find(cell2mat(arrayfun(@(x) strcmp(obj.condition{x},condition),1:length(obj.condition),'UniformOutput',false)));
        end

        cond_data = obj.trial_data(cond_id);
        
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AVERAGE TRIAL DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function input_d_ave=get_average(obj,input_d)
        % Averages trial data signal. Returns cell array of trial timing
        % tables with the average signal value for each event (e.g., word)
        % and signal type (e.g., unipolar, bipolar)
        %
        % NOTE - it does not operate on obj.trial_data.
        input_d_ave = input_d;

        % go through trials of input
        for k=1:size(input_d,1)
            B = input_d{k};

            [keys,strings,values] = obj.get_columns(B);
            
            values_ave = varfun( @(x) cellfun(@(y) nanmean(y,2), x,'uni',false), values,'OutputFormat','table');
            values_ave.Properties.VariableNames = values.Properties.VariableNames;
            input_d_ave{k} = [keys,strings,values_ave];

        end

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT EVENTS FROM TRIAL DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function output_d=get_value(obj,input_d,varargin)
        % Extracts events with label ops.key (e.g., 'word') from trial
        % data. Returns cell array of trial timing tables with only 
        % events of interest.
        %
        % NOTE - this does not operate directly on obj.trial_data.

        p =inputParser();
        addParameter(p, 'key', 'word');
        addParameter(p, 'type', 'match'); % match or contain
        parse(p, varargin{:});
        ops = p.Results;

        if strcmp(ops.type,'match')
            func = @(x,y) ismember(x,y);
        else,
            func = @(x,y) contains(x,y);
        end

        output_d = input_d;

        % go through trials of input
        for k = 1:size(input_d,1)
            B = input_d{k};
            output_d{k} = B(func(B.key,ops.key),:);
        end

    end
 

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT COLUMNS FROM TRIAL DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [keys,strings,values]=get_columns(obj,B)
        % Pulls out columns of provided trial data table.

        keys = B(:,ismember(B.Properties.VariableNames,'key'));
        strings = B(:,ismember(B.Properties.VariableNames,'string'));
        values = B(:,~ismember(B.Properties.VariableNames,{'key','string'}));
        
    end

        
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMBINE TRIAL CONDITIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function output_d=combine_trial_cond(obj,input_d)
        % Combines separate trial tables (e.g., for a condition) into one 
        % table.
        
        all_keys = cellfun(@(x) x.key,input_d,'uni',false);

        % make sure the keys are the same for all trials 
        [X,Y] = ndgrid(1:numel(all_keys));
        Z = tril(true(numel(all_keys)),-1);
        assert(all(arrayfun(@(x,y) isequal(all_keys{x},all_keys{y}),X(Z),Y(Z))));
            
        % function to combine trial data tables
        func = @(x) cell2mat(reshape(vertcat(x),1,[]));

        output_d = table();

        % go through keys in trial data table
        for k=1:numel(all_keys{1})
            
            each_key = all_keys{1}{k};
            temp = (cellfun(@(x) x(ismember(x.key,each_key),:),input_d,'uni',false));
            cond_tbl = vertcat(temp{:});
            [~,strings,values] = obj.get_columns(cond_tbl);
            values_comb = cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
            temp_table = cell2table(horzcat(each_key,{strings.string},values_comb),'VariableNames',cond_tbl.Properties.VariableNames);
            
            output_d = [output_d; temp_table];

        end 
            
    end

        
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AVERAGE SIGNAL FOR CONDITION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [output_tbl,cond_table]=get_ave_cond_trial(obj,varargin)
        % Averages signal for given condition across all words.    

        p = inputParser();
        addParameter(p,'words',1:12);
        addParameter(p,'condition',[]);
        addParameter(p,'keep_trials',[])
        parse(p, varargin{:});
        ops = p.Results;
            
        func = @(x) cell2mat(permute(x,[3,2,1])); % format : electrode*trial_id*state/word
            
        % get trials condition
        if ops.condition
            condition_flag = ops.condition;
            cond_data = obj.get_cond_resp(condition_flag,'keep_trials',ops.keep_trials);
        else 
            if isempty(obj.trial_data)
                obj.make_trials();
            end
            condition_flag = 'all';
            cond_data = obj.trial_data;
        end
        
        cond_data_ave = obj.get_average(cond_data);
        word_data = obj.get_value(cond_data_ave,'key','word','type','contain');
            
        B = obj.combine_trial_cond(word_data);
        [keys,strings,values] = obj.get_columns(B);
        values_comb = cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
        cond_table = cell2table(horzcat(condition_flag,{strings.string},values_comb),'VariableNames',B.Properties.VariableNames);

        % select how many words are selected : default 1:12 and average over words 
        % create the comparision between the two condition
        func_1 = @(x) x{1}(:,:,ops.words); % format : electrode*trial_id
        func_2 = @(x) nanmean(x,3); % format : electrode*trial_id
        
        B = cond_table;
        [keys,strings,values] = obj.get_columns(B);
        condition_ave = cellfun(@(X) func_2(func_1(values.(X))),values.Properties.VariableNames,'uni',false);
        
        output_tbl = cell2table(horzcat(condition_flag,strings.string,condition_ave),'VariableNames',B.Properties.VariableNames);
        
    end 

    
end


end

