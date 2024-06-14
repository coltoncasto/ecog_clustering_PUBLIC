function [values,values_sem]=get_timecourses(obj,varargin)
    % TODO - description

    p = inputParser();
    addParameter(p,'words',[]);
    addParameter(p,'condition',[]);
    addParameter(p,'signalType','unipolar') % 'bipolar'
    addParameter(p,'useLangElecs',true)
    addParameter(p,'split',[]) % empty 'odd' or 'even'
    parse(p, varargin{:});
    ops = p.Results;

    if strcmp(ops.signalType,'unipolar')
        data = obj.elec_data_zs_dec;
        sig_elecs = logical(obj.s_vs_n_sig.elec_data_zs_dec{1,1});
    elseif strcmp(ops.signalType,'bipolar')
        data = obj.bip_elec_data_zs_dec;
        sig_elecs = logical(obj.s_vs_n_sig.bip_elec_data_zs_dec{1,1});
    end

    starts_n_stops = cellfun(@(x) [x.start(2), x.end(length(ops.words)+1)],obj.trial_timing_dec(:,1),'UniformOutput',false);

    if ops.useLangElecs
        values = cellfun(@(x) data(sig_elecs,x(1):x(2)),starts_n_stops,'UniformOutput',false);
    else
        values = cellfun(@(x) data((~sig_elecs & obj.elec_ch_valid),x(1):x(2)),starts_n_stops,'UniformOutput',false);
    end

    values = values(strcmp(obj.trial_type,ops.condition));

    % split the data if specified
    if strcmp(ops.split,'odd')
        values = values(1:2:length(values)); 
    elseif strcmp(ops.split,'even')
        values = values(2:2:length(values));
    end

    % cut trials to same number of samples (some variabi)
    trial_lengths = cell2mat(cellfun(@(x) size(x,2),values,'UniformOutput',false));
    min_length = min(trial_lengths);
    for i=1:length(values)
        curr_trial = values{i,1};
        if size(curr_trial,2)>min_length
            values{i,1} = curr_trial(:,1:min_length);
        end
    end
    
    % average over trials
    values_done = zeros(size(values{1,1},1),size(values{1,1},2));
    values_sem = zeros(size(values{1,1},1),size(values{1,1},2));
    for i=1:size(values{1,1},1)
        tmp = zeros(length(values),size(values{1,1},2));
        for j=1:length(values)
            tmp(j,:) = values{j,1}(i,:);
        end
        values_done(i,:) = mean(tmp,1);
        values_sem(i,:) = std(tmp,1)/sqrt(size(tmp,1));
    end
    values = values_done;

end