function decodeCondition(varargin)
    % assumes that clusters are in MATCHED order across solutions being compared
    % TODO - make work with different file naming convention for experiments

    p = inputParser();
    addRequired(p,'solution'); % folder name of solution 1
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc'
    addParameter(p,'elecType','langElecs');
    addParameter(p,'k',3); % k value for solutions to compare
    addParameter(p,'singleK',[]); % cluster to perform decoding for
    addParameter(p,'srate',60);
    addParameter(p,'windowSize',100); % ms
    addParameter(p,'slideBy',100); % ms 
    addParameter(p,'iterations',20); % number of decoders to train on real data to find best performance
    addParameter(p,'kfold',10); % number of folds in cross validation
    addParameter(p,'decoders','SW'); % binary classifiers to build, see line 126
    addParameter(p,'singleDecoder',true); % whether comparing multiple conditions or just one (if single, assume convention is 'cond1cond2')
    addParameter(p,'useWandJ',true); % MITSWJNTask only
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'averageTrials',true); % use average condition responses
    addParameter(p,'doCummulative',true);
    addParameter(p,'permutations',1000); % number of perumations for constructing null distribution
    addParameter(p,'threshold',0.05); % significance threshold for permutation test
    addParameter(p,'doPlot',false);
    addParameter(p,'loadData',false); % whether to load saved data instead of training classifiers again
    addParameter(p,'words',8);
    addParameter(p,'filterByReliability',false);
    addParameter(p,'reliabThresh',0.1);
    addParameter(p,'signalType','unipolar'); % or 'bipolar'
    parse(p, varargin{:});
    ops = p.Results % print for log files

    % --------------------
    % --- INITIALIZE --- %
    % --------------------

    % solution 1
    [CLUSTER_PATH,SAVE_PATH] = initialize(ops.solution);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'decoding' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'decoding' filesep];
    if ~exist(PLOT_PATH,'dir'), mkdir(PLOT_PATH); end
    if ~exist(PLOT_PATH2,'dir'), mkdir(PLOT_PATH2); end

    % loading solution 1
    expt = split(ops.solution,'_'); expt = expt{1}; % TODO - doesn't work with 'both'
    load([SAVE_PATH 'data' filesep expt '_' ops.elecType '_data_for_clustering.mat']); X = all_X;
    all_X_table = readtable([SAVE_PATH 'clustering' filesep expt '_' ops.elecType '_cluster_assignments.csv']);
    
    % filter by reliability if desired
    if ops.filterByReliability
        if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
            cond_string = '_SWJN';
        elseif strcmp(ops.experiment,'MITSWJNTask') % only 2 conds
            cond_string = '_SN';
        else
            cond_string = '';
        end
        if strcmp(ops.signalType,'bipolar')
            signal_string = '_bipolar';
        else
            signal_string = '';
        end

        RELIAB_PATH = [CLUSTER_PATH 'output/_reliability/reliability/'];
        load([RELIAB_PATH expt cond_string '_' ops.elecType signal_string '_reliability.mat']); % corrs

        X = X(corrs>ops.reliabThresh,:);
        reliab_string = ['_filtered_by_reliability_threshold_' strrep(sprintf('%0.2f',ops.reliabThresh),'0.','')];
    else
        reliab_string = '';
    end

    % get average cluster response
    eval(strcat("assignments=all_X_table.k",num2str(ops.k),";"));
    if ops.filterByReliability
        IDX = assignments(corrs>ops.reliabThresh); 
    else
        IDX = assignments; 
    end     
    C = zeros(ops.k,size(X,2));
    for kk=1:ops.k, C(kk,:) = mean(X(IDX==kk,:),1); end

    % load all trials if desired
    if ~ops.averageTrials
        load([SAVE_PATH 'data' filesep expt '_' ops.elecType '_data_for_clustering_all_trials.mat']); X_all_trials = all_X;
        all_X_table_all_trials = readtable([SAVE_PATH 'data' filesep expt '_' ops.elecType '_labels_for_clustering_all_trials.csv']);
        % find electrode idxs for averaging after shuffling trials
        elec_idxs = zeros(size(X,1),2); % starts in column 1, ends in column 2
        idx_start = 1; elec_num = 1;
        prev_name = all_X_table_all_trials.channel_name{1};
        for i=2:size(all_X_table_all_trials,1)
            curr_name = all_X_table_all_trials.channel_name{i};
            if ~strcmp(prev_name,curr_name) || (i==size(all_X_table_all_trials,1))
                elec_idxs(elec_num,2) = i;
                elec_idxs(elec_num,1) = idx_start;
                idx_start = i+1;
                elec_num = elec_num + 1; 
            end
            prev_name = curr_name;
        end
    end

    % file naming
    if ~ops.averageTrials
        trial_string = '_all_trials';
    else
        trial_string = '';
    end
    if ops.doCummulative
        cummulative_string = '_training_cummulative';
    else
        cummulative_string = '_training_by_timepoint';
    end    
    sig_string = split(num2str(ops.threshold),'.');
    sig_string = sig_string{2};


    % ---------------------------------------
    % --- DECODE CONDTIIONS OF INTEREST --- %
    % ---------------------------------------

    % number of conditions and binary classifiers to constuct
    if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
        nConds = 4;
        cond_order = {'S','W','J','N'};
        if ops.singleDecoder
            decoders = {[ops.decoders(1) ' ' ops.decoders(2)]};
            shapes = {'-'};
            titles = {strrep(decoders{1},' ',' vs ')};
        else % want to compare multiple conditions
            if strcmp(ops.decoders,'timescales')
                decoders = {'S W','W N'};
                shapes = {'-','--'};
                titles = {'S vs W','W vs N'};
            elseif strcmp(ops.decoders,'timescales_revision')
                decoders = {'S W','W N','J N'};
                shapes = {'-','-','-'};
                titles = {'S vs W','W vs N','J vs N'};
            elseif strcmp(ops.decoders,'all')
                decoders = {'S W','S J','S N','W J','W N','J N'};
                shapes = {'-o','-+','-*','-.','-<','->'}; % haven't tested
                titles = {'S vs W','','S vs J','','S vs N','','W vs J','','W vs N','','J vs N',''};
            elseif strcmp(ops.decoders,'test')
                decoders = {'S W'};
                shapes = {'-'};
                titles = {'S vs W'};
            end
        end
    else % only S and N
        nConds = 2;
        cond_order = {'S','N'};
        decoders = {'S N'};
        titles = {'S vs. N'};
    end
    legend_labels = cell(1,ops.k*2);
    kk = 1;
    for k=1:2:ops.k*2
        legend_labels{k} = ['Cluster #' num2str(kk)];
        legend_labels{k+1} = '';
        kk = kk+1;
    end

    % select matrix to perform permutations on
    if ~ops.averageTrials
        X_use = X_all_trials;
    else
        X_use = X;
    end
    cond_starts = 1:size(X_use,2)/nConds:size(X_use,2);
    cond_ends = size(X_use,2)/nConds:size(X_use,2)/nConds:size(X_use,2)+1;

    % variables for binning of neural signal
    trial_len = size(X,2)/nConds;
    sliding_samples = ops.slideBy/1000 / (1/ops.srate);
    window_samples = ops.windowSize/1000 / (1/ops.srate);
    window_starts = 1:sliding_samples:trial_len-window_samples+1;

    % initialize empty 4d matrix for accuracies
    accuracies = nan(ops.kfold,length(window_starts),ops.k,length(decoders));
    null_accuracies = nan(ops.permutations,ops.kfold,length(window_starts),ops.k,length(decoders));
    tsum_null_thresholds = nan(ops.k,length(decoders));

    % values of k to train decoders for
    if ~isempty(ops.singleK)
        minK = ops.singleK;
        maxK = ops.singleK;
    else
        minK = 1;
        maxK = ops.k;
    end

    if ~ops.loadData
        % train separate decoder for each cluster
        parpool;
        fprintf(1,'\nCONSTRUCTING DECODERS ... \n');
        for k=minK:maxK

            % extract cluster 
            if ops.averageTrials
                Xc_use = X_use(IDX==k,:);
            else
                Xc_use = [];
                cluster_chans = find(IDX==k);
                for c=1:length(cluster_chans)
                    chan_idx = cluster_chans(c);
                    Xc_use = [Xc_use; X_use(elec_idxs(chan_idx,1):elec_idxs(chan_idx,2),:)];
                end
            end

            % train separate decoder for each binary classification of interest
            tic
            for d=1:length(decoders)

                % conditions
                conds = split(decoders{d},' ');
                cond_nums = cell2mat(cellfun(@(x) find(strcmp(cond_order,x)),conds,'UniformOutput',false));
                
                % build data to provide to decoder
                Xn = nan(size(Xc_use,1)*length(conds),size(Xc_use,2)/nConds); % nTrials*nConds x nSamples 1 trial
                labels = []; % y labels to use for classification
                idxs = 1:size(Xc_use,1);
                for ii=1:length(conds)
                    labels = [labels; repmat(conds(ii),size(Xc_use,1),1)];
                    i = cond_nums(ii);
                    Xn(idxs,:) = Xc_use(:,cond_starts(i):cond_ends(i));
                    idxs = idxs+size(Xc_use,1);
                end

                % bin neural signal (look ops.slideBy to the left and right)
                % if ops.slideBy == ops.windowSize then windows are non-overlapping
                Xn_binned = nan(size(Xn,1),length(window_starts));
                for w=1:length(window_starts)
                    start = window_starts(w);
                    Xn_binned(:,w) = mean(Xn(:,start:start+window_samples-1),2);
                end

                % --- REAL DATA --- 
                % train separate decoder for each time point (can be cummulative or not)
                fprintf(1,'Training %s classifier on real data from Cluster %d ...\n\n',titles{d},k)
                for w=1:size(Xn_binned,2)
                    if ops.doCummulative
                        X_decode = Xn_binned(:,1:w);
                    else
                        X_decode = Xn_binned(:,w);
                    end
                    best_accs = zeros(ops.kfold,1); % for storing accuracies from best model
                    for iter=1:ops.iterations % take best decoder from iterations
                        model = fitclinear(X_decode,labels,'Learner','logistic','Solver','sgd','CrossVal','on','KFold',ops.kfold);
                        [label,score]=kfoldPredict(model);
                        fold_assignments = model.ModelParameters.Generator.UseObsForIter;
                        curr_accs = nan(ops.kfold,1);
                        for fold=1:size(fold_assignments,2)
                            idxs = ~fold_assignments(:,fold);
                            curr_accs(fold,1) = sum(cell2mat(cellfun(@(x,y) strcmp(x,y),label(idxs),labels(idxs),...
                                                            'UniformOutput',false)))/length(labels(idxs));
                        end
                        % replace accuracies if current model is better than the best so far
                        if mean(curr_accs,1)>mean(best_accs,1)
                            accuracies(:,w,k,d) = curr_accs;
                            best_accs = curr_accs;
                        end
                    end
                end

                % --- NULL DATA --- 
                % train separate decoder for each time point (can be cummulative or not)
                for w=1:size(Xn_binned,2)
                    tic
                    fprintf(1,'Performing permutations for timepoint %d/%d ...',w,size(Xn_binned,2))
                    if ops.doCummulative
                        X_decode = Xn_binned(:,1:w);
                    else
                        X_decode = Xn_binned(:,w);
                    end
                    % shuffle label to construct null distribtution
                    % b = ProgressBar(ops.permutations,'IsParallel', true,'Title', 'Parallel 1');
                    % b.setup([], [], []);
                    parfor perm=1:ops.permutations
                        best_accs = zeros(ops.kfold,1); % for storing accuracies from best model
                        for iter=1:ops.iterations % take best decoder from iterations
                            model = fitclinear(X_decode,labels(randperm(length(labels))),'Learner','logistic','Solver','sgd','CrossVal','on','KFold',ops.kfold);
                            [label,score]=kfoldPredict(model);
                            fold_assignments = model.ModelParameters.Generator.UseObsForIter;
                            curr_accs = nan(ops.kfold,1);
                            for fold=1:size(fold_assignments,2)
                                idxs = ~fold_assignments(:,fold);
                                curr_accs(fold,1) = sum(cell2mat(cellfun(@(x,y) strcmp(x,y),label(idxs),labels(idxs),...
                                                                'UniformOutput',false)))/length(labels(idxs));
                            end
                            % replace accuracies if current model is better than the best so far
                            if mean(curr_accs,1)>mean(best_accs,1)
                                null_accuracies(perm,:,w,k,d) = curr_accs;
                                best_accs = curr_accs;
                            end
                        end
                        % updateParallel();
                    end
                    % b.release();
                    fprintf(1,'\tDone in %0.4f seconds\n',toc)
                end

                % --- COMPUTE NULL DISTRIBUTION ---
                tsum_distribution = [];
                for perm=1:ops.permutations
                    tvals_curr = squeeze(null_accuracies(perm,:,:,k,d));
                    tvals_curr = (mean(tvals_curr,1)-0.5)./(std(tvals_curr,0,1)/sqrt(ops.kfold));
                    tvals_curr_sig = tvals_curr>tinv(1-ops.threshold,ops.kfold);
                    
                    % sum t-values per cluster
                    dx = diff(tvals_curr_sig);
                    isnew = (dx~=0);
                    idx = cumsum([1 isnew]); % array where each cluster takes value from 1-nClusters
                    tsum_curr = splitapply(@sum,tvals_curr,idx); % t-sum for all clusters (including non sig clusters)
                    % only want sum from sig clusters
                    if tvals_curr_sig(1)==1 % when sig cluster is first
                        tsum_curr = tsum_curr(1:2:max(idx));
                    else % when sig cluster is NOT first
                        tsum_curr = tsum_curr(2:2:max(idx));
                    end
                    tsum_distribution = [tsum_distribution, tsum_curr];
                end

                % store null distribution threshold for computing significance later
                tsum_null_thresholds(k,d) = prctile(tsum_distribution,(1-ops.threshold)*100);

                % save decoding accuracies
                DEC_PATH = [SAVE_PATH 'decoding' filesep];
                if ~exist(DEC_PATH,'dir'), mkdir(DEC_PATH); end
                filename = [DEC_PATH expt '_' ops.elecType '_K=' num2str(ops.k) '_condition_decoding' trial_string cummulative_string '_' strrep(decoders{d},' ','_') '_sig' sig_string '_k=' num2str(k) reliab_string '.mat'];
                eval(['accuracies_' strrep(decoders{d},' ','_') '_k' num2str(k) '=squeeze(accuracies(:,:,k,d));']);
                save(filename,['accuracies_' strrep(decoders{d},' ','_') '_k' num2str(k)],'-v7.3');
                filename = [DEC_PATH expt '_' ops.elecType '_K=' num2str(ops.k) '_condition_decoding' trial_string cummulative_string '_' strrep(decoders{d},' ','_') '_sig' sig_string '_k=' num2str(k) reliab_string '_null_distribution.mat'];
                eval(['null_accuracies_' strrep(decoders{d},' ','_') '_k' num2str(k) '=squeeze(null_accuracies(:,:,:,k,d));'])
                save(filename,['null_accuracies_' strrep(decoders{d},' ','_') '_k' num2str(k)],'-v7.3');
                filename = [DEC_PATH expt '_' ops.elecType '_K=' num2str(ops.k) '_condition_decoding' trial_string cummulative_string '_' strrep(decoders{d},' ','_') '_sig' sig_string '_k=' num2str(k) reliab_string '_null_distribution_thresholds.mat'];
                eval(['tsum_null_thresholds_' strrep(decoders{d},' ','_') '_k' num2str(k) '=tsum_null_thresholds(k,d);'])
                save(filename,['tsum_null_thresholds_' strrep(decoders{d},' ','_') '_k' num2str(k)],'-v7.3');

                fprintf(1,'Done with Cluster %d %s decoder in %0.4f seconds\n',k,titles{d},toc);
                tic;

            end

        end
    else
        fprintf(1,'\nLOADING DECODERS ... \n');
        accuracies = nan(ops.kfold,length(window_starts),ops.k,length(decoders));
        null_accuracies = nan(ops.permutations,ops.kfold,length(window_starts),ops.k,length(decoders));
        tsum_null_thresholds = nan(ops.k,length(decoders));
        for k=minK:maxK
            if ops.filterByReliability && (ops.reliabThresh==0.3) && (k==3)
                continue
            end
            for d=1:length(decoders)
                DEC_PATH = [SAVE_PATH 'decoding' filesep];
                filename = [DEC_PATH expt '_' ops.elecType '_K=' num2str(ops.k) '_condition_decoding' trial_string cummulative_string '_' strrep(decoders{d},' ','_') '_sig' sig_string '_k=' num2str(k) reliab_string '.mat'];
                load(filename); % accuracies_{decoders{d}}_k{k}
                eval(['accuracies(:,:,k,d) = accuracies_' strrep(decoders{d},' ','_') '_k' num2str(k) ';']);
                filename = [DEC_PATH expt '_' ops.elecType '_K=' num2str(ops.k) '_condition_decoding' trial_string cummulative_string '_' strrep(decoders{d},' ','_') '_sig' sig_string '_k=' num2str(k) reliab_string '_null_distribution.mat'];
                load(filename); % null_accuracies_{decoders{d}}_k{k}
                eval(['null_accuracies(:,:,:,k,d) = null_accuracies_' strrep(decoders{d},' ','_') '_k' num2str(k) ';']);
                filename = [DEC_PATH expt '_' ops.elecType '_K=' num2str(ops.k) '_condition_decoding' trial_string cummulative_string '_' strrep(decoders{d},' ','_') '_sig' sig_string '_k=' num2str(k) reliab_string '_null_distribution_thresholds.mat'];
                load(filename); % tsum_null_thresholds_{decoders{d}}_k{k}
                eval(['tsum_null_thresholds(k,d) = tsum_null_thresholds_' strrep(decoders{d},' ','_') '_k' num2str(k) ';']);
            end
        end
    end

    
    % --------------------------------------
    % --- PLOTTING ACCURACY BY CLUSTER --- %
    % --------------------------------------

    if ops.doPlot
        pos = [0 0 600 400*length(decoders)];

        close all
        h = ERPfigure; set(h,'Position',pos,'visible',ops.isPlotVisible); hold on
        set(gca,'fontsize',16); box off;

        colors = hsv(ops.k);
        x = (ops.slideBy/1000):(ops.slideBy/1000):(trial_len*(1/ops.srate)); % s

        % histograms
        for d=1:length(decoders)

            h(d) = subplot(length(decoders),1,d); hold on

            for k=1:ops.k
                if ops.filterByReliability && (ops.reliabThresh==0.3) && (k==3)
                    continue
                end
                avg_acc = squeeze(mean(accuracies(:,:,k,d),1));
                sem_acc = squeeze(std(accuracies(:,:,k,d),0,1)) / sqrt(size(accuracies(:,:,k,d),1));
                plot(x,avg_acc,shapes{d},'color',colors(k,:),'linewidth',2);
                patch([x fliplr(x)], [avg_acc-sem_acc fliplr(avg_acc+sem_acc)],...
                        colors(k,:),'linestyle','none','FaceAlpha',0.2);
            end
            ylim([0.3 1]); set(gca,'ytick',0.4:0.2:1); xlim([0 x(end)]);
            plot([x(1) x(end)],[0.5 0.5],':k','linewidth',2); % chance

            for iii=1:ops.words-1 % repeat for last condition
                word_onset = ((trial_len/ops.words)/ops.srate)*iii;
                plot([word_onset word_onset],[0.3 1],'color','#D3D3D3','linewidth',1);
            end
            
            if d==length(decoders)
                xlabel({' ','Time (s)'},'fontsize',18,'fontweight','bold');
            end
            ylabel({' ','Classification Accuracy',' '},'fontsize',18,'fontweight','bold');

            % evaluate significance
            height = 0.38; % y value of significance stars in plot
            for k=1:ops.k
                if ops.filterByReliability && (ops.reliabThresh==0.3) && (k==3)
                    continue
                end
                tvals_curr = squeeze(accuracies(:,:,k,d));
                tvals_curr = (mean(tvals_curr,1)-0.5)./(std(tvals_curr,0,1)/sqrt(ops.kfold));
                tvals_curr_sig = tvals_curr>tinv(1-ops.threshold,ops.kfold);
                    
                % sum t-values per cluster
                dx = diff(tvals_curr_sig);
                isnew = (dx~=0);
                if ~all(isnew==0) % check that there is a significant cluster
                    idx = cumsum([1 isnew]); % array the size of tvals_curr where each cluster takes value from 1-nClusters
                    tsum_curr = splitapply(@sum,tvals_curr,idx); % t-sum for all clusters (including non sig clusters)
                    % only want info from sig clusters (not '0' clusters)
                    if tvals_curr_sig(1)==1 % when sig cluster is first
                        sig_clusts = 1:2:max(idx);
                    else % when sig cluster is NOT first
                        sig_clusts = 2:2:max(idx);
                    end
                    tsum_curr = tsum_curr(sig_clusts);
                    idx_idx = 1; idx_sig = [];
                    for clust=sig_clusts
                        idx_sig{idx_idx} = find(idx==clust);
                        idx_idx = idx_idx+1;
                    end

                    % plot signifiance stars 
                    for clust=1:length(tsum_curr)
                        if tsum_curr(clust)>tsum_null_thresholds(k,d)
                            idx_clust = idx_sig{clust};
                            s = scatter(x(idx_clust),repmat(height,1,length(idx_clust)),10,colors(k,:),'*');
                        end
                    end
                    height = height-0.02;
                else % no significant cluster
                    continue;
                end
            end

            title(titles{d},'fontsize',22);
            legend(legend_labels,'Location','northwest','fontsize',14,'box','off');

        end
        
        % save png
        saveas(gcf,[PLOT_PATH2 expt '_' ops.elecType '_K=' num2str(ops.k) '_condition_decoding_' ops.decoders trial_string cummulative_string '_sig' sig_string reliab_string '.png']);

        % save pdf 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH expt '_' ops.elecType '_K=' num2str(ops.k) '_condition_decoding_' ops.decoders trial_string cummulative_string '_sig' sig_string reliab_string '.pdf'],'pdf')
    
    end
    
end