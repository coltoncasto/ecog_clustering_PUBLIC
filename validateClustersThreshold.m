function validateClustersThreshold(varargin)
    % assumes that clusters are in MATCHED order across solutions being compared

    p = inputParser();
    addRequired(p,'solution'); % folder name of solution 1
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc'
    addParameter(p,'elecType','langElecs');
    addParameter(p,'k',3); % k value for solutions to compare
    addParameter(p,'iterations',1000);
    addParameter(p,'method','kmedoids'); % can be 'assign'
    addParameter(p,'distance','correlation');
    addParameter(p,'repeats',100);
    addParameter(p,'useWandJ',true); % MITSWJNTask only
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'averageTrials',true); % use average condition responses
    parse(p, varargin{:});
    ops = p.Results;

    % --------------------
    % --- INITIALIZE --- %
    % --------------------

    % solution 1
    [~,SAVE_PATH] = initialize(ops.solution);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'clustering' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'clustering' filesep];
    if ~exist(PLOT_PATH,'dir'), mkdir(PLOT_PATH); end
    if ~exist(PLOT_PATH2,'dir'), mkdir(PLOT_PATH2); end

    % loading solution 1
    expt = split(ops.solution,'_'); expt = expt{1}; % TODO - doesn't work with 'both'
    load([SAVE_PATH 'data' filesep expt '_' ops.elecType '_data_for_clustering.mat']); X = all_X;
    all_X_table = readtable([SAVE_PATH 'clustering' filesep expt '_' ops.elecType '_cluster_assignments.csv']);
    % get average cluster response
    eval(strcat("assignments=all_X_table.k",num2str(ops.k),";"));
    IDX = assignments; C = zeros(ops.k,size(X,2));
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


    % ------------------------------------------
    % --- SIGNIFICANCE W/ TRIALS SCRAMBLED --- %
    % ------------------------------------------

    fprintf(1,'\nCALCULATING SIGNIFICANCE W/ TRIALS SCRAMBLED \n');

    if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
        nConds = 4;
    else % only S and N
        nConds = 2;
    end

    % select matrix to perform permutations on
    if ~ops.averageTrials
        X_use = X_all_trials;
    else
        X_use = X;
    end

    totNchans = size(X_use,1);
    dists_same_cluster_scr = nan(ops.k,ops.iterations); % k x iters
        
    tic
    for iter = 1:ops.iterations
        
        % randomly shuffle trials
        Xn = nan(size(X_use));
        cond_starts = 1:size(X_use,2)/nConds:size(X_use,2);
        cond_ends = size(X_use,2)/nConds:size(X_use,2)/nConds:size(X_use,2)+1;
        for i=1:nConds
            Xn(:,cond_starts(i):cond_ends(i)) = X_use(randperm(totNchans),cond_starts(i):cond_ends(i));
        end

        % take mean of shuffled trials if using all trials
        if ~ops.averageTrials
            Xn_new = nan(size(X));
            for elec=1:size(elec_idxs,1)
                Xn_new(elec,:) = mean(Xn(elec_idxs(elec,1):elec_idxs(elec,2),:),1);
            end
            Xn = Xn_new;
        end

        % --- CLUSTERING --- %
        if strcmp(ops.method,'assign')
            IDX = zeros(size(Xn,1),1);
            for kk=1:size(Xn,1)
                dists = cell2mat(arrayfun(@(x) pdist2(Xn(kk,:),C1(x,:),ops.distance),1:ops.k,'UniformOutput',false));
                [~,IDX(kk,1)] = min(dists);
            end
        elseif strcmp(ops.method,'kmeans')
            [IDX,Cn,SUMD,D] = kmeans(Xn,ops.k,'Replicates',ops.repeats,'Distance',ops.distance);
        elseif strcmp(ops.method,'kmedoids')
            [IDX,Cn,SUMD,D] = kmedoids(Xn,ops.k,'Replicates',ops.repeats,'Distance',ops.distance);
        end
            
        % ---- EVALUATE RESULTS ---- %
        if ~strcmp(ops.method,'assign')
            [IDX,Cn,SUMD,D] = reorderClusters(Xn,C,IDX,Cn,SUMD,D,ops.distance);
        end
                    
        % calculate mean of cluster (in case using kmedoids)
        Cn_means = zeros(ops.k,size(Xn,2));
        for kk=1:ops.k
            Cn_means(kk,:) = mean(Xn(IDX==kk,:),1);
        end

        % calculate distance between this solution and the full solution (by cluster)
        dists_curr = pdist2(C,Cn_means,ops.distance);
        dists_same_cluster_scr(:,iter) = ones(size(diag(dists_curr)))-diag(dists_curr);

        if mod(iter,50)==0
            fprintf(1,'Done with 50 iterations in %0.4f seconds\n',toc);
            tic;
        end

    end

    if ~ops.averageTrials
        trial_string = '_all_trials';
    else
        trial_string = '';
    end

    % save null distribution
    COMP_PATH = [SAVE_PATH 'clustering' filesep];
    if ~exist(COMP_PATH,'dir'), mkdir(COMP_PATH); end
    filename = [COMP_PATH expt '_' ops.elecType '_K=' num2str(ops.k) '_null_distribution' trial_string '.mat'];
    null = dists_same_cluster_scr;
    save(filename,'null','-v7.3');


    % ------------------------------------------------
    % --- PLOTTING NULL DISTRIBUTIONS BY CLUSTER --- %
    % ------------------------------------------------

    if ops.k==3
        pos = [0 0 400 600];
    elseif ops.k==5
        pos = [0 0 400 900];
    elseif ops.k==2
        pos = [0 0 400 450];
    end

    k = ops.k;

    close all
    h = ERPfigure; set(h,'Position',pos,'visible',ops.isPlotVisible); hold on
    set(gca,'fontsize',16); box off;

    colors = hsv(ops.k);

    % histograms
    for i=1:ops.k

        h(i) = subplot(ops.k,1,i); hold on
        histogram(null(i,:),0:.05:1,'FaceColor',colors(i,:),'FaceAlpha',.7,...
                  'EdgeColor',colors(i,:),'EdgeAlpha',.7,'Normalization','probability');
        ylim([0 0.6]); set(gca,'ytick',0:.2:.6);
        
        prc95 = prctile(null(i,:),95)
        plot([prc95 prc95],[0 1],'--','Color',colors(i,:),'linewidth',3);
        text(0.1,0.45,['sig thresh (< 0.05) = ' num2str(round(prc95,2))],'Color',colors(i,:),'fontsize',14);

        legend(['Cluster #' num2str(i)],'Location','northwest','fontsize',14,'box','off');
        if i==ops.k
            xlabel({' ','Correlation'},'fontsize',18,'fontweight','bold');
        end
        if (i==3) && (ops.k==5)
            ylabel({' ','Probability',' '},'fontsize',18,'fontweight','bold');
        elseif (i==2) && (ops.k==3)
            ylabel({' ','Probability',' '},'fontsize',18,'fontweight','bold');
        end
    end
    
    % save png
    saveas(gcf,[PLOT_PATH2 expt '_' ops.elecType '_K=' num2str(ops.k) '_null_distribution' trial_string '.png']);

    % save pdf 
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH expt '_' ops.elecType '_K=' num2str(ops.k) '_null_distribution' trial_string '.pdf'],'pdf')

end