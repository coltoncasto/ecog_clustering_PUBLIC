function validateClusters(varargin)
    % TODO - make more general to work with other values of k

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'method','kmedoids');
    addParameter(p,'distance','correlation');
    addParameter(p,'k',3);
    addParameter(p,'cutoff',0.9);
    addParameter(p,'repeats',100);
    addParameter(p,'srate',60);
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs',true);
    addParameter(p,'omitIterations',100); % how many iterations run per each specific number of channels we omit
    addParameter(p,'omitChannelsBy',5); % number of channels in between each dropout
    addParameter(p,'useWandJ',true); % MITSWJNTask only
    addParameter(p,'threshold',[]);
    addParameter(p,'averageTrials',true); % use thresholds from shuffling average condition responses
    parse(p, varargin{:});
    ops = p.Results;

    % --- INITIALIZE --- %

    % paths
    [~,SAVE_PATH] = initialize(ops.saveName);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'validation' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'validation' filesep];
    if ~exist(PLOT_PATH,'dir'), mkdir(PLOT_PATH); end
    if ~exist(PLOT_PATH2,'dir'), mkdir(PLOT_PATH2); end

    % file naming 
    DATA_PATH = [SAVE_PATH 'data' filesep];
    if ops.useLangElecs, elecType = 'langElecs'; else, elecType = 'nonLangElecs'; end
    if ops.useWandJ, cond_string = 'SWJN'; else, cond_string = 'SN'; end
    if strcmp(ops.experiment,'both')
        expt_string = 'bothMITSWJNTaskandMITLangloc';
    else
        expt_string = ops.experiment;
    end
    if ~ops.averageTrials
        trial_string = '_all_trials';
    else
        trial_string = '';
    end

    % load in data matrix 
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering.mat']); % all_X
    X = all_X;

    % load in data labels
    all_X_table = readtable([SAVE_PATH 'clustering' filesep expt_string '_' elecType '_cluster_assignments.csv']);

    % load cluster assignments 
    eval(strcat("assignments=all_X_table.k",num2str(ops.k),";"));
    IDX_all = assignments;

    % get average cluster response
    C = zeros(ops.k,size(X,2));
    for kk=1:ops.k
        C(kk,:) = mean(X(IDX_all==kk,:),1);
    end
    C_all = C; % don't want to get overwritten

    % other params
    k = ops.k;
    srate = ops.srate; 
    t = 0:1/srate:((size(X,2)-1)/srate);
    if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
        nConds = 4;
    else % only S and N
        nConds = 2;
    end

    % -----------------------------------------
    % --- CLUSTER LANG WITHOUT N ELECTRODES --- %
    % -----------------------------------------

    fprintf(1,'\n > CLUSTERING LANG CHANS WITH OMMITTED CHANNELS \n');
    
    totNchans = size(X,1);

    rng('shuffle')
    omitNchannels = [0 flip(totNchans-k:-1*ops.omitChannelsBy:1)];
    dists_same_cluster = nan(k,ops.omitIterations,length(omitNchannels)); % k x iters x nChansRemoved
    dists_btwn_cluster = nan(sum(nonzeros(tril(ones(k),-1))),ops.omitIterations,length(omitNchannels)); % k x iters x nChansRemoved
    iN=0;
    for omitN = omitNchannels
        iN=iN+1;
        ticomit = tic;
        
        for iter = 1:ops.omitIterations
            
            % randomly select a group of totNchans-omitN channels
            currNchannels = totNchans - omitN;
            Xn = X(randperm(totNchans, currNchannels),:);

            % --- CLUSTERING --- %
            if strcmp(ops.method,'kmeans')
                [IDX,Cn,SUMD,D] = kmeans(Xn,k,'Replicates',ops.repeats,'Distance',ops.distance);
            elseif strcmp(ops.method,'kmedoids')
                [IDX,Cn,SUMD,D] = kmedoids(Xn,k,'Replicates',ops.repeats,'Distance',ops.distance);
            end
        
            % ---- EVALUATE RESULTS ---- %
            [IDX,Cn,SUMD,D] = reorderClusters(Xn,C_all,IDX,Cn,SUMD,D,ops.distance);
                
            % calculate mean of cluster (in case using kmedoids)
            Cn_means = zeros(k,size(Xn,2));
            for kk=1:k
                Cn_means(kk,:) = mean(Xn(IDX==kk,:),1);
            end

            % calculate distance between this solution and the full solution (by cluster)
            % calculate distance between this solution and the full solution (by cluster)
            dists_curr = pdist2(C_all,Cn_means,ops.distance);
            dists_same_cluster(:,iter,iN) = ones(size(diag(dists_curr)))-diag(dists_curr);
            dists_curr = pdist2(Cn_means,Cn_means,ops.distance);
            dists_btwn_cluster(:,iter,iN) = ones(size(nonzeros(tril(dists_curr,-1))'))-nonzeros(tril(dists_curr,-1))';

        end

        fprintf(1,['Done ' num2str(ops.omitIterations) ' iterations\tomiting ' num2str(omitN) ' electrodes in ' num2str(toc(ticomit)) ' \n'])
        fprintf(1,['\t\t\tremaining ' num2str(length(omitNchannels)-iN) ' out of ' num2str(length(omitNchannels)) ' omitions \n'])

    end
  

    % --- PLOT --- %

    h = ERPfigure; set(h,'Position',[0 0 800 600],'visible',ops.isPlotVisible);
    colors = hsv(k);
    if k==3
        sig_y_values = [0.2 0.15 0.1];
        colors_btwn = [mean(colors([1 2],:),1); mean(colors([1 3],:),1); mean(colors([2 3],:),1);];
    elseif k==2 
        colors_btwn = [mean(colors([1 2],:),1)];
    end
    
    % save files with omitted channel analysis
    VAL_PATH = [SAVE_PATH 'validation' filesep];
    if ~exist(VAL_PATH,'dir'), mkdir(VAL_PATH); end
    dists_to_full = dists_same_cluster;
    filename = [VAL_PATH ops.experiment '_' elecType '_K=' num2str(k) '_dropout_validation_dists_to_full' trial_string '.mat'];
    save(filename,'dists_to_full','-v7.3');
    dists_within_omitted = dists_btwn_cluster;
    filename = [VAL_PATH ops.experiment '_' elecType '_K=' num2str(k) '_dropout_validation_dists_within_omitted' trial_string '.mat'];
    save(filename,'dists_within_omitted','-v7.3');

    % same cluster to same cluster (how well do the reduced solutions match the full solution)
    subplot(2,1,1); 
    dists_same_mean = squeeze(mean(dists_same_cluster,2)); % k x nChannelsRemoved
    dists_same_sem = squeeze(std(dists_same_cluster,[],2)) / sqrt(size(dists_same_cluster,2));
    for kk=1:k
        plot(omitNchannels,dists_same_mean(kk,:),'color',colors(kk,:),'linewidth',2); hold on;
        patch([omitNchannels fliplr(omitNchannels)], [dists_same_mean(kk,:)+dists_same_sem(kk,:) fliplr(dists_same_mean(kk,:)-dists_same_sem(kk,:))],...
                colors(kk,:),'linestyle','none','FaceAlpha',0.2);
        if ops.threshold
            sig_idxs = dists_same_mean(kk,:) > ops.threshold;
            y_values = ones(1,sum(sig_idxs))*sig_y_values(kk);
            x_values = omitNchannels(sig_idxs);
            scatter(x_values,y_values,[],colors(kk,:),'*');
        end
    end
    if ops.threshold
        plot([0 max(omitNchannels)],[ops.threshold ops.threshold],":k",'linewidth',2);
    end
    if k==3
        if ops.threshold
            legend({'C1 vs C1','','','C2 vs C2','','','C3 vs C3','','',''},'Location','eastoutside','box','off');
        else 
            legend({'C1 vs C1','','C2 vs C2','','C3 vs C3',''},'Location','eastoutside','box','off');
        end
    elseif k==2 
        legend({'C1 vs C1','','C2 vs C2',''},'Location','eastoutside','box','off');
    end
    set(gca,'fontsize',12,'box','off');
    xlabel('Number of omitted electrodes');
    ylabel('Correlation with full solution','fontsize',12);
    ylim([0 1]); xlim([0 max(omitNchannels)]);

    % different cluster to differnt cluster (how well does the *relationship* between clusters in the reduced solutions match the full solution)
    subplot(2,1,2); 
    dists_btwn_mean = squeeze(mean(dists_btwn_cluster,2)); % nRelations x nChannelsRemoved
    dists_btwn_sem = squeeze(std(dists_btwn_cluster,[],2)) / sqrt(size(dists_btwn_cluster,2));
    if k==2
        dists_btwn_mean = dists_btwn_mean';
        dists_btwn_sem = dists_btwn_sem';
    end
    for kk=1:size(dists_btwn_mean,1)
        plot(omitNchannels,dists_btwn_mean(kk,:),'--','color',colors_btwn(kk,:),'linewidth',2); hold on;
        patch([omitNchannels fliplr(omitNchannels)], [dists_btwn_mean(kk,:)+dists_btwn_sem(kk,:) fliplr(dists_btwn_mean(kk,:)-dists_btwn_sem(kk,:))],...
                colors_btwn(kk,:),'linestyle','none','FaceAlpha',0.2);
    end
    if k==3
        legend({'C1 vs C2','','C1 vs C3','','C2 vs C3',''},'Location','eastoutside','box','off');
    elseif k==2
        legend({'C1 vs C2',''},'Location','eastoutside','box','off');
    end
    set(gca,'fontsize',12,'box','off');
    xlabel('Number of omitted electrodes');
    ylabel('Correlation within ommitted solution','fontsize',12);
    ylim([0 1]); xlim([0 max(omitNchannels)]);


    sgtitle({'Cluster Validation Holding Out Channels', [num2str(ops.omitIterations) ' iterations']})
    
    % save png 
    saveas(gcf,[PLOT_PATH2 ops.experiment '_' elecType '_K=' num2str(k) '_dropout_validation' trial_string '.png'],'png');

    % save pdf 
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH ops.experiment '_' elecType '_K=' num2str(k) '_dropout_validation' trial_string '.pdf'],'pdf')

end