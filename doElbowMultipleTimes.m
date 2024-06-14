function doElbowMultipleTimes(varargin)

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); 
    addParameter(p,'method','kmedoids');
    addParameter(p,'distance','correlation');
    addParameter(p,'doElbow',false);
    addParameter(p,'minK',1);
    addParameter(p,'maxK',10);
    addParameter(p,'maxElbowK',10);
    addParameter(p,'cutoff',0.9);
    addParameter(p,'repeats',100);
    addParameter(p,'srate',60);
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs',true);
    addParameter(p,'useWandJ',false); % MITSWJNTask only
    addParameter(p,'split',[]); % 'odd' or 'even'
    addParameter(p,'colors',[]); % when specified will only work with one value of k
    addParameter(p,'filterByReliability',true);
    addParameter(p,'reliabThresh',0:0.05:0.4);
    addParameter(p,'words',8);
    addParameter(p,'signalType','unipolar'); % or 'bipolar'
    parse(p, varargin{:});
    ops = p.Results;

    % --- INITIALIZE ---

    % paths
    [CLUSTER_PATH,SAVE_PATH] = initialize(ops.saveName);
    if ops.filterByReliability
        PLOT_PATH = [SAVE_PATH 'plots' filesep 'clustering' filesep 'filteredByReliability' filesep];
        PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'clustering' filesep 'filteredByReliability' filesep];
        reliab_dir = [filesep 'filteredByReliability'];
    else 
        PLOT_PATH = [SAVE_PATH 'plots' filesep 'clustering' filesep];
        PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'clustering' filesep];
        reliab_dir = '';
    end
    if ~exist(PLOT_PATH,'dir')
        mkdir(PLOT_PATH);
    end
    if ~exist(PLOT_PATH2,'dir')
        mkdir(PLOT_PATH2);
    end

    % file naming 
    DATA_PATH = [SAVE_PATH 'data' filesep];
    if ops.useLangElecs, elecType = 'langElecs'; else, elecType = 'nonLangElecs'; end
    if ops.split, split_string = ['_' ops.split]; else, split_string = ''; end
    if strcmp(ops.experiment,'both')
        expt_string = 'bothMITSWJNTaskandMITLangloc';
    else
        expt_string = ops.experiment;
    end
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

    % other params
    if strcmp(ops.experiment,'MITLangloc')
        ylims1 = [0.2 1];
        ylims2 = [0 0.4];
    elseif strcmp(ops.experiment,'MITSWJNTask')
        ylims1 = [0.4 1];
        ylims2 = [0 0.2];
    end
    ops.reliabThresh = [-1, ops.reliabThresh(2:end)]; % MAY NEED TO ADJUST FOR FUTURE
    max_alpha = 0.9; min_alpha = 0.1;
    alphas = min_alpha:(max_alpha-min_alpha)/(length(ops.reliabThresh)-1):max_alpha;

    reliab_string = ['_filtered_by_reliability_threshold_' strrep(sprintf('%0.2f',min(ops.reliabThresh(1))),'0.','') ...
                     '_' strrep(sprintf('%0.2f',max(ops.reliabThresh(end))),'0.','')];

    % --- FINDING OPTIMAL K --- %

    for i=1:length(ops.reliabThresh)
        tic;
        % load in data matrix
        load([DATA_PATH expt_string '_' elecType '_data_for_clustering' split_string '.mat']); % all_X
        X = all_X;

        % filter by reliability if desired
        RELIAB_PATH = [CLUSTER_PATH 'output/_reliability/reliability/'];
        load([RELIAB_PATH expt_string cond_string '_' elecType signal_string '_reliability.mat']); % corrs

        X = X(corrs>ops.reliabThresh(i),:);

        % search for optimal k
        if strcmp(ops.method,'kmeans')
            [~,~,~,~,~,d{i},Var{i},~] = kmeans_opt(X,ops.maxElbowK,ops.cutoff,ops.repeats,ops.distance);
        elseif strcmp(ops.method,'kmedoids')
            [~,~,~,~,~,d{i},Var{i},~] = kmedoids_opt(X,ops.maxElbowK,ops.cutoff,ops.repeats,ops.distance);
        end
        fprintf(1,['\nDone testing ' num2str(ops.maxElbowK) ' values of K with ' num2str(ops.repeats) ' repeats in ' num2str(toc) ' s\n'])

    end

    % --- PLOTTING --- % 
         
    % initialize plot
    h = ERPfigure; set(h,'Position',[10 10 1100 550],'visible',ops.isPlotVisible) 
    
    % plotting results of search for optimal k
    c = [0.8500 0.3250 0.0980];

    % distance
    subplot(1,2,1); hold on;
    for i=1:length(ops.reliabThresh)
        d_curr = d{i};
        patchline(1:length(d_curr),d_curr/d_curr(1),'linestyle','-','edgecolor',c,'edgealpha',alphas(i),'linewidth',4) 
        scatter(1:length(d_curr),d_curr/d_curr(1),100,c,'filled','MarkerFaceAlpha',alphas(i),'MarkerEdgeAlpha',alphas(i));
        set(gca,'fontsize',20,'box','off'); title('Variance Unexplained','fontsize',28);
        xlabel({' ','Number of Clusters (K)',' '},'fontsize',22); xlim([1 length(d_curr)]); xticks(1:length(d_curr));
        if ops.filterByReliability, ylim(ylims1); end
    end
    
    % change in distance (TODO)
    subplot(1,2,2); hold on;
    for i=length(ops.reliabThresh):-1:1
        d_curr = d{i};
        Var_curr = Var{i};
        bar(1:length(d_curr)-1,Var_curr/d_curr(1),'facecolor',c,'facealpha',alphas(i));
        set(gca,'fontsize',20,'box','off'); title('\Delta Variance Unexplained','fontsize',28);
        xlabel({' ','K',' '},'fontsize',22); xlim([0 ops.maxElbowK-1]); set(gca,'XTick',0:ops.maxElbowK-1,'XTickLabel',strsplit(num2str(1:ops.maxElbowK)));
        if ops.filterByReliability, ylim(ylims2); end
    end

    hs = sgtitle({' ','Search for Optimal K',' '}); set(hs,'fontsize',32,'fontweight','bold');
    
    % save png
    saveas(gcf,[PLOT_PATH2 expt_string '_' elecType '_search_for_optimalK_K=' num2str(ops.maxElbowK)  split_string reliab_string '.png']);

    % save pdf
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH expt_string '_' elecType '_search_for_optimalK_K=' num2str(ops.maxElbowK) split_string reliab_string '.pdf'],'pdf');


end