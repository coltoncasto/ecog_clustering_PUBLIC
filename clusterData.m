function clusterData(varargin)

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'method','kmedoids');
    addParameter(p,'distance','correlation');
    addParameter(p,'doElbow',false);
    addParameter(p,'minK',3);
    addParameter(p,'maxK',3);
    addParameter(p,'maxElbowK',10);
    addParameter(p,'cutoff',0.9);
    addParameter(p,'repeats',100);
    addParameter(p,'srate',60);
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs',true);
    addParameter(p,'useWandJ',false); % MITSWJNTask only
    addParameter(p,'split',[]); % 'odd' or 'even'
    addParameter(p,'template',[]); % filename of template to use in output/_templates/ folder
    addParameter(p,'colors',[]); % when specified will only work with one value of k
    addParameter(p,'filterByReliability',false);
    addParameter(p,'reliabThresh',0.1);
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
    
    % load in data matrix
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering' split_string '.mat']); % all_X
    X = all_X;

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
        load([RELIAB_PATH expt_string cond_string '_' elecType signal_string '_reliability.mat']); % corrs

        X = X(corrs>ops.reliabThresh,:);
        reliab_string = ['_filtered_by_reliability_threshold_' strrep(sprintf('%0.2f',ops.reliabThresh),'0.','')];

        if strcmp(ops.experiment,'MITLangloc')
            ylims1 = [0.2 1];
            ylims2 = [0 0.4];
        elseif strcmp(ops.experiment,'MITSWJNTask')
            ylims1 = [0.4 1];
            ylims2 = [0 0.2];
        end
    else
        reliab_string = '';
    end

    % load in data labels
    save_table_filename = [SAVE_PATH 'clustering' reliab_dir filesep expt_string '_' elecType '_cluster_assignments' split_string reliab_string '.csv'];
    if exist(save_table_filename,'file')
        all_X_table = readtable(save_table_filename);
    else % start from scratch
        all_X_table = readtable([DATA_PATH expt_string '_' elecType '_labels_for_clustering' split_string '.csv']);
    end
    if ops.filterByReliability & ~exist(save_table_filename,'file')
        all_X_table = all_X_table(corrs>ops.reliabThresh,:);
    end

    % other params
    srate = ops.srate; 
    t = (1:size(X,2))/srate;
    if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
        nConds = 4;
        heatmap_title = {' ','Concatenated Timecourses',' ', ...
            'S               W                J              N'};
        cluster_offset = -1;
    else % only S and N
        nConds = 2;
        heatmap_title = {' ','Concatenated Timecourses',' ', ...
            'S                                  N'};
        cluster_offset = -0.4;
    end
    t_per_cond = t(length(t)/nConds);
    

    % ------------------
    % --- CLUSTERING ---
    % ------------------


    
    % --- FINDING OPTIMAL K ---

    if ops.doElbow
    
        tic;
        if strcmp(ops.method,'kmeans');
            [IDX,C,SUMD,D,K,d,Var,PC] = kmeans_opt(X,ops.maxElbowK,ops.cutoff,ops.repeats,ops.distance);
        elseif strcmp(ops.method,'kmedoids');
            [IDX,C,SUMD,D,K,d,Var,PC] = kmedoids_opt(X,ops.maxElbowK,ops.cutoff,ops.repeats,ops.distance);
        end
        fprintf(1,['\nDone testing ' num2str(ops.maxElbowK) ' values of K with ' num2str(ops.repeats) ' repeats in ' num2str(toc) ' s\n'])

        % plotting results of search for optimal k
        h = ERPfigure; set(h,'Position',[10 10 1100 550],'visible',ops.isPlotVisible) 
        c = [0.8500 0.3250 0.0980];

        % distance
        subplot(1,2,1); 
        plot(1:length(d),d/d(1),'-o','color',c,'markersize',12,'markeredgecolor',c,'linewidth',3); 
        set(gca,'fontsize',20,'box','off'); title('Variance Unexplained','fontsize',28);
        xlabel({' ','Number of Clusters (K)',' '},'fontsize',22); xlim([1 length(d)]); xticks(1:length(d));
        if ops.filterByReliability, ylim(ylims1); end

        % change in distance
        subplot(1,2,2); 
        bar(1:length(d)-1,Var/d(1),'facecolor',c,'facealpha',0.75);
        set(gca,'fontsize',20,'box','off'); title('\Delta Variance Unexplained','fontsize',28);
        xlabel({' ','K',' '},'fontsize',22); xlim([0 ops.maxElbowK-1]); set(gca,'XTick',0:ops.maxElbowK-1,'XTickLabel',strsplit(num2str(1:ops.maxElbowK)));
        if ops.filterByReliability, ylim(ylims2); end

        hs = sgtitle({' ','Search for Optimal K',' '}); set(hs,'fontsize',32,'fontweight','bold');

        % save png
        saveas(gcf,[PLOT_PATH2 expt_string '_' elecType '_search_for_optimalK_K=' num2str(ops.maxElbowK)  split_string reliab_string '.png']);

        % save pdf
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH expt_string '_' elecType '_search_for_optimalK_K=' num2str(ops.maxElbowK) split_string reliab_string '.pdf'],'pdf');

    end


    % --- CLUSTERING WITH ALL VALUES OF K ---

    % parameters for all plots
    n = 12; % number of columns in subplot

    % stitch idxs for concatenation
    length_trial = size(X,2)/nConds;
    stitch_idxs = length_trial:length_trial:size(X,2);
    length_word = length_trial/ops.words;

    % stitch idxs for subjects
    unique_subs = unique(all_X_table.subject);
    subject_idxs = zeros(size(unique_subs,1),1);
    for i=1:length(unique_subs)
        idxs = find(cell2mat(cellfun(@(x) strcmp(x,unique_subs{i}),all_X_table.subject,'UniformOutput',false)));
        subject_idxs(i) = idxs(1);
    end
    subject_idxs = subject_idxs(2:end);

    % xticks 
    xlocslabels = repmat({'','1','2','3'},1,nConds);
    xlocs = repmat([0 1 2 3],1,nConds);
    u = 1;
    for i=5:4:length(xlocs)-1
        xlocs(:,i:i+3) = xlocs(:,i:i+3)+(u*t_per_cond);
        u = u+1; 
    end

    % ylimits
    if strcmp(ops.distance,'correlation') & strcmp(ops.method,'kmeans')
        ylims = [-0.2 0.2];
    elseif contains(ops.saveName,'zscored_by_condition')
        ylims = [-2 2];
    else
        ylims = [0 1];
    end

    for k=ops.minK:ops.maxK

        tic
        if strcmp(ops.method,'kmeans')
            [IDX,C,SUMD,D] = kmeans(X,k,'Replicates',ops.repeats,'Distance',ops.distance);
        elseif strcmp(ops.method,'kmedoids')
            [IDX,C,SUMD,D] = kmedoids(X,k,'Replicates',ops.repeats,'Distance',ops.distance);
        end
        fprintf(1,['Done clustering K=' num2str(k) ' with ' num2str(ops.repeats) ' repeats in ' num2str(toc) ' s\n'])
        
        % reorder clusters to match template order (if specified)
        if ops.template & (ops.minK==ops.maxK)
            C_curr = C; % don't want to overwrite it
            load([CLUSTER_PATH 'output' filesep '_templates' filesep ops.template]); % C
            [IDX,C,SUMD,D] = reorderClusters(X,C,IDX,C_curr,SUMD,D,ops.distance);
        end

        % save clustering output
        if ~exist([SAVE_PATH,'clustering' reliab_dir],'dir'), mkdir([SAVE_PATH,'clustering' reliab_dir]); end
        filename = [SAVE_PATH 'clustering' reliab_dir filesep expt_string '_' elecType '_clusters_K=' num2str(k) split_string reliab_string '.mat'];
        save(filename,'C','-v7.3');
        % append cluster assignments to all_X_table (to be saved at the end)
        eval(strcat("all_X_table.k",num2str(k),"=IDX;"));


        % sort heatmap (X) and IDX due to cluster assignment within each subject
        ii=1;
        for i=1:length(subject_idxs)
            Xs{i}=X(ii:subject_idxs(i)-1,:);
            cluster_ind{i}=IDX(ii:subject_idxs(i)-1);
            ii=subject_idxs(i);
        end
        Xs{i+1}=X(ii:end,:);
        cluster_ind{i+1}=IDX(ii:end);   
        
        XSorted = nan(size(X));
        IDXsorted = nan(size(IDX));
        ii=1;
        for i=1:length(Xs)
            [indSorted, I] = sort(cluster_ind{i});
            XsSorted{i} = Xs{i}(I,:);
            
            if i<length(Xs)
                XSorted(ii:subject_idxs(i)-1,:) = XsSorted{i};
                IDXsorted(ii:subject_idxs(i)-1) = indSorted;
                ii=subject_idxs(i);
            else
                XSorted(ii:end,:)=XsSorted{i};
                IDXsorted(ii:end) = indSorted;
            end
        end
        

        %%%% PLOT %%%%

        h = ERPfigure; set(h,'Position',[0 0 2400 1200],'visible',ops.isPlotVisible)

        % cluster assignments
        ax(1) = subplot(k,n,[6:n:k*n]); imagesc(IDXsorted); hold on;
        set(gca,'XTick',[],'YTick',[])
        if ~isempty(ops.colors)
            colors = ops.colors;
        else
            colors = hsv(k); 
        end
        clustColors = colormap(ax(1),colors);
        pos = get(ax(1),'position'); pos(1) = pos(1)-0.005; pos(3) = pos(3)*0.3; pos(4) = pos(4)*0.92;
        set(ax(1),'position',pos);

        % heatmap of timecourses
        ax(2) = subplot(k,n,[1:n:k*n,5:n:k*n]); imagesc(t,1:size(X,1),XSorted,'CDataMapping','scaled',ylims); hold on;
        pos = get(ax(2),'position'); pos(4) = pos(4)*0.92; set(ax(2),'position',pos);
        for i=1:length(stitch_idxs)-1
            plot([t(stitch_idxs(i)) t(stitch_idxs(i))],[1 size(X,1)],'w','linewidth',4);
        end
        set(gca,'YTick',[])
        set(gca,'Xtick',xlocs,'XTickLabels',xlocslabels,'fontsize',26);
        title(heatmap_title,'fontsize',38,'fontweight','bold')
        gr = colormap(ax(2),'gray');
        colormap(ax(2),gr(70:end,:));
        xlabel({' ','Time (seconds relative to trial start)'},'fontsize',34);
        ylabel({'All Electrodes (sorted)',' '},'fontsize',38);
        
        % centriods of clusters
        for i=1:k
            ax(2+k) = subplot(k,n,[(i*n-(n/3)):(i*n)]); plot(t,C(i,:),'k','linewidth',1.5); hold on;
            ylim(ylims); xlim([0 t(size(X,2))]); set(gca,'Ytick',ylims);
            set(gca,'Xtick',xlocs,'XTickLabel',[],'fontsize',26,'box','off');
            for ii=1:length(stitch_idxs)-1
                plot([t(stitch_idxs(ii)) t(stitch_idxs(ii))],ylims,'--k','linewidth',2); % condition boundary lines
                for iii=1:ops.words-1
                    word_onset = t(stitch_idxs(ii)-length_trial+1)+t(length_word)*iii;
                    plot([word_onset word_onset],ylims,'color','#D3D3D3','linewidth',1); % word lines
                end
            end
            for iii=1:ops.words-1 % repeat for last condition
                word_onset = t(stitch_idxs(end)-length_trial+1)+t(length_word)*iii;
                plot([word_onset word_onset],ylims,'color','#D3D3D3','linewidth',1);
            end
            if i==k
                set(gca,'Xtick',xlocs,'XTickLabel',xlocslabels,'fontsize',26)
            end
            ylabel({['Cluster ' num2str(i)]},'fontsize',34,'fontweight','bold','color',clustColors(i,:));
            ylh = get(gca,'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
            set(ylh,'Rotation',0,'Position',ylp+[cluster_offset 0 0],'VerticalAlignment','middle','HorizontalAlignment','right');
        end

        % save png 
        saveas(gcf,[PLOT_PATH2 expt_string '_' elecType '_K=' num2str(k) split_string reliab_string '.png']);

        % save pdf 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH expt_string '_' elecType '_K=' num2str(k) split_string reliab_string '.pdf'],'pdf')

    end

    % save table with cluster assignments
    writetable(all_X_table,save_table_filename);
end