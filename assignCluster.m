function assignCluster(varargin)
    % function for assigning clusters to MITLangloc subjects from MITSWJNTask clusters

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'distance','correlation'); % or 'squaredeuclidean'
    addParameter(p,'expt1','MITSWJNTask'); % expriment to use for assignment
    addParameter(p,'expt2','MITLangloc'); % experiment to assign
    addParameter(p,'srate',60);
    addParameter(p,'k',3); % assumes clusters from expt 1 are ordered as desired
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs',true);
    parse(p, varargin{:});
    ops = p.Results;

    % --- INITIALIZE --- %

    % paths
    [~,SAVE_PATH] = initialize(ops.saveName);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'clustering' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'clustering' filesep];
    if ~exist(PLOT_PATH,'dir'), mkdir(PLOT_PATH); end
    if ~exist(PLOT_PATH2,'dir'), mkdir(PLOT_PATH2); end

    % file naming 
    DATA_PATH = [SAVE_PATH 'data' filesep];
    if ops.useLangElecs, elecType = 'langElecs'; else, elecType = 'nonLangElecs'; end

    % load in data matrix and cluster assignments from expt 1
    load([DATA_PATH ops.expt1 '_' elecType '_data_for_clustering.mat']); % all_X
    X1 = all_X;
    CLUSTER_PATH = [SAVE_PATH 'clustering' filesep];
    all_X_table1 = readtable([SAVE_PATH 'clustering' filesep ops.expt1 '_' elecType '_cluster_assignments.csv']);

    % load in data matrix and data labels from expt 2
    load([DATA_PATH ops.expt2 '_' elecType '_data_for_clustering.mat']); % all_X
    all_X_table = readtable([DATA_PATH ops.expt2 '_' elecType '_labels_for_clustering.csv']);
    X = all_X;

    % other params
    k = ops.k;
    srate = ops.srate; 
    t = 0:1/srate:((size(X,2)-1)/srate);
    nConds = 2;
    heatmap_title = {' ','Concatenated Timecourses',' ', ...
            'S                                  N'};
    cluster_offset = -0.4;
    t_per_cond = t(length(t)/nConds);

    % saving suffix 
    if contains(ops.saveName,'medoids')
        suffix = '_clusterMean';
    else
        suffix = '';
    end

    
    % --------------------------
    % --- CLUSTER ASSIGNMENT --- %
    % --------------------------

    % --- MEAN OF EXPT1 CLUSTER CENTERS --- %

    % load cluster assignments
    eval(strcat("assignments=all_X_table1.k",num2str(k),";"));
    IDX1 = assignments;

    % get average expt1 cluster response
    C1 = zeros(k,size(X1,2));
    for kk=1:k
        C1(kk,:) = mean(X1(IDX1==kk,:),1);
    end

    % --- ASSIGNING EXPT2 CHANNELS TO EXPT1 CLUSTERS --- %

    IDX = zeros(size(X,1),1);
    for kk=1:size(X,1)
        dists = cell2mat(arrayfun(@(x) pdist2(X(kk,:),C1(x,:),ops.distance),1:k,'UniformOutput',false));
        [~,IDX(kk,1)] = min(dists);
    end

    % append cluster assignments to all_X_table (to be saved at the end)
    eval(strcat("all_X_table.k",num2str(k),"=IDX;"));

    % get average langloc cluster response
    C = zeros(k,size(X,2));
    for kk=1:k
        C(kk,:) = mean(X(IDX==kk,:),1);
    end

    
    % --- PLOTTING CLUSTER ASSIGNMENTS --- %

    % parameters for all plots
    n = 12; % number of columns in subplot

    % stitch idxs for concatenation
    length_trial = size(X,2)/nConds;
    stitch_idxs = length_trial:length_trial:size(X,2);

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

    ylims = [0 1];

    % sort heatmap (X) and IDX due to cluster assignment within each subject
    ii=1;
    for i=1:length(subject_idxs)
        Xs{i}=X(ii:subject_idxs(i),:);
        cluster_ind{i}=IDX(ii:subject_idxs(i));
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
            XSorted(ii:subject_idxs(i),:) = XsSorted{i};
            IDXsorted(ii:subject_idxs(i)) = indSorted;
            ii=subject_idxs(i);
        else
            XSorted(ii:end,:)=XsSorted{i};
            IDXsorted(ii:end) = indSorted;
        end
    end
        

    %%%% --- PLOT --- %%%%

    h = ERPfigure; set(h,'Position',[0 0 2400 1200],'visible',ops.isPlotVisible)

    % cluster assignments
    ax(1) = subplot(k,n,[6:n:k*n]); image(IDXsorted); hold on;
    set(gca,'XTick',[],'YTick',[])
    clustColors = colormap(ax(1),hsv(k));
    pos = get(ax(1),'position'); pos(1) = pos(1)-0.005; pos(3) = pos(3)*0.3; pos(4) = pos(4)*0.92;
    set(ax(1),'position',pos);

    % heatmap of timecourses
    ax(2) = subplot(k,n,[1:n:k*n,5:n:k*n]); imagesc(t,1:size(X,1),XSorted,'CDataMapping','scaled'); hold on;
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
        set(gca,'Xtick',xlocs,'XTickLabel',[],'fontsize',26,'box','off')
        for ii=1:length(stitch_idxs)-1
            plot([t(stitch_idxs(ii)) t(stitch_idxs(ii))],ylims,'--k','linewidth',2);
        end
        if i==k
            set(gca,'Xtick',xlocs,'XTickLabel',xlocslabels,'fontsize',26)
        end
        ylabel({['Cluster ' num2str(i)]},'fontsize',34,'fontweight','bold','color',clustColors(i,:));
        ylh = get(gca,'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
        set(ylh,'Rotation',0,'Position',ylp+[cluster_offset 0 0],'VerticalAlignment','middle','HorizontalAlignment','right');
    end

    % save png 
    saveas(gcf,[PLOT_PATH2 ops.expt2 '_' elecType '_K=' num2str(k) suffix '.png']);

    % save pdf 
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH ops.expt2 '_' elecType '_K=' num2str(k) suffix '.pdf'],'pdf')

    % save table with cluster assignments
    filename = [SAVE_PATH 'clustering' filesep ops.expt2 '_' elecType '_cluster_assignments.csv'];
    writetable(all_X_table,filename);

end