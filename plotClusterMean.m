function plotClusterMean(varargin)
    % method for visualize kmedoids cluster centers
    % (typically represented by an exemplar electrode)
    % code is duplicated from cluster.m
    %%%%%%%%%%%%%%% code review - DONE (Tamar June 12 2023)

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'minK',1);%currently only k=3 is supported
    addParameter(p,'maxK',10);
    addParameter(p,'srate',60);
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs',true);
    addParameter(p,'useWandJ',false); % MITSWJNTask only
    addParameter(p,'showVar',true); % use varplot to show the variance
    addParameter(p,'showAllElect',false); % plot all electrodes behind the mean
    addParameter(p,'split',[]); % 'odd' or 'even'
    addParameter(p,'colors',[]);
    addParameter(p,'words',8);
    addParameter(p,'filterByReliability',false);
    addParameter(p,'reliabThresh',0.1);
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
        reliab_string = ['_filtered_by_reliability_threshold_' strrep(sprintf('%0.2f',ops.reliabThresh),'0.','')];
    else 
        PLOT_PATH = [SAVE_PATH 'plots' filesep 'clustering' filesep];
        PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'clustering' filesep];
        reliab_dir = '';
        reliab_string = '';
    end
    if ~exist(PLOT_PATH,'dir'), mkdir(PLOT_PATH); end
    if ~exist(PLOT_PATH2,'dir'), mkdir(PLOT_PATH2); end

    % file naming 
    DATA_PATH = [SAVE_PATH 'data' filesep];
    if ops.useLangElecs, elecType = 'langElecs'; else, elecType = 'nonLangElecs'; end
    if ops.split, split_string = ['_' ops.split]; else, split_string = ''; end
    if strcmp(ops.experiment,'both')
        expt_string = 'bothMITSWJNTaskandMITLangloc';
    else
        expt_string = ops.experiment;
    end
    if ops.showAllElect
        all_string = '_all_electrodes';
    else
        all_string = '';
    end

    % load in data matrix
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering' split_string '.mat']); % all_X\
    X = all_X;
%     % load in all trials matrix
%     load([DATA_PATH expt_string '_' elecType '_data_for_clustering_all_trials' split_string '.mat']); % all_X\
%     Xtr = all_X;

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
    end

    % load in data labels
    all_X_table = readtable([SAVE_PATH 'clustering' reliab_dir filesep expt_string '_' elecType '_cluster_assignments' split_string reliab_string '.csv']);
    
%     % load in all trials data labels
%     all_Xtr_table = readtable([DATA_PATH expt_string '_' elecType '_labels_for_clustering_all_trials' split_string '.csv']);
%     nTrials = nan(size(X,1),1);
%     change_channel_idxs = find(diff(all_Xtr_table.channel_number));
%     nTrials(1) = change_channel_idxs(1);
%     nTrials(2:end-1) = diff(change_channel_idxs);
%     nTrials(end) = length(all_Xtr_table.channel_number)-change_channel_idxs(end);

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
    if contains(ops.saveName,'zscored_by_condition')
        ylims = [-2 2];
    else
        ylims = [0 1];
    end

    for k=ops.minK:ops.maxK

        % load cluster assignments
        eval(strcat("assignments=all_X_table.k",num2str(k),";"));
        IDX = assignments;

        % get average cluster response
        C = zeros(k,size(X,2));
        for kk=1:k
            C(kk,:) = mean(X(IDX==kk,:),1);
        end
        %keep all electrodes in each cluster:
        Call = cell(k,1);
        for kk=1:k
            Call{kk,1} = X(IDX==kk,:);
        end

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
        % colorbar
        
        % centriods of clusters
        for i=1:k
            ax(2+k) = subplot(k,n,[(i*n-(n/3)):(i*n)]); 
            if ops.showAllElect
                plot(t,Call{i,:}','Color',[0.5 0.5 0.5],'linewidth',0.25);hold on
            end
            if ops.showVar
                varplot(t,Call{i,:}','ci',0.99,'k','linewidth',2); hold on;
            else
                plot(t,C(i,:),'k','linewidth',2); hold on;
            end
            
            ylim(ylims); xlim([0 t(size(X,2))]); set(gca,'Ytick',ylims);
            set(gca,'Xtick',xlocs,'XTickLabel',[],'fontsize',26,'box','off')
            for ii=1:length(stitch_idxs)-1
                plot([t(stitch_idxs(ii)) t(stitch_idxs(ii))],ylims,'--k','linewidth',2);
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
        saveas(gcf,[PLOT_PATH2 expt_string '_' elecType '_K=' num2str(k) split_string reliab_string '_clusterMean' all_string '.png']);

        % save pdf 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH expt_string '_' elecType '_K=' num2str(k) split_string reliab_string '_clusterMean' all_string '.pdf'],'pdf')

    end

end