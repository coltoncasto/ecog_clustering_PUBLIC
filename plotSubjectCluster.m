function plotSubjectCluster(varargin)
    % method for plotting cluster assignments by subject
    % does NOT re-cluster

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'k',3); % only does it for 1 value of k to save space
    addParameter(p,'srate',60);
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs',true);
    addParameter(p,'useWandJ',false); % MITSWJNTask only
    addParameter(p,'split',[]);
    addParameter(p,'colors',[]);
    addParameter(p,'words',8);
    addParameter(p,'showVar',true); % use varplot to show the variance
    addParameter(p,'showAllElect',false); % plot all electrodes behind the mean
    parse(p, varargin{:});
    ops = p.Results;


    % --- INITIALIZE ---

    % paths
    [~,SAVE_PATH] = initialize(ops.saveName);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'clustering' filesep 'bySubject' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'clustering' filesep 'bySubject' filesep];
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
    elseif strcmp(ops.experiment,'langloc')
        expt_string = 'MGHlangloc';
    else
        expt_string = ops.experiment;
    end

    % load in data matrix
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering' split_string '.mat']); % all_X
    X = all_X;

    % load in data labels
    all_X_table = readtable([SAVE_PATH 'clustering' filesep expt_string '_' elecType '_cluster_assignments' split_string '.csv']);

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

    % saving suffix 
    if contains(ops.saveName,'medoids')
        suffix = '_clusterMean';
    else
        suffix = '';
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
    if contains(ops.saveName,'zscored_by_condition')
        ylims = [-2 2];
    else
        ylims = [0 1];
    end

    k = ops.k;
    if ~isempty(ops.colors)
        colors = ops.colors;
    else
        colors = hsv(k); 
    end

    for sub=1:length(unique_subs)
        close all
        subject = unique_subs{sub};

        % extract subject's data from matrix and table
        sub_idxs = find(strcmp(all_X_table.subject,subject));
        Xsub = X(sub_idxs,:);
        all_X_table_sub = all_X_table(sub_idxs,:);

        % load cluster assignments
        eval(strcat("assignments=all_X_table_sub.k",num2str(k),";"));
        IDX = assignments;

        % get average cluster response
        C = zeros(k,size(Xsub,2));
        for kk=1:k
            C(kk,:) = mean(Xsub(IDX==kk,:),1);
        end
        % keep all electrodes in each cluster:
        Call = cell(k,1);
        for kk=1:k
            Call{kk,1} = Xsub(IDX==kk,:);
        end

        % sort heatmap (X) and IDX due to cluster assignment within each subject
        [IDXsorted, I] = sort(IDX);
        XSorted = X(I,:);


        %%%% PLOT %%%%

        h = ERPfigure; set(h,'Position',[0 0 2400 1200],'visible',ops.isPlotVisible)

        % cluster assignments
        ax(1) = subplot(k,n,[6:n:k*n]); image(IDXsorted); hold on;
        set(gca,'XTick',[],'YTick',[])
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
            ax(2+k) = subplot(k,n,[(i*n-(n/3)):(i*n)]); 
            if ops.showAllElect
                plot(t,Call{i,:}','Color',[0.5 0.5 0.5],'linewidth',0.25);hold on
            end
            if ops.showVar & (size(Call{i,:},1)>1)
                varplot(t,Call{i,:}','ci',0.99,'k','linewidth',2); hold on;
            else
                plot(t,C(i,:),'k','linewidth',2); hold on;
            end
            ylim(ylims); xlim([0 t(size(X,2))]); set(gca,'Ytick',ylims);
            set(gca,'Xtick',xlocs,'XTickLabel',[],'fontsize',26,'box','off')
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
            ylabel({['Cluster ' num2str(i)]},'fontsize',34,'fontweight','bold','color',colors(i,:));
            ylh = get(gca,'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
            set(ylh,'Rotation',0,'Position',ylp+[cluster_offset 0 0],'VerticalAlignment','middle','HorizontalAlignment','right');
        end

        % save png 
        saveas(gcf,[PLOT_PATH2 subject '_' expt_string '_' elecType '_K=' num2str(k) split_string suffix '.png']);

        % save pdf 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH subject '_' expt_string '_' elecType '_K=' num2str(k) split_string suffix '.pdf'],'pdf')

    end

end