function plot_channels_clustDist(varargin)

    %% initialize

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'k',3);%clustering solution for k clusters
    addParameter(p,'whichk',2);%which cluster to plot out of k
    addParameter(p,'distTo','medoids');%options: 'medoids', 'means'
    addParameter(p,'plotsec',15);%how many seconds to plot on x axis
    addParameter(p,'nselChans',30);%how many signals to display
    addParameter(p,'colorby','order');%option are: 'order' , 'trw' , 'reliability'
    addParameter(p,'transparency','reliability');%option are: 'off', 'reliability'
    addParameter(p,'dispText','trw');%option are: 'order' , 'trw'
    addParameter(p,'dispSub',false);%option are: true, false
    addParameter(p,'savePlot',true);%
    addParameter(p,'distance','correlation');
    addParameter(p,'srate',60);
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'whichKernel','gaussian_wide'); % alternatives: 'square', 'cosine'
    addParameter(p,'useLangElecs',true);
    addParameter(p,'useWandJ',true); % MITSWJNTask only
    addParameter(p,'split',[]); % 'odd' or 'even'
    addParameter(p,'template',[]); % filename of template to use in output/_templates/ folder
    addParameter(p,'colors',[]); % when specified will only work with one value of k
    addParameter(p,'plotind',false); % plot individual channels
    
    parse(p, varargin{:});
    ops = p.Results;

    % --- INITIALIZE ---

    % paths
    [CLUSTER_PATH,SAVE_PATH] = initialize(ops.saveName);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'channels' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'channels' filesep];
    if ~exist(PLOT_PATH,'dir')
        mkdir(PLOT_PATH);
    end

    if ~exist(PLOT_PATH2,'dir')
        mkdir(PLOT_PATH2);
    end
    TRW_PATH = [SAVE_PATH filesep 'trw'];
    RELIABILITY_PATH = [SAVE_PATH filesep 'reliability'];
    CLUSTERING_PATH = [SAVE_PATH 'clustering' filesep];

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
    save_table_filename = [CLUSTERING_PATH filesep expt_string '_' elecType '_cluster_assignments' split_string '.csv'];
    if exist(save_table_filename,'file')
        all_X_table = readtable(save_table_filename);
    else % start from scratch
        all_X_table = readtable([DATA_PATH expt_string '_' elecType '_labels_for_clustering' split_string '.csv']);
    end

    % load cluster assignments
    eval(strcat("assignments=all_X_table.k",num2str(ops.k),";"));
    all_X_table.k3 = categorical(all_X_table.k3);
    IDX_all = assignments;

    %load medoids
    load([CLUSTERING_PATH filesep expt_string '_' elecType '_clusters_K=' num2str(ops.k) '.mat'])%loaded as C
    medoids = C;

    %add fitted trw to the table
    trw_matfile = [TRW_PATH filesep ops.experiment '_langElecs_receptive_window_lengths_words_kernel_' ops.whichKernel '.mat'];
    load([trw_matfile])%loaded as trws

    all_X_table=[all_X_table, table(trws)];

    % load subject number map
    sub_num_map = subjectNumberMap();

    %add reliability to the table
    switch ops.experiment
        case 'MITSWJNTask'
            if ops.useWandJ
                reliability_matfile = [RELIABILITY_PATH filesep ops.experiment '_SWJN_langElecs_reliability.mat'];
                nConds = 4;
            else
                reliability_matfile = [RELIABILITY_PATH filesep ops.experiment '_SN_langElecs_reliability.mat'];
                nConds = 2;
            end
        case 'MITLangloc'
            reliability_matfile = [RELIABILITY_PATH filesep ops.experiment '_langElecs_reliability.mat'];
            nConds = 2;
    end
    load([reliability_matfile])%loaded as corrs
    reliability = corrs;
    all_X_table=[all_X_table, table(reliability)];

    % get average cluster response
    if contains(ops.saveName,'assignedFrom')
        % load in SWJN SN data matrix
        load([DATA_PATH 'MITSWJNTask_' elecType '_data_for_clustering' split_string '.mat']); % all_X
        X_SWJN = all_X;
        save_table_filename_SWJN = [CLUSTERING_PATH filesep 'MITSWJNTask_' elecType '_cluster_assignments' split_string '.csv'];
        all_X_table_SWJN = readtable(save_table_filename_SWJN);
        eval(strcat("assignments_SWJN=all_X_table_SWJN.k",num2str(ops.k),";"));
        all_X_table_SWJN.k3 = categorical(all_X_table_SWJN.k3);
        IDX_all_SWJN = assignments_SWJN;
        C = zeros(ops.k,size(X_SWJN,2));
        for kk=1:ops.k
            C(kk,:) = mean(X_SWJN(IDX_all_SWJN==kk,:),1);
        end
    else
        C = zeros(ops.k,size(X,2));
        for kk=1:ops.k
            C(kk,:) = mean(X(IDX_all==kk,:),1);
        end
    end

    switch ops.distTo
        case 'medoids'
            C_all = medoids;
        case 'means'
            C_all = C; % don't want to get overwritten
    end

    % other params
    k = ops.k;
    srate = ops.srate;
    t = 0:1/srate:((size(X,2)-1)/srate);
    t_per_cond = t(length(t)/nConds);

    if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
        nConds = 4;
    else % only S and N
        nConds = 2;
    end

    %% prep data

    ic = ops.whichk;
    c = C_all(ic,:);
    Xc = X(assignments==ic,:);
    nChans = size(Xc,1);

    Xc_table = all_X_table(assignments==ic,:);
    dist = diag(1-corr(Xc',repmat(c,[size(Xc,1),1])') );
    [ds,idx] = sort(dist);
    Xs = Xc(idx,:);
    Xs_table = Xc_table(idx,:);
    %%%%flip up down channel order to plot the channels most similar to the
    %cluster means at the top:
    %Xs = Xs(end:-1:1,:);
    %Xs_table = Xs_table(end:-1:1,:);
    if ops.nselChans <=nChans
        nselChans = ops.nselChans;
    else
        nselChans = nChans;
    end
    Xss = Xs(1:nselChans,:);
    Xss_table = Xs_table(1:nselChans,:);


    %% plot

    D = Xss';
    channel_labels = Xss_table.channel_name;
    clean_channels = 1:nselChans;

    valid_channels = ones(1,nselChans);

    t_len = ops.plotsec;
    curr_sample_freq = 60;
    downsample = false;
    decimation_freq=300;

    stitch_index=1;

    save = false;

    % ------------------------------
    % FORMAT & NORMALIZE SIGNAL
    % ------------------------------
    x_norm_cell = [];
    for k=1:length(stitch_index) % number of separate data files with signal
        
        if k == length(stitch_index) % signal for file stops at end of matrix
            stop = size(D,1);
        else % signal for file stops before stitch index of next file
            stop = stitch_index(k+1)-1;
        end
        
        D_ = D(stitch_index(k):stop,:);
        
        % formatting
        x_cell = mat2cell(D_',ones(1,size(D_,2))); % make signal from each electrode a cell
        
        % downsampling
        if downsample %not supported
            decimation_factor = curr_sample_freq/decimation_freq;
            x_cell = cellfun(@(x) downsample(x,decimation_factor),x_cell,'uni',false);
            curr_sample_freq = decimation_freq;
        end
        
        % normalizing
        min_max = mean(cell2mat(cellfun(@(y) prctile(y,[0 100]),x_cell,'UniformOutput',false))); % mean 5th and 95th percentiles of ALL electrodes
        x_norm_cell_ = cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false); % normalize
        x_norm_cell_ = arrayfun(@(x) x_norm_cell_{x}*valid_channels(x),1:size(x_norm_cell_,1),'uni',false)'; % set signal of noisy channels to 0
        
        if length(stitch_index) > 1
            x_norm_cell = [x_norm_cell, x_norm_cell_];
        else
            x_norm_cell = x_norm_cell_;
        end
        
    end

    % combine separate x_norm_cell columns
    if length(stitch_index) > 1
        for k=1:length(obj.stitch_index)-1
            x_norm_cell(:,k+1) = arrayfun(@(x) {[x_norm_cell{x,k}, x_norm_cell{x,k+1}]},[1:size(x_norm_cell,1)])';
        end
        x_norm_cell = x_norm_cell(:,length(stitch_index));
    end
    assert(size(x_norm_cell,2)==1,'x_norm_cell not in the correct format');


    % ------------
    % PLOT SIGNAL
    % ------------
    t_length = t_len*curr_sample_freq;

    switch ops.colorby
        case 'order'
            %col_inf=inferno(floor(.8*size(x_norm_cell,1)));
            %col_vir=viridis(floor(.8*size(x_norm_cell,1)));
            %colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];
            colors_more = inferno(floor(1.2*size(x_norm_cell,1)));
            colors = colors_more(1:size(x_norm_cell,1),:);
        case 'trw'
            
            TRWs = Xss_table.trws;
            map = jet;TRWs_linspace = linspace(1,8,length(map));
            mapname = 'jet';
            colors = nan(size(x_norm_cell,1),3);
            for ii=1:size(x_norm_cell,1)
                [m, j] = min(abs(TRWs_linspace - TRWs(ii)));
                colors(ii,:) = map(j,:);
            end
            titleString = 'TRW (words)';
            
            barlabels = {'0','2','4','6','8'};
        case 'reliability'
            
            rels = abs(Xss_table.reliability);
            map = cool;rels_linspace = linspace(0,1,length(map));
            mapname = 'cool';
            colors = nan(size(x_norm_cell,1),3);
            for ii=1:size(x_norm_cell,1)
                [m, j] = min(abs(rels_linspace - rels(ii)));
                colors(ii,:) = map(j,:);
            end
            titleString = {'Reliability','corr coef'};
            %barlabels = {num2str(min(rels),1),num2str((max(rels)+min(rels))*0.25,1),num2str((max(rels)+min(rels))*0.5,1),num2str((max(rels)+min(rels))*0.75,1),num2str(max(rels),1)};
            barlabels = {'0','0.25','0.5','0.75','1'};
    end


    hf=figure;
    height = 2000*nselChans./90;
    switch nConds
        case 4
            width = 900;
        case 2
            width = 700;
    end

    set(gcf,'position',[1000,1000,width,height],'visible',ops.isPlotVisible); % dimensions for my external monitor
    ax = axes('position',[.1,.05,.7,.9]);

    hold on
    time_stamps = [1:size(x_norm_cell{1},2)]/curr_sample_freq;

    hold on
    H = arrayfun(@(x) plot(time_stamps,x_norm_cell{x}+x,'color',colors(x,:),'linewidth',2,'DisplayName',strrep(sprintf('ch %d, tag %s',x,channel_labels{x}),'_',' ')),1:size(x_norm_cell,1));
    H0 = arrayfun(@(x) plot(time_stamps,zeros(size(time_stamps))+x,'color','k','linewidth',1,'LineStyle',':','HandleVisibility','off'),1:size(x_norm_cell,1));

    switch ops.transparency
        case 'none'
            
        case 'reliability'
            
            rels = Xss_table.reliability;rels(rels<0)=0;
            maxrel = 50.^0.92;
            alphas = 50.^rels;alphas=alphas./maxrel;
            for i=1:numel(H)
                c = get(H(i), 'Color');
                set(H(i), 'Color', [c alphas(i)]);
            end
    end


    set(gcf,'doublebuffer','on');
    set(ax,'ytick',[1:size(x_norm_cell,1)]);
    set(ax,'yticklabel','');

    switch ops.dispText
        case 'order'
            arrayfun(@(x) text(0,x+0.5,num2str(x),'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',20),clean_channels);
        case 'trw'
            TRWs = Xss_table.trws;
            arrayfun(@(x) text(0,x+0.5,num2str(TRWs(x),2),'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',17),clean_channels);
            text(0,clean_channels(end) + 2,{'TRW', '(words)'},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',25);
    end

    switch ops.dispSub
        case true
            sub = cellfun(@(x) sub_num_map(x),Xss_table.subject,'uniformoutput',false);
            arrayfun(@(x) text(ops.plotsec,x+0.5,['P' num2str(sub{x})],'Color',colors(x,:),'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',17),clean_channels);
        case false
    end

    set(ax,'ylim',[0,size(x_norm_cell,1)+2]);
    ax.XAxis.TickLength = [0.005,0.01];
    ax.YAxis.TickLength = [0.005,0.01];
    set(ax,'xlim',[0 t_len]);
    pos = get(ax,'position');
    Newpos = [pos(1) pos(2)-0.1 pos(3) 0.05];
    xmax=max(time_stamps);

    %%%% To make this a scrollable GUI
    %S = ['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(t_len) '])'];
    %h = uicontrol('style','slider','units','normalized','position',Newpos,'callback',S,'min',0,'max',xmax-t_len);
    %datacursormode on
    %waitfor(findobj('type','figure','number',1));

    %hl=legend('fontsize',16,'location','eastoutside');

    colormap(mapname)
    hc=colorbar('Ticks',[0:0.25:1],'TickLabels',barlabels,'fontsize',20);
    set(hc,'Position',[0.87,0.05,0.02,0.9])
    ht = get(hc,'Title');
    set(ht ,'String',titleString);
    set(gca,'fontsize',20)

    title({['Cluster ' num2str(ic) ', ' num2str(nselChans) ' electrodes, ' ops.experiment], ['sorted by dist from '  ops.distTo]},'fontsize',25)

    % xticks
    xlocslabels = repmat({'','1','2','3'},1,nConds);
    xlocs = repmat([0 1 2 3],1,nConds);
    u = 1;
    for i=5:4:length(xlocs)-1
        xlocs(:,i:i+3) = xlocs(:,i:i+3)+(u*t_per_cond);
        u = u+1;
    end
    set(gca,'Xtick',xlocs,'XTickLabels',xlocslabels,'fontsize',26);

    % stitch idxs for concatenation
    length_trial = size(X,2)/nConds;
    stitch_idxs = length_trial:length_trial:size(X,2);

    % plot white vertical lines to separate conds
    for i=1:length(stitch_idxs)-1
        plot([t(stitch_idxs(i)) t(stitch_idxs(i))],[1 size(X,1)],'w','linewidth',4);
    end

    savename = [ops.experiment '_channels_cluster' num2str(ic) '_sortedby_distance_colorby_' ops.colorby '_alpha_' ops.transparency '_dispText_' ops.dispText '_nChans_' num2str(nselChans) ];

    saveas(gcf,[PLOT_PATH filesep savename],'fig')
    saveas(gcf,[PLOT_PATH2 filesep savename],'png')

    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH filesep savename],'pdf')

end
