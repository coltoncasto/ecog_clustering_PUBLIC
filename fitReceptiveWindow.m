function fitReceptiveWindow(varargin)

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'useWandJ',true); % MITSWJNTask only
    addParameter(p,'distance','correlation'); 
    addParameter(p,'k',3);
    addParameter(p,'whichKernel','gaussian_wide'); % alternatives: 'square', 'cosine'
    addParameter(p,'numberOfWindows',1); % times the number of words
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs',true);
    addParameter(p,'wordIncrements',false);
    addParameter(p,'nWords',8); % can be 12 for MITLangloc when in right folder
    addParameter(p,'condition','S');
    parse(p, varargin{:});
    ops = p.Results;

    % --- INITIALIZE ---

    % paths
    [CLUSTER_PATH,SAVE_PATH] = initialize(ops.saveName);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'trw' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'trw' filesep];
    if ~exist(PLOT_PATH,'dir')
        mkdir(PLOT_PATH);
    end
    if ~exist(PLOT_PATH2,'dir')
        mkdir(PLOT_PATH2);
    end

    % file naming conventions
    DATA_PATH = [SAVE_PATH 'data' filesep];
    if ops.useLangElecs
        elecType = 'langElecs';
    else
        elecType = 'nonLangElecs';
    end
    if strcmp(ops.experiment,'both')
        expt_string = 'bothMITSWJNTaskandMITLangloc';
    else
        expt_string = ops.experiment;
    end

    % load in data matrix
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering.mat']); % all_X
    X = all_X;

    % load in all trial data matrix
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering_all_trials.mat']); % all_X
    Xtr = all_X;

    % load in data labels
    all_X_table = readtable([SAVE_PATH 'clustering' filesep expt_string '_' elecType '_cluster_assignments.csv']);
    eval(strcat("assignments=all_X_table.k",num2str(ops.k),";"));

    % load in all trials data labels
    all_Xtr_table = readtable([DATA_PATH expt_string '_' elecType '_labels_for_clustering_all_trials.csv'],'Delimiter',',');

    % extract unique subject/electrode identified
    elec_id = arrayfun(@(x) strcat(all_X_table.subject{x},'_',num2str(all_X_table.channel_number(x)),'_',all_X_table.channel_name{x}),1:size(all_X_table,1),'uniformoutput',false)';
    elec_id_tr = arrayfun(@(x) strcat(all_Xtr_table.subject{x},'_',num2str(all_Xtr_table.channel_number(x)),'_',all_Xtr_table.channel_name{x}),1:size(all_Xtr_table,1),'uniformoutput',false)';
    
    % experiment-specific info
    if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
        conds = {'S','W','J','N'};
        nConds = 4;
        nWords = 8;
        cond_string = '_SWJN';
    elseif ~ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask') % only S and N
        conds = {'S','N'};
        nConds = 2;
        nWords = 8;
        cond_string = '_SN';
    else % MITLangloc
        conds = {'S','N'};
        nConds = 2;
        nWords = ops.nWords;
        cond_string = '';
    end

    % load in reliability values (for dot size)
    RELIAB_PATH = [CLUSTER_PATH 'output/_reliability/reliability/'];
    load([RELIAB_PATH expt_string cond_string '_' elecType '_reliability.mat']); % corrs
    corrs_for_lm = corrs; corrs_for_lm(corrs_for_lm<0) = 0;

    % -----------------------------------
    % --- FIT RECEPTIVE WINDOW SIZE --- %
    % -----------------------------------

    % extract only sentence condition
    cond_starts = 1:size(X,2)/nConds:size(X,2);
    cond_ends = size(X,2)/nConds:size(X,2)/nConds:size(X,2)+1;
    X = X(:,cond_starts(1):cond_ends(1));
    Xtr = Xtr(:,cond_starts(1):cond_ends(1));
    
    % construct stim train
    train = zeros(1,size(X,2));
    word_starts = 1:size(X,2)/nWords:size(X,2);
    train(word_starts) = 1;

    % construct kernels to convolve with stim train
    kernels = {};
    if ops.wordIncrements
        window_sizes = (size(X,2)/nWords):(size(X,2)/nWords):(size(X,2)*ops.numberOfWindows); % in samples
        word_inc_string = '_word_increments';
    else
        window_sizes = (size(X,2)/nWords/3):(size(X,2)*ops.numberOfWindows); % in samples
        word_inc_string = '';
    end
    convs = zeros(length(window_sizes),size(X,2));
    switch ops.whichKernel
        case 'gaussian_wide'
            %old:
            % for i=1:length(window_sizes) 
            %     w = window_sizes(i);
            %     points = 1:w;
            %     kernels{i} = normpdf(points,mean(points),length(points)/2);
            %     kernels{i} = (kernels{i}-min(kernels{i}))*...
            %                        (1/(max(kernels{i})-min(kernels{i})));
            %     kernel_x_stim = conv(train,kernels{i});
            %     convs(i,:) = kernel_x_stim(1:size(X,2));
            % end
            for i=1:length(window_sizes) 
                w = window_sizes(i);
                points = 1:w;
                kernels{i} = normpdf(points,mean(points),length(points)/2);
                kernels{i} = (kernels{i}-min(kernels{i}))*...
                                   (1/(max(kernels{i})-min(kernels{i})));
                kernel_x_stim = conv(train,kernels{i});
                convs(i,:) = kernel_x_stim(1:size(X,2));
            end
        case 'gaussian_narrow'
            for i=1:length(window_sizes) 
                w = window_sizes(i);
                points = 1:w;
                kernels{i} = normpdf(points,mean(points),length(points)/16);
                kernels{i} = (kernels{i}-min(kernels{i}))*...
                                   (1/(max(kernels{i})-min(kernels{i})));
                kernel_x_stim = conv(train,kernels{i});
                convs(i,:) = kernel_x_stim(1:size(X,2));
            end
        case 'gaussian_medium'
            for i=1:length(window_sizes) 
                w = window_sizes(i);
                points = 1:w;
                kernels{i} = normpdf(points,mean(points),length(points)/8);
                kernels{i} = (kernels{i}-min(kernels{i}))*...
                                   (1/(max(kernels{i})-min(kernels{i})));
                kernel_x_stim = conv(train,kernels{i});
                convs(i,:) = kernel_x_stim(1:size(X,2));
            end
        case 'square'
            for i=1:length(window_sizes) 
                w = window_sizes(i);
                points = 1:w;
                kernels{i} = ones(1,length(points)); kernels{i}(1) = 0; kernels{i}(end) = 0;
                kernel_x_stim = conv(train,kernels{i});
                convs(i,:) = kernel_x_stim(1:size(X,2));
            end
        case 'cosine'
            for i=1:length(window_sizes) 
                w = window_sizes(i);
                points = 1:w;
                x4cos = linspace(-pi/2,pi/2,length(points));
                kernels{i} = cos(x4cos);
                kernel_x_stim = conv(train,kernels{i});
                convs(i,:) = kernel_x_stim(1:size(X,2));
            end
        case 'linear_asym'
            for i=1:length(window_sizes) 
                w = window_sizes(i);
                points = 1:w;
                kernels{i} = linspace(0,1,length(points));
                kernel_x_stim = conv(train,kernels{i});
                convs(i,:) = kernel_x_stim(1:size(X,2));
            end
        case 'evoked'
           for i=1:length(window_sizes) 
                w = window_sizes(i);
                points = 1:w;
                x4cos = linspace(-pi/2,pi/2,length(points));
                kernels{i} = cos(x4cos);
                kernel_x_stim = conv(train,kernels{i});
                convs(i,:) = kernel_x_stim(1:size(X,2));
                
                evoked_result = zeros(size(kernel_x_stim));
                onsets = find(train);
                %i
                for ii=1:length(onsets)
                    evoked_result(onsets(ii):onsets(ii)+length(points)-1) = evoked_result(onsets(ii):onsets(ii)+length(points)-1) + kernels{i};
                end
                evoked_result = evoked_result(1:size(X,2));
                % figure
                % subplot(2,1,1)
                % plot(convs(i,:))
                % title('Convolution')
                % set(gca,'fontsize',16,'xlim',[0 size(X,2)])
                % hold on
                % subplot(2,1,2)
                % plot(evoked_result)
                % title('Evoked')
                % set(gca,'fontsize',16,'xlim',[0 size(X,2)])
                if ~isequal(convs(i,:),evoked_result)
                    warning('evoked is not equal to conv')
                end
            end     
    end

    % calculate distance between simulated kernels and signal from channels
    dists = pdist2(X,convs,ops.distance); % chans x kernels
    if strcmp(ops.distance,'correlation')
        dists = 1-dists;%this is the correlation
    end
    [M,trws] = max(dists,[],2); % in words
    X_sim = convs(trws,:);
    if ~ops.wordIncrements
        trws = (trws+(size(X,2)/nWords/3)-1)/(size(X,2)/nWords);
    end

    % save distance matrix and TRWs
    TRW_PATH = [SAVE_PATH 'trw' filesep];
    if ~exist(TRW_PATH,'dir'), mkdir(TRW_PATH); end
    filename = [TRW_PATH expt_string '_' elecType '_receptive_window_distances' word_inc_string '_kernel_' ops.whichKernel '.mat'];
    save(filename,'dists','-v7.3');
    filename = [TRW_PATH expt_string '_' elecType '_receptive_window_lengths_words' word_inc_string '_kernel_' ops.whichKernel '.mat'];
    save(filename,'trws','-v7.3');


    % ------------------
    % --- PLOTTING --- %
    % ------------------

    % other parameters
    k = ops.k;
    cluster_legend_names = cell(1,k);
    for kk=1:k
        cluster_legend_names{kk} = ['Cluster #' num2str(kk)];
    end
    colors = hsv(k);
    
    
    % --- EXAMPLE PLOT --- %
    %close all
    h = ERPfigure; set(h,'Position',[10 10 700 1050],'visible',ops.isPlotVisible);
    set(h,'defaultAxesColorOrder',[0,0,0;0,0,0]);
    
    widths_to_plot = [1,2,3,4,6,8];
    colors_grey = [0.05,0.05,0.05;...
                   0.20,0.20,0.20;...
                   0.35,0.35,0.35;...
                   0.50,0.50,0.50;...
                   0.65,0.65,0.65;...
                   0.80,0.80,0.80];
    
    % intergration windows and stim train
    ax = subplot(length(widths_to_plot)+3,1,1:2);
    set(gca,'fontsize',16); box off;
    x_values = (1:length(train))/60;
    
    yyaxis left
    plot(ax,x_values,train,'k','linewidth',2); hold on;
    legend_names = {'stim'};
    set(ax,'box','off','ytick',[],'yticklabels',[]);
    
    yyaxis right
    for w=1:length(widths_to_plot)
        width = widths_to_plot(w);
        y = kernels{(window_sizes/(size(X,2)/nWords))==width};
        plot(ax,x_values(1:min(size(X,2),length(y))),y(1:min(size(X,2),length(y))),...
            'color',colors_grey(w,:),'linestyle','--','marker','none','linewidth',2);
        legend_names = [legend_names,{[num2str(width) ' words']}];
    end
    xlabel('Time (seconds)','fontsize',16); xlim([x_values(1) x_values(end)]);
    legend(legend_names,'location','eastoutside','box','off','location','eastoutside','fontsize',16);
    set(ax,'box','off','ytick',[],'yticklabels',[]);
    title({'Window of Integration',''},'fontsize',16);

    % example of convolution
    for w=1:length(widths_to_plot)
        width = widths_to_plot(w);
        
        ax = subplot(length(widths_to_plot)+3,1,w+3); hold on;
        set(gca,'fontsize',16); box off;
        
        y = convs((window_sizes/(size(X,2)/nWords))==width,:);
        y = (y-min(y))/(max(y)-min(y)); % normalize
        plot(ax,x_values,y,'color',colors_grey(w,:),'linestyle','-','linewidth',2);
        set(ax,'box','off','ytick',[],'yticklabels',[]);
        xlim([x_values(1) x_values(end)]);
        
        onsets=x_values(logical(train));
        for o=1:length(onsets)
            plot([onsets(o) onsets(o)],[0 1],':k','linewidth',1,'displayname','');
        end
        legend_name = {[num2str(width) ' words']};
        legend(legend_name,'location','eastoutside','box','off','location','eastoutside','fontsize',16);
        
        if w==length(widths_to_plot)
            xlabel('Time (seconds)','fontsize',16); 
        else
            set(gca,'xticklabels',[]);
        end
        
        
    end
    
    % save png
    saveas(gcf,[PLOT_PATH2 expt_string '_' elecType '_convolution_example' '_kernel_' ops.whichKernel '.png']);

    % save pdf
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH expt_string '_' elecType '_convolution_example' '_kernel_' ops.whichKernel '.pdf'],'pdf');
    
    
    % --- DISTRIBUTION OF RECEPTIVE WINDOWS --- %
    %close all
    h = ERPfigure; set(h,'Position',[10 10 350 600],'visible',ops.isPlotVisible); hold on
    set(gca,'fontsize',16); box off;
    colorsv = [0.3,0.3,0.3; colors]; % append grey to beginning of colors
    
    % violin plots
    x_locs = -k:1:0;
    v = cell(1,k+1);
    h(1) = subplot(1,1,1); hold on;
    for i=1:(k+1)
        
        % all channels first, then by cluster
        if i > 1
            j = i-1;
            cluster_mask = (assignments==j);
            data = trws(cluster_mask);
        else
            data = trws;
        end
        
        % remove outliers
        data_outliers = data(isoutlier(data,'quartiles'));
        data = data(~isoutlier(data,'quartiles'));
        
        v{1,i} = violinplot(data,1,'Bandwidth',(max(data)-min(data))/10,'ShowMean',true,...
                            'BoxColor',[0.3 0.3 0.3],'ShowData',true,'ViolinColor',colorsv(i,:));
        if i~=(k+1)
            v{1,i} = moveViolin(v{1,i},x_locs(i));
        end
        v{1,i} = invisibleViolin(v{1,i});
        v{1,i}.MeanPlot.Color = 'k';
        v{1,i}.MeanPlot.LineWidth = 4;
        
        % plot data outliers 
        scatter(repmat(round(mean(v{1,i}.ScatterPlot.XData),0),1,length(data_outliers)),data_outliers,80,'k','x','linewidth',2);
    end
    
    set(gca,'XTick',x_locs+1,'xticklabel',{'All Channels','Cluster 1','Cluster 2','Cluster 3'},'fontsize',16); 
    hold on; legend hide;
    if ops.wordIncrements
        ylim([0.5 nWords*ops.numberOfWindows]);
    else
        ylim([0 nWords*ops.numberOfWindows]);
    end
    ylabel({'Window of Integration',''});
    
    % save png
    saveas(gcf,[PLOT_PATH2 expt_string '_' elecType '_trw_distribution_' ops.distance '_K=' num2str(k) word_inc_string '_kernel_' ops.whichKernel '.png']);

    % save pdf
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH expt_string '_' elecType '_trw_distribution_' ops.distance '_K=' num2str(k) word_inc_string '_kernel_' ops.whichKernel '.pdf'],'pdf')
    
    
    % --- TRW VS CORRELATION --- % 
    % closest distance as a function of window size 
    % (dot size by reliability)

    %close all
    h = ERPfigure; set(h,'Position',[10 10 900 700],'visible',ops.isPlotVisible); hold on
    set(gca,'fontsize',16); box off;
    ax = subplot(1,1,1);
    
    % plot words
    onsets = 1:nWords;
    for o=1:length(onsets)
        plot(ax,[onsets(o) onsets(o)],[-0.2,1],'-','color',[0.9,0.9,0.9],'linewidth',1);
    end
    
    colors_all = colors(assignments,:);
    sizes_all = 80*(50.^(abs(corrs)));
    % sizes_all = 50*(2.^(abs(corrs)));
    % sizes_all = 1800*abs(corrs);
    hh = scatter(ax,trws,M,sizes_all,colors_all,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
    xlabel({' ','Window of Integration',' '},'fontsize',16,'fontweight','bold');
    ylabel({' ','Correlation',' '},'fontsize',16,'fontweight','bold');
    % legend(cluster_legend_names,'location','eastoutside','box','off','location','eastoutside','fontsize',24);
    set(ax,'box','off');
    if ops.wordIncrements
        xlim([0.5 nWords*ops.numberOfWindows]);
    else
        xlim([0 nWords*ops.numberOfWindows]);
    end
    
    % make 'legend' with marker sizes
    sizes_legend = 80*(50.^([0.05,0.3,0.6,0.9]));
    hh_legend = scatter(ax,[4,5,6,7],[0,0,0,0],sizes_legend,[0.1,0.1,0.1],...
                        'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);

    % save png
    saveas(gcf,[PLOT_PATH2 expt_string '_' elecType '_trw_vs_' ops.distance '_K=' num2str(k) '_size_by_reliability' word_inc_string  '_kernel_' ops.whichKernel '.png']);

    % save pdf
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH expt_string '_' elecType '_trw_vs_' ops.distance '_K=' num2str(k) '_size_by_reliability' word_inc_string '_kernel_' ops.whichKernel '.pdf'],'pdf');

end