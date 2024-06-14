function reduceDimensionality(varargin)

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'method','PCA'); % or 'tSNE' or 'dist'
    addParameter(p,'distance','correlation'); % for tSNE and dist visualizations
    addParameter(p,'minK',1);
    addParameter(p,'maxK',10);
    addParameter(p,'elecSize',20);% 60 is a good size when shading by reliability
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs',true);
    addParameter(p,'useWandJ',true);
    addParameter(p,'shadeByReliability',false); % set alpha by split half reliability
    addParameter(p,'colors',[]); % when specified will only work with one value of k
    addParameter(p,'signalType','unipolar');
    parse(p, varargin{:});
    ops = p.Results;

    % --- INITIALIZE ---

    % paths
    [CLUSTER_PATH,SAVE_PATH] = initialize(ops.saveName);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'dimensionality' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'dimensionality' filesep];
    if ~exist(PLOT_PATH,'dir')
        mkdir(PLOT_PATH);
    end
    if ~exist(PLOT_PATH2,'dir')
        mkdir(PLOT_PATH2);
    end

    % load in data matrix
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

    % load in data
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering.mat']); % all_X
    X = all_X;

    % load in data labels
    all_X_table = readtable([SAVE_PATH 'clustering' filesep expt_string '_' elecType '_cluster_assignments.csv']);

    % alpha values 
    if ops.shadeByReliability
        if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
            cond_string = '_SWJN';
        elseif strcmp(ops.experiment,'MITSWJNTask')% only 2 conds
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
        corrs(corrs<0) = 0; 
        corrs = corrs/max(corrs);
        
        legend_alphas = [0.05,0.3,0.6,0.9];
        legend_alphas = legend_alphas/max(corrs);

        reliab_string = '_shaded_by_reliability';

    else
        corrs = ones(size(all_X_table,1),1);
        reliab_string = '';
    end


    % -------------------------------
    % --- REDUCE DIMENSIONALITY --- %
    % -------------------------------


    % --- PRINCIPAL COMPONENT ANALYSIS (PCA) ---
    if strcmp(ops.method,'PCA')
        [coeff,score,latent,tsquared,explained,mu] = pca(X);
        % Z = X * coeff(:,1:2); % project X onto the first 2 principal components
        Z = score(:,1:2);
        dist_string = '';
    elseif strcmp(ops.method,'dist')
        D = squareform(pdist(X,ops.distance)); % visualize clusters on eigenvectors of distance space
        [V,Deig] = eig(D);
        Z = D * V(:,1:3);
        dist_string = ['_pairwise_' ops.distance];
    elseif strcmp(ops.method,'tSNE')
        Z = tsne(X,'Distance',ops.distance);
        dist_string = ['_' ops.distance];
    end

    % ------------------
    % --- PLOTTING --- %
    % ------------------

    % CHANNELS PLOTTED ON FIRST TWO PRINCIPAL COMPONENTS
    for k=ops.minK:ops.maxK
        eval(strcat("assignments=all_X_table.k",num2str(k),";"));
        cluster_legend_names = cell(1,k);
        for kk=1:k
            cluster_legend_names{kk} = ['Cluster #' num2str(kk)];
        end

        % colors for plotting
        if ~isempty(ops.colors)
            colors = ops.colors;
        else
            colors = hsv(k); 
        end
        black = zeros(size(assignments,1),3);

        close all
        h = ERPfigure; set(h,'Position',[10 10 700 660],'visible',ops.isPlotVisible) 
        hold on

        % scatter plot
        ax = subplot(1,1,1);
        if ops.shadeByReliability
            hh = scatter(ax,Z(:,1),Z(:,2),ops.elecSize,colors(assignments,:),'filled');
            hh.AlphaData = corrs; hh.MarkerFaceAlpha = 'flat';
            % alpha legend
            % dim1_min = min(Z(:,1));
            % dim2_min = min(Z(:,2));
            % hhl = scatter(ax,[dim1_min,dim1_min,dim1_min,dim1_min],...
            %                  [dim2_min,dim2_min,dim2_min,dim2_min],...
            %                  ops.elecSize,[0,0,0;0,0,0;0,0,0;0,0,0;],'filled');
            % hhl.AlphaData = legend_alphas; hhl.MarkerFaceAlpha = 'flat';
        else
            hh = gscatter(ax,Z(:,1),Z(:,2),assignments,colors,[],ops.elecSize);
        end
        xlabel({' ','Principal Component #1',' '},'fontsize',16,'fontweight','bold');
        ylabel({' ','Principal Component #2',' '},'fontsize',16,'fontweight','bold');
        % TODO - fix legend for shading by reliability
        % legend(cluster_legend_names,'location','eastoutside','box','off','location','eastoutside','fontsize',24);
        set(ax,'box','off','xtick',0,'ytick',0);

        % save png
        saveas(gcf,[PLOT_PATH2 expt_string '_' elecType '_' ops.method dist_string '_K=' num2str(k) reliab_string '.png']);

        % save pdf
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH expt_string '_' elecType '_' ops.method dist_string '_K=' num2str(k) reliab_string '.pdf'],'pdf')

    end

end