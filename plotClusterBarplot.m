function plotClusterBarplot(varargin)
    % method for visualize kmedoids cluster centers
    % (typically represented by an exemplar electrode)
    % code is duplicated from cluster.m

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'minK',3);
    addParameter(p,'maxK',3);
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs',true);
    addParameter(p,'useWandJ',false); % MITSWJNTask only
    addParameter(p,'colors',[]); % when specified will only work with one value of k
    addParameter(p,'addScatter',true)
    parse(p, varargin{:});
    ops = p.Results;


    % --- INITIALIZE --- %

    % paths
    [~,SAVE_PATH] = initialize(ops.saveName);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'clustering' filesep 'barplots' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'clustering' filesep 'barplots' filesep];
    if ~exist(PLOT_PATH,'dir'), mkdir(PLOT_PATH); end
    if ~exist(PLOT_PATH2,'dir'), mkdir(PLOT_PATH2); end

    % file naming 
    DATA_PATH = [SAVE_PATH 'data' filesep];
    if ops.useLangElecs, elecType = 'langElecs'; else, elecType = 'nonLangElecs'; end
    if strcmp(ops.experiment,'both')
        expt_string = 'bothMITSWJNTaskandMITLangloc';
    elseif strcmp(ops.experiment,'langloc')
        expt_string = 'MGHlangloc';
    else
        expt_string = ops.experiment;
    end
    if ops.addScatter, scatter_string = '_withScatter'; else, scatter_string = ''; end

    % load in data matrix
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering.mat']); % all_X
    X = all_X;

    % load in data labels
    all_X_table = readtable([SAVE_PATH 'clustering' filesep expt_string '_' elecType '_cluster_assignments.csv']);

    % variable definition
    if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
        conditions = {'S','W','J','N'};
        labels = {'Sentences','Words','Jabberwocky','Nonwords'};
        colors = {[0.15 0.15 0.15],[0.35 0.35 0.35],[0.55 0.55 0.55],[0.75 0.75 0.75]};
    else
        conditions = {'S','N'};
        labels = {'Sentences','Nonwords'};
        colors = {[0.15 0.15 0.15],[0.75 0.75 0.75]};
    end
    nConds = length(conditions);
    trial_len_samples = size(X,2)/nConds;
    trial_stitch_idxs = 1:trial_len_samples:size(X,2);


    % --- SAVE CONDITION MEANS TO CLUSTERING FILE --- %

    C_byChan = zeros(size(X,1),nConds);
    for cc=1:nConds
        window = trial_stitch_idxs(cc):trial_stitch_idxs(cc)+trial_len_samples-1;
        C_byChan(:,cc) = mean(X(:,window),2);
        eval(strcat("all_X_table.",conditions(cc),"=C_byChan(:,cc);"));
    end
    % save table with condition means
    writetable(all_X_table,[SAVE_PATH 'clustering' filesep expt_string '_' elecType '_cluster_assignments.csv']);

    
    % --- CLUSTERING WITH ALL VALUES OF K --- %

    for k=ops.minK:ops.maxK

        if ~isempty(ops.colors)
            ccolors = ops.colors;
        else
            ccolors = hsv(k); 
        end

        % load cluster assignments
        eval(strcat("assignments=all_X_table.k",num2str(k),";"));

        % get average cluster response
        C = zeros(k,nConds);
        Csem = zeros(k,nConds);
        Cscatter = cell(k,nConds);
        for kk=1:k
            for cc=1:nConds
                window = trial_stitch_idxs(cc):trial_stitch_idxs(cc)+trial_len_samples-1;
                C(kk,cc) = mean(mean(X(assignments==kk,window),2),1);
                Csem(kk,cc) = std(mean(X(assignments==kk,window),2),[],1) / sqrt(sum(assignments==kk)); % over electrodes
                Cscatter(kk,cc) = {mean(X(assignments==kk,window),2)};
            end
        end

        %%%% PLOT %%%%

        if k==3
            % pos = [0 0 430 200];
            pos = [0 0 800 400];
        elseif k==5
            pos = [0 0 700 200];
        elseif k==2
            % pos = [0 0 250 200];
            pos = [0 0 500 400];
        elseif k==4
            pos = [0 0 580 200];
        end

        h = ERPfigure; set(h,'Position',pos,'visible',ops.isPlotVisible)

        % loop through clusters
        for i=1:k

            hh(i) = subplot(1,k,i); hold on

            % loop through conditions
            for j=1:nConds
                % bar plot
                h(j) = bar(j,C(i,j),'displayname',labels{j});
                set(h(j), 'FaceColor', colors{j});
                if ops.addScatter
                    try
                        assert(max(Cscatter{i,j})<0.8)
                    catch
                        max(Cscatter{i,j})
                    end
                    if j==1
                        swarmchart(ones(length(Cscatter{i,j}),1)*j,Cscatter{i,j},20,'k','filled','MarkerEdgeColor',[0.3,0.3,0.3],'XJitter','density','XJitterWidth',0.3);
                    else
                        swarmchart(ones(length(Cscatter{i,j}),1)*j,Cscatter{i,j},20,colors{j},'filled','MarkerEdgeColor','k','XJitter','density','XJitterWidth',0.3);
                    end
                end
                if j==1
                    errorbar(j,C(i,j),Csem(i,j),'color',[0.3,0.3,0.3],'LineWidth',2,'Capsize',8);
                else
                    errorbar(j,C(i,j),Csem(i,j),'k','LineWidth',2,'Capsize',8);
                end
            end

            % other plotting parameters
            ylim([0,0.8]); 
            set(gca,'FontSize',12,'Box','off');
            subtitle({'';['    Cluster #' num2str(i)]},'FontSize',12,'fontweight','bold','Color',ccolors(i,:));
            if i==1
                set(gca,'xtick',[]);
                set(gca,'ytick',[0 0.5]);
            elseif i==median(1:k)
                set(gca,'xtick',1:nConds,'xticklabels',conditions);
                set(gca,'ytick',[]);
            else
                set(gca,'xtick',[]);
                set(gca,'ytick',[]);
            end

            % adjusting position of plots
            h_pos=get(hh(i),'Position'); 
            if i/median(1:k) < 1
                set(hh(i),'Position',[h_pos(1)-0.04 h_pos(2) h_pos(3)*1.2 h_pos(4)]);
            elseif i==median(1:k)
                set(hh(i),'Position',[h_pos(1)-0.02 h_pos(2) h_pos(3)*1.2 h_pos(4)]);
            elseif i/median(1:k) > 1
                set(hh(i),'Position',[h_pos(1)+0 h_pos(2) h_pos(3)*1.2 h_pos(4)]);
            end

        end
        
        % save png 
        saveas(gcf,[PLOT_PATH2 expt_string '_' elecType '_K=' num2str(k) '_barplot' scatter_string '.png']);

        % save pdf 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH expt_string '_' elecType '_K=' num2str(k) '_barplot' scatter_string '.pdf'],'pdf')

    end

end