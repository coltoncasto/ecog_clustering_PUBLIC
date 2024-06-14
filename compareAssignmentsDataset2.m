function compareAssignmentsDataset2(varargin)
%%%% This function was created for the 2nd NHB revision to compare the
%%%% "winner-take-all" method to the full clustering for Dataset 2
%%%% (only works for Dataset 2, k=3)
    
    p = inputParser();
    addParameter(p,'solution1','MITLangloc_SN_8words_kmedoids_correlation_FIGURES'); % folder name of solution 1
    addParameter(p,'solution2','MITLangloc_assignedFrom_MITSWJNTask_SN_kmedoids_correlation_FIGURES'); % folder name of solution 2
    addParameter(p,'experiment1','MITLangloc'); % 'MITSWJNTask' or 'MITLangloc'
    addParameter(p,'elecType1','langElecs');
    addParameter(p,'elecType2','langElecs');
    addParameter(p,'k1',3); % k value for solution 1
    addParameter(p,'k2',3); % k value for solution 2
    addParameter(p,'split1',[]); % '_odd' or '_even'
    addParameter(p,'split2',[]); % '_odd' or '_even'
    addParameter(p,'useWandJ',false); % MITSWJNTask only
    addParameter(p,'isPlotVisible',false);
    parse(p, varargin{:});
    ops = p.Results;

    % --------------------
    % --- INITIALIZE --- %
    % --------------------

    % --- SOLUTION 1 --- %

    [~,SAVE_PATH1] = initialize(ops.solution1);
    PLOT_PATH1 = [SAVE_PATH1 'plots' filesep 'comparisons' filesep];
    PLOT_PATH21 = [SAVE_PATH1 'plots' filesep 'pngs' filesep 'comparisons' filesep];
    if ~exist(PLOT_PATH1,'dir'), mkdir(PLOT_PATH1); end
    if ~exist(PLOT_PATH21,'dir'), mkdir(PLOT_PATH21); end

    % loading solution 1
    expt1 = split(ops.solution1,'_'); expt1 = expt1{1}; % TODO - doesn't work with 'both'
    load([SAVE_PATH1 'data' filesep expt1 '_' ops.elecType1 '_data_for_clustering' ops.split1 '.mat']); X1 = all_X;
    all_X_table1 = readtable([SAVE_PATH1 'clustering' filesep expt1 '_' ops.elecType1 '_cluster_assignments' ops.split1 '.csv']);
    % get average cluster response
    eval(strcat("assignments=all_X_table1.k",num2str(ops.k1),";"));
    IDX1 = assignments; C1 = zeros(ops.k1,size(X1,2));
    for kk=1:ops.k1, C1(kk,:) = mean(X1(IDX1==kk,:),1); end

    % --- SOLUTION 2 --- %
    
    [~,SAVE_PATH2] = initialize(ops.solution2);
    PLOT_PATH2 = [SAVE_PATH2 'plots' filesep 'comparisons' filesep];
    PLOT_PATH22 = [SAVE_PATH2 'plots' filesep 'pngs' filesep 'comparisons' filesep];
    if ~exist(PLOT_PATH2,'dir'), mkdir(PLOT_PATH2); end
    if ~exist(PLOT_PATH22,'dir'), mkdir(PLOT_PATH22); end

    % loading solution 2
    expt2 = split(ops.solution2,'_'); expt2 = expt2{1}; % TODO - doesn't work with 'both'
    load([SAVE_PATH2 'data' filesep expt2 '_' ops.elecType2 '_data_for_clustering' ops.split2 '.mat']); X2 = all_X;
    all_X_table2 = readtable([SAVE_PATH2 'clustering' filesep expt2 '_' ops.elecType2 '_cluster_assignments' ops.split2 '.csv']);
    % get average cluster response
    eval(strcat("assignments=all_X_table1.k",num2str(ops.k2),";"));
    IDX2 = assignments; C2 = zeros(ops.k2,size(X2,2));
    for kk=1:ops.k2, C2(kk,:) = mean(X2(IDX2==kk,:),1); end

    T_WTA=all_X_table2;
    T=all_X_table1;

    Trans= nan(3,3);
    whichElects = cell(3,3);

    for iclust = 1:3
        for iwta = 1:3
            Trans(iclust,iwta) = numel(find(T.k3==iclust & T_WTA.k3==iwta));
            whichElects{iclust,iwta} = find(T.k3==iclust & T_WTA.k3==iwta); 
        end
    end

    hf=figure;
    set(hf,'Position',[100 100 450 400],'visible',ops.isPlotVisible)
    imagesc([1:3],[1:3],Trans)
    cm=gray(100);
    colormap(cm(1:100,:))
    set(gca,'XTick',[1:3],'YTick',[1:3],'fontsize',20)
    xlabel('WTA assignment')
    ylabel('Clustering assignment')
    for iclust = 1:3
        for iwta = 1:3
            text(iwta,iclust,num2str(Trans(iclust,iwta)),'FontSize',22)
        end
    end

    % save png
    saveas(gcf,[PLOT_PATH21 ops.solution1 ops.split1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 ops.split2 '_K=' num2str(ops.k2) '_n_electrodes_assigned.png']);

    % save pdf 
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH1 ops.solution1 ops.split1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 ops.split2 '_K=' num2str(ops.k2) '_n_electrodes_assigned.pdf'],'pdf')
         
    % save png
    saveas(gcf,[PLOT_PATH22 ops.solution2 ops.split2 '_K=' num2str(ops.k2) '_vs_' ops.solution1 ops.split1 '_K=' num2str(ops.k1) '_n_electrodes_assigned.png']);

    % save pdf 
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH2 ops.solution2 ops.split2 '_K=' num2str(ops.k2) '_vs_' ops.solution1 ops.split1 '_K=' num2str(ops.k1) '_n_electrodes_assigned.pdf'],'pdf')
        

    %% NOTE - THINGS ARE HARDCODED
    all_X = X1;
    hf=figure;
    set(hf,'Position',[100 100 450 400],'visible',ops.isPlotVisible)

    % other params
    srate = 60; 
    t = (1:size(all_X,2))/srate;
    
    nConds = 2;
    heatmap_title = {' ','Concatenated Timecourses',' ', ...
            'S                                  N'};
    cluster_offset = -0.4;
    t_per_cond = t(length(t)/nConds);

    % xticks 
    nConds=2;
    xlocslabels = repmat({'','1','2','3'},1,nConds);
    xlocs = repmat([0 1 2 3],1,nConds);
    u = 1;
    for i=5:4:length(xlocs)-1
        xlocs(:,i:i+3) = xlocs(:,i:i+3)+(u*t_per_cond);
        u = u+1; 
    end

    % stitch idxs for concatenation
    length_trial = size(all_X,2)/nConds;
    stitch_idxs = length_trial:length_trial:size(all_X,2);
    length_word = length_trial/8;

    i=0;
    for iclust = 1:3
        for iwta = 1:3
            i=i+1;
            subplot(3,3,i)
            if Trans(iclust,iwta) >1
                varplot(t,all_X(whichElects{iclust,iwta},:)','ci',0.99,'linewidth',0.75,'Color','k')
            else
                plot(t,all_X(whichElects{iclust,iwta},:),'linewidth',0.75,'Color','k')
            end
            ylim([0,1])
            xlim([0 t(end)])
            set(gca,'Xtick',xlocs,'XTickLabels',xlocslabels);
            hold on
            for ii=1:length(stitch_idxs)-1
                plot([t(stitch_idxs(ii)) t(stitch_idxs(ii))],[0 1],'--k','linewidth',1);
            end
            for ii=1:length(stitch_idxs)
                for iii=1:8-1
                    word_onset = t(stitch_idxs(ii)-length_trial+1)+t(length_word)*iii;
                    plot([word_onset word_onset],[0 1],'color',[0.8 0.8 0.8],'linewidth',0.25); % word lines
                end
            end
        end
    end

    % save png
    saveas(gcf,[PLOT_PATH21 ops.solution1 ops.split1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 ops.split2 '_K=' num2str(ops.k2) '_electrodes_signals.png']);

    % save pdf 
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH1 ops.solution1 ops.split1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 ops.split2 '_K=' num2str(ops.k2) '_electrodes_signals.pdf'],'pdf')
         
    % save png
    saveas(gcf,[PLOT_PATH22 ops.solution2 ops.split2 '_K=' num2str(ops.k2) '_vs_' ops.solution1 ops.split1 '_K=' num2str(ops.k1) '_electrodes_signals.png']);

    % save pdf 
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH2 ops.solution2 ops.split2 '_K=' num2str(ops.k2) '_vs_' ops.solution1 ops.split1 '_K=' num2str(ops.k1) '_electrodes_signals.pdf'],'pdf')
        
end
