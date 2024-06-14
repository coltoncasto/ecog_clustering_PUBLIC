function calculatePartialCorrelation(varargin)
    % Calculates partial correlation between electrodes and cluster medoids

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'k',3);
    addParameter(p,'split',[]);
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'whichKernel','gaussian_wide'); % alternatives: 'square', 'cosine'
    addParameter(p,'useLangElecs',true);
    addParameter(p,'useWandJ',true); % MITSWJNTask only
    addParameter(p,'reliabThresholds',[-1]); % if array will repeat for multiple reliab threshs
    addParameter(p,'srate',60);
    addParameter(p,'words',8);
    parse(p, varargin{:});
    ops = p.Results;

    % --- INITIALIZE --- %

    % paths
    [CLUSTER_PATH,SAVE_PATH] = initialize(ops.saveName);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'partial' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'partial' filesep];
    if ~exist(PLOT_PATH,'dir'), mkdir(PLOT_PATH); end
    if ~exist(PLOT_PATH2,'dir'), mkdir(PLOT_PATH2); end


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

    % --- LOAD --- %

    % load in averaged matrix
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering' split_string '.mat']); % all_X
    X = all_X;

    % load in medoids
    MED_PATH = [SAVE_PATH 'clustering' filesep];
    eval(strcat("medoid_file = [MED_PATH expt_string '_' elecType '_clusters_K=' num2str(ops.k) '.mat'];"));
    load(medoid_file); % C
        
    % load in cluster assignments
    all_X_table = readtable([SAVE_PATH 'clustering' filesep expt_string '_' elecType '_cluster_assignments' split_string '.csv']);
    eval(strcat("assignments=all_X_table.k",num2str(ops.k),";"));
    IDX = assignments;
    
    % extract unique subject/electrode identified
    elec_id = arrayfun(@(x) strcat(all_X_table.subject{x},'_',num2str(all_X_table.channel_number(x)),'_',all_X_table.channel_name{x}),1:size(all_X_table,1),'uniformoutput',false);

    % experiment-specific info
    if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
        nConds = 4;
        cond_string = '_SWJN';
    elseif ~ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask') % only S and N
        nConds = 2;
        cond_string = '_SN';
    else % MITLangloc
        nConds = 2;
        cond_string = '';
    end

    % load in reliability values 
    RELIAB_PATH = [CLUSTER_PATH 'output/_reliability/reliability/'];
    load([RELIAB_PATH expt_string cond_string '_' elecType '_reliability.mat']); % corrs
    
    % load trw values
    TRW_PATH = [SAVE_PATH filesep 'trw'];
    load([TRW_PATH filesep ops.experiment '_' elecType '_receptive_window_lengths_words_kernel_' ops.whichKernel '.mat']); % trws


    % ------------------------------
    % --- PARTIAL CORRELATIONS --- %
    % ------------------------------

    analyses = {'partial','normal'};
    for a=1:length(analyses)
        for k=1:ops.k
            for r=1:length(ops.reliabThresholds)
                
                % filter by reliability
                if length(ops.reliabThresholds)==1
                    curr_thresh = ops.reliabThresholds;
                else 
                    curr_thresh = ops.reliabThresholds(r);
                end
                X_to_use = X(corrs>curr_thresh,:);
                IDX_to_use = IDX(corrs>curr_thresh);

                % partial correlation
                k_idxs = zeros(ops.k,1); k_idxs(k) = 1;
                k_idxs = logical(k_idxs);
                if strcmp(analyses{a},'partial')
                    pcs = partialcorr([C(k_idxs,:)', X_to_use'],C(~k_idxs,:)');
                    part_string = '_partial';
                elseif strcmp(analyses{a},'normal')
                    pcs = corrcoef([C(k_idxs,:)', X_to_use']);
                    part_string = '';
                end
                pcs = pcs(2:end,1);

                
                % --- PLOTTING ---

                close all
                h = ERPfigure; set(h,'Position',[0 0 800 700],'visible',ops.isPlotVisible)

                % histogram 
                yd = 0.8; % y displacement
                colors = hsv(ops.k);
                for kk=1:ops.k
                    corrs_curr = pcs(IDX_to_use==kk);
                    histogram(corrs_curr,12,'FaceColor',colors(kk,:),'FaceAlpha',0.4); hold on;

                    % dashed line with mean corr
                    mean_curr = mean(corrs_curr,1,'omitnan');
                    yl = ylim;
                    plot([mean_curr mean_curr],[yl(1) yl(2)],'--','Color',colors(kk,:),'linewidth',3);
                    text(mean_curr+0.03,yl(2)*yd,['mean = ' num2str(round(mean_curr,4))],'Color',colors(kk,:),'fontsize',11);
                    yd = yd - 0.1;
                end

                % other plotting params
                set(gca,'fontsize',14); box off; xlim([-0.5 1]);
                xlabel({' ','Correlation'},'fontsize',18,'fontweight','bold');
                ylabel({'Electrodes',' '},'fontsize',18,'fontweight','bold');

                % file naming
                if curr_thresh==-1
                    reliab_string = '';
                else
                    reliab_string = ['_reliability_threshold_' strrep(sprintf('%0.2f',curr_thresh),'0.','')];
                end

                % save png 
                saveas(gcf,[PLOT_PATH2 expt_string '_' elecType part_string '_corrs_K=' num2str(k) reliab_string '.png']);

                % save pdf 
                set(gcf,'Units','inches');
                screenposition = get(gcf,'Position');
                set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
                saveas(gcf,[PLOT_PATH expt_string '_' elecType part_string '_corrs_K=' num2str(k) reliab_string '.pdf'],'pdf')

                
                % --- PLOTTING INIDIVIDUAL ELECS --- 
                if strcmp(analyses{a},'partial') && (k~=3) && (curr_thresh==-1)
                    if k==1
                        other_k = 2;
                    elseif k==2
                        other_k = 1;
                    end
                    
                    % electrodes that show mixed profile
                    eois = find((corrs>0.3) & (IDX_to_use==other_k) & (pcs>0.2));
                    for e=1:length(eois)
                        curr_e = eois(e);
                        
                        close all
                        h = ERPfigure; set(h,'Position',[0 0 800 200],'visible',ops.isPlotVisible); hold on;
                        
                        % plotting params
                        colors = hsv(ops.k);
                        ylims = [0 1];
                        t = (1:size(X,2))/ops.srate;
                        t_per_cond = t(length(t)/nConds);
                        
                        % stitch idxs for concatenation
                        length_trial = size(X,2)/nConds;
                        stitch_idxs = length_trial:length_trial:size(X,2);
                        length_word = length_trial/ops.words;
                        
                        % word and condition onset lines
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

                        % electrode timecourse
                        % varplot(t,X(curr_e,:),'ci',0.99,'k','linewidth',2); hold on;
                        plot(t,C(2,:),'color',colors(2,:),'linewidth',1.15); hold on;
                        plot(t,C(1,:),'color',colors(1,:),'linewidth',1.15); hold on;
                        plot(t,X(curr_e,:),'k','linewidth',1.5); hold on;

                        % xticks 
                        xlocslabels = repmat({'','1','2','3'},1,nConds);
                        xlocs = repmat([0 1 2 3],1,nConds);
                        u = 1;
                        for i=5:4:length(xlocs)-1
                            xlocs(:,i:i+3) = xlocs(:,i:i+3)+(u*t_per_cond);
                            u = u+1; 
                        end
                        
                        % other plotting params
                        ylim(ylims); xlim([0 t(size(X,2))]); set(gca,'Ytick',ylims);
                        set(gca,'Xtick',xlocs,'XTickLabel',xlocslabels,'fontsize',14,'box','off')
                        title({['Cluster ' num2str(other_k) ' ' strrep(elec_id{curr_e},'_',' ')],...
                               ['partial corr w Cluster ' num2str(k) ' medoid: ' num2str(round(pcs(curr_e),4)) ...
                                '; reliability: ' num2str(round(corrs(curr_e),4)) ...
                                '; TRW: ' num2str(round(trws(curr_e),2))]},'fontsize',15);
                        
                        IND_PATH = [PLOT_PATH filesep 'individual_timecourses' filesep];
                        if ~exist(IND_PATH,'dir'), mkdir(IND_PATH); end
                        IND_PATH2 = [PLOT_PATH2 filesep 'individual_timecourses' filesep];
                        if ~exist(IND_PATH2,'dir'), mkdir(IND_PATH2); end
                        
                        % save png 
                        saveas(gcf,[IND_PATH2 expt_string '_' elecType '_Cluster' num2str(other_k) '_elec_w_high_partial_correlation_to_Cluster' num2str(k) '_medoid_' elec_id{curr_e} '.png']);

                        % save pdf 
                        set(gcf,'Units','inches');
                        screenposition = get(gcf,'Position');
                        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
                        saveas(gcf,[IND_PATH expt_string '_' elecType '_Cluster' num2str(other_k) '_elec_w_high_partial_correlation_to_Cluster' num2str(k) '_medoid_' elec_id{curr_e} '.pdf'],'pdf')

                    end
                                      
                    % repeat for the average
                    close all
                    h = ERPfigure; set(h,'Position',[0 0 800 200],'visible',ops.isPlotVisible); hold on;
              
                    % word and condition onset lines
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

                    % electrode timecourse
                    plot(t,C(2,:),'color',colors(2,:),'linewidth',1.15); hold on;
                    plot(t,C(1,:),'color',colors(1,:),'linewidth',1.15); hold on;
                    varplot(t,X(eois,:)','ci',0.99,'k','linewidth',1.5); hold on;
                        
                    % other plotting params
                    ylim(ylims); xlim([0 t(size(X,2))]); set(gca,'Ytick',ylims);
                    set(gca,'Xtick',xlocs,'XTickLabel',xlocslabels,'fontsize',14,'box','off')
                    title({['Cluster ' num2str(other_k) ' high partial corr w Cluster ' num2str(k) ' medoid'],...
                            ['(n=' num2str(length(eois)) ')']},'fontsize',15);
                        
                    % save png 
                    saveas(gcf,[IND_PATH2 expt_string '_' elecType '_Cluster' num2str(other_k) '_elec_w_high_partial_correlation_to_Cluster' num2str(k) '_medoid_average_n=' num2str(length(eois)) '.png']);

                    % save pdf 
                    set(gcf,'Units','inches');
                    screenposition = get(gcf,'Position');
                    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
                    saveas(gcf,[IND_PATH expt_string '_' elecType '_Cluster' num2str(other_k) '_elec_w_high_partial_correlation_to_Cluster' num2str(k) '_medoid_average_n=' num2str(length(eois)) '.pdf'],'pdf')
  
                end
            end

        end

    end

end