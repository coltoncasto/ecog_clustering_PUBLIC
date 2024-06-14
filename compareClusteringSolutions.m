function compareClusteringSolutions(varargin)
    % assumes that clusters are in MATCHED order across solutions being compared

    p = inputParser();
    addRequired(p,'solution1'); % folder name of solution 1
    addRequired(p,'solution2'); % folder name of solution 2, if comparing stability, include here not solution 1
    addParameter(p,'experiment1','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc'
    addParameter(p,'useWandJ',true); % MITSWJNTask only
    addParameter(p,'elecType1','langElecs');
    addParameter(p,'elecType2','langElecs');
    addParameter(p,'k1',3); % k value for solution 1
    addParameter(p,'k2',3); % k value for solution 2
    addParameter(p,'method','kmedoids'); % can be 'assign'
    addParameter(p,'distance','correlation');
    addParameter(p,'repeats',100);
    addParameter(p,'compareAssignments',false); % only for the same set of channels (not across experiments!)
    addParameter(p,'doHeatmap',true);
    addParameter(p,'averageTrials',true); % use average condition responses
    addParameter(p,'matchedClusters',[1,3]); % clusters from solution 1 that match the clusters from solution 2, to be used when k2<k1 (the order specified matters)
    addParameter(p,'doSignificance',false);
    addParameter(p,'iterations',1000);
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
    load([SAVE_PATH1 'data' filesep expt1 '_' ops.elecType1 '_data_for_clustering.mat']); X1 = all_X;
    all_X_table1 = readtable([SAVE_PATH1 'clustering' filesep expt1 '_' ops.elecType1 '_cluster_assignments.csv']);
    % get average cluster response
    eval(strcat("assignments=all_X_table1.k",num2str(ops.k1),";"));
    IDX1 = assignments; C1 = zeros(ops.k1,size(X1,2));
    for kk=1:ops.k1, C1(kk,:) = mean(X1(IDX1==kk,:),1); end

    % --- SOLUTION 3 --- % (if required)

    expt3 = [];
    split_string2 = '';
    if contains(ops.solution2,'STABILITY') % need to look at both odd even splits, now three comparisons
        split_string3 = '_odd';
        [~,SAVE_PATH3] = initialize(ops.solution2);
        expt3 = split(ops.solution2,'_'); expt3 = expt3{1}; % TODO - doesn't work with 'both'
        load([SAVE_PATH3 'data' filesep expt3 '_' ops.elecType2 '_data_for_clustering' split_string3 '.mat']); X3 = all_X;
        all_X_table3 = readtable([SAVE_PATH3 'clustering' filesep expt3 '_' ops.elecType2 '_cluster_assignments' split_string3 '.csv']);
        % get average cluster response
        eval(strcat("assignments=all_X_table3.k",num2str(ops.k2),";"));
        IDX3 = assignments; C3 = zeros(ops.k2,size(X3,2));
        for kk=1:ops.k2, C3(kk,:) = mean(X3(IDX3==kk,:),1); end
        split_string2 = '_even';

        % load all trials if desired
        if ~ops.averageTrials
            load([SAVE_PATH3 'data' filesep expt3 '_' ops.elecType2 '_data_for_clustering' split_string3 '_all_trials.mat']); X_all_trials3 = all_X;
            all_X_table_all_trials3 = readtable([SAVE_PATH3 'data' filesep expt3 '_' ops.elecType2 '_labels_for_clustering' split_string3 '_all_trials.csv']);
            % find electrode idxs for averaging after shuffling trials
            elec_idxs3 = zeros(size(X3,1),2); % starts in column 1, ends in column 2
            idx_start = 1; elec_num = 1;
            prev_name = all_X_table_all_trials3.channel_name{1};
            for i=2:size(all_X_table_all_trials3,1)
                curr_name = all_X_table_all_trials3.channel_name{i};
                if ~strcmp(prev_name,curr_name) || (i==size(all_X_table_all_trials3,1))
                    elec_idxs3(elec_num,2) = i;
                    elec_idxs3(elec_num,1) = idx_start;
                    idx_start = i+1;
                    elec_num = elec_num + 1; 
                end
                prev_name = curr_name;
            end
        end
    end

    % --- SOLUTION 2 --- %

    [~,SAVE_PATH2] = initialize(ops.solution2);
    PLOT_PATH2 = [SAVE_PATH2 'plots' filesep 'comparisons' filesep];
    PLOT_PATH22 = [SAVE_PATH2 'plots' filesep 'pngs' filesep 'comparisons' filesep];
    if ~exist(PLOT_PATH2,'dir'), mkdir(PLOT_PATH2); end
    if ~exist(PLOT_PATH22,'dir'), mkdir(PLOT_PATH22); end

    % loading solution 2
    expt2 = split(ops.solution2,'_'); expt2 = expt2{1}; % TODO - doesn't work with 'both'
    load([SAVE_PATH2 'data' filesep expt2 '_' ops.elecType2 '_data_for_clustering' split_string2 '.mat']); X2 = all_X;
    all_X_table2 = readtable([SAVE_PATH2 'clustering' filesep expt2 '_' ops.elecType2 '_cluster_assignments' split_string2 '.csv']);
    % get average cluster response
    eval(strcat("assignments=all_X_table2.k",num2str(ops.k2),";"));
    IDX2 = assignments; C2 = zeros(ops.k2,size(X2,2));
    for kk=1:ops.k2, C2(kk,:) = mean(X2(IDX2==kk,:),1); end

    % load all trials if desired
    if ~ops.averageTrials
        load([SAVE_PATH2 'data' filesep expt2 '_' ops.elecType2 '_data_for_clustering' split_string2 '_all_trials.mat']); X_all_trials2 = all_X;
        all_X_table_all_trials2 = readtable([SAVE_PATH2 'data' filesep expt2 '_' ops.elecType2 '_labels_for_clustering' split_string2 '_all_trials.csv']);
        % find electrode idxs for averaging after shuffling trials
        elec_idxs2 = zeros(size(X2,1),2); % starts in column 1, ends in column 2
        idx_start = 1; elec_num = 1;
        prev_name = all_X_table_all_trials2.channel_name{1};
        for i=2:size(all_X_table_all_trials2,1)
            curr_name = all_X_table_all_trials2.channel_name{i};
            if ~strcmp(prev_name,curr_name) || (i==size(all_X_table_all_trials2,1))
                elec_idxs2(elec_num,2) = i;
                elec_idxs2(elec_num,1) = idx_start;
                idx_start = i+1;
                elec_num = elec_num + 1; 
            end
            prev_name = curr_name;
        end
    end

    % ------------------------------------------------
    % --- COMPARE SOLUTIONS (CORRELATION MATRIX) --- %
    % ------------------------------------------------

    if ops.doHeatmap

        if expt3
            fullC = [C1; C2; C3];
            axis_label =  'FULL vs EVEN vs ODD';
            position = [0 0 800 800];
        else
            fullC = [C1; C2];
            axis_label = {strrep(ops.solution1,'_',' '),'vs',strrep(ops.solution2,'_',' ')};
            position = [0 0 550 550];
        end
        corrs = corrcoef(fullC');

        % plot 
        close all
        h = ERPfigure; set(h,'Position',position,'visible',ops.isPlotVisible)

        h = heatmap(round(corrs,2),'colormap',gray); 
        set(gca,'fontsize', 14);
        h.XDisplayLabels(1:size(fullC,1)) = {''};
        h.YDisplayLabels(1:size(fullC,1)) = {''};        
        h.XLabel = axis_label;
        % addTitle(h,['k=' num2str(ops.k)]);
        
        % save png
        saveas(gcf,[PLOT_PATH21 ops.solution1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 '_K=' num2str(ops.k2) '.png']);

        % save pdf 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH1 ops.solution1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 '_K=' num2str(ops.k2) '.pdf'],'pdf')
            
        % save png
        saveas(gcf,[PLOT_PATH22 ops.solution2 '_K=' num2str(ops.k2) '_vs_' ops.solution1 '_K=' num2str(ops.k1) '.png']);

        % save pdf 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH2 ops.solution2 '_K=' num2str(ops.k2) '_vs_' ops.solution1 '_K=' num2str(ops.k1) '.pdf'],'pdf')

    end

    % ----------------------------------------
    % --- COMPARE SOLUTION (ASSIGNMENTS) --- %
    % ----------------------------------------

    if ops.compareAssignments 
        if expt3, perc_same_solution3 = round(sum(IDX1==IDX3)/length(IDX1),3); end
        perc_same_solution2 = round(sum(IDX1==IDX2)/length(IDX1),3);
        
        for s=1:2
            % write output to text file for saving
            eval(strcat("COMP_PATH = [SAVE_PATH",num2str(s)," 'comparisons' filesep];"));
            if ~exist(COMP_PATH,'dir'), mkdir(COMP_PATH); end
            if s==1
                filename = [COMP_PATH ops.solution1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 '_K=' num2str(ops.k2) '.txt'];
            elseif s==2
                filename = [COMP_PATH ops.solution2 '_K=' num2str(ops.k2) '_vs_' ops.solution1 '_K=' num2str(ops.k1) '.txt'];
            end
            fileID = fopen(filename,'w');
            fprintf(fileID,'Percent of clusters assigned to the same clusters in %s as in %s : \t%5.3f\n\n',...
                                        ops.solution1,...
                                        [ops.solution2 split_string2],...
                                        perc_same_solution2);
            if expt3
                fprintf(fileID,'Percent of clusters assigned to the same clusters in %s as in %s : \t%5.3f',...
                                        ops.solution1,...
                                        [ops.solution2 split_string3],...
                                        perc_same_solution3);
            end
            fclose(fileID);
        end
    end


    % ------------------------------------------
    % --- SIGNIFICANCE W/ TRIALS SCRAMBLED --- %
    % ------------------------------------------

    % file naming convention
    if ~ops.averageTrials
        trial_string = '_all_trials';
    else
        trial_string = '';
    end

    % --- SOLUTION 2 --- %

    if ops.doHeatmap & ops.doSignificance
        
        fprintf(1,'\nCALCULATING SIGNIFICANCE W/ TRIALS SCRAMBLED FOR SOLUTION 2\n');

        if ops.useWandJ && strcmp(ops.experiment1,'MITSWJNTask')
            nConds = 4;
        else % only S and N
            nConds = 2;
        end

        % select matrix to perform permutations on
        if ~ops.averageTrials
            X_use2 = X_all_trials2;
        else
            X_use2 = X2;
        end

        % chose matching clusters when k2<k1
        if ops.k2<ops.k1
            C1_temp = C1(ops.matchedClusters,:);
        elseif ops.k2==ops.k1
            C1_temp = C1;
        else
            error('Code not designed for k2>k1 at this time');
        end

        totNchans = size(X_use2,1);
        dists_same_cluster_scr = nan(ops.k2,ops.iterations); % k x iters
              
        tic;
        for iter = 1:ops.iterations

            % randomly shuffle trials
            Xn = nan(size(X_use2));
            cond_starts = 1:size(X_use2,2)/nConds:size(X_use2,2);
            cond_ends = size(X_use2,2)/nConds:size(X2,2)/nConds:size(X_use2,2)+1;
            for i=1:nConds
                Xn(:,cond_starts(i):cond_ends(i)) = X_use2(randperm(totNchans),cond_starts(i):cond_ends(i));
            end

            % take mean of shuffled trials if using all trials
            if ~ops.averageTrials
                Xn_new = nan(size(X2));
                for elec=1:size(elec_idxs2,1)
                    Xn_new(elec,:) = mean(Xn(elec_idxs2(elec,1):elec_idxs2(elec,2),:),1);
                end
                Xn = Xn_new;
            end

            % --- CLUSTERING --- %
            if strcmp(ops.method,'assign')
                IDX = zeros(size(Xn,1),1);
                for kk=1:size(Xn,1)
                    dists = cell2mat(arrayfun(@(x) pdist2(Xn(kk,:),C1(x,:),ops.distance),1:ops.k2,'UniformOutput',false));
                    [~,IDX(kk,1)] = min(dists);
                end
            elseif strcmp(ops.method,'kmeans')
                [IDX,Cn,SUMD,D] = kmeans(Xn,ops.k2,'Replicates',ops.repeats,'Distance',ops.distance);
            elseif strcmp(ops.method,'kmedoids')
                [IDX,Cn,SUMD,D] = kmedoids(Xn,ops.k2,'Replicates',ops.repeats,'Distance',ops.distance);
            end
                
            % ---- EVALUATE RESULTS ---- %
            if ~strcmp(ops.method,'assign')
                [IDX,Cn,SUMD,D] = reorderClusters(Xn,C1_temp,IDX,Cn,SUMD,D,ops.distance);
            end
                        
            % calculate mean of cluster (in case using kmedoids)
            Cn_means = zeros(ops.k2,size(Xn,2));
            for kk=1:ops.k2
                Cn_means(kk,:) = mean(Xn(IDX==kk,:),1);
            end

            % calculate distance between this solution and the full solution (by cluster)
            dists_curr = pdist2(C1,Cn_means,ops.distance);
            dists_same_cluster_scr(:,iter) = ones(size(diag(dists_curr)))-diag(dists_curr);

            if mod(iter,50)==0
                fprintf(1,'Done with 50 iterations in %0.4f seconds\n',toc);
                tic;
            end

        end

        % save null distribution
        COMP_PATH = [SAVE_PATH1 'comparisons' filesep];
        if ~exist(COMP_PATH,'dir'), mkdir(COMP_PATH); end
        filename = [COMP_PATH ops.solution1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 split_string2 '_K=' num2str(ops.k2) '_null_distribution' trial_string '.mat'];
        null = dists_same_cluster_scr;
        save(filename,'null','-v7.3');
        null2 = null;

        COMP_PATH = [SAVE_PATH2 'comparisons' filesep];
        if ~exist(COMP_PATH,'dir'), mkdir(COMP_PATH); end
        filename = [COMP_PATH ops.solution2 split_string2 '_K=' num2str(ops.k2) '_vs_' ops.solution1 '_K=' num2str(ops.k1) '_null_distribution' trial_string '.mat'];
        save(filename,'null','-v7.3');

        % --- SOLUTION 3 --- % (if applicable)

        if expt3 

            fprintf(1,'\nCALCULATING SIGNIFICANCE W/ TRIALS SCRAMBLED FOR SOLUTION 3\n');

            % select matrix to perform permutations on
            if ~ops.averageTrials
                X_use3 = X_all_trials3;
            else
                X_use3 = X3;
            end

            totNchans = size(X_use3,1);
            dists_same_cluster_scr = nan(ops.k2,ops.iterations); % k x iters
               
            tic;
            for iter = 1:ops.iterations
        
                % randomly shuffle trials
                Xn = nan(size(X_use3));
                cond_starts = 1:size(X_use3,2)/nConds:size(X_use3,2);
                cond_ends = size(X_use3,2)/nConds:size(X_use3,2)/nConds:size(X_use3,2)+1;
                for i=1:nConds
                    Xn(:,cond_starts(i):cond_ends(i)) = X_use3(randperm(totNchans),cond_starts(i):cond_ends(i));
                end

                % take mean of shuffled trials if using all trials
                if ~ops.averageTrials
                    Xn_new = nan(size(X3));
                    for elec=1:size(elec_idxs3,1)
                        Xn_new(elec,:) = mean(Xn(elec_idxs3(elec,1):elec_idxs3(elec,2),:),1);
                    end
                    Xn = Xn_new;
                end
        
                % --- CLUSTERING --- %
                if strcmp(ops.method,'assign')
                    IDX = zeros(size(Xn,1),1);
                    for kk=1:size(Xn,1)
                        dists = cell2mat(arrayfun(@(x) pdist2(Xn(kk,:),C1(x,:),ops.distance),1:ops.k2,'UniformOutput',false));
                        [~,IDX(kk,1)] = min(dists);
                    end
                elseif strcmp(ops.method,'kmeans')
                    [IDX,Cn,SUMD,D] = kmeans(Xn,ops.k2,'Replicates',ops.repeats,'Distance',ops.distance);
                elseif strcmp(ops.method,'kmedoids')
                    [IDX,Cn,SUMD,D] = kmedoids(Xn,ops.k2,'Replicates',ops.repeats,'Distance',ops.distance);
                end
                    
                % ---- EVALUATE RESULTS ---- %
                [IDX,Cn,SUMD,D] = reorderClusters(Xn,C1_temp,IDX,Cn,SUMD,D,ops.distance);
                            
                % calculate mean of cluster (in case using kmedoids)
                Cn_means = zeros(ops.k2,size(Xn,2));
                for kk=1:ops.k2
                    Cn_means(kk,:) = mean(Xn(IDX==kk,:),1);
                end
        
                % calculate distance between this solution and the full solution (by cluster)
                dists_curr = pdist2(C1,Cn_means,ops.distance);
                dists_same_cluster_scr(:,iter) = ones(size(diag(dists_curr)))-diag(dists_curr);

                if mod(iter,50)==0
                    fprintf(1,'Done with 50 iterations in %0.4f seconds\n',toc);
                    tic;
                end
        
            end
        
            % save null distribution
            COMP_PATH = [SAVE_PATH1 'comparisons' filesep];
            filename = [COMP_PATH ops.solution1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 split_string3 '_K=' num2str(ops.k2) '_null_distribution' trial_string '.mat'];
            null = dists_same_cluster_scr;
            save(filename,'null','-v7.3');
            null3 = null;

            COMP_PATH = [SAVE_PATH2 'comparisons' filesep];
            filename = [COMP_PATH ops.solution2 split_string3 '_K=' num2str(ops.k2) '_vs_' ops.solution1 '_K=' num2str(ops.k1) '_null_distribution' trial_string '.mat'];
            save(filename,'null','-v7.3');

        end

        % ------------------------------------------------
        % --- PLOTTING NULL DISTRIBUTIONS BY CLUSTER --- %
        % ------------------------------------------------

        if ops.k2==3
            pos = [0 0 400 600];
        elseif ops.k2==5
            pos = [0 0 400 900];
        elseif ops.k2==2
            pos = [0 0 400 450];
        end

        if ops.k1==ops.k2
            k = ops.k1;
            real_values2 = [];
            real_values3 = []; % may go unused
            if expt3
                for kk=1:ops.k1
                    real_values2 = [real_values2, corrs(kk,kk+k)];
                    real_values3 = [real_values3, corrs(kk,kk+2*k)];
                end
            else
                for kk=1:ops.k1
                    real_values2 = [real_values2, corrs(kk,kk+k)];
                end
            end
        elseif ops.k2<ops.k1 
            real_values2 = []; % no third solution
            for kk=1:ops.k2
                row = ops.matchedClusters(kk);
                real_values2 = [real_values2, corrs(row,ops.k1+kk)];
            end
        end

        close all
        h = ERPfigure; set(h,'Position',pos,'visible',ops.isPlotVisible); hold on
        set(gca,'fontsize',16); box off;

        colors = hsv(ops.k1);
        if ops.k2<ops.k1
            colors = colors(ops.matchedClusters,:);
        end

        % histograms
        for i=1:ops.k2

            h(i) = subplot(ops.k2,1,i); hold on
            histogram(null2(i,:),0:.05:1,'FaceColor',colors(i,:),'FaceAlpha',.7,...
                    'EdgeColor',colors(i,:),'EdgeAlpha',.7,'Normalization','probability');
            ylim([0 0.6]); set(gca,'ytick',0:.2:.6);
            
            plot([real_values2(i) real_values2(i)],[0 1],'--','Color',colors(i,:),'linewidth',3);
            text(0.1,0.45,['uncorrupted corr = ' num2str(round(real_values2(i),2))],'Color',colors(i,:),'fontsize',14);
            pvalue = sum(null2(i,:)>real_values2(i))/ops.iterations;
            if pvalue == 0 
                pstring = '< 0.001';
            else 
                pstring = num2str(round(pvalue,5));
            end
            text(0.1,0.38,['p = ' pstring],'Color',colors(i,:),'fontsize',14);

            legend(['Cluster #' num2str(i)],'Location','northwest','fontsize',14,'box','off');
            if i==ops.k2
                xlabel({' ','Correlation'},'fontsize',18,'fontweight','bold');
            end
            if (i==3) && (ops.k2==5)
                ylabel({' ','Probability',' '},'fontsize',18,'fontweight','bold');
            elseif (i==2) && (ops.k2==3)
                ylabel({' ','Probability',' '},'fontsize',18,'fontweight','bold');
            elseif (i==2) && (ops.k2==2)
                ylabel({' ','Probability',' '},'fontsize',18,'fontweight','bold');
            end
        end
        

        % save png
        saveas(gcf,[PLOT_PATH21 ops.solution1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 split_string2 '_K=' num2str(ops.k2) '_null_distribution' trial_string '.png']);

        % save pdf 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH1 ops.solution1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 split_string2 '_K=' num2str(ops.k2) '_null_distribution' trial_string '.pdf'],'pdf')
            
        % save png
        saveas(gcf,[PLOT_PATH22 ops.solution2 split_string2 '_K=' num2str(ops.k2) '_vs_' ops.solution1 '_K=' num2str(ops.k1) '_null_distribution' trial_string '.png']);

        % save pdf 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATH2 ops.solution2 split_string2 '_K=' num2str(ops.k2) '_vs_' ops.solution1 '_K=' num2str(ops.k1) '_null_distribution' trial_string '.pdf'],'pdf')


        % repeat for solution 3
        if expt3
            close all
            h = ERPfigure; set(h,'Position',pos,'visible',ops.isPlotVisible); hold on
            set(gca,'fontsize',16); box off;

            colors = hsv(ops.k1);
            if ops.k2<ops.k1
                colors = colors(ops.matchedClusters,:);
            end

            % histograms
            for i=1:ops.k2

                h(i) = subplot(ops.k2,1,i); hold on
                histogram(null3(i,:),0:.05:1,'FaceColor',colors(i,:),'FaceAlpha',.7,...
                'EdgeColor',colors(i,:),'EdgeAlpha',.7,'Normalization','probability');
                ylim([0 0.6]); set(gca,'ytick',0:.2:.6);
                
                plot([real_values3(i) real_values3(i)],[0 1],'--','Color',colors(i,:),'linewidth',3);
                text(0.1,0.45,['uncorrupted corr = ' num2str(round(real_values3(i),2))],'Color',colors(i,:),'fontsize',14);
                pvalue = sum(null3(i,:)>real_values3(i))/ops.iterations;
                if pvalue == 0 
                    pstring = '< 0.001';
                else 
                    pstring = num2str(round(pvalue,4));
                end
                text(0.1,0.38,['p = ' pstring],'Color',colors(i,:),'fontsize',14);

                legend(['Cluster #' num2str(i)],'Location','northwest','fontsize',14,'box','off');
                if i==ops.k2
                    xlabel({' ','Correlation'},'fontsize',18,'fontweight','bold');
                end
                if (i==3) && (ops.k2==5)
                    ylabel({' ','Probability',' '},'fontsize',18,'fontweight','bold');
                elseif (i==2) && (ops.k2==3)
                    ylabel({' ','Probability',' '},'fontsize',18,'fontweight','bold');
                elseif (i==2) && (ops.k2==2)
                    ylabel({' ','Probability',' '},'fontsize',18,'fontweight','bold');
                end
            end
            

            % save png
            saveas(gcf,[PLOT_PATH21 ops.solution1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 split_string3 '_K=' num2str(ops.k2) '_null_distribution' trial_string '.png']);

            % save pdf 
            set(gcf,'Units','inches');
            screenposition = get(gcf,'Position');
            set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
            saveas(gcf,[PLOT_PATH1 ops.solution1 '_K=' num2str(ops.k1) '_vs_' ops.solution2 split_string3 '_K=' num2str(ops.k2) '_null_distribution' trial_string '.pdf'],'pdf')
                
            % save png
            saveas(gcf,[PLOT_PATH22 ops.solution2 split_string3 '_K=' num2str(ops.k2) '_vs_' ops.solution1 '_K=' num2str(ops.k1) '_null_distribution' trial_string '.png']);

            % save pdf 
            set(gcf,'Units','inches');
            screenposition = get(gcf,'Position');
            set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
            saveas(gcf,[PLOT_PATH2 ops.solution2 split_string3 '_K=' num2str(ops.k2) '_vs_' ops.solution1 '_K=' num2str(ops.k1) '_null_distribution' trial_string '.pdf'],'pdf')

        end

    end


end