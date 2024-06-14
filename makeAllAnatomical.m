function makeAllAnatomical(varargin)
% makes images or gifs of all electrodes across subjects

% ---------------
%% --- SETUP ----
% ---------------
p = inputParser();
addRequired(p,'saveName');
addParameter(p,'experiment','MITSWJNTask');
addParameter(p,'k',3); % clustering solution for k clusters
addParameter(p,'doGIF',false);
addParameter(p,'angle',270);
addParameter(p,'alpha',0.6);
addParameter(p,'elecSize',100);
addParameter(p,'color','r');
addParameter(p,'doLang',false);
addParameter(p,'doBoth',false);
addParameter(p,'doCluster',[]); % value of k 
addParameter(p,'excludeCluster',[]);
addParameter(p,'useWandJ',true);
addParameter(p,'shadeByReliability',false); % set alpha by split half reliability
addParameter(p,'shadeByClusterDistance',false); % set alpha by dist to cluster center
addParameter(p,'doTRW',false);
addParameter(p,'doClusterDistance',false); % plot colors due to cluster distance. Make 3 plots one for each cluster (k=3)
addParameter(p,'doPartialCorrelations',false); % plot colors due to *partial correlation* cluster distance. Make 3 plots one for each cluster (k=3)
addParameter(p,'AllChans',true); % only relevant if 'doClusterDistance/partialCorrelations'=true. If false then plots just assigned electrodes by cluster. If true plots all
addParameter(p,'colormapPerCluster',false); % only relevant if 'doClusterDistance'=true or 'doPartialCorrelations'=true
addParameter(p,'isPlotVisible',false);
addParameter(p,'whichKernel','gaussian_wide'); % alternatives: 'square', 'cosine'
addParameter(p,'colors',[]); % when specified will only work with one value of k
parse(p, varargin{:});
ops = p.Results;

% experiments
expt = ops.experiment;
if strcmp(expt,'both')
    nSummaries = 2;
    expt = {'MITSWJNTask','MITLangloc'};
    expt_suffix = 'bothMITSWJNTaskandMITLangloc';
else
    nSummaries = 1;
    expt_suffix = expt;
    expt = {expt};
end

% viewing angle(s)
if ops.doGIF
    n = ops.angle:1:ops.angle+360;
else
    n = ops.angle;
end

% colormap
remap = [];
cmap = 'Gray';

% update lang to be true if 'doBoth' specified
if ops.doBoth
    ops.doLang = true;
end

% initalize paths
[CLUSTER_PATH,SAVE_PATH] = initialize(ops.saveName);
CLUSTERING_PATH = [SAVE_PATH 'clustering' filesep];


% electrode type for clustering
if ops.doLang
    elec_type_cluster = 'langElecs';
elseif ops.doCluster
    elec_type_cluster = 'nonLangElecs';
end

all_X_table = readtable([SAVE_PATH 'clustering' filesep expt_suffix '_' elec_type_cluster '_cluster_assignments.csv']);

% load cluster assignments 
if ops.doCluster
    eval(strcat("assignments=all_X_table.k",num2str(ops.doCluster),";"));

    % load in data (might go unused)
    DATA_PATH = [SAVE_PATH 'data' filesep];
    load([DATA_PATH expt_suffix '_' elec_type_cluster '_data_for_clustering.mat']); % all_X
    X = all_X;

    % colors for cluster 
    if ~isempty(ops.colors)
        colors = ops.colors;
    else
        colors = hsv(ops.doCluster);
    end
    shapes = {'o','s','d','<','>','p'};
    cluster_legend_names = cell(1,ops.doCluster);
    for k=1:ops.doCluster
        cluster_legend_names{k} = ['Cluster #' num2str(k)];
    end
    if ops.excludeCluster
        cluster_legend_names(ops.excludeCluster) = [];
    end

end

% load trw values
if ops.doTRW
    load([SAVE_PATH 'trw' filesep expt_suffix '_' elec_type_cluster '_receptive_window_lengths_words_kernel_' ops.whichKernel '.mat']); %trws
    trws = uint8(round(trws,1)*10); % because double comparisons in matlab cause issues
    reference = min(trws):1:max(trws);
    color_map = spring(length(reference));
end

if ops.doClusterDistance || ops.doPartialCorrelations 
    nplots = 3;
    % get average cluster response
    eval(strcat("assignments=all_X_table.k3;"));
    
    % load in data (might go unused)
    DATA_PATH = [SAVE_PATH 'data' filesep];
    load([DATA_PATH expt_suffix '_' elec_type_cluster '_data_for_clustering.mat']); % all_X
    X = all_X;

    IDX = assignments; Cm = zeros(3,size(X,2));
    for kk=1:3, Cm(kk,:) = mean(X(IDX==kk,:),1); end
   
    %load medoids
    load([CLUSTERING_PATH filesep expt_suffix '_' elec_type_cluster '_clusters_K=' num2str(ops.k) '.mat'])%loaded as C
    medoids = C;
else
    nplots = 1;
end

% alpha values 
if ops.shadeByReliability
    RELIAB_PATH = [CLUSTER_PATH 'output/_reliability/reliability/'];
    if ops.useWandJ & strcmp(ops.experiment,'MITSWJNTask')
        cond_string = 'SWJN';
    else % only 2 conds
        cond_string = 'SN';
    end
    load([RELIAB_PATH expt_suffix '_' cond_string '_' elec_type_cluster '_reliability.mat']); % corrs
    %corrs = abs(corrs);
    corrs(corrs<0) = 0; 
    corrs = corrs/max(corrs);
    corrs = corrs*0.7 + 0.3;
    alpha_string = '_shaded_by_reliability';

elseif ops.shadeByClusterDistance
    % get average cluster response
    eval(strcat("assignments=all_X_table.k",num2str(ops.doCluster),";"));
    IDX = assignments; C = zeros(ops.doCluster,size(X,2));
    for kk=1:ops.doCluster, C(kk,:) = mean(X(IDX==kk,:),1); end

    % calculate distance from every chan to cluster center (correlation distance)
    corrs = nan(size(all_X_table,1),1);
    for c = 1:size(corrs,1)
        corrs(c,1) = 1 - pdist2(X(c,:),C(IDX(c),:),'correlation');
    end
    corrs = abs(corrs);%TODO = zero/keep the negative values? 
    alpha_string = '_shaded_by_cluster_distance';

else
    corrs = ones(size(all_X_table,1),1);
    alpha_string = '';
end

% ------------------
%% --- PLOTTING ----
% ------------------
%close all

for ik=1:nplots

    fig = figure('Position',[0 0 1000 700],'visible',ops.isPlotVisible);
    
    if ops.colormapPerCluster
        if ops.doClusterDistance
            reference = min(dist_corrs(assignments==ik)):0.01:max(dist_corrs(assignments==ik));
            color_map = winter(length(reference));
        elseif ops.doPartialCorrelations
            reference = -0.5:0.01:1;
            color_map = cool(length(reference));
        end
    end


    % store images at each viewing angle
    nAngles = length(n);
    for i = 1:nAngles
        
        % reset cluster legend handles
        if ops.doCluster
            cluster_legend_handles = cell(1,ops.doCluster);
        end
        
        % --- MNI BRAIN ---
        filename = ['anatomical' filesep 'cvs_avg35_inMNI152.mat'];
        MNI_mat = load(filename);
        plot3DModel(gca,MNI_mat.cortex,remap,ops.alpha);
        colormap(cmap)
        hold on
        
        
        % --- ELECTRODES FROM ALL SUBJECTS ---
        % go through experiments
        for j=1:nSummaries
            
            % load summary statistics file
            load(['anatomical' filesep 'summary_statistics_' expt{j} '.mat']); % all_summary_stats
            
            % go through subjects in experiment
            subjects = all_summary_stats.subject;
            for is=1:length(subjects)
                
                if strcmp(subjects{is},'MCJ014')
                    continue
                end
                
                filename = ['anatomical' filesep subjects{is} '_MNI_brain.mat'];
                load(filename); % vera_mat_minimal
                vera_mat = vera_mat_minimal;
                
                [chans,~,langChans] = channelMappingAnatomical(subjects{is},...
                                                             all_summary_stats,...
                                                             vera_mat,...
                                                             CLUSTER_PATH);
    
                if ops.doBoth % do both lang AND all electrodes
                    
                    chans_use = cell2mat(chans);
                    scatter3(vera_mat.tala.electrodes(chans_use,1),...
                        vera_mat.tala.electrodes(chans_use,2),...
                        vera_mat.tala.electrodes(chans_use,3),...
                        ops.elecSize,'k','filled');
                    
                    scatter3(vera_mat.tala.electrodes(langChans,1),...
                        vera_mat.tala.electrodes(langChans,2),...
                        vera_mat.tala.electrodes(langChans,3),...
                        ops.elecSize,'MarkerEdgeColor','k',...
                        'MarkerFaceColor','r');
    
                elseif ops.doLang & isempty(ops.doCluster) & ~ops.doTRW & ~ops.doClusterDistance & ~ops.doPartialCorrelations
                    
                    scatter3(vera_mat.tala.electrodes(langChans,1),...
                        vera_mat.tala.electrodes(langChans,2),...
                        vera_mat.tala.electrodes(langChans,3),...
                        ops.elecSize,'MarkerEdgeColor','k',...
                        'MarkerFaceColor',ops.color);
                
                elseif ops.doCluster & ~ops.doTRW & ~ops.doClusterDistance & ~ops.doPartialCorrelations % plot cluster assignments for specified k
    
                    subject_idxs = cell2mat(cellfun(@(x) strcmp(subjects{is},x),all_X_table.subject,'UniformOutput',false));
                    subject_assignments = assignments(subject_idxs);
                    subject_corrs = corrs(subject_idxs);
                    assert(length(subject_assignments)==length(langChans));
    
                    for m=1:length(subject_assignments)
                        cluster = subject_assignments(m);
                        alpha = subject_corrs(m);
                        if ops.excludeCluster & (cluster==ops.excludeCluster)
                            continue
                        else
                            h = scatter3(vera_mat.tala.electrodes(langChans(m),1),...
                                vera_mat.tala.electrodes(langChans(m),2),...
                                vera_mat.tala.electrodes(langChans(m),3),...
                                ops.elecSize,'o','MarkerEdgeColor','k',...
                                'MarkerFaceColor',colors(cluster,:),...
                                'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
                            
                            % fill in legned handles until full
                            if isempty(cluster_legend_handles{cluster})
                                cluster_legend_handles{cluster} = h;
                            end
                        end
                    end
                    title({'All Subjects',' ',' '},'fontsize',20);
                    
                elseif ops.doTRW
    
                    subject_idxs = cell2mat(cellfun(@(x) strcmp(subjects{is},x),all_X_table.subject,'UniformOutput',false));
                    subject_trws = trws(subject_idxs);
                    assert(length(subject_trws)==length(langChans));%
                    subject_corrs = corrs(subject_idxs);%for transparency
    
                    for m=1:length(subject_trws)
                        trw = subject_trws(m);
                        alpha = subject_corrs(m);
                        h = scatter3(vera_mat.tala.electrodes(langChans(m),1),...
                                vera_mat.tala.electrodes(langChans(m),2),...
                                vera_mat.tala.electrodes(langChans(m),3),...
                                ops.elecSize,'o','MarkerEdgeColor','k',...
                                'MarkerFaceColor',color_map(reference==trw,:),...
                                'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
                    end
                    title({'All Subjects','Estimated TRWs',' '},'fontsize',20);

                elseif ops.doClusterDistance

                    % calculate distance from every chan to cluster center (correlation distance)
                    dist_corrs = nan(size(all_X_table,1),1);
                    for c = 1:size(dist_corrs,1)
                        dist_corrs(c,1) = 1 - pdist2(X(c,:),C(ik,:),'correlation');%distance from cluster ik
                    end
                    dist_corrs = abs(dist_corrs);
                   
                    reference = linspace(min(dist_corrs),max(dist_corrs),100);
                    color_map = cool(length(reference));


                    subject_idxs = cell2mat(cellfun(@(x) strcmp(subjects{is},x),all_X_table.subject,'UniformOutput',false));
                    if ops.AllChans
                        subject_cluster_idx = subject_idxs;
                    else
                        subject_cluster_idx = subject_idxs & assignments==ik;
                    end

                    subject_dist_corrs = dist_corrs(subject_cluster_idx);
                    langChans_cluster = langChans(subject_cluster_idx(subject_idxs));
                    subject_corrs = corrs(subject_cluster_idx);%for transparency

                    for m=1:length(subject_dist_corrs)
                        D = subject_dist_corrs(m);
                        [~, idx] = min(abs(reference-D));

                        alpha = subject_corrs(m);

                        h = scatter3(vera_mat.tala.electrodes(langChans_cluster(m),1),...
                                vera_mat.tala.electrodes(langChans_cluster(m),2),...
                                vera_mat.tala.electrodes(langChans_cluster(m),3),...
                                ops.elecSize,'o','MarkerEdgeColor','k',...
                                'MarkerFaceColor',color_map(idx,:),...
                                'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
                    end
                    title({'All Subjects',['Cluster ' num2str(ik)],'Correlation with medoid'},'fontsize',20);

                elseif ops.doPartialCorrelations
                  
                    % partial correlation
                    k_idxs = zeros(ops.k,1); k_idxs(ik) = 1;
                    k_idxs = logical(k_idxs);
                    pcs = partialcorr([C(k_idxs,:)', X'],C(~k_idxs,:)');
                    pcs = pcs(2:end,1);
                      
                     reference = linspace(min(pcs),max(pcs),100);
                    color_map = cool(length(reference));

                    subject_idxs = cell2mat(cellfun(@(x) strcmp(subjects{is},x),all_X_table.subject,'UniformOutput',false));
                    if ops.AllChans
                        subject_cluster_idx = subject_idxs;
                    else
                        subject_cluster_idx = subject_idxs & assignments==ik;
                    end

                    subject_dist_corrs = pcs(subject_cluster_idx);
                    langChans_cluster = langChans(subject_cluster_idx(subject_idxs));
                    subject_corrs = corrs(subject_cluster_idx);%for transparency

                    for m=1:length(subject_dist_corrs)
                        D = subject_dist_corrs(m);
                        [~, idx] = min(abs(reference-D));

                        alpha = subject_corrs(m);

                        h = scatter3(vera_mat.tala.electrodes(langChans_cluster(m),1),...
                                vera_mat.tala.electrodes(langChans_cluster(m),2),...
                                vera_mat.tala.electrodes(langChans_cluster(m),3),...
                                ops.elecSize,'o','MarkerEdgeColor','k',...
                                'MarkerFaceColor',color_map(idx,:),...
                                'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
                    end
                    title({'All Subjects',['Cluster ' num2str(ik)],'Partial correlation with medoid'},'fontsize',20);

                end
            end
        end
    
        % add legend
        if ops.doCluster
            string_for_legend = '[';
            for idx=1:length(cluster_legend_handles)
                string_for_legend = [string_for_legend 'cluster_legend_handles{' num2str(idx) '} '];
            end
            string_for_legend = [string_for_legend(1:end-1) ']'];
            eval(strcat("legend(",string_for_legend,",cluster_legend_names,'box','off',", ...
                                "'location','eastoutside','fontsize',26);"));
        end
        
        
        view(n(i),0);
        set(fig,'Color','w');
    
        frame = getframe(fig);
        im{i} = frame2im(frame);
        % clf;
        if ops.doTRW
            ax_new = axes('visible','off');
            colormap(ax_new,spring)
            hc = colorbar;
            trw_map = double(reference)./10;
            ticks = find(ismember(trw_map,[1 2 3 4 5 6 7 8]))./numel(trw_map);
            set(hc,'position',[0.9, 0.1, 0.02,0.8],'fontsize',16,'Ticks',ticks,'TickLabels',{'1','2','3','4','5','6','7','8'})
            title(hc,'TRW (s)','fontsize',20)
            %set(hc,'position',[0.9, 0.1, 0.02,0.8],'fontsize',20,'Ticks',[0,0.2,0.4,0.6,0.8,1],'TickLabels',strsplit(num2str(linspace(trw_map(1),trw_map(end),6))))
        
        elseif ops.doClusterDistance
            ax_new = axes('visible','off');
            colormap(ax_new,cool)
            hc = colorbar;
            dist_map = reference;
            tick_values = [min(reference), max(reference)];
            tick_idx = nan(size(tick_values));
            for it = 1:length(tick_values)
                [~, tick_idx(it)] = min(abs(reference-tick_values(it)));
            end
            tick_idx = tick_idx/numel(reference);
            set(hc,'position',[0.9, 0.1, 0.02,0.8],'fontsize',16,'Ticks',tick_idx,'TickLabels',strsplit(num2str(tick_values,2)))
            title(hc,{'Corr. with medoid'},'fontsize',20)
            %set(hc,'position',[0.9, 0.1, 0.02,0.8],'fontsize',20,'Ticks',[0,0.2,0.4,0.6,0.8,1],'TickLabels',strsplit(num2str(linspace(trw_map(1),trw_map(end),6))))
   
        elseif ops.doPartialCorrelations
            ax_new = axes('visible','off');
            colormap(ax_new,cool)
            hc = colorbar;
            dist_map = reference;
            tick_values = [min(reference), max(reference)];
            tick_idx = nan(size(tick_values));
            for it = 1:length(tick_values)
                [~, tick_idx(it)] = min(abs(reference-tick_values(it)));
            end
            tick_idx = tick_idx/numel(reference);
            set(hc,'position',[0.9, 0.1, 0.02,0.8],'fontsize',16,'Ticks',tick_idx,'TickLabels',strsplit(num2str(tick_values,2)))
            title(hc,{'Partial Corr.','with medoid'},'fontsize',20)
            %set(hc,'position',[0.9, 0.1, 0.02,0.8],'fontsize',20,'Ticks',[0,0.2,0.4,0.6,0.8,1],'TickLabels',strsplit(num2str(linspace(trw_map(1),trw_map(end),6))))
   
        end
    end
    
    %% SAVE
    % output directory
    saveDir = [SAVE_PATH filesep 'plots' filesep 'anatomical' filesep];
    saveDir2 = [SAVE_PATH filesep 'plots' filesep 'pngs' filesep 'anatomical' filesep];
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    if ~exist(saveDir2, 'dir')
        mkdir(saveDir2);
    end
    
    % componenets of output filename
    % GIF vs PDF
    if ops.doGIF
        file_suffix = '.gif';
    else
        file_suffix = ['_angle' num2str(ops.angle)];
    end
    % electrode type
    if ~ops.doBoth && ops.doLang 
        elec_type = '_langElecs';
    elseif ops.doCluster 
        elec_type = ['_' elec_type_cluster];
    else
        elec_type = '_allElecs';
    end
    % clustering info
    if ops.doCluster
        cluster_string = ['_K=' num2str(ops.doCluster) ''];
    else
        cluster_string = '';
    end
    % cluster to exclude
    if ops.excludeCluster
        exclude_string = ['_excludedCluster' num2str(ops.excludeCluster)];
    else
        exclude_string= '';
    end
    % trw info
    if ops.doTRW
        trw_string = '_coloredByTRW';
    else
        trw_string = '';
    end
    % clustDist info
    if ops.doClusterDistance
        if ops.AllChans
            dist_string = ['_dist2MedoidCluster' num2str(ik) '_allchans'];
        else
            dist_string = ['_dist2MedoidCluster' num2str(ik)];
        end
    else
        dist_string = '';
    end
    if ops.doPartialCorrelations
        if ops.AllChans
            part_string = ['_partCorr2MedoidCluster' num2str(ik) '_allchans'];
        else
            part_string = ['_partCorr2MedoidCluster' num2str(ik)];
        end
    else
        part_string = '';
    end
    filename = [saveDir expt_suffix elec_type cluster_string exclude_string trw_string dist_string part_string alpha_string file_suffix];
    filename2 = [saveDir2 expt_suffix elec_type cluster_string exclude_string trw_string dist_string part_string alpha_string file_suffix];
    
    % save images to gif
    if ops.doGIF
        for i = 1:nAngles
            [A,map] = rgb2ind(im{i},256);
            if i == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',.05);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.05);
            end
        end
        
    % save single image
    else
    
        % png
        assert(length(im)==1,'Images created from more than one viewing angle.');
        [A,map] = rgb2ind(im{i},256);
        imwrite(A,map,[filename2 '.png']);
        
        % pdf
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[filename '.pdf'],'pdf');
        
    end

end
%close all


end

