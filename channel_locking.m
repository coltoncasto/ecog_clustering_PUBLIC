function channel_locking(varargin)

    p = inputParser();
    addRequired(p,'saveName');
    addParameter(p,'experiment','MITSWJNTask'); % 'MITSWJNTask' or 'MITLangloc' or 'both'
    addParameter(p,'k',3);
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs',true);
    addParameter(p,'useWandJ',false); % MITSWJNTask only
    addParameter(p,'calcFlag',0); % run fitting of sine
    addParameter(p,'srate',60);
    parse(p, varargin{:});
    ops = p.Results;

    %% INITIALIZE
    [CLUSTER_PATH,SAVE_PATH] = initialize(ops.saveName);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'locking' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'locking' filesep];
    if ~exist(PLOT_PATH,'dir')
        mkdir(PLOT_PATH);
    end
    if ~exist(PLOT_PATH2,'dir')
        mkdir(PLOT_PATH2);
    end

    % experiment-specific info
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
    if ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask')
        conditions = {'S','W','J','N'};
        nConds = 4;
        nWords = 8;
        cond_string = '_SWJN';
        colors = {[0.15 0.15 0.15],[0.35 0.35 0.35],[0.55 0.55 0.55],[0.75 0.75 0.75]};
    elseif ~ops.useWandJ && strcmp(ops.experiment,'MITSWJNTask') % only S and N
        conditions = {'S','N'};
        nConds = 2;
        nWords = 8;
        cond_string = '_SN';
        colors = {[0.15 0.15 0.15],[0.75 0.75 0.75]};
    else % MITLangloc
        conditions = {'S','N'};
        nConds = 2;
        nWords = 8;
        cond_string = '';
        colors = {[0.15 0.15 0.15],[0.75 0.75 0.75]};
    end

    % load in data matrix
    DATA_PATH = [SAVE_PATH 'data' filesep];
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering.mat']); % all_X
    X = all_X;

    % load in all trial data matrix
    load([DATA_PATH expt_string '_' elecType '_data_for_clustering_all_trials.mat']); % all_X
    X_all = all_X;

    % load in odd and even data matrix
    stab_dir = split(ops.saveName,'_');
    stab_dir = join(stab_dir(1:end-1),'_');
    stab_dir = [stab_dir{1} '_STABILITY_FIGURES'];
    DATA_PATH_STAB = [CLUSTER_PATH 'output' filesep stab_dir filesep 'data' filesep];
    load([DATA_PATH_STAB expt_string '_' elecType '_data_for_clustering_odd.mat']); % all_X
    Xodd = all_X;
    load([DATA_PATH_STAB expt_string '_' elecType '_data_for_clustering_even.mat']); % all_X
    Xeven = all_X;

    % load in data labels
    Tl = readtable([SAVE_PATH 'clustering' filesep expt_string '_' elecType '_cluster_assignments.csv']);

    % load in all trials data labels
    Tl_all = readtable([DATA_PATH expt_string '_' elecType '_labels_for_clustering_all_trials.csv'],'Delimiter',',');

    %% run calculation
    if ops.calcFlag
        fs=ops.srate;
        [nChan,nSamp] = size(X);
        condLen = nSamp/nConds;
        stimLen = condLen/nWords;
        freq_stim_samp = 1/stimLen;


        condLen_s = condLen/fs;
        stimLen_s = condLen_s/nWords;
        freq_stim = 1/stimLen_s;

        chanchange=(find(diff(Tl_all.channel_number)));

        rhos = nan(nChan,nConds);%channels x conditions
        ps = nan(nChan,nConds);%channels x conditions

        fromrow=1;
        for ichan = 1:nChan

            if ichan<nChan
                torow = chanchange(ichan);
            else
                torow = height(Tl_all);
            end
            Tl_all_chan = Tl_all(fromrow:torow,:);
            ntrials = height(Tl_all_chan);

            for icond = 1:nConds
                ticnow=tic;
                condition = conditions{icond};
                data = X_all(fromrow:torow,1+condLen*(icond-1):condLen*(icond));
                data_mean = mean(data);
                data_odd = data(1:2:end,:);
                data_even = data(2:2:end,:);

                %fit a sine wave onto the odd trials:
                mean_odd = mean(data_odd);
                mean_even = mean(data_even);   

                t=0:1/fs:(length(data_mean)-1)/fs;
                sinewave = 0.5*sin(2*pi*freq_stim_samp*fs*t);

                amplitude=0.5;freq=freq_stim;phase=pi;offset=amplitude/2;

                initialparameter=[amplitude,phase,offset];
                mx=@(initialparameter)fita(initialparameter,freq,t,mean_odd);
                outputparameters=fminsearch(mx,initialparameter);

                sine=outputparameters(1)*sin(2*pi*freq*t+outputparameters(2))+outputparameters(3);
                [rho,p]=corr(mean_even',sine');

                    %plot(x,even,'.','linewidth',200);
                    %hold on
                    %plot(x,sine,'r');
                    %axis([min(x) max(x) min(even) max(even)]);
                    %title(['Corr Rho=' num2str(rho) ' p=' num2str(p)],'fontsize',20)

                rhos(ichan,icond)=rho;
                ps(ichan,icond)=p;
                disp(['Done iChan= ' num2str(ichan) ', icond= ' num2str(icond) ' in ' num2str(toc(ticnow)) ' sec']);
            end
            fromrow=torow+1;
        end
        % save
        LOCK_PATH = [SAVE_PATH 'locking' filesep];
        if ~exist(LOCK_PATH,'dir'), mkdir(LOCK_PATH); end
        filename = [LOCK_PATH expt_string '_' elecType '_locking_fitted_sine_values.mat'];
        save(filename,'rhos','ps');
        
    end
    
    %% load locking
    LOCK_PATH = [SAVE_PATH 'locking' filesep];
    filename = [LOCK_PATH expt_string '_' elecType '_locking_fitted_sine_values.mat'];
    load(filename); % rhos, ps

    %% save as a table
    if nConds==2
        T1=table;
        locking_S = rhos(:,1);
        locking_N = rhos(:,4);
        T1=[table(locking_S),table(locking_N)];
        Tl = [Tl,T1];
        writetable(Tl,[LOCK_PATH expt_string '_' elecType '_chan_assignments_wLocking.csv'])
        names = {'locking_S','locking_N'};
    elseif nConds==4
        T1=table;
        locking_S = rhos(:,1);
        locking_W = rhos(:,2);
        locking_J = rhos(:,3);
        locking_N = rhos(:,4);
        T1=[table(locking_S),table(locking_W),table(locking_J),table(locking_N)];
        Tl = [Tl,T1];
        writetable(Tl,[LOCK_PATH expt_string '_' elecType '_chan_assignments_wLocking.csv']) 
        names = {'locking_S','locking_W','locking_J','locking_N'};
    end

    %% plot per cluster
    clusters = ops.k;
    means = nan(clusters,length(conditions));
    errors = nan(clusters,length(conditions));
    Cscatter = cell(clusters,length(conditions));
    for iclust = 1:clusters
        for icond = 1:length(conditions)
            eval(strcat('rows = find(Tl.k',num2str(clusters),'==iclust);'));
            means(iclust,icond) = mean(Tl(rows,:).(names{icond}));
            errors(iclust,icond) = std(Tl(rows,:).(names{icond}))/sqrt(length(rows));
            Cscatter{iclust,icond} = Tl(rows,:).(names{icond});
        end
    end

    figure('Position',[100 100 600 500],'visible',ops.isPlotVisible)
    if (nConds==4) && (clusters==3)
        xs = [1 2 3 4 6 7 8 9 11 12 13 14];
    elseif (nConds==2) && (clusters==3)
        xs = [1 2 4 5 7 8];
    end

    ii=0;
    for i=1:size(means,1)
        for j=1:size(means,2)
            ii=ii+1;
            h(i,j) = bar(xs(ii),means(i,j),'DisplayName',conditions{j});
            set(h(i,j),'FaceColor',colors{j});
            hold on
            if j==1
                swarmchart(ones(length(Cscatter{i,j}),1)*xs(ii),Cscatter{i,j},20,'k','filled','MarkerEdgeColor',[0.3,0.3,0.3],'XJitter','density','XJitterWidth',0.3,'HandleVisibility','off','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
            else
                swarmchart(ones(length(Cscatter{i,j}),1)*xs(ii),Cscatter{i,j},20,colors{j},'filled','MarkerEdgeColor','k','XJitter','density','XJitterWidth',0.3,'HandleVisibility','off','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
            end
            if j==1
                errorbar(xs(ii),means(i,j),errors(i,j),'color',[0.3 0.3 0.3],'LineWidth',2,'CapSize',8,'HandleVisibility','off');
            else
                errorbar(xs(ii),means(i,j),errors(i,j),'k','LineWidth',2,'CapSize',8,'HandleVisibility','off');
            end
        end
    end

    set(gca,'fontsize',20)
    set(gca,'ytick',[-0.5 0 0.5 1])
    set(gca,'xtick',[2.5 7.5 12.5],'xticklabels',{'1','2','3'})

    ylabel('Correlation to sinusoidal stimulus train','FontSize',20)
    xlabel('Cluster','FontSize',16)
    legend(conditions,'location','nw')
    title('Locking of channels to stimulus')
    
    % save png
    saveas(gcf,[PLOT_PATH2 expt_string '_' elecType '_locking.png']);

    % save pdf
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH expt_string '_' elecType '_locking.pdf'],'pdf');
    
    
    %% stats -- ONLY WORKS FOR K==3 and SWJN 4 CONDITIONS CURRENTLY
    if (nConds==4) && (ops.k==3)
        Tlme = table;
        i=0;
        for ii=1:height(Tl)
            for icond=1:nConds
                i=i+1;
                temp=table;
                temp.subject = Tl.subject(ii);
                temp.channel_number = Tl.channel_number(ii);
                temp.channel_name = Tl.channel_name(ii);
                temp.k3 = Tl.k3(ii);
                temp.condition = conditions{icond};
                temp.locking = Tl.(names{icond})(ii);
                Tlme=[Tlme;temp];
                clear temp
            end
        end
        Tlme.k3=categorical(Tlme.k3);

        formula = 'locking ~ k3*condition + (k3|subject) + (condition|subject)';

        lme_wint = fitlme(Tlme,formula,'FitMethod','REML');
        ANOVA = anova(lme_wint,'DFMethod','satterthwaite');
        Tanova=anova2table(ANOVA);
        writetable(Tanova,[LOCK_PATH expt_string '_' elecType '_ANOVA_locking_withInt.csv']);

        [beta,betanames,stats] = fixedEffects(lme_wint,'DFMethod','satterthwaite');
        Tstats=stats2table(stats);
        Tstats.Name = {'Intercept: C1, Condition S';'C2';'C3';'Condition W';'Condition J';'Condition N';'C2 * Condition W';'C3 * Condition W';'C2 * Condition J';'C3 * Condition J';'C2 * Condition N';'C3 * Condition N'};
        writetable([Tstats],[LOCK_PATH expt_string '_' elecType '_Stats_locking_withInt.csv']);


        Tlme2 = Tlme;
        Tlme2.k3 = reordercats(Tlme2.k3,[2,3,1]);
        lme2_wint = fitlme(Tlme2,formula,'FitMethod','REML');
        ANOVA = anova(lme2_wint,'DFMethod','satterthwaite')
        [beta2,betanames2,stats2] = fixedEffects(lme2_wint,'DFMethod','satterthwaite');
        Tstats2=stats2table(stats2);
        Tstats2.Name = {'Intercept: C2, Condition S';'C3';'C1';'Condition W';'Condition J';'Condition N';'C3 * Condition W';'C1 * Condition W';'C3 * Condition J';'C1 * Condition J';'C3 * Condition N';'C1 * Condition N'};

        Tlme3 = Tlme;
        Tlme3.k3 = reordercats(Tlme3.k3,[3,1,2]);
        lme3_wint = fitlme(Tlme3,formula,'FitMethod','REML');
        ANOVA = anova(lme3_wint,'DFMethod','satterthwaite')
        [beta3,betanames3,stats3] = fixedEffects(lme3_wint,'DFMethod','satterthwaite');
        Tstats3=stats2table(stats3);
        Tstats3.Name = {'Intercept: C3, Condition S';'C1';'C2';'Condition W';'Condition J';'Condition N';'C1 * Condition W';'C2 * Condition W';'C1 * Condition J';'C2 * Condition J';'C1 * Condition N';'C2 * Condition N'};

        TstatsAll = [Tstats; Tstats2; Tstats3];
        writetable([TstatsAll],[LOCK_PATH expt_string '_' elecType '_StatsAll_locking_withInt.csv']);
        
    end

end