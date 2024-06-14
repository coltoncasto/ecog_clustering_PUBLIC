function calculateReliability(varargin)
% Calculates reliability and compares two groups of electrodes(no need to 
% specify saveName). Assumes that _odd and _even data files are in the 
% output/_reliabilty folder for the given experiments. If they are not, 
% see formatData.m for how to make these files. The second group of 
% channels should probably be the larger set (plotted first).

    p = inputParser();
    addParameter(p,'saveNameAdditional',[]); % additional folder to save output in 
    addParameter(p,'expt1','MITSWJNTask');
    addParameter(p,'expt2','MITSWJNTask');
    addParameter(p,'isPlotVisible',false);
    addParameter(p,'useLangElecs1',true);
    addParameter(p,'useLangElecs2',true);
    addParameter(p,'useWandJ',true) % MITSWJNTask only
    addParameter(p,'signalType','unipolar'); % or 'bipolar'
    parse(p, varargin{:});
    ops = p.Results;

    % --- INITIALIZE --- %

    ops.saveName = '_reliability';

    % paths
    [~,SAVE_PATH] = initialize(ops.saveName);
    PLOT_PATH = [SAVE_PATH 'plots' filesep 'reliability' filesep];
    PLOT_PATH2 = [SAVE_PATH 'plots' filesep 'pngs' filesep 'reliability' filesep];
    if ~exist(PLOT_PATH,'dir'), mkdir(PLOT_PATH); end
    if ~exist(PLOT_PATH2,'dir'), mkdir(PLOT_PATH2); end

    if ops.saveNameAdditional
        [~,SAVE_PATHa] = initialize(ops.saveNameAdditional);
        PLOT_PATHa = [SAVE_PATHa 'plots' filesep 'reliability' filesep];
        PLOT_PATH2a = [SAVE_PATHa 'plots' filesep 'pngs' filesep 'reliability' filesep];
        if ~exist(PLOT_PATHa,'dir'), mkdir(PLOT_PATHa); end
        if ~exist(PLOT_PATH2a,'dir'), mkdir(PLOT_PATH2a); end
    end


    % --- LOAD GROUP #1 --- %

    DATA_PATH = [SAVE_PATH 'data' filesep];
    if ops.useWandJ & strcmp(ops.expt1,'MITSWJNTask')
        cond_string = '_SWJN'; 
    elseif strcmp(ops.expt1,'MITSWJNTask')
        cond_string = '_SN';
    else
        cond_string = '';
    end
    if ops.useLangElecs1, elecType1 = 'langElecs'; else, elecType1 = 'nonLangElecs'; end
    if strcmp(ops.expt1,'both')
        expt_string1 = 'bothMITSWJNTaskandMITLangloc';
    elseif strcmp(ops.expt1,'langloc')
        expt_string1 = 'MGHlangloc';
    else
        expt_string1 = ops.expt1;
    end
    if strcmp(ops.signalType,'bipolar')
        signal_string1 = '_bipolar';
    else
        signal_string1 = '';
    end

    % load data files
    load([DATA_PATH expt_string1 cond_string '_' elecType1 signal_string1 '_data_for_clustering_odd.mat']); % all_X
    X1_odd = all_X;
    load([DATA_PATH expt_string1 cond_string '_' elecType1 signal_string1 '_data_for_clustering_even.mat']); % all_X
    X1_even = all_X;
    assert(size(X1_odd,1)==size(X1_even,1));
    assert(size(X1_odd,2)==size(X1_even,2));

    % load in data labels
    all_X_table1 = readtable([DATA_PATH expt_string1 cond_string '_' elecType1 signal_string1 '_labels_for_clustering_odd.csv']);

    % --- LOAD GROUP #2 --- %

    if ops.useWandJ & strcmp(ops.expt2,'MITSWJNTask')
        cond_string = '_SWJN'; 
    elseif strcmp(ops.expt2,'MITSWJNTask')
        cond_string = '_SN';
    else
        cond_string = '';
    end
    if ops.useLangElecs2, elecType2 = 'langElecs'; else, elecType2 = 'nonLangElecs'; end
    if strcmp(ops.expt2,'both')
        expt_string2 = 'bothMITSWJNTaskandMITLangloc';
    elseif strcmp(ops.expt2,'langloc')
        expt_string2 = 'MGHlangloc';
    else
        expt_string2 = ops.expt2;
    end
    if strcmp(ops.signalType,'bipolar')
        signal_string2 = '_bipolar';
    else
        signal_string2 = '';
    end


    % load data files
    load([DATA_PATH expt_string2 cond_string '_' elecType2 signal_string2 '_data_for_clustering_odd.mat']); % all_X
    X2_odd = all_X;
    load([DATA_PATH expt_string2 cond_string '_' elecType2 signal_string2 '_data_for_clustering_even.mat']); % all_X
    X2_even = all_X;
    assert(size(X2_odd,1)==size(X2_even,1));
    assert(size(X2_odd,2)==size(X2_even,2));

    % load in data labels
    all_X_table2 = readtable([DATA_PATH expt_string2 cond_string '_' elecType2 signal_string2 '_labels_for_clustering_odd.csv']);


    % -------------------
    % --- RELIABILITY --- %
    % -------------------

    % --- GROUP #1 RELIABILITY --- %

    corrs1 = zeros(size(X1_odd,1),1);
    for i=1:size(X1_odd,1)
        corr = corrcoef(X1_odd(i,:),X1_even(i,:)); % 4x4 matrix
        corrs1(i,1) = corr(1,2); % just want corr of A and B
    end

    % save reliabilities
    if ~exist([SAVE_PATH,'reliability'],'dir'), mkdir([SAVE_PATH,'reliability']); end
    corrs = corrs1;
    filename = [SAVE_PATH 'reliability' filesep expt_string1 cond_string '_' elecType1 signal_string1 '_reliability.mat'];
    save(filename,'corrs','-v7.3');

    if ops.saveNameAdditional
        if ~exist([SAVE_PATHa,'reliability'],'dir'), mkdir([SAVE_PATHa,'reliability']); end
        filename = [SAVE_PATHa 'reliability' filesep expt_string1 cond_string '_' elecType1 signal_string1 '_reliability.mat'];
        save(filename,'corrs','-v7.3');
    end

    % --- GROUP #2 RELIABILITY --- %

    corrs2 = zeros(size(X2_odd,1),1);
    for i=1:size(X2_odd,1)
        corr = corrcoef(X2_odd(i,:),X2_even(i,:)); % 4x4 matrix
        corrs2(i,1) = corr(1,2); % just want corr of A and B
    end

    % save reliabilities
    corrs = corrs2;
    filename = [SAVE_PATH 'reliability' filesep expt_string2 cond_string '_' elecType2 signal_string2 '_reliability.mat'];
    save(filename,'corrs','-v7.3');

    if ops.saveNameAdditional
        if ~exist([SAVE_PATHa,'reliability'],'dir'), mkdir([SAVE_PATHa,'reliability']); end
        filename = [SAVE_PATHa 'reliability' filesep expt_string2 cond_string '_' elecType2 signal_string2 '_reliability.mat'];
        save(filename,'corrs','-v7.3');
    end

    % --- PLOTTING DISTRIBUTIONS --- %

    close all
    h = ERPfigure; set(h,'Position',[0 0 800 700],'visible',ops.isPlotVisible)

    % histogram 
    color2 = [0 0 0]; color1 = [0.8 0.2 0.2];
    histogram(corrs2,20,'FaceColor',color2,'FaceAlpha',0.4); hold on;
    histogram(corrs1,20,'FaceColor',color1,'FaceAlpha',0.4); 
    set(gca,'fontsize',14); box off; xlim([-0.5 1]);

    % dashed lines with mean corr
    yl = ylim;
    mean2 = mean(corrs2,1);
    plot([mean2 mean2],[yl(1) yl(2)],'--','Color',color2,'linewidth',3);
    mean1 = mean(corrs1,1);
    plot([mean1 mean1],[yl(1) yl(2)],'--','Color',color1,'linewidth',3);

    % label means
    text(mean2+0.03,yl(2)*0.8,['mean = ' num2str(round(mean2,2))],'Color',color2,'fontsize',11);
    text(mean1+0.03,yl(2)*0.65,['mean = ' num2str(round(mean1,2))],'Color',color1,'fontsize',11);

    % other plotting params
    xlabel({' ','Correlation'},'fontsize',18,'fontweight','bold');
    ylabel({'Channels',' '},'fontsize',18,'fontweight','bold');
    legend({[expt_string2 ' ' elecType2],[expt_string1 ' ' elecType1]},'Location','northeast','fontsize',20,'box','off');
    
    % save png 
    saveas(gcf,[PLOT_PATH2 expt_string1 '_' elecType1 signal_string1 '_' expt_string2 '_' elecType2 signal_string2 cond_string '_odd_even_split.png']);

    % save pdf 
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    saveas(gcf,[PLOT_PATH expt_string1 '_' elecType1 signal_string1 '_' expt_string2 '_' elecType2 signal_string2 cond_string '_odd_even_split.pdf'],'pdf')

    if ops.saveNameAdditional

        % save png 
        saveas(gcf,[PLOT_PATH2a expt_string1 '_' elecType1 signal_string1 '_' expt_string2 '_' elecType2 signal_string2 cond_string '_odd_even_split.png']);

        % save pdf 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        saveas(gcf,[PLOT_PATHa expt_string1 '_' elecType1 signal_string1 '_' expt_string2 '_' elecType2 signal_string2 cond_string '_odd_even_split.pdf'],'pdf')

    end

end