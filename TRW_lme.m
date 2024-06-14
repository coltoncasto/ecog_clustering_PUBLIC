%% swjn
% paths
saveName = 'MITSWJNTask_SWJN_kmedoids_correlation_FIGURES';
[~,SAVE_PATH] = initialize(saveName);

randTag = 1;
Dataset = 1;
hemi = 'LH';

conditions = {'S','W','J','N'};
ResultsFolder = [SAVE_PATH 'trw'];
trw_mat = ['MITSWJNTask_langElecs_receptive_window_lengths_words_kernel_gaussian_wide.mat'];
load([ResultsFolder filesep trw_mat])
tablename = 'MITSWJNTask_langElecs_cluster_assignments.csv';
T=readtable([SAVE_PATH filesep 'clustering' filesep tablename]);
T=[T, table(trws)];
T.k3=categorical(T.k3);

if randTag
    formula = 'trws ~ k3 + (k3|subject)';
    nameTag = 'wRand';
else
    formula = 'trws ~ k3';
    nameTag = 'noRand';
end
lme = fitlme(T,formula,'FitMethod','REML');
[beta,betanames,stats] = fixedEffects(lme,'DFMethod','satterthwaite');
Tstats=stats2table(stats);
Tstats.Name = {'TRW C1 - reference';'TRW C2 - relative to C1';'TRW C3 - relative to C1'};

T2 = T;
T2.k3 = reordercats(T2.k3,[2,3,1]);
lme2 = fitlme(T2,formula,'FitMethod','REML');
[beta2,betanames2,stats2] = fixedEffects(lme2,'DFMethod','satterthwaite');
Tstats2=stats2table(stats2);
Tstats2.Name = {'TRW C2 - reference';'TRW C3 - relative to C2';'TRW C1 - relative to C2'};

T3 = T;
T3.k3 = reordercats(T3.k3,[3,1,2]);
lme3 = fitlme(T3,formula,'FitMethod','REML');
[beta3,betanames3,stats3] = fixedEffects(lme3,'DFMethod','satterthwaite');
Tstats3=stats2table(stats3);
Tstats3.Name = {'TRW C3 - reference';'TRW C1 - relative to C3';'TRW C2 - relative to C3'};

writetable([Tstats;Tstats2;Tstats3],[ResultsFolder filesep 'predictTRW_allLMEs_Dataset' num2str(Dataset) '_' nameTag '.csv'])


%% swjn by presentation rate
saveName = 'MITSWJNTask_SWJN_kmedoids_correlation_FIGURES';
[~,SAVE_PATH] = initialize(saveName);

rates = {'450ms','700ms'};
for r=1:length(rates)
    ResultsFolder = [SAVE_PATH 'trw' filesep 'byTime'];
    trw_mat = ['MITSWJNTask_langElecs_receptive_window_lengths_kernel_gaussian_wide_time_' rates{r} '.mat'];
    load([ResultsFolder filesep trw_mat])
    tablename = ['MITSWJNTask_langElecs_receptive_window_lengths_kernel_gaussian_wide_time_' rates{r} '.csv'];
    T=readtable([ResultsFolder filesep tablename]);
    T=[T, table(trws)];
    T.k3=categorical(T.k3);
    
    formula = 'trws ~ k3 + (1 + k3|subject)';

    lme = fitlme(T,formula,'FitMethod','REML');
    [beta,betanames,stats] = fixedEffects(lme,'DFMethod','satterthwaite');
    Tstats=stats2table(stats);
    Tstats.Name = {'TRW C1 - reference';'TRW C2 - relative to C1';'TRW C3 - relative to C1'};
    

    T2 = T;
    T2.k3 = reordercats(T2.k3,[2,3,1]);
    lme2 = fitlme(T2,formula,'FitMethod','REML');
    [beta2,betanames2,stats2] = fixedEffects(lme2,'DFMethod','satterthwaite');
    Tstats2=stats2table(stats2);
    Tstats2.Name = {'TRW C2 - reference';'TRW C3 - relative to C2';'TRW C1 - relative to C2'};

    T3 = T;
    T3.k3 = reordercats(T3.k3,[3,1,2]);
    lme3 = fitlme(T3,formula,'FitMethod','REML');
    [beta3,betanames3,stats3] = fixedEffects(lme3,'DFMethod','satterthwaite');
    Tstats3=stats2table(stats3);
    Tstats3.Name = {'TRW C3 - reference';'TRW C1 - relative to C3';'TRW C2 - relative to C3'};
    writetable([Tstats;Tstats2;Tstats3],[ResultsFolder filesep 'predictTRW_allLMEs_gaussian_wide_Dataset' num2str(Dataset) '_' hemi '_' rates{r} '.csv'])

end


%% langloc 8 words
saveName = 'MITLangloc_SN_8words_kmedoids_correlation_FIGURES';
[~,SAVE_PATH] = initialize(saveName);

Dataset = 2;
stimlen = '8words';
conditions = {'S','N'};
ResultsFolder = [SAVE_PATH 'trw'];
trw_mat = 'MITLangloc_langElecs_receptive_window_lengths_words_kernel_gaussian_wide.mat';
load([ResultsFolder filesep trw_mat])
tablename = 'MITLangloc_langElecs_cluster_assignments.csv';
T=readtable([SAVE_PATH filesep 'clustering' filesep tablename]);
T=[T, table(trws)];
T.k3=categorical(T.k3);

formula = 'trws ~ k3 + (1 + k3|subject)';
lme = fitlme(T,formula,'FitMethod','REML');
[beta,betanames,stats] = fixedEffects(lme,'DFMethod','satterthwaite');
Tstats=stats2table(stats);
Tstats.Name = {'TRW C1 - reference';'TRW C2 - relative to C1';'TRW C3 - relative to C1'};

T2 = T;
T2.k3 = reordercats(T2.k3,[2,3,1]);
lme2 = fitlme(T2,formula,'FitMethod','REML');
[beta2,betanames2,stats2] = fixedEffects(lme2,'DFMethod','satterthwaite');
Tstats2=stats2table(stats2);
Tstats2.Name = {'TRW C2 - reference';'TRW C3 - relative to C2';'TRW C1 - relative to C2'};

T3 = T;
T3.k3 = reordercats(T3.k3,[3,1,2]);
lme3 = fitlme(T3,formula,'FitMethod','REML');
[beta3,betanames3,stats3] = fixedEffects(lme3,'DFMethod','satterthwaite');
Tstats3=stats2table(stats3);
Tstats3.Name = {'TRW C3 - reference';'TRW C1 - relative to C3';'TRW C2 - relative to C3'};

writetable([Tstats;Tstats2;Tstats3],[ResultsFolder filesep 'predictTRW_allLMEs_Dataset' num2str(Dataset) '_' stimlen '.csv'])


%% langloc assignment by correlation
saveName = 'MITLangloc_assignedFrom_MITSWJNTask_SN_kmedoids_correlation_FIGURES';
[~,SAVE_PATH] = initialize(saveName);

Dataset = 2;
conditions = {'S','N'};
ResultsFolder = [SAVE_PATH 'trw'];
trw_mat = 'MITLangloc_langElecs_receptive_window_lengths_words_kernel_gaussian_wide.mat';
load([ResultsFolder filesep trw_mat])
tablename = 'MITLangloc_langElecs_cluster_assignments.csv';
T=readtable([SAVE_PATH filesep 'clustering' filesep tablename]);
T=[T, table(trws)];
T.k3=categorical(T.k3);

formula = 'trws ~ k3 + (1 + k3|subject)';
lme = fitlme(T,formula,'FitMethod','REML');
[beta,betanames,stats] = fixedEffects(lme,'DFMethod','satterthwaite');
Tstats=stats2table(stats);
Tstats.Name = {'TRW C1 - reference';'TRW C2 - relative to C1';'TRW C3 - relative to C1'};

T2 = T;
T2.k3 = reordercats(T2.k3,[2,3,1]);
lme2 = fitlme(T2,formula,'FitMethod','REML');
[beta2,betanames2,stats2] = fixedEffects(lme2,'DFMethod','satterthwaite');
Tstats2=stats2table(stats2);
Tstats2.Name = {'TRW C2 - reference';'TRW C3 - relative to C2';'TRW C1 - relative to C2'};

T3 = T;
T3.k3 = reordercats(T3.k3,[3,1,2]);
lme3 = fitlme(T3,formula,'FitMethod','REML');
[beta3,betanames3,stats3] = fixedEffects(lme3,'DFMethod','satterthwaite');
Tstats3=stats2table(stats3);
Tstats3.Name = {'TRW C3 - reference';'TRW C1 - relative to C3';'TRW C2 - relative to C3'};

writetable([Tstats;Tstats2;Tstats3],[ResultsFolder filesep 'predictTRW_allLMEs_Dataset' num2str(Dataset) '_assignmentByCorrelation.csv'])
