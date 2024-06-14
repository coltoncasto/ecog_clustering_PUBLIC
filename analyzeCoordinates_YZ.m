% currently does not work as a general function

%% Definitions
DataFolder = 'output/_coordinates';
ResultsFolder = [DataFolder '/results/'];
if ~exist(ResultsFolder,'dir'), mkdir(ResultsFolder); end
filename1='MITSWJNTask_SWJN_langElecs_topographical_coordinates.csv';
filename2 = 'MITLangloc_langElecs_topographical_coordinates.csv';    
isPlotVisible = false;
          
Datasets = [2,1];
hemis = {'LH','RH'};
lobes = {'temp','front','all'};


RELIAB_PATH = 'output/_reliability/reliability/';
addpath(genpath('utils/'))

for d=1:2
    Dataset = Datasets(d);
    for ih=1:2
        hemi = hemis{ih};
        if strcmp(hemi,'RH') && (Dataset==1)
            continue
        end
        for il=1:3
            lobe = lobes{il};
            
            switch Dataset
                case 1
                    filename = filename1;
                    expt_suffix = 'MITSWJNTask';
                    cond_string = 'SWJN_';
                    elec_type_cluster = 'langElecs';
                case 2
                    filename = filename2;
                    expt_suffix = 'MITLangloc';
                    cond_string = '';
                    elec_type_cluster = 'langElecs';
            end
            
            %% Load data - a single dataset
            T = readtable([DataFolder filesep 'coordinates/' filename]);
            T.k3 = categorical(T.k3);
            T.r = sqrt(T.x.^2+T.y.^2+T.z.^2);
            colors = [1 0 0;...
                      0 1 0;...
                      0 0 1];
            
            load([RELIAB_PATH expt_suffix '_' cond_string elec_type_cluster '_reliability.mat']); % corrs
            corrs(corrs<0) = 0; 
            corrs = corrs/max(corrs);

            if ~isequal(height(T),size(corrs,1))
                error('number of channels is not equal in table and reliability structures')
            end
            switch hemi
                case 'LH'
                    corrs(T.x>0)=[];
                    T(T.x>0,:)=[];
                case 'RH'
                    corrs(T.x<0)=[];
                    T(T.x<0,:)=[];
            end
            
            switch lobe
                case 'temp'
                    corrs=corrs(strcmp(T.lobe,'temporal'));
                    T=T(strcmp(T.lobe,'temporal'),:);
                case 'front'
                    corrs=corrs(strcmp(T.lobe,'frontal'));
                    T=T(strcmp(T.lobe,'frontal'),:);
                case 'all'
            end
            
            T.catlobe = ones(size(T.lobe));
            T.catlobe(strcmp(T.lobe,'temporal'))=-1;

            

            %% LME
            
            formula_x = 'x ~ k3 + (1 + k3|subject)';
            formula_y = 'y ~ k3 + (1 + k3|subject)';
            formula_z = 'z ~ k3 + (1 + k3|subject)';
            %formula_r = 'r ~ k3 + (1 + k3|subject)';
            formula_lobe = 'catlobe ~ k3 + (1 + k3|subject)';
            
            if strcmp(lobe,'all')
                lme_lobe = fitlme(T,formula_lobe,'FitMethod','REML');

                [beta,betanames,stats] = fixedEffects(lme_lobe,'DFMethod','satterthwaite');
                Tstats=stats2table(stats);
                Tstats.Name = {'whichlobe C1 - reference';'whichlobe C2 - relative to C1';'whichlobe C3 - relative to C1'};

                T2 = T;
                T2.k3 = reordercats(T2.k3,[2,3,1]);
                lme2_lobe = fitlme(T2,formula_lobe,'FitMethod','REML');
                [beta2,betanames2,stats2] = fixedEffects(lme2_lobe,'DFMethod','satterthwaite');
                Tstats2=stats2table(stats2);
                Tstats2.Name = {'whichlobe C2 - reference';'whichlobe C3 - relative to C2';'whichlobe C1 - relative to C2'};

                T3 = T;
                T3.k3 = reordercats(T3.k3,[3,1,2]);
                lme3_lobe = fitlme(T3,formula_lobe,'FitMethod','REML');
                [beta3,betanames3,stats3] = fixedEffects(lme3_lobe,'DFMethod','satterthwaite');
                Tstats3=stats2table(stats3);
                Tstats3.Name = {'whichlobe C3 - reference';'whichlobe C1 - relative to C3';'whichlobe C2 - relative to C3'};
                
                Tallstats = [Tstats;Tstats2;Tstats3];

                 Tdataset_coordinate = table;
                 Tdataset_coordinate.Dataset = cell(size(Tallstats.Name));
                 Tdataset_coordinate.Dataset(:) = {['Dataset ' num2str(Dataset)]}; 
                 Tdataset_coordinate.Predict = cell(size(Tallstats.Name));
                 Tdataset_coordinate.Predict(:) = {'lobe'};
                 Tallstats=[Tdataset_coordinate, Tallstats];
                 writetable(Tallstats,[ResultsFolder filesep 'predictLobe_Dataset' num2str(Dataset) '_' hemi '_' lobe '.csv']);

            end
            
            lme_y = fitlme(T,formula_y,'FitMethod','REML');
            lme_y_w = fitlme(T,formula_y,'FitMethod','REML','Weights',corrs);

            [~,~,stats] = fixedEffects(lme_y,'DFMethod','satterthwaite');
            [~,~,stats_w] = fixedEffects(lme_y_w,'DFMethod','satterthwaite');

            Tstats=stats2table(stats);
            Tstats.Name = {'y C1 - reference';'y C2 - relative to C1';'y C3 - relative to C1'};
            Tstats_w=stats2table(stats_w);
            Tstats_w.Name = {'y C1 - reference';'y C2 - relative to C1';'y C3 - relative to C1'};

            T2 = T;
            T2.k3 = reordercats(T2.k3,[2,3,1]);
            lme2_y = fitlme(T2,formula_y,'FitMethod','REML');
            lme2_y_w = fitlme(T2,formula_y,'FitMethod','REML','Weights',corrs);
            [~,~,stats2] = fixedEffects(lme2_y,'DFMethod','satterthwaite');
            [~,~,stats2_w] = fixedEffects(lme2_y_w,'DFMethod','satterthwaite');
            Tstats2=stats2table(stats2);
            Tstats2.Name = {'y C2 - reference';'y C3 - relative to C2';'y C1 - relative to C2'};
            Tstats2_w=stats2table(stats2_w);
            Tstats2_w.Name = {'y C2 - reference';'y C3 - relative to C2';'y C1 - relative to C2'};

            T3 = T;
            T3.k3 = reordercats(T3.k3,[3,1,2]);
            lme3_y = fitlme(T3,formula_y,'FitMethod','REML');
            lme3_y_w = fitlme(T3,formula_y,'FitMethod','REML','Weights',corrs);

            [~,~,stats3] = fixedEffects(lme3_y,'DFMethod','satterthwaite');
            [~,~,stats3_w] = fixedEffects(lme3_y_w,'DFMethod','satterthwaite');

            Tstats3=stats2table(stats3);
            Tstats3.Name = {'y C3 - reference';'y C1 - relative to C3';'y C2 - relative to C3'};
            Tstats3_w=stats2table(stats3_w);
            Tstats3_w.Name = {'y C3 - reference';'y C1 - relative to C3';'y C2 - relative to C3'};
           
            Tallstats = [Tstats;Tstats2;Tstats3];
            Tallstats_w = [Tstats_w;Tstats2_w;Tstats3_w];


             Tdataset_coordinate = table;
             Tdataset_coordinate.Dataset = cell(size(Tallstats.Name));
             Tdataset_coordinate.Dataset(:) = {['Dataset ' num2str(Dataset)]}; 
             Tdataset_coordinate.Predict = cell(size(Tallstats.Name));
             Tdataset_coordinate.Predict(:) = {'y'};
             Tallstats_y=[Tdataset_coordinate, Tallstats];

             Tdataset_coordinate_w = table;
             Tdataset_coordinate_w.Dataset = cell(size(Tallstats_w.Name));
             Tdataset_coordinate_w.Dataset(:) = {['Dataset ' num2str(Dataset)]}; 
             Tdataset_coordinate_w.Predict = cell(size(Tallstats_w.Name));
             Tdataset_coordinate_w.Predict(:) = {'y'};
             Tallstats_y_w=[Tdataset_coordinate_w, Tallstats_w];           

            lme_z = fitlme(T,formula_z,'FitMethod','REML');
            lme_z_w = fitlme(T,formula_z,'FitMethod','REML','Weights',corrs);

            [~,~,stats] = fixedEffects(lme_z,'DFMethod','satterthwaite');
            Tstats=stats2table(stats);
            Tstats.Name = {'z C1 - reference';'z C2 - relative to C1';'z C3 - relative to C1'};
            [~,~,stats_w] = fixedEffects(lme_z_w,'DFMethod','satterthwaite');
            Tstats_w=stats2table(stats_w);
            Tstats_w.Name = {'z C1 - reference';'z C2 - relative to C1';'z C3 - relative to C1'};

            T2 = T;
            T2.k3 = reordercats(T2.k3,[2,3,1]);
            
            lme2_z = fitlme(T2,formula_z,'FitMethod','REML');
            [beta2,betanames2,stats2] = fixedEffects(lme2_z,'DFMethod','satterthwaite');
            Tstats2=stats2table(stats2);
            Tstats2.Name = {'z C2 - reference';'z C3 - relative to C2';'z C1 - relative to C2'};
            
            lme2_z_w = fitlme(T2,formula_z,'FitMethod','REML','Weights',corrs);
            [~,~,stats2_w] = fixedEffects(lme2_z_w,'DFMethod','satterthwaite');
            Tstats2_w=stats2table(stats2_w);
            Tstats2_w.Name = {'z C2 - reference';'z C3 - relative to C2';'z C1 - relative to C2'};
            
            T3 = T;
            T3.k3 = reordercats(T3.k3,[3,1,2]);
            lme3_z = fitlme(T3,formula_z,'FitMethod','REML');
            [~,~,stats3] = fixedEffects(lme3_z,'DFMethod','satterthwaite');
            Tstats3=stats2table(stats3);
            Tstats3.Name = {'z C3 - reference';'z C1 - relative to C3';'z C2 - relative to C3'};
            
            
            lme3_z_w = fitlme(T3,formula_z,'FitMethod','REML','Weights',corrs);
            [~,~,stats3_w] = fixedEffects(lme3_z_w,'DFMethod','satterthwaite');
            Tstats3=stats2table(stats3_w);
            Tstats3_w.Name = {'z C3 - reference';'z C1 - relative to C3';'z C2 - relative to C3'};
            
            Tallstats = [Tstats;Tstats2;Tstats3];
            Tallstats_w = [Tstats_w;Tstats2_w;Tstats3_w];

             Tdataset_coordinate = table;
             Tdataset_coordinate.Dataset = cell(size(Tallstats.Name));
             Tdataset_coordinate.Dataset(:) = {['Dataset ' num2str(Dataset)]}; 
             Tdataset_coordinate.Predict = cell(size(Tallstats.Name));
             Tdataset_coordinate.Predict(:) = {'z'};
             
             Tdataset_coordinate_w = table;
             Tdataset_coordinate_w.Dataset = cell(size(Tallstats_w.Name));
             Tdataset_coordinate_w.Dataset(:) = {['Dataset ' num2str(Dataset)]}; 
             Tdataset_coordinate_w.Predict = cell(size(Tallstats_w.Name));
             Tdataset_coordinate_w.Predict(:) = {'z'};
             
             Tallstats_z=[Tdataset_coordinate, Tallstats];
             Tallstats_z_w=[Tdataset_coordinate_w, Tallstats_w];         

            %lme_r = fitlme(T,formula_r);
            Tstats = [Tallstats_y; Tallstats_z];
            Tstats_w = [Tallstats_y_w; Tallstats_z_w];
            
            writetable(Tstats,[ResultsFolder filesep 'coordStatsYZ_Dataset' num2str(Dataset) '_' hemi '_' lobe '.csv'])
            writetable(Tstats_w,[ResultsFolder filesep 'coordStatsYZ_weightedByReliability_Dataset' num2str(Dataset) '_' hemi '_' lobe '.csv'])

            
            %% prep data
            
            measures = {'x','y','z','r'};
            means = nan(2,3);
            errs = nan(2,3);
            
            %figure
            i=0;
            for im=2:3
                i=i+1;
               % subplot(2,2,i)
                coord1=T.(measures{im})(T.k3=="1");err1=std(coord1)/sqrt(length(coord1));
                means(i,1)=mean(coord1);errs(i,1)=err1;
                coord2=T.(measures{im})(T.k3=="2");err2=std(coord2)/sqrt(length(coord2));
                means(i,2)=mean(coord2);errs(i,2)=err2;
                coord3=T.(measures{im})(T.k3=="3");err3=std(coord3)/sqrt(length(coord3));
                means(i,3)=mean(coord3);errs(i,3)=err3;
                % v = violinplot(coord1,1,'ShowMean',true);
                % hold on
                % v = violinplot(coord2,2,'ShowMean',true);
            end
            
            %% bar plot
            figure
            set(gcf,'visible',isPlotVisible);
            [h, hE]=barwitherr(errs,means);
            for ib=1:length(h)
                h(ib).FaceColor=colors(ib,:);
                h(ib).LineWidth=1;
            end
            xticklabels(measures(2:3))
            
            set(gca,'fontsize',20)
            legend({'Cluster 1','Cluster 2','Cluster 3'})
            ylabel('MNI coordinate values')
            title(['Dataset' num2str(Dataset) ' ' hemi ' ' lobe])
            
            %% violin plot
            dxs = [0 1 2 4 5 6 8 9 10];
            figure
            set(gcf,'Position',[10,10,400,500],'visible',isPlotVisible)
            v=cell(3,3);
            j=0;
            for i=2:3 % measures x y z
                for ic=1:3 %clusters 1 2 3
                    j=j+1;
                    data = T.(measures{i})(T.k3==num2str(ic));
                    v{i,ic} = violinplot(data,1,'Bandwidth',(max(data)-min(data))/10,'ShowMean',true,'BoxColor',[0 0 0],'ShowData',true,'ViolinColor',colors(ic,:));
                   %v=violinplot(data,1,'BoxColor',[0 0 0],'ShowMean',true,'ViolinAlpha',0.8,'ViolinColor',[0.9 0.9 0.9]);
                   %v{im,j}.ScatterPlot.SizeData = 100;
                   %v{i,ic}.ViolinColor=colors(ic,:);
                   v{i,ic}=moveViolin(v{i,ic},dxs(j));
                   v{i,ic} = invisibleViolin(v{i,ic});
                   v{i,ic}.MeanPlot.Color='k';
                   v{i,ic}.MeanPlot.LineWidth = 5;
                end
            end
            set(gca,'XTick',[2,6,10],'xticklabel',{'X','Y','Z'},'ytick',[-80 -40 0 40 80],'fontsize',20)
            hold on
            plot(get(gca,'xlim'),[0 0],'-.k')
            ylabel('MNI Coordinate Value')
            legend hide
            title(['Dataset' num2str(Dataset) ' ' hemi ' ' lobe ' Channels'])
            saveas(gcf,[ResultsFolder 'violinCoordYZ_Clusters_Dataset' num2str(Dataset) '_' hemi '_' lobe],'pdf')
            
            %% prep data per participant:
            subjects = unique(T.subject);
            Ts=table;
            for is=1:length(subjects)
                for ik=1:3
                    Tk=T(strcmp(T.subject,subjects{is}) & T.k3==categorical(ik),:);
                    switch Dataset
                        case 1
                            tmp=T(1,[1,4,8:11]);
                        case 2
                            tmp=T(1,[1,4,8:11]);
                    end
                    tmp.subject = subjects(is);
                    tmp.k3 = categorical(ik);
                    tmp.x=mean(Tk(strcmp(Tk.subject,subjects{is}),:).x);
                    tmp.y=mean(Tk(strcmp(Tk.subject,subjects{is}),:).y);
                    tmp.z=mean(Tk(strcmp(Tk.subject,subjects{is}),:).z);
                    tmp.r=mean(Tk(strcmp(Tk.subject,subjects{is}),:).r);
                    Ts=[Ts;tmp];
                end
            end
            %% violin plot - by participant
            figure
            set(gcf,'Position',[10,10,400,400],'visible',isPlotVisible)
            switch Dataset
                case 1
                    clusters=[1,2,3];
                    dxs = [0 1 2 4 5 6];
                    xtickloc = [2,6,10];
                case 2
                    clusters=[1,3];
                    dxs = [0 1 3 4];
                    xtickloc = [1.5 4.5];
            end
            v=cell(3,3);
            j=0;
            for i=2:3 % measures x y z
                for ic=clusters %clusters 1 2 3
                    j=j+1;
                    data = Ts.(measures{i})(Ts.k3==num2str(ic));
                    v{i,ic} = violinplot(data,1,'Bandwidth',(max(data)-min(data))/10,'ShowMean',true,'BoxColor',[0 0 0],'ShowData',true,'ViolinColor',colors(ic,:));
                    %v=violinplot(data,1,'BoxColor',[0 0 0],'ShowMean',true,'ViolinAlpha',0.8,'ViolinColor',[0.9 0.9 0.9]);
                    %v{im,j}.ScatterPlot.SizeData = 100;
                    %v{i,ic}.ViolinColor=colors(ic,:);
                    v{i,ic} = moveViolin(v{i,ic},dxs(j));
                    v{i,ic} = invisibleViolin(v{i,ic});
                    v{i,ic}.MeanPlot.Color='k';
                    v{i,ic}.MeanPlot.LineWidth = 2;
            
            %         if ic>1
            %             xi=v{i,ic}.ScatterPlot.XData;
            %             yi=v{i,ic}.ScatterPlot.YData;
            %             xip=v{i,ic-1}.ScatterPlot.XData;
            %             yip=v{i,ic-1}.ScatterPlot.YData;
            %             for is=1:length(xi)
            %                 line([xi(is) xip(is)],[yi(is) yip(is)],'Color','k')
            %             end
            %         end
                end
            end
            set(gca,'XTick',xtickloc,'xticklabel',{'Y','Z'},'ytick',[-80 -40 0 40 80],'fontsize',20)
            hold on
            plot(get(gca,'xlim'),[0 0],'-.k')
            ylabel('MNI Coordinate Value')
            legend hide
            title(['Dataset' num2str(Dataset) ' ' hemi ' ' lobe ' Participants'])
            saveas(gcf,[ResultsFolder 'violinCoordYZ_ClustersPerSubject_Dataset' num2str(Dataset) '_' hemi '_' lobe],'pdf')
            
        end 
    end
end