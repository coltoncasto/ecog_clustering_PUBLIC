%% Definitions (currently not a general function)
saveName = 'MITSWJNTask_SWJN_kmedoids_correlation_FIGURES';
[CLUSTER_PATH,SAVE_PATH] = initialize(saveName);
DataFolder = [CLUSTER_PATH 'stimuli' filesep];
NgramFolder = [DataFolder 'ngram' filesep];
ResultsFolder = [SAVE_PATH 'plots' filesep 'stimuli' filesep];
if ~exist(ResultsFolder,'dir')
    mkdir(ResultsFolder);
end

sows = {'sentences','words'};
%% Read data

Tw = readtable([DataFolder 'stimulus_words_capitalized.csv']);
Ts = readtable([DataFolder 'stimulus_sentences_capitalized.csv']);

Tw = readtable([DataFolder 'stimulus_words_lowercase.csv']);
Ts = readtable([DataFolder 'stimulus_sentences_lowercase.csv']);
%% create Ngrams

for n=1:4
    for isow=1:2
        nametag = sows{isow};
    
        switch nametag
            case 'sentences'
                Tt = Ts;
            case 'words'
                Tt = Tw;
        end
    
        T = table;
        j=0;%row count
        for i=1:height(Tt)
            for ii=1:(8-n+1)% starting point of ngram, running along a single line/stimulus
                j=j+1;
                for in = 1:n% within-ngram index 
                    T(j,in) = Tt(i,ii+in-1);
                end
            end
        end
        writetable(T,[NgramFolder 'stimulus_' nametag '_ngram_n' num2str(n) '.csv'])
    end
end

%% %%%%%%% Run python script to scrape google ngram data and save.
% (google_ngram_SWstim.ipynb)

%% load and arrange results
sent_means = nan(8,1);
sent_stes = nan(8,1);
word_means = nan(8,1);
word_stes = nan(8,1);

for n=1:4
    for isow=1:2
        nametag = sows{isow};
        opts = detectImportOptions([NgramFolder 'frequency_' nametag '_ngram_n' num2str(n) '.csv'],'NumHeaderLines', 1);
        opts.VariableNamesLine = 1;
        T = readtable([NgramFolder 'frequency_' nametag '_ngram_n' num2str(n) '.csv'],opts);
        
        w = size(T,2);
        for i=1:height(T)
            if ~isnan(T{i,3})
                T{i,w+1}=mean(T{i,3:w});
            else
                T{i,w+1}=0;
            end

        end
        switch nametag
            case 'sentences'
                sent_means(n) = mean(T{:,w+1});
                sent_stes(n) = std(T{:,w+1})/sqrt(height(T));
            case 'words'
                word_means(n) = mean(T{:,w+1});
                word_stes(n) = std(T{:,w+1})/sqrt(height(T));
        end
        
    end
end
%% histograms:
h=ERPfigure;
set(h,'Position',[10 10 600 700])

j=0;
for n=1:4
    for isow=1:2
        j=j+1;
        nametag = sows{isow};
        opts = detectImportOptions([NgramFolder 'frequency_' nametag '_ngram_n' num2str(n) '.csv'],'NumHeaderLines', 1);
        opts.VariableNamesLine = 1;
        T = readtable([NgramFolder 'frequency_' nametag '_ngram_n' num2str(n) '.csv'],opts);
        
        w = size(T,2);
        for i=1:height(T)
            if ~isnan(T{i,3})
                T{i,w+1}=mean(T{i,3:w});
            else
                T{i,w+1}=0;
            end
        end
        T.Properties.VariableNames{'Var13'} = 'MeanFreq';
        
        if isow==2
            xlim = get(gca,'xlim');
            subplot(4,2,j)
            %hold on
            histogram(T.MeanFreq)
            %set(gca,'xlim',xlim)
        else
            subplot(4,2,j)
            histogram(T.MeanFreq)
        end  
        if n==1
            title(nametag,'fontsize',16)
        end
        set(gca,'fontsize',16)
        ylabel(['N = ' num2str(n)])

        hold on
        ylim = get(gca,'ylim');
        line([mean(T.MeanFreq) mean(T.MeanFreq)],[ylim(1) ylim(2)],'Color','red','linewidth',2)
        disp([num2str(n) ' ' nametag])
        disp(mean(T.MeanFreq))
        
        switch nametag
            case 'sentences'
                disp(sent_means(n))
                m = sent_means(n);
                ste = sent_stes(n);
            case 'words'
                disp(word_means(n))
                m = word_means(n);
                ste = word_stes(n);
        end 
        text(m + 10*ste,200,[num2str(m,2) ' +- ' num2str(ste,2)],'fontsize',16)
        set(gca,'ylim',ylim)
    end
end
saveas(gcf,[ResultsFolder 'Ngram_freq_distributions'],'png')
%% plot
figure
ratio=sent_means(1:4)./word_means(1:4);
ste_ratio = ratio.*sqrt((sent_stes(1:4)./sent_means(1:4)).^2 + (word_stes(1:4)./word_means(1:4)).^2);
plot(ratio,'.-','linewidth',2,'MarkerSize',40)
%errorbar(ratio,ste_ratio,'.-','linewidth',2,'MarkerSize',40)
hold on
set(gca,'fontsize',16,'xtick',[1 2 3 4])
xlabel('TRW / Ngram length (words)')
ylabel('ratio of Ngram frequencies P(sent)/P(word) ')
title('SENT to WORD frequency ratio increases with TRW/Ngram length')

figure
plot(sent_means(1:4)./word_means(1:4),'.-','linewidth',2,'MarkerSize',40)
set(gca,'fontsize',16,'xtick',[1 2 3 4])
xlabel('TRW / Ngram size')
set(gca,'yScale','log')
xlabel('TRW / Ngram length (words)')
ylabel('ratio of Ngram frequencies P(sent)/P(word), log scale')
title('SENT to WORD frequency ratio increases with TRW/Ngram length')

figure
%plot(sent_means(1:4),'.-','linewidth',2,'MarkerSize',40)
errorbar(sent_means(1:4),sent_stes(1:4),'.-','linewidth',2,'MarkerSize',40)
hold on
%plot(word_means(1:4),'r.-','linewidth',2,'MarkerSize',40)
errorbar(word_means(1:4),word_stes(1:4),'r.-','linewidth',2,'MarkerSize',40)
legend('senrences','wordlists')
set(gca,'fontsize',14,'xtick',[1 2 3 4])
set(gca,'yScale','log')
xlabel('TRW / Ngram length (words)')
ylabel('google Ngram frequencies, log scale')
title({'SENT and WORD Ngram frequencies','become less similar with TRW/Ngram length'})
legend({'sentences','wordlists'},'fontsize',16)
saveas(gcf,[ResultsFolder 'SENT_vs_WORD_Ngram_freq_log'],'png')
saveas(gcf,[ResultsFolder 'SENT_vs_WORD_Ngram_freq_log'],'pdf')

figure
%plot(sent_means(1:4),'.-','linewidth',2,'MarkerSize',40)
errorbar(sent_means(1:4),sent_stes(1:4),'.-','linewidth',2,'MarkerSize',40)
hold on
%plot(word_means(1:4),'r.-','linewidth',2,'MarkerSize',40)
errorbar(word_means(1:4),word_stes(1:4),'r.-','linewidth',2,'MarkerSize',40)
legend('senrences','wordlists')
set(gca,'fontsize',14,'xtick',[1 2 3 4])
xlabel('TRW / Ngram length (words)')
ylabel('google Ngram frequencies')
title({'SENT and WORD Ngram frequencies','become less similar with TRW/Ngram length'})
legend({'senrences','wordlists'},'fontsize',16)
saveas(gcf,[ResultsFolder 'SENT_vs_WORD_Ngram_freq'],'png')
