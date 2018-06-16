%PlotERP_MWP.m - 2018 - Kyle Mathewson and Sayeed Devraj-Kizuk
%set of code to analyze the Moral word ERP project (https://osf.io/5jmze/)
%Load in the data and create ERPs comparing the different conditions

addpath('subfunctions')
addpath('lib');

ccc
eeglab
%uses EEGLAB13_3_2b

%% load in files specified in settings
load('4conds_rej3_Settings.mat')
anal.segments = 'on'; %load the EEG segments
Load_data_MWP(exp,anal)

%% Plot ERPs of your epochs
electrodes = {EEG.chanlocs.labels};
ERP_effect_name = {'P2','N2','P3','LPP'};
ERP_effect_timerange = {[200 250];[250 350];[350 600];[600 800]};
data_out = [];

% The code below uses a nested loo'p to determine which segmented dataset corresponds to the right argument in data_out
% e.g. if you have 5 sets, 20 participants, and 4 events, for i_set ==
% 2 and i_part == 1 and i_event == 1, the code uses the data from set (2-1)*4*20 + (1-1)*20 + 1 == set 81
for i_set = 1:nsets
    for eegset = 1:nevents
        for i_part = 1:nparts
            data_out(:,:,eegset,i_part,i_set) = mean(ALLEEG((i_set-1)*nevents*nparts + (i_part-1)*(nevents) + eegset).data,3);
        end
    end
end


%% Count how many trials were retained after artifact rejection
for i_set = 1:nsets
    for eegset = 1:nevents
        for i_part = 1:nparts
            trial_count(eegset,i_part,i_set) = length(ALLEEG((i_set-1)*nevents*nparts + (i_part-1)*(nevents) + eegset).epoch);
        end
    end
end
totalTrials = squeeze(sum(sum(trial_count,1),3));
justCorrect = squeeze(trial_count(:,:,2));
justCorrectMean = mean(justCorrect,2);
justCorrectSD = std(justCorrect,[],2);

%prints stats for paper in command window
fprintf(['An average of ' num2str(mean(totalTrials)) ' (SD = ' num2str(std(totalTrials)) ') of the original 500 trials remained for analysis after artifact rejection. \n']);
for i=1:4
    fprintf(['An average of ' num2str(justCorrectMean(i)) ' (SD = ' num2str(justCorrectSD(i)) ') of the original correct ' exp.event_names{1,i} 'trials remained for analysis after artifact rejection. \n']);
end

%% Loop through each time window and plot topographies
channels = [7 7 7 8];

for i_set = 2
    for i_effect = 1:length(ERP_effect_name)
        
        i_chan = channels(i_effect);
        
        %% Figure 3A % 4A
        % this is the averaged ERP data. We take the mean across the selected channels (dim1) and the participants (dim4).
        
        if i_effect > 2 %since first three effects are both at Cz only need one plot
            
            erpdata = squeeze(mean(mean(data_out(i_chan,:,:,:,i_set),1),4));
            erpsems = squeeze(SEMws(data_out(i_chan,:,:,:,i_set),3,4));
            
            figure;
            boundedline(EEG.times(801:end),erpdata(801:end,1),erpsems(801:end,1),'r',EEG.times(801:end),erpdata(801:end,2),erpsems(801:end,2),'g',EEG.times(801:end),erpdata(801:end,3),erpsems(801:end,3),'b',EEG.times(801:end),erpdata(801:end,4),erpsems(801:end,4),'c');
            xlim([-200 1000])
            ylim([-2 11]);     set(gca,'ydir','reverse');
            title([ERP_effect_name{i_effect} ' at electrode ' electrodes{i_chan}]);
            
            legend(exp.event_names{i_set,:});
            line([0 0],[-11 11],'color','k');
            line([ERP_effect_timerange{1}(1) ERP_effect_timerange{1}(1)],[-10 11],'color','k');
            line([ERP_effect_timerange{1}(2) ERP_effect_timerange{1}(2)],[-10 11],'color','k');
            line([ERP_effect_timerange{2}(2) ERP_effect_timerange{2}(2)],[-10 11],'color','k');
            line([ERP_effect_timerange{3}(2) ERP_effect_timerange{3}(2)],[-10 11],'color','k');
            line([ERP_effect_timerange{4}(2) ERP_effect_timerange{4}(2)],[-10 11],'color','k');
            line([exp.epochslims*1000],[0 0],'color','k');
            savefig(['AllFour_' ERP_effect_name{i_effect} ' at electrode ' electrodes{i_chan}]);
            
        end
        %% MoralWord NonmoralWord MoralNonword NonmoralNonword
        WordMinusNonword = squeeze(data_out(:,:,[1 2],:,i_set)-data_out(:,:,[3 4],:,i_set));
        MoralMinusNonmoral = squeeze(data_out(:,:,[1 3],:,i_set)-data_out(:,:,[2 4],:,i_set)); %exclude non-words
        
        %% Figure 3B
        % Words vs nonwords
        conds = {'Words','Non-words'};
        time_window = find(EEG.times>ERP_effect_timerange{i_effect}(1),1): find(EEG.times>ERP_effect_timerange{i_effect}(2),1);
        figure('Color',[1 1 1]);
        set(gca,'Color',[1 1 1]);
        erpeffect = squeeze(mean(mean(mean(WordMinusNonword(:,time_window,:,:),4),3),2))';
        erpeffect(1,16:18) = NaN;
        topoplot(erpeffect,'EOG-electrode-locs.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4]  )
        title([ERP_effect_name{i_effect} ' - Words Minus Nonwords']);
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
        savefig([ERP_effect_name{i_effect} ' - Words Minus Nonwords']);
        
        
        %% Figure 4B
        % Moral vs. Non-moral
        conds = {'Moral words','Non-moral words'};
        time_window = find(EEG.times>ERP_effect_timerange{i_effect}(1),1):find(EEG.times>ERP_effect_timerange{i_effect}(2),1);
        figure('Color',[1 1 1]);
        set(gca,'Color',[1 1 1]);
        erpeffect = squeeze(mean(mean(MoralMinusNonmoral(:,time_window,1,:),4),2))';
        erpeffect(1,16:18) = NaN;
        topoplot(erpeffect,'EOG-electrode-locs.ced', 'whitebk','on','plotrad',.6,'maplimits',[-1 1]  )
        title([ERP_effect_name{i_effect} ' - Moral words - Non-moral words']);
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
        
        savefig([ERP_effect_name{i_effect} ' - Moral words - Non-moral words']);
        
        
        %% Figure 3C & 4C
        if i_effect > 2
            % this is the averaged ERP data. We take the mean across the selected channels (dim1) and the participants (dim4).
            erpdata = squeeze(mean(mean(WordMinusNonword(i_chan,:,:,:),3),4));
            erpsems = squeeze(std(mean(WordMinusNonword(i_chan,:,:,:),3),[],4)./sqrt(nparts));
            erpdata2 = squeeze(mean(mean(MoralMinusNonmoral(i_chan,:,1,:),3),4));
            erpstds2 = squeeze(std(mean(MoralMinusNonmoral(i_chan,:,1,:),3),[],4)./sqrt(nparts));
            
            figure; boundedline(EEG.times(801:end),erpdata(801:end),erpsems(801:end),'b',EEG.times(801:end),erpdata2(801:end),erpstds2(801:end),'r');
            xlim([-200 1000])
            ylim([-2 5]);     set(gca,'ydir','reverse');
            title([ERP_effect_name{i_effect} ' at electrode ' electrodes{i_chan}])
            
            legend('Word - Nonword','Moral words - Non-moral words','Location','SouthEast'); legend BOXOFF
            line([0 0],[-10 10],'color','k');
            line([exp.epochslims*1000],[0 0],'color','k');
            line([ERP_effect_timerange{1}(1) ERP_effect_timerange{1}(1)],[-10 11],'color','k');
            line([ERP_effect_timerange{1}(2) ERP_effect_timerange{1}(2)],[-10 11],'color','k');
            line([ERP_effect_timerange{2}(2) ERP_effect_timerange{2}(2)],[-10 11],'color','k');
            line([ERP_effect_timerange{3}(2) ERP_effect_timerange{3}(2)],[-10 11],'color','k');
            line([ERP_effect_timerange{4}(2) ERP_effect_timerange{4}(2)],[-10 11],'color','k');
            
            savefig(['Difference_' ERP_effect_name{i_effect} ' at electrode ' electrodes{i_chan}]);
        end
    end
end



