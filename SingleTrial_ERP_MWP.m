%SingleTrial_ERP_MWP.m - 2018 - Kyle Mathewson and Sayeed Devraj-Kizuk
%set of code to analyze the Moral word ERP project (https://osf.io/5jmze/)
%load in datasets and pull out single trial erp averages
%save csv file with data for analysis

addpath('subfunctions')
ccc
eeglab

%% load in settings and segments
load('4conds_rej3_Settings.mat')
anal.segments = 'on'; %load the EEG segments?
Load_data_MWP(exp,anal)

%% Plot ERPs of your epochs
% An ERP is a plot of the raw, usually filtered, data, at one or multiple electrodes. It doesn't use time-frequency data.
% We make ERPs if we have segmented datasets that we want to compare across conditions.

% In this case, you are using all your electrodes, not a partset.
electrodes = {EEG.chanlocs.labels};
% Type "electrodes" into the command line. This will show you which number to use for i_chan

ERP_effect_name = {'P2','N2','P3','LPP'};
ERP_effect_timerange = {[200 250];[250 350];[350 600];[600 800]};
n_effect = length(ERP_effect_name);

clear all_trials

% The code below uses a nested loop to determine which segmented dataset corresponds to the right argument in data_out
% e.g. if you have 5 sets, 20 participants, and 4 events, for i_set ==
% 2 and i_part == 1 and i_event == 1, the code uses the data from set (2-1)*4*20 + (1-1)*20 + 1 == set 81
for i_set = 1:nsets
    for eegset = 1:nevents
        for i_part = 1:nparts
            alltrials{eegset,i_part,i_set}= ALLEEG((i_set-1)*nevents*nparts + (i_part-1)*(nevents) + eegset).data;
            alltrigs = [{ALLEEG((i_set-1)*nevents*nparts + (i_part-1)*(nevents) + eegset).urevent.type}];
            alltrigs(strcmp(alltrigs,'S210')==1 | strcmp(alltrigs,'S220')==1 | strcmp(alltrigs,'S230')==1 | strcmp(alltrigs,'S240')==1) = {'Target'};
            alltrigs(strcmp(alltrigs,'Target')==0) = {'NonTarget'};
            
            for target = 1:500
                alltrigs(find(strcmp('Target',alltrigs),1)) = {num2str(target)};
            end
            
            if (i_set-1)*nevents*nparts + (i_part-1)*(nevents)+ eegset == 100 %this guy has just 1 epoch, need to use event field
                triggernumber = [{ALLEEG((i_set-1)*nevents*nparts + (i_part-1)*(nevents) + eegset).event.urevent}];
                triggertype = [{ALLEEG((i_set-1)*nevents*nparts + (i_part-1)*(nevents) + eegset).event.type}];
            else
                triggernumber = [ALLEEG((i_set-1)*nevents*nparts + (i_part-1)*(nevents) + eegset).epoch.eventurevent];
                triggertype = [ALLEEG((i_set-1)*nevents*nparts + (i_part-1)*(nevents) + eegset).epoch.eventtype];               
            end
 
            triggernumber = cell2mat(triggernumber);
            targetnumber = alltrigs(triggernumber(strcmp(num2str(exp.events{i_set,eegset}),triggertype)));
            S = sprintf('%s ', targetnumber{:});
            
            trialnum{eegset,i_part,i_set} = sscanf(S, '%f');
        end
    end
end

%% loop through ang get the average erp for each effect
for i_set = 1:nsets
    for i_effect = 1:n_effect %ERP components
        time_window = find(EEG.times>ERP_effect_timerange{i_effect}(1),1):find(EEG.times>ERP_effect_timerange{i_effect}(2),1);
        for eegset = 1:nevents %experiment conditions
            for i_part = 1:nparts
                submeans{i_set,eegset,i_part,i_effect} = squeeze(mean(alltrials{eegset,i_part,i_set}(:,time_window,:),2));
            end
        end
    end
end


%% load in beh and fix a few that were wrong

for i_part = 1:nparts
    load(['M:\Data\Moral\BEH\' exp.participants{i_part} '_MWP.mat']);
    
    if strcmp(exp.participants{i_part},'001') || strcmp(exp.participants{i_part},'005') || strcmp(exp.participants{i_part},'006') %wrong buttons
        subject.type_num = floor((word_order-1)/125)+1;  %recompute the correct number
        subject.type_corr_resp = ceil(subject.type_num/2);  %and correct response
        subject.correct = ~abs(subject.coded_answer-subject.type_corr_resp);  %now compute the correct answers
    end
    
    if strcmp(exp.participants{i_part},'017')  %used the number pad
        subject.answer(subject.answer == 97) = 49;
        subject.coded_answer(subject.answer == 49) = 1;
        subject.answer(subject.answer == 101) = 53;
        subject.coded_answer(subject.answer == 53) = 2;
        subject.type_corr_resp = ceil(subject.type_num/2);  %and correct response
        subject.correct = ~abs(subject.coded_answer-subject.type_corr_resp);  %now compute the correct answers
    end
    
    if strcmp(exp.participants{i_part},'018') %reversed keys
        subject.correct = ~subject.correct;
    end
    
    if strcmp(exp.participants{i_part},'050') %first trials no eeg
        subject.correct(1:5) = [];
        subject.type_num(1:5) = [];
        subject.word(1:5) = [];
        subject.answer(1:5) = [];
        subject.rt(1:5) = [];
    end
    
    names = {'MoralWord','NonMoralWord','MoralNonWord','NonMoralNonWord'};
    
    %get some of the behavoiural info
    for i_set = 1:nsets
        for eegset = 1:nevents
            unrejectedtrials = zeros(1,500);
            if strcmp(exp.participants{i_part},'050')==1
                 unrejectedtrials = zeros(1,495);
            end
            alltrialnums = sort(trialnum{eegset,i_part,i_set});
            unrejectedtrials(alltrialnums) = 1;
            
            wordspresented{eegset,i_part,i_set} = subject.word(unrejectedtrials == 1 & subject.correct == i_set-1 & subject.type_num==eegset);
            buttonpress{eegset,i_part,i_set} = subject.answer(unrejectedtrials == 1 & subject.correct == i_set-1 & subject.type_num==eegset);
            reactiontime{eegset,i_part,i_set} = subject.rt(unrejectedtrials == 1 & subject.correct == i_set-1 & subject.type_num==eegset);
        
        end
    end
end

%% Get the columns to paste to excel
Window0Pz = []; Window0Cz = []; Window1Pz = []; Window1Cz = []; Window2Pz = []; Window2Cz = []; Window3Pz = []; Window3Cz = [];  
SubjectNumber = []; TrialNumber = []; WordCategory = []; LetterString = []; ButtonPress = []; Accuracy = []; 
WordType = []; WordCategory = []; RT = []; originalorder = 0;

Categories =  {'Moral', 'Nonmoral', 'Moral', 'Nonmoral'};
Types = {'Word','Word','Nonword','Nonword'};
ChanCz = 8;  %moral vs nonmoral at Cz
ChanPz = 7;  %words - nonwords at Pz

for i_set = 1:nsets
    for i_part = 1:nparts
        for eegset = 1:nevents
            numberoftrials = length(submeans{i_set,eegset,i_part,1}(ChanPz,:)'); %could use any channel
            
            Window0Pz = [Window0Pz ; submeans{i_set,eegset,i_part,1}(ChanPz,:)']; %P2
            Window0Cz = [Window0Cz ; submeans{i_set,eegset,i_part,1}(ChanCz,:)'];
            
            Window1Pz = [Window1Pz ; submeans{i_set,eegset,i_part,2}(ChanPz,:)']; %N2
            Window1Cz = [Window1Cz ; submeans{i_set,eegset,i_part,2}(ChanCz,:)'];
            
            Window2Pz = [Window2Pz ; submeans{i_set,eegset,i_part,3}(ChanPz,:)']; %p3
            Window2Cz = [Window2Cz ; submeans{i_set,eegset,i_part,3}(ChanCz,:)'];
            
            Window3Pz = [Window3Pz ; submeans{i_set,eegset,i_part,4}(ChanPz,:)']; %LPP
            Window3Cz = [Window3Cz ; submeans{i_set,eegset,i_part,4}(ChanCz,:)'];
            
            SubjectNumber = [SubjectNumber; repmat(str2num(exp.participants{i_part}),[1 numberoftrials])'];
            TrialNumber = [TrialNumber; trialnum{eegset,i_part,i_set}];
            WordType = [WordType; repmat(Types(eegset),[1 numberoftrials])'];
            WordCategory = [WordCategory; repmat(Categories(eegset),[1 numberoftrials])'];
            LetterString = [LetterString; wordspresented{eegset,i_part,i_set}'];
            ButtonPress = [ButtonPress; buttonpress{eegset,i_part,i_set}'];
            Accuracy = [Accuracy; repmat(exp.setname(i_set),[1 numberoftrials])'];
            RT = [RT; reactiontime{eegset,i_part,i_set}'];
            originalorder = [originalorder; (originalorder(end)+1:originalorder(end)+numberoftrials)'];
        end
    end
end

%writes out the data to xlsx to save for GMM, where the still 
%NEED TO BE SORTED by subjectnumber (A) and trialnumber(B) in xcell and
%then pasted into GEE_Data_ERP_MWP.xlsx final three columns
xlswrite('testout_ERP.xlsx',[SubjectNumber TrialNumber Window0Pz Window0Cz Window1Pz Window1Cz Window2Pz Window2Cz Window3Pz Window3Cz])

