%SingleTrial_MWP.m - 2018 - Kyle Mathewson and Sayeed Devraj-Kizuk
%set of code to analyze the Moral word ERP project (https://osf.io/5jmze/)
%load in behavioural datasets and put everything in a big spreadsheet
%save csv file with data for analysis

addpath('subfunctions')
ccc
eeglab

%% load in settings and segments
load('4conds_rej3_Settings.mat')
anal.segments = 'on'; %load the EEG segments?
Load_data_MWP(exp,anal)

%% Behavioural data
nsets = 3;
for i_part = 1:nparts
    load(['M:\Data\Moral\BEH\' exp.participants{i_part} '_MWP.mat']);
    
    if strcmp(exp.participants{i_part},'001') || strcmp(exp.participants{i_part},'005') || strcmp(exp.participants{i_part},'006') %wrong button
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
    
    if strcmp(exp.participants{i_part},'050') %eeg didn't start for first couple trials
        subject.correct(1:5) = [];
        subject.type_num(1:5) = [];
        subject.word(1:5) = [];
        subject.answer(1:5) = [];
        subject.rt(1:5) = [];
    end
    
    names = {'MoralWord','NonMoralWord','MoralNonWord','NonMoralNonWord'};
    answertype = [ 0 1 3];
    for i_set = 1:nsets
        for eegset = 1:nevents
            unrejectedtrials = ones(1,500);
            if strcmp(exp.participants{i_part},'050')==1
                 unrejectedtrials = ones(1,495);
            end
            
            wordspresented{eegset,i_part,i_set} = subject.word(unrejectedtrials == 1 & subject.correct == answertype(i_set) & subject.type_num==eegset);
            buttonpress{eegset,i_part,i_set} = subject.answer(unrejectedtrials == 1 & subject.correct == answertype(i_set) & subject.type_num==eegset);
            reactiontime{eegset,i_part,i_set} = subject.rt(unrejectedtrials == 1 & subject.correct == answertype(i_set) & subject.type_num==eegset);
            trialnum{eegset,i_part,i_set} = find(unrejectedtrials == 1 & subject.correct == answertype(i_set) & subject.type_num==eegset);
        end
    end
end

%% Get the columns to paste to excel
Window1Pz = []; Window1Fz = []; Window2Pz = []; Window2Fz = []; Window3Pz = []; Window3Fz = []; SubjectNumber = []; TrialNumber = []; WordCategory = []; LetterString = []; ButtonPress = []; Accuracy = []; WordType = []; WordCategory = []; RT = []; originalorder = 0;
Categories =  {'Moral', 'Nonmoral', 'Moral', 'Nonmoral'};
Types = {'Word','Word','Nonword','Nonword'};
exp.setname = {'Incorrect';'Correct';'No Response'};
for i_part = 1:nparts
    for i_set = 1:nsets
        for eegset = 1:nevents
            numberoftrials = length(trialnum{eegset,i_part,i_set}');
            SubjectNumber = [SubjectNumber; repmat(exp.participants(i_part),[1 numberoftrials])'];
            TrialNumber = [TrialNumber; trialnum{eegset,i_part,i_set}'];
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
%then pasted into all columns of GEE_Data_MWP.xlsx

xlswrite('testout.xlsx',[SubjectNumber num2cell(TrialNumber) WordType WordCategory LetterString num2cell(ButtonPress) Accuracy num2cell(RT)])

