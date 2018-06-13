function Preprocessing_MWP(exp)

%Preprocessing_MWP.m - 2018 - Kyle Mathewson and Sayeed Devraj-Kizuk
%set of code to analyze the Moral word ERP project (https://osf.io/5jmze/)
%preprocess the raw EEG data

try
   
    nparts = length(exp.participants);
    nsets = length(exp.setname);
    
    if any(size(exp.event_names) ~= size(exp.events))
        repfactor = size(exp.events)./size(exp.event_names);
        exp.event_names = repmat(exp.event_names, repfactor);
    end
    
    if isempty(exp.epochs) == 1
        exp.epochs = cellstr(num2str(cell2mat( reshape(exp.events,1,size(exp.events,1)*size(exp.events,2)) )'))';
    end
    
    %initialize EEGLAB
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    %subject numbers to analyze
    
    %% Load data and channel locations
    if strcmp('on',exp.segments) == 1
        for i_part = 1:nparts
            sprintf(['Processing Participant ' num2str(exp.participants{i_part})])
            
            %% load a data file
            EEG = pop_loadbv(exp.pathname, [exp.participants{i_part} '_' exp.name '.vhdr']);
            
            %load behaviour
            temp = load(['M:\Data\Moral\BEH\' exp.participants{i_part} '_MWP.mat'])
            if strcmp(exp.participants{i_part},'001') || strcmp(exp.participants{i_part},'005') || strcmp(exp.participants{i_part},'006')
                temp.subject.type_num = floor((temp.word_order-1)/125)+1;  %recompute the correct number
                temp.subject.type_corr_resp = ceil(temp.subject.type_num/2);  %and correct response
                temp.subject.correct = ~abs(temp.subject.coded_answer-temp.subject.type_corr_resp);  %now compute the correct answers
            end
            
            % load channel information
            EEG=pop_chanedit(EEG, 'load',{exp.electrode_locs 'filetype' 'autodetect'});
            
            %% Filter the data with low pass of 30
            if strcmp('on',exp.filter) == 1
                EEG = pop_eegfilt( EEG, exp.highpass, exp.lowpass, [], 0);
            end
            
            %% arithmetically rereference to linked mastoid (M1 + M2)/2
            for x=exp.brainelecs
                EEG.data(x,:) = (EEG.data(x,:)-((EEG.data(exp.refelec,:))*.5));
            end
            
            %change markers so they can be used by the gratton_emcp script
            allevents = length(EEG.event);
            for i_event = 2:allevents %skip the first
                EEG.event(i_event).type = num2str(str2num(EEG.event(i_event).type(2:end)));
            end
            
            %% The triggers are early
            [EEG] = VpixxEarlyTriggerFix(EEG);
            
            %% Extract epochs of data time locked to event
            %Extract data time locked to targets and remove all other events
            EEG = pop_epoch( EEG, exp.epochs, exp.epochslims, 'newname', [exp.participants{i_part} '_epochs'], 'epochinfo', 'yes');
            
            %replace with triggers that are correct
            if strcmp(exp.participants{i_part},'001') || strcmp(exp.participants{i_part},'005') || strcmp(exp.participants{i_part},'006')
                for i_epoch = 1:length(EEG.epoch)
                    EEG.epoch(i_epoch).eventtype([EEG.epoch(i_epoch).eventlatency{:}] == 0) = { num2str(200+(temp.subject.type_num(i_epoch)*10)) };
                    EEG.event( EEG.epoch(i_epoch).event([EEG.epoch(i_epoch).eventlatency{:}] == 0)).type = 200+(temp.subject.type_num(i_epoch)*10);
                end
            end
            
            if strcmp(exp.participants{i_part},'017')  %used the number pad
                temp.subject.answer(temp.subject.answer == 97) = 49;
                temp.subject.coded_answer(temp.subject.answer == 49) = 1;
                temp.subject.answer(temp.subject.answer == 101) = 53;
                temp.subject.coded_answer(temp.subject.answer == 53) = 2;
                temp.subject.type_corr_resp = ceil(temp.subject.type_num/2);  %and correct response
                temp.subject.correct = ~abs(temp.subject.coded_answer-temp.subject.type_corr_resp);  %now compute the correct answers
            end
            
            if strcmp(exp.participants{i_part},'018') %reversed keys
                temp.subject.correct = ~temp.subject.correct;
            end
            
            if strcmp(exp.participants{i_part},'050') %first trials did not record
                temp.subject.correct(1:6) = [];
                temp.subject.type_num(1:6) = [];
            end
            
            %replace triggers with triggers indicating accuracy %1 of correct, 0 if incorrect
            for i_epoch = 1:length(EEG.epoch)
                EEG.epoch(i_epoch).eventtype([EEG.epoch(i_epoch).eventlatency{:}] == 0) = { num2str(200+(temp.subject.type_num(i_epoch)*10)+temp.subject.correct(i_epoch)) };   %replace with 241 if correct, 240 incorrect
                EEG.event( EEG.epoch(i_epoch).event([EEG.epoch(i_epoch).eventlatency{:}] == 0)).type = 200+(temp.subject.type_num(i_epoch)*10)+temp.subject.correct(i_epoch);
            end
            
            %subtract baseline
            EEG = pop_rmbase( EEG, exp.epochbaseline);
            %% Artifact Rejection, Correction, then 2nd Rejection
            
            %Artifact rejection, trials with range >exp.preocularthresh uV
            if isempty(exp.preocularthresh) == 0
                rejtrial = struct([]);
                EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],exp.preocularthresh(1),exp.preocularthresh(2),EEG.xmin,EEG.xmax,0,1);
                rejtrial(i_part,1).ids = find(EEG.reject.rejthresh==1);
            end
            
            %EMCP occular correction
            temp_ocular = EEG.data(end-1:end,:,:); %to save the EYE data for after
            EEG = gratton_emcp(EEG,exp.selection_cards,{'VEOG'},{'HEOG'}); %this assumes the eye channels are called this
            EEG.emcp.table %this prints out the regression coefficients
            EEG.data(end-1:end,:,:) = temp_ocular; %replace the eye data
            
            %Artifact rejection, trials with range >exp.postocularthresh uV
            if isempty(exp.postocularthresh) == 0
                EEG = pop_rmbase( EEG, exp.epochbaseline); %baseline again since this changed it
                EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],exp.postocularthresh(1),exp.postocularthresh(2),EEG.xmin,EEG.xmax,0,1);
                rejtrial(i_part,2).ids = find(EEG.reject.rejthresh==1);
            end
            
            %% Event coding
            %here we rename the event codes to the event names in the wrapper so stimuli can be easily identified later
            for i_trial = 1:EEG.trials
                for i_set = 1:nsets
                    nevents = length(exp.events(i_set,:));
                    for i_event = 1:nevents
                        nperevent = length(exp.events{i_set,i_event});
                        for j_event = 1:nperevent
                            EEG.epoch(i_trial).eventcode(find(strcmp(num2str(exp.events{i_set,i_event}(j_event)),EEG.epoch(i_trial).eventtype))) = exp.event_names(i_event);
                        end
                    end
                end
            end
            
            %here we collect the triggers for each trial that matched an event
            output = struct([]);
            output(i_part).target(EEG.trials) = 0;
            output(i_part).targetlatency(EEG.trials) = 0;
            for i_trial = 1:EEG.trials
                for i_event = 1:nevents
                    if isempty(find(strcmp(exp.event_names(i_event),EEG.epoch(i_trial).eventcode)== 1,1)) == 0
                        output(i_part).target(i_trial) = str2num( EEG.epoch(i_trial).eventtype{ find(strcmp(exp.event_names(i_event),EEG.epoch(i_trial).eventcode)== 1,1) } );
                        output(i_part).targetlatency(i_trial) = EEG.epoch(i_trial).eventlatency{ find(strcmp(exp.event_names(i_event),EEG.epoch(i_trial).eventcode)==1,1) };
                    end
                end
            end
            
            %here we make the trial type and latency a NaN if no trial matched an event
            output(i_part).target(output(i_part).target == 0) = NaN;
            output(i_part).targetlatency(output(i_part).targetlatency == 0) = NaN;
            
            %here we replace rejected trials with NaNs as well
            if isempty(exp.postocularthresh) == 0
                for i_trial = rejtrial(i_part,2).ids
                    output(i_part).target = [output(i_part).target(1:i_trial-1), NaN, output(i_part).target(i_trial:end)];
                    output(i_part).targetlatency = [output(i_part).targetlatency(1:i_trial-1), NaN, output(i_part).targetlatency(i_trial:end)];
                end
            end
            
            if isempty(exp.preocularthresh) == 0
                for i_trial = rejtrial(i_part,1).ids
                    output(i_part).target = [output(i_part).target(1:i_trial-1), NaN, output(i_part).target(i_trial:end)];
                    output(i_part).targetlatency = [output(i_part).targetlatency(1:i_trial-1), NaN, output(i_part).targetlatency(i_trial:end)];
                end
            end
            %here we save a variable into the EEG file with the original trial numbers, with NaNs for rejected trials as well as trials that matched no events
            EEG.replacedtrials.target = output(i_part).target;
            EEG.replacedtrials.targetlatency = output(i_part).targetlatency;
            
            for i_set = 1:nsets
                sprintf(exp.setname{i_set})
                
                if ~exist([exp.pathname '\' exp.settings '\'])
                    mkdir([exp.pathname '\' exp.settings '\']);
                end
                
                if ~exist([exp.pathname '\' exp.settings '\Segments\' exp.setname{i_set} '\'])
                    mkdir([exp.pathname '\' exp.settings '\Segments\' exp.setname{i_set} '\']);
                end
                %% Select individual events and Save
                EEG.exp = exp;
                setEEG = EEG;   %replace the stored data with this new set
                for i_event = 1:nevents
                    filename = [exp.name '_' exp.settings '_' exp.participants{i_part} '_' exp.setname{i_set} '_' exp.event_names{i_set,i_event}]
                    eventEEG = pop_selectevent(setEEG, 'type', exp.events{i_set,i_event}, 'deleteevents','off','deleteepochs','on','invertepochs','off');
                    eventEEG = pop_editset(eventEEG, 'setname', filename );
                    eventEEG = pop_saveset(eventEEG, 'filename',['Segments_' filename '.set'],'filepath',[exp.pathname  '\' exp.settings '\Segments\' exp.setname{i_set} '\']);
                end
            end
        end
    end
    
    
    
catch ME
    A = who;
    for i = 1:length(A)
        assignin('base', A{i}, eval(A{i}));
    end
    throw(ME)
end
end