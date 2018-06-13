%% MoralWordPopout.m 
% Dr. Kyle Mathewson - University of Alberta
% Moral Word Popout Task
% uses Psychtoolbox Version 3.0.14 
% uses parallelport toolbox found here: https://github.com/kylemath/MathewsonMatlabTools/tree/master/ParallelPorts
% any mention of quest throughout refers to a practice session that did not
% change the target (not a staircasing), a constant target luminance was
% used for each subject

KbName('UnifyKeyNames');
clear all
close all
warning off MATLAB:DeprecatedLogicalAPI
Screen('Preference', 'SkipSyncTests', 1); %for test on an LCD
Priority(2);
rng('shuffle'); %seed the random number generator
 
% Input participant's number, name and date of testing
% Output file will be named after the inputs
Info.number = input('Participant Number:','s');
if isempty(Info.number)
    Info.number = 0; %for testing
end
Info.date = datestr(now,30); % 'dd-mmm-yyyy HH:MM:SS' 


% Create text file to write data in
Filename = ['EEG_data\Subject_' Info.number '_data_' Info.date '.mat'];


%% Open the main window and get dimensions
white=[255,255,255];  %WhiteIndex(window);
black=[0,0,0];   %BlackIndex(window);

grey= (white+black)/2; %background colour

%load the window up
%fMRI projector at NYU is 1024x768, at 60hz. 
%EEG monitor at UofA is 1200x800?? at 120hz. 

screenNumber = max(Screen('Screens')); % Get the maximum screen number i.e. get an external screen if avaliable
[window,rect]=Screen(screenNumber ,'OpenWindow',grey(1));
HideCursor;     
v_res = rect(4);
h_res = rect(3);
v_center = v_res/2;
h_center = h_res/2;
fixation = [h_center-10 v_center-10];

% Get presentation timing information
refresh = Screen('GetFlipInterval', window); % Get flip refresh rate
slack = refresh/2; % Divide by 2 to get slack
Info.refresh = round(1/refresh); 

 % Select specific text font, style and size:
text_size = 20;
% Screen('TextFont',window, 'Courier New');
Screen('TextSize',window, text_size);
%  Screen('TextStyle', window, 1+2); 

 %% Set up parallel port to zero out everything
 
%output triggers in Recorder
        
    %1,2,3,4 by type of word (1-Moral word, 2-Nonmoral word; 3-Moral non-word; 4-Nonmoral non-word)
    
    % Moral word trial: 1-Fixation; 11-Word onset; 21-Mask onset; 31-Response; 41-Correct Response
    %Quest 1-word, 2-nonword
    % Quest word trial: 101-Fixation; 111-Word onset; 121-Mask onset; 131-Response; 
 
 
%initialize the inpoutx64 low-level I/O driver
config_io;
%optional step: verify that the inpoutx64 driver was successfully installed
global cogent;
if( cogent.io.status ~= 0 )
   error('inp/outp installation failed');
end
%write a value to the default LPT1 printer output port (at 0x378)
address_eeg = hex2dec('B010');
outp(address_eeg,0);  %set pins to zero  

%% but use the monitor top left pixel to trigger
trigger_size = [0 0 1 1]; %use [0 0 1 1] for eeg, 50 50 for photodiode

%% Import data from spreadsheet
load('Wordlists.mat')
Words = [MoralWords; NonMoralWords; MoralNonWords; NonMoralNonWords];
Masks = cell(size(Words));
for i_Masks = 1:length(Words)
    mask_word = [];
    for i_mask = 1:length(Words{i_Masks})
        mask_word = [mask_word '&'];
    end
    Masks{i_Masks} = mask_word;
end
Quest_Words = [QuestNonMoralWords; QuestNonMoralNonWords];
word_order = randperm(500);
quest_word_order = randperm(100);
word_types = {'MoralWord'; 'NonMoralWord'; 'MoralNonWord'; 'NonMoralNonWord'};
Types = zeros(size(Words));
for i_Types = 1:length(Words)
    Types(i_Types) = floor((word_order(i_Types)-1)/125)+1; %type of word (1-Moral word, 2-Nonmoral word; 3-Moral non-word; 4-Nonmoral non-word)
end
Answers = ceil(Types/2); %correct answer (1 - word, 2 - nonword);

%% Variables to adjust
%trial numbers    
n_blocks = 10; %10how many blocks
trials_per_block = 50; %50round(n_trials/n_blocks);
n_trials = n_blocks*trials_per_block; %500how many trials overall
Info.fmri_face_order = randi(2); %order of the face manipulation is randomized

%quest settings
n_quest_trials = 45; %45
n_quests = 3; % 3how many quests to run
quest_trials_per_block = round(n_quest_trials/n_quests);

%random fixation lags before target (400 -700 ms)
lags =  48+randi(84-48,1,500); %one for every trial
n_lags = length(lags); %number of unique lags

%response keys
key_word = KbName('1!'); % one '1!' or left arrow 'right'
key_nonword = KbName('5%'); % 5 '5%' or right arrow 'left'
key_mri = KbName('`'); 

%response deadline
resp_timelim = 1.5; %2seconds

%target
vis_thresh = 108; %target RGB points darker than grey background
targ_length = 2; %in refresh cycles; each refresh is 1000msec / 60 Hz = 16.66 msec
targ_grey = grey-vis_thresh; %stim_grey - vis_thresh; %target colour 128-90=38
fixation_offset = 35; %vertical offset above fixation of target and mask

%SOA and mask
mask_thresh = 50; %points darker than background
SOA = 6; %refresh target onset to mask onset (50 ms optimal/16.66 msec = 3 cycles 
ISI = SOA - targ_length; %target OFFSET to mask onset 
mask_length = 3;%refresh cycles of mask
mask_grey = grey-mask_thresh; %mask colourIt no 
wrongRight={'Incorrect','Correct'};



%% Set up the offScreen windows and stimuli

left_arrow = [h_center-10 v_center-30; h_center+10 v_center-40; h_center+10 v_center-20];   %create the attention cues
right_arrow = [h_center+10 v_center-30; h_center-10 v_center-40; h_center-10 v_center-20];

%setup the blank screen
blank=Screen(window,'OpenoffScreenWindow',grey);
    Screen(blank, 'FillRect',black,trigger_size); %eeg trigger for start of trial

%setup the fixation screen
fixate=Screen(window,'OpenoffScreenWindow',grey);
    Screen('TextSize',fixate, text_size);
    DrawFormattedText(fixate,'+','center','center',0);  %Print the fixation,
    Screen(fixate, 'FillRect',black,trigger_size); %eeg trigger

%setup the target and mask backgrounds
target=Screen(window,'OpenoffScreenWindow',grey);
mask=Screen(window,'OpenoffScreenWindow',grey);

    

%% Instructions
DrawFormattedText(window,'+','center','center',0);  %Print the fixation,
DrawFormattedText(window,'Please keep your eyes fixed on the central cross the entire time.', 'center',fixation(2)-100,0); 
Screen('FillRect',window,black,trigger_size); %eeg trigger
Screen('Flip', window,[],0); %flip it to the screen
KbWait; %wait for subject to press button

DrawFormattedText(window,'example','center',fixation(2)-fixation_offset,0);  %Print the example word
DrawFormattedText(window,'+','center','center',0);  %Print the fixation,
DrawFormattedText(window,'On each trial decide if you see a word or a non-word at fixation before the &&&&&&''s.','center',fixation(2)-110,0);  
DrawFormattedText(window,'Press 1 if you see a word, and press 5 if you see a non-word','center',fixation(2)-85,0); 
Screen('FillRect',window,black,trigger_size); %eeg trigger
Screen('Flip', window,[],0);
WaitSecs(1);
KbWait; 

DrawFormattedText(window,'+','center','center',0);  %Print the fixation,
DrawFormattedText(window,'Let the experimenter know if you have any questions. Press any key to begin the practice session.','center',fixation(2)-100,0);
Screen('FillRect',window,black,trigger_size); %eeg trigger
Screen('Flip', window,[],0); %
WaitSecs(1);
KbWait;

%wait a bit before it starts
Screen('CopyWindow',blank ,window,rect,rect);
Screen(window, 'Flip');
WaitSecs(2);

% 
%% Quest Staircase
for i_quest = 1:n_quests  
    
    
    %% Loop for trials
    for i_trial = ((i_quest-1)*quest_trials_per_block)+1:(i_quest*quest_trials_per_block)
        
        %current word and mask and lag
        current_lag = lags(randi(n_lags)); %For quest just pick randomly from the list with replacement
        current_word = Quest_Words{quest_word_order(i_trial)};
        current_type = floor(quest_word_order(i_trial)/50)+1; %type of word (1-Nonmoral word; 2-Nonmoral non-word)
        mask_word = [];
        for i_mask = 1:length(current_word)
            mask_word = [mask_word '&'];
        end       
 
   

        if i_trial == 1 %make practice really easy
            targ_grey = 0;
        end
        


        %draw the target
        target=Screen(window,'OpenoffScreenWindow',grey);
        Screen('TextSize',target, text_size);
        DrawFormattedText(target,current_word,'center',fixation(2)-fixation_offset,targ_grey); %adjust targ_grey in the quest loop
        DrawFormattedText(target,'+','center','center',0);  %Print the fixation,
        Screen(target, 'FillRect',Vpixx2Vamp(110+current_type),trigger_size); %eeg trigger 111 (quest nonmoral word) or 112 (quest nonmoral nonword)
        
        %Setup the mask
        mask=Screen(window,'OpenoffScreenWindow',grey);
        Screen('TextSize',mask, text_size);
        DrawFormattedText(mask,mask_word,'center',fixation(2)-fixation_offset,mask_grey);
        DrawFormattedText(mask,'+','center','center',0);  %Print the fixation,
        Screen(mask, 'FillRect',Vpixx2Vamp(120+current_type),trigger_size); %eeg trigger 121,122
   
        %trial 
        
        %put up fixation
        Screen('CopyWindow',fixate ,window,rect,rect);
        Screen(fixate, 'FillRect',Vpixx2Vamp(100+current_type),trigger_size); %eeg trigger  101,102
        onsets.fix = Screen(window, 'Flip', 0); %already waited on the previous trial
        Screen('CopyWindow',fixate ,window,rect,rect);
        Screen(fixate, 'FillRect',black,trigger_size); %eeg trigger  101,102
        Screen(window, 'Flip', 0); %already waited on the previous trial

        % present the Target
        Screen('CopyWindow',target ,window,rect,rect);
        onsets.target = Screen(window, 'Flip', onsets.fix + current_lag*refresh - slack);

        % blank Inter stimulus interval
        Screen('CopyWindow',fixate ,window,rect,rect);
        Screen(fixate, 'FillRect',black,trigger_size); %eeg trigger  101,102
        onsets.ISI = Screen(window, 'Flip', onsets.target + targ_length*refresh - slack);

        % present the mask
        Screen('CopyWindow',mask ,window,rect,rect);
        onsets.mask = Screen(window, 'Flip', onsets.ISI + ISI*refresh - slack);

        % Response period
        Screen('CopyWindow',blank ,window,rect,rect);
        Screen(blank, 'FillRect',black,trigger_size); %eeg trigger  101,102
        onsets.response_period = Screen(window, 'Flip', onsets.mask + mask_length*refresh - slack);
         
        t1 = GetSecs;
        keyIsDown = 0;
%         while  ~keyIsDown 
        key_pressed = 0;
        while GetSecs - t1 < resp_timelim
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyIsDown
                Screen('CopyWindow',blank ,window,rect,rect);
                Screen(blank, 'FillRect',Vpixx2Vamp(130+current_type),trigger_size); %eeg trigger 40
                onsets.response = Screen(window, 'Flip', 0);
                Screen('CopyWindow',blank ,window,rect,rect);
                Screen(blank, 'FillRect',black,trigger_size);
                Screen(window, 'Flip', 0);
                
                key_pressed = 1;
                response = find(keyCode>0,1);   %1 is 49, 5 is 53 %37 is left arrow, %39 is right arrow
                wrong_key = 0;
                break
            end
        end
        WaitSecs( resp_timelim - (GetSecs - t1)); %wait the rest of time
        Screen('CopyWindow',blank ,window,rect,rect);
        Screen(blank, 'FillRect',black,trigger_size); 
        onsets.response_period = Screen(window, 'Flip', 0);
        
        if key_pressed == 1;
             %keep a log of the subject answers
            if response == key_word
                quest_subject_answer(i_trial) = 1; %word 
            elseif response == key_nonword
                quest_subject_answer(i_trial) = 2; %nonword
            else
                wrong_key = 1;
            end
            if wrong_key == 0
                quest_subject_correct(i_trial) = ~(abs(quest_subject_answer(i_trial)- current_type)); %1 of correct, 0 if incorrect
          
            
            else %if an incorrect button pressed %if you have a second monitor you can watch how they do
                fprintf('Trial %3d at %5.2f is %s\n',i_trial,targ_grey,'wrong key');
            end
        else %if there was no response
            fprintf('Trial %3d at %5.2f is %s\n',i_trial,targ_grey,'unanswered');
        end

    end %i trial  
    
    %% Wait for the subject to move onto the next block, or end the practice.   

   if i_quest == n_quests %last block
        DrawFormattedText(window,'+','center','center',0);  %Print the fixation,
        DrawFormattedText(window,['You have now completed all ' num2str(n_quests) ' practice blocks. Press any key to move on to the experiment.'],'center',fixation(2)-100,0);
        Screen('FillRect',window,black,trigger_size); %eeg trigger
        Screen('Flip', window,[],0); %
        WaitSecs(1);
        KbWait ;
  
         
    end

    

end %iblock



targ_grey = grey-vis_thresh; %note we don't use a quest derived value we use the same value for each subject and quest just acts as a practice



%wait a bit before experiment starts
Screen('CopyWindow',blank ,window,rect,rect);
Screen(window, 'Flip');
WaitSecs(2);
lag_order = [];


%initialize some stuff
t1 = 0;
keyIsDown = 0;
key_pressed = 0;
response = 0;   
response_time = 0;
wrong_key = 0;

subject.coded_answer = zeros(1,n_trials);
subject.answer = zeros(1,n_trials);
subject.rt = zeros(1,n_trials);
subject.correct = zeros(1,n_trials);
subject.type = cell(1,n_trials);
subject.type_num = zeros(1,n_trials);
subject.block = zeros(1,n_trials);
subject.lags = zeros(1,n_trials);
subject.word = cell(1,n_trials);

    onsets.blockstart = zeros(1,n_blocks);
onsets.preload = zeros(1,n_trials);
onsets.fix = zeros(1,n_trials);
onsets.lag = zeros(1,n_trials);
onsets.target = zeros(1,n_trials);
onsets.ISI = zeros(1,n_trials);
onsets.mask = zeros(1,n_trials);
onsets.response_period = zeros(1,n_trials);
onsets.response_button = zeros(1,n_trials);
onsets.posttrial = zeros(1,n_trials);
onsets.posttrialend = zeros(1,n_trials);
    onsets.blockend = zeros(1,n_blocks);

target_triggers = [210 220 230 240];
    
%% Real Experiment
for i_block = 1:n_blocks
   
    %create the random order for this block, and adds them to the list
%     lag_order = [lag_order randperm(n_lags)]; %order used in the task
 
       
    
    %% trigger start of block
    
    if i_block == 1;   %first block
        DrawFormattedText(window,'+','center','center',0);  %Print the fixation,
        if Info.refresh == 120; %MRI study
            message = 'Let the experimenter know if you have any questions. Press any key to begin the first block of the experiment.';
        else
            message = 'Let the experimenter know if you have any questions and are ready to start. Experimenter: Start scanner now.';
        end
        DrawFormattedText(window, message ,'center',fixation(2)-100,0);
        
    else %all other blocks
        
        DrawFormattedText(window,'+','center','center',0);  %Print the fixation,
        if i_block-1 == 1
            suffix = '';
        else
            suffix = 's';
        end
        if Info.refresh == 120; %MRI study
            message = 'Press any key to begin the next block of the experiment.';
        else
            message = 'Experimenter: Start scanner now.';
        end
        DrawFormattedText(window,['You have now completed ' num2str(i_block-1) ' block' suffix ' out of ' num2str(n_blocks) '. ' message],'center',fixation(2)-100,0);
        
    end
    
    Screen('FillRect',window,black,trigger_size); %eeg trigger
    Screen('Flip', window); %
    WaitSecs(1);
    KbWait;
    
    
    
    %wait a bit before the block starts
    Screen('CopyWindow',blank ,window,rect,rect);
    onsets.blockstart(i_block) = Screen(window, 'Flip');
    WaitSecs(2);
    
    %% Loop for trials
    for i_trial = ((i_block-1)*trials_per_block)+1:(i_block*trials_per_block)
        
        
        onsets.preload(i_trial) = GetSecs; %%%%
        
        %current word and mask and lag
        current_lag = lags(i_trial);
        current_word = Words{word_order(i_trial)};
        current_mask = Masks{word_order(i_trial)};
        current_type = Types(i_trial); %type of word (1-Moral word, 2-Nonmoral word; 3-Moral non-word; 4-Nonmoral non-word)
        current_ans = Answers(i_trial); %1 word, 2 nonword     
              
        
                
        %Setup the target word
        target=Screen(window,'OpenoffScreenWindow',grey);
        Screen('TextSize',target, text_size);
        DrawFormattedText(target,current_word,'center',fixation(2)-fixation_offset,targ_grey);
        DrawFormattedText(target,'+','center','center',0);  %Print the fixation,
        Screen(target, 'FillRect',Vpixx2Vamp(target_triggers(current_type)),trigger_size); %eeg trigger 11, 12, 13, 14   (depending on trial type)    
        
        %Setup the mask of &'s
        mask=Screen(window,'OpenoffScreenWindow',grey);
        Screen('TextSize',mask, text_size);
        DrawFormattedText(mask,current_mask,'center',fixation(2)-fixation_offset,mask_grey);        
        DrawFormattedText(mask,'+','center','center',0);  %Print the fixation,
        Screen(mask, 'FillRect',Vpixx2Vamp(20+current_type),trigger_size); %eeg trigger  21, 22, 23, 24 (depending on trial type)
       
        %% trial 
    
        
        %put up fixation
        Screen('CopyWindow',fixate ,window,rect,rect);
        Screen(fixate, 'FillRect',Vpixx2Vamp(current_type),trigger_size); %eeg trigger  1 2 3 4 (depending on trial type)
        onsets.fix(i_trial) = Screen(window, 'Flip', 0); %already waited on the previous trial
        Screen('CopyWindow',fixate ,window,rect,rect);
        Screen(fixate, 'FillRect',black,trigger_size); %eeg trigger  1 2 3 4 (depending on trial type)
        Screen(window, 'Flip', 0); %already waited on the previous trial
            
        % present the Target
        Screen('CopyWindow',fixate ,window,rect,rect);
        Screen(fixate, 'FillRect',black,trigger_size); %eeg trigger  1 2 3 4 (depending on trial type)
        Screen(window, 'Flip', onsets.fix(i_trial) + current_lag*refresh - slack); %already waited on the previous trial        
        Screen('CopyWindow',target ,window,rect,rect);
        onsets.target(i_trial) = Screen(window, 'Flip', 0); %so every lag just got one refresh added to it

        % blank Inter stimulus interval
        Screen('CopyWindow',fixate ,window,rect,rect);
        Screen(fixate, 'FillRect',black,trigger_size); %eeg trigger reset
        onsets.ISI(i_trial) = Screen(window, 'Flip', onsets.target(i_trial) + targ_length*refresh - slack);

        % present the mask
        Screen('CopyWindow',mask ,window,rect,rect);
        onsets.mask(i_trial) = Screen(window, 'Flip', onsets.ISI(i_trial) + ISI*refresh - slack);

        % Response period
        Screen('CopyWindow',blank ,window,rect,rect);
        Screen(blank, 'FillRect',black,trigger_size); %eeg trigger reset
        onsets.response_period(i_trial) = Screen(window, 'Flip', onsets.mask(i_trial) + mask_length*refresh - slack);
        
        t1 = GetSecs;
        keyIsDown = 0;
        key_pressed = 0;
        while GetSecs - t1 < resp_timelim
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyIsDown
                Screen('CopyWindow',blank ,window,rect,rect);
                Screen(blank, 'FillRect',Vpixx2Vamp(30 + current_type),trigger_size); %eeg trigger 31,32,33,34 just response moment for now
                onsets.response_button(i_trial) = Screen(window, 'Flip', 0);
                Screen('CopyWindow',blank ,window,rect,rect);
                Screen(blank, 'FillRect',black,trigger_size); %eeg trigger  1 2 3 4 (depending on trial type)
                Screen(window, 'Flip', 0); %already waited on the previous trial
                
                key_pressed = 1;
                response = find(keyCode>0,1);   %1 is 49, 5 is 53 %37 is left arrow, %39 is right arrow
                response_time = secs;
                wrong_key = 0;
                
                break %stop waiting for key
            end
        end
        WaitSecs( resp_timelim - (GetSecs - t1)); %wait the rest of time
        Screen('CopyWindow',blank ,window,rect,rect);
        Screen(blank, 'FillRect',black,trigger_size); %eeg trigger reset
        onsets.posttrial(i_trial) = Screen(window, 'Flip', 0);

        if key_pressed == 1 %if they pressed something
            %keep a log of the subject answers
            if response == key_word
                subject.coded_answer(i_trial) = 1; %word
            elseif response == key_nonword
                subject.coded_answer(i_trial) = 2; %nonword
            else
                subject.coded_answer(i_trial) = 99; %wrong key
                wrong_key = 1;
            end
            if wrong_key == 0 % if right key
                subject.correct(i_trial) = ~(abs(subject.coded_answer(i_trial)- current_ans)); %1 of correct, 0 if incorrect
                
            else %if an incorrect button pressed
                subject.correct(i_trial) = 3; %incorrect button
            end
            %keep a log of the subject answers 
            subject.answer(i_trial) = response;
            subject.rt(i_trial) = response_time-t1; %compute response time and log        
        else %if there was no response
            subject.coded_answer(i_trial) = 0; %no response
            subject.correct(i_trial) = 3; %no response
            onsets.response_button(i_trial) = 0;
            subject.answer(i_trial) = 0;
            subject.rt(i_trial) = 0; %compute response time and log
        end
              
        subject.type{i_trial} = word_types{current_type}; %(1-Moral word, 2-Nonmoral word; 3-Moral non-word; 4-Nonmoral non-word)
        subject.type_num(i_trial) = current_type;
        subject.block(i_trial) = i_block;
        subject.lags(i_trial) = current_lag;
        subject.word{i_trial} = current_word;
 
        onsets.posttrialend(i_trial) = GetSecs; %%%
        
    end %i trial

    
    
    %% Wait for the subject to move onto the next block, or end the experiment.   

   
    %wait a bit before the next block starts or experiment ends
    Screen('CopyWindow',blank ,window,rect,rect);
    Screen(window, 'Flip');
    WaitSecs(1); %six second ramp down
    
    %timestamp for end of block
    onsets.blockend(i_block) = GetSecs;
    
     if i_block == n_blocks %last block
            DrawFormattedText(window,'+','center','center',0);  %Print the fixation,
            DrawFormattedText(window,['You have now completed all ' num2str(n_blocks) ' blocks. Thank you for your participation.'],'center',fixation(2)-100,0);
            Screen('FillRect',window,black,trigger_size); %eeg trigger
            Screen('Flip', window,[],0); %
            WaitSecs(1);
            KbWait;
     end
     
end %iblock


%% Save the data and close the window, exit
Screen('Close', window);
ShowCursor;
% manip_responses = [resp1 resp2 resp3 resp4];
save(Filename,'subject', 'onsets', 'Info','word_order','quest_word_order','lags','targ_grey');


%% compute some results
for i_cond = 1:4
    results(i_cond) = length(find(subject.correct(subject.type_num == i_cond) == 1))/length(subject.correct(subject.type_num == i_cond));
end

for i_cond = 1:2
    resultscom(i_cond) = (length(find(subject.correct(subject.type_num == i_cond+2) == 1))+ length(find(subject.correct(subject.type_num == i_cond) == 1)))  / (length(subject.correct(subject.type_num == i_cond+2)) + length(subject.correct(subject.type_num == i_cond))) ;
end

figure; subplot(2,1,1); barweb(results, ones(1,4)*0,1,[],[],[],'Proportion Correct','gray',[],{'Moral Words';'Non-moral Words';'Moral Non-words'; 'Non-moral Non-words'},[],'axis');
subplot(2,1,2); barweb(resultscom, ones(1,2)*0,1,[],[],[],'Proportion Correct','gray',[],{'Moral';'Non-moral'},[],'axis');
