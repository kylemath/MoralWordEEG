%Analysis_Wrapper_MWP.m - 2018 - Kyle Mathewson and Sayeed Devraj-Kizuk
%set of code to analyze the Moral word ERP project (https://osf.io/5jmze/)


addpath('subfunctions')
ccc %clear and close everything


%% Settings for loading the raw data
%uses EEGLAB13_3_2b
exp.pathname = 'M:\Data\MoralP\EEG';
exp.electrode_locs = 'M:\Data\MoralP\EEG\electrode_locations.ced';

%Datafiles must be in the format exp_participant, e.g. EEGexp_001.vhdr
exp.name = 'MWP';
exp.participants = {'001','002','003','005','006','007','008','009','010','012','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030','031','032','033','035','036','037','038','039','040','041','042','043','044','045','046','047','048','049','050','052','053','054'}; %all subjects
% exp.participants = {'001','002'}; %for testing

% removed data:
% Subject 4  - Data bug - labels don't match words
% Subject 11 - didn't answer
% Subject 34 - too many of one answer
% Subject 51 - too many of one answer


% Choose how to organize your datasets. 
% The sets are different trial types that contain comparable events. 
% Examples of different sets: 10Hz/12Hz stimulation on different trials; emotional vs non-emotional images on different trials.
% Each row is a different set. 
% The events are stimuli within the trial. 
% Examples of different events. Missing vs present targets, valid vs invalid targets, detected vs undetected responses.
% Each column is a different event. 
% You can collect multiple triggers into one event with square brackets [].

exp.events = {[210],[220],[230],[240];...
              [211],[221],[231],[241]};    %must be matrix (sets x events)
          
exp.setname = {'Incorrect';'Correct'}; %name the rows
exp.event_names = {'MoralWord','NonmoralWord','MoralNonword','NonmoralNonword'}; %name the columns

% Each item in the exp.events matrix will become a seperate dataset, including only those epochs referenced by the events in that item. 
%e.g. 3 rows x 4 columns == 12 datasets/participant

% The settings will be saved as a new folder. It lets you save multiple datasets with different preprocessing parameters.
exp.settings = '4conds_rej3';
 


%% Preprocessing Settings
%segments settings
exp.segments = 'on'; %Do you want to make new epoched datasets? Set to "off" if you are only changing the tf settings.
%Where are your electrodes? (.ced file)

%% Filter the data?
exp.filter = 'on';
exp.lowpass = 30;
exp.highpass = 0;

%% Re-referencing the data
exp.refelec = 16; %which electrode do you want to re-reference to?
exp.brainelecs = [1:15]; %list of every electrode collecting brain data (exclude mastoid reference, EOGs, HR, EMG, etc.

%% Epoching the data
%Choose what to epoch to. The default [] uses every event listed above.
%Alternatively, you can epoch to something else in the format {'TRIG'}. Use triggers which are at a consistent time point in every trial.
exp.epochs = {210,220,230,240}; 
exp.epochslims = [-1 1]; %in seconds; epoched trigger is 0 e.g. [-1 2]
exp.epochbaseline = [-200 0]; %remove the for each epoched set, in ms. e.g. [-200 0] 

%% Artifact rejection. 
% Choose the threshold to reject trials. More lenient threshold followed by an (optional) stricter threshold 
exp.preocularthresh = [-1000 1000]; %First happens before the ocular correction.
exp.postocularthresh = [-500 500]; %Second happens after. Leave blank [] to skip
 
%% Blink Correction
%the Blink Correction wants dissimilar events (different erps) seperated by commas and similar events (similar erps) seperated with spaces. See 'help gratton_emcp'
exp.selection_cards = {'210 211','220 221','230 231','240 241'};
%%%%


%% Save your pipeline settings
save([exp.settings '_Settings'],'exp') %save these settings as a .mat file. This will help you remember what settings were used in each dataset




%% Run Preprocessing
Preprocessing_MWP(exp) %comment out if you're only doing analysis



































