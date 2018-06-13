function Load_data_MWP(exp,anal)

%Load_data_MWP.m - 2018 - Kyle Mathewson and Sayeed Devraj-Kizuk
%set of code to analyze the Moral word ERP project (https://osf.io/5jmze/)
%load in datasets

sprintf(exp.settings)
try
    %Some lines here to make the variables work in different cases
    
    nparts = length(exp.participants);
    nsets = length(exp.setname);
    
    if any(size(exp.event_names) ~= size(exp.events))
        repfactor = size(exp.events)./size(exp.event_names);
        exp.event_names = repmat(exp.event_names, repfactor);
    end
 
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
    %% Load the data
    %The main loop loops over events, then participants, then sets.
    for i_set = 1:nsets
        exp.setname(i_set)
        tic
        for i_part = 1:nparts
            sprintf(['Loading Participant ' num2str(exp.participants{i_part}) '...' ])
            
            
            nevents = length(exp.events(i_set,:));
            for i_event = 1:nevents
                filename = [exp.name '_' exp.settings '_' exp.participants{i_part} '_' exp.setname{i_set} '_' exp.event_names{i_set,i_event}]
                
                
                % Load the EEGLAB datasets, if needed.
                if strcmp('on',anal.segments) == 1
                    try
                        EEG = pop_loadset('filename',['Segments_' filename '.set'],'filepath',[exp.pathname  '\' exp.settings '\Segments\' exp.setname{i_set} '\']);
                        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                    catch
                        %                     WaitSecs(.5)
                        pause(.5)
                        EEG = pop_loadset('filename',['Segments_' filename '.set'],'filepath',[exp.pathname  '\' exp.settings '\Segments\' exp.setname{i_set} '\']);
                        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                    end
                end
                
            end
        end
        toc
    end
    
    A = who;
    for i = 1:length(A)
        assignin('base', A{i}, eval(A{i}));
    end
    
catch ME
    A = who;
    for i = 1:length(A)
        assignin('base', A{i}, eval(A{i}));
    end
    throw(ME)
end
end