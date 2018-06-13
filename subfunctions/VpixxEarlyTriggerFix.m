function [EEG] = VpixxEarlyTriggerFix(EEG)

    %adds 6 ms to every trigger time in order to account for the delay between
    %the trigger and the pixel onset in the top left corner (could always add
    %more time to account for the delay until the stimulus pixel is
    %illuminated)

    for i_event = 1:length(EEG.event); 
        EEG.event(i_event).latency = EEG.event(i_event).latency+6;
    end


end