%
% gratton_emcp() - Run the eye movement correction procedure developed
%                  by Gratton, Coles, and Donchin (1983) on an epoched eeglab
%                  dataset.
% Usage:
%   >> OUT_EEG = gratton_emcp( EEG, selection_cards, eog_channel_pairs )
%
% Required Inputs:
%
%   EEG                 - input epoched EEG dataset
%   selection_cards     - [cell array of strings] specifying averaging bins (see below)
%   eog_channel_pairs   - [cell array of strings] one or more EOG channels or pairs of channels,
%                                                 the first of which records blinks (see below)
%
% Output%
%
%   OUT_EEG             - corrected EEG dataset, including new structure EEG.emcp with diagnostic information
%
%  Examples:
%
%  1. Performs traditional Gratton emcp with one vertical bipolar and one horizontal bipolar EOG channel:
%
%    selection_cards =  {'1 2','3 4', '5 6','7 8'};
%    EEG = gratton_emcp(EEG,selection_cards,{'VEOG'},{'HEOG'});
%
%  2. Performs traditional Gratton emcp with vertical and horizontal channels derived from monopolar recordings:
%
%    selection_cards =  {'1 2','3 4', '5 6','7 8'};
%    EEG = gratton_emcp(EEG,selection_cards,{'LVEOGUP','LVEOGLO'},{'HEOGL','HEOGR'});
%
%  3. Regression has three EOG channels: left vertical, horizontal, right vertical:
%
%    selection_cards =  {'1 2','3 4', '5 6','7 8'};
%    EEG = gratton_emcp(EEG,selection_cards,{'LVEOGUP','LVEOGLO'},{'HEOGL','HEOGR'},{'RVEOGUP','RVEOGLO'});
%
%  4. First channel in the regression is vertical EOG consisting of the average of left and right vertical channels:
%
%    selection_cards =  {'1 2','3 4', '5 6','7 8'};
%    EEG = gratton_emcp(EEG,selection_cards,{'LVEOGUP RVEOGUP','LVEOGLO REOGLO'},{'HEOGL','HEOGR'});
%
%  5. Correct for a single EOG channel (must be a bipolar vertical channel, corrects only blinks and vertical EOG)
%
%    selection_cards =  {'1 2','3 4', '5 6','7 8'};
%    EEG = gratton_emcp(EEG,selection_cards,{'VEOG'});
%
%
%  A structure EEG.emcp containing diagnostic information is added to output EEG dataset. Example:
%
%          bintotals: [0 0 80 80]       - number matching each selection card
%              table: [91x40 char]      - easy-to-read table of diagnostic information
%          markblink: [512x321 logical] - which data points (pts x epochs) are blinks
%      blink_factors: [72x2 double]     - propagation factors for blinks
%    saccade_factors: [72x2 double]     - propagation factors for saccades
%     regress_blinks: 11773             - number of blink points in regression
%   regress_saccades: 70147             - number of saccade points in regression
% regress_propblinks: 0.1437            - blink proportion in regression
%     overall_blinks: 23943             - number of blink points in all data
%   overall_saccades: 140409            - number of saccade points in all data
% overall_propblinks: 0.1457            - blink proportion in all data
%
% Author:
%
%   Bill Gehring, University of Michigan
%   Version 0.1, July 2007 - Initial version
%   Version 0.2, July 2007 - Added capability for multiple, average, and bipolar EOG channels
%   Version 0.3, July 2007 - Add checks for zero blinks and zero saccades,
%                          - Performance tuning, including bsxfun for ver >= R2007a
%
%   Based on the algorithm reported in:
%
%	Gratton, G., Coles, M. G. H., & Donchin, E. (1983).  A new method for
%	off-line removal of ocular artifact.  Electroencephalography and
%	Clinical Neurophysiology, 55, 468-484.  (Reports the original algorithm
%	for one channel of vertical EOG).
%
%	Miller, G. A., Gratton, G., & Yee, C. M. (1988).  Generalized
%	implementation of an eye movement correction procedure.
%	Psychophysiology, 25, 241-243. (Reports a distribution of this program
%	for the RT11 operating system FORTRAN IV that corrects for vertical EOG
%	and horizontal EOG.)
%
%   Derived from the original C code developed at the University of
%   Illinois Cognitive Psychophysiology Laboratory, available at
%   http://www.umich.edu/~wgehring/emcp2001.zip.
%
%  Notes:
%
%  What is a selection card?  Selection cards (yes, they used to be paper
%  cards) identify epochs according to time-locking events (i.e.,
%  EEG.epoch.eventtypes that have latency of 0). Epochs identified by a
%  selection card are used to create an average waveform, and those
%  averages are subtracted from each single trial before correction factors
%  are computed.  The aim of this is so that the correction factors are not
%  affected by event-related activity (real ERPs).  In practice, you want
%  to choose averages that will differ meaningfully, so if you have an
%  oddball task where trial type 1 is frequent and trial type 2 is
%  infrequent, you should have separate selection cards for trial types 1
%  and 2.  Before the regression is performed, all EEG epochs with trial
%  type 1 will have the frequent average subtracted.
%
%  In the example shown above, four averages are computed. The syntax is
%  that trial types between single quotes are OR'd together, so the first
%  selection card will match any epoch with the time-locking event 1 or 2.
%
%  Epochs that do not match a selection card are not used in calculating
%  the regression, but those data are corrected (using the correction
%  factors used on the other data).  This is useful if you have separate
%  stimulus-locked and response-locked epochs from the same trial. You can
%  specify selection cards only for the stimulus-locked epochs, which
%  ensures that each data point is only used once in the regression.
%
%  Assumptions of the procedure:
%
%  1) The EEG dataset must have at least one channel that can be used to derive
%     vertical EOG.
%
%  2) The epochs identified by the selection cards should be
%     non-overlapping, otherwise data points will be used more than once in
%     computing the regression.
%
%  3) EEG data are already scaled in microVolts.
%
%  4) The datafile is epoched.
%
%  Differences from original Gratton procedure:
%
%  In the case where each epoch matches one selection card and there is
%  one vertical and ine horizontal EOG channel, the output from this
%  function should match the emcpx C code output from the Cognitive
%  Psychophysiology Laboratory at Illinois. (see
%  http://www.umich.edu/~wgehring/emcp2001.zip).
%
%  "Trash data" are data epochs that do not match one of the selection
%  cards. All epochs in the data file are corrected, even those that don't
%  match a selection card. If there are such trash data in the input EEG
%  data, the output of this function will be slightly different from the
%  Illinois C program as the C code used trash data epochs to compute the
%  variance but left the trash data out of other computations.  Greg
%  Miller's FORTRAN code (which more closely corresponds to Gratton's
%  original FORTRAN code) used the trash data in all steps of computing the
%  regression. Thus, it is unclear whether the use of the trash data in the
%  C code for one step of the procedure but not for others was a desired
%  feature.  In any case, the trash data play no role in computing the
%  correction factors in the present function. Modern epoched datasets can
%  be time-locked to any event (e.g., stimulus, response), so computing a
%  single average across such hetereogeneous events and subtracting it from
%  single trials does not make sense. Leaving the trash data out also
%  allows the desirable feature in which certain data are excluded from
%  computing the correction factors.
%
%  The function does not do the following things that the original code
%  did:
%
%  1. The original Gratton code checked for data out of scale and flat
%     lines and discarded those trials.
%
%  2. In some cases it tried to reconstruct the blink on the VEOG channel
%     when data on the channel went out of scale.
%
%     Because of these changes, artifact rejection and other data
%     screening should take place prior to emcp, because the procedure is
%     sensitive to bad data.
%
function EEG = gratton_emcp(EEG, selection_cards, varargin)
% % % %
tic;
starttime = cputime;
if nargin < 3
    errstring = cat(2,'gratton_emcp(): syntax error\n',...
        'syntax: EEG = gratton_emcp(EEG, selection_cards, eog_channel_pairs)\n',...
        'example: EEG = gratton_emcp(EEG,{''1 2'',''3 4''},{''LVEOGUP'',''LVEOGLO''},{''HEOGR'',''HEOGL''})');
    error(errstring,'')
end
if EEG.trials == 1
    error('gratton_emcp(): EEG dataset must be epoched');
end
disp('gratton_emcp(): Correcting data for eye movement artifacts')
for i = 1:length(varargin)
    chanpair{i} = varargin{i};
end

% % % % for testing as script:
% % % % [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% % % % EEG = pop_loadset('pfk0805.bdf_read_12bin.set');
% % % % selection_cards =  {'1 2','3 4','5 6','7 8'};
% % % % chanpair{1} = {'LVEOGUP','LVEOGLO'};
% % % % chanpair{2} = {'HEOGL','HEOGR'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Determine which channels to use for EOG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bsxyes = exist('bsxfun','builtin');
bsxyes = 0
[electrodes{1:size(EEG.chanlocs,2)}] = deal(EEG.chanlocs.labels);
electrodes = electrodes';

EEGdata = zeros(EEG.nbchan + length(chanpair),EEG.pnts,EEG.trials);
table_label=cell(1,length(chanpair));
for i = 1:length(chanpair)
    if length(chanpair{i}) > 2
        error('too many channel labels')
        return
    end
    for chans=1:length(chanpair{i})
        [LABEL1, LABEL2] = strtok(chanpair{i}{chans});  % parse the string
        if ~isempty(strtrim(LABEL2))  %if there was > one label in the string, take the average
            [s,r] = strtok(strtrim(LABEL2));
            if ~isempty(r)
                error('too many channel labels')
                return
            end
            chan1 = find(strcmpi(electrodes,LABEL1));
            chan2 = find(strcmpi(electrodes,strtrim(LABEL2)));
            % create subtracted channels
            EOGCHAN{i}{chans} = (EEG.data(chan1,:,:)+ EEG.data(chan2,:,:))/2;
            EOGLABEL{i}{chans} = ['mean(' LABEL1 ',' strtrim(LABEL2) ')'];
        else
            if(~isempty(strtok(LABEL2)))
                error('too many channel labels')
                return
            else
                chan1 = find(strcmpi(electrodes,LABEL1));
                EOGCHAN{i}{chans} = EEG.data(chan1,:,:);
            end
            EOGLABEL{i}{chans} = LABEL1;
        end
    end
    if length(chanpair{i}) == 2   %subtract if there are two
        EOGCHAN{i}{1} = EOGCHAN{i}{1} - EOGCHAN{i}{2};
        table_label{i} = [strtrim(EOGLABEL{i}{1}) '-' strtrim(EOGLABEL{i}{2})];
    else  %channel data already derived
        table_label{i} = strtrim(EOGLABEL{i}{1});
    end
    EEGdata(i,:,:) = EOGCHAN{i}{1};
end
EEGdata(length(chanpair)+1:end,:,:) = EEG.data;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Identify blinks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute blink template information based on epoched data
parameters.trialLength = EEG.pnts;
blinkCriteria = 14.0*1.0;  %  slope
msTen = round((10.0/(1000/EEG.srate))+0.5);  % number of sample points in 10 ms
thirdOfBlink = msTen*7;  %third of blink is 70 ms
if ~mod(thirdOfBlink,2)
    thirdOfBlink = thirdOfBlink + 1;  %force odd
end
middleOfBlink = round((thirdOfBlink+1)/2);  %point # for template center
lengthOfWind = thirdOfBlink*3.0;  % length of blink window 210 ms
initBlinkScan = thirdOfBlink + middleOfBlink;  %first point for scan
endBlinkScan = EEG.pnts - initBlinkScan - thirdOfBlink + 1;  %last point to scan
windVariance = 2.0;
sumall = zeros(3,100);

covariance = zeros(EEG.pnts,EEG.trials);
markBlink = false(EEG.pnts,EEG.trials);
if bsxyes
    for i = 1:EEG.trials  %fastest method)
        vec2scan = [initBlinkScan:endBlinkScan]'; %32:460
        templateidx=[1:lengthOfWind]-(thirdOfBlink+middleOfBlink); %-31:31
        % yes this template isn't exactly symmetrical, but that's the way
        % it was in the original
        blinkcov=[-1 * ones(1,thirdOfBlink+1) 2 * ones(1,thirdOfBlink) -1 * ones(1,thirdOfBlink-1)]; % -1 -1 2 2 -1 -1    = 1 x 63
        idxvec = bsxfun(@plus,templateidx,vec2scan);  %vectors to compare to template 1 - scanpts
        multvec = bsxfun(@times,reshape(EEGdata(1,idxvec,i)',[size(idxvec)]),blinkcov);  %mutiply data x covariance vector
        covariance = sum(multvec,2)/lengthOfWind;
        slope=covariance/(windVariance*windVariance);  %429 pts 1:429
        blinkmiddleidx = vec2scan(abs(slope)>blinkCriteria); %have to use vec2scan indices, not slope indices
        middlevec =-middleOfBlink+1:middleOfBlink-1;
        blinkidx = unique(bsxfun(@plus,middlevec,blinkmiddleidx));
        markBlink(blinkidx,i) = true;
    end
    %         elseif old == 3 %completely vectorized is slower than 2, true for lots of trials?
    %             k = false(EEG.pnts,EEG.trials);
    %             k(initBlinkScan:endBlinkScan,1:EEG.trials) = true;
    %             vec2scan = find(k); % get the linear indices  (32:420) column; (1 x 137709)
    %             EOGvec=EEGdata(1,:,:);
    %             %%%% maybe use sub2ind somehow?
    %             templateidx=[1:lengthOfWind]-(thirdOfBlink+middleOfBlink); %-31:31  = 1 x 63
    %             % yes this template isn't exactly symmetrical, but that's the way
    %             % it was in the original
    %             blinkcov=[-1 * ones(1,thirdOfBlink+1) 2 * ones(1,thirdOfBlink) -1 * ones(1,thirdOfBlink-1)]; % -1 -1 2 2 -1 -1    = 1 x 63
    %             idxvec = bsxfun(@plus,templateidx,vec2scan);  % (137709 x 63) = vectors to compare to template 1 - scanpts  = 429 x 63
    %             multvec = bsxfun(@times,EOGvec(idxvec),blinkcov);  % 429 x 63
    %             covariance = sum(multvec,2)/lengthOfWind;  % 429
    %             slope=covariance/(windVariance*windVariance);  %429 pts 1:429
    %             blinkmiddleidx = vec2scan(abs(slope)>blinkCriteria); %have to use scanvec indices, not slope indices  ~79
    %             middlevec =-middleOfBlink+1:middleOfBlink-1;  %1 x 21
    %             blinkidx = unique(bsxfun(@plus,middlevec,blinkmiddleidx));  %111
    %             markBlink(blinkidx) = true;  %512 x 31
    %             % end
    %             % %
    %
else
    for pt = initBlinkScan:endBlinkScan
        % performance note: computing these on the fly as indices doesn't save
        % time, initBlinkscan = 32, thirdofblink = 21
        start1 = pt - initBlinkScan+1;  %1
        end1 = start1 + thirdOfBlink;   %22
        start2 = end1 + 1; %23
        end2 = start1 + 2 * thirdOfBlink; %43
        start3 = end2+1;  %44
        end3 = pt + initBlinkScan-1;  %63
        covariance(pt,:) = covariance(pt,:) - squeeze(sum(EEGdata(1,start1:end1,:),2))';
        covariance(pt,:) = covariance(pt,:) + 2 * squeeze(sum(EEGdata(1,start2:end2,:),2))';
        covariance(pt,:) = covariance(pt,:) - squeeze(sum(EEGdata(1,start3:end3,:),2))';
    end

    covariance = covariance/lengthOfWind;
    slope = covariance/(windVariance*windVariance);
    for pt = initBlinkScan:endBlinkScan
        for tr = 1:EEG.trials
            if(abs(slope(pt,tr))>blinkCriteria)
                markBlink(pt-middleOfBlink+1:pt+middleOfBlink-1,tr) = true;
            end
        end
    end
end
%%%%%%%%%%%%%%%
%
% test for zero blink condition
%%%markBlink = markBlink * 0.0;
if bsxyes
    EEGdata = bsxfun(@minus,EEGdata,mean(EEGdata,2));
else
    for tr = 1:EEG.trials

        for ch=1:EEG.nbchan+length(chanpair)
            % subtract mean of each waveform from the wavefrom
            EEGdata(ch,:,tr) = EEGdata(ch,:,tr) - mean(EEGdata(ch,:,tr));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Match epochs to selection cards
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binned_indices = [];
for cards=1:length(selection_cards)
    indices =  eeg_getepochevent(EEG, strread(selection_cards{cards},'%s'));
%     indices =  eeg_getepochevent(EEG, {'X  1'});
    subindices = find(indices == 0);
    % used later on to determine which bin a trial falls under:
    EEG.emcp.bintotals(cards) = length(subindices);
    if ~isempty(subindices) % ~isempty(subindices)
        binned_indices = cat(2,binned_indices,subindices);
        %
        %  subtract the appropriate mean from each trial
        %
        if bsxyes %512 datapoints ~.30 sec
            EEGavg = mean(EEGdata(:,:,subindices),3);
            EEGdata(:,:,subindices) = bsxfun(@minus,EEGdata(:,:,subindices),EEGavg);
        else  %512 datapoints ~.26 sec
            EEGavg = mean(EEGdata(:,:,subindices),3); %faster outside of the loop!
            for idx = 1:length(subindices)  % do this for EOG?
                EEGdata(:,:,subindices(idx)) = EEGdata(:,:,subindices(idx)) - EEGavg;
            end
        end
    end
end
if isempty(binned_indices)
    warning('gratton_emcp(): No trial matches a selection card')
    EEG.emcp.table = [];
    return
end

binned_indices = sort(binned_indices);
rmarkBlink = markBlink(:,binned_indices);
rblinkpts = numel(find(rmarkBlink));
rsaccpts = numel(find(~rmarkBlink));

disp(['gratton_emcp(): Regression blink points: ' num2str(rblinkpts)])
disp(['gratton_emcp(): Regression saccade points: ' num2str(rsaccpts)])
disp(['gratton_emcp(): Regression blink proportion: ' num2str(rblinkpts/(rblinkpts+rsaccpts))]);

blinkpts = numel(find(markBlink));
saccpts = numel(find(~markBlink));
disp(['gratton_emcp(): Corrected blink points: ' num2str(blinkpts)])
disp(['gratton_emcp(): Corrected saccade points: ' num2str(saccpts)])
disp(['gratton_emcp(): Corrected blink proportion: ' num2str(blinkpts/(blinkpts+saccpts))]);
%
%  Identify unusual conditions with zero blinks or saccades
%
%  Because the beginning and end of the waveform will not be tagged as blink,
%  rsaccpts and saccpts should always both be > 0
%
warn=[];
if blinkpts == 0     % subject didn't blink
    warn{1} = sprintf('*** Warning:');
    warn{2} = sprintf('*** No blinks found in data.');
    warn{3} = sprintf('*** All data treated as saccades.');
    for i = 1:length(warn)
        disp(['gratton_emcp(): ' char(warn{i})]);
    end
elseif blinkpts > 0 && rblinkpts == 0    % subject only blinked during non-selected trials
    warn{1} = sprintf('*** Warning:');
    warn{2} = sprintf('*** No blinks found in data selected for regression.');
    warn{3} = sprintf('*** All data treated as saccades.');
    for i = 1:length(warn)
        disp(['gratton_emcp(): ' char(warn{i})]);
    end
    markBlink(:)=false;
    blinkpts = numel(find(markBlink));
    saccpts = numel(find(~markBlink));    
end
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Compute regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take the subset of the data used for computing regression
%
subdata = EEGdata(:,:,binned_indices);
%
% Subtract mean of blinks or saccades used for regression
%
if ~bsxyes

    for ch=1:EEG.nbchan+length(chanpair)
        if rblinkpts > 0
            subdata(ch,rmarkBlink) = subdata(ch,rmarkBlink)- mean(subdata(ch,rmarkBlink));
        end
        subdata(ch,~rmarkBlink) = subdata(ch,~rmarkBlink)- mean(subdata(ch,~rmarkBlink));
    end

    PropBlink = zeros(EEG.nbchan+length(chanpair),length(chanpair));
    PropSaccade = zeros(EEG.nbchan+length(chanpair),length(chanpair));

    xvecB = subdata(1:length(chanpair),rmarkBlink)';
    xvecS = subdata(1:length(chanpair),~rmarkBlink)';
    if rblinkpts > 0
        for ch = 1:EEG.nbchan+length(chanpair)
            yvec = subdata(ch,rmarkBlink)';
            %  Perform regression - no intercept term
            %  Propagation factors are unstandardized coefficients
            XTXI = inv(xvecB' * xvecB);
            PropBlink(ch,:) = XTXI * xvecB' * yvec ;
        end
    end
    for ch = 1:EEG.nbchan+length(chanpair)
        yvec = subdata(ch,~rmarkBlink)';
        XTXI = inv(xvecS' * xvecS);
        %        whos XTXI xvecS yvec;
        %        whos PropSaccade
        PropSaccade(ch,:) = XTXI * xvecS' * yvec ;
    end
else
    PropBlink = zeros(EEG.nbchan+length(chanpair),length(chanpair));
    PropSaccade = zeros(EEG.nbchan+length(chanpair),length(chanpair));
    if rblinkpts > 0
        subdata(:,rmarkBlink) = bsxfun(@minus,subdata(:,rmarkBlink),mean(subdata(:,rmarkBlink),2));
        xvecB = subdata(1:length(chanpair),rmarkBlink)';
        XTXI = inv(xvecB' * xvecB);
        for ch=1:EEG.nbchan+length(chanpair)
            yvec = subdata(ch,rmarkBlink)';
            %  Perform regression - no intercept term
            %  Propagation factors are unstandardized coefficients
            %  XTXI = 2 x 2
            %  xvecB' = 2 x 11773
            %  yvec = 11773 x 1s
            PropBlink(ch,:) = XTXI * xvecB' * yvec ;
        end
    end
    subdata(:,~rmarkBlink) = bsxfun(@minus,subdata(:,~rmarkBlink),mean(subdata(:,~rmarkBlink),2));
    xvecS = subdata(1:length(chanpair),~rmarkBlink)';
    XTXI = inv(xvecS' * xvecS);
    for ch=1:EEG.nbchan+length(chanpair)
        % yvec = subdata(ch,~rmarkBlink)';
        %        whos XTXI xvecS yvec;
        %        whos PropSaccade
        PropSaccade(ch,:) = XTXI * xvecS' * subdata(ch,~rmarkBlink)';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Create output table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
key=cell(1,length(chanpair));
for i = 1:length(chanpair)
    if i == 1
        disp(['gratton_emcp(): Vertical EOG data = ' table_label{1}]);
        key{1} =  ['Vertical = ' table_label{1}];
    else
        disp(['gratton_emcp(): EOG Channel #' num2str(i) ' data = ' table_label{i}]);
        key{i} =  ['EOG' num2str(i) ' = ' table_label{i}];
    end
end


emcplabels = [];
for i = 1:length(chanpair)
    if i == 1
        emcplabels = {'Vertical'};
    else
        emcplabels{i} = ['EOG' num2str(i)];
    end
end
emcplabels = strjust(char([ emcplabels, {EEG.chanlocs.labels}]'),'right');
labelwidth = size(emcplabels,2);

chanheader = [];
for i = 1:length(chanpair)
    if i == 1
        chanheader = strjust(sprintf('%8s','Vert'),'center');
    else
        chanheader = [ chanheader strjust(sprintf('%8s',['EOG' num2str(i)]),'center')];
    end;
end

chanheader = [chanheader chanheader]; %blinks + saccades
datawidth = size(chanheader,2);
tablewidth = max(labelwidth + size(chanheader,2),size(char(key),2));
if tablewidth<120
    tablewidth = 120;
end
filler = tablewidth-size(emcplabels,2)-datawidth;

blk = strjust(sprintf(['%' num2str(8*length(chanpair)) 's'],'Blinks'),'center');
sac = strjust(sprintf(['%' num2str(8*length(chanpair)) 's'],'Saccades'),'center');

EEG.emcp.table(1,:) = strjust(sprintf(['%' num2str(tablewidth) 's'],'Eye Movement Correction Procedure'),'left');
EEG.emcp.table = [EEG.emcp.table; sprintf(blanks(tablewidth))];
EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],['Regression blink points: ' num2str(rblinkpts)]),'left')];
EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],['Regression saccade points: ' num2str(rsaccpts)]),'left')];
EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],['Regression blink proportion: ' num2str(rblinkpts/(rblinkpts+rsaccpts))]),'left')];
EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],['Corrected blink points: ' num2str(blinkpts)]),'left')];
EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],['Corrected saccade points: ' num2str(saccpts)]),'left')];
EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],['Corrected blink proportion: ' num2str(blinkpts/(blinkpts+saccpts))]),'left')];
EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],['Epochs per bin: ' num2str(EEG.emcp.bintotals)]),'left')];
EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],['Total epochs: ' num2str(EEG.trials)]),'left')];
if ~isempty(warn)
    EEG.emcp.table = [EEG.emcp.table; sprintf(blanks(tablewidth))];
    for i = 1:length(warn)
        EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],char(warn{i})),'left')];
    end
end
EEG.emcp.table = [EEG.emcp.table; sprintf(blanks(tablewidth))];
EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],'Propagation Factors'),'left')];
EEG.emcp.table = [EEG.emcp.table; sprintf(blanks(tablewidth))];
for i=1:length(chanpair)
    EEG.emcp.table = [EEG.emcp.table; strjust(sprintf(['%' num2str(tablewidth) 's'],key{i}),'left')];
end
EEG.emcp.table = [EEG.emcp.table; sprintf(blanks(tablewidth))];
EEG.emcp.table = [EEG.emcp.table; [blanks(labelwidth) strjust(sprintf(['%' num2str(datawidth) 's'],[blk sac]),'center') blanks(filler)]];
EEG.emcp.table = [EEG.emcp.table; [blanks(labelwidth) strjust(sprintf(['%' num2str(datawidth) 's'],chanheader),'center') blanks(filler)]];

tabline =  char(zeros(EEG.nbchan+length(chanpair),datawidth+filler));
for i = 1:EEG.nbchan+length(chanpair)
    tmp = '';
    for j = 1:length(chanpair)
        tmp = [tmp sprintf('%8.1f',PropBlink(i,j))];
    end
    for j = 1:length(chanpair)
        tmp = [tmp sprintf('%8.1f',PropSaccade(i,j))];
    end
    tabline(i,:) = [tmp blanks(filler)];
end
EEG.emcp.table = [EEG.emcp.table; [emcplabels tabline] ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Correct the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstitute EEGdata so it no longer has bin mean subtracted
% also, this will enable correction of non-selected data
for i = 1:length(chanpair)
    for chans=1:length(chanpair{i})
        [LABEL1, LABEL2] = strtok(chanpair{i}{chans});  % parse the string
        if ~isempty(strtrim(LABEL2))  %if there was > one label in the string, take the average
            %%%%%%   [s,r] = strtok(strtrim(LABEL2));
            chan1 = find(strcmpi(electrodes,LABEL1));
            chan2 = find(strcmpi(electrodes,strtrim(LABEL2)));
            % create average channels
            EOGCHAN{i}{chans} = (EEG.data(chan1,:,:)+ EEG.data(chan2,:,:))/2;
        else
            chan1 = find(strcmpi(electrodes,LABEL1));
            EOGCHAN{i}{chans} = EEG.data(chan1,:,:);
        end
    end
    if length(chanpair{i}) == 2   %subtract if there are two
        EOGCHAN{i}{1} = EOGCHAN{i}{1} - EOGCHAN{i}{2};
    end
    EEGdata(i,:,:) = EOGCHAN{i}{1};
end
EEGdata(length(chanpair)+1:end,:,:) = EEG.data;
%
% Subtract epoch mean from ALL trials incl non-binned
%
if bsxyes
    EEGdata = bsxfun(@minus,EEGdata,mean(EEGdata,2));
else
    for tr = 1:EEG.trials
        for ch=1:EEG.nbchan+length(chanpair)
            % subtract mean of each waveform from the wavefrom
            EEGdata(ch,:,tr) = EEGdata(ch,:,tr) - mean(EEGdata(ch,:,tr));
        end
    end
end
%
% Compute mean blink and saccade values used in computing the regression
%

if bsxyes
    subdata=EEGdata(:,:,binned_indices);
    EEGdata(length(chanpair)+1:end,:,:) = EEG.data;
    if rblinkpts > 0
        EEGblinkmean = mean(subdata(:,rmarkBlink),2);  %zero if no blinks
        deltaBlink = bsxfun(@minus,EEGdata(:,markBlink),EEGblinkmean);% 
        avgAdjustmentB = zeros(size(deltaBlink));

        for i = 1:length(chanpair)
            avgAdjustmentB = avgAdjustmentB + bsxfun(@times,deltaBlink,PropBlink(:,i)); %
        end
        EEGdata(:,markBlink) = EEGdata(:,markBlink) - (bsxfun(@plus,avgAdjustmentB,EEGblinkmean));% 
    end
    EEGsaccademean = mean(subdata(:,~rmarkBlink),2);
    deltaSaccade = bsxfun(@minus,EEGdata(:,~markBlink),EEGsaccademean); %  
    avgAdjustmentS = zeros(size(deltaSaccade));
    for i = 1:length(chanpair)
        avgAdjustmentS = avgAdjustmentS + bsxfun(@times,deltaSaccade,PropSaccade(:,i)); %  .* repmat(PropSaccade(:,i),1,saccpts);
    end
    EEGdata(:,~markBlink) = EEGdata(:,~markBlink) - (bsxfun(@plus,avgAdjustmentS,EEGsaccademean)); % 
    %  stick original data in here so that matrix dimensions match up. the
    %  first chanpair channels of EEGdata are junk, don't need it
    EEG.data = EEGdata(length(chanpair)+1:end,:,:);
else
    subdata=EEGdata(:,:,binned_indices);
    if(rblinkpts > 0)
        for ch = 1:EEG.nbchan+length(chanpair)
            EEGblinkmean(ch) = mean(subdata(ch,rmarkBlink));
        end

        for i = 1:length(chanpair)
            deltaBlink{i} = EEGdata(i,markBlink) - EEGblinkmean(i);
        end
        for ch = 1:EEG.nbchan
            avgAdjustmentB = zeros(1,blinkpts);
            corch = ch+length(chanpair);
            for i = 1:length(chanpair)
                avgAdjustmentB = avgAdjustmentB+ PropBlink(corch,i) * deltaBlink{i};
            end
            EEG.data(ch,markBlink) = EEG.data(ch,markBlink) - (avgAdjustmentB + EEGblinkmean(corch));
        end
    end
    for ch = 1:EEG.nbchan+length(chanpair)
        EEGsaccademean(ch) = mean(subdata(ch,~rmarkBlink));
    end
    for i = 1:length(chanpair)
        deltaSaccade{i} = EEGdata(i,~markBlink) - EEGsaccademean(i);
    end
    for ch = 1:EEG.nbchan
        avgAdjustmentS = zeros(1,saccpts);
        corch = ch+length(chanpair);
        for i = 1:length(chanpair)
            avgAdjustmentS = avgAdjustmentS + PropSaccade(corch,i) * deltaSaccade{i};
        end
        EEG.data(ch,~markBlink) = EEG.data(ch,~markBlink) - (avgAdjustmentS + EEGsaccademean(corch));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Create EEG.emcp with info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG.emcp.markblink = markBlink;
EEG.emcp.blink_factors = PropBlink;
EEG.emcp.saccade_factors = PropSaccade;
EEG.emcp.regress_blinks = rblinkpts;
EEG.emcp.regress_saccades = rsaccpts;
EEG.emcp.regress_propblinks = rblinkpts/(rblinkpts+rsaccpts);
EEG.emcp.overall_blinks = blinkpts;
EEG.emcp.overall_saccades = saccpts;
EEG.emcp.overall_propblinks = blinkpts/(blinkpts+saccpts);
disp(['gratton_emcp(): Number of epochs per bin: ' num2str(EEG.emcp.bintotals) ', total of ' num2str(EEG.trials) ' epochs corrected'])
disp(['gratton_emcp(): Finished in ' num2str(toc,'%.2f') 's real time, ' num2str(cputime-starttime,'%.2f') 's cpu time']);
disp('gratton_emcp(): See EEG.emcp.table for table of correction factors')
% can use this code to show what points are identified as blinks
%EEG.data(1,markBlink) = 100;
%EEG.data(1,~markBlink) = 0;