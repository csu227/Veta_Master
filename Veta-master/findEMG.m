function findEMG(filename)
% Authors: Nick Jackson (njackson@uoregon.edu) & Ian Greenhouse
% By Cazmon & Jessica  

% 
% This function identifies various EMG/TMS events in prerecorded EMG data.
% 
% The user defines appropriate data parameters either at the command line
% or in the first section of the code. Analysis parameters are defined in 
% the second section, and should be adjusted according to the data.
% 
% Event metrics will be added to the trials table depending on the type of
% data collected. The following shows the type of event detection and the
% associated metrics that are outputed:
% 
% If EMG (electromyography) burst detection is on:
% 
% ch#_EMGburst_onset = beginning of EMG burst
% ch#_EMGburst_offset = end of EMG burst
% ch#_EMG_RT = reaction time relative to photodiode event
% ch#_EMGburst_area = area under the EMG burst curve
% 
% If MEP (motor evoked potential) detection is on:
% 
% artloc = TMS artefact location
% ch#_MEP_time = MEP onset
% ch#_MEP_latency = time from TMS artefact to MEP time
% ch#_MEP_duration = length of MEP
% ch#_RMS_preMEP = root mean square tolerance, determines if pre-MEP noise
% may contaminate resting MEP measurements
% ch#_MEP_amplitude = peak-to-peak amplitude of MEP
% ch#_MEP_area = area under the MEP curve
% 
% If CSP (cortical silent period) detection is on:
% 
% ch#_CSP_onset = beginning of CSP
% ch#_CSP_offset = end of CSP
% 
% output:
% Once the analysis code has run, a UI will open and prompt the user to enter
% a file name. The default file name is the original file name appended with
% "preprocessed". The saved file will contain the following:
%     parameters: struct
%         analysis parameters
%     subject : struct
%         generated from recordEMG
%     trials : trials
%         updated trials table with addition columns that contain calculated
%         metrics itemized above
%% define data parameters
%index 101 corresponds to time 0(.1s if starting at 0 and not negatives)
%which is time of stimulus
%%start mep search at .115 and end at .165

time_matlab = [0:.001:.399]';
time_signal = [-.1:.001:0.299]';

time = [time_matlab time_signal];

use_command_line = 1;%toggle bool to suit parameter input preferences. IF this is set to 0 it wont ask user to define parameters
if ~use_command_line
%     parameters.EMG = 1; % Detect EMG bursts: 0 = no, 1 = yes
%     parameters.EMG_burst_channels = [1 2];
    parameters.MEP = 1; % Detect MEPs: 0 = no, 1 = yes
    parameters.artchan_index = 3;
    parameters.MEP_channels = [1];
    
 %Silent Period   
    
    parameters.CSP = 1; % Detect CSP: 0 = no, 1 = yes
     parameters.CSP_channels = [0];
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     parameters.MEP_std_or_chngpts = 1; % 0 = std of baseline, 1 = findchangepts
end
%% define analysis parameters
%creates a structure for just parameters. Can set a structure by making a
%variable with anything after a period

%edit these to suit your analysis needs
parameters.sampling_rate = 1000; % samples per second (Hz) &chnaged from 5000 to 100 for our study 
% parameters.emg_burst_threshold = .3; % raw threshold in V to consider for EMG
% parameters.emg_onset_std_threshold = 2; % number of std to consider for EMG burst onsets/offsets
% parameters.tms_artefact_threshold = .0001; % raw threshold magnitude in V to consider for TMS artefact

%define MEP search range
parameters.min_TMS_to_MEP_latency = .015; % number of secs after TMS to begin MEP onset detection  % lower_limit_MEP_window
parameters.MEP_window_post_artefact = .065; % time in s after TMS to measure MEP in seconds % upper_limit_MEP_window

parameters.pre_TMS_reference_window = .05;  % time before TMS to serve as reference baseline for MEP onset



parameters.MEP_onset_std_threshold = .5; % number of std to consider for MEP onsets


% parameters.RMS_preMEP_EMG_tolerance = .05; % root mean square EMG tolerance for including MEP

% parameters for removing TMS artefact & MEP when detecting EMG in same channel as MEP
%parameters.time_prior_to_TMS_artefact = .005; % time in s prior to TMS artefact to set to zero %Changed to .1s or 100ms for time before stimulus in our files



parameters.end_of_MEP_relative_to_TMS = .1; % time in s of end of MEP relative to TMS artefact


%% Parse input
if (nargin < 1)
    % open file with finder/file explore
    [filename, pathname] = uigetfile;
    File = fullfile(pathname, filename);
else
    File = fullfile(pwd, filename);
end
load(File);
%% set number of channels
trials.trial_accept(:,1) = 1;
chs = ["ch1","ch2","ch3","ch4","ch5","ch6","ch7","ch8"];
parameters.num_channels = sum(contains(trials.Properties.VariableNames,chs));
%% Toggles EMG and TMS functionality
if use_command_line

    
    parameters.MEP = input('Do you want to detect motor evoked potentials (MEPs)? yes(1) or no(0): ');
    parameters.CSP = input('Do you want to detect cortical silent period (CSP) epochs? yes(1) or no(0): ');

    if parameters.MEP | parameters.CSP %  "|" = "=="
        trials.artloc(:,1) = 0;
        trials.preTMS_period_start(:,1) = 0;
        
%Commented out the question for artifact channel        
%         if ~isfield(parameters,'artchan_index')
%             parameters.artchan_index = input('Enter TMS artefact channel #: ');
%         end
    end
    
    if parameters.MEP & ~isfield(parameters, 'MEP_channels')
        parameters.MEP_channels = input('Enter MEP channels (e.g. [2] or [1 3 5]): ');
    end
       if parameters.CSP & ~isfield(parameters,'CSP_channels')
        parameters.CSP_channels = input('Enter CSP channels (e.g. [2] or [1 3 5]): ');
    end
    
    parameters.MEP_std_or_chngpts = 1; % 0 = std of baseline, 1 = findchangepts
end
%% find photodiode event %% replace for pinchmeter later
% if any(strcmp('photodiode', trials.Properties.VariableNames))
%     trials = findDiode(trials,parameters);
% else
%     trials.stim_onset(:,1) = zeros;
% end
%% Find TMS, MEP, and EMG events
trials = findEvents(trials,parameters); %findEvents() function is nested at the bottom 
%% save file
outfile=[File(1:end-4),'_preprocessed'];
uisave({'trials','subject','parameters'},outfile);
end
%% HELPER FUNCTIONS
%% Find Diode
% function trials=findDiode(trials,parameters)
% % detects large changes in photodiode signal
%     for i=1:height(trials)
%         %% photodiode event
%         diff_diode = diff(trials.photodiode{i});
%         [max_diode_value,max_diode_index] = max(diff_diode);
%         trials.stim_onset(i,1) = max_diode_index/parameters.sampling_rate; % location for GUI
%     end
% end

%% Find TMS, MEP, CSP, and EMG
function trials = findEvents(trials,parameters) % parameterize these

%% initialize trials columns

% commenting out artchan index
% if parameters.artchan_index
%     trials.artloc(:,1) = 0;
% end

if parameters.MEP
    for chan = 1:length(parameters.MEP_channels)
        trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_onset'])(:,1) = 0; % adding to structure for trials. middle portion just is calling a previous structure to get value for 'ch' and adding outcome variables
        trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_latency'])(:,1) = 0;                                                                                                                                                                 
       % trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_offset'])(:,1) = 0;
        trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_area'])(:,1) = 0;
    end
end

if parameters.CSP
    for chan = 1:length(parameters.CSP_channels)
        trials.(['ch', num2str(parameters.CSP_channels(chan)), '_CSP_onset'])(:,1) = 0;
        trials.(['ch', num2str(parameters.CSP_channels(chan)), '_CSP_offset'])(:,1) = 0;
    end
end

%% identify MEP and non-MEP channels
%dont think we need this since all channels are MEP Channels
% if parameters.MEP & parameters.EMG
%     non_MEP_channels = parameters.EMG_burst_channels(parameters.EMG_burst_channels ~= parameters.MEP_channels); % do not want to detect MEPs as EMG bursts, so ignore MEP channel
% elseif parameters.EMG
%     non_MEP_channels = parameters.EMG_burst_channels;
% end

%% sweep loop
for i = 1:height(trials)    
    %% find TMS artefact and MEP
    if parameters.MEP
        for chan = 1:length(parameters.MEP_channels)
            MEPchannel = trials.(['ch', num2str(parameters.MEP_channels(chan))]){i,1}; %pulls first sweet from first channel

% Commented out the use of artifact channel. Setting an index location. Currently 101 (time = 0)          
            
%             if parameters.artchan_index
%                 artchannel = trials.(['ch', num2str(parameters.artchan_index)]){i,1};
%                 [artefact_value, TMS_artefact_sample_index] = max(artchannel); %max(abs(artchannel));
%             else
%                 TMS_artefact_sample_index = find(MEPchannel > parameters.tms_artefact_threshold,1);
%                 artefact_value = abs(MEPchannel(TMS_artefact_sample_index));
%             end                                    


%Test to try a set value for the artefact index in this case where to start
%searching for MEP
            TMS_artefact_sample_index = 115;

            %set range to look for preMEP EMG activity to calculate RMS
            baseline_lower_bound = 41;% index 41 = cursor 3 %TMS_artefact_sample_index - (parameters.pre_TMS_reference_window * parameters.sampling_rate); 1 is start of trial
            baseline_upper_bound = 91; % index 91 = cursor 4 %TMS_artefact_sample_index; 100 for end of 100 ms

            %redefine range if it extends beyond lower x limit
            if baseline_lower_bound < 0
                baseline_lower_bound = 1;
            end
            %taking our rms use
            preTMS_reference_data = MEPchannel(baseline_lower_bound:baseline_upper_bound); %pulls out baseline
            %RMS_of_preMEP_window = rms(preTMS_reference_data);
            
            trials.(['ch',num2str(parameters.MEP_channels(chan)),'baselineAVG_preMEP'])(i,1) = mean(preTMS_reference_data); %stores baseline in structure
            BA = trials.(['ch',num2str(parameters.MEP_channels(chan)),'baselineAVG_preMEP'])(i,1);
            %use this later with baseline peak to peak to compare with mep
            %peak to peak to accept or reject
            % reject trial if RMS is above tolerance threshold
%             if RMS_of_preMEP_window < parameters.RMS_preMEP_EMG_tolerance  %determines whether to 
%                 trials.trial_accept(i,1)=1;
%             else
%                 trials.trial_accept(i,1)=0;
%             end



 % Commented out more artifact variables           
%             if artefact_value > abs(parameters.tms_artefact_threshold) % TMS artefact must exceed a threshold to be classified as an artefact
%                 trials.artloc(i,1) = TMS_artefact_sample_index/parameters.sampling_rate; %artefact location scaled for visualization
 % changed the if statement to look for MEP to look no matter what     
 
             if    baseline_upper_bound == 91
                %define MEP search range
                %eveyrthing is already off by 1 so need to adjust by 1 and
                %then another 1 to account for values at cursors
                
                % index 116 = cursor 2 
                % index 166 = cursor 3 
                lower_limit_MEP_window = 116; % one point after the cursors to align with how spike finds the max value %* parameters.sampling_rate; %changing mep search range lower limit to 115 (15ms) %TMS_artefact_sample_index + (parameters.min_TMS_to_MEP_latency * parameters.sampling_rate);
                upper_limit_MEP_window = 166; % one point before the cursors to align with how spike finds the max value %* parameters.sampling_rate; %changing mep search range lower limit to 165 (65ms) %TMS_artefact_sample_index + (parameters.MEP_window_post_artefact * parameters.sampling_rate);
                if lower_limit_MEP_window>length(MEPchannel)
                    lower_limit_MEP_window = 1;
                end
                if upper_limit_MEP_window>length(MEPchannel)
                    upper_limit_MEP_window = length(MEPchannel);
                end
                MEPsearchrange = MEPchannel(lower_limit_MEP_window+1:upper_limit_MEP_window-1);
                
                % detect MEP onset and offset point;
                
                %use con_find
                %here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
               % MEP_onset_from_TMS = find(MEPsearchrange > parameters.MEP_onset_std_threshold * std(abs(preTMS_reference_data)),1); % first value that exceeds std threshold within rectified MEP search range
                ipoints = findchangepts(MEPsearchrange, 'MaxNumChanges', 10, 'Statistic', 'mean'); % fewer change points may suffice
                
                
                time = .001:.001:.4;
                
                
%               [upperbound,lowerbound,MCD] = MCD_Find(preTMS_reference_data);
                

                [upperbound,lowerbound,MCD] = MCD_Find(preTMS_reference_data);

                
                
               [MEP_onset_time,MEP_onset_index] = CON_Finder(MEPsearchrange,time,upperbound,'U',5,1);     %CON_Finder(EMG_wave,time,Threshold,direction,varargin)%,n,start,direction)
               
               MEP_onset_index = MEP_onset_index(1)+lower_limit_MEP_window;
                  
                if ~isnan(MEP_onset_index) %if we found onset then look for offset
                   
                    trials.trial_accept(i,1) = 1; % if first criteria worked trial accept (1)
                    offsearchrange =  MEPchannel(MEP_onset_index:end);
                  
                    
                    
                    [MEP_offset_time,MEP_offset_index] = CON_Finder(offsearchrange,time,upperbound,'D',5,1);   %CON_Finder(EMG_wave,time,Threshold,direction,varargin)%,n,start,direction)
                     MEP_offset_index = MEP_offset_index(1)+MEP_onset_index;
                
                     if isnan(MEP_offset_index)
                               MEP_offset_index = ipoints(end) + lower_limit_MEP_window;
                     else
                                MEP_offset_index = MEP_onset_index;
                     end
                
                     
                     
                else  %if our first criteria didnt find onset then use this next and find offset
                        %MEP_offset_index = NaN; 
                        trials.trial_accept(i,1) = 0; %if first criteria didnt work do not accept trial (0)
                        
                        MEP_onset_from_TMS = find(MEPsearchrange > parameters.MEP_onset_std_threshold * std(abs(preTMS_reference_data)),1); % first value that exceeds std threshold within rectified MEP search range
                         ipoints = findchangepts(MEPsearchrange, 'MaxNumChanges', 10, 'Statistic', 'mean'); % fewer change points may suffice
                
                            if parameters.MEP_std_or_chngpts & ipoints
                                 MEP_onset_index = ipoints(1) + lower_limit_MEP_window; % use findchangepts value
                            else
                                 MEP_onset_index = MEP_onset_from_TMS + lower_limit_MEP_window; % use num std of baseline
                            end
                
                            if ipoints
                                 MEP_offset_index = ipoints(end) + lower_limit_MEP_window;
                            else
                                 MEP_offset_index = MEP_onset_index;
                            end
                end
                  
%                 MEP_onset_index = find((MEPsearchrange < upperbound),1);
%                 
                
%                 
%                                
%                  MEP_offset_index = find((offsearchrange < lowerbound),1);
                
                
                %select between two options for determining MEP onset
%                 if parameters.MEP_std_or_chngpts & ipoints
%                     MEP_onset_index = ipoints(1) + lower_limit_MEP_window; % use findchangepts value
%                 else
%                     MEP_onset_index = MEP_onset_from_TMS + lower_limit_MEP_window; % use num std of baseline
%                 end
%                 
%                 if ipoints
%                     MEP_offset_index = ipoints(end) + lower_limit_MEP_window;
%                 else
%                     MEP_offset_index = MEP_onset_index;
%                 end
%                 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                % define lower limit of pre-TMS artefact reference window
                % for determining threshold for MEP.
                %preTMS_reference_window_lower_limit = TMS_artefact_sample_index - (parameters.pre_TMS_reference_window * parameters.sampling_rate);
                
                %redefine range if preTMS reference range extends beyond lower x limit
               % if preTMS_reference_window_lower_limit < 1
                %    preTMS_reference_window_lower_limit = 1;
               % end
               trials.preTMS_period_start(i,1) = baseline_lower_bound/parameters.sampling_rate;
                
                %redefine range if it extends beyond upper x limit
     
                if MEP_offset_index > length(trials.ch1{1,1})
                    MEP_offset_index=length(trials.ch1{1,1})-1;
                end
                
                
                %look only in range after artefact
                %preTMS_MEP_reference_data = MEPchannel(preTMS_reference_window_lower_limit:TMS_artefact_sample_index);
                %MEPsearchrange = MEPchannel(lower_limit_MEP_window:upper_limit_MEP_window);
                [max_MEP_value,MEP_max_sample_point] = max(MEPsearchrange);  % fix alignment to match cursors
                [min_MEP_value,MEP_min_sample_point] = min(MEPsearchrange);
                
                
                %for MEP max min it does not use values at cursors in
                %signal but for baseline it looks like it does use those
                %values so I dont need +- 1
                Baselinesearchrange = MEPchannel(baseline_lower_bound:baseline_upper_bound);
                
                [max_Baseline_value,Baseline_max_sample_point] = max(Baselinesearchrange);  % fix alignment to match cursors
                [min_Baseline_value,Baseline_min_sample_point] = min(Baselinesearchrange);
                
                
                
                % MEP AREA formula
                if ~isnan(MEP_onset_index)
                MEParea = trapz(MEPchannel(MEP_onset_index:MEP_offset_index));
                else
                    MEParea = NaN;
                end
                
                % orginal method of MEP area %sum(abs(MEPchannel(MEP_onset_index:MEP_offset_index)))/(MEP_max_sample_point - MEP_min_sample_point);
                
                % identify MEP onset
                %if ~isempty(MEP_onset_index)
                 if ~isnan(MEP_onset_index)
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_onset'])(i,1) = MEP_onset_index/parameters.sampling_rate;                   
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_offset'])(i,1) = (MEP_offset_index/parameters.sampling_rate);% - trials.artloc(i,1);
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_latency'])(i,1) =
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_duration'])(i,1) = trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_offset'])(i,1) - trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_time'])(i,1);  %changed the formula from 'MEP_latency to 'MEP_time'
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_amplitude'])(i,1) = max_MEP_value - min_MEP_value;
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_Baseline_amplitude'])(i,1) = max_Baseline_value - min_Baseline_value;
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_Peak-Peak_Subtraction'])(i,1) = trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_amplitude'])(i,1) - trials.(['ch', num2str(parameters.MEP_channels(chan)), '_Baseline_amplitude'])(i,1);
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_Trial_Accept_200'])(i,1) = (trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_Peak-Peak_Subtraction'])(i,1) > 200);
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_Trial_Accept_50'])(i,1) = (trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_Peak-Peak_Subtraction'])(i,1) > 50);
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_area'])(i,1) = MEParea;
                    
                 else
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_time'])(i,1) = NaN;                   
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_offset'])(i,1) = NaN;
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_duration'])(i,1) = NaN;
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_amplitude'])(i,1) = NaN;
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_Baseline_amplitude'])(i,1) = NaN;
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_Peak-Peak_Subtraction'])(i,1) = NaN;
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_Trial_Accept_200'])(i,1) = NaN;
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_Trial_Accept_50'])(i,1) = NaN;
                    trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_area'])(i,1) = MEParea;
                    
                end
                
            end
        end
    end    
    
     %% find CSP
     %START HERE TODAY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if parameters.CSP
        for chan = 1:length(parameters.CSP_channels) %for # of channels
                            MEPchannel = trials.(['ch', num2str(parameters.MEP_channels(chan))]){i,1}; %pulls first sweep from first channel
    preTMS_reference_data = MEPchannel(baseline_lower_bound:baseline_upper_bound); %pulls out baseline
           
   [upperbound,lowerbound,MCD] = MCD_Find(preTMS_reference_data);
            
    CSPsearchrange = MEPchannel(MEP_onset_index:end);
    
    [CSP_onset_time,CSP_onset_index] = CON_Finder(CSPsearchrange,time,lowerbound,'D',5,1);
    CSP_onset_index = CSP_onset_index + MEP_onset_index;
   % offsearchrange =  MEPchannel(CSP_onset_index:end);
    
%     for i =  CSP_onset_index:length(MEPchannel)-10
%         csp_test = MEPchannel(i:i+10);
%         proof = csp_test > BA;
%         if sum(proof => 5)
%             CSP_offset_index = i;
%         break
%         end
%     end
    
    
    
            csp_signal = trials.(['ch', num2str(parameters.CSP_channels(chan))]){i,1};
            CSP_start_search_range =  csp_signal(MEP_offset_index+1:end);
            
            
             if ~isnan(CSP_start_index)
                    offsearchrange =  MEPchannel(CSP_onset_index:end);
                  
                    [CSP_offset_time,CSP_offset_index] = CON_Finder(offsearchrange,time,lowerbound,'U',5,1);
                     CSP_offset_index = CSP_offset_index+MEP_onset_index;
                else 
                       CSP_start_index = MEP_offset_index + find(CSP_start_search_range < lowerbound,1); %probably create a function similar to find but for consecutive values
            
                        CSP_end_search_range = csp_signal(CSP_start_index+1:end);
            
                            try
                          CSP_end_index = MEP_offset_index + find(CSP_end_search_range > upperbound,1);
                            catch
                      CSP_end_index = 400;
                             end
             end
            
            
            
            
            
            
            
           
            
            if CSP_start_index
                trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_onset'])(i,1) = CSP_start_index/parameters.sampling_rate; % DOUBLE CHECK TO SEE IF THIS GIVES THE RIGHT TIME VALUE
                try
                    trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_offset'])(i,1) = CSP_end_index/parameters.sampling_rate;  % DOUBLE CHECK TO SEE IF THIS GIVES THE RIGHT TIME VALUE
                catch
                    trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_offset'])(i,1) = NaN;%/parameters.sampling_rate;
                end
            elseif ~CSP_start_index %no silent period
                trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_onset'])(i,1) = NaN;%/parameters.sampling_rate;
                trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_offset'])(i,1) = NaN;%/parameters.sampling_rate;
            end
            
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%             if trials.artloc(i,1) %if artc loc has a non zero value
%                 csp_signal = trials.(['ch', num2str(parameters.CSP_channels(chan))]){i,1};
%                 [IUPPER, ILOWER, UPPERSUM] = cusum(abs(csp_signal));
%                 
%                 
%                 if sum(ismember(parameters.CSP_channels, parameters.MEP_channels))
%                     trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_onset'])(i,1) = trials.artloc(i,1) + ...
%                                                                                                 trials.(['ch', num2str(parameters.MEP_channels(chan)), '_MEP_offset'])(i,1);
%                     csp_start_loc = trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_onset'])(i,1) * parameters.sampling_rate;
%                 else                    
%                     [peak1, csp_start_loc] = findpeaks(UPPERSUM, 'MinPeakProminence', 9, 'NPeaks', 1);
%                     trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_onset'])(i,1) = csp_start_loc/parameters.sampling_rate;
%                 end
%                 if csp_start_loc
%                     [peak2, csp_end_position] = findpeaks(-1*UPPERSUM(round(csp_start_loc):end), 'MinPeakProminence', 5, 'NPeaks', 1); % default 5, lower may suffice
%                     csp_end_loc = csp_end_position + csp_start_loc;
%                 end
%                 if csp_end_loc
%                     trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_offset'])(i,1) = csp_end_loc/parameters.sampling_rate;
%                 end
%             else %else if art loc is zero
%                 artchannel = trials.(['ch', num2str(parameters.artchan_index)]){i,1}; %not using
%                 [artefact_value, TMS_artefact_sample_index] = max(abs(artchannel));
%                 
%                 if artefact_value > abs(parameters.tms_artefact_threshold) % TMS artefact must exceed a threshold to be classified as an artefact
%                     trials.artloc(i,1) = TMS_artefact_sample_index/parameters.sampling_rate;%artefact location scaled for visualization
%  %%%%%%%%%%%%%%%%%%%%%%%%%%This is the only part that actually looks for csp
%                     
%                     csp_signal = trials.(['ch', num2str(parameters.CSP_channels(chan))]){i,1};
%                     [IUPPER, ILOWER, UPPERSUM] = cusum(abs(csp_signal));
%                     [peak1, csp_start_loc] = findpeaks(UPPERSUM,'MinPeakProminence', 9, 'NPeaks', 1);
%                     if csp_start_loc
%                         [peak2, csp_end_position] = findpeaks(-1*UPPERSUM(csp_start_loc:end), 'MinPeakProminence', 5, 'NPeaks', 1); % default 5, lower may suffice
%                         csp_end_loc = csp_end_position + csp_start_loc;
%                     end
%                     
%                     if csp_start_loc & csp_end_loc
%                         trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_onset'])(i,1) = csp_start_loc/parameters.sampling_rate;
%                         trials.(['ch', num2str(parameters.CSP_channels(chan)) '_CSP_offset'])(i,1) = csp_end_loc/parameters.sampling_rate;
%                     end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%                 end
%             end            
%         end
%     end
 end % end trial loop

end % end findEvent function


% 
% function [upperbound,lowerbound] = MCD(x)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% MCD = mean(abs(diff(x)));
% avg = mean(x);
% 
% upperbound = avg + MCD * 2.66;
% lowerbound = avg - MCD * 2.66;
% end