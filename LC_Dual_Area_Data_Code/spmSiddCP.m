function spmSiddCP(func)

%  Sidd: Modified from: spmSiddOdd.
%  This is the spm file for task 740.
%
%  Extract and organize data for new task - reward based change point.
%  Task - visually guided saccade task with 2 TG's.
%  No Cue other than FP OFF!
%  These is NOT a reaction time task!
%  "Correct" targets are rewarded as: 20-80, 80-20 or 40-40.
%
%  061016 - Sidd - Started.
%  063016 - Sidd - tweaked and updated.
%  110716 - Sidd - Changed start and end offset from -1, -1 to 0, 0.
% 
%  Added the following ECode (Sidd: 070116):
%  "choice" - this is the TG that the mnky chose.
%  "trialType" - regular or boost rew + BEEP.
%  "rew_or_not" - Reward is probabilistic - so need to check if a given choice was rewarded or not, even if "correct".
%  fix_off (is the fixation point off time)

%  Origins:
%
%  Modified from:
%  function spmRKcmd(func)

%  File describing how to interpret data for FIRA
%  Use input argument "func" to switch functions:
%  'init'      ... fills FIRA.spm
%  'trial'     ... parse current trial (in FIRA.trial)
%  'cleanup'   ... whatever you need to do at the end

%  Returns:
%  nada, but can change FIRA

%  Origin notes:
%  created by jig 10/21/04
%  modified from Sharath's spm724rg.m
%  LD 2-21-2007
%  modified from Long's spmLDdots.m
%  RK 3-19-2007

try
    global FIRA EC
    declareEC_ecodes_RKLCPupil;
    ScoreList = [EC.CORRECTCD EC.WRONGCD EC.NOCHCD EC.BRFIXCD];
    
    % if called with empty FIRA.spm, fill it and return
    
    if strcmp(func, 'init')
        
        % useful ecode markers for tagged ecodes
        cb = 7000;
        cm = 6000;
        cx = 8000;
        
        % FIRA.spm is a struct with fields corresponding to "dataTypes"
        % that have values that are argument lists to the given dataType
        % constructor method (empty means use defaults).
        % (or, if you prefer, an nx2 cell array, with the first column
        % corresponding to the fieldnames, the second to the arglist.)
        %
        % Thus, to create FIRA with dataType "ecodes" including a list of
        % computed values called "one" and "two", use:
        % spm.ecodes = {'compute', {'one' 'value'; 'two' 'value}}
        %
        % buildFIRA_init then converts these descriptors into actual dataType
        % objects used to parse data. Note that the arglists to buildFIRA_init
        % can override the choices of dataTypes here (so you can use the same
        % base spm file to create different FIRAs)
        %
        % see FIRA/build/dataTypes for possible dataTypes
        %     FIRA/build/buildFIRA_init.m for details of how the spm struct is
        %       used
        
        FIRA.spm = struct( ...
            'trial',  {{ ...
            'startcd', 1005,                             ...
            'startcd_offset', -1,                    ... % offset to start code
            'endcd_offset', 0,                      ... % offset to end code
            'anycd',   [EC.CORRECTCD EC.WRONGCD EC.NOCHCD EC.BRFIXCD], ...
            'allcd',   [] ...
            'textra',  0, ...
            }}, ...
            'ecodes', {{  ...
            'extract', {   ...  % values to extract & save in FIRA{2}: <name> <type> <method> <args>
            'fp_on'       'time'      @getEC_time         {EC.FPONCD      1}; ...
            'eyefix'       'time'      @getEC_time         {EC.EYINWD      1}; ...
            'tg1_on'      'time'      @getEC_time         {EC.TRGC1CD     1}; ...
            'tg2_on'      'time'      @getEC_time         {EC.TRGC2CD     1}; ...
            'tgt_off'       'time'      @getEC_time        {EC.TARGOFFCD   1}; ...
            % 'RT_task'  'time'    @getEC_time         {EC.QOFFCD   1}; ... % Used in (short-lived) oddball task.
            'fp_off'        'time'      @getEC_time        {EC.FPOFFCD     1}; ...
            'targ_aq'      'time'      @getEC_time        {EC.TRGACQUIRECD 1}; ...
            'beep_on'     'time'      @getEC_time        {EC.BEEPON       1}; ...
            'stim_time'   'time'      @getEC_time        {EC.ELESTM       1}; ...
            'OLscore'      'value'     @getEC_fromList   {ScoreList}; ...
            'fp_x'           'value'     @getEC_tagged     {EC.I_FIXXCD    cb    cm    cx  0.1};   ...
            'fp_y'           'value'     @getEC_tagged     {EC.I_FIXYCD    cb    cm    cx  0.1};   ...
            't1_x'           'value'     @getEC_tagged     {EC.I_TRG1XCD   cb    cm    cx  0.1};   ...
            't1_y'           'value'     @getEC_tagged     {EC.I_TRG1YCD   cb    cm    cx  0.1};   ...
            't2_x'           'value'     @getEC_tagged     {EC.I_TRG2XCD   cb    cm    cx  0.1};   ...
            't2_y'           'value'     @getEC_tagged     {EC.I_TRG2YCD   cb    cm    cx  0.1};   ...
            'taskid'         'id'        @getEC_tagged      {EC.I_TASKIDCD  cb    cm    cx  1};     ...
            'trialid'         'id'        @getEC_tagged      {EC.I_TRIALIDCD cb    cm    cx  1};     ...
            'seed_base' 'value'     @getEC_tagged      {EC.I_DTVARCD   cb    cm    cx  1};     ...
            }, ...
            'compute', { ...                % names/types of fields to compute & save later in FIRA{2}
            'correct'       'id'; ...          % offline score: 0=no saccade, -1=nc, -2=brfix
            'fix_time'      'value' ; ...   % time to attain fixation from beginning of trial (first fp_on)
            'choice'         'value' ; ...
            'RT'              'value' ; ...
            'sac1_time'    'time'  ; ...
            'sac1_length'  'value' ; ...
            'sac1_vmax'    'value' ; ...
            'sac1_endx'    'value' ; ...
            'sac1_endy'    'value' ; ...
            'sac2_time'     'time'  ; ...
            'sac2_length'  'value' ; ...
            'sac2_vmax'    'value' ; ...
            'sac2_endx'     'value' ; ...
            'sac2_endy'     'value' ; ...
            'sac3_time',    'value' ; ...
            'sac3_length'   'value' ; ...
            'sac3_vmax'    'value' ; ...
            'rew1_time',    'time'  ; ...
            'rew2_time',    'time'  ; ...
            'rew3_time',    'time'  ; ...
            'analog_flag'    'id'; ...
            'trialType'        'value'; ...
            'rew_or_not'        'value'; ...
            % 'corrTG'           'value'; ...
            }, ...
            'tmp', { ...  % values to extract & use temporarily: <name> <type> <method> <args>
            }}}, ...
            'spikes',  [], ...
            'analog', [] );
        
        %% Parse trials:
        
    elseif strcmp(func, 'trial')
        
        % Get trial stuffings:
        
        % get this trial index
        tti = FIRA.raw.trial.good_count;
        
        % What kinda task (1, 2, 3...etc;  3 is the CP task we want).
        % These should all be 3 only...
        task_id = getFIRA_ec(tti, {'taskid'});
        
        % What kinda trial/block? For TG1:TG2 reward probability:
        % This decides the block type...
        % (300 is 80:20; 302 is 20:80; 304 is 40:40).
        trial_id = getFIRA_ec(tti, {'trialid'});
        
        % Get FP_OFF time:
        fpOff  = getFIRA_ec(tti, 'fp_off'); % FP Off.
        
        %         % Use the following to decide which was the randomly chosen "correct" target:
        %         tg1on  = getFIRA_ec(tti, 'tg1_on'); % "Correct" was tg1.
        %         tg2on  = getFIRA_ec(tti, 'tg2_on'); % "Correct" was tg2.
        
        % Now get to testing the trial outcome:
        
        if ~isnan(fpOff) % Good trials must have FP off...
            
            % Idiot check:
            if trial_id==300 || trial_id==302 || trial_id ==304
                
                % Got eyes?
                aNames = FIRA.analog.name;
                hidx =  find(strncmp('HEye',aNames,4) == 1);
                vidx =  find(strncmp('VEye',aNames,4) == 1);
                HEye = FIRA.analog.data(tti,hidx).values;
                VEye = FIRA.analog.data(tti,vidx).values;
                
                % Parse saccades (uses "getFIRA_saccadesCP"):
                
                %  getFIRA_saccadesCP
                %  Returns an nx8 vectors, rows are saccades found, columns are:
                %  1 .. latency
                %  2 .. duration (ms)
                %  3 .. maximum speed
                %  4 .. average speed
                %  5 .. end point x
                %  6 .. end point y
                %  7 .. distance raw  (the actual length of the saccade, calculated sample-by-sample)
                %  8 .. distance vect (the length of the saccade calculated as the linear vector)
                
                setFIRA_ec(tti, 'fp_x', 0);
                setFIRA_ec(tti, 'fp_y', 0);
                fpx = getFIRA_ec(tti, {'fp_x'});
                fpy = getFIRA_ec(tti, {'fp_y'});
                [sacs, bf,ft] = getFIRA_saccadesCP(tti, 3, true , HEye, VEye, 'fp_off', 'fp_x', 'fp_y', 6000);
                
                % Never really used - remove...
                %                 % set fix_start time
                %                 if ~isempty(ft)
                %                     setFIRA_ec(tti, 'fix_start', getFIRA_ec(tti, 'fp_off') - ft);
                %                 end
                
                % Check saccades
                % Broken Fix:
                if bf
                    % broken fixation ... just punt
                    setFIRA_ec(tti, 'correct', -2);
                    % No Saccades:
                elseif size(sacs, 1)==0
                    
                    % No saccade is error for all CP trials.
                    if task_id == 3 % CP task, all trials shd have saccs.
                        setFIRA_ec(tti, 'correct', 0);
                        % else
                        % setFIRA_ec(tti, 'correct', 1);
                    end
                    
                else
                    % One Saccade:
                    % Got at least one saccade
                    % Set parameters, saccade #1
                    setFIRA_ec(tti, 'RT',          sacs(1,1));
                    if isfinite(sacs(1,1))
                        setFIRA_ec(tti, 'sac1_time',   sacs(1,1)+fpOff);
                    end
                    if ~isfinite(sacs(1,1))
                        setFIRA_ec(tti, 'sac1_time',   sacs(end)+fpOff);
                    end
                    setFIRA_ec(tti, 'sac1_length', sacs(1,8));
                    setFIRA_ec(tti, 'sac1_vmax',   sacs(1,3));
                    setFIRA_ec(tti, 'sac1_endx',   sacs(1,5));
                    setFIRA_ec(tti, 'sac1_endy',   sacs(1,6));
                    
                    % Set parameters, saccade #2
                    if size(sacs,1) >= 2
                        setFIRA_ec(tti, 'sac2_time',   sacs(2,1)+fpOff);
                        setFIRA_ec(tti, 'sac2_length', sacs(2,8));
                        setFIRA_ec(tti, 'sac2_vmax',   sacs(2,3));
                        setFIRA_ec(tti, 'sac2_endx',   sacs(2,5));
                        setFIRA_ec(tti, 'sac2_endy',   sacs(2,6));
                        
                        if size(sacs, 1) >= 3
                            setFIRA_ec(tti, 'sac3_time',   sacs(3,1)+fpOff);
                            setFIRA_ec(tti, 'sac3_length', sacs(3,8));
                            setFIRA_ec(tti, 'sac3_vmax',   sacs(3,3));
                        end
                    end
                    
                    % For CP task:
                    % Check if saccade is within 2 deg of either target.
                    
                    % Check if saccade occured within 100-1000 ms after fpoff.
                    % FOR NOW - this is turned off!!! 070116 - Sidd.
                    
                    sacLim = 5; % Use 5 degrees for training (Sidd: 041116).
                    
                    % What was mnky's choice?
                    tg1X = getFIRA_ec(tti, 't1_x'); tg1Y = getFIRA_ec(tti, 't1_y');
                    tg2X = getFIRA_ec(tti, 't2_x'); tg2Y = getFIRA_ec(tti, 't2_y');
                    
                    % Get distance of eye from target:
                    tdist1 = sqrt((tg1X-sacs(1,5)).^2 + (tg1Y-sacs(1,6)).^2);
                    tdist2 = sqrt((tg2X-sacs(1,5)).^2 + (tg2Y-sacs(1,6)).^2);
                    
                    % Set the choice ECode:
                    % ie which TG did the mnky pick?
                    % 1 for t1; 2 for t2; 0 neither.
                    if tdist1 < sacLim % && sacs(1,1) >= 100 && sacs(1,1) <= 1000.
                        setFIRA_ec(tti, 'choice', 1);
                    end
                    if tdist2 < sacLim % && sacs(1,1) >= 100 && sacs(1,1) <= 1000.
                        setFIRA_ec(tti, 'choice', 2);
                    end
                    if tdist1 > sacLim && tdist2 > sacLim % && sacs(1,1) >= 100 && sacs(1,1) <= 1000.
                        setFIRA_ec(tti, 'choice', 0);
                    end
                    
                end % End of the distance check loop.
            end % Idiot check for trial ID.
        end % Check FP OFF.
        
        % end
        
        % **************************************
        % **************************************

        % Testing this 10/27/16 - Sidd
        
        % Reward delivery times, if rewarded:
      
        %         if isfield(FIRA, 'dio') && ~isempty(FIRA.dio{tti})
        %             % rts = find(FIRA.dio{tti}(:,3)==2,3);
        %             rts = find(FIRA.dio{tti}(:,3)==3,4); % 3 and 4 are the only values in here; Sidd: 070116.
        %             ri = getFIRA_ec(tti, 'rew1_time');
        %             for rr = 1:length(rts)
        %                 FIRA.ecodes.data(tti,ri+rr-1) = FIRA.dio{tti}(rts(rr),1);
        %             end
        %         end
        
        if ~isempty(FIRA.dio{tti})
            rewTime = FIRA.dio{tti}(1,1);
            setFIRA_ec(tti, 'rew1_time',rewTime); % Trial has regular reward.
        end
        if isempty(FIRA.dio{tti})
            setFIRA_ec(tti, 'rew1_time',nan); % Trial has regular reward.
        end
        
 
        % *************************************
        % *************************************
        
        % Set "rew_or_not" ECode and also check for BEEP + reward boost:
        if isfield(FIRA, 'dio') && ~isempty(FIRA.dio{tti})
            
            setFIRA_ec(tti, 'rew_or_not', 1);
            
            di = FIRA.dio{tti};
            if size(di,1) > 4
                extraBeep = 1;
            end
            if size(di,1) == 4
                extraBeep = nan;
            end
            
            % Set "trialType" ECode for BEEP + reward boost:
            if isnan(extraBeep)
                setFIRA_ec(tti, 'trialType',1); % Trial has regular reward.
            end
            if ~isnan(extraBeep)
                setFIRA_ec(tti, 'trialType',2); % Task has BEEP + boost reward.
            end
        end
        
        if isfield(FIRA, 'dio') && isempty(FIRA.dio{tti})
            setFIRA_ec(tti, 'rew_or_not', 0);
        end
        
        
        %% CHECK to see if the trial is a continuous recording!
        %% Some data seem continuous before toggle was set, these data
        %% are probably correct SSD or VGS trials which went 'til end
        %% anyway
        if FIRA.analog.data(tti,3).length < getFIRA_ec(tti, {'trial_end'})
            setFIRA_ec(tti, {'analog_flag'}, 0)
        elseif FIRA.analog.data(tti,3).length >= getFIRA_ec(tti, {'trial_end'})
            setFIRA_ec(tti, {'analog_flag'}, 1)
        end
        
        % %cleanup
    else
        
    end
    
catch
    evalin('base', 'e=lasterror')
end



