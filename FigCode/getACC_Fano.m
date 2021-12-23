function [ACC_sp_] = getACC_Fano(siteData,mm,minNumTrials)
%
% function [pair_dat_] = getLC_Pair_rsc_xcor(siteData)
%
% For ACC data.
%
% Fixation task ONLY.
% Calculates spike counts for different epochs relative to beep.

% ********************************************
% Summary of standard data structure:

%  siteData{1}: trialsxcols matrix, cols are:
%   1 ... fix start time wrt fixation on
%   2 ... fix end time wrt fix start time (fix duration)
%   3 ... reported correct
%   4 ... beep on time (when appropriate), wrt to fix start time
%   5 ... trial begin time, wrt to fix start time
%   6 ... trial end time, wrt to fix start time
%   7 ... trial wrt time (cpu clock)
%   8 ... LFP index corresponding to fix start time (coded above)
%   9 ... ELESTM on time (when appropriate), wrt fix start time

%  siteData{2}: Analog:
%   dim1: trial
%   dim2: sample
%   dim3: 1 = x, 2 = y, 3 = z-pupil, 4 = corrected z-pupil, 5 = pupil slope
%   [remember first sample is considered time=0 (i.e., wrt fix start time)]
%     = eyedat(Lgood,:,:);
%
%  siteData{3}: spikes, re-coded wrt fix start time

%  siteData{4}: LFP
%  Should now be 9 channels of LFP: one from LC and 8 from ACC.

%  siteData{5}: pupil events
%   1. trial number
%   2. start time of event (wrt fix start time)
%   3. end time of event (wrt fix start time)
%   4. magnitude at start of event (raw z-score)
%   5. magnitude at end of event (raw z-score)
%   6. magnitude at start of event (corrected z-score)
%   7. magnitude at end of event (corrected z-score)
%   8. time of subsequent max slope
%   9. magnitude of subsequent max slope (corrected z/sample)

%  siteData{6}: microsaccades
%   1. trial number
%   2. start time of event (wrt fix start time)
%   3. duration of event (wrt fix start time)
%   4. maximum velocity (deg/ms)
%   5. magnitude of microsaccade event (deg)
%   6. onset time wrt phase of associated pupil event (fraction)
%   7. magnitude of associated pupil event

%  siteData{7}: Spike and analog signal channels
%   1. spike channel numbers
%   2. Analog channel names (LFP's, eye signals, eeg, pulse-ox)

% Sidd: Added the two additional cells below (062116).

%  siteData{8}: EEG

%  siteData{9}: Pulse-Ox

% ********************************************

%%

% bs = linspace(100,1000,10); % 032220b.
bsiz = 100; % April 2021.
backTime = -1000; fwdTime = 1000; % 032320b.
tmin   = backTime; tmax   = fwdTime; xBin = []; nb = []; nbs = length(bsiz);
tsize  = bsiz; tstep  = tsize/2;
xBin  = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
nb = size(xBin,1);

winSiz = length(backTime:fwdTime)/1000;

%%

bsiz2 = 200; % April 2021.
backTime2 = -1000; fwdTime2 = 1000; % 032320b.
tmin2   = backTime2; tmax2   = fwdTime2; xBin2 = []; nb2 = []; nbs2 = length(bsiz2);
tsize2  = bsiz2; tstep2  = tsize2/2;
xBin2 = [(tmin2:tstep2:tmax2-tsize2)' (tmin2+tsize2:tstep2:tmax2)'];
nb2 = size(xBin2,1);

winSiz2 = length(backTime2:fwdTime2)/1000;

%%

ccgbSize = 1;
% binsX = -1100:ccgbSize:1100; % Bins for creating binary spike train for 1st 1 sec of stable fixation after beep.
binsX3 = -1000:ccgbSize:1000; % Bins for creating binary spike train for 1st 1 sec of stable fixation after beep.

%% Find no beep trials ...
LnoBeep    = isnan(siteData{1}(:,4)); % Only non-beep trials
FnoBeep    = find(LnoBeep);
num_noBeep = length(FnoBeep);

LBeep    = ~isnan(siteData{1}(:,4)); % Only non-beep trials
FBeep    = find(LBeep);
num_Beep = length(FBeep);

beepTime = siteData{1}(LBeep,4);

%% Get pupil data ...

% .... or not.

%%  Get unit information for this session:

sp_ch = siteData{7}{1};     % All spike channel id's.
ACC_ch = find(sp_ch>8000);   % ACC channel id's.

% How many units in LC?
n_ACC = length(ACC_ch);

% n_ACC








%% Get timing information for this session:

% fst = siteData{1}(LnoBeep,1); % Fixation start time (wrt fix ON).

% Fixation duration.
fd = siteData{1}(LnoBeep,2);
fd_beep = siteData{1}(LBeep,2);

%% Get LC binary spike trains and trial spike counts.

nL_count = cell(1,n_ACC);
nL_count1 = cell(1,n_ACC);
nL_count2 = cell(1,n_ACC);
nL_count3 = cell(1,n_ACC);

% Get LC trial spikes..
sps_ACC = siteData{3}(LnoBeep,ACC_ch); % LC spikes.
sps_ACC_beep = siteData{3}(LBeep,ACC_ch); % LC spikes.

% Trial indices for fixation duration > 3sec.
gti = find(fd>=2000);
% gti_beep = find(beepTime>=2000 & fd_beep-beepTime>=1000);
gti_beep = find(beepTime>=1000 & fd_beep-beepTime>=1000);
beepTime = beepTime(gti_beep);

%%

% Get the spike times for each trial and create binary spike trains.
for ai = 1:n_ACC
    ats = sps_ACC(gti,ai);
    nL = nans(1,length(gti));
    for ti = 1:length(gti)
        ts = []; ts = ats{ti}; % Spike times were recoded wrt mnky fixation start.
        ts_fix = ts(ts>1001 & ts<=2000); % Get spikes from the 1st 2sec of fixation.
        nL(ti) = length(ts_fix);
    end % Trials loop.
    
    if ti>=minNumTrials
        nL_count{ai} = nL;
    end
    if ti<minNumTrials
        nL_count{ai} = [];
    end
    
end

%%

binsX = xBin;
binsX2 = xBin2;

FF_unit = nans(n_ACC,nb2);
% FF_unit = nans(n_ACC,2001);
nA = cell(1,n_ACC);
nA1 = cell(1,n_ACC);
for ai = 1:n_ACC
    ats = sps_ACC_beep(gti_beep,ai);
    nL_bef = nans(1,length(gti_beep));
    nL_phasic = nans(1,length(gti_beep));
    nL_aft_phasic = nans(1,length(gti_beep));
    bs = nans(length(gti_beep),length(binsX)); % Preallocate.
    bs2 = nans(length(gti_beep),length(binsX2)); % Preallocate.
    
    for ti = 1:length(gti_beep)
        
        ts = [];
        ts = ats{ti}-beepTime(ti); % Spike times were recoded wrt mnky fixation start.
        ts_fix = ts(ts>=-1000 & ts<= 1000); % Get spikes from the 1sec after beep.
        
        % 1ms bins for 2D plots!
        bsTemp = hist(ts_fix,binsX3); % Spikes binned at 1msec.
        bs_1(ti,:) = bsTemp;
        
        for bi = 1:nb % Iterate over bins.
            bsTemp = sum((ts_fix > binsX(bi,1) & ts_fix < binsX(bi,2)));
            bs(ti,bi) = bsTemp;
        end
        
        for bi = 1:nb2 % Iterate over bins.
            bsTemp2 = sum((ts_fix > binsX2(bi,1) & ts_fix < binsX2(bi,2)));
            bs2(ti,bi) = bsTemp2;
        end
        
        ts_fix_bef = ts(ts>=-1000 & ts <0); % Get spikes from before beep period.
        ts_fix_phasic = ts(ts>=0 & ts <200); % Get spikes from after beep period.
        ts_fix_aft_phasic = ts(ts>=100 & ts <1100); % Get spikes after beep period, but from after phasic response.
        
        nL_bef(ti) = length(ts_fix_bef);
        nL_phasic(ti) = length(ts_fix_phasic);
        nL_aft_phasic(ti) = length(ts_fix_aft_phasic);
        
    end % Beep trials loop.
    
    nA{ai} = bs; % Binned ACC spikes for this unit - ie - binary spike trains.
    nA1{ai} = bs_1; % Binned ACC spikes for this unit - ie - binary spike trains.
    
    % pause
    
    numsp = sum(bs,2);
    g_ind = find(numsp>1);
    bs = bs(g_ind,:);
    if length(g_ind)>=minNumTrials
        nL_count1{ai} = nL_bef;
        nL_count2{ai} = nL_phasic;
        nL_count3{ai} = nL_aft_phasic;
        
        % msc = bs;
        msc = bs2;
        %         msc = movsum(bs,2,2);
        mean_sc = nanmean(msc);
        var_sc = nanvar(msc);
        unit_FF = nans(1,length(mean_sc));
        finite_ind = find(mean_sc>0 &var_sc>0);
        if length(finite_ind)>minNumTrials
            unit_FF(finite_ind) = var_sc(finite_ind)./mean_sc(finite_ind);
        end
        % infinite_ind = find(mean_sc==0);
        % unit_FF(infinite_ind) = nan;
        
        % Don't smooth here.
        FF_unit(ai,:) = unit_FF;
        
    end
    
    if ti<minNumTrials
        nL_count1{ai} = [];
        nL_count2{ai} = [];
        nL_count3{ai} = [];
    end
    
end % ACC neuron loop.

%%

ACC_sp_ = {nL_count nL_count1 nL_count2 nL_count3 nA FF_unit nA1};


% % Test.
% taxis = -1000:1000;
% meanFF = nanmean(FF_unit);
% seFF = nanse(FF_unit,1);
% figure; hold on;
% errorBarFilled(taxis,meanFF,seFF,[.5 .5 .5]);
% plot(taxis,meanFF,'k-','linewidth',2)

% % pause
% % close
