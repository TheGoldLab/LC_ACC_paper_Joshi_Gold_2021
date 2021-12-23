function [beep_rsc_] = getLC_DualArea_FakeBeep_binned_rsc(siteData,numTrialsLim,fi)
%

% Mod: 041220: Sidd
% Mod: 022019: Sidd
% Streamlined so only considers case when LC has phasic with post phasic inhibition.
% Added spike rate and Fano calculation.
%
% Mod: 123018 - Sidd.
% For low and high LC spiking trials, calculate pairwise spike count correlation for ACC pairs.
% function [beep_rsc] = getLC_DualArea_beep_rsc(siteData)
%
% Fixation task, BEEP trials.
%
% Origin: 062516 - Sidd
% Calculate pairwise spike-train cross-correllograms (CCG's) between ACC neurons, conditioned on LC activity.
% CCG calculation based on Bair et al and Kohn et al.
% Each session is divided (median split) into high and low LC spike count trials.
% This median split is the same as done for the spike-count correlation calculations.
%
% Mod: 100218 - Sidd
% Now look at BEEP trials. Session beep trials are grouped based on whether the beep evoked a phasic response or not.
% Mod: 120318: now also write out the psth's.
% Calls get_ccg to calculate the CCG's.

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
%% Find no beep trials ...
LNoBeep = isnan(siteData{1}(:,4)); % Only non-beep trials.
FNoBeep    = find(LNoBeep);
num_noBeep = length(FNoBeep);

%% Set up bins ...

% Create binsizes.
tt_time = 1000;
bs = linspace(100,tt_time,10); % 021720.

% backTime = 1001; fwdTime = 2100; % We're looking at the 1sec to 2.1sec window of stable fixation.
% tmin = 1; tmax = tt_time+100;
tmin = 1; tmax = tt_time+1;
xBinA = []; nbA = [];
nbs = length(bs);
for i = 1:nbs
    tsize  = bs(i); tstep  = 0.5*tsize;
    xBinA{i}  = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
    %     if i == 10
    %         xBinA{i}  = [1 1000];
    %     end
    % tax    = mean(binsX,2);
    nbA{i} = size(xBinA{i},1);
end

% backTime = 1001; fwdTime = 2100; % We're looking at the 1sec to 2.1sec window of stable fixation.
% tmin = -tt_time-100; tmax = -1;
tmin = -tt_time-1; tmax = -1;
xBinB = []; nbB = [];
nbs = length(bs);
for i = 1:nbs
    tsize  = bs(i); tstep  = 0.5*tsize;
    xBinB{i}  = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
    %     if i == 10
    %         xBinB{i}  = [-1000 -1];
    %     end
    % tax    = mean(binsX,2);
    nbB{i} = size(xBinB{i},1);
end

%%  Get unit information for this session:

sp_ch = siteData{7}{1};     % All spike channel id's.
LC_ch = find(sp_ch<5000);   % LC channel id's.
ACC_ch = find(sp_ch>8000); % ACC channel id's.

% How many units in LC and ACC?
n_LC = length(LC_ch);
n_ACC = length(ACC_ch);

% How many pairs in ACC?
ACC_pairs = nchoosek(ACC_ch,2);
n_ap = size(ACC_pairs,1);

%% Get timing information for this session:

% fst = siteData{1}(LBeep,1); % Fixation start time (wrt fix ON).
fd = siteData{1}(LNoBeep,2); % Fixation duration.

%% Get ACC binned spike counts.

sps_LC = siteData{3}(LNoBeep,LC_ch); % LC spikes.
sps_ACC = siteData{3}(LNoBeep,ACC_ch); % ACC spikes.

% Want to drop 1st second and then get 1 sec of ACC spikes after "beep".
gti = find(fd >2100);
ngt = length(gti);
fd = fd(gti); % Stable fixation duration.

ccgbSize = 1;
binsX = -tt_time:ccgbSize:tt_time; % Bins for creating binary spike train for 1st 1 sec of stable fixation after beep.

%%

% Collect spikes from window of size winSiz - needed for calculating mean spike rate.
winSiz = (2*tt_time)/1000;

% Get trial spike trains for ACC units.
% Get the spike times for each trial and create binary spike trains.
% spCountBef = cell(1,n_ACC);
% spCountAft = cell(1,n_ACC);
nA_binnedB = [];
nA_binnedA = [];
for ai = 1:n_ACC
    ats = sps_ACC(gti,ai);
    for bsi = 1:nbs % Iterate over binsizes.
        binsXB = xBinB{bsi}; % Get bins for this binsize
        binsXA = xBinA{bsi}; % Get bins for this binsize
        btsB = nans(length(gti),size(binsXB,1));
        btsA = nans(length(gti),size(binsXA,1));
        nts = [];
        nbef = []; naft = [];
        ngt = length(gti);
        % rand_beep_times = 1000*rand(1,100000);
        for ti = 1:ngt
            bTime = fd(ti) - 1100; % Create an arbitrary "beep" time.
            ts = [];
            ts = ats{ti}-bTime; % Spike times were recoded wrt mnky fixation start.
            % ts_fix = ts(ts>-1000 & ts<= 1000); % Get spikes from the 1sec after beep.
            ts_fixB = ts(ts>=-tt_time & ts<0); % Get spikes from the 1sec after beep.
            ts_fixA = ts(ts>0 & ts<=tt_time); % Get spikes from the 1sec after beep.
            ntsB(ti) = length(ts_fixB);
            ntsA(ti) = length(ts_fixA);
            % bsTemp = hist(ts_fix,binsX); % Spikes binned at 1msec.
            % wwi = find(bsTemp>0);
            % bsTemp(wwi)=1;
            % bs(ti,:) = bsTemp;
            for bi = 1:nbB{bsi} % Iterate over bins.
                % if binsXB(bi,1) >= tst(ti) &  binsXB(bi,2) <= tet(ti)
                btsB(ti,bi) = sum((ts_fixB > binsXB(bi,1) & ts_fixB <= binsXB(bi,2)));
                % end
            end
            for bi = 1:nbA{bsi} % Iterate over bins.
                % if binsXA(bi,1) >= tst(ti) &  binsXA(bi,2) <= tet(ti)
                btsA(ti,bi) = sum((ts_fixA > binsXA(bi,1) & ts_fixA <= binsXA(bi,2)));
                % end
            end % Iterate over bins.
        end % Trials loop.
        nA_binnedB{ai,bsi} = btsB;
        nA_binnedA{ai,bsi} = btsA;
        nA_rate(ai) = (1/winSiz)*nanmedian(ntsA+ntsB); % 042821.
        
    end % Counting bins loop.
    %     spCountBef{ai} = nbef;
    %     spCountAft{ai} = naft;
    
    % nA{ai} = bs; % Binned ACC spikes for this unit - ie - binary spike trains.
end

%% LC: Define FAKE phasic/no-phasic trials.

nL = cell(1,n_LC);
phasicU = cell(1,n_LC);
phasicU2 = cell(1,n_LC);
phasicU3 = cell(1,n_LC);
phasicU4 = cell(1,n_LC);
nts_u = [];

for nu = 1:n_LC
    prebeepspikes = nans(1,ngt);
    phasicYesNo = nans(1,ngt);
    phasicYesNo2 = nans(1,ngt);
    phasicYesNo3 = nans(1,ngt);
    phasicYesNo4 = nans(1,ngt);
    gzi = []; nts = [];
    ats = sps_LC(gti,nu);
    bs = nans(ngt,length(binsX)); % Preallocate.
    for ti = 1:ngt
        bTime = fd(ti) - 1100;
        ts = ats{ti} - bTime; % Spike times were recoded wrt mnky fixation start.
        ts_fix = ts(ts>-1000 & ts < 1000);
        
        ts_fix_postBeep = (1000/200)*length(ts(ts>200 & ts<=400)); % Get spikes from the PBI epoch.
        ts_fix_postBeep_all = (1000/1000)*length(ts(ts>0 & ts<=1000)); % Get spikes from PB epoch.
        ts_fix_Beep = (1000/200)*length(ts(ts> 0 & ts <=200)); % Get spikes from the phasic epoch.
        ts_fix_preBeep = (1000/1000)*length(ts(ts>-1000 & ts < 0)); % Get spikes from the pre-beep epoch.
        
        % Classify beep responses (or lack of response).
        if ts_fix_preBeep<ts_fix_Beep && ts_fix_postBeep<=ts_fix_preBeep % Phasic response and post-beep inhibition.
            phasicYesNo(ti) = 1;
            
        end
        if 1*ts_fix_preBeep<ts_fix_Beep % Phasic response but NO post-beep inhibition.
            phasicYesNo2(ti) = 1;
            
        end
        if ts_fix_preBeep>=ts_fix_Beep % NO Phasic response.
            phasicYesNo3(ti) = 1;
            
        end
        if ts_fix_preBeep>=2*ts_fix_Beep % NO Phasic response.
            phasicYesNo4(ti) = 1;
            
        end
        
        bs(ti,:) = hist(ts_fix,binsX); % Spikes binned at 1msec.
        
        % nts(ti) = length(ts_fix);
    end % Trials loop.
    
    phasicU{nu} = phasicYesNo;
    phasicU2{nu} = phasicYesNo2;
    phasicU3{nu} = phasicYesNo3;
    phasicU4{nu} = phasicYesNo4;
    
    
    nL{nu} = bs; % Binned LC spikes for this unit - ie - binary spike trains.
    
end

%% Data check, debug.

% pause

% Check to see if "phasic" responses are being identified.
% figure; hold on;
% for uii = 1:size(nL{1},1)
%
%     trial_sp = nL{1}(uii,:);
%     plot(binsX,trial_sp);
%     title(num2str(phasicU{1}(uii)));
%     %     xlim([-1000 1000]);
%     pause
%     clf
% end
%
% spRate = nanmean(nL{1});
% plot(spRate,'k.-');
%
% pause

% *********************************************
% % Debugging only: Look at the data, in case wonky...
%
% tz = [];
% for ai = 1:n_ACC
%     ty = nA{ai};
%     tz = [tz smooth(nanmean(ty),10)];
% end
%
% tz = tz';
% mz = nanmean(tz);
% plot(1000*smooth(mz,10));
%
% % Looked'ed at the data.
% *********************************************

%% Calculate ACC pairwise spike-count correlations. Case: Phasic with post-beep inhibition.

bsi = size(nbB,2);
udat = cell(1,n_LC);
if ~isempty(phasicU)
    for li = 1:n_LC
        accept_ai = 0; out_rsc_binned = cell(1,nbs);
        if ~isempty(phasicU{li})
            lhi = phasicU{li};
            % nts =nans(1,length(gti)); % Preallocate.
            lo_sp_ind = find(isnan(lhi));
            hi_sp_ind = find(lhi==1);
            % if length(lo_sp_ind)>=3 && length(hi_sp_ind)>=3
            %             if length(hi_sp_ind)>=3
            if length(hi_sp_ind)>=numTrialsLim
                for nbi = 1:bsi % Iterate over number of binsizes used (9) - drop the 1msec bin size!
                    out_rscB_NP = nans(n_ap,5); out_rscA_NP = nans(n_ap,5);
                    out_rscB_P = nans(n_ap,5); out_rscA_P = nans(n_ap,5);
                    out_rscB_All = nans(n_ap,5); out_rscAll_P = nans(n_ap,5);
                    mscB = nans(n_ACC,1); mscA = nans(n_ACC,1);
                    vscB = nans(n_ACC,1); vscA = nans(n_ACC,1);
                    % pB_a = nans(n_ap,nns); pA_a = nans(n_ap,nns); p_allB_a = nans(n_ap,nns); p_allA_a = nans(n_ap,nns);
                    
                    for jji = 1:n_ACC
                        spCountB = nA_binnedB{jji,nbi}(hi_sp_ind,:); spCountB = spCountB(:); mscB(jji) = nanmean(spCountB); vscB(jji) = nanvar(spCountB);
                        spCountA = nA_binnedA{jji,nbi}(hi_sp_ind,:); spCountA = spCountA(:); mscA(jji) = nanmean(spCountA); vscA(jji) = nanvar(spCountA);
                    end
                    % For this LC unit, iterate over all ACC pairs:
                    for ai = 1:n_ap
                        a1 = find(ACC_pairs(ai,1) == ACC_ch);
                        a2 = find(ACC_pairs(ai,2) == ACC_ch);
                        
                        if (nA_rate(a1) >= 1 & nA_rate(a2) >= 1) % 042821
                            
                            % FAKE Beep trials with "phasic" response.
                            
                            apcount1 = nA_binnedB{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedB{a2,nbi}(hi_sp_ind,:);
                            % rsc_result = get_pwc(reshape(apcount1.',1,[]),reshape(apcount2.',1,[]),numTrialsLim);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),numTrialsLim);
                            out_rscB_P(ai,:) = [rsc_result.RSC rsc_result.pval fi li ai];
                            
                            apcount1 = nA_binnedA{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedA{a2,nbi}(hi_sp_ind,:);
                            % rsc_result = get_pwc(reshape(apcount1.',1,[]),reshape(apcount2.',1,[]),numTrialsLim);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),numTrialsLim);
                            out_rscA_P(ai,:) = [rsc_result.RSC rsc_result.pval fi li ai];
                            
                            % ALL FAKE beep trials.
                            apcount1 = nA_binnedB{a1,nbi}; apcount2 = nA_binnedB{a2,nbi};
                            % rsc_result = get_pwc(reshape(apcount1.',1,[]),reshape(apcount2.',1,[]),numTrialsLim);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),numTrialsLim);
                            out_rscB_All(ai,:) = [rsc_result.RSC rsc_result.pval fi li ai];
                            
                            apcount1 = nA_binnedA{a1,nbi}; apcount2 = nA_binnedA{a2,nbi};
                            % rsc_result = get_pwc(reshape(apcount1.',1,[]),reshape(apcount2.',1,[]),numTrialsLim);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),numTrialsLim);
                            out_rscA_All(ai,:) = [rsc_result.RSC rsc_result.pval fi li ai];
                            
                        end % 042821
                        
                        if (nA_rate(a1) < 1 | nA_rate(a2) < 1) % 042821
                            
                            out_rscB_P(ai,:) = [nan nan fi li ai];
                            out_rscA_P(ai,:) = [nan nan fi li ai];
                            out_rscB_All(ai,:) = [nan nan fi li ai];
                            out_rscA_All(ai,:) = [nan nan fi li ai];
                            
                        end % 042821
                        
                    end % Pairs.
                    
                    % out_rsc_binned{nbi} = {out_rscB_NP out_rscA_NP out_rscB_P out_rscA_P out_rscB_All out_rscA_All};
                    
                    out_rsc_binned{nbi} = {out_rscB_P out_rscA_P out_rscB_All out_rscA_All  mscB vscB mscA vscA};
                end % Bins.
                udat{li} = out_rsc_binned;
            end % Check for enough trials.
            % Return empty if ANY of the LC spike conditions are not satisfied.
            % if length(lo_sp_ind)<3 && length(hi_sp_ind)<3
            if length(hi_sp_ind)<numTrialsLim
                udat{li} = [];
            end
        end
    end
    
    if isempty(phasicU{li})
        udat{li} = [];
    end
    
end % All LC units for this session.

%% Calculate ACC pairwise spike-count correlations. Case: Phasic with NO post-beep inhibition.

udat2 = cell(1,n_LC);

%% Calculate ACC pairwise spike-count correlations. Case: Phasic with NO post-beep inhibition.

bsi = size(nbB,2);
udat3 = cell(1,n_LC);
if ~isempty(phasicU3)
    for li = 1:n_LC
        accept_ai = 0; out_rsc_binned = cell(1,nbs);
        if ~isempty(phasicU3{li})
            lhi = phasicU3{li};
            % nts =nans(1,length(gti)); % Preallocate.
            lo_sp_ind = find(isnan(lhi));
            hi_sp_ind = find(lhi==1);
            % if length(lo_sp_ind)>=3 && length(hi_sp_ind)>=3
            %             if length(hi_sp_ind)>=3
            if length(hi_sp_ind)>=numTrialsLim
                for nbi = 1:bsi % Iterate over number of binsizes used (9) - drop the 1msec bin size!
                    out_rscB_NP = nans(n_ap,5); out_rscA_NP = nans(n_ap,5);
                    out_rscB_P = nans(n_ap,5); out_rscA_P = nans(n_ap,5);
                    out_rscB_All = nans(n_ap,5); out_rscA_All = nans(n_ap,5);
                    mscB = nans(n_ACC,1); mscA = nans(n_ACC,1);
                    vscB = nans(n_ACC,1); vscA = nans(n_ACC,1);
                    
                    for jji = 1:n_ACC
                        spCountB = nA_binnedB{jji,nbi}(hi_sp_ind,:); spCountB = spCountB(:); mscB(jji) = nanmean(spCountB); vscB(jji) = nanvar(spCountB);
                        spCountA = nA_binnedA{jji,nbi}(hi_sp_ind,:); spCountA = spCountA(:); mscA(jji) = nanmean(spCountA); vscA(jji) = nanvar(spCountA);
                    end
                    
                    % For this LC unit, iterate over all ACC pairs:
                    for ai = 1:n_ap
                        a1 = find(ACC_pairs(ai,1) == ACC_ch);
                        a2 = find(ACC_pairs(ai,2) == ACC_ch);
                        
                        if (nA_rate(a1) >= 1 & nA_rate(a2) >= 1) % 042821
                            
                            % FAKE Beep trials with "no phasic" response.
                            
                            apcount1 = nA_binnedB{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedB{a2,nbi}(hi_sp_ind,:);
                            % rsc_result = get_pwc(reshape(apcount1.',1,[]),reshape(apcount2.',1,[]),numTrialsLim);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),numTrialsLim);
                            out_rscB_P(ai,:) = [rsc_result.RSC rsc_result.pval fi li ai];
                            
                            apcount1 = nA_binnedA{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedA{a2,nbi}(hi_sp_ind,:);
                            % rsc_result = get_pwc(reshape(apcount1.',1,[]),reshape(apcount2.',1,[]),numTrialsLim);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),numTrialsLim);
                            out_rscA_P(ai,:) = [rsc_result.RSC rsc_result.pval fi li ai];
                            
                            % ALL beep trials.
                            apcount1 = nA_binnedB{a1,nbi}; apcount2 = nA_binnedB{a2,nbi};
                            % rsc_result = get_pwc(reshape(apcount1.',1,[]),reshape(apcount2.',1,[]),numTrialsLim);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),numTrialsLim);
                            out_rscB_All(ai,:) = [rsc_result.RSC rsc_result.pval fi li ai];
                            
                            apcount1 = nA_binnedA{a1,nbi}; apcount2 = nA_binnedA{a2,nbi};
                            % rsc_result = get_pwc(reshape(apcount1.',1,[]),reshape(apcount2.',1,[]),numTrialsLim);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),numTrialsLim);
                            out_rscA_All(ai,:) = [rsc_result.RSC rsc_result.pval fi li ai];
                            
                        end
                        
                        if (nA_rate(a1) < 1 | nA_rate(a2) < 1) % 042821
                            
                            out_rscB_P(ai,:) = [nan nan fi li ai];
                            out_rscA_P(ai,:) = [nan nan fi li ai];
                            out_rscB_All(ai,:) = [nan nan fi li ai];
                            out_rscA_All(ai,:) = [nan nan fi li ai];
                            
                        end % 042821
                        
                    end % Pairs.
                    
                    % out_rsc_binned{nbi} = {out_rscB_NP out_rscA_NP out_rscB_P out_rscA_P};
                    out_rsc_binned{nbi} = {out_rscB_P out_rscA_P out_rscB_All out_rscA_All mscB vscB mscA vscA};
                end % Bins.
                udat3{li} = out_rsc_binned;
            end % Check for enough trials.
            % Return empty if ANY of the LC spike conditions are not satisfied.
            % if length(lo_sp_ind)<3 && length(hi_sp_ind)<3
            % if length(hi_sp_ind)<numTrialsLim
            if length(hi_sp_ind)<numTrialsLim
                udat3{li} = [];
            end
        end
    end
    
    if isempty(phasicU3{li})
        udat3{li} = [];
    end
    
end % All LC units for this session.

%% Calculate ACC pairwise spike-count correlations. Case: Phasic with NO post-beep inhibition.

udat4 = cell(1,n_LC);

%%

% Collect results for saving:
beep_rsc_ = {udat nL phasicU phasicU2 phasicU3 phasicU4 udat2 udat3 udat4};

% ALL DONE!

%%

% figure;
% cg = udat{1}{1}-udat{1}{2}; plot(smooth(nanmean(cg),100),'k'); hold on; pause;
% cg = udat{1}{3}-udat{1}{4}; plot(smooth(nanmean(cg),100),'r'); pause; close all;


%% Troubleshooting and data check figures:

% ********************************************************

% % Debugging only:
% %
% % Plots to check analysis thus far...
% %
% % sp1 = al_hi1(n_hi,:);
% % sp2 = al_hi2(n_hi,:);
% smW = 100;
% % Preallocate.
% smLo = nans(size(all_Lo_ccg));
% smHi = nans(size(all_Hi_ccg));
% smLo_s = nans(size(all_Lo_ccg));
% smHi_s = nans(size(all_Hi_ccg));
% smLo_su = nans(size(all_Lo_ccg));
% smHi_su = nans(size(all_Hi_ccg));
% smLoSCu = nans(size(all_Lo_ccg));
% smHiSCu = nans(size(all_Hi_ccg));
% smLou = nans(size(all_Lo_ccg));
% smHiu = nans(size(all_Hi_ccg));
% smLoSC = nans(size(all_Lo_ccg));
% smHiSC = nans(size(all_Hi_ccg));
%
% lcount = 0; hcount = 0;
% for ii = 1:size(all_Lo_ccg,1)
%     if ~isnan(sum(all_Lo_ccg(ii,:))) & ~isnan(sum(all_Lo_ccg_s(ii,:)))
%         lcount = lcount+1;
%         smLo(ii,:) = smooth(all_Lo_ccg(ii,:),smW);
%         smLo_s(ii,:) = smooth(all_Lo_ccg_s(ii,:),smW);
%         smLou(ii,:) = smooth(all_Lo_ccg_u(ii,:),smW);
%         smLoSC(ii,:) = smooth(all_Lo_ccg(ii,:)-all_Lo_ccg_s(ii,:),smW);
%         smLo_su(ii,:) = smooth(all_Lo_ccg_s_u(ii,:),smW);
%         smLoSCu(ii,:) = smooth(all_Lo_ccg_u(ii,:)-all_Lo_ccg_s_u(ii,:),smW);
%     end
%     if ~isnan(sum(all_Hi_ccg(ii,:))) & ~isnan(sum(all_Hi_ccg_s(ii,:)))
%         hcount = hcount+1;
%         smHi(ii,:) = smooth(all_Hi_ccg(ii,:),smW);
%         smHi_s(ii,:) = smooth(all_Hi_ccg_s(ii,:),smW);
%         smHiu(ii,:) = smooth(all_Hi_ccg_u(ii,:),smW);
%         smHiSC(ii,:) = smooth(all_Hi_ccg(ii,:)-all_Hi_ccg_s(ii,:),smW);
%         smHi_su(ii,:) = smooth(all_Hi_ccg_s_u(ii,:),smW);
%         smHiSCu(ii,:) = smooth(all_Hi_ccg_u(ii,:)-all_Hi_ccg_s_u(ii,:),smW);
%     end
% end
%
% figure; hold on;
% subplot(3,1,1); hold on; plot(nanmean(smLo),'k'); plot(nanmean(smHi),'r'); ylim([-0.1 0.15]);
% subplot(3,1,2); hold on; plot(nanmean(smLo_s),'k'); plot(nanmean(smHi_s),'r'); ylim([-0.1 0.15]);
% subplot(3,1,3); hold on; plot(nanmean(smLoSC),'k'); plot(nanmean(smHiSC),'r'); ylim([-0.1 0.15]);
%
% figure; hold on;
% subplot(3,1,1); hold on; plot(nanmean(smLou),'k'); plot(nanmean(smHiu),'r'); ylim([0 0.0001]);
% subplot(3,1,2); hold on; plot(nanmean(smLo_su),'k'); plot(nanmean(smHi_su),'r'); ylim([0 0.0001]);
% subplot(3,1,3); hold on; plot(nanmean(smLoSCu),'k'); plot(nanmean(smHiSCu),'r'); ylim([0 0.0001]);

% % ********************************************************


