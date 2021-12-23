function [beep_rsc_] = getLC_DualArea_Beep_binned_pupil_rsc3(siteData,nlim,fi)

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
%% Find beep trials ...
LBeep    = ~isnan(siteData{1}(:,4)); % Only beep trials
beepTime = siteData{1}(LBeep,4);
FBeep    = find(LBeep);

%% Get pupil data ...

% default -- pupil diameter
pdat = siteData{2}(FBeep,:,4); % trial/sample
% pupil slope
% pdatS = siteData{2}(FBeep,:,5); % trial/sample

% num_samps = size(pdat,2);

%% Set up bins ...

% Create binsizes.
tt_time = 1000;
% bs = linspace(100,tt_time,10); % Before 043021.
bs = linspace(200,tt_time,5); % 043021.

% backTime = 1001; fwdTime = 2100; % We're looking at the 1sec to 2.1sec window of stable fixation.
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
fd = siteData{1}(LBeep,2); % Fixation duration.

%% Get ACC binned spike counts.

sps_LC = siteData{3}(LBeep,LC_ch); % LC spikes.
sps_ACC = siteData{3}(LBeep,ACC_ch); % ACC spikes.

% Want 1 sec of ACC spikes after beep.
% gti = find(beepTime>500 & fd-beepTime>=tt_time+100);
% gti = find(beepTime>1000 & fd-beepTime>=tt_time+100);
gti = find(beepTime>1000 & fd-beepTime>=tt_time+1);
ngt = length(gti);
fd = fd(gti);
beepTime = beepTime(gti);
ccgbSize = 1;
binsX = -tt_time:ccgbSize:tt_time; % Bins for creating binary spike train for 1st 1 sec of stable fixation after beep.

%% Get mean pupil size for good no beep trials.

% Test plot.
% figure; hold on;
beep_pupil = nans(ngt,1500);
t_pdat = pdat(gti,:); % All good beep trials.
% t_pdatS = pdatS(gti,:); % All good beep trials.

% tax = -500:1000;
peak_pupil = nans(1,ngt);
all_beep_pupil = nans(1,ngt);

for ti=1:ngt
    
    bTime = ceil(beepTime(ti)); % Get beep time.
    
    %         baseline_pupil(ti) = nanmean(t_pdat(ti,bTime-150:bTime-75));
    %         peak_pupil(ti) = max(t_pdat(ti,bTime+100:bTime+1000)) - nanmean(t_pdat(ti,bTime-150:bTime-75));
    %         evoked_pupil = max(t_pdat(ti,bTime+100:bTime+1000)) - nanmean(t_pdat(ti,bTime-150:bTime-75));
    
    % peak_pupil(ti) = max(t_pdat(ti,bTime+100:bTime+1000)) - nanmean(t_pdat(ti,bTime-100:bTime));
    baseline_pupil(ti) = nanmean(t_pdat(ti,bTime-200:bTime));
    evoked_pupil = max(t_pdat(ti,bTime+100:bTime+1000)) - nanmean(t_pdat(ti,bTime-200:bTime));
    
    all_beep_pupil(ti) = evoked_pupil;
    %     peak_pupil(ti) = evoked_pupil;
    
    if evoked_pupil > 0
        peak_pupil(ti) = evoked_pupil;
    end
    
    if evoked_pupil <= 0
        peak_pupil(ti) = nan;
    end
    
    %     plot(tax,t_pdat(ti,bTime-500:bTime+1000),'k-');
    %     %     hold on;
    %     plot([0 0],[-2 6],'k--');
    %     piii(ti) = 700+find(t_pdat(ti,bTime+200:bTime+1000) == peak_pupil(ti));
    %     plot(tax(piii(ti)),peak_pupil(ti),'r.','markersize',12);
    %     %     pause
    %     %     clf
    
end

% Do regression on max pupil vs baseline pupil and get residuals.
xd1 = [ones(ngt,1) baseline_pupil'];
yd1 = peak_pupil';
[b,bint,r_pupil,rint,stats]=regress(yd1,xd1); % r_pupil is the residual vector.

% % Test plots.
% figure; hold on;
% plot(baseline_pupil,peak_pupil,'k.');
% plot(baseline_pupil,r_pupil,'r.');
% x_plot = linspace(min(baseline_pupil),max(baseline_pupil),10);
% y_plot = b(1) + b(2)*x_plot;
% plot(x_plot,y_plot,'g-'); % Plot linear regression line.
%
% figure; hold on;
% plot(peak_pupil,r_pupil,'k.')
% plot(peak_pupil-baseline_pupil,r_pupil,'ro')
% plot(peak_pupil-baseline_pupil,peak_pupil,'g^')

%%

% Collect spikes from window of size winSiz - needed for calculating mean spike rate.
% winSiz = length(binsX)/1000;
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
        for ti = 1:ngt
            bTime = beepTime(ti);
            ts = [];
            ts = ats{ti}-bTime; % Spike times were recoded wrt mnky fixation start.
            ts_fixB = ts(ts>=-tt_time & ts<0); % Get spikes from the 1sec after beep.
            ts_fixA = ts(ts>0 & ts<=tt_time); % Get spikes from the 1sec after beep.
            ntsA(ti) = length(ts_fixB);
            ntsB(ti) = length(ts_fixA);
            for bi = 1:nbB{bsi} % Iterate over bins.
                btsB(ti,bi) = sum((ts_fixB > binsXB(bi,1) & ts_fixB < binsXB(bi,2)));
                % end
            end
            for bi = 1:nbA{bsi} % Iterate over bins.
                btsA(ti,bi) = sum((ts_fixA > binsXA(bi,1) & ts_fixA < binsXA(bi,2)));
            end % Iterate over bins.
        end % Trials loop.
        nA_binnedB{ai,bsi} = btsB;
        nA_binnedA{ai,bsi} = btsA;
        nA_rate(ai) = (1/winSiz)*nanmedian(ntsA+ntsB); % 043021.
    end % Counting bins loop.
end

%% LC: Define phasic/no-phasic trials.

nL = cell(1,n_LC); phasicU = cell(1,n_LC); phasicU2 = cell(1,n_LC); phasicU3 = cell(1,n_LC); phasicU4 = cell(1,n_LC); phasicRespU = cell(1,n_LC); nts_u = [];
for nu = 1:n_LC
    phasicYesNo = nans(1,length(gti)); phasicYesNo2 = nans(1,length(gti)); phasicYesNo3 = nans(1,length(gti)); phasicYesNo4 = nans(1,length(gti)); phasicResp = nans(1,length(gti));
    gzi = []; nts = [];
    ats = sps_LC(gti,nu);
    bs = nans(length(gti),length(binsX)); % Preallocate.
    for ti = 1:ngt
        bTime = beepTime(ti);
        ts = ats{ti} - bTime; % Spike times were recoded wrt mnky fixation start.
        ts_fix = ts(ts > -1000 & ts < 1000);
        
        ts_fix_postBeep = (1000/200)*length(ts(ts>200 & ts<=400)); % Get spikes from the PBI epoch.
        ts_fix_postBeep_all = (1000/1000)*length(ts(ts>0 & ts<=1000)); % Get spikes from PB epoch.
        ts_fix_Beep = (1000/200)*length(ts(ts> 0 & ts <=200)); % Get spikes from the phasic epoch.
        ts_fix_preBeep = (1000/1000)*length(ts(ts>-1000 & ts <=0)); % Get spikes from the pre-beep epoch.
        
        % Classify beep responses (or lack of response).
        if ts_fix_preBeep<ts_fix_Beep && ts_fix_postBeep<=ts_fix_preBeep; % Phasic response and post-beep inhibition. % 041320.
            phasicYesNo(ti) = 1;
        end
        if ts_fix_preBeep<ts_fix_Beep % Phasic response but NO post-beep inhibition.
            phasicYesNo2(ti) = 1;
        end
        if ts_fix_preBeep>=ts_fix_Beep % NO Phasic response.
            phasicYesNo3(ti) = 1;
        end
        if ts_fix_preBeep>=2*ts_fix_Beep % NO Phasic response.
            phasicYesNo4(ti) = 1;
        end
        phasicResp(ti) = ts_fix_Beep;
        bs(ti,:) = hist(ts_fix,binsX); % Spikes binned at 1msec.
    end % Trials loop.
    
    phasicU{nu} = phasicYesNo; phasicU2{nu} = phasicYesNo2; phasicU3{nu} = phasicYesNo3; phasicU4{nu} = phasicYesNo4;
    phasicRespU{nu} = phasicResp;
    nL{nu} = bs; % Binned LC spikes for this unit - ie - binary spike trains.
    
end % LC units loop.

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

%% Calculate ACC pairwise spike-count correlations. Case: LC Phasic with post-beep inhibition.

bsi = size(nbB,2);
pp_lo_phasic = []; pp_hi_phasic = [];

udat = cell(1,n_LC);
if ~isempty(phasicU)
    for li = 1:n_LC
        accept_ai = 0; out_rsc_binned = cell(1,nbs);
        if ~isempty(phasicU{li})
            lhi = phasicU{li};
            asi = 1:length(gti);
            nts =nans(1,length(gti)); % Preallocate.
            
            %             lo_sp_ind = find(isnan(lhi));
            hi_sp_ind = find(lhi==1 & isfinite(peak_pupil));
            % hi_sp_ind = find(lhi==1);
            
            % evoked_pupil = peak_pupil;
            % evoked_pupil = r_pupil;
            
            r_pupil_hi = peak_pupil(hi_sp_ind); % ONLY phasic LC trials.
            % r_pupil_hi = evoked_pupil; hi_sp_ind = 1:length(r_pupil_hi); % ALL trials.
            
            r_vals = prctile(r_pupil_hi,[25 50 75]);
            
            %             % Quartile split.
            %             r_lo_i = find(r_pupil_hi<=r_vals(1));
            %             r_hi_i = find(r_pupil_hi>r_vals(3));
            
            % Median split.
            r_lo_i = find(r_pupil_hi<r_vals(2));
            r_hi_i = find(r_pupil_hi>=r_vals(2));
            
            pp_lo_phasic = peak_pupil(r_lo_i);
            pp_hi_phasic = peak_pupil(r_hi_i);
            
            if ~isempty(pp_lo_phasic) & ~isempty(pp_hi_phasic)
                
                if length(r_lo_i)>=nlim & length(r_hi_i)>=nlim
                    for nbi = 1:bsi % Iterate over number of binsizes used (9) - drop the 1msec bin size!
                        out_rscB_LO = nans(n_ap,2); out_rscA_LO = nans(n_ap,2); out_rscB_HI = nans(n_ap,2); out_rscA_HI = nans(n_ap,2);
                        
                        for jji = 1:n_ACC
                            spCountB = nA_binnedB{jji,nbi}(hi_sp_ind,:); spCountB = spCountB(:); mscB(jji) = nanmean(spCountB); vscB(jji) = nanvar(spCountB);
                            spCountA = nA_binnedA{jji,nbi}(hi_sp_ind,:); spCountA = spCountA(:); mscA(jji) = nanmean(spCountA); vscA(jji) = nanvar(spCountA);
                        end
                        % For this LC unit, iterate over all ACC pairs:
                        for ai = 1:n_ap
                            a1 = find(ACC_pairs(ai,1) == ACC_ch);
                            a2 = find(ACC_pairs(ai,2) == ACC_ch);
                            
                            % Spike count correlation.
                            % if (nA_rate(a1) >= 1 & nA_rate(a2) >= 1) % 042821
                            
                            % Phasic response Beep trials, pupil based median split - LOW evoked pupil.
                            
                            apcount1 = nA_binnedB{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedB{a2,nbi}(hi_sp_ind,:);
                            apcount1 = apcount1(r_lo_i,:); apcount2 = apcount2(r_lo_i,:);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                            out_rscB_LO(ai,:) = [rsc_result.RSC rsc_result.pval];
                            
                            apcount1 = nA_binnedA{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedA{a2,nbi}(hi_sp_ind,:);
                            apcount1 = apcount1(r_lo_i,:); apcount2 = apcount2(r_lo_i,:);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                            out_rscA_LO(ai,:) = [rsc_result.RSC rsc_result.pval];
                            
                            % Phasic response Beep trials, pupil based median split - HIGH evoked pupil.
                            apcount1 = nA_binnedB{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedB{a2,nbi}(hi_sp_ind,:);
                            apcount1 = apcount1(r_hi_i,:); apcount2 = apcount2(r_hi_i,:);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                            out_rscB_HI(ai,:) = [rsc_result.RSC rsc_result.pval];
                            
                            apcount1 = nA_binnedA{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedA{a2,nbi}(hi_sp_ind,:);
                            apcount1 = apcount1(r_hi_i,:); apcount2 = apcount2(r_hi_i,:);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                            out_rscA_HI(ai,:) = [rsc_result.RSC rsc_result.pval];
                            
                            % end
                            %                             if (nA_rate(a1) < 1 | nA_rate(a2) < 1) % 042821
                            %                                 out_rscB_LO(ai,:) = [nan nan]; out_rscA_LO(ai,:) = [nan nan]; out_rscB_HI(ai,:) = [nan nan]; out_rscA_HI(ai,:) = [nan nan];
                            %                             end % 042821
                        end % Pairs.
                        %                         out_rsc_binned{nbi} = {out_rscB_P out_rscA_P out_rscB_All out_rscA_All mscB vscB mscA vscA out_rscB_LO out_rscA_LO out_rscB_HI out_rscA_HI pp_lo_phasic(:) pp_hi_phasic(:)};
                        out_rsc_binned{nbi} = {out_rscB_LO out_rscA_LO out_rscB_HI out_rscA_HI pp_lo_phasic(:) pp_hi_phasic(:)};
                    end % Bins.
                    udat{li} = out_rsc_binned;
                end % Check for enough trials.
            end % Check that there are lo and hi pupil trials.
            % Return empty if ANY of the LC spike conditions are not satisfied.
            if length(hi_sp_ind)<nlim
                udat{li} = [];
            end
            if isempty(pp_lo_phasic) | isempty(pp_hi_phasic)
                udat{li} = [];
            end
        end
    end
    if isempty(phasicU{li})
        udat{li} = [];
    end
end % All LC units for this session.

%% NOW - calculate ACC pairwise spike-count correlations. Case: LC No Phasic.

udat3 = cell(1,n_LC);

bsi = size(nbB,2);

if ~isempty(phasicU3)
    for li = 1:n_LC
        accept_ai = 0; out_rsc_binned = cell(1,nbs);
        if ~isempty(phasicU3{li})
            lhi = phasicU3{li};
            nts =nans(1,length(gti)); % Preallocate.
            
            %             lo_sp_ind = find(isnan(lhi));
            hi_sp_ind = find(lhi==1 & isfinite(peak_pupil));
            % hi_sp_ind = find(lhi==1);
            
            % evoked_pupil = peak_pupil;
            %             %             evoked_pupil = r_pupil;
            
            r_pupil_hi = peak_pupil(hi_sp_ind); % ONLY phasic LC trials.
            % r_pupil_hi = evoked_pupil; hi_sp_ind = 1:length(r_pupil_hi); % ALL trials.
            
            r_vals = prctile(r_pupil_hi,[25 50 75]);
            
            %             % Quartile split.
            %             r_lo_i = find(r_pupil_hi<=r_vals(1));
            %             r_hi_i = find(r_pupil_hi>r_vals(3));
            
            % Median split.
            r_lo_i = find(r_pupil_hi<r_vals(2));
            r_hi_i = find(r_pupil_hi>=r_vals(2));
            
            pp_lo_nophasic = peak_pupil(r_lo_i);
            pp_hi_nophasic = peak_pupil(r_hi_i);
            
            if length(r_lo_i)>=nlim & length(r_hi_i)>=nlim
                for nbi = 1:bsi % Iterate over number of binsizes used (9) - drop the 1msec bin size!
                    
                    out_rscB_LO = nans(n_ap,2); out_rscA_LO = nans(n_ap,2); out_rscB_HI = nans(n_ap,2); out_rscA_HI = nans(n_ap,2);
                    
                    for jji = 1:n_ACC
                        spCountB = nA_binnedB{jji,nbi}(hi_sp_ind,:); spCountB = spCountB(:); mscB(jji) = nanmean(spCountB); vscB(jji) = nanvar(spCountB);
                        spCountA = nA_binnedA{jji,nbi}(hi_sp_ind,:); spCountA = spCountA(:); mscA(jji) = nanmean(spCountA); vscA(jji) = nanvar(spCountA);
                    end
                    
                    % For this LC unit, iterate over all ACC pairs:
                    for ai = 1:n_ap
                        a1 = find(ACC_pairs(ai,1) == ACC_ch);
                        a2 = find(ACC_pairs(ai,2) == ACC_ch);
                        
                        % if (nA_rate(a1) >= 1 & nA_rate(a2) >= 1) % 042821
                        
                        % Phasic response Beep trials, pupil based median split - LOW evoked pupil.
                        
                        apcount1 = nA_binnedB{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedB{a2,nbi}(hi_sp_ind,:);
                        apcount1 = apcount1(r_lo_i,:); apcount2 = apcount2(r_lo_i,:);
                        rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                        out_rscB_LO(ai,:) = [rsc_result.RSC rsc_result.pval];
                        
                        apcount1 = nA_binnedA{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedA{a2,nbi}(hi_sp_ind,:);
                        apcount1 = apcount1(r_lo_i,:); apcount2 = apcount2(r_lo_i,:);
                        rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                        out_rscA_LO(ai,:) = [rsc_result.RSC rsc_result.pval];
                        
                        % Phasic response Beep trials, pupil based median split - HIGH evoked pupil.
                        apcount1 = nA_binnedB{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedB{a2,nbi}(hi_sp_ind,:);
                        apcount1 = apcount1(r_hi_i,:); apcount2 = apcount2(r_hi_i,:);
                        rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                        out_rscB_HI(ai,:) = [rsc_result.RSC rsc_result.pval];
                        
                        apcount1 = nA_binnedA{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedA{a2,nbi}(hi_sp_ind,:);
                        apcount1 = apcount1(r_hi_i,:); apcount2 = apcount2(r_hi_i,:);
                        rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                        out_rscA_HI(ai,:) = [rsc_result.RSC rsc_result.pval];
                        
                        % end
                        %                         if (nA_rate(a1) < 1 | nA_rate(a2) < 1) % 042821
                        %                             out_rscB_LO(ai,:) = [nan nan]; out_rscA_LO(ai,:) = [nan nan]; out_rscB_HI(ai,:) = [nan nan]; out_rscA_HI(ai,:) = [nan nan];
                        %                         end % 042821
                    end % Pairs.
                    % out_rsc_binned{nbi} = {out_rscB_NP out_rscA_NP out_rscB_P out_rscA_P};
                    %                     out_rsc_binned{nbi} = {out_rscB_P out_rscA_P out_rscB_All out_rscA_All mscB vscB mscA vscA out_rscB_LO out_rscA_LO out_rscB_HI out_rscA_HI pp_lo_nophasic(:) pp_hi_nophasic(:)};
                    out_rsc_binned{nbi} = {out_rscB_LO out_rscA_LO out_rscB_HI out_rscA_HI pp_lo_nophasic(:) pp_hi_nophasic(:)};
                end % Bins.
                udat3{li} = out_rsc_binned;
            end % Check for enough trials.
            % Return empty if ANY of the LC spike conditions are not satisfied.
            if length(hi_sp_ind)<nlim
                udat3{li} = [];
            end
        end
    end
    if isempty(phasicU3{li})
        udat3{li} = [];
    end
end % All LC units for this session.


% *****************************
% *****************************
% *****************************

%% Calculate ACC pairwise spike-count correlations. Case: Phasic with post-beep inhibition.

udatB = cell(1,n_LC);
if ~isempty(phasicU)
    for li = 1:n_LC
        accept_ai = 0; out_rsc_binned = cell(1,nbs);
        if ~isempty(phasicU{li})
            lhi = phasicU{li};
            asi = 1:length(gti);
            nts =nans(1,length(gti)); % Preallocate.
            
            evoked_pupil = baseline_pupil;
            % % evoked_pupil = r_pupil;
            
            % lo_sp_ind = find(isnan(lhi));
            hi_sp_ind = find(lhi==1 & isfinite(evoked_pupil));
            % hi_sp_ind = find(lhi==1);
            
            r_pupil_hi = evoked_pupil(hi_sp_ind); % ONLY phasic LC trials.
            % r_pupil_hi = evoked_pupil; hi_sp_ind = 1:length(r_pupil_hi); % ALL trials.
            
            r_vals = prctile(r_pupil_hi,[25 50 75]);
            
            %             % Quartile split.
            %             r_lo_i = find(r_pupil_hi<=r_vals(1));
            %             r_hi_i = find(r_pupil_hi>r_vals(3));
            
            % Median split.
            r_lo_i = find(r_pupil_hi<r_vals(2));
            r_hi_i = find(r_pupil_hi>=r_vals(2));
            
            pp_lo_phasic = evoked_pupil(r_lo_i);
            pp_hi_phasic = evoked_pupil(r_hi_i);
            
            if ~isempty(pp_lo_phasic) & ~isempty(pp_hi_phasic)
                
                if length(r_lo_i)>=nlim & length(r_hi_i)>=nlim
                    for nbi = 1:bsi % Iterate over number of binsizes used (9) - drop the 1msec bin size!
                        out_rscB_LO = nans(n_ap,2); out_rscA_LO = nans(n_ap,2); out_rscB_HI = nans(n_ap,2); out_rscA_HI = nans(n_ap,2);
                        
                        for jji = 1:n_ACC
                            spCountB = nA_binnedB{jji,nbi}(hi_sp_ind,:); spCountB = spCountB(:); mscB(jji) = nanmean(spCountB); vscB(jji) = nanvar(spCountB);
                            spCountA = nA_binnedA{jji,nbi}(hi_sp_ind,:); spCountA = spCountA(:); mscA(jji) = nanmean(spCountA); vscA(jji) = nanvar(spCountA);
                        end
                        % For this LC unit, iterate over all ACC pairs:
                        for ai = 1:n_ap
                            a1 = find(ACC_pairs(ai,1) == ACC_ch);
                            a2 = find(ACC_pairs(ai,2) == ACC_ch);
                            
                            % Spike count correlation.
                            % if (nA_rate(a1) >= 1 & nA_rate(a2) >= 1) % 042821
                            
                            % Phasic response Beep trials, pupil based median split - LOW evoked pupil.
                            
                            apcount1 = nA_binnedB{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedB{a2,nbi}(hi_sp_ind,:);
                            apcount1 = apcount1(r_lo_i,:); apcount2 = apcount2(r_lo_i,:);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                            out_rscB_LO(ai,:) = [rsc_result.RSC rsc_result.pval];
                            
                            apcount1 = nA_binnedA{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedA{a2,nbi}(hi_sp_ind,:);
                            apcount1 = apcount1(r_lo_i,:); apcount2 = apcount2(r_lo_i,:);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                            out_rscA_LO(ai,:) = [rsc_result.RSC rsc_result.pval];
                            
                            % Phasic response Beep trials, pupil based median split - HIGH evoked pupil.
                            apcount1 = nA_binnedB{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedB{a2,nbi}(hi_sp_ind,:);
                            apcount1 = apcount1(r_hi_i,:); apcount2 = apcount2(r_hi_i,:);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                            out_rscB_HI(ai,:) = [rsc_result.RSC rsc_result.pval];
                            
                            apcount1 = nA_binnedA{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedA{a2,nbi}(hi_sp_ind,:);
                            apcount1 = apcount1(r_hi_i,:); apcount2 = apcount2(r_hi_i,:);
                            rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                            out_rscA_HI(ai,:) = [rsc_result.RSC rsc_result.pval];
                            
                            % end
                            %                             if (nA_rate(a1) < 1 | nA_rate(a2) < 1) % 042821
                            %                                 out_rscB_LO(ai,:) = [nan nan]; out_rscA_LO(ai,:) = [nan nan]; out_rscB_HI(ai,:) = [nan nan]; out_rscA_HI(ai,:) = [nan nan];
                            %                             end % 042821
                        end % Pairs.
                        %                         out_rsc_binned{nbi} = {out_rscB_P out_rscA_P out_rscB_All out_rscA_All mscB vscB mscA vscA out_rscB_LO out_rscA_LO out_rscB_HI out_rscA_HI pp_lo_phasic(:) pp_hi_phasic(:)};
                        out_rsc_binned{nbi} = {out_rscB_LO out_rscA_LO out_rscB_HI out_rscA_HI pp_lo_phasic(:) pp_hi_phasic(:)};
                    end % Bins.
                    udatB{li} = out_rsc_binned;
                end % Check for enough trials.
            end % Check that there are lo and hi pupil trials.
            % Return empty if ANY of the LC spike conditions are not satisfied.
            if length(hi_sp_ind)<nlim
                udatB{li} = [];
            end
            if isempty(pp_lo_phasic) | isempty(pp_hi_phasic)
                udatB{li} = [];
            end
        end
    end
    if isempty(phasicU{li})
        udatB{li} = [];
    end
end % All LC units for this session.

%% NOW - calculate ACC pairwise spike-count correlations. Case: LC No Phasic.

udatB3 = cell(1,n_LC);

bsi = size(nbB,2);

if ~isempty(phasicU3)
    for li = 1:n_LC
        accept_ai = 0; out_rsc_binned = cell(1,nbs);
        if ~isempty(phasicU3{li})
            lhi = phasicU3{li};
            nts =nans(1,length(gti)); % Preallocate.
            
            evoked_pupil = baseline_pupil;
            % evoked_pupil = r_pupil;
            
            % lo_sp_ind = find(isnan(lhi));
            hi_sp_ind = find(lhi==1 & isfinite(evoked_pupil));
            % hi_sp_ind = find(lhi==1);
            
            r_pupil_hi = evoked_pupil(hi_sp_ind); % ONLY phasic LC trials.
            % r_pupil_hi = evoked_pupil; hi_sp_ind = 1:length(r_pupil_hi); % ALL trials.
            
            r_vals = prctile(r_pupil_hi,[25 50 75]);
            
            %             % Quartile split.
            %             r_lo_i = find(r_pupil_hi<=r_vals(1));
            %             r_hi_i = find(r_pupil_hi>r_vals(3));
            
            % Median split.
            r_lo_i = find(r_pupil_hi<r_vals(2));
            r_hi_i = find(r_pupil_hi>=r_vals(2));
            
            pp_lo_nophasic = evoked_pupil(r_lo_i);
            pp_hi_nophasic = evoked_pupil(r_hi_i);
            
            if length(r_lo_i)>=nlim & length(r_hi_i)>=nlim
                for nbi = 1:bsi % Iterate over number of binsizes used (9) - drop the 1msec bin size!
                    
                    out_rscB_LO = nans(n_ap,2); out_rscA_LO = nans(n_ap,2); out_rscB_HI = nans(n_ap,2); out_rscA_HI = nans(n_ap,2);
                    
                    for jji = 1:n_ACC
                        spCountB = nA_binnedB{jji,nbi}(hi_sp_ind,:); spCountB = spCountB(:); mscB(jji) = nanmean(spCountB); vscB(jji) = nanvar(spCountB);
                        spCountA = nA_binnedA{jji,nbi}(hi_sp_ind,:); spCountA = spCountA(:); mscA(jji) = nanmean(spCountA); vscA(jji) = nanvar(spCountA);
                    end
                    
                    % For this LC unit, iterate over all ACC pairs:
                    for ai = 1:n_ap
                        a1 = find(ACC_pairs(ai,1) == ACC_ch);
                        a2 = find(ACC_pairs(ai,2) == ACC_ch);
                        
                        % if (nA_rate(a1) >= 1 & nA_rate(a2) >= 1) % 042821
                        
                        % Phasic response Beep trials, pupil based median split - LOW evoked pupil.
                        
                        apcount1 = nA_binnedB{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedB{a2,nbi}(hi_sp_ind,:);
                        apcount1 = apcount1(r_lo_i,:); apcount2 = apcount2(r_lo_i,:);
                        rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                        out_rscB_LO(ai,:) = [rsc_result.RSC rsc_result.pval];
                        
                        apcount1 = nA_binnedA{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedA{a2,nbi}(hi_sp_ind,:);
                        apcount1 = apcount1(r_lo_i,:); apcount2 = apcount2(r_lo_i,:);
                        rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                        out_rscA_LO(ai,:) = [rsc_result.RSC rsc_result.pval];
                        
                        % Phasic response Beep trials, pupil based median split - HIGH evoked pupil.
                        apcount1 = nA_binnedB{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedB{a2,nbi}(hi_sp_ind,:);
                        apcount1 = apcount1(r_hi_i,:); apcount2 = apcount2(r_hi_i,:);
                        rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                        out_rscB_HI(ai,:) = [rsc_result.RSC rsc_result.pval];
                        
                        apcount1 = nA_binnedA{a1,nbi}(hi_sp_ind,:); apcount2 = nA_binnedA{a2,nbi}(hi_sp_ind,:);
                        apcount1 = apcount1(r_hi_i,:); apcount2 = apcount2(r_hi_i,:);
                        rsc_result = get_pwc(apcount1(:),apcount2(:),nlim);
                        out_rscA_HI(ai,:) = [rsc_result.RSC rsc_result.pval];
                        
                        % end
                        %                         if (nA_rate(a1) < 1 | nA_rate(a2) < 1) % 042821
                        %                             out_rscB_LO(ai,:) = [nan nan]; out_rscA_LO(ai,:) = [nan nan]; out_rscB_HI(ai,:) = [nan nan]; out_rscA_HI(ai,:) = [nan nan];
                        %                         end % 042821
                    end % Pairs.
                    % out_rsc_binned{nbi} = {out_rscB_NP out_rscA_NP out_rscB_P out_rscA_P};
                    %                     out_rsc_binned{nbi} = {out_rscB_P out_rscA_P out_rscB_All out_rscA_All mscB vscB mscA vscA out_rscB_LO out_rscA_LO out_rscB_HI out_rscA_HI pp_lo_nophasic(:) pp_hi_nophasic(:)};
                    out_rsc_binned{nbi} = {out_rscB_LO out_rscA_LO out_rscB_HI out_rscA_HI pp_lo_nophasic(:) pp_hi_nophasic(:)};
                end % Bins.
                udatB3{li} = out_rsc_binned;
            end % Check for enough trials.
            % Return empty if ANY of the LC spike conditions are not satisfied.
            if length(hi_sp_ind)<nlim
                udatB3{li} = [];
            end
        end
    end
    if isempty(phasicU3{li})
        udatB3{li} = [];
    end
end % All LC units for this session.

%%

for li = 1:n_LC
    
    if (~isempty(udat{li}))
        % if (~isempty(udat{li}) & ~isempty(udat3{li}))
        % if (isempty(udat{li}{nbs}{9}) | isempty(udat3{li}{nbs}{9}))
        if (isempty(udat{li}{nbs}{1}))
            % if (isempty(udat{li}) | isempty(udat3{li}))
            udat{li} = [];
            udat3{li} = [];
        end
    end
    
    if (isempty(udat{li}))
        % if (isempty(udat{li}) | isempty(udat3{li}))
        udat{li} = [];
        udat3{li} = [];
    end
end

% Collect results for saving:
% beep_rsc_ = {udat nL phasicU phasicU2 phasicU3 phasicU4 udat2 udat3 udat4 phasicRespU r_pupil};

% beep_rsc_ = {udat nL phasicU phasicU2 phasicU3 phasicU4 udat2 udat3 udat4 phasicRespU peak_pupil r_pupil};
%                     1        2          3           4             5             6        7        8       9             10              11          12

% beep_rsc_ = {udat phasicU phasicU3 udat3 phasicRespU peak_pupil r_pupil udatB udatB3 all_beep_pupil};
beep_rsc_ = {udat phasicU phasicU3 udat3 phasicRespU peak_pupil r_pupil udatB udatB3 all_beep_pupil baseline_pupil}; % Added baseline_pupil 11/24/21 after submitting ms.


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


