function [ccg_dat] = getLC_DualArea_Covariation(siteData)

% 090621
% Modified from "getLC_DualArea_sts"
% function [sts_dat_] = getLC_DualArea_sts(siteData)
%
% Calculates LC spike triggered ACC spiking and ACC pairwise correlations.

% 022221 - Sidd - Post-review changes
% 072617 - Sidd - Improved binning loop.
% 042917 - Sidd - Tweaked SfN version of code and re-analzyed.
% 062516 - Sidd
% 070216 - Sidd
% 110816 - Sidd
% 111016 - Sidd

% ********************************************
% Summary of standardized data structure:

%  siteData{1}: trialsxcols matrix, cols are:
%   1   ... fix start time wrt fixation on
%   2   ... fix end time wrt fix start time (fix duration)
%   3   ... reported correct
%   4   ... beep on time (when appropriate), wrt to fix start time
%   5   ... trial begin time, wrt to fix start time
%   6   ... trial end time, wrt to fix start time
%   7   ... trial wrt time (cpu clock)
%   8   ... LFP index corresponding to fix start time (coded above)
%   9   ... ELESTM on time (when appropriate), wrt fix start time

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

% ************************************************************************************

% Create 0.5 sec bins.

backTime1 = -1500; fwdTime1 = 0; % 091921b
% tmin = backTime1; tmax = fwdTime1; tsize  = 500; tstep = 50; % 101021 NOT great
tmin = backTime1; tmax = fwdTime1; tsize  = 250; tstep = 50; % 100821a
% tmin = backTime1; tmax = fwdTime1; tsize  = 500; tstep = 125; % Upto 100821
% tmin = backTime1; tmax = fwdTime1; tsize  = 250; tstep = 125;
binsXB = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
taxB = mean(binsXB,2);
nb = size(binsXB,1);

backTime2 = 0; fwdTime2 = 1500; % 091921b
% tmin = backTime2; tmax = fwdTime2; tsize  = 500; tstep = 50; % 101021 - NOT great
tmin = backTime2; tmax = fwdTime2; tsize  = 250; tstep = 50; % 100821
% tmin = backTime2; tmax = fwdTime2; tsize  = 500; tstep = 125; % Upto 100821
% tmin = backTime2; tmax = fwdTime2; tsize  = 250; tstep = 125;
binsXA = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
taxF = mean(binsXA,2);
nb = size(binsXA,1);

% tax2 = -(1050-tstep):tstep:(1050-tstep); % 100821
tax2 = -(1300-tstep):tstep:(1300-tstep); % 100821
% tax2 = -(1125-tstep):tstep:(1125-tstep); % Upto 100821
% tax2 = -(taxF(end)-tstep):tstep:(taxF(end)-tstep); % j

% tax2 = -(1125-tstep):tstep:(1125-tstep);
% backTime = -1500; fwdTime = 2500; % 500e,f,g,h,i
% backTime = 1000; fwdTime = 2000; % 091921a
% tmin = backTime; tmax = fwdTime; tsize  = 500; tstep = 100; % 500e.
% tmin = backTime; tmax = fwdTime; tsize  = 500; tstep = 250; % 500f,g,h,i.
% tmin = backTime; tmax = fwdTime; tsize  = 500; tstep = 250; % 500c... NOISY rsc with smaller bins.
% tax2 = -(875-tstep):tstep:(875-tstep); % 091921a
% tax2 = -(3250-tstep):tstep:(3250-tstep); % 091921a
% tax2 = -(tax(end)-tstep):tstep:(tax(end)-tstep); % j
% tax2 = -(3625-tstep):tstep:(3625-tstep); % 500e
% tax2 = -(3750-tstep):tstep:(3750-tstep); % 500e,f,g,h,i
% tax2 = -(2100-tstep):tstep:(2100-tstep);


% ************************************************************************************

%% Find no beep trials ...

LnoBeep = isnan(siteData{1}(:,4)); % Only non-beep trials
FnoBeep = find(LnoBeep);
num_noBeep = length(FnoBeep);

%%  Get unit information for this session:
% Note that with we are using channel 4 for LC and channels 9-16 for ACC.

sp_ch = siteData{7}{1};     % All spike channel id's.
LC_ch = find(sp_ch<5000);   % LC channel id's.
ACC_ch = find(sp_ch>8000); % ACC channel id's.

% How many units in LC?
n_LC = length(LC_ch);
n_ACC = length(ACC_ch);

% How many pairs in ACC?
ACC_pairs = nchoosek(ACC_ch,2);
n_ap = size(ACC_pairs,1);

numLim = 10;

%% Get timing information for this session:

fst = siteData{1}(LnoBeep,1); % Fixation start time (wrt fix ON).
fet = siteData{1}(LnoBeep,2); % Fixation end time (wrt FP ON).
fd = siteData{1}(LnoBeep,2); % Fixation duration.

tst = siteData{1}(LnoBeep,5); % Trial start time.
tet = siteData{1}(LnoBeep,6); % Trial end time.

% ************************************************************************************

%% Calculate LC spike triggered ACC spiking:

sps_LC = siteData{3}(LnoBeep,LC_ch); % LC spikes.
sps_ACC = siteData{3}(LnoBeep,ACC_ch); % ACC spikes.

% gti = find(fd>=1500);
gti = find(fd>=1500 & tst<=-1500);
fd=fd(gti);
tst = tst(gti);
tet = tet(gti);
ngt = length(gti);

%%

% pdata = siteData{2}(LnoBeep,:,5);
% pdata = pdata(gti,:);

%% ******************************

% % Get ACC firing rate in trial.
% bnACC = nans(n_ACC,nb);
% bnACCZ = nans(n_ACC,nb);
% for ai = 1:n_ACC % Iterate over ACC units.
%     bsACC = nans(1000,nb); % Initialize with nans... mmmmm nans!
%     ciii = 0;
%     for ti = 1:ngt;   % Iterate over trials.
%         ts_acc = sps_ACC{gti(ti),ai}; % ACC trial spikes.
%         % ts_acc_trial = ts_acc(ts_acc>=-tback & ts_acc<=fd(ti));
%         ts_acc_trial = ts_acc(ts_acc>=backTime & ts_acc<=fwdTime);
%         % ts_acc_trial = ts_acc; % 090621, Sidd.
%         if ~isempty(ts_acc_trial)
%             %if tst(ti)<=backTime & tet(ti)>=fwdTime
%                 ciii = ciii + 1;
%                 %                 si = find(binsX(:,1)>=tst(ti)); si = si(1);
%                 %                 ei = find(binsX(:,2)<=tet(ti)); ei = ei(end);
%                 for bi = 1:nb
%                     bsACC(ciii,bi) = sum(ts_acc_trial>=binsX(bi,1) & ts_acc_trial<=binsX(bi,2)); % Binned ACC spikes for this trial, for this LC spike.
%                 end % Binning loop.
%             % end % Check start and end bins.
%         end % Check for ACC spikes.
%     end % Trials for this unit.
%     bsACC(ciii+1:end,:) = [];
%     % Firing rate.
%     bsACC = (1000/tsize)*bsACC;
%     % Z-scored firing rate.
%     szz = size(bsACC); bsACCZ = reshape(zscore(bsACC(:)),szz);
%
%     bnACC(ai,:) = nanmean(bsACC);
%     bnACCZ(ai,:) = nanmean(bsACCZ);
%
% end % ACC units.
%
% % % Data lookit. Comment out when functionalize.
% % figure;
% % mmn = nanmean(bnACC); gbb = find(isfinite(mmn));
% % subplot(2,1,1); plot(tax(gbb),nanmean(bnACC(:,gbb)),'k-'); hold on;
% % errorBarFilled(tax(gbb)',nanmean(bnACC(:,gbb)),nanse(bnACC(:,gbb),1),[0.7 0.7 0.7]);
% %
% % mmn = nanmean(bnACCZ); gbb = find(isfinite(mmn));
% % subplot(2,1,2); plot(tax(gbb),nanmean(bnACCZ(:,gbb)),'k-'); hold on;
% % errorBarFilled(tax(gbb)',nanmean(bnACCZ(:,gbb)),nanse(bnACCZ(:,gbb),1),[0.7 0.7 0.7]);
%
% %% Get LC firing rate in trial.
% bnLC = nans(n_LC,nb);
% bnLCZ = nans(n_LC,nb);
% for ai = 1:n_LC % Iterate over ACC units.
%     bsLC = nans(1000,nb); % Initialize with nans... mmmmm nans!
%     ciii = 0;
%     for ti = 1:ngt;   % Iterate over trials.
%         ts_lc = sps_LC{gti(ti),ai}; % ACC trial spikes.
%         % ts_acc_trial = ts_acc(ts_acc>=-tback & ts_acc<=fd(ti));
%         ts_lc_trial = ts_lc; % 090621, Sidd.
%         if ~isempty(ts_lc_trial)
%             if tst(ti)<=backTime & tet(ti)>=fwdTime
%                 ciii = ciii + 1;
%                 for bi = 1:nb
%                     bsLC(ciii,bi) = sum(ts_lc_trial>=binsX(bi,1) & ts_lc_trial<=binsX(bi,2)); % Binned ACC spikes for this trial, for this LC spike.
%                 end % Binning loop.
%             end % Check start and end bins.
%         end % Check for ACC spikes.
%     end % Trials for this unit.
%     bsLC(ciii+1:end,:) = [];
%     % Firing rate.
%     bsLC = (1000/tsize)*bsLC;
%     % Z-scored firing rate.
%     szz = size(bsLC); bsLCZ = reshape(zscore(bsLC(:)),szz);
%
%     bnLC(ai,:) = nanmean(bsLC);
%     bnLCZ(ai,:) = nanmean(bsLCZ);
%
% end % LC units.
%
% % % Data lookit. Comment out when functionalize.
% % figure;
% % mmn = nanmean(bnLC); gbb = find(isfinite(mmn));
% % subplot(2,1,1); plot(tax(gbb),nanmean(bnLC(:,gbb)),'k-'); hold on;
% % errorBarFilled(tax(gbb)',nanmean(bnLC(:,gbb)),nanse(bnLC(:,gbb),1),[0.7 0.7 0.7]);
% %
% % mmn = nanmean(bnLCZ); gbb = find(isfinite(mmn));
% % subplot(2,1,2); plot(tax(gbb),nanmean(bnLCZ(:,gbb)),'k-'); hold on;
% % errorBarFilled(tax(gbb)',nanmean(bnLCZ(:,gbb)),nanse(bnLCZ(:,gbb),1),[0.7 0.7 0.7]);


%% LC YesNo

% n_spikes_trial_lcu = cell(1,n_LC);
% yesNo_spikes_trial_lcu = cell(1,n_LC);
% frate = cell(1,n_LC);
%
% for li = 1:n_LC % Iterate over LC units.
%     n_spikes_trial_lc = nans(ngt,1);
%     lc_spike_trial_yesNo = nans(ngt,1);
%     frate_unit = nans(ngt,1);
%     for ti = 1:ngt;   % Iterate over trials.
%         ts_lc = sps_LC{gti(ti),li}; % LC trial spikes.
%         n_spikes_trial_lc(ti) = length(ts_lc(ts_lc>=0 & ts_lc<=fd(ti)));
%         lc_spike_yesNo(ti) = isempty(ts_lc(ts_lc>=0 & ts_lc<=fd(ti)));
%         frate_unit(ti) = (1000/(fd(ti)))*length(ts_lc(ts_lc>=0 & ts_lc<=fd(ti)));
%     end
%     frate{li} = frate_unit;
%     n_spikes_trial_lcu{li} = n_spikes_trial_lc;
%     yesNo_spikes_trial_lcu{li} = lc_spike_yesNo;
% end

%% Get binned spikes for ACC units and calculate trial FR, var and FF.

ci = 0;
nA_allB = cell(1,n_ACC); sta_ACCB = nans(1000,nb);
nA_allA = cell(1,n_ACC); sta_ACCA = nans(1000,nb);
for ai = 1:n_ACC % Iterate over ACC units.
    bsb = nans(ngt,nb); % Initialize with nans... mmmmm nans!
    bsa = nans(ngt,nb); % Initialize with nans... mmmmm nans!
    for ti = 1:ngt;   % Iterate over trials.
        ts_acc = sps_ACC{gti(ti),ai}; % ACC trial spikes.
        ts_acc_trial = ts_acc(ts_acc>=backTime1 & ts_acc<=fwdTime2);
        % ts_acc_trial = ts_acc; % 090621, Sidd.
        if ~isempty(ts_acc_trial)
            for bi = 1:nb
                bsb(ti,bi) = sum(ts_acc_trial>=binsXB(bi,1) & ts_acc_trial<=binsXB(bi,2)); % Binned ACC spikes for this trial, for this LC spike.
                bsa(ti,bi) = sum(ts_acc_trial>=binsXA(bi,1) & ts_acc_trial<=binsXA(bi,2)); % Binned ACC spikes for this trial, for this LC spike.
            end % Binning loop.
        end % Check for empty spike structure.
    end % All trials for this LC unit, this ACC unit.
    
    % if (size(bsb,1)>=10 & size(bsa,1)>=10)
    ci = ci+1;
    sta_ACCB(ci,:) = (1000/tsize)*nanmean(bsb);
    sta_ACCA(ci,:) = (1000/tsize)*nanmean(bsa);
    % end
    nA_allB{ai} = bsb; % LC aligned spikes for each ACC unit - there is one cell for each ACC unit; each cell has n_LC cells.
    nA_allA{ai} = bsa; % LC aligned spikes for each ACC unit - there is one cell for each ACC unit; each cell has n_LC cells.
end % ACC unit loop.

% Cleanup.
sta_ACCB(ci+1:end,:) = [];
sta_ACCA(ci+1:end,:) = [];

%% Get binned spikes for LC units and calculate trial FR, var and FF.

ci = 0;
zeib = cell(1,2); nzeib = cell(1,2);
zeia = cell(1,2); nzeia = cell(1,2);
sta_LCB = nans(10,nb);
sta_LCA = nans(10,nb);
for ai = 1:n_LC % Iterate over LC units.
    bsb = nans(ngt,nb); % Initialize with nans... mmmmm nans!
    bsa = nans(ngt,nb); % Initialize with nans... mmmmm nans!
    for ti = 1:ngt  % Iterate over trials.
        ts_lc = sps_LC{gti(ti),ai}; % LC trial spikes.
        ts_lc_trial = ts_lc(ts_lc>=backTime1 & ts_lc<=fwdTime2);
        if ~isempty(ts_lc_trial)
            for bi = 1:nb
                bsb(ti,bi) = sum(ts_lc_trial>=binsXB(bi,1) & ts_lc_trial<=binsXB(bi,2)); % Binned ACC spikes for this trial, for this LC spike.
                bsa(ti,bi) = sum(ts_lc_trial>=binsXA(bi,1) & ts_lc_trial<=binsXA(bi,2)); % Binned ACC spikes for this trial, for this LC spike.
            end % Binning loop.
        end % Check for empty spike structure.
    end % All trials for this LC unit, this ACC unit.
    
    % nL{ai} = bs; % LC aligned spikes for each ACC unit - there is one cell for each ACC unit; each cell has n_LC cells.
    zeib{ai} = find(~isfinite(sum(bsb,2)));
    nzeib{ai} = find(isfinite(sum(bsb,2)));
    zeia{ai} = find(~isfinite(sum(bsa,2)));
    nzeia{ai} = find(isfinite(sum(bsa,2)));
    ci = ci+1;
    sta_LCB(ci,:) = (1000/tsize)*nanmean(bsb(nzeib{ai},:));
    sta_LCA(ci,:) = (1000/tsize)*nanmean(bsa(nzeia{ai},:));
end % ACC unit loop.

% Cleanup.
sta_LCB(ci+1:end,:) = [];
sta_LCA(ci+1:end,:) = [];

%% Get unconditioned ACC spike counts for rsc calculation.

% gti = find(fd>=2100);
% ngt = length(gti);

nAA = [];

for ai = 1:n_ACC
    ats = sps_ACC(gti,ai);
    nts = [];
    for ti = 1:ngt
        ts = ats{ti};
        ts_fix = ts(ts>1001 & ts <=2100); % Drop 1sec of fixation at beginning of stable fixation period.
        nts(ti) = length(ts_fix);
    end % Iterate over trials.
    nAA{ai} = nts;
end % Iterate over ACC units.

%% Troubleshooting plots:
%
% figure; hold on;
% plot(tax,nanmean(sta_ACC_LC-sta_ACC_LC_R),'k-');
%
% figure; hold on;
% nn = size(sta_ACC_LC,1)
% for i = 1:nn
%     rr = sta_ACC_LC(i,:)
%     plot(tax,nanmean(rr),'k-');
% end

%% Now calculate the pairwise correlations in ACC.

msr = nans(n_ap,2);
allcount = 0;
% Regular trials.
for li = 1:n_LC
    for ai = 1:n_ap % Iterate over pairs:
        
        a1 = find(ACC_pairs(ai,1) == ACC_ch);
        a2 = find(ACC_pairs(ai,2) == ACC_ch);
        
        % Binned ACC spikes, aligned to start of stable fixation.
        ap1B = nA_allB{a1}; % Get spike count for 1st unit in pair.
        ap2B = nA_allB{a2}; % Get spike count for 2nd unit in pair.
        ap1A = nA_allA{a1}; % Get spike count for 1st unit in pair.
        ap2A = nA_allA{a2}; % Get spike count for 2nd unit in pair.
        ap1All = nAA{a1}; % Get spike count for 1st unit in pair.
        ap2All = nAA{a2}; % Get spike count for 2nd unit in pair.
        
        % msr(ai,:) = [(1000/(backTime+fwdTime))/nanmean(nanmean(ap1)) (1000/(backTime+fwdTime))/nanmean(nanmean(ap2))];
        out_rsc_nzb = nans(nb,1);
        out_rsc_zb = nans(nb,1);
        out_rsc_nza = nans(nb,1);
        out_rsc_za = nans(nb,1);
        
        if (~isempty(ap1B) & ~isempty(ap2B))
            for kkk = 1:nb
                
                % ALL LC nonzero trials - before start of stable fixation.
                al_1 = ap1B(nzeib{li},kkk); al_2 = ap2B(nzeib{li},kkk);
                ci = find(~isnan(al_1) & ~isnan(al_2));
                if length(ci)>=5
                    al_1 = al_1(ci); % ACC pair unit 1.
                    al_2 = al_2(ci); % ACC pair unit 2.
                    rsc_out = get_pwc(al_1,al_2,numLim);
                    out_rsc_nzb(kkk) = rsc_out.RSC;
                else
                    out_rsc_nzb(kkk) = nan;
                end
                
                % ALL LC zero trials - before start of stable fixation.
                al_1 = ap1B(zeib{li},kkk); al_2 = ap2B(zeib{li},kkk);
                ci = find(~isnan(al_1) & ~isnan(al_2));
                if length(ci)>=5
                    al_1 = al_1(ci); % ACC pair unit 1.
                    al_2 = al_2(ci); % ACC pair unit 2.
                    rsc_out = get_pwc(al_1,al_2,numLim);
                    out_rsc_zb(kkk) = rsc_out.RSC;
                else
                    out_rsc_zb(kkk) = nan;
                end
                
            end % Bins loop.
        end % Check for empty ninned ACC spikes.
        
        if ~isempty(ap1A) && ~isempty(ap2A)
            for kkk = 1:nb
                % ALL LC nonzero trials - during stable fixation.
                
                al_1 = ap1A(nzeia{li},kkk); al_2 = ap2A(nzeib{li},kkk);
                ci = find(~isnan(al_1) & ~isnan(al_2));
                if length(ci)>=5
                    al_1 = al_1(ci); % ACC pair unit 1.
                    al_2 = al_2(ci); % ACC pair unit 2.
                    rsc_out = get_pwc(al_1,al_2,numLim);
                    out_rsc_nza(kkk) = rsc_out.RSC;
                else
                    out_rsc_nza(kkk) = nan;
                end
                
                % ALL LC zero trials - during stable fixation..
                al_1 = ap1A(zeia{li},kkk); al_2 = ap2A(zeia{li},kkk);
                ci = find(~isnan(al_1) & ~isnan(al_2));
                if length(ci)>=5
                    al_1 = al_1(ci); % ACC pair unit 1.
                    al_2 = al_2(ci); % ACC pair unit 2.
                    rsc_out = get_pwc(al_1,al_2,numLim);
                    out_rsc_za(kkk) = rsc_out.RSC;
                else
                    out_rsc_za(kkk) = nan;
                end
                
            end % Bins loop.
        end % Check for empty ninned ACC spikes.
        
        % Call "get_pwc" to calculate the pairwise correlations in spike counts.
        rsc_all = get_pwc(reshape(ap1All.',1,[]),reshape(ap2All.',1,[]),numLim);
        
        % out_rsc_uA{ai} = [rsc_all.RSC rsc_all.pval];
        out_rsc_uA{ai} = rsc_all.RSC;
        out_rsc_uNZB{ai} = out_rsc_nzb;
        out_rsc_uNZA{ai} = out_rsc_nza;
        
        out_rsc_uZB{ai} = out_rsc_zb;
        out_rsc_uZA{ai} = out_rsc_za;
        
    end % ACC Pairs loop.
    % pause
    out_rsc_all{li} = {out_rsc_uNZB out_rsc_uZB out_rsc_uNZA out_rsc_uZA out_rsc_uA};
end

%%
% Triangle function to normalize for number of bins as per Bair et al (2001).
tBack = nb;
TT = (1/tBack)*(tBack-(abs(-(tBack-1):(tBack-1))));

%% Calculate the cross-correlation between LC spiking and ACC rsc.

out_FRL_AB = nans(n_LC,nb);
out_FRL_AA = nans(n_LC,nb);
out_FRA_AB = nans(n_ACC,nb);
out_FRA_AA = nans(n_ACC,nb);

for li = 1:n_LC
    
    out_ccgB = nans(n_ap,2*nb-1); out_rsc_ABnz = nans(n_ap,nb); out_rsc_ABz = nans(n_ap,nb);
    out_ccgA = nans(n_ap,2*nb-1); out_rsc_AAnz = nans(n_ap,nb); out_rsc_AAz = nans(n_ap,nb);
    out_ccgAB = nans(n_ACC,2*nb-1); out_ccgAA = nans(n_ACC,2*nb-1);
    
    LC_spikesB = sta_LCB(li,:)'; out_FRL_AB(li,:) = LC_spikesB;
    LC_spikesA = sta_LCA(li,:)'; out_FRL_AA(li,:) = LC_spikesA;
    
    for ai = 1:n_ACC
        ACC_spikesB = sta_ACCB(ai,:)'; ACC_spikesA = sta_ACCA(ai,:)';
        raw_ccg = xcorr(ACC_spikesB,LC_spikesB,'unbiased'); out_ccgAB(ai,:) = raw_ccg';
        raw_ccg = xcorr(ACC_spikesA,LC_spikesA,'unbiased'); out_ccgAA(ai,:) = raw_ccg';
        if li == 1
            out_FRA_AB(ai,:) = ACC_spikesB; out_FRA_AA(ai,:) = ACC_spikesA;
        end
    end
    
    for jj = 1:n_ap
        
        ACC_rsc = out_rsc_all{li}{1}{jj};
        %         raw_ccg = xcorr(ACC_rsc,LC_spikesB,'unbiased');
        %         out_ccgB(jj,:) = raw_ccg';
        out_rsc_ABnz(jj,:) = ACC_rsc;
        out_rsc_ABz(jj,:)= out_rsc_all{li}{2}{jj};
        
        ACC_rsc = out_rsc_all{li}{3}{jj};
        %         raw_ccg = xcorr(ACC_rsc,LC_spikesA,'unbiased');
        %         out_ccgA(jj,:) = raw_ccg';
        out_rsc_AAnz(jj,:) = ACC_rsc;
        out_rsc_AAz(jj,:)= out_rsc_all{li}{4}{jj};
        
    end
    out_ccgB = []; out_ccgA = [];
    out_ccg_all{li} = {out_ccgB out_ccgA out_FRL_AB out_FRL_AA out_rsc_ABnz out_rsc_AAnz out_rsc_ABz out_rsc_AAz out_FRA_AB out_FRA_AA out_ccgAB out_ccgAA};
end


%%

% % Test plots.
% figure; hold on;
% mean_ccg = nanmedian(out_ccg);
% std_ccg = nanstd(out_ccg,1);
% plot(tax2,mean_ccg,'k-');
% errorBarFilled(tax2,mean_ccg,std_ccg,[0.7 0.7 0.7]);

%% Write 'em out:

% ccg_dat = {out_ccg out_rsc_A out_FR_A tax tax2};
% ccg_dat = {out_rsc_all sta_LC tax tax2 TT out_ccg out_rsc_A out_FR_A};

ccg_dat = {out_ccg_all taxB taxF tax2 TT};
% out_ccg_all{li} = {out_ccgB out_ccgA out_FR_AB out_FR_AA out_rsc_ABnz out_rsc_AAnz out_rsc_ABz out_rsc_AAz};

% out_rsc_all{li} = {out_rsc_uNZ out_rsc_uZ out_rsc_uA};
% out_ccg_all{li} = {out_ccg out_FR_A out_rsc_Anz out_rsc_Az};










