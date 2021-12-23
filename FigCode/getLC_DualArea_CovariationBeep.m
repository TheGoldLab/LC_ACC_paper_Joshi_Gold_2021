function [ccg_dat] = getLC_DualArea_CovariationBeep(siteData)

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

backTime = 0; fwdTime = 1000; % j
% tmin = backTime; tmax = fwdTime; tsize  = 500; tstep = 125; % 500f,g,h,i.
% tmin = backTime; tmax = fwdTime; tsize  = 250; tstep = 125; % Upto 100821.
tmin = backTime; tmax = fwdTime; tsize  = 250; tstep = 50;
binsX = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
taxF = mean(binsX,2);
nb = size(binsX,1);
tax2F = -(800-tstep):tstep:(800-tstep);
% tax2F = -(taxF(end)-tstep):tstep:(taxF(end)-tstep); % j
% tax2F = tax2F(2:end-1);
% tax2A= -(taxF(end)-tstep):tstep:(taxF(end)-tstep); % j

backTime2 = 1000; fwdTime2 = 0; % j
% tmin = -backTime2; tmax = fwdTime2; tsize  = 500; tstep = 125; % 500f,g,h,i.
% tmin = -backTime2; tmax = fwdTime2; tsize  = 250; tstep = 125; % Upto 100821.
tmin = -backTime2; tmax = fwdTime2; tsize  = 250; tstep = 50;
binsX2 = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
taxB = mean(binsX2,2);
nb = size(binsX,1);
% tax2B = -(taxB(end)-tstep):tstep:(taxB(end)-tstep); % j

% tax2 = -(3625-tstep):tstep:(3625-tstep); % 500e
% tax2 = -(3750-tstep):tstep:(3750-tstep); % 500e,f,g,h,i


% ************************************************************************************

%% Find no beep trials ...

LBeep = ~isnan(siteData{1}(:,4)); % Only beep trials
LNoBeep = isnan(siteData{1}(:,4)); % Only non-beep trials
FBeep = find(LBeep);
num_Beep = length(FBeep);

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

numLim = 5;

%% Get timing information for this session:

fst = siteData{1}(LBeep,1); % Fixation start time (wrt fix ON).
fet = siteData{1}(LBeep,2); % Fixation end time (wrt FP ON).
fd = siteData{1}(LBeep,2); % Fixation duration.

tst = siteData{1}(LBeep,5); % Trial start time.
tet = siteData{1}(LBeep,6); % Trial end time.

bTime = siteData{1}(LBeep,4);

% ************************************************************************************

%% Calculate LC spike triggered ACC spiking:

sps_LC = siteData{3}(LBeep,LC_ch); % LC spikes.
sps_ACC = siteData{3}(LBeep,ACC_ch); % ACC spikes.

% gti1 = find((fd-bTime)>=1000);
% gti2 = find((fd-bTime)>=1100);
% gti1'
% gti2'

% gti = find(fd>=fwdTime);
gti = find((fd-bTime)>=fwdTime);

fd=fd(gti);
tst = tst(gti);
tet = tet(gti);
bTime = bTime(gti);
ci = 0;
ngt = length(gti);

%% LC: Define phasic/no-phasic trials.

phasicU = cell(1,n_LC);
phasicU3 = cell(1,n_LC);

for nu = 1:n_LC
    prebeepspikes = nans(1,length(gti));
    phasicYesNo = nans(1,length(gti));
    phasicYesNo3 = nans(1,length(gti));
    
    gzi = []; nts = [];
    ats = sps_LC(gti,nu);
    bs = nans(length(gti),length(binsX)); % Preallocate.
    for ti = 1:ngt
        ts = ats{ti} - bTime(ti); % Spike times were recoded wrt mnky fixation start.
        ts_fix = ts(ts>-1000 & ts < 1000);
        
        ts_fix_postBeep = (1000/200)*length(ts(ts > 200 & ts <= 400)); % Get spikes from the PBI epoch.
        %         ts_fix_postBeep_all = (1000/1000)*length(ts(ts > 0 & ts <= 1000)); % Get spikes from PB epoch.
        ts_fix_Beep = (1000/200)*length(ts(ts > 0 & ts <= 200)); % Get spikes from the phasic epoch.
        ts_fix_preBeep = (1000/1000)*length(ts(ts > -1000 & ts < 0)); % Get spikes from the pre-beep epoch.
        
        % Classify beep responses (or lack of response).
        % if ts_fix_preBeep<ts_fix_Beep && ts_fix_postBeep<=ts_fix_preBeep % Phasic response and post-beep inhibition. 041020.
        if ts_fix_preBeep<ts_fix_Beep && ts_fix_postBeep<=ts_fix_preBeep % Phasic response and post-beep inhibition. 041020.
            phasicYesNo(ti) = 1;
            prebeepspikes(ti) = ts_fix_preBeep;
        end
        if ts_fix_preBeep>=ts_fix_Beep % NO Phasic response.
            phasicYesNo3(ti) = 1;
            prebeepspikes(ti) = ts_fix_preBeep;
        end
    end % Trials loop.
    
    phasicU{nu} = phasicYesNo;
    phasicU3{nu} = phasicYesNo3;
    
end

%% Get binned spikes for ACC units and calculate trial FR, var and FF.

nA_allA = cell(1,n_ACC);
nA_allB = cell(1,n_ACC);
% sta_ACC = nans(1000,nb); FF_ACC = nans(1000,nb); var_out_ACC = nans(1000,nb); mean_out_ACC = nans(1000,nb); ci = 0;
for ai = 1:n_ACC % Iterate over ACC units.
    bsB = nans(ngt,nb); % Initialize with nans... mmmmm nans!
    bsA = nans(ngt,nb); % Initialize with nans... mmmmm nans!
    for ti = 1:ngt;   % Iterate over trials.
        ts_acc = sps_ACC{gti(ti),ai}; % ACC trial spikes.
        ts_acc_beep = ts_acc - bTime(ti);
        if ~isempty(ts_acc)
            for bi = 1:nb
                bsA(ti,bi) = sum(ts_acc_beep>=binsX(bi,1) & ts_acc_beep<=binsX(bi,2)); % Binned ACC spikes for this trial, for this LC spike.
                bsB(ti,bi) = sum(ts_acc_beep>=binsX2(bi,1) & ts_acc_beep<=binsX2(bi,2)); % Binned ACC spikes for this trial, for this LC spike.
            end % Binning loop.
        end % Check for empty spike structure.
    end % All trials for this LC unit, this ACC unit.
    
    % if size(bsA,1)>=5 & size(bsB,1)>=5
    nAA = bsA; % LC aligned spikes for this ACC unit. There is one cell for each LC unit.
    nAB = bsB; % LC aligned spikes for this ACC unit. There is one cell for each LC unit.
    %     end
    nA_allA{ai} = nAA; % LC aligned spikes for each ACC unit - there is one cell for each ACC unit; each cell has n_LC cells.
    nA_allB{ai} = nAB; % LC aligned spikes for each ACC unit - there is one cell for each ACC unit; each cell has n_LC cells.
end % ACC unit loop.

%% Get binned spikes for LC units and calculate trial FR, var and FF.

nL_allA = cell(1,n_LC);
nL_allB = cell(1,n_LC);
sta_LCAp = nans(1000,nb);
sta_LCBp = nans(1000,nb);
sta_LCAnp = nans(1000,nb);
sta_LCBnp = nans(1000,nb);
% FF_LC = nans(1000,nb); var_out_LC = nans(1000,nb); mean_out_LC = nans(1000,nb);
ci = 0;
for ai = 1:n_LC % Iterate over LC units.
    bsB = nans(ngt,nb); % Initialize with nans... mmmmm nans!
    bsA = nans(ngt,nb); % Initialize with nans... mmmmm nans!
    for ti = 1:ngt;   % Iterate over trials.
        ts_lc = sps_LC{gti(ti),ai}; % LC trial spikes.
        ts_lc_beep = ts_lc - bTime(ti);
        if ~isempty(ts_lc)
            for bi = 1:nb
                bsA(ti,bi) = sum(ts_lc_beep>=binsX(bi,1) & ts_lc_beep<=binsX(bi,2)); % Binned ACC spikes for this trial, for this LC spike.
                bsB(ti,bi) = sum(ts_lc_beep>=binsX2(bi,1) & ts_lc_beep<=binsX2(bi,2)); % Binned ACC spikes for this trial, for this LC spike.
            end % Binning loop.
        end % Check for empty spike structure.
    end % All trials for this LC unit, this ACC unit.
    
    phasic_index = find(phasicU{ai}==1); np = nansum(phasic_index);
    no_phasic_index = find(phasicU3{ai}==1); nnp = nansum(no_phasic_index);
    if np>=5 & nnp>=5
        ci = ci+1;
        nLA = bsA; % LC aligned spikes for this ACC unit. There is one cell for each LC unit.
        nLB = bsB; % LC aligned spikes for this ACC unit. There is one cell for each LC unit.
        sta_LCAp(ci,:) = (1000/tsize)*nanmean(bsA(phasic_index,:));
        sta_LCBp(ci,:) = (1000/tsize)*nanmean(bsB(phasic_index,:));
        sta_LCAnp(ci,:) = (1000/tsize)*nanmean(bsA(no_phasic_index,:));
        sta_LCBnp(ci,:) = (1000/tsize)*nanmean(bsB(no_phasic_index,:));
    end
    %     nL_allA{ai} = nLA; % LC aligned spikes for each ACC unit - there is one cell for each ACC unit; each cell has n_LC cells.
    %     nL_allB{ai} = nLB; % LC aligned spikes for each ACC unit - there is one cell for each ACC unit; each cell has n_LC cells.
end % ACC unit loop.

% Cleanup.
sta_LCAp(ci+1:end,:) = [];
sta_LCBp(ci+1:end,:) = [];
sta_LCAnp(ci+1:end,:) = [];
sta_LCBnp(ci+1:end,:) = [];

% FF_LC(ci+1:end,:) = [];
% var_out_LC(ci+1:end,:) = [];
% mean_out_LC(ci+1:end,:) = [];
%
% % pause
%
% % figure(2); hold on;
% % subplot(4,1,1); plot(nanmean(sta_LC),'k-');
% % subplot(4,1,2); plot(nanmean(FF_LC),'k-');
% % subplot(4,1,3); plot(nanmean(var_out_LC),'k-');
% % subplot(4,1,4); plot(nanmean(mean_out_LC),'k-');

%% Get unconditioned ACC spike counts for rsc calculation.

fd2 = siteData{1}(LNoBeep,2); % Fixation duration.
gti2 = find(fd2>=2100);
ngt2 = length(gti2);

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

%% Now calculate pairwise correlations in ACC.

msr = nans(n_ap,2);
allcount = 0;
% Regular trials:
for li = 1:n_LC
    phasic_index = find(phasicU{li}==1); np = nansum(phasic_index);
    no_phasic_index = find(phasicU3{li}==1); nnp = nansum(no_phasic_index);
    
    for ai = 1:n_ap % Iterate over pairs:
        
        a1 = find(ACC_pairs(ai,1) == ACC_ch);
        a2 = find(ACC_pairs(ai,2) == ACC_ch);
        
        % Binned ACC spikes, aligned to start of stable fixation.
        ap1A = nA_allA{a1}; % Get spike count for 1st unit in pair.
        ap2A = nA_allA{a2}; % Get spike count for 2nd unit in pair.
        ap1B = nA_allB{a1}; % Get spike count for 1st unit in pair.
        ap2B = nA_allB{a2}; % Get spike count for 2nd unit in pair.
        ap1All = nAA{a1}; % Get spike count for 1st unit in pair.
        ap2All = nAA{a2}; % Get spike count for 2nd unit in pair.
        
        % msr(ai,:) = [(1000/(backTime+fwdTime))/nanmean(nanmean(ap1)) (1000/(backTime+fwdTime))/nanmean(nanmean(ap2))];
        
        out_rscAp = nans(nb,1);
        out_rscBp = nans(nb,1);
        out_rscAnp = nans(nb,1);
        out_rscBnp = nans(nb,1);
        
        if ~isempty(ap1A) && ~isempty(ap2A) && ~isempty(ap1B) && ~isempty(ap2B)
            
            for kkk = 1:nb
                
                % ALL trials, LC phasic.
                al_1Ap = ap1A(phasic_index,kkk); al_2Ap = ap2A(phasic_index,kkk);
                al_1Bp = ap1B(phasic_index,kkk); al_2Bp = ap2B(phasic_index,kkk);
                cip = find(~isnan(al_1Ap) & ~isnan(al_2Ap) & ~isnan(al_1Bp) & ~isnan(al_2Bp));
                if length(cip)>=5
                    
                    al_1Ap = al_1Ap(cip); % ACC pair unit 1.
                    al_2Ap = al_2Ap(cip); % ACC pair unit 2.
                    rsc_outAp = get_pwc(al_1Ap,al_2Ap,numLim);
                    out_rscAp(kkk) = rsc_outAp.RSC;
                    
                    al_1Bp = al_1Bp(cip); % ACC pair unit 1.
                    al_2Bp = al_2Bp(cip); % ACC pair unit 2.
                    rsc_outBp = get_pwc(al_1Bp,al_2Bp,numLim);
                    out_rscBp(kkk) = rsc_outBp.RSC;
                    
                else
                    out_rscAp(kkk) = nan;
                    out_rscBp(kkk) = nan;
                end
                
                % ALL trials, LC no phasic.
                al_1Anp = ap1A(no_phasic_index,kkk); al_2Anp = ap2A(no_phasic_index,kkk);
                al_1Bnp = ap1B(no_phasic_index,kkk); al_2Bnp = ap2B(no_phasic_index,kkk);
                cinp = find(~isnan(al_1Anp) & ~isnan(al_2Anp) & ~isnan(al_1Bnp) & ~isnan(al_2Bnp));
                if length(cinp)>=5
                    
                    al_1Anp = al_1Anp(cinp); % ACC pair unit 1.
                    al_2Anp = al_2Anp(cinp); % ACC pair unit 2.
                    rsc_outAnp = get_pwc(al_1Anp,al_2Anp,numLim);
                    out_rscAnp(kkk) = rsc_outAnp.RSC;
                    
                    al_1Bnp = al_1Bnp(cinp); % ACC pair unit 1.
                    al_2Bnp = al_2Bnp(cinp); % ACC pair unit 2.
                    rsc_outBnp = get_pwc(al_1Bnp,al_2Bnp,numLim);
                    out_rscBnp(kkk) = rsc_outBnp.RSC;
                    
                else
                    out_rscAnp(kkk) = nan;
                    out_rscBnp(kkk) = nan;
                end
                
            end % bins loop
        end
        
        % Call "get_pwc" to calculate the pairwise correlations in spike counts.
        rsc_all = get_pwc(reshape(ap1All.',1,[]),reshape(ap2All.',1,[]),numLim);
        
        % out_rsc_uA{ai} = [rsc_all.RSC rsc_all.pval];
        out_rsc_uAll{ai} = rsc_all.RSC;
        out_rsc_uAp{ai} = out_rscAp;
        out_rsc_uBp{ai} = out_rscBp;
        out_rsc_uAnp{ai} = out_rscAnp;
        out_rsc_uBnp{ai} = out_rscBnp;
        
    end % ACC Pairs loop.
    
    out_rsc_allAp{li} = out_rsc_uAp;
    out_rsc_allBp{li} = out_rsc_uBp;
    out_rsc_allAnp{li} = out_rsc_uAnp;
    out_rsc_allBnp{li} = out_rsc_uBnp;
    out_rsc_uA_all = out_rsc_uAll;
    
end % LC unit loop.

% pause

%%
% Triangle function to normalize for number of bins as per Bair et al (2001).
tBack = nb;
TT = (1/tBack)*(tBack-(abs(-(tBack-1):(tBack-1))));


%% Calculate the cross-correlation between LC spiking and ACC rsc.

for li = 1:n_LC
    
    out_ccgBp = nans(n_ap,2*nb-1);
    out_ccgAp = nans(n_ap,2*nb-1);
    out_rsc_Ap = nans(n_ap,nb);
    out_rsc_Bp = nans(n_ap,nb);
    
    out_ccgBnp = nans(n_ap,2*nb-1);
    out_ccgAnp = nans(n_ap,2*nb-1);
    out_rsc_Anp = nans(n_ap,nb);
    out_rsc_Bnp = nans(n_ap,nb);
    
    out_FR_Ap = nans(n_LC,nb);
    out_FR_Bp = nans(n_LC,nb);
    out_FR_Anp = nans(n_LC,nb);
    out_FR_Bnp = nans(n_LC,nb);
    
    LC_spikesBp = sta_LCBp(li,:)'; LC_spikesAp = sta_LCAp(li,:)';
    LC_spikesBnp = sta_LCBnp(li,:)'; LC_spikesAnp = sta_LCAnp(li,:)';
    
    % z-score.
    % ZSp= zscore([LC_spikesBp; LC_spikesAp]);
    % ZSnp= zscore([LC_spikesBnp; LC_spikesAnp]);
    
    % Don't z-score.
    ZSp= [LC_spikesBp; LC_spikesAp];
    ZSnp= [LC_spikesBnp; LC_spikesAnp];
    
    out_FR_Bp(li,:) = ZSp(1:nb);
    out_FR_Ap(li,:) = ZSp(nb+1:end);
    
    out_FR_Bnp(li,:) = ZSnp(1:nb);
    out_FR_Anp(li,:) = ZSnp(nb+1:end);
    
    for jj = 1:n_ap
        
        ACC_rscAp = out_rsc_allAp{li}{jj};
        ACC_rscBp = out_rsc_allBp{li}{jj};
        ACC_rscAnp = out_rsc_allAnp{li}{jj};
        ACC_rscBnp = out_rsc_allBnp{li}{jj};
        
        % Calculate CCG and ACG for LC zero spiking condition.
        
        %         msr_ze1 = msr(jj,1); % Mean spike rate for ACC pair 1, low LC spiking.
        %         msr_ze2 = msr(jj,2); % Mean spike rate for ACC pair 1, low LC spiking.
        
        % raw_ccg = xcorr(zscore(ACC_rsc),zscore(LC_spikes));
        
        raw_ccgBp = xcorr((ACC_rscAp),ZSp(1:nb),'unbiased');
        raw_ccgAp = xcorr((ACC_rscBp),ZSp(nb+1:end),'unbiased');
        raw_ccgBnp = xcorr((ACC_rscAnp),ZSnp(1:nb),'unbiased');
        raw_ccgAnp = xcorr((ACC_rscBnp),ZSnp(nb+1:end),'unbiased');
        
        out_ccgBp(jj,:) = raw_ccgBp';
        out_ccgAp(jj,:) = raw_ccgAp';
        out_ccgBnp(jj,:) = raw_ccgBnp';
        out_ccgAnp(jj,:) = raw_ccgAnp';
        
        %         out_ccgBp(jj,:) = raw_ccgBp'./TT;
        %         out_ccgAp(jj,:) = raw_ccgAp'./TT;
        %         out_ccgBnp(jj,:) = raw_ccgBnp'./TT;
        %         out_ccgAnp(jj,:) = raw_ccgAnp'./TT;
        
        out_rsc_ABp(jj,:) = ACC_rscBp;
        out_rsc_AAp(jj,:) = ACC_rscAp;
        out_rsc_ABnp(jj,:) = ACC_rscBnp;
        out_rsc_AAnp(jj,:) = ACC_rscAnp;
        
    end
    
    LC_dat{li} = {out_ccgBp out_ccgAp out_ccgBnp out_ccgAnp out_rsc_ABp out_rsc_AAp out_rsc_ABnp out_rsc_AAnp out_FR_Bp out_FR_Ap out_FR_Bnp out_FR_Anp};
    
end



% % Test plots.
% figure; hold on;
% mean_ccg = nanmedian(out_ccg);
% std_ccg = nanstd(out_ccg,1);
% plot(tax2,mean_ccg,'k-');
% errorBarFilled(tax2,mean_ccg,std_ccg,[0.7 0.7 0.7]);

%% Write 'em out:

% ccg_dat = {out_ccg out_rsc_A out_FR_A tax tax2};
ccg_dat = {LC_dat taxF taxB tax2F TT};
% LC_dat = {out_ccgBp out_ccgAp out_ccgBnp out_ccgAnp out_rsc_ABp out_rsc_AAp out_rsc_ABnp out_rsc_AAnp out_FR_Bp out_FR_Ap out_FR_Bnp out_FR_Anp};











