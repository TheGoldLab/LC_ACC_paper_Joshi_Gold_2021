function [pwc_] = getACC_DualArea_rsc2(siteData,numLim,fNum)
%
% function [pwc_] = getACC_DualArea_rsc2(siteData)
%
% Calculate pairwise trial spike-count correlations between LC neurons, conditioned on ACC activity.
%
% Origin: 062516 - Sidd.
% History:
% Mod: 120219: Latest iteration for Fig 2 supplementary.
% Mod: 101618: adapted from "getLC_DualArea_rsc".
% Mod: 082518 - Sidd - now uses analysis with LC spiking divided into 3 classes - zero spikes and median split on remainder.
% Mod - tweaked code for rsc analysis on binned spike counts... output cells are slightly different now.
% Mod - 060817 - Added analysis for measuring the effect of window size on spike count and on variability (FF and rsc).
% Mod - 050817 - Wrote out spike rates for mean-matched plots.
% Mod - 050317 - Wrote out ISI's for CV calculation.
% Mod - 071516 - Wrote out LC conditioned ACC spike counts.
% Mod - 101416 - Fixed the duplication of count vectors in data structure (it was previously writing these out twice).
%
% Fixation task.
% "rsc" is the term used by Bair et al and Kohn et al for pairwise spike-count (sc) correlations (r).
%
% Each session is divided (median split) into high and low LC spike count trials.
%
% Calls "get_pwc" to do the pairwise correlation calculation.


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

% ************************************************************************************
% ************************************************************************************

%% Find no beep trials ...
LnoBeep    = isnan(siteData{1}(:,4)); % Only non-beep trials
FnoBeep    = find(LnoBeep);
num_noBeep = length(FnoBeep);

%% Set up bins ...

% bs = linspace(100,1000,10); % 032220b.
bs = linspace(200,1000,5); % April 2021.
backTime = 1001; fwdTime = 2100; % 032320b.
tmin   = backTime; tmax   = fwdTime; xBin = []; nb = []; nbs = length(bs);
for i = 1:nbs
    tsize  = bs(i); tstep  = tsize/2;
    xBin{i}  = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
    nb{i} = size(xBin{i},1);
end

winSiz = length(backTime:fwdTime)/1000;

%% Get pupil data ...

% % Sidd- 062316:
% % Simplified for now - use slope only where applicable...

% if nargin < 2 || isempty(ptype) || strcmp(ptype, 'pupil')
%     % default -- pupil diameter
% pdat = siteData{2}(FnoBeep,:,4); % trial/sample
% else
%     % pupil slope
% pdat = siteData{2}(FnoBeep,:,5); % trial/sample
% end
% num_samps = size(pdat,2);

%%  Get unit information for this session:
% Note that with we are using channel 4 for LC and channels 9-16 for ACC.

sp_ch = siteData{7}{1};      % All spike channel id's.
LC_ch = find(sp_ch<5000);   % LC channel id's.
ACC_ch = find(sp_ch>8000); % ACC channel id's.

% How many units in LC and ACC?
n_LC = length(LC_ch);
n_ACC = length(ACC_ch);

if n_LC>1
    % How many pairs in LC?
    LC_pairs = nchoosek(LC_ch,2);
    n_ap = size(LC_pairs,1);
    
    % ************************************************************************************
    
    %% Get timing information for this session:
    
    fst = siteData{1}(LnoBeep,1); % Fixation start time (wrt fix ON).
    fet = siteData{1}(LnoBeep,2); % Fixation end time (wrt FP ON).
    fd = siteData{1}(LnoBeep,2); % Fixation duration.
    
    tst = siteData{1}(LnoBeep,5); % Trial start time.
    tet = siteData{1}(LnoBeep,6); % Trial end time.
    
    % ************************************************************************************
    
    %%
    
    % For this, we first need the spike count per trial per session.
    
    sps_LC = siteData{3}(LnoBeep,LC_ch); % LC spikes.
    sps_ACC = siteData{3}(LnoBeep,ACC_ch); % ACC spikes.
    
    % **********************************
    nL = [];
    nL_binned = [];
    out_rsc_LCU = [];
    % **********************************
    
    %% Get trial spike counts for LC units.
    
    gti = find(fd>=2100);
    ngt = length(gti);
    nA_rate = nans(1,n_LC);
    for ai = 1:n_LC
        ats = sps_LC(gti,ai);
        for bsi = 1:nbs % Iterate over binsizes.
            binsX = xBin{bsi}; % Get bins for this binsize
            bts = zeros(ngt,min(size((binsX))));
            nts = [];
            for ti = 1:ngt
                ts = ats{ti};
                ts_fix = ts(ts>backTime & ts <=fwdTime); % Drop 1sec of fixation at beginning of stable fixation period.
                nts(ti) = length(ts_fix);
                for bi = 1:nb{bsi} % Iterate over bins.
                    % if binsX(bi,1) >= tst(ti) &  binsX(bi,2) <= tet(ti)
                    bts(ti,bi) = sum((ts_fix > binsX(bi,1) & ts_fix < binsX(bi,2)));
                    % end
                end % Iterate over bins.
                nL_binned{ai,bsi} = bts; % These are the binned spike counts - for 10 binsizes, equally logspaced.
            end % Iterate over trials.
        end % Iterate over bin sizes.
        nL{ai} = nts;
        
        nA_rate(ai) = (1/winSiz)*nanmedian(nts); % 042821.
    end % Iterate over LC units.
    
    lio = 0;
    
    % % numShuf = 10;
    % numShuf = 1;
    
    % pause
    
    %% ACC trials.
    
    numShuf = max([ngt 50]); % Number of shuffles to do.
    numShuf = 10; % Number of shuffles to do.
    ACC_LowHigh = cell(1,n_ACC);
    ntsA = nans(n_ACC,ngt);
    
    for nu = 1:n_ACC % Loop through ACC units.
        gzi = []; % nts = [];
        us = sps_ACC(gti,nu); % Get spikes for good trials only
        for ii = 1:ngt % Iterate over trials.
            spdata = us{ii}; % Get trial spikes.
            ntsA(nu,ii) = length(spdata(spdata>1001 & spdata <= 2100)); % Drop 0.2 sec of fixation at beginning of stable fixation period.
        end
    end
    
    % Get the mean spike count across units, per trial.
    % This will provide the basis for ACC low or high, etc...
    nts = nanmean(ntsA);
    
    % nts_prc = quantile(nts,2); Before 042420
    nts_prc = nanmedian(nts); % 042420
    
    %     lo_sp_ind = find(nts<=nts_prc(1));
    %     mid_sp_ind = find(nts>nts_prc(1) & nts<=nts_prc(2));
    %     hi_sp_ind = find(nts>nts_prc(2));
    lo_sp_ind = find(nts<=nts_prc);
    hi_sp_ind = find(nts>nts_prc);
    
    % If not three quantiles...
    if isempty(intersect(lo_sp_ind,hi_sp_ind))
        % lhi = nans(1,length(nts)); lhi(lo_sp_ind)=0; lhi(mid_sp_ind) = 1; lhi(hi_sp_ind) = 2;
        lhi = nans(1,length(nts)); lhi(lo_sp_ind)=0; lhi(hi_sp_ind) = 1;
        ACC_LowHigh = lhi;
    end
    if ~isempty(intersect(lo_sp_ind,hi_sp_ind))
        % if ~isempty(intersect(lo_sp_ind,hi_sp_ind))
        ACC_LowHigh = [];
    end
    
    nts_u = nts;
    
    %% NOW - calculate LC pairwise spike-count correlations:
    
    
    
    
    if ~isempty(ACC_LowHigh)
        out_rsc_binned = cell(1,nbs);
        
        lo_sp_ind = find(ACC_LowHigh == 0);
        hi_sp_ind = find(ACC_LowHigh == 1);
        
        asi = 1:ngt;
        
        % For this ACC session:
        nu = size(nL_binned,1); bsi = size(nb,2);
        
        for nbi = 1:bsi % Iterate over number of binsizes used.
            % Preallocate, OR, initialize as empty.
            out_rsc = [];
            p_loP = nans(n_ap,100); p_hiP = nans(n_ap,100); p_midP = nans(n_ap,100); p_allP = nans(n_ap,100);
            out_rate = [];
            
            
            for ai = 1:n_ap % For this binsize, iterate over all LC pairs.
                
                a1 = find(LC_pairs(ai,1) == LC_ch);
                a2 = find(LC_pairs(ai,2) == LC_ch);
                
                ap1 = nL_binned{a1,nbi}; % Get spike count for 1st unit in pair.
                ap2 = nL_binned{a2,nbi}; % Get spike count for 2nd unit in pair.
                
                al_lo1 = ap1(lo_sp_ind,:); % LC pair unit 1; LC Spikes from trials corr to low ACC spiking.
                al_lo2 = ap2(lo_sp_ind,:); % LC pair unit 2; LC Spikes from trials corr to low ACC spiking.
                
                al_hi1 = ap1(hi_sp_ind,:); % ACC pair unit 1; Spikes from trials corr to high ACC spiking.
                al_hi2 = ap2(hi_sp_ind,:); % ACC pair unit 2; Spikes from trials corr to high ACC spiking.
                
                % Call "get_pwc" to calculate the pairwise correlations in spike counts.
                
                if (nA_rate(a1) >= 0 & nA_rate(a2) >= 0) % 051621. Criterion is different from ACC because of low LC fr's.
                    rsc_lo =get_pwc(al_lo1(:),al_lo2(:),numLim);
                    rsc_hi =get_pwc(al_hi1(:),al_hi2(:),numLim);
                    rsc_all = get_pwc(ap1(:),ap2(:),numLim);
                    
                    % Create the output matrix:
                    out_rsc = [out_rsc; rsc_lo.RSC rsc_hi.RSC rsc_all.RSC rsc_lo.pval rsc_hi.pval rsc_all.pval];
                end
                
                out_rate{ai} = (1000/bs(nbi))*[nanmean(al_lo1(:)); nanmean(al_lo2(:)); nanmean(al_hi1(:)); nanmean(al_hi2(:)); nanmean(ap1(:)); nanmean(ap2(:))];
                
                if (nA_rate(a1) < 1 | nA_rate(a2) < 1)
                    out_rsc = [out_rsc; nan nan nan nan nan nan];
                end
                
                % **********************************************
                % NO shuffled calculation: % 051621 - Sidd.
                
                % Now do a shuffled calculation: % Added 032017 - Sidd.
                % Modified on 031119 to do 100 x 100 shuffles per pair % - Sidd.
                % lo_r_shuf = []; hi_r_shuf = []; all_r_shuf = [];
                
                p_lo = nans(1,100); p_hi = nans(1,100); p_all = nans(1,100);
                p_loP(ai,:) = p_lo;
                p_hiP(ai,:) = p_hi;
                p_allP(ai,:) = p_all;
                
            end % Loop over pairs
            
            % **********************************
            % CV calculation in archived file if needed.
            % **********************************
            out_rsc_binned{nbi} = {out_rsc out_rate p_loP p_hiP p_allP};
        end % Loop over number of binsizes tested.
        lio = lio+1;
        out_rsc_LCU{lio} = out_rsc_binned;
        
        if lio>0
            dd1 = LC_ch(1)-1;
            pwc_ = {out_rsc_LCU nL nL_binned  ACC_LowHigh nts_u  LC_pairs-dd1};
        end
        
        if lio==0
            pwc_ = [];
        end
        
    end
    
    if isempty(ACC_LowHigh)
        pwc_ = [];
    end
    
end

if n_LC <= 1
    pwc_ = [];
end

