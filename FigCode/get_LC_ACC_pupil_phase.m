function [pwc_] = get_LC_ACC_pupil_phase(siteData,nlim,fi,mmnum)
%
% Fixation task.

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

%%

% tBack = 2000;
% TT = (1/2000)*(tBack-(abs(-(tBack-1):(tBack-1))));

tBack = 1000;
TT = (1/1000)*(tBack-(abs(-(tBack-1):(tBack-1))));

%%  Get unit information for this session:
% Note that with we are using channel 4 for LC and channels 9-16 for ACC.

sp_ch = siteData{7}{1};      % All spike channel id's.

n_ACC = 0;
n_LC = 0;

ccg_0_out = [];
ccg_180_out = [];
ccg_90_out = [];
ccg_270_out = [];


% NOTE: Remember that the 2016 LC pupil data set was sorted with MU as the last channel - so drop this for now for Ci and Oz old data!!!
if mmnum >=3
    LC_ch = find(sp_ch>0);   % LC channel id's.
    LC_ch = LC_ch(1:end-1);
    n_LC = length(LC_ch);
    
    if n_LC > 1
        % How many pairs in LC?
        LC_pairs = nchoosek(LC_ch,2);
        n_lp = size(LC_pairs,1);
    end
end


if mmnum <3
    % How many units in LC?
    LC_ch = find(sp_ch<5000);   % LC channel id's.
    ACC_ch = find(sp_ch>8000); % ACC channel id's.
    n_LC = length(LC_ch);
    n_ACC = length(ACC_ch);
    
    if n_LC>1
        % How many pairs in LC?
        LC_pairs = nchoosek(LC_ch,2);
        n_lp = size(LC_pairs,1);
    end
    
    if n_ACC>1
        % How many pairs in ACC?
        ACC_pairs = nchoosek(ACC_ch,2);
        n_ap = size(ACC_pairs,1);
    end
end

%% Get timing information for this session:

fst = siteData{1}(LnoBeep,1); % Fixation start time (wrt fix ON).
fet = siteData{1}(LnoBeep,2); % Fixation end time (wrt FP ON).
fd = siteData{1}(LnoBeep,2); % Fixation duration.

tst = siteData{1}(LnoBeep,5); % Trial start time.
tet = siteData{1}(LnoBeep,6); % Trial end time.

%% Calculate ACC correlations for trials corresponding to zero LC spiking and nonzero LC spiking.

% For this, we first need the spike count per trial per session.

sps_LC = siteData{3}(LnoBeep,LC_ch); % LC spikes.
if mmnum<3
    sps_ACC = siteData{3}(LnoBeep,ACC_ch); % ACC spikes.
end

% **********************************
nA = [];
nA_binned = [];
out_rsc_ACCU = [];
% **********************************

ngt = length(FnoBeep);

%% Pupil phase.

pdat = siteData{2}(LnoBeep,:,4); % trial/sample
evt_dat = siteData{5}; % Pupil events.
gei = evt_dat(:,1); % Trial numbers for events.
% est = evt_dat(:,2); % Event start time.
% eet = evt_dat(:,3); % Event end time.
% estm = evt_dat(:,6); % Event start magnitude.
% eetm = evt_dat(:,7); % Event end magnitude.
edt = evt_dat(:,8); % Time of max slope.
edtm = evt_dat(:,9); % Magnitude of max slope.

evt_out_time = [];
evt_out_sign = [];
gti2 = [];
outphase = [];
out_time_phase = [];

% figureNumber = 1; num = 1; wid = 17.6; hts = [2]; cols = {2 2 2 2 2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 2, 2, [12], 'Joshi & Gold, 2018', true,1,figureNumber); set(axs,'Units','normalized');
% movegui(fig_,[1,1]); % Move figure so it is visible on screen.

plt_index = 0;

for kk = 1:num_noBeep
    if fd(kk)>=1500
        tpdat=pdat(kk,:); gi = find(~isfinite(tpdat)); gi = 1:gi-1; plotP = tpdat(gi);
        ti = FnoBeep(kk); gpi=[]; gpi = find(gei==ti);
        t_edt = edt(gpi); t_edtm = edtm(gpi);
        
        if max(t_edt)<max(gi)
            tti = find(t_edt >=1500);
            if length(tti)>=3
                evt_time = t_edt(tti); evt_time = evt_time(1:3);
                evt_sign = sign(t_edtm(tti)); evt_sign = evt_sign(1:3);
                %                 evt_time = t_edt(tti); evt_time = evt_time(end-2:end);
                %                 evt_sign = sign(t_edtm(tti)); evt_sign = evt_sign(end-2:end);
                %             if sum(evt_sign) ==2 | sum(evt_sign) == -2
                %                 evt_out_time = [evt_out_time; evt_time];
                %                 evt_out_sign = [evt_out_sign; evt_sign];
                %             end
                
                % Test plots.
                tax = evt_time(1):evt_time(3);
                
                % Find peak/trough.
                if evt_sign(1) > 0
                    pii = evt_time(1):evt_time(2); a1 = length(pii);
                    pii1 = find(plotP(pii) == max(plotP(pii))); pii1 = pii1(1);
                    pii = evt_time(2):evt_time(3);
                    pii2 = -1+a1+find(plotP(pii) == min(plotP(pii))); pii2 = pii2(1);
                    % trial_phase = 0:15:255;
                    trial_phase = 0:15:345;
                end
                
                if evt_sign(1) < 0
                    pii = evt_time(1):evt_time(2); a1 = length(pii);
                    pii1 = find(plotP(pii) == min(plotP(pii))); pii1 = pii1(1);
                    pii = evt_time(2):evt_time(3);
                    pii2 = -1+a1+find(plotP(pii) == max(plotP(pii))); pii2 = pii2(1);
                    % trial_phase = [180:15:255 0:15:165];
                    trial_phase = [180:15:345 0:15:165];
                end
                
                
                x1 = evt_time(1):tax(pii1);
                x2 = tax(pii1):evt_time(2);
                x3 = evt_time(2):tax(pii2);
                x4 = tax(pii2):evt_time(3);
                
                cyc_lim = 50;
                if length(x1) > cyc_lim & length(x2) > cyc_lim & length(x3) > cyc_lim & length(x4) > cyc_lim
                    
                    gti2 = [gti2; kk];
                    
                    xx1 = linspace(evt_time(1),tax(pii1),6);
                    xx2 = linspace(1+tax(pii1),evt_time(2),6);
                    xx3 = linspace(1+evt_time(2),tax(pii2),6);
                    xx4 = linspace(1+tax(pii2),evt_time(3),6);
                    
                    interp_y1 = spline(x1,plotP(x1),xx1);
                    interp_y2 = spline(x2,plotP(x2),xx2);
                    interp_y3 = spline(x3,plotP(x3),xx3);
                    interp_y4 = spline(x4,plotP(x4),xx4);
                    
                    %                 % ***********************
                    %                 % Plot things to test code - COMMENT OUT WHEN RUNNING AS FUNCTION!!!.
                    %                 plt_index = plt_index+1;
                    %
                    %                 axes(axs(plt_index)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
                    %                 plot(tax,plotP(tax),'k-');
                    %                 plot(evt_time,plotP(evt_time),'bx');
                    %
                    %                 plot(tax(pii1),plotP(tax(pii1)),'r*');
                    %                 plot(tax(pii2),plotP(tax(pii2)),'r*');
                    %
                    %                 ylabel('pupil size');
                    %                 xlabel('time (msec)');
                    %                 title(num2str(kk));
                    %                 xlim([tax(1) tax(end)]);
                    %
                    %                 plt_index = plt_index+1;
                    %                 axes(axs(plt_index)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
                    %                 plot(x1,plotP(x1),'g-');
                    %                 plot(x2,plotP(x2),'m-');
                    %                 plot(x3,plotP(x3),'c-');
                    %                 plot(x4,plotP(x4),'r-');
                    %
                    %                 plot(xx1,interp_y1,'ko');
                    %                 plot(xx2,interp_y2,'ko');
                    %                 plot(xx3,interp_y3,'ko');
                    %                 plot(xx4,interp_y4,'ko');
                    %
                    %                 ylabel('pupil size');
                    %                 xlabel('time (msec)');
                    %                 xlim([tax(1) tax(end)]);
                    %
                    %                 % Plotted things.
                    %                 % ***********************
                    
                    outphase = [outphase; trial_phase];
                    % out_time_phase = [out_time_phase; xx1 xx2(2:end) xx3(2:end)];
                    out_time_phase = [out_time_phase; xx1 xx2 xx3 xx4];
                    %                     out_time_phase = [out_time_phase; xx1 xx2 xx3];
                    
                end
                
                %             pause;
                %             clf;
                
            end
        end
    end
end

%%

ngt = length(gti2);


%% Get trial spike counts for ACC units.

backTime = 1001; fwdTime = 2100; % 032320b.
winSiz = length(backTime:fwdTime)/1000;

nA_rate = nans(1,n_ACC);
for ai = 1:n_ACC
    ats = sps_ACC(gti2,ai);
    nts = [];
    for ti = 1:ngt
        ts = ats{ti};
        ts_fix = ts(ts>1001 & ts <=2100); % Drop 1sec of fixation at beginning of stable fixation period.
        nts(ti) = length(ts_fix);
    end % Iterate over trials.
    nA{ai} = nts;
    nA_rate(ai) = (1/winSiz)*nanmedian(nts); % 042821.
end % Iterate over ACC units.

%% LC trials.

for ai = 1:n_LC
    ats = sps_LC(gti2,ai);
    nts = [];
    for ti = 1:ngt
        ts = ats{ti};
        ts_fix = ts(ts>1001 & ts <=2100); % Drop 1sec of fixation at beginning of stable fixation period.
        nts(ti) = length(ts_fix);
    end % Iterate over trials.
    nL{ai} = nts;
end % Iterate over ACC units.

%% LC. Calculate PETH's re: pupil phase.

% backTime = -1000; fwdTime = 500; % We're looking at the 1sec to 2.1sec window of stable fixation.
% tmin   = backTime; tmax   = fwdTime; xBin = []; nb = []; tsize  = 200; tstep  = 100;
% xBin = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)']; nb = size(xBin,1);

backTime = -1000; fwdTime = 500; % We're looking at the 1sec to 2.1sec window of stable fixation.
tmin   = backTime; tmax   = fwdTime; xBin = []; nb = []; tsize  = 500; tstep  = 100;
xBin = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)']; nb = size(xBin,1);

pwc_all_lc = [];
pwc_all_acc = [];
unit_dat = [];
unit_dat_Shuf = [];
unit_dat_A = [];
unit_dat_A_Shuf = [];
all_LC = [];
all_LC_bs = []
all_ACC = [];
all_ACC_bs = []
pspike_LC_out = [];
pspike_ACC_out = [];
rsc_LC_out = [];
rsc_ACC_out = [];

ngt;

if ngt > 0
    
    for nu = 1:n_LC % Loop through LC units.
        us = sps_LC(gti2,nu); % Get spikes for good trials only.
        trial_dat = nans(length(trial_phase),size(xBin,1),ngt);
        for ti = 1:ngt % Iterate over trials.
            spdata = us{ti}; % Get trial spikes.
            t_phase = out_time_phase(ti,:);
            p_phase = outphase(ti,:);
            
            bts = zeros(length(p_phase),size(xBin,1));
            for jj = 1:length(t_phase)
                sp_phase = spdata-t_phase(jj);
                % Create PETH.
                for bi = 1:nb % Iterate over bins.
                    bts(jj,bi) = (1000/tsize)*sum((sp_phase > xBin(bi,1) & sp_phase < xBin(bi,2)));
                end
                % If selected pupil period starts with constriction, sort the phases so it starts with zero.
                if p_phase(1) ~=0
                    bts = [bts(13:end,:); bts(1:12,:)];
                end
                %                 pause
            end % Phases for this trial.
            
            % bts = reshape(zscore(bts(:)),size(bts));
            % bts = zscore(bts')';
            trial_dat(:,:,ti) = bts;
            
        end
        unit_dat{nu} = nanmean(trial_dat,3);
        all_LC{nu} = reshape(trial_dat(:)>0,size(trial_dat));
        all_LC_bs{nu} = trial_dat;
    end
    
    
    % % Plot things.
    % pax = 0:15:270;
    % tax = nanmean(xBin,2);
    % plot_matrix(unit_dat{1}',tax,pax,'n');
    
    
    %% LC. Calculate PETH's re: pupil phase - SHUFFLED!
    
    %     for nu = 1:n_LC % Loop through LC units.
    %         us = sps_LC(gti2,nu); % Get spikes for good trials only.
    %         trial_dat = nans(length(trial_phase),size(xBin,1),ngt);
    %         for rti = 1:500 % Iterate over number of shuffled trials.
    %             ti1 = randi(ngt);
    %             ti2 = randi(ngt);
    %             spdata = us{ti1}; % Get trial spikes.
    %             t_phase = out_time_phase(ti2,:);
    %             p_phase = outphase(ti2,:);
    %
    %             bts = zeros(length(t_phase),size(xBin,1));
    %             for jj = 1:length(t_phase)
    %                 sp_phase = spdata-t_phase(jj);
    %                 % Create PETH.
    %                 for bi = 1:nb % Iterate over bins.
    %                     bts(jj,bi) = (1000/tsize)*sum((sp_phase > xBin(bi,1) & sp_phase < xBin(bi,2)));
    %                 end
    %                 % If selected pupil period starts with constriction, sort the phases so it starts with zero.
    %                 if p_phase(1) ~=0
    %                     bts = [bts(13:end,:); bts(1:12,:)];
    %                 end
    %             end % Phases for this trial.
    %
    %             % bts = reshape(zscore(bts(:)),size(bts));
    %             % bts = zscore(bts')';
    %             trial_dat(:,:,rti) = bts;
    %
    %         end
    %         unit_dat_Shuf{nu} = nanmean(trial_dat,3);
    %     end
    
    
    hist_edges = -500:500;
    % hist_edges = -1000:1000;
    
    if mmnum<3
        %% ACC. Calculate PETH's re: pupil phase.
        
        for nu = 1:n_ACC % Loop through ACC units.
            bs_ACC_0 = nans(ngt,length(hist_edges)-1);
            bs_ACC_90 = nans(ngt,length(hist_edges)-1);
            bs_ACC_180 = nans(ngt,length(hist_edges)-1);
            bs_ACC_270 = nans(ngt,length(hist_edges)-1);
            us = sps_ACC(gti2,nu); % Get spikes for good trials only.
            trial_dat = nans(length(trial_phase),size(xBin,1),ngt);
            for ti = 1:ngt % Iterate over trials.
                spdata = us{ti}; % Get trial spikes.
                t_phase = out_time_phase(ti,:);
                p_phase = outphase(ti,:);
                
                bts = zeros(length(t_phase),size(xBin,1));
                for jj = 1:length(t_phase)
                    sp_phase = spdata-t_phase(jj);
                    % Create PETH.
                    for bi = 1:nb % Iterate over bins.
                        bts(jj,bi) = sum((sp_phase > xBin(bi,1) & sp_phase < xBin(bi,2)));
                    end
                    % If selected pupil period starts with constriction, sort the phases so it starts with zero.
                    if p_phase(1) ~=0
                        bts = [bts(13:end,:); bts(1:12,:)];
                    end
                    % end % Phases for this trial.
                    
                    % For CCG's, get 1msec binned spikes for zero and 180 phases only.
                    if p_phase(jj) == 0
                        sp_phase(sp_phase<-1000)=[];
                        sp_phase(sp_phase>1000)=[];
                        bs_ACC_0(ti,:) = histcounts(sp_phase,hist_edges);
                    end
                    if p_phase(jj) == 90
                        sp_phase(sp_phase<-1000)=[];
                        sp_phase(sp_phase>1000)=[];
                        bs_ACC_90(ti,:) = histcounts(sp_phase,hist_edges);
                    end
                    if p_phase(jj) == 180
                        sp_phase(sp_phase<-1000)=[];
                        sp_phase(sp_phase>1000)=[];
                        bs_ACC_180(ti,:) = histcounts(sp_phase,hist_edges);
                    end
                    if p_phase(jj) == 270
                        sp_phase(sp_phase<-1000)=[];
                        sp_phase(sp_phase>1000)=[];
                        bs_ACC_270(ti,:) = histcounts(sp_phase,hist_edges);
                    end
                    
                end % Phases for this trial.
                
                % pause
                
                % bts = reshape(zscore(bts(:)),size(bts));
                % bts = zscore(bts);
                trial_dat(:,:,ti) = bts;
                
            end
            bs_0{nu} = bs_ACC_0;
            bs_90{nu} = bs_ACC_90;
            bs_180{nu} = bs_ACC_180;
            bs_270{nu} = bs_ACC_270;
            unit_dat_A{nu} = nanmean(trial_dat,3);
            all_ACC{nu} = reshape(trial_dat(:)>0,size(trial_dat));
            all_ACC_bs{nu} = trial_dat;
        end
        
        % % Plot things.
        % pax = 0:15:270;
        % tax = nanmean(xBin,2);
        % plot_matrix(unit_dat{1}',tax,pax,'n');
        
        
        %% ACC. Calculate PETH's re: pupil phase - SHUFFLED!
        
        %         for nu = 1:n_ACC % Loop through ACC units.
        %             us = sps_ACC(gti2,nu); % Get spikes for good trials only.
        %             trial_dat = nans(length(trial_phase),size(xBin,1),ngt);
        %             for rti = 1:500 % Iterate over number of shuffled trials.
        %                 ti1 = randi(ngt);
        %                 ti2 = randi(ngt);
        %                 spdata = us{ti1}; % Get trial spikes.
        %                 t_phase = out_time_phase(ti2,:);
        %                 p_phase = outphase(ti2,:);
        %
        %                 bts = zeros(length(t_phase),size(xBin,1));
        %                 for jj = 1:length(t_phase)
        %                     sp_phase = spdata-t_phase(jj);
        %                     % Create PETH.
        %                     for bi = 1:nb % Iterate over bins.
        %                         bts(jj,bi) = sum((sp_phase > xBin(bi,1) & sp_phase < xBin(bi,2)));
        %                     end
        %                     % If selected pupil period starts with constriction, sort the phases so it starts with zero.
        %                     if p_phase(1) ~=0
        %                         bts = [bts(13:end,:); bts(1:12,:)];
        %                     end
        %                 end % Phases for this trial.
        %
        %                 % bts = reshape(zscore(bts(:)),size(bts));
        %                 % bts = zscore(bts);
        %                 trial_dat(:,:,rti) = bts;
        %
        %             end
        %             unit_dat_A_Shuf{nu} = nanmean(trial_dat,3);
        %         end
        
        
    end
    
    
    %% Probability of coinicident spiking in these bins.
    
    
    if n_LC>1
        rsc_LC_out = nans(size(outphase,2),nb,n_lp);
        pwc_all_lc = [];
        for pin = 1:n_lp
            pwc_pb = [];
            p1 = all_LC{find(LC_pairs(pin,1) == LC_ch)};
            p2 = all_LC{find(LC_pairs(pin,2) == LC_ch)};
            pspike_LC = zeros(size(all_LC{1}(:,:,1)));
            for ttii = 1:ngt
                pspike_LC = pspike_LC + p1(:,:,ttii).*p2(:,:,ttii);
            end
            
            p1 = find(LC_pairs(pin,1) == LC_ch);
            p2 = find(LC_pairs(pin,2) == LC_ch);
            
            pwc_all_lc_dat = get_pwc(nL{p1},nL{p2},10);
            pwc_all_lc(pin) = pwc_all_lc_dat.RSC;
            
            for phase_i = 1:size(outphase,2)
                for time_i = 1:nb
                    sp1 = squeeze(all_LC_bs{p1}(phase_i,time_i,:)); % pause
                    sp2 = squeeze(all_LC_bs{p2}(phase_i,time_i,:));
                    rsc_pb = get_pwc(sp1,sp2,10);
                    pwc_pb(phase_i,time_i) = rsc_pb.RSC;
                end
            end
            pspike_LC_out{pin} = pspike_LC/ngt;
            rsc_LC_out(:,:,pin) = pwc_pb;
        end
    end
    
    % pause
    
    if n_ACC > 1
        rsc_ACC_out = nans(size(outphase,2),nb,n_ap);
        pwc_all_acc = [];
        if mmnum<3
            for pin = 1:n_ap
                pwc_pb = [];
                p1 = all_ACC{find(ACC_pairs(pin,1) == ACC_ch)};
                p2 = all_ACC{find(ACC_pairs(pin,2) == ACC_ch)};
                pspike_ACC = zeros(size(all_ACC{1}(:,:,1)));
                
                for ttii = 1:ngt
                    pspike_ACC = pspike_ACC + p1(:,:,ttii).*p2(:,:,ttii);
                end
                
                p1 = find(ACC_pairs(pin,1) == ACC_ch);
                p2 = find(ACC_pairs(pin,2) == ACC_ch);
                
                if (nA_rate(p1) >= 1 & nA_rate(p2) >= 1) % 042821
                    
                    pwc_all_acc_dat = get_pwc(nA{p1},nA{p2},10);
                    pwc_all_acc(pin) = pwc_all_acc_dat.RSC;
                    
                    for phase_i = 1:size(outphase,2)
                        for time_i = 1:nb
                            sp1 = squeeze(all_ACC_bs{p1}(phase_i,time_i,:));
                            sp2 = squeeze(all_ACC_bs{p2}(phase_i,time_i,:));
                            rsc_pb = get_pwc(sp1,sp2,10);
                            pwc_pb(phase_i,time_i) = rsc_pb.RSC;
                        end
                    end
                    
                    pspike_ACC_out{pin} = pspike_ACC/ngt;
                    rsc_ACC_out(:,:,pin) = pwc_pb;
                end
                
                if (nA_rate(p1) < 1 | nA_rate(p2) < 1)
                    pspike_ACC_out{pin} = nan;
                    rsc_ACC_out(:,:,pin) = nan;
                end
                
            end % Pairs loop.
        end % Mnky check.
    end % Pair check.
    
end

%% Write stuff out.

pwc_ = {unit_dat unit_dat_Shuf unit_dat_A unit_dat_A_Shuf pspike_LC_out pspike_ACC_out rsc_LC_out rsc_ACC_out pwc_all_lc pwc_all_acc};

% sumP = zeros(size(all_ACC{1}(:,:,1)));
% for ki = 1:n_ap
%     sumP = sumP+pspike_ACC_out{kk};
% end
% probA = sumP/n_ap;
% imagesc(probA)




