function [pwc_] = getLC_DualArea_rsc_3(siteData,numLim,fi)

% function [pwc_] = getLC_DualArea_rsc2(siteData,numLim,fi)
%
% Calculate pairwise trial spike-count correlations between ACC neurons, conditioned on LC activity.
%
% Origin: 062516 - Sidd.
% History:
% Mod: 032521 - Sidd - new shuffle analysis for reviewer response.
% Mod: 112319 - Sidd - wrote out LC_NZ.
% Mod: 011219 - Now calculate based on spike/no spike rather than zero/median split.
% Mod: 082518 - Sidd - now uses analysis with LC spiking divided into 3 classes - zero spikes and median split on remainder.
% Mod - tweaked code for rsc analysis on binned spike counts... output cells are slightly different now.
% Mod - 060817 - Added analysis for measuring the effect of window size on spike count and on variability (FF and rsc).
% % Mod - 050817 - Wrote out spike rates for mean-matched plots.
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

% How many units in LC?
n_LC = length(LC_ch);
n_ACC = length(ACC_ch);

% How many pairs in ACC?
ACC_pairs = nchoosek(ACC_ch,2);
n_ap = size(ACC_pairs,1);

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
nA = [];
nA_binned = [];
out_rsc_LCU = [];
% **********************************

%% Get trial spike counts for ACC units.

gti = find(fd>=2100);
ngt = length(gti);
nA_rate = nans(1,n_ACC);
for ai = 1:n_ACC
    ats = sps_ACC(gti,ai);
    for bsi = 1:nbs % Iterate over binsizes.
        binsX = xBin{bsi}; % Get bins for this binsize
        bts = zeros(ngt,size(binsX,1));
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
            nA_binned{ai,bsi} = bts; % These are the binned spike counts - for 10 binsizes, equally logspaced.
        end % Iterate over trials.
    end % Iterate over bin sizes.
    nA{ai} = nts;
    %     nA_rate(ai) = (1/winSiz)*nanmean(nts); % Before 042821.
    nA_rate(ai) = (1/winSiz)*nanmedian(nts); % 042821.
end % Iterate over ACC units.

% lio = 0;

% % numShuf = 10;
% numShuf = 1;

% pause

%% LC trials.

LC_LowHigh = cell(1,n_LC);
LC_NZ = cell(1,n_LC);

for nu = 1:n_LC % Loop through LC units.
    lhi = [];
    gzi = []; nts = [];
    us = sps_LC(gti,nu); % Get spikes for good trials only
    for ii = 1:ngt % Iterate over trials.
        spdata = us{ii}; % Get trial spikes.
        nts(ii) = length(spdata(spdata>backTime & spdata <= fwdTime)); % Drop 0.2 sec of fixation at beginning of stable fixation period.
    end
    
    % Get percentiles in case needed later.
    ze_sp_ind = find(nts==0);
    one_sp_ind = find(nts==1);
    two_sp_ind = find(nts==2);
    three_sp_ind = find(nts==3);
    many_sp_ind = find(nts>3);
    nz_sp_ind = find(nts>0);
    
    lhi(ze_sp_ind)=0;
    lhi(one_sp_ind) = 1;
    lhi(two_sp_ind) = 2;
    lhi(three_sp_ind) = 3;
    lhi(many_sp_ind) = 4;
    % lhi(nz_sp_ind) = 9;
    
    LC_LowHigh{nu} = lhi;
    LC_NZ{nu} = nz_sp_ind;
    
    nts_u{nu} = nts;
    
end

% pause

%% NOW - calculate ACC pairwise spike-count correlations:

numShuf = 1; % 032521.
pair_id = [];

if ~isempty(LC_LowHigh)
    for li = 1:n_LC
        out_rsc_binned = cell(1,nbs);
        ze_sp_ind = find(LC_LowHigh{li}==0);
        nz_sp_ind = LC_NZ{li};
        
        % Added 080719 after Neuron review+rejection.
        % LC spike conditions: 1, 2, 3 or >3 LC spikes in each trial.
        one_sp_ind = find(LC_LowHigh{li}==1);
        two_sp_ind = find(LC_LowHigh{li}==2);
        three_sp_ind = find(LC_LowHigh{li}==3);
        many_sp_ind = find(LC_LowHigh{li}==4);
        
        asi = 1:ngt;
        
        % For this LC unit:
        nu = size(nA_binned,1); bsi = size(nb,2);
        
        for nbi = 1:bsi % Iterate over number of binsizes used.
            % Preallocate, OR, initialize as empty.
            out_rsc = [];
            p_zeP = nans(1,n_ap); p_nzP = nans(1,n_ap); p_allP = nans(1,n_ap);
            p_zeP_all = nans(n_ap,numShuf); p_nzP_all = nans(n_ap,numShuf); p_allP_all = nans(n_ap,numShuf);
            out_rate = [];
            out_FF = [];
            
            for ai = 1:n_ap % For this binsize, iterate over all ACC pairs.
                
                if nbi == 1
                    pair_id = [pair_id; li ai];
                end
                
                a1 = find(ACC_pairs(ai,1) == ACC_ch);
                a2 = find(ACC_pairs(ai,2) == ACC_ch);
                
                ap1 = nA_binned{a1,nbi}; % Get spike count for 1st unit in pair.
                ap2 = nA_binned{a2,nbi}; % Get spike count for 2nd unit in pair.
                
                al_ze1 = ap1(ze_sp_ind,:); % ACC pair unit 1; Spikes from ACC trials corr to zero LC spiking.
                al_ze2 = ap2(ze_sp_ind,:); % ACC pair unit 2; Spikes from ACC trials corr to zero LC spiking.
                % rsc_ze = get_pwc(reshape(al_ze1.',1,[]),reshape(al_ze2.',1,[]),numLim);
                
                
                al_nz1 = ap1(nz_sp_ind,:); % ACC pair unit 1; Spikes from ACC trials corr to nonzero LC spiking.
                al_nz2 = ap2(nz_sp_ind,:); % ACC pair unit 2; Spikes from ACC trials corr to nonzero LC spiking.
                % rsc_nz = get_pwc(reshape(al_nz1.',1,[]),reshape(al_nz2.',1,[]),numLim);
                
                al_11 = ap1(one_sp_ind,:); % ACC pair unit 1; Spikes from ACC trials corr to nonzero LC spiking.
                al_12 = ap2(one_sp_ind,:); % ACC pair unit 2; Spikes from ACC trials corr to nonzero LC spiking.
                
                al_21 = ap1(two_sp_ind,:); % ACC pair unit 1; Spikes from ACC trials corr to nonzero LC spiking.
                al_22 = ap2(two_sp_ind,:); % ACC pair unit 2; Spikes from ACC trials corr to nonzero LC spiking.
                
                al_31 = ap1(three_sp_ind,:); % ACC pair unit 1; Spikes from ACC trials corr to nonzero LC spiking.
                al_32 = ap2(three_sp_ind,:); % ACC pair unit 2; Spikes from ACC trials corr to nonzero LC spiking.
                
                al_41 = ap1(many_sp_ind,:); % ACC pair unit 1; Spikes from ACC trials corr to nonzero LC spiking.
                al_42 = ap2(many_sp_ind,:); % ACC pair unit 2; Spikes from ACC trials corr to nonzero LC spiking.
                
                if (nA_rate(a1) >= 1 & nA_rate(a2) >= 1) % 042821
                    
                    rsc_ze = get_pwc(al_ze1(:),al_ze2(:),numLim);
                    rsc_nz = get_pwc(al_nz1(:),al_nz2(:),numLim);
                    rsc_1 = get_pwc(al_11(:),al_12(:),numLim);
                    rsc_2 = get_pwc(al_21(:),al_22(:),numLim);
                    rsc_3 = get_pwc(al_31(:),al_32(:),numLim);
                    rsc_many = get_pwc(al_41(:),al_42(:),numLim);
                    rsc_all = get_pwc(ap1(:),ap2(:),numLim);
                    
                    % Create the output matrix:
                    out_rsc = [out_rsc; rsc_ze.RSC rsc_nz.RSC rsc_all.RSC rsc_ze.pval rsc_nz.pval rsc_all.pval fi li ai rsc_1.RSC rsc_2.RSC rsc_3.RSC rsc_many.RSC rsc_1.pval rsc_2.pval rsc_3.pval rsc_many.pval];
                    %                                 1              2               3              4               5               6      7 8 9        10           11            12                13            14            15           16                 17
                end
                
                out_rate{ai} = (1000/bs(nbi))*[nanmean(al_ze1(:)); nanmean(al_ze2(:)); nanmean(al_nz1(:)); nanmean(al_nz2(:)); nanmean(ap1(:)); nanmean(ap2(:))];
                out_FF{ai} = [nanvar(al_ze1(:))/nanmean(al_ze1(:)); nanvar(al_ze2(:))/nanmean(al_ze2(:)); nanvar(al_nz1(:))/nanmean(al_nz1(:)); nanvar(al_nz2(:))/nanmean(al_nz2(:)); nanvar(ap1(:))/nanmean(ap1(:)); nanvar(ap2(:))/nanmean(ap2(:))];
                
                if (nA_rate(a1) < 1 | nA_rate(a2) < 1)
                    out_rsc = [out_rsc; nan nan nan nan nan nan fi li ai nan nan nan nan nan nan nan nan];
                end
                
                % **********************************************
                % Now do a shuffled calculation: % Added 032017 - Sidd.
                % Modified on 030819 to do 100 x 100 shuffles per pair % - Sidd.
                
                p_ze = nans(1,numShuf);
                p_nz = nans(1,numShuf);
                p_all = nans(1,numShuf);
                
                for psi = 1:numShuf
                    
                    if length(ze_sp_ind)>numLim & length(nz_sp_ind)>numLim & length(asi)>numLim
                        % Zero.
                        ri1 = randi(length(ze_sp_ind),1000,1); ri2 = randi(length(ze_sp_ind),1000,1);
                        drop_ind = intersect(ri1,ri2);
                        ri1(drop_ind) = []; ri2(drop_ind) = [];
                        ri1 = ri1(1:length(ze_sp_ind)); ri2 = ri2(1:length(ze_sp_ind));
                        
                        ze_sp_ind_shuf = ze_sp_ind(ri1); ze_sp_ind_shuf2 = ze_sp_ind(ri2);
                        al_ze1s = ap1(ze_sp_ind_shuf,:); % ACC pair unit 1; Spikes from ACC trials corr to zero LC spiking.
                        al_ze2s = ap2(ze_sp_ind_shuf2,:); % ACC pair unit 2; Spikes from ACC trials corr to zero LC spiking.
                        
                        % Non-zero.
                        ri1 = randi(length(nz_sp_ind),1000,1); ri2 = randi(length(nz_sp_ind),1000,1);
                        drop_ind = intersect(ri1,ri2);
                        ri1(drop_ind) = []; ri2(drop_ind) = [];
                        ri1 = ri1(1:length(nz_sp_ind)); ri2 = ri2(1:length(nz_sp_ind));
                        
                        nz_sp_ind_shuf = nz_sp_ind(ri1); nz_sp_ind_shuf2 = nz_sp_ind(ri2);
                        al_nz1s = ap1(nz_sp_ind_shuf,:); % ACC pair unit 1; Spikes from ACC trials corr to zero LC spiking.
                        al_nz2s = ap2(nz_sp_ind_shuf2,:); % ACC pair unit 2; Spikes from ACC trials corr to zero LC spiking.
                        
                        % All.
                        ri1 = randi(length(asi),1000,1); ri2 = randi(length(asi),1000,1);
                        drop_ind = intersect(ri1,ri2);
                        ri1(drop_ind) = []; ri2(drop_ind) = [];
                        ri1 = ri1(1:length(asi)); ri2 = ri2(1:length(asi));
                        
                        sp_ind_shuf = asi(ri1); sp_ind_shuf2 = asi(ri2);
                        ap1_shuf = ap1(sp_ind_shuf,:); % ACC pair unit 1; Spikes from all trials.
                        ap2_shuf = ap2(sp_ind_shuf2,:); % ACC pair unit 2; Spikes from all trials.
                        
                        % Call "get_pwc" to calculate the pairwise correlations in spike counts.
                        rsc_ze_shuf =get_pwc(al_ze1s(:),al_ze2s(:),numLim);
                        rsc_nz_shuf =get_pwc(al_nz1s(:),al_nz2s(:),numLim);
                        rsc_all_shuf = get_pwc(ap1_shuf(:),ap2_shuf(:),numLim);
                        
                        % Create the output matrix:
                        % lo_r_shuf(jj) = rsc_lo_shuf.RSC;
                        % mid_r_shuf(jj) = rsc_mid_shuf.RSC;
                        % hi_r_shuf(jj) = rsc_hi_shuf.RSC;
                        % all_r_shuf(jj) = rsc_all_shuf.RSC;
                        
                        %                     ze_r_shuf = [ze_r_shuf rsc_ze_shuf.RSC];
                        %                     nz_r_shuf = [nz_r_shuf rsc_nz_shuf.RSC];
                        %                     all_r_shuf = [all_r_shuf rsc_all_shuf.RSC];
                        
                        % Shuffled calculation is done.
                        % *****************************
                        
                        % pause
                        
                        %                     p_ze = [p_ze ze_r_shuf];
                        %                     p_nz = [p_nz nz_r_shuf];
                        %                     p_all = [p_all all_r_shuf];
                        
                        p_ze(psi) = rsc_ze_shuf.RSC;
                        p_nz(psi) = rsc_nz_shuf.RSC;
                        p_all(psi) = rsc_all_shuf.RSC;
                        
                    end
                end % Loop over shuffles for this pair.
                
                p_zeP(ai) = nanmedian(p_ze);
                p_nzP(ai) = nanmedian(p_nz);
                p_allP(ai) = nanmedian(p_all);
                
                p_zeP_all(ai,:) = p_ze;
                p_nzP_all(ai,:) = p_nz;
                p_allP_all(ai,:) = p_all;
                
                
            end % Loop over pairs
            
            % **********************************
            % CV calculation in archived file if needed.
            % **********************************
            out_rsc_binned{nbi} = {out_rsc out_rate p_zeP p_nzP p_allP p_zeP_all p_nzP_all p_allP_all out_FF};
            
            
            
        end % Loop over number of binsizes tested.
        
        %             pause
        
        out_rsc_LCU{li} = out_rsc_binned;
        % end % Check for #trials in each -tile.
    end % LC unit loop.
    dd1 = ACC_ch(1)-1;
    pwc_ = {out_rsc_LCU nA nA_binned LC_LowHigh nts_u ACC_pairs-dd1 LC_NZ};
end

if isempty(LC_LowHigh)
    pwc_ = [];
end

%%

% % Find example.
% figureNumber = 333; num = 3333; wid = 14; hts = [6 6]; cols = {2 1}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,3,4, [12], 'Joshi & Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
% movegui(fig_,[1,1]); % Move figure so it is visible on screen.
%
% xbs = bs;
%
% for jli = 1:n_LC
%     % jili=2
%
%     lhi = LC_LowHigh{jli};
%     low_ind = find(lhi==0);
%     high_ind = find(lhi>0);
%
%     axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
%
%     c1='k';
%
%     ACC_pair_id = ACC_pairs-dd1;
%
%     for pind = 1:length(ACC_pair_id)
%         pind=44
%
%         %% Plot rsc for this pair.
%
%         axes(axs(3)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax); hold on;
%         r_sc_z = []; r_sc_nz = [];
%         for bindex = 1:nbs
%             r_sc_z(bindex) = (out_rsc_LCU{jli}{bindex}{1}(pind,1));
%             r_sc_nz(bindex) = (out_rsc_LCU{jli}{bindex}{1}(pind,2));
%         end
%
%         plot(xbs,r_sc_z,'ms','linewidth',1,'markerfacecolor','m','markeredgecolor','none','markersize',5);
%         plot(xbs,r_sc_nz,'kd','linewidth',1,'markersize',5);
%         plot([50 1050],[0 0],'k--','linewidth',0.5);
%         legend({'LC zero','LC nonzero'},'autoupdate','off','fontsize',8,'location','northwest','fontname','arial'); legend('boxoff');
%         xlabel('binsize','fontsize',10,'fontname','arial');
%         ylabel('rsc','fontsize',10,'fontname','arial');
%         xlim([100 1100]); % ylim([-0.1 0.52]);
%         title(['fileindex=',num2str(fi),' lcindex = ',num2str(jli),' pairindex=',num2str(pind)]);
%
%         %% Spike count scatter plots for this pair - small binsize.
%
%         uind = ACC_pair_id(pind,:)
%
%         axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
%
%         spCount1 = nA_binned{uind(1),1}; spCount1 = spCount1(low_ind,:); spCount1 = spCount1(:);
%         spCount2 = nA_binned{uind(2),1}; spCount2 = spCount2(low_ind,:); spCount2 = spCount2(:);
%         plot(spCount1,spCount2,'ms','markerfacecolor','m','markeredgecolor','none','markersize',5);
%
%         spCount3 = nA_binned{uind(1),1}; spCount3 = spCount3(high_ind,:); spCount3 = spCount3(:);
%         spCount4 = nA_binned{uind(2),1}; spCount4 = spCount4(high_ind,:); spCount4 = spCount4(:);
%         plot(spCount3,spCount4,'kd','markerfacecolor','k','markeredgecolor','none','markersize',5);
%
%         legend('LC zero','LC nonzero','autoupdate','off','fontsize',7,'location','best','fontname','arial'); % legend('boxoff');
%
%         if ~isempty(spCount1) & ~isempty(spCount2)
%             br = fitlm(spCount1,spCount2);
%             fitParams = br.Coefficients.Estimate
%             xv = linspace(min(spCount1),max(spCount1),10);
%             fitY = fitParams(1)+fitParams(2)*xv;
%             plot(xv,fitY,'m-','linewidth',0.5);
%             arrr = corrcoef(spCount1,spCount2);
%             % text(5,8,strcat('rsc ze=',num2str(arrr(1,2))),'color','m','fontsize',9,'fontname','arial');
%             axis([0 13 0 13]);
%             xlabel('spike count, unit 1','fontsize',10,'fontname','arial');
%             ylabel('spike count, unit 2','fontsize',10,'fontname','arial');
%         end
%
%         if ~isempty(spCount1) & ~isempty(spCount2)
%             br = fitlm(spCount3,spCount4);
%             fitParams = br.Coefficients.Estimate
%             xv = linspace(min(spCount3),max(spCount3),10);
%             fitY = fitParams(1)+fitParams(2)*xv;
%             plot(xv,fitY,'k-','linewidth',0.5);
%             arrr = corrcoef(spCount3,spCount4);
%             % text(5,10,strcat('rsc ze=',num2str(arrr(1,2))),'color','k','fontsize',9,'fontname','arial');
%             %         axis([0 11 0 11]);
%             xlabel('spike count, unit 1','fontsize',10,'fontname','arial');
%             ylabel('spike count, unit 2','fontsize',10,'fontname','arial');
%         end
%
%
%         %% Spike count scatter plots for this pair - large binsize.
%
%         axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
%
%         spCount5 = nA_binned{uind(1),5}; spCount5 = spCount5(low_ind,:); spCount5 = spCount5(:);
%         spCount6 = nA_binned{uind(2),5}; spCount6 = spCount6(low_ind,:); spCount6 = spCount6(:);
%         if ~isempty(spCount5) & ~isempty(spCount6)
%             plot(spCount5,spCount6,'ms','markerfacecolor','m','markeredgecolor','none','markersize',5);
%             br = fitlm(spCount5,spCount6);
%             fitParams = br.Coefficients.Estimate
%             xv = linspace(min(spCount5),max(spCount5),10);
%             fitY = fitParams(1)+fitParams(2)*xv;
%             plot(xv,fitY,'m-','linewidth',0.5);
%             arrr = corrcoef(spCount5,spCount6);
%             % text(5,8,strcat('rsc ze=',num2str(arrr(1,2))),'color','m','fontsize',9,'fontname','arial');
%             %         axis([0 13 0 13]);
%             xlabel('spike count, unit 1','fontsize',10,'fontname','arial');
%             ylabel('spike count, unit 2','fontsize',10,'fontname','arial');
%         end
%
%         spCount7 = nA_binned{uind(1),5}; spCount7 = spCount7(high_ind,:); spCount7 = spCount7(:);
%         spCount8 = nA_binned{uind(2),5}; spCount8 = spCount8(high_ind,:); spCount8 = spCount8(:);
%         if ~isempty(spCount7) & ~isempty(spCount8)
%             plot(spCount7,spCount8,'kd','markerfacecolor','k','markeredgecolor','none','markersize',5);
%             br = fitlm(spCount7,spCount8);
%             fitParams = br.Coefficients.Estimate
%             xv = linspace(min(spCount7),max(spCount7),10);
%             fitY = fitParams(1)+fitParams(2)*xv;
%             plot(xv,fitY,'k-','linewidth',0.5);
%             arrr = corrcoef(spCount7,spCount8);
%             % text(5,10,strcat('rsc ze=',num2str(arrr(1,2))),'color','k','fontsize',9,'fontname','arial');
%             %         axis([0 13 0 13]);
%             xlabel('spike count, unit 1','fontsize',10,'fontname','arial');
%             ylabel('spike count, unit 2','fontsize',10,'fontname','arial');
%         end
%
%         pause;
%
%     end
%
% end
%
% close all;
%
% cd C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC;
% save exampleFig3_28_1_44 r_sc_z r_sc_nz spCount1 spCount2 spCount3 spCount4 spCount5 spCount6 spCount7 spCount8;

