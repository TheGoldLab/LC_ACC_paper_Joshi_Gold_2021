% Fig2_LC_ACC_Basics_ShuffOnly.m
%
% Shuffle analysis only - set # shuffles appropriately in "getLC_DualArea_rsc_3.m" if reanalyzing from scratch.
% Run on a fast computer with plenty of memory.
%
% *****************************
% HISTORY.
% *****************************
% Mod from: Fig4_LC_ACC_Dual_Recording_Firing_Rate.m
% 012820: Changed marker types and added bootstrap analysis for errors.
% 120219: Modified to add plot of variance of ACC  spike counts; also, create LC conditioned ACC figure.
% 080719: Modified to use "getLC_DualArea_rsc_2" - to address reviewer concern
% Instead of LC zero/nonzero split, we will try dividing trials based on spike counts of 1,2,3,4.... see if that is possible.
%
% 011219: Modified - now split trials based on spike/no spike rather than zero/median split.
% 011819: Finalizing figure layouts. Ha ha. Finalizing....
% 121218. Modified from: LC_Dual_Recording_rsc.m and ACC_Dual_Recording_rsc.m
%
% Origin: Fig2_LC_ACC_Basics; 062116 - Sidd
%
% LC spiking divided into 2 classes: zero and non zero spike trials.
% ACC spiking divided into 2 classes: low (<25th prctile of spike rates across trials per session) and high spike trials.
%
% *****************************
% INFO.
% *****************************
% "rsc" is the term used by Bair, Kohn et al for pairwise spike-count (sc) correlations (r).
% Dual area recording - Analysis script; FIXATION task.
% Calls "getLC_DualArea_rsc" and "getACC_DualArea_rsc" to do the heavy lifting.
%
% Calculates pairwise correlations (pwc) between ACC neurons.
% Plots distributions of spike rates in ACC, conditioned on zero, low or high LC counts.
% Saves pwc conditioned on LC/ACC low vs high spiking regimes.
% Plots variability of trial spike counts and trial ISI's using Fano Factor (CV dropped for now).

%% Setup stuff.

% Clear everything?
clear; clear all;

monks = {'Sprout','Cicero'}; % Add mnks as needed.
sites = {'LC_ACC_Fixation'}; % Only LC for now.
base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.
nMonks = length(monks); nSites = length(sites);

% Don't reanalyze data?
reanalyzeLC = false;
reanalyzeACC = false;
% Reanalyze data?
% reanalyzeLC = true;
% reanalyzeACC = true;

%% Analysis Section: Calculate ACC pairwise spike count correlations (rsc) and ACC spiking statistics conditioned on LC spiking.

% Loops through monkeys and sessions.
% Uses function "getLC_DualArea_rsc.m" to do the heavy lifting.
% Saves results for plotting later.

%% LC linked ACC measurements.

min_trials = 5;
% min_trials = 10; % 021020.

if reanalyzeLC
    mnkRes = [];
    for mm = 1:nMonks % Loop over monkeys (wheeeee....!).
        sitRes = [];
        for ss = 1:nSites  % Loop over brainstem sites.
            inDir= strcat([base_dir,'\',monks{mm},'\',sites{ss},'\clean']); % Create dir name for input (clean) files.
            cd(inDir);  dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
            fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
            rsc_dat = [];
            if ~isempty(fnames)
                nf = length(fnames);
                for ff = 1:nf
                    % ff = 15% For troubleshooting.
                    load(fnames{ff});
                    % pause % For troubleshooting.
                    rsc_dat{ff} = getLC_DualArea_rsc_3(siteData,min_trials,ff); % Added analyses for window size effects.
                    fprintf('File %d%sof%s%d\n', ff,' ',' ',nf);
                    % disp(sprintf('File %d%sof%s%d', ff,' ',' ',nf)); % Older versions of MATLAB.
                end % Sessions.
            end % Check for empty file list.
            sitRes{ss} = rsc_dat;
        end % Sites.
        mnkRes{mm} = sitRes;
    end % Mnks.
    
    clear sitRes;
    
    % Save analysis:
    % savedir = strcat([base_dir,'\Results\Results_2019\LC_ACC']); cd(savedir);
    % savedir = strcat([base_dir,'\Results\Results_2020\LC_ACC']); cd(savedir);
    savedir = strcat([base_dir,'\Results\Results_2021\LC_ACC']); cd(savedir);
    
    % save LC_ACC_Results_rsc_032320 mnkRes; % With new way of dividing trials into zero, mid, high.
    % save LC_ACC_Results_rsc_042020 mnkRes; % With new way of dividing trials into zero, mid, high.
    % save LC_ACC_Results_rsc_042320 mnkRes; % With new way of dividing trials into zero, mid, high.
    % save LC_ACC_Results_rsc_071820b mnkRes; % With new way of dividing trials into zero, mid, high.
    % save LC_ACC_Results_rsc_100120 mnkRes; % With new way of dividing trials into zero, mid, high.
    save LC_ACC_Results_rsc_032521 mnkRes; % With new way of dividing trials into zero, mid, high.
    
    clear mnkRes;
end % If reanalyze LC loop.

%% Plotting Section - Load results file created above, calculate summary statistics and plot stuff.

% Don't organize results?
orgYesNo = false;
% Plot results?
orgYesNo = true;

if orgYesNo
    
    cd C:\Sidd\PostDoc2\SidCode\General;
    load('sidd_colors');
    c1 = cmap_M(110,:);
    c2 = cmap_M(170,:);
    
    %% LC linked ACC results.
    
    clear; clear all;
    
    base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC'; % Base directory for brain area.
    savedir = base_dir;
    % inFile = 'LC_ACC_Results_rsc_032521';
    % inFile = 'LC_ACC_Results_rsc_041321';
    inFile = 'LC_ACC_Results_rsc_050321'; % Ran this on faster laptop, therefore not in run section above (Line 100).
    
    cd(savedir); load(inFile);
    
    monks = {'Sprout', 'Cicero'}; % Add mnks as needed...
    sites = {'LC_ACC_Fixation'}; % For now LC only.
    nMonks = length(monks); nSites = length(sites);
    
    nshuft = 1000; % Number of shuffles used.
    nbin = 5;
    % **********************************************
    
    for mm = 1:nMonks % Loop over mnks....
        allShufDatZ = nans(nbin,10000,nshuft);
        allShufDatNZ = nans(nbin,10000,nshuft);
        allShufDatA = nans(nbin,10000,nshuft);
        mnkDat = mnkRes{mm};
        
        si=1; ei = 0;
        
        for ss = 1:nSites % Loop over brainstem sites.... LC and/or IC.
            
            sitDat = mnkDat{ss};
            nnuu = 0; nnAA = 0; % Housekeeping - how many LC and ACC units did we have?
            
            if ~isempty(sitDat)
                for ap = 1:length(sitDat) % Loop over sessions.
                    dat = sitDat{ap}; % Get session data.
                    if ~isempty(dat)
                        allZ = []; allH = []; allA = [];
                        nnAA = nnAA + length(dat{2});
                        
                        % *********************************************
                        % Get pair indices.
                        pair_ind = dat{6};
                        numPairs = size(pair_ind,1);
                        for sp = 1:length(dat{1}) % Loop over LC units in this session (can be 1 or more).
                            
                            LC_lohi = dat{4}{sp};
                            LC_nz = dat{7}{sp}; % Added: 112319.
                            hii = LC_nz;
                            zei = find(LC_lohi==0);
                            
                            if ~isempty(zei) && ~isempty(hii)
                                
                                % Processed data structure:
                                % Cell1:
                                % out_rsc_binned{nbi} = {out_rsc out_rate p_lo p_hi p_all lo_sp_ind hi_sp_ind nts LC_LowHigh};
                                % 1. out_rsc.
                                % 2. out_rate.
                                % 3. p_lo
                                % 4. p_hi
                                % 5. p_all
                                % 6. lo_sp_ind
                                % 7. hi_sp_ind
                                % 8. nts
                                % 9. LC_LowHigh
                                
                                % Cell2: nA
                                
                                % Cell3: nA_binned
                                
                                nnuu = nnuu + 1;
                                dat0 = dat{1}{sp};
                                
                                dat1_LC_spikes = dat{5}{sp}; % LC trial spike counts.
                                LC_sp_hi = dat1_LC_spikes; % 100219; get all spike counts for distribution.
                                
                                nbin = length(dat0);
                                
                                si = 1+ei;
                                ei = si+numPairs-1;
                                for bhi= 1:nbin % Iterate over bins.
                                    
                                    dat1_pze = dat0{bhi}{6}; % Zero LC. Pairwise spike-count correlation results - shuffle tests.
                                    dat1_phi = dat0{bhi}{7}; % High LC. Pairwise spike-count correlation results - shuffle tests.
                                    dat1_pall = dat0{bhi}{8}; % All trials. Pairwise spike-count correlation results - shuffle tests.
                                    
                                    allShufDatZ(bhi,si:ei,:) = dat1_pze;
                                    allShufDatNZ(bhi,si:ei,:) = dat1_phi;
                                    allShufDatA(bhi,si:ei,:) = dat1_pall;
                                    
                                end % Bin loop.
                            end % Check if LC zero and nonzero indices are real.
                        end % LC units.
                    end % Data struct empty or not?
                end % Sessions loop.
            end % Site data struct empty?
        end % Sites.
        
        allShufDatZ(:,ei+1:end,:) = [];
        allShufDatNZ(:,ei+1:end,:) = [];
        allShufDatA(:,ei+1:end,:) = [];
        
        mnkZ{mm} = allShufDatZ;
        mnkNZ{mm} = allShufDatNZ;
        mnkA{mm} = allShufDatA;
        
    end % Monks.
end % ifOrgYesNo.

%% 050421.

for mm = 1:2
    ddZ = mnkZ{mm};
    ddNZ = mnkNZ{mm};
    ddA = mnkA{mm};
    
    all_pvals = nans(3,nbin,nshuft);
    diff_all_rsc = nans(3,nbin,nshuft);
    sigi = nans(3,nbin,nshuft);
    median_diff_rsc = nans(3,nbin);
    median_diff_rsc_iqr = nans(3,nbin,2);
    median_diff_rsc_sig1 = nans(1,5);
    median_diff_rsc_sig2 = nans(1,5);
    median_diff_rsc_sig3 = nans(1,5);
    
    for bi=1:nbin
        for shufInd = 1:nshuft
            bdZ = squeeze(ddZ(bi,:,shufInd)); % LC zero.
            bdNZ = squeeze(ddNZ(bi,:,shufInd)); % LC non-zero.
            bdA = squeeze(ddA(bi,:,shufInd)); % Not conditioned on LC.
            
            terc_vals = quantile(bdA,2); % Get terciles from unconditioned ACC rsc's.
            t1 = find(bdA<=terc_vals(1)); t2 = find(bdA>terc_vals(1) & bdA<terc_vals(2)); t3 = find(bdA>=terc_vals(2));
            diffVal1 = bdNZ(t1)-bdZ(t1); diffVal2 = bdNZ(t2)-bdZ(t2); diffVal3 = bdNZ(t3)-bdZ(t3);
            diffVal1 = diffVal1(isfinite(diffVal1));
            diffVal2 = diffVal2(isfinite(diffVal2));
            diffVal3 = diffVal3(isfinite(diffVal3));
            
            all_pvals(1,bi,shufInd) = signrank(diffVal1); diff_all_rsc(1,bi,shufInd) = nanmedian(diffVal1);
            all_pvals(2,bi,shufInd) = signrank(diffVal2); diff_all_rsc(2,bi,shufInd) = nanmedian(diffVal2);
            all_pvals(3,bi,shufInd) = signrank(diffVal3); diff_all_rsc(3,bi,shufInd) = nanmedian(diffVal3);
            
        end
        
        sigi(:,bi,:) = all_pvals(:,bi,:) < 0.05;
        
        median_diff_rsc(:,bi) = nanmedian(diff_all_rsc(:,bi,:),3);
        median_diff_rsc_sig1(:,bi) = signrank(squeeze(diff_all_rsc(1,bi,:)));
        median_diff_rsc_sig2(:,bi) = signrank(squeeze(diff_all_rsc(2,bi,:)));
        median_diff_rsc_sig3(:,bi) = signrank(squeeze(diff_all_rsc(3,bi,:)));
        
        [ci_1 ci_m1]= bootci(2000,{@median,diff_all_rsc(1,bi,:)},'alpha',0.05,'type','per');
        [ci_2 ci_m2]= bootci(2000,{@median,diff_all_rsc(2,bi,:)},'alpha',0.05,'type','per');
        [ci_3 ci_m3]= bootci(2000,{@median,diff_all_rsc(3,bi,:)},'alpha',0.05,'type','per');
        
        median_diff_rsc_iqr(1,bi,:) = ci_1;
        median_diff_rsc_iqr(2,bi,:) = ci_2;
        median_diff_rsc_iqr(3,bi,:) = ci_3;
        
    end
    mnkSig{mm} = sigi;
    mnk_rsc{mm} = median_diff_rsc;
    mnk_rsc_iqr{mm} = median_diff_rsc_iqr;
    mnk_rsc_sig{mm} = {median_diff_rsc_sig1 median_diff_rsc_sig2 median_diff_rsc_sig3};
end

% Now calculate the probability of getting significant correlations in adjacent bins.
sig_propA = 0; sig_propB = 0;
for si = 1:nshuft
    sig_propA = sig_propA + sum(sum(mnkSig{1}(:,:,si)));
    sig_propB = sig_propB + sum(sum(mnkSig{2}(:,:,si)));
end
pA1 = (1/(15*1000))*(sig_propA);
pB1 = (1/(15*1000))*(sig_propB);

sig_propA = 0; sig_propB = 0;
for ti = 1:3
    for si = 1:nshuft
        sig_propA = sig_propA + (mnkSig{1}(ti,1,si) * mnkSig{1}(ti,2,si)) + (mnkSig{1}(ti,2,si) * mnkSig{1}(ti,3,si)) + (mnkSig{1}(ti,3,si) * mnkSig{1}(ti,4,si)) + (mnkSig{1}(ti,4,si) * mnkSig{1}(ti,5,si));
        sig_propB = sig_propB + (mnkSig{2}(ti,1,si) * mnkSig{2}(ti,2,si)) + (mnkSig{2}(ti,2,si) * mnkSig{2}(ti,3,si)) + (mnkSig{2}(ti,3,si) * mnkSig{2}(ti,4,si)) + (mnkSig{2}(ti,4,si) * mnkSig{2}(ti,5,si));
    end
end
pA2 = (1/(15*1000))*(sig_propA);
pB2 = (1/(15*1000))*(sig_propB);

sig_propA = 0; sig_propB = 0;
for ti = 1:3
    for si = 1:nshuft
        sig_propA = sig_propA + (mnkSig{1}(ti,1,si) * mnkSig{1}(ti,2,si) * mnkSig{1}(ti,3,si)) + (mnkSig{1}(ti,2,si) * mnkSig{1}(ti,3,si) * mnkSig{1}(ti,4,si)) + (mnkSig{1}(ti,3,si) * mnkSig{1}(ti,4,si) * mnkSig{1}(ti,5,si));
        sig_propB = sig_propB + (mnkSig{2}(ti,1,si) * mnkSig{2}(ti,2,si) * mnkSig{2}(ti,3,si)) + (mnkSig{2}(ti,2,si) * mnkSig{2}(ti,3,si) * mnkSig{2}(ti,4,si)) + (mnkSig{2}(ti,3,si) * mnkSig{2}(ti,4,si) * mnkSig{2}(ti,5,si));
    end
end
pA3 = (1/(15*1000))*(sig_propA);
pB3 = (1/(15*1000))*(sig_propB);

sig_propA = 0; sig_propB = 0;
for ti = 1:3
    for si = 1:nshuft
        sig_propA = sig_propA + (mnkSig{1}(ti,1,si) * mnkSig{1}(ti,2,si) * mnkSig{1}(ti,3,si) * mnkSig{1}(ti,4,si)) + (mnkSig{1}(ti,2,si) * mnkSig{1}(ti,3,si) * mnkSig{1}(ti,4,si) * mnkSig{1}(ti,5,si));
        sig_propB = sig_propB + (mnkSig{2}(ti,1,si) * mnkSig{2}(ti,2,si) * mnkSig{2}(ti,3,si) * mnkSig{2}(ti,4,si)) + (mnkSig{2}(ti,2,si) * mnkSig{2}(ti,3,si) * mnkSig{2}(ti,4,si) * mnkSig{2}(ti,5,si));
    end
end
pA4 = (1/(15*1000))*(sig_propA);
pB4 = (1/(15*1000))*(sig_propB);

sig_propA = 0; sig_propB = 0;
for ti = 1:3
    for si = 1:nshuft
        sig_propA = sig_propA + (mnkSig{1}(ti,1,si) * mnkSig{1}(ti,2,si) * mnkSig{1}(ti,3,si) * mnkSig{1}(ti,4,si) * mnkSig{1}(ti,5,si));
        sig_propB = sig_propB + (mnkSig{2}(ti,1,si) * mnkSig{2}(ti,2,si) * mnkSig{2}(ti,3,si) * mnkSig{2}(ti,4,si) * mnkSig{2}(ti,5,si));
    end
end
pA5 = (1/(15*1000))*(sig_propA);
pB5 = (1/(15*1000))*(sig_propB);

pa_All_s = [[pA1 pB1]; [pA2 pB2]; [pA3 pB3]; [pA4 pB4]; [pA5 pB5]];

%% Setup figure:
figureNumber = 41; num = 41; wid = 17.6; hts = [7]; cols = {1 3}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 3, 1, [12], 'Joshi and Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
movegui(fig_,[1,1]); % Move figure so it is visible on screen.

%%

colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};

axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

yd1L = nans(5,1); yd2L = nans(5,1);
yd1 = pa_All_s(:,1); yd1L(yd1>0) = log(yd1(yd1>0));
yd2 = pa_All_s(:,2); yd2L(yd2>0) = log(yd2(yd2>0));

xAx = 1:5;
bar(xAx-0.15,yd1,0.30,'facecolor','none','edgecolor',colors{2},'linewidth',1);
bar(xAx+0.15,yd2,0.30,'facecolor','none','edgecolor',colors{3},'linewidth',1);
% axis([0.5 5.5 0 0.06]);
xlabel('number of bins'); ylabel('probability of p<0.05');
text(1100,0.065,'Monkey Sp','fontweight','bold');
title('Shuffled trials ACC r_s_c analysis ','fontname','arial','fontsize',12);

yd1LA = nans(5,1); yd2LA = nans(5,2);
aa = load('propSig_ACC_051521')';
yd1LA = nans(5,1); yd2LA = nans(5,1);
yd1A = aa.pa_All(:,1); yd1LA(yd1A>0) = log(yd1A(yd1A>0)); yd1A = round(yd1A,2);
yd2A = aa.pa_All(:,2); yd2LA(yd2A>0) = log(yd2A(yd2A>0)); yd2A = round(yd2A,2);

for kind = 1:5
    
    %     text(xAx(kind)-0.25,yd1(kind)+0.005,num2str(yd1A(kind)),'color',colors{2});
    %     text(xAx(kind)+0.1,yd2(kind)+0.005,num2str(yd2A(kind)),'color',colors{3});
    
    plot(xAx(kind)-0.25,yd1A(kind),'o','markeredgecolor','none','markerfacecolor',colors{2},'markersize',8);
    plot(xAx(kind)+0.1,yd2A(kind),'o','markeredgecolor','none','markerfacecolor',colors{3},'markersize',8);
    
end

% bar(xAx-0.15,yd1A,0.25,'facecolor','none','edgecolor',colors{2},'linewidth',1);
% bar(xAx+0.15,yd2A,0.25,'facecolor','none','edgecolor',colors{3},'linewidth',1);
axis([0.5 5.5 0 0.4]);

% pause
%%

xbs = linspace(200,1000,nbin); % 032320.
xjit = 30;

for ti = 1:3
    
    axes(axs(ti+1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    for bi = 1:5
        
        mdat1 = mnk_rsc{1}(ti,bi);
        mdat2 = mnk_rsc{2}(ti,bi);
        
        mdat1_cia = mnk_rsc_iqr{1}(ti,bi,1); mdat1_cib = mnk_rsc_iqr{1}(ti,bi,2);
        mdat2_cia = mnk_rsc_iqr{1}(ti,bi,1); mdat1_cib = mnk_rsc_iqr{1}(ti,bi,2);
        
        mdat1_cia = mnk_rsc_iqr{2}(ti,bi,1); mdat1_cib = mnk_rsc_iqr{2}(ti,bi,2);
        mdat2_cia = mnk_rsc_iqr{2}(ti,bi,1); mdat1_cib = mnk_rsc_iqr{2}(ti,bi,2);
        
        bar(xbs(bi)-xjit,mdat1,2*xjit,'facecolor',colors{2},'edgecolor','none');
        bar(xbs(bi)+xjit,mdat2,2*xjit,'facecolor',colors{4},'edgecolor','none');
        
        plot([xbs(bi)-xjit xbs(bi)-xjit],[mdat1_cia mdat1_cib],'k-');
        plot([xbs(bi)+xjit xbs(bi)+xjit],[mdat2_cia mdat1_cib],'k-');
        
        
    end
    %     axis([100 1100 -0.0115 0.005]);
    axis([100 1100 -0.115 0.05]);
    
    if ti == 1
        ylabel('\DeltaACC rsc re:LCzero','fontname','arial','fontsize',10);
    end
    xlabel('binsize (msec)','fontname','arial','fontsize',10);
    
end


%% Pre 050421 Plots.

% for mm = 1:2
%     ddZ = mnkZ{mm};
%     ddNZ = mnkNZ{mm};
%     ddA = mnkA{mm};
%
%     all_pvals = nans(nbin,nshuft);
%     diff_all_rsc = nans(nbin,nshuft);
%     sigi = nans(nbin,nshuft);
%     median_diff_rsc = nans(1,nbin);
%     median_diff_rsc_iqr = nans(nbin,2);
%     for bi=1:nbin
%         for shufInd = 1:nshuft
%             bdZ = squeeze(ddZ(bi,:,shufInd));
%             bdNZ = squeeze(ddNZ(bi,:,shufInd));
%             bdA = squeeze(ddA(bi,:,shufInd));
%             gi = find(isfinite(bdZ) & isfinite(bdNZ) & abs(bdZ)>0);
%             diffVal = nans(1,length(gi));
%             diffVal = bdNZ(gi)-bdZ(gi);
%             all_pvals(bi,shufInd) = signrank(diffVal);
%             diff_all_rsc(bi,shufInd) = nanmedian(diffVal);
%         end
%         sigi(bi,:) = all_pvals(bi,:) < 0.05;
%         median_diff_rsc(bi) = nanmedian(diff_all_rsc(bi,:));
%         median_diff_rsc_iqr(bi,:) = prctile(diff_all_rsc(bi,:),[25 75])';
%     end
%     mnkSig{mm} = sigi;
%     mnk_rsc{mm} = median_diff_rsc;
%     mnk_rsc_iqr{mm} = median_diff_rsc_iqr;
% end

%% Set up plots.

% Don't plot results?
% plotYesNo = false;
% Plot results?
% plotYesNo = true;

% if plotYesNo
%     %% Plot without terciles.
%
%     % Setup figure:
%     figureNumber = 1; num = 1; wid = 17.6; hts = [7]; cols = {2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 3, 3, [12], 'Joshi and Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
%     movegui(fig_,[1,1]); % Move figure so it is visible on screen.
%
%     %%
%
%     xbs = linspace(200,1000,nbin); % 032320.
%
%     propSig1 = nansum(mnkSig{1},2)/nshuft;
%     xMax1 = max(propSig1); xMin1 = min(propSig1);
%     propSig2 = nansum(mnkSig{2},2)/nshuft;
%     xMax2 = max(propSig1); xMin2 = min(propSig1);
%
%     minY1 = min(mnk_rsc_iqr{1}(:,1)); maxY1 = max(mnk_rsc_iqr{1}(:,2));
%     minY2 = min(mnk_rsc_iqr{2}(:,1)); maxY2 = max(mnk_rsc_iqr{2}(:,2));
%     minY = min([minY1 minY2]);
%     maxY = max([maxY1 maxY2]);
%     minY = -0.015;
%     maxY = 0.015;
%
%
%     axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
%     bar(xbs,propSig1,'edgecolor','none','facecolor','k');
%     axis([100 1100 0 0.065]);
%     xlabel('bin size (ms)'); ylabel('proportion significant');
%     text(1100,0.065,'Monkey Sp','fontweight','bold');
%
%     axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
%     plot(xbs,mnk_rsc{1},'ko-');
%     plot(xbs,mnk_rsc_iqr{1}(:,1),'k--');
%     plot(xbs,mnk_rsc_iqr{1}(:,2),'k--');
%     plot([100 1100],[0 0],'k:');
%     axis([100 1100 minY maxY]);
%     xlabel('bin size (ms)'); ylabel({'% change in ACC rsc';'(LC_h_i_g_h re: LC_l_o_w)'});
%
%     axes(axs(3)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
%     bar(xbs,propSig2,'edgecolor','none','facecolor','k');
%     axis([100 1100 0 0.065]);
%     xlabel('bin size (ms)'); ylabel('proportion significant');
%     text(1100,0.065,'Monkey Ci','fontweight','bold');
%
%     axes(axs(4)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
%     plot(xbs,mnk_rsc{2},'ko-');
%     plot(xbs,mnk_rsc_iqr{2}(:,1),'k--');
%     plot(xbs,mnk_rsc_iqr{2}(:,2),'k--');
%     plot([100 1100],[0 0],'k:');
%     axis([100 1100 minY maxY]);
%     xlabel('bin size (ms)'); ylabel({'% change in ACC rsc';'(LC_h_i_g_h re: LC_l_o_w)'});
%
%
%
%     %%
%
%
%     %     legend({'zero','> 0'},'autoupdate','off','fontsize',8); legend('boxoff');
%     %     % plot([prcHi(1) prcHi(3)],[0.3 0.3],'-','linewidth',2,'color',c2);
%     %     plot([prcHi(1) prcHi(3)],[0.3 0.3],'k-','linewidth',2);
%     %     plot([prcHi(2)],[0.3],'kd','markersize',10,'linewidth',2);
%     %     xlabel('LC firing rate','fontsize',10);
%     %     ylabel('proportion of trials','fontsize',10);
%     %     % title('LC','fontsize',12);
%     %     % if jjm == 1 title('LC spike counts'); end
%     %     n_neurons = nOut(1)+nOut(2);
%     %     text(7,0.05,strcat('# neurons=',num2str(n_neurons)),'fontsize',8);
%     %     text(7,0.15,strcat('# trials =',num2str(NNN)),'fontsize',8);
%     %     axis([-1 xMax+1 -0.01 0.36]);
%     %     title('LC','fontsize',10);
%     %     set(gca,'fontname','arial');
%
%     %     pause
% end

