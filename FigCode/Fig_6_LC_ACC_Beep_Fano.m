% Fig_6_LC_ACC_Beep_Fano.m
%
% INFO: Plot beep driven spike rates and Fano factors.
% 
% Origin: LC_IC_ACC_Fano.m
% 120118: Combined LC_Fano, IC_Fano, ACC_Fano.
%
% FIXATION task ONLY.
%
% Calls "getLC_Fano" and "getACC_Fano" to do the heavy lifting.
%

%% Setup stuff:

% Clear everything?
clear; clear all;

%% Analysis Section: Calculate pairwise spike count correlations (rsc) conditioned on LC (or IC) spiking:

% Loops through monkeys and sessions.
% Uses function "getLC_DualArea_rsc.m" to do the heavy lifting.

% Don't reanalyze data?
reanalyzeLC = false;
% % Reanalyze data?
% reanalyzeLC = true;

if reanalyzeLC
    
    clear; clear all;
    
    monks = {'Sprout','Cicero1','Cicero2','Oz'}; % Add mnks as needed.
    % sites = {'LC_ACC_Fixation'}; % Only LC for now.
    
    % Base directory for brain area.
    base_dir = {
        'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Sprout\LC_ACC_Fixation\clean';
        'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Cicero\LC_ACC_Fixation\clean';
        'C:\Sidd\PostDoc2\Data\ALL_PUPIL_DATA\LC\Cicero\clean'
        'C:\Sidd\PostDoc2\Data\ALL_PUPIL_DATA\LC\Oz\clean'
        };
    
    nMonks = length(monks);
    % nSites = length(sites);
    
    mnkRes = [];
    binSiz = 1; % In case want to look at coincidences per bin for bins> 1msec.
    minTrials = 5;
    
    for mm = 1:nMonks % Loop over monkeys (wheeeee....!).
        
        inDir= base_dir{mm};
        cd(inDir);  dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
        fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
        pair_dat = [];
        
        if ~isempty(fnames)
            nf = length(fnames);
            % for ff = 1:10 % For testing.
            for ff = 1:nf % Loop over sessions for this mnky.
                % ff=4
                load(fnames{ff});
                % pause
                pair_dat{ff} = getLC_Fano(siteData,mm,minTrials); % Added analyses for window size effects.
                disp(sprintf('File %d%sof%s%d', ff,' ',' ',nf));
            end % Sessions for this mnky.
        end
        mnkRes{mm} = pair_dat;
    end % Mnks.
    
    % Save analysis:
    savedir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC\';
    % savedir = strcat([base_dir,'\Results\Results_June_2017']);
    cd(savedir);
    
    % save LC_ACC_Results_rsc_082618b mnkRes;
    % save LC_Fano_Results_121718 mnkRes;
    % save LC_Fano_Results_010119 mnkRes;
    % save LC_Fano_Results_120819 mnkRes;
    % save LC_Fano_Results_121119 mnkRes;
    %     save LC_Fano_Results_042620 mnkRes;
    % save LC_Fano_Results_092420 mnkRes;
    save LC_Fano_Results_051821 mnkRes;
    
end % If reanalyze loop.

%% ACC.

% Loops through monkeys and sessions.
% Uses function "getACC_DualArea_rsc.m" to do the heavy lifting.

% Don't reanalyze data?
reanalyzeACC = false;
% Reanalyze data?
% reanalyzeACC = true;

if reanalyzeACC
    
    clear; clear all;
    
    monks = {'Sprout','Cicero2'}; % Add mnks as needed.
    % sites = {'LC_ACC_Fixation'}; % Only LC for now.
    
    % Base directory for brain area.
    base_dir = {
        'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Sprout\LC_ACC_Fixation\clean';
        'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Cicero\LC_ACC_Fixation\clean';
        };
    
    nMonks = length(monks);
    % nSites = length(sites);
    
    mnkRes = [];
    binSiz = 1; % In case want to look at coincidences per bin for bins> 1msec.
    minTrials = 5;
    
    for mm = 1:nMonks % Loop over monkeys (wheeeee....!).
        
        inDir= base_dir{mm};
        cd(inDir);  dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
        fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
        pair_dat = [];
        
        if ~isempty(fnames)
            nf = length(fnames);
            % for ff = 1:10 % For testing.
            for ff = 1:nf % Loop over sessions for this mnky.
                % ff=20
                load(fnames{ff});
                % pause
                pair_dat{ff} = getACC_Fano(siteData,mm,minTrials); % Added analyses for window size effects.
                disp(sprintf('File %d%sof%s%d', ff,' ',' ',nf));
            end % Sessions for this mnky.
        end
        mnkRes{mm} = pair_dat;
    end % Mnks.
    
    % Save analysis:
    savedir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC\';
    % savedir = strcat([base_dir,'\Results\Results_June_2017']);
    cd(savedir);
    
    % save LC_ACC_Results_rsc_082618b mnkRes;
    % save ACC_Fano_Results_121718 mnkRes;
    %     save ACC_Fano_Results_010119 mnkRes;
    % save ACC_Fano_Results_120619 mnkRes;
    % save ACC_Fano_Results_121119 mnkRes;
    % save ACC_Fano_Results_042620 mnkRes;
    save ACC_Fano_Results_051821 mnkRes;
    
end % If reanalyze loop.

% ****************************************************

%% Plotting Section - Load results file created above, calculate summary statistics and plot stuff:

% Don't plot results?
plotYesNo = false;
% Plot results?
plotYesNo = true;

%%


if plotYesNo
    % Clear everything?
    clear; clear all;
    
    %     % Load colormap.
    %     cd C:\Sidd\PostDoc2\SidCode\General;
    %     load('sidd_colors');
    %     c1 = cmap_M(110,:);
    %     c2 = cmap_M(170,:);
    
    %% LC.
    
    base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC\';
    
    % inFile = 'LC_Fano_Results_092420';
    inFile = 'LC_Fano_Results_051821';
    
    savedir = base_dir;
    cd(savedir); load(inFile);
    monks = {'Sprout', 'Cicero1','Cicero2','Oz'}; % Add mnks as needed...
    nMonks = length(monks);
    
    % **********************************************
    % Initialize data matrices for plot data.
    mnkData = [];
    
    for mm = 1:nMonks % Loop over mnks....
        
        mnkDat = mnkRes{mm}; % Get mnky data from saved analysis.
        
        FF_unit_dat = [];
        unitMean1 = []; unitVar1 = [];
        unitMean2 = []; unitVar2 = [];
        unitMean3 = []; unitVar3 = [];
        unitMean4 = []; unitVar4 = [];
        sess_id = [];
        bsi = 0; binned_sp = [];
        binned_sp1 = [];
        
        for sesi = 1:length(mnkDat)
            sesDat = mnkDat{sesi};
            if ~isempty(sesDat)
                
                %                 {nL_count nL_count1 nL_count2 nL_count3 nA FF_unit nA1};
                
                % pause
                FF_unit_dat = [FF_unit_dat; sesDat{6}];
                
                for ui = 1:size(sesDat{1},2)
                    bsi = bsi+1;
                    
                    spCount1 = sesDat{1}{ui}; % No beep.
                    spCount2 = sesDat{2}{ui}; % Before beep.
                    spCount3 = sesDat{3}{ui}; % Phasic.
                    spCount4 = sesDat{4}{ui}; % After phasic.
                    
                    % if (nanmean(spCount1)>0 && nanmean(spCount2)>0 && nanmean(spCount3)>0 && nanmean(spCount4) > 0)
                    if (nanmean(spCount2)>0 && nanmean(spCount3)>0 && nanmean(spCount4) > 0)
                        if (nanvar(spCount2)>0 && nanvar(spCount3)>0 && nanvar(spCount4) > 0)
                            unitMean1 = [unitMean1 nanmean(spCount1)]; unitVar1 = [unitVar1 nanvar(spCount1)];
                            unitMean2 = [unitMean2 nanmean(spCount2)]; unitVar2 = [unitVar2 nanvar(spCount2)];
                            unitMean3 = [unitMean3 nanmean(spCount3)]; unitVar3 = [unitVar3 nanvar(spCount3)];
                            unitMean4 = [unitMean4 nanmean(spCount4)]; unitVar4 = [unitVar4 nanvar(spCount4)];
                            sess_id = [sess_id sesi];
                        end
                    end
                    binned_sp{bsi} = sesDat{5}{ui};
                    binned_sp1{bsi} = sesDat{7}{ui};
                end
                mnkData{mm} = {unitMean1 unitMean2 unitMean3 unitMean4 unitVar1 unitVar2 unitVar3 unitVar4 binned_sp sess_id FF_unit_dat binned_sp1};
            end
        end
    end
    
    %% ACC.
    
    base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC\';
    
    % inFile = 'ACC_Fano_Results_042620';
    inFile = 'ACC_Fano_Results_051821';
    
    savedir = base_dir;
    cd(savedir); load(inFile);
    monks = {'Sprout','Cicero2'}; % Add mnks as needed...
    nMonks = length(monks);
    
    % **********************************************
    % Initialize data matrices for plot data.
    mnkData3 = [];
    
    for mm = 1:nMonks % Loop over mnks....
        
        mnkDat = mnkRes{mm}; % Get mnky data from saved analysis.
        
        FF_unit_dat = [];
        unitMean1 = []; unitVar1 = [];
        unitMean2 = []; unitVar2 = [];
        unitMean3 = []; unitVar3 = [];
        unitMean4 = []; unitVar4 = [];
        sess_id = [];
        bsi = 0; binned_sp = [];
        binned_sp1 = [];
        
        for sesi = 1:length(mnkDat)
            sesDat = mnkDat{sesi};
            if ~isempty(sesDat)
                
                % pause
                FF_unit_dat = [FF_unit_dat; sesDat{6}];
                
                for ui = 1:size(sesDat{1},2)
                    
                    bsi = bsi+1;
                    
                    spCount1 = sesDat{1}{ui};
                    spCount2 = sesDat{2}{ui};
                    spCount3 = sesDat{3}{ui};
                    spCount4 = sesDat{4}{ui};
                    
                    
                    if (nanmean(spCount1) > 0 && nanmean(spCount2) > 0 && nanmean(spCount3) > 0 && nanmean(spCount4) > 0)
                        if (nanvar(spCount2)>0 && nanvar(spCount3)>0 && nanvar(spCount4) > 0)
                            unitMean1 = [unitMean1 nanmean(spCount1)]; unitVar1 = [unitVar1 nanvar(spCount1)];
                            unitMean2 = [unitMean2 nanmean(spCount2)]; unitVar2 = [unitVar2 nanvar(spCount2)];
                            unitMean3 = [unitMean3 nanmean(spCount3)]; unitVar3 = [unitVar3 nanvar(spCount3)];
                            unitMean4 = [unitMean4 nanmean(spCount4)]; unitVar4 = [unitVar4 nanvar(spCount4)];
                            sess_id = [sess_id sesi];
                        end
                    end
                    binned_sp{bsi} = sesDat{5}{ui};
                    binned_sp1{bsi} = sesDat{7}{ui};
                end
                mnkData3{mm} = {unitMean1 unitMean2 unitMean3 unitMean4 unitVar1 unitVar2 unitVar3 unitVar4 binned_sp sess_id FF_unit_dat binned_sp1};
            end
        end
    end
    
    % pause
    
    %% Setup figure:
    
    figureNumber = 7; num = 7; wid = 15; hts = [5]; cols = {2 2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 2, 2, [12], 'Joshi and Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    % bs = linspace(100,1000,10); % 032220b.
    bsiz = 100; % April 2021.
    backTime = -1000; fwdTime = 1000; % 032320b.
    tmin   = backTime; tmax   = fwdTime; xBin = []; nb = []; nbs = length(bsiz);
    tsize  = bsiz; tstep  = tsize/2;
    xBin  = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
    nb = size(xBin,1);
    
    bsiz2 = 200; % April 2021.
    backTime2 = -1000; fwdTime2 = 1000; % 032320b.
    tmin2   = backTime2; tmax2   = fwdTime2; xBin2 = []; nb2 = []; nbs2 = length(bsiz2);
    tsize2  = bsiz2; tstep2  = tsize2/2;
    xBin2 = [(tmin2:tstep2:tmax2-tsize2)' (tmin2+tsize2:tstep2:tmax2)'];
    nb2 = size(xBin2,1);
    
    xAx1 = -1000:1000;
    winSiz = length(backTime:fwdTime)/1000;
    xAx = nanmean(xBin,2);
    winSiz2 = length(backTime2:fwdTime2)/1000;
    xAx2 = nanmean(xBin2,2);
    
    %     pause
    %% Plot LC example PSTH.
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    % binned_spikes1 = mnkData{1}{9}{23}; % Binned at 100 ms.
    binned_spikes1 = mnkData{1}{12}{23}; % Binned at 1 ms.
    
    %     phasic_sp = mnkData{1}{3}; phasic_resp = 5*nanmean(phasic_sp);
    
    ntrials = size(binned_spikes1,1);
    % binned_mean = (1000/bsiz)*nanmean(binned_spikes1);
    binned_mean = nanmean(binned_spikes1);
    maxFR = (1000/25)*max(binned_mean);
    minFR = (1000/25)*min(binned_mean);
    binned_mean = smooth(binned_mean,25);
    psth = ntrials*binned_mean/max(binned_mean);    % binned_se = smooth(nanse(1000*binned_spikes1{yy},1),100);
    % plot(xAx,psth,'-','linewidth',1,'color','k'); % errorBarFilled(xAx,binned_mean',binned_se',[0.5 0.5 0.5]);
    plot(xAx1,psth,'-','linewidth',1,'color','k'); % errorBarFilled(xAx,binned_mean',binned_se',[0.5 0.5 0.5]);
    for kk = 1:ntrials
        trial_spikes = binned_spikes1(kk,:);
        a1 = find(trial_spikes>0);
        for jj = 1:length(a1)
            % plot([xAx(a1(jj)) xAx(a1(jj))],[(kk-0.4) (kk+0.4)],'k-','color',[.5 .5 .5],'linewidth',0.5);
            plot([xAx1(a1(jj)) xAx1(a1(jj))],[(kk-0.4) (kk+0.4)],'k-','color',[.5 .5 .5],'linewidth',0.5);
        end
    end
    text(200,32,['max fr = ',num2str(maxFR),'sp/sec'],'fontsize',8,'fontname','arial');
    text(200,28,['min fr = ',num2str(minFR),'sp/sec'],'fontsize',8,'fontname','arial');
    ylabel('firing rate (sp/s)','fontsize',10,'fontname','arial');
    xlabel('time re:beep (ms)','fontsize',10,'fontname','arial');
    plot([0 0],[0 kk],'k--','linewidth',1); plot([xAx(1) xAx(end)],[0 0],'k--','linewidth',1);
    axis([-1000 1000 -1 kk+1]); title('LC', 'fontsize',10,'fontname','arial'); set(gca,'fontname','arial');
    
    %     pause;
    
    %% Population PSTH.
    
    % Averaged psth's across all monkeys.
    msr1 = []; ntLim = 3;
    % binned_spikes1 = mnkData{1}{9};
    binned_spikes1 = mnkData{1}{12}; bsiz = 1;
    for qq = 1:length(binned_spikes1)
        if size (binned_spikes1{qq},1) > ntLim
            msr1 = [msr1; (1000/bsiz)*nanmean(binned_spikes1{qq})];
        end
    end
    msr2 = [];
    % binned_spikes2 = mnkData{2}{9};
    binned_spikes2 = mnkData{2}{12}; bsiz = 1;
    for qq = 1:length(binned_spikes2)
        if size (binned_spikes2{qq},1) > ntLim
            msr2 = [msr2; (1000/bsiz)*nanmean(binned_spikes2{qq})];
        end
    end
    msr3 = [];
    % binned_spikes3 = mnkData{3}{9};
    binned_spikes3 = mnkData{3}{12}; bsiz = 1;
    for qq = 1:length(binned_spikes3)
        if size (binned_spikes3{qq},1) > ntLim
            msr3 = [msr3; (1000/bsiz)*nanmean(binned_spikes3{qq})];
        end
    end
    msr4 = [];
    % binned_spikes4 = mnkData{4}{9};
    binned_spikes4 = mnkData{4}{12}; bsiz = 1;
    for qq = 1:length(binned_spikes4)
        if size (binned_spikes4{qq},1) > ntLim
            msr4 = [msr4; (1000/bsiz)*nanmean(binned_spikes4{qq})];
        end
    end
    
    phasic_sp = (1000/200)*[mnkData{1}{3} mnkData{2}{3} mnkData{3}{3} mnkData{4}{3}];
    phasic_raw = nanmean(phasic_sp);
    rawRate = phasic_raw;
    phasic_raw_se = nanse(phasic_sp,2);
    rawSe = phasic_raw_se;
    
    msr = [zscore(msr1')';zscore(msr2')';zscore(msr3')';zscore(msr4')'];
    meanRate = nanmean(msr);
    seRate = nanse(msr,1);
    
    meanRate = smooth(nanmean(msr),25);
    seRate = smooth(nanse(msr,1),25);
    
    axes(axs(3)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    plot(xAx1,meanRate,'-','color','k');
    errorBarFilled(xAx1,meanRate',seRate','k');
    
    n = size(msr1,1)+size(msr2,1)+size(msr3,1)+size(msr4,1);
    text(-900,1.2,['n=',num2str(n)],'fontsize',8);
    text(-900,1.0,['avg phasic=',num2str(rawRate),'+-',num2str(rawSe),'sp/sec'],'fontsize',8,'fontname','arial');
    xlabel('time re:beep (ms)','fontsize',10,'fontname','arial');
    ylabel([{'zscored';'firing rate'}],'fontsize',10,'fontname','arial');
    plot([0 0],[-0.16 1.1],'k--','linewidth',1); plot([xAx(1) xAx(end)],[0 0],'k--','linewidth',1);
    axis([-1000 1000 -0.16 1.1]); set(gca,'fontname','arial');
    
    %% Fano factor.
    
    axes(axs(5)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    FF1 = mnkData{1}{11}; FF2 = mnkData{2}{11}; FF3 = mnkData{3}{11}; FF4 = mnkData{4}{11};
    
    %     FFA = [FF1;FF2]; % 051821.
    FFA = [FF1;FF2;FF3;FF4];
    
    mean_FFA = nanmean(FFA);
    se_FFA = nanse(FFA,1);
    
    errorBarFilled(xAx2',mean_FFA,se_FFA,[0.5 0.5 0.5]);
    plot(xAx2,mean_FFA,'k-','linewidth',1);
    plot([0 0],[-0.92 1.3],'k--'); axis([-1000 1000 0.92 1.25]);
    xlabel('time re:beep (ms)','fontsize',10,'fontname','arial');
    ylabel('Fano factor','fontsize',10,'fontname','arial');
    set(gca,'fontname','arial');
    
    %% Next:  ACC.
    
    % Plot ACC example PSTH.
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    %     % *******************************************************************************************
    % %     for kii = 1:length(mnkData3{1}{9})
    %
    % % kii=324
    %     % figure;
    % %         hold on;
    %         % *******************************************************************************************
    
    kii = 183;
    
    % binned_spikes1 = mnkData{1}{9}{kii}; % Binned at 100 ms.
    binned_spikes1 = mnkData3{1}{12}{kii}; % Binned at 1 ms.
    
    %     phasic_sp = mnkData{1}{3}; phasic_resp = 5*nanmean(phasic_sp);
    
    ntrials = size(binned_spikes1,1);
    % binned_mean = (1000/bsiz)*nanmean(binned_spikes1);
    binned_mean = nanmean(binned_spikes1);
    maxFR = (1000/25)*max(binned_mean);
    minFR = (1000/25)*min(binned_mean);
    binned_mean = smooth(binned_mean,25);
    psth = ntrials*binned_mean/max(binned_mean);    % binned_se = smooth(nanse(1000*binned_spikes1{yy},1),100);
    % plot(xAx,psth,'-','linewidth',1,'color','k'); % errorBarFilled(xAx,binned_mean',binned_se',[0.5 0.5 0.5]);
    plot(xAx1,psth,'-','linewidth',1,'color','k'); % errorBarFilled(xAx,binned_mean',binned_se',[0.5 0.5 0.5]);
    for kk = 1:ntrials
        trial_spikes = binned_spikes1(kk,:);
        a1 = find(trial_spikes>0);
        for jj = 1:length(a1)
            % plot([xAx(a1(jj)) xAx(a1(jj))],[(kk-0.4) (kk+0.4)],'k-','color',[.5 .5 .5],'linewidth',0.5);
            plot([xAx1(a1(jj)) xAx1(a1(jj))],[(kk-0.4) (kk+0.4)],'k-','color',[.5 .5 .5],'linewidth',0.5);
        end
    end
    text(200,55,['max fr = ',num2str(maxFR),'sp/sec'],'fontsize',8,'fontname','arial');
    text(200,50,['min fr = ',num2str(minFR),'sp/sec'],'fontsize',8,'fontname','arial');
    ylabel('firing rate (sp/s)','fontsize',10,'fontname','arial');
    xlabel('time re:beep (ms)','fontsize',10,'fontname','arial');
    plot([0 0],[0 kk],'k--','linewidth',1); plot([xAx(1) xAx(end)],[0 0],'k--','linewidth',1);
    axis([-1000 1000 -1 kk+5]); title('ACC', 'fontsize',10,'fontname','arial'); set(gca,'fontname','arial');
    
    %% Population PSTH.
    
    % Averaged psth's across all monkeys.
    msr1 = [];
    % binned_spikes1 = mnkData3{1}{9};
    binned_spikes1 = mnkData3{1}{12}; bsiz = 1;
    for qq = 1:length(binned_spikes1)
        if size (binned_spikes1{qq},1) > 5
            msr1 = [msr1; (1000/bsiz)*nanmean(binned_spikes1{qq})];
        end
    end
    msr2 = [];
    % binned_spikes2 = mnkData3{2}{9};
    binned_spikes2 = mnkData3{2}{12}; bsiz = 1;
    for qq = 1:length(binned_spikes2)
        if size (binned_spikes2{qq},1) > 5
            msr2 = [msr2; (1000/bsiz)*nanmean(binned_spikes2{qq})];
        end
    end
    
    phasic_sp = (1000/200)*[mnkData{1}{3} mnkData{2}{3} mnkData{3}{3} mnkData{4}{3}];
    phasic_raw = nanmean(phasic_sp);
    rawRate = phasic_raw;
    phasic_raw_se = nanse(phasic_sp,2);
    rawSe = phasic_raw_se;
    
    msr = [zscore(msr1')';zscore(msr2')'];
    meanRate = nanmean(msr);
    seRate = nanse(msr,1);
    
    meanRate = smooth(nanmean(msr),25);
    seRate = smooth(nanse(msr,1),25);
    
    axes(axs(4)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    plot(xAx1,meanRate,'-','color','k');
    errorBarFilled(xAx1,meanRate',seRate','k');
    
    n = size(msr1,1)+size(msr2,1);
    text(-900,0.11,['n=',num2str(n)],'fontsize',8,'fontname','arial');
    text(-900,0.08,['avg phasic=',num2str(rawRate),'+-',num2str(rawSe),'sp/sec'],'fontsize',8);
    xlabel('time re:beep (ms)','fontsize',10,'fontname','arial');
    ylabel([{'zscored';'firing rate'}],'fontsize',10,'fontname','arial');
    plot([0 0],[-0.1 0.14],'k--'); plot([xAx(1) xAx(end)],[0 0],'k--');
    axis([-1000 1000 -0.1 0.14]); set(gca,'fontname','arial');
    
    %% Fano factor.
    
    axes(axs(6)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    FF1 = mnkData3{1}{11}; FF2 = mnkData3{2}{11};
    FFA = [FF1;FF2];
    
    %         gi = find(isfinite(sum(FFA,2)));
    %         FFA_f = FFA(gi,:);
    
    FFA_f = FFA;
    
    %     for ki = 1:size(FFA_f,2)
    %         FFA_f(ki,:) = smooth(FFA_f(ki,:),20);
    %     end
    
    mean_FFA = nanmean(FFA_f);
    se_FFA = nanse(FFA_f,1);
    
    errorBarFilled(xAx2',mean_FFA,se_FFA,[0.5 0.5 0.5]);
    plot(xAx2,mean_FFA,'k-','linewidth',2);
    plot([0 0],[1.37 1.66],'k--'); axis([-1000 1000 1.37 1.66]);
    xlabel('time re:beep (ms)','fontsize',10,'fontname','arial');
    ylabel('Fano factor','fontsize',10,'fontname','arial');
    set(gca,'fontname','arial');
    
    
    %% SUPPLEMENTARY FIGURE.
    
    figureNumber = 71; num = 71; wid = 15; hts = [10]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 2, 2, [12], 'Joshi and Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    %%
    
    % Averaged psth's across all monkeys.
    msr1A = [];
    binned_spikes1 = mnkData3{1}{12};
    for qq = 1:length(binned_spikes1)
        if size (binned_spikes1{qq},1) > 5
            msr1A = [msr1A; nanmean(binned_spikes1{qq})];
        end
    end
    msr2A = [];
    binned_spikes2 = mnkData3{2}{12};
    for qq = 1:length(binned_spikes2)
        if size (binned_spikes2{qq},1) > 5
            msr2A = [msr2A; nanmean(binned_spikes2{qq})];
        end
    end
    
    % LC neuron example psth showing eopchs.
    msr1AS=nans(size(msr1A));
    for kki = 1:size(msr1A,1)
        msr1AS(kki,:) = smooth(msr1A(kki,:),50);
    end
    
    msr2AS=nans(size(msr2A));
    for kki = 1:size(msr2A,1)
        msr2AS(kki,:) = smooth(msr2A(kki,:),50);
    end
    
    msr1 = msr1AS;
    msr2 = msr2AS;
    
    msr1Z = zscore(msr1')';
    di = find(sum(msr1Z,2)==0 | ~isfinite(sum(msr1Z,2)));
    msr1Z(di,:) = [];
    msr1_phasic=nanmean(msr1Z(:,1025:1100),2);
    
    msr2Z = zscore(msr2')';
    di = find(sum(msr2Z,2)==0 | ~isfinite(sum(msr2Z,2)));
    msr2Z(di,:) = [];
    msr2_phasic=nanmean(msr2Z(:,1010:1110),2);
    
    [pk pindex1] = sort(msr1_phasic);
    [pk pindex2] = sort(msr2_phasic);
    
    msr1_sort = msr1Z(pindex1,751:1251);
    msr2_sort = msr2Z(pindex2,751:1251);
    
    xx = -250:250;
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    yy = 1:size(msr1,1);
    colormap(gray);
    imagesc(xx,yy,msr1_sort); colorbar;
    plot([0 0],[1 max(yy)],'r--','linewidth',0.5);
    axis([min(xx) max(xx) min(yy) max(yy)]);
    xlabel('time re:beep','fontname','arial','fontsize',10);
    ylabel('Unit #','fontname','arial','fontsize',10);
    title('Monkey Sp','fontname','arial','fontsize',12);
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    yy = 1:size(msr2,1);
    colormap(gray);
    imagesc(xx,yy,msr2_sort); colorbar;
    plot([0 0],[1 max(yy)],'r--','linewidth',0.5);
    axis([min(xx) max(xx) min(yy) max(yy)]);
    xlabel('time re:beep','fontname','arial','fontsize',10);
    ylabel('Unit #','fontname','arial','fontsize',10);
    title('Monkey Ci','fontname','arial','fontsize',12);
    
    
end % Plot Yes/No.
