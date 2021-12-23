% Fig_9_LC_ACC_Covariation.m
%
% INFO: Compare LC spiking and ACC rsc time courses.
% NO BEEP trials only.
% Previous version: LC (ACC) STA of ACC (LC) spiking, Fano and pairwise spike count correlation (rsc).
% Dual area recording - figure script; FIXATION task.

% For paper ... Combined LC_Dual_Recording_sts.m and ACC_Dual_Recording_sts.m
%
% OLD: Plots ACC (LC) spike rate aligned to LC (ACC) spiking.
% OLD: Plots ACC (LC) pairwise spike count correlations (rsc); ACC (LC) pairs with ACC (LC) unit spikes aligned to LC (ACC) spikes.
%
% Origin of "LC_Dual_Recording_sts.m": 111016 - Sidd.
% Modified: 062116, 070216, 072518, 092618.
% Mod: 112918: Cleaning and finalizing - Sidd.
% 010519: "finalizing" .... ha ha.
% Calls "getLC_DualArea_Covariation" to do the heavy lifting.

%% Setup stuff:

% Clear everything?
clear; clear all;

%% Setup names and directories:

monks = {'Sprout','Cicero'}; % Add mnks as needed.
sites = {'LC_ACC_Fixation'}; % Add sites as needed.
base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.
nMonks = length(monks); nSites = length(sites);

%% LC STA of ACC spiking.
% Loop over clean datafiles and do the math.

% Don't reanalyze data?
reanalyze = false;
% % Reanalyze data?
% reanalyze = true;

mnkRes = []; % Analyzed data structure.
if reanalyze
    for mm = 1:nMonks % Loop over monkeys.
        sitRes = [];
        for ss = 1:nSites % Loop over brainstem sites.
            % Get directory and files.
            inDir= strcat([base_dir,'\',monks{mm},'\',sites{ss},'\clean']); % Create dir name for input (clean) files.
            cd(inDir);  dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
            fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
            CP_dat = []; % This monkey,this site data structure.
            if ~isempty(fnames)
                nf = length(fnames); % Number of data files to analyze.
                for ff = 1:nf
                    load(fnames{ff}); % Load "clean" data file.
                    % pause
                    CP_dat{ff} = getLC_DualArea_Covariation(siteData); % Call "getLC_DualArea_Covariation" to do heavy lifting.
                    disp(sprintf('File %d%sof%s%d', ff,' ',' ',nf)); % Display progress on command window.
                end
            end
            sitRes{ss} = CP_dat; % Collect this site analyses.
        end % Sites loop.
        mnkRes{mm} = sitRes; % Collect this monkey analyses.
    end % Mnks loop.
    
    % **************************************************************
    % savedir = strcat([base_dir,'\Results\Results_2018\LC_ACC\']);
    savedir = strcat([base_dir,'\Results\Results_2021\LC_ACC\']);
    cd(savedir);
    
    % save LC_ACC_CCG_092021b mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    % save LC_ACC_CCG_092121b mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    % save LC_ACC_CCG_092321a mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    % save LC_ACC_CCG_100721 mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    %     save LC_ACC_CCG_101021b mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    save LC_ACC_CCG_101121 mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    
    % **************************************************************
end

%% Plot results for LC spike triggered ACC spiking.

% % Don't plot results.
plotYesNo = false;
% % Plot results.
plotYesNo = true;

if plotYesNo
    
    %%
    clear; clear all;
    
    %%
    
    % inFile = 'LC_ACC_CCG_090621_500b';
    % inFile = 'LC_ACC_CCG_090721_500i';
    % inFile = 'LC_ACC_CCG_091221_500';
    % inFile = 'LC_ACC_CCG_091921d';
    % inFile = 'LC_ACC_CCG_092021b';
    % inFile = 'LC_ACC_CCG_092121b';
    % inFile = 'LC_ACC_CCG_092321a';
    % inFile = 'LC_ACC_CCG_100721';
    % inFile = 'LC_ACC_CCG_101021b';
    inFile = 'LC_ACC_CCG_101121';
    
    smW = 1;
    base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.
    savedir = strcat([base_dir,'\Results\Results_2021\LC_ACC\']);
    cd(savedir);
    load(inFile);
    
    % **********************************************
    
    monks = {'Sprout','Cicero'}; % Add mnks as needed...
    sites = {'LC_ACC_Fixation'}; % Add sites as needed...
    nMonks = length(monks); nSites = length(sites);
    
    dat_trial_FRAB = cell(1,2); dat_trial_FRAA = cell(1,2);
    dat_trial_rscBz = cell(1,2); dat_trial_rscAz = cell(1,2);
    dat_trial_rscBnz = cell(1,2); dat_trial_rscAnz = cell(1,2);
    dat_trial_ccgB = cell(1,2); dat_trial_ccgA = cell(1,2);
    dat_trial_ccgB2 = cell(1,2); dat_trial_ccgA2 = cell(1,2);
    dat_trial_FRAB2 = cell(1,2); dat_trial_FRAA2 = cell(1,2);
    
    % Data structure.
    % ccg_dat = {out_ccg_all taxB taxF tax2 TT};
    % out_ccg_all{li} = {out_ccgB out_ccgA out_FR_AB out_FR_AA out_rsc_ABnz out_rsc_AAnz out_rsc_ABz out_rsc_AAz};
    
    % {out_ccgB out_ccgA out_FRL_AB out_FRL_AA out_rsc_ABnz out_rsc_AAnz out_rsc_ABz out_rsc_AAz out_FRA_AB out_FRA_AA out_ccgAB out_ccgAA};
    
    % First - get the processed data:
    for mm = 1:nMonks % Loop over mnks.
        mnkDat = mnkRes{mm};
        for ss = 1:nSites % Loop over brainstem sites.
            sitDat = mnkDat{ss};
            for ap = 1:length(sitDat) % Iterate over sessions.
                for li = 1:length(sitDat{ap}{1}) % Iterate over LC units.
                    
                    dat_trial_ccgB{mm} = [dat_trial_ccgB{mm}; sitDat{ap}{1}{li}{1}];
                    dat_trial_ccgA{mm} = [dat_trial_ccgA{mm}; sitDat{ap}{1}{li}{2}];
                    dat_trial_FRAB{mm} = [dat_trial_FRAB{mm}; sitDat{ap}{1}{li}{3}];
                    dat_trial_FRAA{mm} = [dat_trial_FRAA{mm}; sitDat{ap}{1}{li}{4}];
                    dat_trial_rscBnz{mm} = [dat_trial_rscBnz{mm}; sitDat{ap}{1}{li}{5}];
                    dat_trial_rscAnz{mm} = [dat_trial_rscAnz{mm}; sitDat{ap}{1}{li}{6}];
                    dat_trial_rscBz{mm} = [dat_trial_rscBz{mm}; sitDat{ap}{1}{li}{7}];
                    dat_trial_rscAz{mm} = [dat_trial_rscAz{mm}; sitDat{ap}{1}{li}{8}];
                    
                    dat_trial_ccgB2{mm} = [dat_trial_ccgB2{mm}; sitDat{ap}{1}{li}{11}];
                    dat_trial_ccgA2{mm} = [dat_trial_ccgA2{mm}; sitDat{ap}{1}{li}{12}];
                    dat_trial_FRAB2{mm} = [dat_trial_FRAB2{mm}; sitDat{ap}{1}{li}{9}];
                    dat_trial_FRAA2{mm} = [dat_trial_FRAA2{mm}; sitDat{ap}{1}{li}{10}];
                    
                end
            end % Session loop.
        end % Sites loop.
    end % Mnks loop.
    
    %%
    taxB =  sitDat{ap}{2};
    taxF =  sitDat{ap}{3};
    tax2 =  sitDat{ap}{4}; % tax2 = tax2(2:end-1);
    % TT =  sitDat{ap}{4};
    
    %     out_ccg = nans(1000,length(tax2));
    %     dat_rsc_all = [dat_trial_rsc{1};dat_trial_rsc{2}]; nrsc = size(dat_rsc_all,1);
    %     dat_FR_all = [dat_trial_FR{1};dat_trial_FR{2}]; nFR = size(dat_FR_all,1);
    %     for tt = 1:1000
    %         raw_ccg = xcorr(zscore(dat_rsc_all(1+floor(rand(1)*nrsc),:)),zscore(dat_FR_all(1+floor(rand(1)*nFR),:)));
    %         out_ccg(tt,:) = raw_ccg./TT;
    %     end
    
    
    %% Setup figure:
    figureNumber = 10; num = 10; wid = 17.6; hts = [8]; cols = {3}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1, 1, [12], '', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    %%
    
    %     pause
    
    %% LC FR
    
    peth1 = dat_trial_FRAB{1}; peth2 = dat_trial_FRAB{2}; pethAllB = [peth1;peth2];
    peth1 = dat_trial_FRAA{1}; peth2 = dat_trial_FRAA{2}; pethAllA = [peth1;peth2];
    
    nb = size(pethAllB,2);
    erbarsB = nans(2,nb);
    erbarsA = nans(2,nb);
    
    for uu = 1:nb
        
        datTemp1 = pethAllB(:,uu);
        datTemp1 = datTemp1(~isnan(datTemp1));
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
        erbarsB(:,uu) = ci_hi_temp1;
        
        datTemp1 = pethAllA(:,uu);
        datTemp1 = datTemp1(~isnan(datTemp1));
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
        erbarsA(:,uu) = ci_hi_temp1;
        
    end
    
    %% Plot.
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    plot(taxB,nanmean(pethAllB),'k-','linewidth',2);
    plot(taxF,nanmean(pethAllA),'k-','linewidth',2);
    plot(taxB,erbarsB(1,:),'k:','linewidth',1); plot(taxB,erbarsB(2,:),'k:','linewidth',1);
    plot(taxF,erbarsA(1,:),'k:','linewidth',1); plot(taxF,erbarsA(2,:),'k:','linewidth',1);
    % legend('LC FR','ACC \it r_s_c','autoupdate','off','fontsize',8); legend('boxoff');
    xlabel('time re:start of stable fixation (ms)');
    ylabel('LC firing rate (sp/s)');
    % plot([tax(1) tax(end)],[0 0],'k--');
    minY = 1.75; maxY = 2.85;
    % minY = 1.7; maxY = 3.4;
    %     minY = 0; maxY = 8;
    axis([taxB(1)-50 taxF(end)+50 minY maxY]);
    plot([0 0],[minY maxY],'k--','linewidth',0.5);
    
    %% ACC FR
    
    %     peth1 = dat_trial_FRAB2{1}; peth2 = dat_trial_FRAB2{2}; pethAllB2 = [peth1;peth2];
    %     peth1 = dat_trial_FRAA2{1}; peth2 = dat_trial_FRAA2{2}; pethAllA2 = [peth1;peth2];
    %
    %     nb = size(pethAllB2,2);
    %     erbarsB2 = nans(2,nb);
    %     erbarsA2 = nans(2,nb);
    %
    %     for uu = 1:nb
    %
    %         datTemp1 = pethAllB2(:,uu);
    %         datTemp1 = datTemp1(~isnan(datTemp1));
    %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %         erbarsB2(:,uu) = ci_hi_temp1;
    %
    %         datTemp1 = pethAllA2(:,uu);
    %         datTemp1 = datTemp1(~isnan(datTemp1));
    %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %         erbarsA2(:,uu) = ci_hi_temp1;
    %
    %     end
    %
    %     axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %     plot(taxB,nanmean(pethAllB2),'k-','linewidth',2);
    %     plot(taxF,nanmean(pethAllA2),'k-','linewidth',2);
    %     plot(taxB,erbarsB2(1,:),'k:','linewidth',1); plot(taxB,erbarsB2(2,:),'k:','linewidth',1);
    %     plot(taxF,erbarsA2(1,:),'k:','linewidth',1); plot(taxF,erbarsA2(2,:),'k:','linewidth',1);
    %     % legend('LC FR','ACC \it r_s_c','autoupdate','off','fontsize',8); legend('boxoff');
    %     % xlabel('time re:start of stable fixation (ms)');
    %     ylabel('ACC firing rate (sp/s)');
    %     % plot([tax(1) tax(end)],[0 0],'k--');
    %     % minY = 1.7; maxY = 2.75;
    %     minY = 2.9; maxY = 3.4;
    %     %     minY = 0; maxY = 8;
    %     axis([taxB(1)-50 taxF(end)+50 minY maxY]);
    %     plot([0 0],[minY maxY],'k--','linewidth',0.5);
    
    
    %% rsc
    
    % pause
    
    peth1 = dat_trial_rscBnz{1}; peth2 = dat_trial_rscBnz{2}; pethAllBNZ = [peth1;peth2];
    peth1 = dat_trial_rscBz{1}; peth2 = dat_trial_rscBz{2}; pethAllBZ = [peth1;peth2];
    peth1 = dat_trial_rscAnz{1}; peth2 = dat_trial_rscAnz{2}; pethAllANZ = [peth1;peth2];
    peth1 = dat_trial_rscAz{1}; peth2 = dat_trial_rscAz{2}; pethAllAZ = [peth1;peth2];
    
    gibnz = find(isfinite(sum(pethAllBNZ,2))); gibz = find(isfinite(sum(pethAllBZ,2)));
    gianz = find(isfinite(sum(pethAllANZ,2))); giaz = find(isfinite(sum(pethAllAZ,2)));
    
    gi = find(isfinite(sum(pethAllBNZ,2)) & isfinite(sum(pethAllBZ,2)) & isfinite(sum(pethAllANZ,2)) & isfinite(sum(pethAllAZ,2)));
    
    nb = size(pethAllBNZ,2);
    erbars_bnz = nans(2,nb); erbars_bz = nans(2,nb); erbars_anz = nans(2,nb); erbars_az = nans(2,nb);
    
    for uu = 1:nb
        
        datTemp1 = pethAllBNZ(gibnz,uu);
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
        erbars_bnz(:,uu) = ci_hi_temp1;
        
        datTemp1 = pethAllBZ(gibz,uu);
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
        erbars_bz(:,uu) = ci_hi_temp1;
        
        datTemp1 = pethAllANZ(gianz,uu);
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
        erbars_anz(:,uu) = ci_hi_temp1;
        
        datTemp1 = pethAllAZ(giaz,uu);
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
        erbars_az(:,uu) = ci_hi_temp1;
        
    end
    
    
    %% Plot rsc.
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    plot(taxB,nanmean(pethAllBZ(gibz,:)),'-','color',[0.5 0.5 0.5],'linewidth',2);
    plot(taxB,nanmean(pethAllBNZ(gibnz,:)),'k-','linewidth',2);
    legend('LC=0','LC>0','autoupdate','off','fontsize',8); legend('boxoff');
    
    plot(taxF,nanmean(pethAllAZ(giaz,:)),'-','color',[0.5 0.5 0.5],'linewidth',2);
    plot(taxF,nanmean(pethAllANZ(gianz,:)),'k-','linewidth',2);
    
    plot(taxB,erbars_bnz(1,:),'k:','linewidth',1); plot(taxB,erbars_bnz(2,:),'k:','linewidth',1);
    plot(taxF,erbars_anz(1,:),'k:','linewidth',1); plot(taxF,erbars_anz(2,:),'k:','linewidth',1);
    
    plot(taxB,erbars_bz(1,:),'--','color',[0.5 0.5 0.5],'linewidth',0.5); plot(taxB,erbars_bz(2,:),'--','color',[0.5 0.5 0.5],'linewidth',0.5);
    plot(taxF,erbars_az(1,:),'--','color',[0.5 0.5 0.5],'linewidth',0.5); plot(taxF,erbars_az(2,:),'--','color',[0.5 0.5 0.5],'linewidth',0.5);
    
    % errorBarFilled(tax',nanmedian(pethAll),nanse(pethAll,1),[0.7 0.7 0.7]);
    % plot(tax,pethAllQ(1,:),'k--'); plot(tax,pethAllQ(2,:),'k--');
    xlabel('time re:start of stable fixation (ms)');
    ylabel('ACC \it r_s_c');
    minY = 0.02; maxY = 0.062;
    axis([taxB(1)-50 taxF(end)+50 minY maxY]);
    % axis([tax(1)-100 tax(end)+100 -0.15 0.15]);
    plot([0 0],[minY maxY],'k--','linewidth',0.5);
    
    % pause
    
    
    %% CCG between LC FR and ACC FR.
    
    %     %     peth1 = dat_trial_ccgB2{1};
    %     %     peth2 = dat_trial_ccgB2{2};
    %     %     pethAllB = [peth1;peth2];
    %     %     giB = find(isfinite(sum(pethAllB,2)));
    %     %     pethAllB = pethAllB(giB,:);
    %     %
    %     %     peth1 = dat_trial_ccgA2{1};
    %     %     peth2 = dat_trial_ccgA2{2};
    %     %     pethAllA = [peth1;peth2];
    %     %     giA = find(isfinite(sum(pethAllA,2)));
    %     %     pethAllA = pethAllA(giA,:);
    %     %
    %     %     nb = size(pethAllB,2);
    %     %
    %     %     erbarsB = nans(2,nb);
    %     %     erbarsA = nans(2,nb);
    %     %
    %     %     for uu = 1:nb
    %     %         datTemp1 = pethAllB(:,uu);
    %     %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %     %         erbarsB(:,uu) = ci_hi_temp1;
    %     %
    %     %         datTemp1 = pethAllA(:,uu);
    %     %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %     %         erbarsA(:,uu) = ci_hi_temp1;
    %     %     end
    %     %
    %     %     axes(axs(4)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %     %     plot(tax2,nanmean(pethAllB),'k--','linewidth',2);
    %     %     plot(tax2,nanmean(pethAllA),'k-','linewidth',2);
    %     %     legend('Before stable fixation','After start of stable fixation','autoupdate','off','fontsize',8); legend('boxoff');
    %     %
    %     %     plot(tax2,erbarsB(1,:),'k:','linewidth',.5);
    %     %     plot(tax2,erbarsB(2,:),'k:','linewidth',.5);
    %     %
    %     %     plot(tax2,erbarsA(1,:),'k:','linewidth',.5);
    %     %     plot(tax2,erbarsA(2,:),'k:','linewidth',.5);
    %     %
    %     %     xlabel({'time lag (ms)';'ACC rsc re: LC spike rate'}); ylabel('CCG');
    %     %     minY = 5; maxY = 7.5;
    %     %     axis([tax2(1)-100 tax2(end)+100 minY maxY]);
    %     %     plot([0 0],[minY maxY],'k--','linewidth',0.5);
    
    %% CCG between LC FR and ACC rsc
    
    %     peth1 = dat_trial_ccgB{1};
    %     peth2 = dat_trial_ccgB{2};
    %     pethAllB = [peth1;peth2];
    %     giB = find(isfinite(sum(pethAllB,2)));
    %     pethAllB = pethAllB(giB,:);
    %
    %     peth1 = dat_trial_ccgA{1};
    %     peth2 = dat_trial_ccgA{2};
    %     pethAllA = [peth1;peth2];
    %     giA = find(isfinite(sum(pethAllA,2)));
    %     pethAllA = pethAllA(giA,:);
    %
    %     nb = size(pethAllB,2);
    %     %     nu = size(pethAllB,1);
    %
    %     %     % For cross-correlations:
    %     %     sigV = 2./(sqrt(length(tax)-abs([-20:1:20])));
    %     %     sigOut = zeros(nu,nb);
    %     %     for uu = 1:nu
    %     %         sigOutTemp = zeros(1,nb);
    %     %         datTemp1 = pethAll(uu,:);
    %     %         sigOutTemp(find(datTemp1>sigV))=1;
    %     %         sigOutTemp(find(isnan(datTemp1)))=nan;
    %     %         sigOut(uu,:) = sigOutTemp;
    %     %     end
    %
    %     %     figure; hold on;
    %     %     plot(tax2,nanmean(sigOut));
    %
    %     %     smDat = nans(nu,nb);
    %     erbarsB = nans(2,nb);
    %     erbarsA = nans(2,nb);
    %
    %     for uu = 1:nb
    %         datTemp1 = pethAllB(:,uu);
    %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %         erbarsB(:,uu) = ci_hi_temp1;
    %
    %         datTemp1 = pethAllA(:,uu);
    %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %         erbarsA(:,uu) = ci_hi_temp1;
    %     end
    %
    %     %% Plot CCG.
    %
    %     %tax2 = [-1750 tax2 1750]
    %     % axes(axs(3)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %     axes(axs(3)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %     plot(tax2,nanmean(pethAllB),'k--','linewidth',2);
    %     plot(tax2,nanmean(pethAllA),'k-','linewidth',2);
    %     legend('Before stable fixation','After start of stable fixation','autoupdate','off','fontsize',8); legend('boxoff');
    %
    %     plot(tax2,erbarsB(1,:),'k:','linewidth',.5);
    %     plot(tax2,erbarsB(2,:),'k:','linewidth',.5);
    %
    %     plot(tax2,erbarsA(1,:),'k:','linewidth',.5);
    %     plot(tax2,erbarsA(2,:),'k:','linewidth',.5);
    %
    %     % errorBarFilled(tax2,nanmedian(pethAll),nanse(pethAll,1),[0.7 0.7 0.7]);
    %     % plot(tax2,pethAllQ(1,:),'k--'); plot(tax2,pethAllQ(2,:),'k--');
    %     % plot([tax2(1) tax2(end)],[0 0],'k--');
    %     xlabel({'time lag (ms)';'ACC rsc re: LC spike rate'}); ylabel('CCG');
    %     % minY = 0.04; maxY = 0.08;
    %     minY = 0.05; maxY = 0.11;
    %     axis([tax2(1)-100 tax2(end)+100 minY maxY]);
    %     plot([0 0],[minY maxY],'k--','linewidth',0.5);
    %
    %     %     sigOutAll = zeros(1,nb);
    %     %     sigOutAll(find(datTemp1>sigV))=1;
    %     %     plot(tax2(find(sigOutAll)),maxY,'k*');
    %     %
    %     %     for uu = 1:length(tax2)
    %     %         % p(uu) = signrank(smDat(:,uu));
    %     %         p(uu) = signrank(pethAll(:,uu));
    %     %         if p(uu) < 0.05
    %     %             plot(tax2(uu),maxY,'k*');
    %     %         end
    %     %     end
    %
    %     %     axis([-1500 2000 0 4.5]);
    %     %     text(1000,3,'median w. IQR');
    %     %     text(2000,5,'ACC trial firing rate','fontweight','bold','fontsize',10);
    %
    %     %     axes(axs(4)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %     %     plot(tax2,nanmean(pethAll),'k-'); errorBarFilled(tax2,nanmean(pethAll),nanse(pethAll,1),[0.7 0.7 0.7]);
    %     %     %     axis([-1500 1500 2.9 3.4]);
    %     %     xlabel('time re:start of stable fixation'); ylabel('firing rate (sp/s)');
    %     %     %     text(1000,3.3,'mean+-sem');
    %
    %     %%
    %
    %     %     axes(axs(4));
    %     %     set(gca,'visible','off');
    %
    %     %% Find peaks and troughs in CCG.
    %
    %     %     ci = 0;
    %     %     ci2 = 0;
    %     %     maxVals = nans(1000000,1);
    %     %     minVals = nans(1000000,1);
    %     %     pind = find(tax2>0)
    %     %     tax2p = tax2(pind);
    %     %     for uu = 1:size(pethAll,1)
    %     %         datTemp1 = pethAll(uu,pind);
    %     %         tmaxVal = tax2p(find((smooth(datTemp1,3))>5));
    %     %         tminVal = tax2p(find((smooth(datTemp1,3))<-5));
    %     %         if ~isempty(tmaxVal)
    %     %             nvals = length(tmaxVal);
    %     %             maxVals(ci+1:ci+nvals) = tmaxVal;
    %     %             ci = ci+nvals;
    %     %         end
    %     %         if ~isempty(tminVal)
    %     %             nvals = length(tminVal);
    %     %             minVals(ci2+1:ci2+nvals) = tminVal;
    %     %             ci2 = ci2+nvals;
    %     %         end
    %     %     end
    %     %
    %     %     maxVals(ci+1:end) = [];
    %     %     minVals(ci2+1:end) = [];
    %     %
    %     %     axes(axs(4)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %     %     bar(bc,nn/sum(nn),'facecolor',[0 0 0],'edgecolor','none');
    %     %     bar(bc,-nn2/sum(nn2),'facecolor',[0.7 0.7 0.7],'edgecolor','none');
    
end








