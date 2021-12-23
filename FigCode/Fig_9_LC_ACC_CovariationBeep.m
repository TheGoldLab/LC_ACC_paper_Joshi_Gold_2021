% Fig_9_LC_ACC_CovariationBeep.m
%
% INFO: Compare LC spiking and ACC rsc time courses.
% BEEP trials only.
% Previous version: LC (ACC) STA of ACC (LC) spiking, Fano and pairwise spike count correlation (rsc).
% Dual area recording - figure script; FIXATION task.

% Origin: 010519.
%
% For paper ... Combined LC_Dual_Recording_sts.m and ACC_Dual_Recording_sts.m
%
% OLD: Plots ACC (LC) spike rate aligned to LC (ACC) spiking.
% OLD: Plots ACC (LC) pairwise spike count correlations (rsc); ACC (LC) pairs with ACC (LC) unit spikes aligned to LC (ACC) spikes.
%
% Origin of "LC_Dual_Recording_sts.m": 111016 - Sidd.
% Modified: 062116, 070216, 072518, 092618.
% Mod: 112918: Cleaning and finalizing - Sidd.
% Calls "getLC_DualArea_sts" to do the heavy lifting.
% 010519: "finalizing" .... ha ha. Nope. Still not funny.

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
% Reanalyze data?
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
                    %                     pause
                    CP_dat{ff} = getLC_DualArea_CovariationBeep(siteData); % Call "getLC_DualArea_sts" to do heavy lifting.
                    % pause
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
    
    % save LC_ACC_CCG_091221_Beep mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    % save LC_ACC_CCG_091921_Beep_d mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    % save LC_ACC_CCG_092021_Beep_a mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    % save LC_ACC_CCG_100821_Beep_a mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    save LC_ACC_CCG_101021_Beep mnkRes; % LC spikes form 1 s to 2.1 sec after stable fixation starts.
    
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
    % inFile = 'LC_ACC_CCG_091221_Beep';
    % inFile = 'LC_ACC_CCG_091821_Beepb';
    % inFile = 'LC_ACC_CCG_091921_Beep_d';
    
    % inFile = 'LC_ACC_CCG_092021_Beep_a';
    inFile = 'LC_ACC_CCG_101021_Beep';
    
    
    smW = 1;
    base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.
    savedir = strcat([base_dir,'\Results\Results_2021\LC_ACC\']);
    cd(savedir);
    load(inFile);
    
    % **********************************************
    
    monks = {'Sprout','Cicero'}; % Add mnks as needed...
    sites = {'LC_ACC_Fixation'}; % Add sites as needed...
    nMonks = length(monks); nSites = length(sites);
    
    dat_trial_FRBp = cell(1,2); dat_trial_FRAp = cell(1,2);
    dat_trial_rscBp = cell(1,2); dat_trial_rscAp = cell(1,2);
    dat_trial_ccgBp = cell(1,2); dat_trial_ccgAp = cell(1,2);
    
    dat_trial_FRBnp = cell(1,2); dat_trial_FRAnp = cell(1,2);
    dat_trial_rscBnp = cell(1,2); dat_trial_rscAnp = cell(1,2);
    dat_trial_ccgBnp = cell(1,2); dat_trial_ccgAnp = cell(1,2);
    
    % First - get the processed data:
    
    % ccg_dat = {LC_dat taxF taxB tax2F TT};
    % LC_dat = {out_ccgBp out_ccgAp out_ccgBnp out_ccgAnp out_rsc_ABp out_rsc_AAp out_rsc_ABnp out_rsc_AAnp out_FR_Bp out_FR_Ap out_FR_Bnp out_FR_Anp};
    
    for mm = 1:nMonks % Loop over mnks.
        mnkDat = mnkRes{mm};
        for ss = 1:nSites % Loop over brainstem sites.
            sitDat = mnkDat{ss};
            for ap = 1:length(sitDat) % Iterate over sessions.
                nLCU = length(sitDat{ap}{1});
                for li = 1:nLCU
                    
                    dat_trial_ccgBp{mm} = [dat_trial_ccgBp{mm}; sitDat{ap}{1}{li}{1}];
                    dat_trial_ccgAp{mm} = [dat_trial_ccgAp{mm}; sitDat{ap}{1}{li}{2}];
                    
                    dat_trial_ccgBnp{mm} = [dat_trial_ccgBnp{mm}; sitDat{ap}{1}{li}{3}];
                    dat_trial_ccgAnp{mm} = [dat_trial_ccgAnp{mm}; sitDat{ap}{1}{li}{4}];
                    
                    dat_trial_rscBp{mm} = [dat_trial_rscBp{mm}; sitDat{ap}{1}{li}{5}];
                    dat_trial_rscAp{mm} = [dat_trial_rscAp{mm}; sitDat{ap}{1}{li}{6}];
                    
                    dat_trial_rscBnp{mm} = [dat_trial_rscBnp{mm}; sitDat{ap}{1}{li}{7}];
                    dat_trial_rscAnp{mm} = [dat_trial_rscAnp{mm}; sitDat{ap}{1}{li}{8}];
                    
                    dat_trial_FRBp{mm} = [dat_trial_FRBp{mm}; sitDat{ap}{1}{li}{9}];
                    dat_trial_FRAp{mm} = [dat_trial_FRAp{mm}; sitDat{ap}{1}{li}{10}];
                    
                    dat_trial_FRBnp{mm} = [dat_trial_FRBnp{mm}; sitDat{ap}{1}{li}{11}];
                    dat_trial_FRAnp{mm} = [dat_trial_FRAnp{mm}; sitDat{ap}{1}{li}{12}];
                    
                    %                 all_rsc = sitDat{ap}{1};
                    %                 for tt = 1:length(all_rsc)
                    %                     dat_trial_rscA{mm} = [dat_trial_rscA{mm}; all_rsc{tt}'];
                    %                 end
                    %                 all_FR = (sitDat{ap}{2});
                    %                 dat_trial_FRA{mm} = [dat_trial_FRA{mm}; all_FR];
                end % LC unit loop.
            end % Session loop.
        end % Sites loop.
    end % Mnks loop.
    
    %%
    % ccg_dat = {LC_dat taxF taxB tax2F TT};
    taxA =  sitDat{1}{2};
    taxB =  sitDat{1}{3};
    tax2A =  sitDat{1}{4}; % tax2 = tax2(2:end-1);
    % tax2B =  sitDat{ap}{8}; % tax2 = tax2(2:end-1);
    TT =  sitDat{1}{5};
    
    %     out_ccg = nans(1000,length(tax2));
    %     dat_rsc_all = [dat_trial_rsc{1};dat_trial_rsc{2}]; nrsc = size(dat_rsc_all,1);
    %     dat_FR_all = [dat_trial_FR{1};dat_trial_FR{2}]; nFR = size(dat_FR_all,1);
    %     for tt = 1:1000
    %         raw_ccg = xcorr(zscore(dat_rsc_all(1+floor(rand(1)*nrsc),:)),zscore(dat_FR_all(1+floor(rand(1)*nFR),:)));
    %         out_ccg(tt,:) = raw_ccg./TT;
    %     end
    
    
    %% Setup figure:
    figureNumber = 216; num = 216; wid = 17.6; hts = [8]; cols = {3 3}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 2, 1, [12], '', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    %%
    
    %     pause
    
    %% FR
    
    peth1Ap = dat_trial_FRAp{1}; peth2Ap = dat_trial_FRAp{2}; pethAllAp = [peth1Ap;peth2Ap];
    peth1Bp = dat_trial_FRBp{1}; peth2Bp = dat_trial_FRBp{2}; pethAllBp = [peth1Bp;peth2Bp];
    
    peth1Anp = dat_trial_FRAnp{1}; peth2Anp = dat_trial_FRAnp{2}; pethAllAnp = [peth1Anp;peth2Anp];
    peth1Bnp = dat_trial_FRBnp{1}; peth2Bnp = dat_trial_FRBnp{2}; pethAllBnp = [peth1Bnp;peth2Bnp];
    
    nb = size(pethAllAp,2);
    erbarsAp = nans(2,nb); erbarsBp = nans(2,nb);
    erbarsAnp = nans(2,nb); erbarsBnp = nans(2,nb);
    % bind = (nb/2):(nb-1);
    %     bind = ((nb/2)-1):(nb);
    %     aind = 1:5;
    
    %     bind = (nb-5):(nb);
    %     aind = 1:5;
    
    bind = find(taxB>-500); % Points just before beep. Essentially -500 to 0 ms re: beep.
    aind = find(taxA<500); % Points just after beep. Essentially 0-500 ms re: beep.
    late_ind = find(taxA>500); % Points for "late" part of post-beep activity. Essentially 500-1000 ms re: beep.
    
    for uu = 1:nb
        
        datTemp1 = pethAllBp(:,uu); datTemp11 = datTemp1(~isnan(datTemp1));
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp11},'alpha',0.05,'type','per');
        erbarsBp(:,uu) = ci_hi_temp1;
        
        datTemp1 = pethAllAp(:,uu); datTemp12 = datTemp1(~isnan(datTemp1));
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp12},'alpha',0.05,'type','per');
        erbarsAp(:,uu) = ci_hi_temp1;
        
        allBp = (pethAllBp(:,bind));
        allBp = allBp(:);
        allBp = allBp(isfinite(allBp));
        p_p(uu) = ranksum(datTemp12,allBp);
        
        datTemp1 = pethAllBnp(:,uu); datTemp13 = datTemp1(~isnan(datTemp1));
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp13},'alpha',0.05,'type','per');
        erbarsBnp(:,uu) = ci_hi_temp1;
        
        datTemp1 = pethAllAnp(:,uu); datTemp14 = datTemp1(~isnan(datTemp1));
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp14},'alpha',0.05,'type','per');
        erbarsAnp(:,uu) = ci_hi_temp1;
        
        allBnp = (pethAllBnp(:,bind));
        allBnp = allBnp(:);
        allBnp = allBnp(isfinite(allBnp));
        p_np(uu) = ranksum(datTemp14,allBnp);
        
        % pB(uu) = ranksum(datTemp11,datTemp13);
        % pA(uu) = ranksum(datTemp12,datTemp14);
        
    end
    
    bef_p_dat = nanmean(pethAllBp(:,bind),2); % gip = find(isfinite(bef_p_dat));
    bef_np_dat = nanmean(pethAllBnp(:,bind),2); % ginp = find(isfinite(bef_np_dat));
    giB = find(isfinite(bef_np_dat) & isfinite(bef_p_dat));
    pB = ranksum(bef_p_dat(giB),bef_np_dat(giB)); % 0.0143
    
    aft_p_dat = nanmean(pethAllAp(:,aind),2); % gip = find(isfinite(bef_p_dat));
    aft_np_dat = nanmean(pethAllAnp(:,aind),2); % ginp = find(isfinite(bef_np_dat));
    giA = find(isfinite(aft_np_dat) & isfinite(aft_p_dat));
    pA = ranksum(aft_p_dat(giA),aft_np_dat(giA)); % 7.1174 10^-22
    
    % pause
    
    %% Phasic response in LC.
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    tax_inc = unique(diff(taxA));
    ydat = nanmean(pethAllAp); xdat = taxA; x_int = [taxA(1):(tax_inc/4):taxA(end)];
    interp_ydat = spline(xdat,ydat,x_int);
    for pii = 1:nb
        if pii==1
            if p_p(pii)<0.05      plot([x_int(1:3)],[interp_ydat(1:3)],'r-','linewidth',2); end
            if p_p(pii)>=0.05    plot([x_int(1:3)],[interp_ydat(1:3)],'k-','linewidth',2); end
        end
        if pii>1 & pii<nb
            si = -2+find(x_int== taxA(pii)); ei = 2+find(x_int == taxA(pii));
            if p_p(pii)<0.05      plot(x_int(si:ei),interp_ydat(si:ei),'r-','linewidth',2); end
            if p_p(pii)>=0.05    plot(x_int(si:ei),interp_ydat(si:ei),'k-','linewidth',2); end
        end
        % pause
        if pii==nb
            if p_p(pii)<0.05      plot([x_int((end-3):end)],[interp_ydat((end-3):end)],'r-','linewidth',2); end
            if p_p(pii)>=0.05    plot([x_int((end-3):end)],[interp_ydat((end-3):end)],'k-','linewidth',2); end
        end
    end
    
    plot(taxB,nanmean(pethAllBp),'k-','linewidth',2);
    % plot(taxA,nanmean(pethAllAp),'k-','linewidth',2);
    % % legend('before','after','autoupdate','off','fontsize',8); legend('boxoff');
    
    plot(taxA,erbarsAp(1,:),'k:','linewidth',1); plot(taxA,erbarsAp(2,:),'k:','linewidth',1);
    plot(taxB,erbarsBp(1,:),'k:','linewidth',1); plot(taxB,erbarsBp(2,:),'k:','linewidth',1);
    
    xlabel('time re:beep (s)'); ylabel('LC firing rate (sp/s)');
    % plot([taxB(1) taxA(end)],[0 0],'k--');
    plot([0 0],[0 8],'k--','linewidth',0.5);
    minY = 0; maxY = 8;
    axis([taxB(1)-100 taxA(end)+100 0 8]);
    
    % sigi = find(p_p<0.05);
    % plot(taxA(sigi),maxY,'k*');
    % bar(taxA(sigi),7*ones(1,length(sigi)),'r','facealpha',0.2,'edgecolor','r','edgealpha',0.2);
    
    %     sigi = find(pB<0.05); plot(taxB(sigi),maxY,'k*');
    %     sigi = find(pA<0.05); plot(taxA(sigi),maxY,'k*');
    
    si = bind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = bind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxB(si)-tax_inc/2 taxB(si)-tax_inc/2 taxB(ei)+tax_inc/2 taxB(ei)+tax_inc/2], [0 maxY maxY 0],'k','facealpha',0.2,'edgecolor','none');
    si = aind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = aind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxA(si)-tax_inc/2 taxA(si)-tax_inc/2 taxA(ei)+tax_inc/2 taxA(ei)+tax_inc/2], [0 maxY maxY 0],'g','facealpha',0.2,'edgecolor','none');
    si = late_ind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = late_ind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxA(si)-tax_inc/2 taxA(si)-tax_inc/2 taxA(ei)+tax_inc/2 taxA(ei)+tax_inc/2], [0 maxY maxY 0],'k','facealpha',0.2,'edgecolor','none');
    
    %% No Phasic response in LC.
    
    axes(axs(4)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    tax_inc = unique(diff(taxA));
    ydat = nanmean(pethAllAnp); xdat = taxA; x_int = [taxA(1):(tax_inc/4):taxA(end)];
    interp_ydat = spline(xdat,ydat,x_int);
    for pii = 1:nb
        if pii==1
            if p_np(pii)<0.05      plot([x_int(1:3)],[interp_ydat(1:3)],'r-','linewidth',2); end
            if p_np(pii)>=0.05    plot([x_int(1:3)],[interp_ydat(1:3)],'k-','linewidth',2); end
        end
        if pii>1 & pii<nb
            si = -2+find(x_int== taxA(pii)); ei = 2+find(x_int == taxA(pii));
            if p_np(pii)<0.05      plot(x_int(si:ei),interp_ydat(si:ei),'r-','linewidth',2); end
            if p_np(pii)>=0.05    plot(x_int(si:ei),interp_ydat(si:ei),'k-','linewidth',2); end
        end
        % pause
        if pii==nb
            if p_np(pii)<0.05      plot([x_int((end-3):end)],[interp_ydat((end-3):end)],'r-','linewidth',2); end
            if p_np(pii)>=0.05    plot([x_int((end-3):end)],[interp_ydat((end-3):end)],'k-','linewidth',2); end
        end
    end
    
    plot(taxB,nanmean(pethAllBnp),'k-','linewidth',2);
    % plot(taxA,nanmean(pethAllAp),'k-','linewidth',2);
    % % legend('before','after','autoupdate','off','fontsize',8); legend('boxoff');
    
    plot(taxA,erbarsAnp(1,:),'k:','linewidth',1); plot(taxA,erbarsAnp(2,:),'k:','linewidth',1);
    plot(taxB,erbarsBnp(1,:),'k:','linewidth',1); plot(taxB,erbarsBnp(2,:),'k:','linewidth',1);
    
    xlabel('time re:beep (s)'); ylabel('LC firing rate (sp/s)');
    % plot([taxB(1) taxA(end)],[0 0],'k--');
    plot([0 0],[0 8],'k--','linewidth',0.5);
    minY = 0; maxY = 8;
    axis([taxB(1)-100 taxA(end)+100 0 8]);
    
    % sigi = find(p_np<0.05);
    % plot(taxA(sigi),maxY,'k*');
    % bar(taxA(sigi),7*ones(1,length(sigi)),'r','facealpha',0.2,'edgecolor','r','edgealpha',0.2);
    
    %     sigi = find(pB<0.05); plot(taxB(sigi),maxY,'k*');
    %     sigi = find(pA<0.05); plot(taxA(sigi),maxY,'k*');
    
    si = bind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = bind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxB(si)-tax_inc/2 taxB(si)-tax_inc/2 taxB(ei)+tax_inc/2 taxB(ei)+tax_inc/2], [0 maxY maxY 0],'k','facealpha',0.2,'edgecolor','none');
    si = aind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = aind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxA(si)-tax_inc/2 taxA(si)-tax_inc/2 taxA(ei)+tax_inc/2 taxA(ei)+tax_inc/2], [0 maxY maxY 0],'g','facealpha',0.2,'edgecolor','none');
    si = late_ind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = late_ind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxA(si)-tax_inc/2 taxA(si)-tax_inc/2 taxA(ei)+tax_inc/2 taxA(ei)+tax_inc/2], [0 maxY maxY 0],'k','facealpha',0.2,'edgecolor','none');
    
    %% ACC rsc.
    
    % pause
    
    peth1Ap = dat_trial_rscAp{1}; peth2Ap = dat_trial_rscAp{2}; pethAllAp = [peth1Ap;peth2Ap];
    peth1Bp = dat_trial_rscBp{1}; peth2Bp = dat_trial_rscBp{2}; pethAllBp = [peth1Bp;peth2Bp];
    
    peth1Anp = dat_trial_rscAnp{1}; peth2Anp = dat_trial_rscAnp{2}; pethAllAnp = [peth1Anp;peth2Anp];
    peth1Bnp = dat_trial_rscBnp{1}; peth2Bnp = dat_trial_rscBnp{2}; pethAllBnp = [peth1Bnp;peth2Bnp];
    
    gi = find(isfinite(sum(pethAllAp,2)) & isfinite(sum(pethAllBp,2)) & isfinite(sum(pethAllAnp,2)) & isfinite(sum(pethAllBnp,2)));
    % gi = find(isfinite(sum(pethAllAp,2)) & isfinite(sum(pethAllBp,2)) & isfinite(sum(pethAllAnp,2)) & isfinite(sum(pethAllBnp,2)));
    
    %     diff_before = [pethAllBp(gi,:)-pethAllBnp(gi,:)];
    %     diff_after = [pethAllAp(gi,:)-pethAllAnp(gi,:)];
    
    nb = size(pethAllAp,2);
    erbarsAp = nans(2,nb); erbarsBp = nans(2,nb);
    erbarsAnp = nans(2,nb); erbarsBnp = nans(2,nb);
    % erbarsDiffB = nans(2,nb); erbarsDiffA = nans(2,nb);
    
    %     bind = (nb-5):(nb);
    %     aind = 1:5;
    %     bind = find(taxB>-500);
    %     aind = find(taxA<500);
    %     late_ind = find(taxA>500);
    
    for uu = 1:nb
        
        datTemp1 = pethAllBp(gi,uu); datTemp11 = datTemp1(~isnan(datTemp1));
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp11},'alpha',0.05,'type','per');
        erbarsBp(:,uu) = ci_hi_temp1;
        
        datTemp1 = pethAllAp(gi,uu); datTemp12 = datTemp1(~isnan(datTemp1));
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp12},'alpha',0.05,'type','per');
        erbarsAp(:,uu) = ci_hi_temp1;
        
        allBp = pethAllBp(:,bind);
        allBp = allBp(:);
        allBp = allBp(isfinite(allBp));
        p_p(uu) = ranksum(datTemp12,allBp);
        
        datTemp1 = pethAllBnp(gi,uu); datTemp13 = datTemp1(~isnan(datTemp1));
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp13},'alpha',0.05,'type','per');
        erbarsBnp(:,uu) = ci_hi_temp1;
        
        datTemp1 = pethAllAnp(gi,uu); datTemp14 = datTemp1(~isnan(datTemp1));
        [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp14},'alpha',0.05,'type','per');
        erbarsAnp(:,uu) = ci_hi_temp1;
        
        allBnp = pethAllBnp(:,bind);
        allBnp = allBnp(isfinite(allBnp));
        allBnp = allBnp(:);
        p_np(uu) = ranksum(datTemp14,allBnp);
        
        % pB(uu) = ranksum(datTemp11,datTemp13);
        % pA(uu) = ranksum(datTemp12,datTemp14);
        
        % differences
        
        %         datTemp1 = diff_before(:,uu); datTemp11 = datTemp1(~isnan(datTemp1));
        %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp11},'alpha',0.05,'type','per');
        %         erbarsDiffB(:,uu) = ci_hi_temp1;
        %
        %         datTemp1 = diff_after(:,uu); datTemp11 = datTemp1(~isnan(datTemp1));
        %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp11},'alpha',0.05,'type','per');
        %         erbarsDiffA(:,uu) = ci_hi_temp1;
        
    end
    
    %     bef_p_dat = nanmean(pethAllBp(:,bind),2); % gip = find(isfinite(bef_p_dat));
    %     bef_np_dat = nanmean(pethAllBnp(:,bind),2); % ginp = find(isfinite(bef_np_dat));
    %     giB = find(isfinite(bef_np_dat) & isfinite(bef_p_dat));
    %     pB = ranksum(bef_p_dat(giB),bef_np_dat(giB)); % 0.2083
    %
    %     aft_p_dat = nanmean(pethAllAp(:,aind),2); % gip = find(isfinite(bef_p_dat));
    %     aft_np_dat = nanmean(pethAllAnp(:,aind),2); % ginp = find(isfinite(bef_np_dat));
    %     giA = find(isfinite(aft_np_dat) & isfinite(aft_p_dat));
    %     pA = ranksum(aft_p_dat(giA),aft_np_dat(giA)); % 0.3639
    
    %     late_ind = find(taxA>500);
    aft_p_dat = pethAllAp(:,late_ind); aft_p_dat = aft_p_dat(:); % gip = find(isfinite(bef_p_dat));
    aft_np_dat = pethAllAnp(:,late_ind); aft_np_dat = aft_np_dat(:)% ginp = find(isfinite(bef_np_dat));
    giA = find(isfinite(aft_np_dat) & isfinite(aft_p_dat));
    pA = ranksum(aft_p_dat(giA),aft_np_dat(giA)); % 0.0350
    
    %     e_ind = find(taxB>-500);
    bef_p_dat = pethAllBp(:,bind); bef_p_dat = bef_p_dat(:); % gip = find(isfinite(bef_p_dat));
    bef_np_dat = pethAllBnp(:,bind); bef_np_dat = bef_np_dat(:)% ginp = find(isfinite(bef_np_dat));
    giB = find(isfinite(bef_np_dat) & isfinite(bef_p_dat));
    pB = ranksum(bef_p_dat(giB),bef_np_dat(giB)); % 0.0.0693
    
    % pause
    
    %% Slope differences, if any, in ACC rsc, for phasic LC re: no phasic LC response to beep.
    
    % n_points = 8;
    n_points = length(aind);
    % xdat = (taxA(1:n_points))/1000;
    xdat = (taxA(aind))/1000;
    XX = [ones(length(xdat),1) xdat];
    ydatp = pethAllAp(:,1:n_points);
    ydatnp = pethAllAnp(:,1:n_points);
    
    gi = find(isfinite(sum(ydatp,2)) & isfinite(sum(ydatnp,2)));
    
    ydatp = ydatp(gi,:);
    ydatnp = ydatnp(gi,:);
    
    nv = size(ydatp,1);
    pp = nans(1,nv);
    pnp = nans(1,nv);
    diffVals = nans(1,nv);
    
    for ti = 1:nv
        pp_temp = XX\ydatp(ti,:)';
        pnp_temp = XX\ydatnp(ti,:)';
        pp(ti) = pp_temp(2);
        pnp(ti) = pnp_temp(2);
        diffVals(ti) = pp_temp(2) - pnp_temp(2);
    end
    
    pslope = ranksum(pp,pnp); % p = 0.0297
    prcVals = prctile(diffVals,[25 50 75]); %  -0.8573    0.0733    0.9439
    diffVals = diffVals(find(diffVals>prcVals(1) & diffVals<prcVals(3)));
    
    %% Plot slopes.
    
    %     figure; hold on;
    %     [nn bb] = hist(diffVals,15);
    %     bar(bb,nn/sum(nn),'facecolor','none','edgecolor','k');
    %     plot([prcVals(2) prcVals(2)],[0 0.2],'r--');
    %
    %     figure; hold on;
    %     plot(pnp,pp,'ko','markersize',4,'markerfacecolor','none','markeredgecolor','k');
    %     minVal = -2; maxVal = 2;
    %     axis([minVal maxVal minVal maxVal]);
    %     plot([0 0],[minVal maxVal],'r--');
    %     plot([minVal maxVal],[0 0],'r--');
    
    
    % pause
    
    
    %% Test plot difference.
    
    % figure; hold on; ax = gca; disableDefaultInteractivity(ax);
    % plot(taxB,nanmean(diff_before),'k-','linewidth',2);
    % plot(taxA,nanmean(diff_after),'r-','linewidth',2);
    % plot(taxB,erbarsDiffB(1,:),'k:','linewidth',1); plot(taxB,erbarsDiffB(2,:),'k:','linewidth',1);
    % plot(taxA,erbarsDiffA(1,:),'k:','linewidth',1); plot(taxA,erbarsDiffA(2,:),'k:','linewidth',1);
    % plot([taxB(1) taxA(end)],[0 0],'k--', 'linewidth',1);
    
    %% For LC phasic, plot ACC rsc.
    
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    tax_inc = unique(diff(taxA));
    ydat = nanmean(pethAllAp(gi,:)); xdat = taxA; x_int = [taxA(1):(tax_inc/4):taxA(end)];
    interp_ydat = spline(xdat,ydat,x_int);
    for pii = 1:nb
        if pii==1
            if p_p(pii)<0.05      plot([x_int(1:3)],[interp_ydat(1:3)],'r-','linewidth',2); end
            if p_p(pii)>=0.05    plot([x_int(1:3)],[interp_ydat(1:3)],'k-','linewidth',2); end
        end
        if pii>1 & pii<nb
            si = -2+find(x_int== taxA(pii)); ei = 2+find(x_int == taxA(pii));
            if p_p(pii)<0.05      plot(x_int(si:ei),interp_ydat(si:ei),'r-','linewidth',2); end
            if p_p(pii)>=0.05    plot(x_int(si:ei),interp_ydat(si:ei),'k-','linewidth',2); end
        end
        % pause
        if pii==nb
            if p_p(pii)<0.05      plot([x_int((end-3):end)],[interp_ydat((end-3):end)],'r-','linewidth',2); end
            if p_p(pii)>=0.05    plot([x_int((end-3):end)],[interp_ydat((end-3):end)],'k-','linewidth',2); end
        end
    end
    
    plot(taxB,nanmean(pethAllBp(gi,:)),'k-','linewidth',2);
    %     plot(taxA,nanmean(pethAllAp(gi,:)),'k-','linewidth',.5);
    % % legend('before','after','autoupdate','off','fontsize',8); legend('boxoff');
    
    plot(taxA,erbarsAp(1,:),'k:','linewidth',1); plot(taxA,erbarsAp(2,:),'k:','linewidth',1);
    plot(taxB,erbarsBp(1,:),'k:','linewidth',1); plot(taxB,erbarsBp(2,:),'k:','linewidth',1);
    
    xlabel('time re:beep (s)'); ylabel('LC firing rate (sp/s)');
    % plot([taxB(1) taxA(end)],[0 0],'k--');
    plot([0 0],[0 maxY],'k--','linewidth',0.5);
    minY = 0; maxY = 0.055;
    axis([taxB(1)-100 taxA(end)+100 0 maxY]);
    
    % sigi = find(p_p<0.05);
    % plot(taxA(sigi),maxY,'k*');
    % bar(taxA(sigi),7*ones(1,length(sigi)),'r','facealpha',0.2,'edgecolor','r','edgealpha',0.2);
    
    %     sigi = find(pB<0.05); plot(taxB(sigi),maxY,'k*');
    %     sigi = find(pA<0.05); plot(taxA(sigi),maxY,'k*');
    
    si = bind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = bind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxB(si)-tax_inc/2 taxB(si)-tax_inc/2 taxB(ei)+tax_inc/2 taxB(ei)+tax_inc/2], [0 maxY maxY 0],'k','facealpha',0.2,'edgecolor','none');
    si = aind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = aind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxA(si)-tax_inc/2 taxA(si)-tax_inc/2 taxA(ei)+tax_inc/2 taxA(ei)+tax_inc/2], [0 maxY maxY 0],'g','facealpha',0.2,'edgecolor','none');
    si = late_ind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = late_ind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxA(si)-tax_inc/2 taxA(si)-tax_inc/2 taxA(ei)+tax_inc/2 taxA(ei)+tax_inc/2], [0 maxY maxY 0],'k','facealpha',0.2,'edgecolor','none');
    
    
    %% For LC no phasic, plot ACC rsc.
    
    axes(axs(5)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    tax_inc = unique(diff(taxA));
    ydat = nanmean(pethAllAnp(gi,:)); xdat = taxA; x_int = [taxA(1):(tax_inc/4):taxA(end)];
    interp_ydat = spline(xdat,ydat,x_int);
    for pii = 1:nb
        if pii==1
            if p_np(pii)<0.05      plot([x_int(1:3)],[interp_ydat(1:3)],'r-','linewidth',2); end
            if p_np(pii)>=0.05    plot([x_int(1:3)],[interp_ydat(1:3)],'k-','linewidth',2); end
        end
        if pii>1 & pii<nb
            si = -2+find(x_int== taxA(pii)); ei = 2+find(x_int == taxA(pii));
            if p_np(pii)<0.05      plot(x_int(si:ei),interp_ydat(si:ei),'r-','linewidth',2); end
            if p_np(pii)>=0.05    plot(x_int(si:ei),interp_ydat(si:ei),'k-','linewidth',2); end
        end
        % pause
        if pii==nb
            if p_np(pii)<0.05      plot([x_int((end-3):end)],[interp_ydat((end-3):end)],'r-','linewidth',2); end
            if p_np(pii)>=0.05    plot([x_int((end-3):end)],[interp_ydat((end-3):end)],'k-','linewidth',2); end
        end
    end
    
    plot(taxB,nanmean(pethAllBnp(gi,:)),'k-','linewidth',2);
    % plot(taxA,nanmean(pethAllAp),'k-','linewidth',2);
    % % legend('before','after','autoupdate','off','fontsize',8); legend('boxoff');
    
    plot(taxA,erbarsAnp(1,:),'k:','linewidth',1); plot(taxA,erbarsAnp(2,:),'k:','linewidth',1);
    plot(taxB,erbarsBnp(1,:),'k:','linewidth',1); plot(taxB,erbarsBnp(2,:),'k:','linewidth',1);
    
    xlabel('time re:beep (s)'); ylabel('LC firing rate (sp/s)');
    % plot([taxB(1) taxA(end)],[0 0],'k--');
    plot([0 0],[0 maxY],'k--','linewidth',0.5);
    minY = 0; maxY = 0.055;
    axis([taxB(1)-100 taxA(end)+100 0 maxY]);
    
    % sigi = find(p_p<0.05);
    % plot(taxA(sigi),maxY,'k*');
    % bar(taxA(sigi),7*ones(1,length(sigi)),'r','facealpha',0.2,'edgecolor','r','edgealpha',0.2);
    
    %     sigi = find(pB<0.05); plot(taxB(sigi),maxY,'k*');
    %     sigi = find(pA<0.05); plot(taxA(sigi),maxY,'k*');
    
    si = bind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = bind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxB(si)-tax_inc/2 taxB(si)-tax_inc/2 taxB(ei)+tax_inc/2 taxB(ei)+tax_inc/2], [0 maxY maxY 0],'k','facealpha',0.2,'edgecolor','none');
    si = aind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = aind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxA(si)-tax_inc/2 taxA(si)-tax_inc/2 taxA(ei)+tax_inc/2 taxA(ei)+tax_inc/2], [0 maxY maxY 0],'g','facealpha',0.2,'edgecolor','none');
    si = late_ind(1); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    ei = late_ind(end); % Don't ask. MATLAB spazzing when I write this in the patch line below...
    patch([taxA(si)-tax_inc/2 taxA(si)-tax_inc/2 taxA(ei)+tax_inc/2 taxA(ei)+tax_inc/2], [0 maxY maxY 0],'k','facealpha',0.2,'edgecolor','none');
    
    
    
    %% CCG
    
    %     % pause
    %
    %     peth1A = dat_trial_ccgAp{1}; peth2A = dat_trial_ccgAp{2}; pethAllAp = [peth1A;peth2A];
    %     peth1B = dat_trial_ccgBp{1}; peth2B = dat_trial_ccgBp{2}; pethAllBp = [peth1B;peth2B];
    %
    %     peth1A = dat_trial_ccgAnp{1}; peth2A = dat_trial_ccgAnp{2}; pethAllAnp = [peth1A;peth2A];
    %     peth1B = dat_trial_ccgBnp{1}; peth2B = dat_trial_ccgBnp{2}; pethAllBnp = [peth1B;peth2B];
    %
    %     gip = find(isfinite(sum(pethAllAp,2)) & isfinite(sum(pethAllBp,2)));
    %     ginp = find(isfinite(sum(pethAllAnp,2)) & isfinite(sum(pethAllBnp,2)));
    %
    %     diff_p = pethAllAp(gip,:) - pethAllBp(gip,:);
    %     diff_np = pethAllAnp(ginp,:) - pethAllBnp(ginp,:);
    %
    %     nb = size(pethAllAp,2); % nu = size(pethAllAp,1);
    %     erbarsAp = nans(2,nb); erbarsBp = nans(2,nb);
    %
    %     erbarsDiffp = nans(2,nb); erbarsDiffnp = nans(2,nb);
    %
    %     for uu = 1:nb
    %
    %         datTemp1 = pethAllAp(gip,uu); % gi = find(~isnan(datTemp1)); datTemp1 = datTemp1(gi);
    %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %         erbarsAp(:,uu) = ci_hi_temp1;
    %
    %         datTemp1 = pethAllBp(gip,uu); % gi = find(~isnan(datTemp1)); datTemp1 = datTemp1(gi);
    %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %         erbarsBp(:,uu) = ci_hi_temp1;
    %
    %         datTemp1 = pethAllAnp(ginp,uu); % gi = find(~isnan(datTemp1)); datTemp1 = datTemp1(gi);
    %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %         erbarsAnp(:,uu) = ci_hi_temp1;
    %
    %         datTemp1 = pethAllBnp(ginp,uu); % gi = find(~isnan(datTemp1)); datTemp1 = datTemp1(gi);
    %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %         erbarsBnp(:,uu) = ci_hi_temp1;
    %
    %         datTemp1 = diff_p(:,uu); % gi = find(~isnan(datTemp1)); datTemp1 = datTemp1(gi);
    %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %         erbarsDiffp(:,uu) = ci_hi_temp1;
    %
    %         datTemp1 = diff_np(:,uu); % gi = find(~isnan(datTemp1)); datTemp1 = datTemp1(gi);
    %         [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@mean,datTemp1},'alpha',0.05,'type','per');
    %         erbarsDiffnp(:,uu) = ci_hi_temp1;
    %
    %     end
    %
    %     %% LC Phasic.
    %
    %     axes(axs(3)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %
    %     plot(tax2A,nanmean(pethAllAp(gip,:)),'r-','linewidth',2);
    %     plot(tax2A,nanmean(pethAllBp(gip,:)),'k-','linewidth',2);
    %     legend('after beep','before beep','autoupdate','off','fontsize',8); legend('boxoff');
    %
    %     plot(tax2A,erbarsAp(1,:),'r--'); plot(tax2A,erbarsAp(2,:),'r--','linewidth',1);
    %     plot(tax2A,erbarsBp(1,:),'k--'); plot(tax2A,erbarsBp(2,:),'k--','linewidth',1);
    %
    %     xlabel('time (ms)'); ylabel('CCG');
    %     minY = 0; maxY = 0.27;
    %     axis([tax2A(1)-100 tax2A(end)+100 minY maxY]);
    %     plot([0 0],[0 0.3],'k--','linewidth',0.5);
    %
    %     %% LC No Phasic.
    %
    %     axes(axs(6)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %
    %     plot(tax2A,nanmean(pethAllAnp(gi,:)),'r-','linewidth',2);
    %     plot(tax2A,nanmean(pethAllBnp(gi,:)),'k-','linewidth',2);
    %
    %     plot(tax2A,erbarsAnp(1,:),'r--'); plot(tax2A,erbarsAnp(2,:),'r--','linewidth',1);
    %     plot(tax2A,erbarsBnp(1,:),'k--'); plot(tax2A,erbarsBnp(2,:),'k--','linewidth',1);
    %
    %     xlabel('time (ms)'); ylabel('CCG');
    %     minY = 0; maxY = 0.27;
    %     axis([tax2A(1)-100 tax2A(end)+100 minY maxY]);
    %     plot([0 0],[0 0.3],'k--','linewidth',0.5);
    %
    %     %% Difference plots.
    %
    %     %     axes(axs(4)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %     %
    %     %     plot(tax2A,nanmean(diff_p),'k-','linewidth',2);
    %     %     plot(tax2A,nanmean(diff_np),'k-','linewidth',0.5);
    %     %
    %     %     plot(tax2A,erbarsDiffp(1,:),'k--','linewidth',2); plot(tax2A,erbarsDiffp(2,:),'k--','linewidth',2);
    %     %     plot(tax2A,erbarsDiffnp(1,:),'k--','linewidth',0.5); plot(tax2A,erbarsDiffnp(2,:),'k--','linewidth',0.5);
    %     %
    %     %     xlabel('time (ms)'); ylabel('CCG');
    %     %     axis([tax2A(1)-200 tax2A(end)+200 -0.06 0.23]);
    %     %     plot([0 0],[-0.06 0.23],'k--','linewidth',0.5);
    %     %
    %     %%
    %
    %     %     for uu = 1:length(tax2)
    %     %         % p(uu) = signrank(smDat(:,uu));
    %     %         p(uu) = signrank(pethAll(:,uu));
    %     %         if p(uu) < 0.05
    %     %             plot(tax2(uu),0.8,'k*');
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








