% Fig5_LC_ACC_noBeep_pupil_phase.m
%
% INFO: Plot LC spiking and ACC rsc relative to pupil phases.
%
% Calls: "get_LC_ACC_pupil_phase.m" to do the heavy lifting.

%% Setup stuff.

% Clear everything?
clear; clear all;

monks = {'Sprout','Cicero'}; % Add mnks as needed.
monks = {'Sprout','Cicero','Cicero','Oz'}; % Add mnks as needed.
sites = {'LC_ACC_Fixation'}; % Only LC for now.
% base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.
base_dir = {'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data' ...
    'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data' ...
    'C:\Sidd\PostDoc2\Data\ALL_PUPIL_DATA\LC\Cicero\clean' ...
    'C:\Sidd\PostDoc2\Data\ALL_PUPIL_DATA\LC\Oz\clean'}

nMonks = length(monks); nSites = length(sites);

% Don't reanalyze data?
reanalyzeLC = false;
% reanalyzeACC = false;
% % Reanalyze data?
% reanalyzeLC = true;
% reanalyzeACC = true;

%% Analysis Section: Calculate ACC pairwise spike count correlations (rsc) and ACC spiking statistics conditioned on LC spiking.

% Loops through monkeys and sessions.
% Uses function "getLC_DualArea_rsc.m" to do the heavy lifting.
% Saves results for plotting later.

%% LC linked ACC measurements.

if reanalyzeLC
    % medianOut{jjm} = [median(loLC) median(hiLC) median(loACC) median(hiACC)];
    %     cd C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2018\;
    mnkRes = [];
    
    % mm=1; ss=1; ff=3; % For troubleshoooting - test any one file.
    for mm = 1:nMonks % Loop over monkeys (wheeeee....!).
        sitRes = [];
        
        for ss = 1:nSites  % Loop over brainstem sites.
            
            if mm<=2
                inDir= strcat([base_dir{mm},'\',monks{mm},'\',sites{ss},'\clean']); % Create dir name for input (clean) files.
            end
            if mm>2
                inDir= base_dir{mm};
            end
            
            cd(inDir);  dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
            fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
            rsc_dat = [];
            if ~isempty(fnames)
                nf = length(fnames);
                for ff = 1:nf
                    % for ff = 1:10 % For troubleshooting.
                    % ff = 26% For troubleshooting.
                    load(fnames{ff});
                    pause % For troubleshooting.
                    % rsc_dat{ff} = getLC_DualArea_rsc_060617(siteData); Used for analyses upto 060617.
                    rsc_dat{ff} = get_LC_ACC_pupil_phase(siteData,5,ff,mm); % Added analyses for window size effects.
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
    
    % savedir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2020\LC_ACC';
    savedir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC';
    cd(savedir);
    
    % save LC_pupil_phase_100120 mnkRes; % With new way of dividing trials into zero, mid, high.
    save LC_pupil_phase_051821 mnkRes; % With new way of dividing trials into zero, mid, high.
    
    clear mnkRes;
    % save(saveFileName,'mnkRes','-v7.3'); % Results structure is LARGE... need to use -v7.3
end % If reanalyze LC loop.

%% ACC linked LC measurements.

% if reanalyzeACC
%     % medianOut{jjm} = [median(loLC) median(hiLC) median(loACC) median(hiACC)];
%     %     cd C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2019\;
%     mnkRes = [];
%
%     % mm=1; ss=1; ff=3; % For troubleshoooting - test any one file.
%     for mm = 1:nMonks % Loop over monkeys (wheeeee....!).
%         sitRes = [];
%
%         for ss = 1:nSites  % Loop over brainstem sites.
%             inDir= strcat([base_dir,'\',monks{mm},'\',sites{ss},'\clean']); % Create dir name for input (clean) files.
%             cd(inDir);  dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
%             fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
%             rsc_dat = [];
%             if ~isempty(fnames)
%                 nf = length(fnames);
%                 for ff = 1:nf
%                     % for ff = 1:10 % For troubleshooting.
%                     % ff = 15
%                     load(fnames{ff});
%                     % pause % For troubleshooting.
%                     rsc_dat{ff} = get_LC_DualArea_rsc_pupil(siteData,5); % Added analyses for window size effects.
%                     disp(sprintf('File %d%sof%s%d', ff,' ',' ',nf));
%                 end
%             end
%             sitRes{ss} = rsc_dat;
%         end % Sites.
%         mnkRes{mm} = sitRes;
%     end % Mnks.
%
%     % Save analysis:
%     savedir = strcat([base_dir,'\Results\Results_2019\LC_ACC']);
%     cd(savedir);
%     % save LC_Results_rsc_pupil_032919b mnkRes; % With new way of dividing trials into zero, mid, high.
%     save LC_Results_rsc_pupil_040519 mnkRes; % With new way of dividing trials into zero, mid, high.
%
% end % If reanalyze ACC loop.

% ****************************************************

%% Load results file created above, organize data for plots.

clear; clear all;

% Don't organize results?
orgYesNo = false;
% Organize results?
orgYesNo = true;

if orgYesNo
    
    cd C:\Sidd\PostDoc2\SidCode\General;
    load('sidd_colors');
    c1 = cmap_M(110,:);
    c2 = cmap_M(170,:);
    
    
    %% Pupil phase linked LC results.
    
    %     backTime = -1000; fwdTime = 500; % We're looking at the 1sec to 2.1sec window of stable fixation.
    %     tmin   = backTime; tmax   = fwdTime; xBin = []; nb = []; tsize  = 200; tstep  = 100;
    %     xBin = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)']; nb = size(xBin,1);
    
    backTime = -1000; fwdTime = 500; % We're looking at the 1sec to 2.1sec window of stable fixation.
    tmin   = backTime; tmax   = fwdTime; xBin = []; nb = []; tsize  = 500; tstep  = 100;
    xBin = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)']; nb = size(xBin,1);
    
    pax = 0:15:345; nphase = length(pax);
    tax = nanmean(xBin,2);  ntime = length(tax);
    
    %     base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2019\LC_ACC'; % Base directory for brain area.
    %     inFile = 'LC_pupil_phase_043019'; % New way of dividing ACC trials into low, mid, high. [200 100]
    % base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2020\LC_ACC'; % Base directory for brain area.
    base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC'; % Base directory for brain area.
    % inFile = 'LC_pupil_phase_081520'; % New way of dividing ACC trials into low, mid, high. [200 100]
    
    %     inFile = 'LC_pupil_phase_100120'; % New way of dividing ACC trials into low, mid, high. [200 100]
    inFile = 'LC_pupil_phase_051821'; % New way of dividing ACC trials into low, mid, high. [200 100]
    
    savedir = base_dir;
    cd(savedir); load(inFile);
    
    % monks = {'Sprout', 'Cicero'}; % Add mnks as needed...
    monks = {'Sprout', 'Cicero','Cicero','Oz'}; % Add mnks as needed...
    sites = {'LC_ACC_Fixation'}; % For now LC only.
    nMonks = length(monks); nSites = length(sites);
    
    m_dat = [];
    for mm = 1:nMonks % Loop over mnks....
        u_dat = []; u_dat_S = []; u_dat_A = []; u_dat_AS = [];
        probLC = []; probACC = [];
        nlp = 0; nap = 0;
        all_index = 0;
        all_index2 = 0;
        nSites = length(sites);
        mnkDat = mnkRes{mm};
        sumLP = zeros(nphase,ntime);
        sumAP = zeros(nphase,ntime);
        rscDat_LC = zeros(nphase,ntime);
        rscDat_ACC = zeros(nphase,ntime);
        pwc_lc = [];
        pwc_acc = [];
        %         ccgDat_0 = [];
        %         ccgDat_180 = [];
        %         ccgDat_90 = [];
        %         ccgDat_270 = [];
        
        for ss = 1:nSites % Loop over brainstem sites.... LC and/or IC.
            sitDat = mnkDat{ss};
            
            if ~isempty(sitDat)
                for ap = 1:length(sitDat) % Loop over sessions.
                    dat = sitDat{ap}; % Get session data.
                    if ~isempty(dat)
                        
                        if ~isempty(dat{9})
                            temp_dat = dat{9}; % temp_dat(find(temp_dat==1))=nan;
                            pwc_lc = [pwc_lc temp_dat];
                        end
                        
                        if ~isempty(dat{10})
                            temp_dat = dat{10}; % temp_dat(find(temp_dat==1))=nan;
                            pwc_acc = [pwc_acc temp_dat];
                        end
                        
                        % pause
                        
                        %                         if ~isempty(dat{9})
                        %                             temp_dat = dat{9}; % temp_dat(find(temp_dat==1))=nan;
                        %                             ccgDat_0 = [ccgDat_0; temp_dat];
                        %                         end
                        %
                        %                         if ~isempty(dat{10})
                        %                             temp_dat = dat{10}; % temp_dat(find(temp_dat==1))=nan;
                        %                             ccgDat_90 = [ccgDat_90; temp_dat];
                        %                         end
                        %
                        %                         if ~isempty(dat{11})
                        %                             temp_dat = dat{11}; % temp_dat(find(temp_dat==1))=nan;
                        %                             ccgDat_180 = [ccgDat_180; temp_dat];
                        %                         end
                        %
                        %                         if ~isempty(dat{12})
                        %                             temp_dat = dat{12}; % temp_dat(find(temp_dat==1))=nan;
                        %                             ccgDat_270 = [ccgDat_270; temp_dat];
                        %                         end
                        
                        
                        if ~isempty(dat{7})
                            temp_dat = dat{7}; temp_dat(find(temp_dat==1))=nan;
                            rscDat_LC = cat(3,rscDat_LC,temp_dat);
                        end
                        if ~isempty(dat{8})
                            temp_dat = dat{8}; temp_dat(find(temp_dat==1))=nan;
                            rscDat_ACC = cat(3,rscDat_ACC,temp_dat);
                        end
                        
                        pspDat_LC = dat{5};
                        np = length(pspDat_LC);
                        for ppi = 1:np
                            sumLP = sumLP+pspDat_LC{ppi};
                        end
                        nlp = nlp+np;
                        
                        nu = length(dat{1});
                        for ui = 1:nu
                            if ~isempty(dat{1}{ui})
                                all_index = all_index + 1;
                                temp_spike_data = dat{1}{ui};
                                % temp_spike_data = zscore(temp_spike_data')';
                                u_dat(:,:,all_index) = temp_spike_data;
                                % u_dat_S(:,:,all_index) = dat{2}{ui};
                            end % Check for empty PETH matrix.
                        end % Units for this session.
                        
                        if mm < 3
                            
                            pspDat_ACC = dat{6};
                            np = length(pspDat_ACC);
                            for ppi = 1:np
                                sumAP = sumAP+pspDat_ACC{ppi};
                            end
                            nap = nap+np;
                            %
                            nu = length(dat{3});
                            for ui = 1:nu
                                if ~isempty(dat{3}{ui})
                                    all_index2 = all_index2 + 1;
                                    temp_spike_data = dat{3}{ui};
                                    % temp_spike_data = zscore(temp_spike_data')';
                                    u_dat_A(:,:,all_index2) = temp_spike_data;
                                    % u_dat_AS(:,:,all_index2) = dat{4}{ui};
                                end % Check for empty PETH matrix.
                            end % Units for this session.
                        end % Check for mm<3.
                        
                    end % Check for empty data structure.
                end % Sessions for this site.
            end % Check for empty site data structure.
        end % Sites for this mnky.
        % m_dat{mm} = {u_dat u_dat_S u_dat_A u_dat_AS};
        
        rscDat_LC(:,:,1) = [];
        rscDat_ACC(:,:,1) = [];
        
        probLC = sumLP/nlp;
        if mm<3
            probACC = sumAP/nap;
        end
        m_dat{mm} = {u_dat u_dat_S u_dat_A u_dat_AS probLC probACC rscDat_LC rscDat_ACC pwc_lc pwc_acc};
        % m_dat{mm} = {u_dat u_dat_S u_dat_A u_dat_AS probLC probACC rscDat_LC rscDat_ACC ccgDat_0 ccgDat_90 ccgDat_180 ccgDat_270};
    end % Monkeys.
    
end % Organize yes/no?


%% Plot.

% Don't plot results?
plotYesNo = false;
% Plot results?
plotYesNo = true;

% pause

pax = 0:15:345;
nphase = length(pax);
tax = nanmean(xBin,2);
ntime = length(tax);

x_ticks = [-800 -600 -400 -200 0 200 400];
y_ticks = [0 90 180 270 360];

% pause

if plotYesNo
    
    %% Setup figure:
    
    figureNumber = 6; num = 6; wid = 15; hts = [8]; cols = {1 1}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 3, 1, [12], 'Joshi and Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    
    %% Plot averaged LC PETH's.
    
    % sm_sigma = 0.1;
    
    % Pupil phase linked LC psth.
    all_dat_LC = (1000/tsize)*cat(3,m_dat{1}{1},m_dat{2}{1},m_dat{3}{1},m_dat{4}{1});
    nu = size(all_dat_LC,3);
    
    % Don't Z-score.
    all_dat_LC_Z = all_dat_LC;
    gi = [];
    
    % ************************************
    % Do Z-score.
    msiz = size(all_dat_LC(:,:,1));
    all_dat_LC_Z = nans(size(all_dat_LC));
    all_dat_LC_ZS = nans(size(all_dat_LC));
    %     figure;
    for ki = 1:nu
        %         if zi(ki) == 0
        this_u = all_dat_LC(:,:,ki);
        
        this_uZ = this_u(:);
        this_uZS = smooth(this_uZ,3);
        
        this_uZ = zscore(this_uZ);
        this_uZS = zscore(this_uZS);
        
        this_uZ = reshape(this_uZ,msiz);
        this_uZS = reshape(this_uZS,msiz);
        all_dat_LC_Z(:,:,ki) = this_uZ;
        all_dat_LC_ZS(:,:,ki) = this_uZS;
        
        %         % ************************************
        %         % Look at data. Done and re-checked: 082420, Sidd.
        %
        %         this_uZ_S = imgaussfilt(this_uZ,sm_sigma);
        %         axes(axs(1)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;
        %         imagesc(tax,pax,this_uZ_S); % hold on;
        %         colorbar('box','off');
        %         title('Spike PETH; LC','fontname','arial','fontsize',12);
        %         xlabel('time re:pupil phase (ms)','fontname','arial','fontsize',10);
        %         ylabel('pupil phase (degrees)','fontname','arial','fontsize',10);
        %         axis([-50+tax(1) 50+tax(end) -7.5+pax(1) 7.5+pax(end)]);
        %         xticks(x_ticks); yticks(y_ticks);
        %         plot([0 0],[pax(1) pax(end)],'k--');
        %
        %         pause;
        %         % clf
        %
        %         % ************************************
        %
        %         end
    end
    % Z-scoring is done.
    % ************************************
    
    sm_sigma = 0.2;
    mean_dat = nanmean(all_dat_LC_Z,3);
    mean_datS = nanmean(all_dat_LC_ZS,3);
    mean_datL = mean_dat;
    mean_datLS = imgaussfilt(mean_datS,sm_sigma);
    
    % pause
    max_out_L = [];
    
    % Look at data for timepoints preceding t=0 for each phase.
    taxN = tax(tax<0);
    for ii = 1:length(taxN)
        pdt = mean_datS(:,ii);
        max_pdt = max(pdt); max_i = find(pdt ==max_pdt); max_out_L(ii) = max_i;
    end
    
    %         axes(axs(1)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;
    %         %     colormap('parula');
    %         imagesc(tax,pax,mean_datLS); colorbar('box','off');
    %         title('Spike PETH; LC','fontname','arial','fontsize',12);
    %         xlabel('time re:pupil phase (ms)','fontname','arial','fontsize',10);
    %         ylabel('pupil phase (degrees)','fontname','arial','fontsize',10);
    %         axis([-50+tax(1) 50+tax(end) -7.5+pax(1) 7.5+pax(end)]);
    %         xticks(x_ticks); yticks(y_ticks);
    %         plot([0 0],[pax(1) pax(end)],'k--');
    
    % Indices for changing psth (ie if successive indices are the same, then it is not changing).
    gi = find(diff(max_out_L)~=0);
    gi = 1:gi(end);
    %     gi = 1:5;
    % gi = [gi gi(end)+1];
    % gi_out = gi;
    
    % plot(tax(gi),pax(max_out_L(gi)),'r.','markersize',20); box('off');
    %     pause
    % set(gca,'fontname','arial');
    
    % Fit a line to the peaks of the average data.
    yd = pax(max_out_L(gi));
    xd = tax(gi);
    
    %     b = polyfit(xd(:),yd(:),1);
    %     slopeL = b(1); intL = b(2); % statsL(kk) = stats(3);
    [b1,bint1,r,rint1,stats1]=regress(yd',[ones(length(xd),1) xd]);
    slopeL = b1(2); intL = b1(1); % statsL(kk) = stats(3);
    
    %     fitdat = intL+slopeL*xd;
    %     plot(xd,fitdat,'r-');
    %     xd_0_L = -intL./slopeL;
    
    % gi = 1:gi(end);
    max_psth = [];
    nu = size(all_dat_LC,3);
    sm_sigma = 0.5;
    for ki = 1:nu
        max_nu = [];
        % this_u = all_dat_LC(:,:,ki);
        this_u = all_dat_LC_ZS(:,:,ki);
        mean_datSm = this_u;
        % mean_datSm = imgaussfilt(this_u,sm_sigma);
        
        % for ji = 1:length(tax)
        for ji = 1:length(gi)
            % for ki = 1:nu
            % this_u = all_dat_LC(:,ji,ki);
            % this_u = smooth(this_u,3);
            
            % max_nu_temp = pax(find(mean_datSm(:,ji) == max(mean_datSm(:,ji))));
            max_nu_temp = pax(find(mean_datSm(:,gi(ji)) == max(mean_datSm(:,gi(ji)))));
            
            % max_nu = [max_nu max_nu_temp(end)];
            max_nu = [max_nu nanmedian(max_nu_temp)];
            % max_nu = [max_nu nanmean(max_nu_temp)];
            
        end
        max_psth = [max_psth; max_nu];
    end
    max_psth = max_psth';
    
    %     axes(axs(2)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;
    %     axi = 1:6;
    %     axes(axs(2)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;
    %     for ii = 1:length(axi)
    %         %        for kk = 1:length(max_psth)
    %             m1= bootci(1000,{@median,max_psth(ii,:)},'alpha',0.05,'type','per');
    %             plot(tax(ii),nanmedian(max_psth(ii,:)),'k.','markersize',6);
    %             plot([tax(ii) tax(ii)],[m1(1) m1(2)],'k-');
    %             %     end
    %     end
    
    % xd = [ones(length(tax(gi)),1) tax(gi)];
    % gi = 1:11;
    
    %     LOOK at each fit... what is the pattern, if any?
    %     figure; hold on;
    % xd = tax(gi);
    
    xd_0_L_All = [];
    slopesL = nans(1,nu);
    intsL = nans(1,nu);
    slopes_pvL = nans(1,nu);
    
    for kk = 1:nu
        yd = max_psth(:,kk);
        % [b,bint,r,rint,stats] = regress(yd',xd);
        %         if sum(yd==0)
        %             zii = zii+1;
        %             zi(kk) = kk;
        %         end
        %         if sum(yd>0)
        
        %         b = polyfit(xd,yd,1);
        %         slopesL(kk) = b(1); intL(kk) = b(2); % statsL(kk) = stats(3);
        %         xd_0_L_All = [xd_0_L_All (-b(2)/b(1))];
        
        [b1,bint1,r,rint1,stats1]=regress(yd,[ones(length(xd),1) xd]);
        slopesL(kk) = b1(2); intsL(kk) = b1(1); % statsL(kk) = stats(3);
        xd_0_L_All = [xd_0_L_All (-b1(1)/b1(2))];
        slopes_pvL(kk) = stats1(3);
        
        %         end
        %         fit_yd = b(1)*xd + b(2);
        %         plot(xd,yd,'ko'); hold on;
        %         plot(xd,fit_yd,'r-');
        %         pause; clf;
    end
    
    % gsi = 1:291;
    gsi = find(slopes_pvL<0.05);
    % gsi = find(xd_0_L_All>-1000 & xd_0_L_All<500);
    % gsi = find(isfinite(slopesL));
    % gsi = find(isfinite(xd_0_L_All));
    pL = signrank(slopesL(gsi));
    
    
    [ci_LC ci_LC_m]= bootci(2000,{@median,slopesL(gsi)},'alpha',0.05,'type','per');
    med_slope_LC = nanmedian(slopesL(gsi))
    slopeL_IQR = prctile(slopesL(gsi),[25 50 75]);
    
    % gsi = find(isfinite(xd_0_L_All));
    [ci_LC_Int ci_LC_Int_m]= bootci(2000,{@median,xd_0_L_All(gsi)},'alpha',0.05,'type','per');
    med_Xint_LC = nanmedian(xd_0_L_All(gsi));
    intXL_IQR = prctile(xd_0_L_All(gsi),[25 50 75]);
    
    [cii_LC cii_LC_m]= bootci(2000,{@median,intsL(gsi)},'alpha',0.05,'type','per');
    med_int_LC = nanmedian(intsL(gsi));
    intL_IQR = prctile(intsL(gsi),[25 50 75]);
    
    mean_datS = nanmean(all_dat_LC_ZS(:,:,gsi),3);
    mean_datS = imgaussfilt(mean_datS,0.5);
    
    axes(axs(1)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;
    %     colormap('parula');
    
    imagesc(tax,pax,mean_datS); colorbar('box','off');
    title('Spike PETH; LC','fontname','arial','fontsize',12);
    xlabel('time re:pupil phase (ms)','fontname','arial','fontsize',10);
    ylabel('pupil phase (degrees)','fontname','arial','fontsize',10);
    axis([-50+tax(1) 50+tax(end) -7.5+pax(1) 7.5+pax(end)]);
    xticks(x_ticks); yticks(y_ticks);
    plot([0 0],[pax(1) pax(end)],'k--');
    
    fitdat = med_int_LC+med_slope_LC*xd;
    plot(xd,fitdat,'r-','linewidth',2);
    % xd_0_L = -med_int_LC./med_slope_LC;
    
    %     n = hist(slopesL);
    %     figure; hold on
    
    %     % x-intercept.
    %     % y_all = slopesL*xd + intL;
    %     xd_0 = -intL./slopesL;
    %
    %     %     %     % Check if senshible...
    %     %     xax = tax(gi);
    %     %     figure; hold on;
    %     s1 = nanmedian(slopesL(intL<=-50));
    %     i1 = nanmedian(intL(intL<=-50));
    %     % s1 = nanmedian(slopesL);
    %     % i1 = nanmedian(intL);
    %     fitdat = i1+s1*xd;
    %     % plot(xax,fitdat,'r-');
    %     % plot([xax(1) xax(end)],[0 0],'k--');
    %     % fitdat = i1+s1*tax(gi);
    %
    %     plot(xd,fitdat,'r-');
    %
    
    %             pause
    
    
    %% Pupil linked ACC rsc data.
    
    ti = find(tax<0);
    taxN = tax(1:ti(end)+1);
    
    nu = size(all_rsc_ACC,3);
    
    gi=[];
    % Pupil phase linked ACC rsc.
    sm_sigma = .5;
    all_rsc_ACC = cat(3,m_dat{1}{8},m_dat{2}{8});
    %     all_rsc_ACC = cat(3,m_dat{1}{7},m_dat{2}{7});
    
    % Don't Z-score.
    all_rsc_ACC_Z = all_rsc_ACC;
    
    % ************************************
    % Do Z-score.
    msiz = size(all_rsc_ACC(:,:,1));
    all_rsc_ACC_Z = nans(size(all_rsc_ACC));
    all_rsc_ACC_ZS = nans(size(all_rsc_ACC));
    %     figure;
    for ki = 1:nu
        %         if zi(ki) == 0
        this_u = all_rsc_ACC(:,:,ki);
        
        this_uZ = this_u(:);
        this_uZS = smooth(this_uZ,2);
        
        % this_uZ = zscore(this_uZ);
        % this_uZS = zscore(this_uZS);
        
        this_uZ = reshape(this_uZ,msiz);
        this_uZS = reshape(this_uZS,msiz);
        
        all_rsc_ACC_Z(:,:,ki) = this_uZ;
        all_rsc_ACC_ZS(:,:,ki) = this_uZS;
        
        %         % ************************************
        %         % Look at data. Done and re-checked: 082420, Sidd.
        %
        %         this_uZ_S = imgaussfilt(this_uZ,sm_sigma);
        %         axes(axs(1)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;
        %         imagesc(tax,pax,this_uZ_S); % hold on;
        %         colorbar('box','off');
        %         title('Spike PETH; LC','fontname','arial','fontsize',12);
        %         xlabel('time re:pupil phase (ms)','fontname','arial','fontsize',10);
        %         ylabel('pupil phase (degrees)','fontname','arial','fontsize',10);
        %         axis([-50+tax(1) 50+tax(end) -7.5+pax(1) 7.5+pax(end)]);
        %         xticks(x_ticks); yticks(y_ticks);
        %         plot([0 0],[pax(1) pax(end)],'k--');
        %
        %         pause;
        %         % clf
        %
        %         % ************************************
        %
        %         end
    end
    % Z-scoring is done.
    % ************************************
    
    %     pause
    
    sm_sigma = .6;
    mdat_rscALL = nanmean(all_rsc_ACC_Z,3);
    mdat_rscALL_A = mdat_rscALL;
    mdat_rscALL_AS = imgaussfilt(mdat_rscALL,sm_sigma);
    
    min_out_A = [];
    
    for ii = 1:length(taxN)
        pdt = mdat_rscALL_A(:,ii);
        min_pdt = min(pdt); min_i = find(pdt == min_pdt); min_out_A(ii) = min_i;
    end
    
    %     axes(axs(2)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;
    %     imagesc(tax,pax,mdat_rscALL_AS);
    %     colorbar('box','off');
    %     title('rsc; ACC','fontname','arial','fontsize',12);
    %     xlabel('time re:pupil phase (ms)','fontname','arial','fontsize',10);
    %     ylabel('pupil phase (degrees)','fontname','arial','fontsize',10);
    %     axis([-50+tax(1) 50+tax(end) -7.5+pax(1) 7.5+pax(end)]);
    %     xticks(x_ticks); yticks(y_ticks);
    %     plot([0 0],[pax(1) pax(end)],'k--');
    
    giT = find(diff(min_out_A)~=0);
    gi = 3:giT(end);
    
    %     % plot(tax(gi),pax(min_out_A(gi)),'d','markersize',6,'markerfacecolor','r','markeredgecolor','none');
    %     box('off');
    %     set(gca,'fontname','arial');
    
    %     gi = intersect(gi,gi_out);
    
    
    
    pax_plot = pax(min_out_A(gi));
    
    %     % pax_plot(find(pax_plot<90)) = pax_plot(find(pax_plot<90))+180;
    %     plot(tax(gi),pax_plot,'d','markersize',6,'markerfacecolor','r','markeredgecolor','none');
    
    yd = pax_plot;
    % yd = pax(min_out_A(gi));
    xd = tax(gi);
    
    %     b = polyfit(xd(:),yd(:),1);
    %     slopeA = b(1); intA = b(2); % statsL(kk) = stats(3);
    
    [b1,bint1,r,rint1,stats1]=regress(yd',[ones(length(xd),1) xd]);
    slopeA = b1(2); intA = b1(1); % statsL(kk) = stats(3);
    
    %     fitdat = intA+slopeA*xd;
    %     plot(xd,fitdat,'r-');
    %     xd_0_A = -intA./slopeA;
    
    max_psth = [];
    %     nu = size(all_rsc_ACC,3);
    % sm_sigma = .1;
    for ki = 1:nu
        max_nu = [];
        this_u = all_rsc_ACC_ZS(:,:,ki);
        mean_datSm = this_u; % DON'T SMOOTH for fits!
        % mean_datSm = imgaussfilt(this_u,sm_sigma);
        for ji = 1:length(gi)
            if isfinite(sum(this_u))
                max_nu_temp = pax(find(mean_datSm(:,gi(ji)) == min(mean_datSm(:,gi(ji)))));
                % max_nu_temp = pax(find(this_u == min(this_u)));
                max_nu = [max_nu nanmedian(max_nu_temp)];
            end
            if ~isfinite(sum(this_u))
                max_nu = [max_nu nan];
            end
        end
        max_psth = [max_psth; max_nu];
    end
    max_psth = max_psth';
    
    nu = size(max_psth,2);
    % gi = gi_out;
    %     gi = [3 4 5 6];
    %     LOOK at each fit... what is the pattern, if any?
    %         figure; hold on;
    
    xd = tax(gi);
    xd_0_A_All = [];
    slopesA = nans(1,nu);
    intsA = nans(1,nu);
    slopes_pvA = nans(1,nu);
    
    for kk = 1:nu
        yd = max_psth(:,kk);
        % [b,bint,r,rint,stats] = regress(yd',xd);
        
        %         b = polyfit(xd,yd,1);
        %         slopesA(kk) = b(1); intL(kk) = b(2); % statsL(kk) = stats(3);
        %         xd_0_A_All = [xd_0_A_All (-b(2)/b(1))];
        
        [b1,bint1,r,rint1,stats1]=regress(yd,[ones(length(xd),1) xd]);
        slopesA(kk) = b1(2); intsA(kk) = b1(1); % statsL(kk) = stats(3);
        xd_0_A_All = [xd_0_A_All (-b1(1)/b1(2))];
        slopes_pvA(kk) = stats1(3);
        
        
        %             fit_yd = b(1)*xd + b(2);
        %             plot(xd,yd,'ko'); hold on;
        %             plot(xd,fit_yd,'r-');
        %             pause; clf;
        
    end
    
    gsiA = find(slopes_pvA<0.05);
    % gsi = find(xd_0_A_All>-1000 & xd_0_A_All<500);
    pA = signrank(slopesA(gsiA));
    
    % gsi = find(isfinite(slopesA));
    [ci_ACC ci_ACC_m]= bootci(2000,{@median,slopesA(gsiA)},'alpha',0.05,'type','per');
    med_slope_ACC = nanmedian(slopesA(gsiA))
    slopeA_IQR = prctile(slopesA(gsiA),[25 50 75]);
    
    % gsi = find(isfinite(xd_0_A_All) & slopes_pv<0.05);
    [ci_ACC_Int ci_ACC_Int_m]= bootci(2000,{@median,xd_0_A_All(gsiA)},'alpha',0.05,'type','per');
    med_Xint_ACC = nanmedian(xd_0_A_All(gsiA));
    intA_IQR = prctile(xd_0_A_All(gsiA),[25 50 75]);
    
    [cii_ACC cii_ACC_m]= bootci(2000,{@median,intsA(gsiA)},'alpha',0.05,'type','per');
    med_int_ACC = nanmedian(intsA(gsiA));
    intA_IQR = prctile(intsA(gsiA),[25 50 75]);
    
    mean_datS = nanmean(all_rsc_ACC_ZS(:,:,gsiA),3);
    mean_datS = imgaussfilt(mean_datS,0.5);
    
    axes(axs(2)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;
    imagesc(tax,pax,mean_datS);
    colorbar('box','off');
    title('rsc; ACC','fontname','arial','fontsize',12);
    xlabel('time re:pupil phase (ms)','fontname','arial','fontsize',10);
    ylabel('pupil phase (degrees)','fontname','arial','fontsize',10);
    axis([-50+tax(1) 50+tax(end) -7.5+pax(1) 7.5+pax(end)]);
    xticks(x_ticks); yticks(y_ticks);
    plot([0 0],[pax(1) pax(end)],'k--');
    giT = find(diff(min_out_A)~=0); gi = [gi giT(end)];
    % plot(tax(gi),pax(min_out_A(gi)),'d','markersize',6,'markerfacecolor','r','markeredgecolor','none');
    box('off');
    set(gca,'fontname','arial');
    
    fitdat =  med_int_ACC+med_slope_ACC*xd;
    plot(xd,fitdat,'r-','linewidth',2);
    % xd_0_A = -intA./slopeA;
    
    
    % x-intercept.
    
    %     % y_all = slopesL*xd + intL;
    %     xd_0_A = -intA./slopesA;
    %
    %     %     %     % Check if senshible...
    %     %     xax = tax(gi);
    %     %     figure; hold on;
    % %     s1 = nanmedian(slopesA(intA<=-50));
    % %     i1 = nanmedian(intA(intA<=-50));
    %     s1 = nanmean(slopesL);
    %     i1 = nanmean(intL);
    %     fitdat = i1+s1*xd;
    %     % plot(xax,fitdat,'r-');
    %     % plot([xax(1) xax(end)],[0 0],'k--');
    %     % fitdat = i1+s1*tax(gi);
    %     %     xd = tax(gi_out);
    %     plot(xd,fitdat,'r-');
    %
    %
    % %     xd = [ones(length(tax(gi)),1) tax(gi)];
    % %     for kk = 1:nu
    % %         yd = max_psth(kk,gi);
    % %         [b,bint,r,rint,stats] = regress(yd',xd);
    % %         slopesA(kk) = b(2); intA(kk) = b(1); statsA(kk) = stats(3);
    % %
    % %         % % Check if senshible...
    % %         % p = polyfit(xd,yd,1);
    % %         % slopesA(kk) = p(1); intA(kk) = p(2); % statsA(kk) = stats(3);
    % %         % fitdat = p(2)+p(1)*tax(1:7);
    % %         % plot(tax(1:7),max_psth(kk,1:7),'ko'); hold on;
    % %         % plot(tax(1:7),fitdat,'r-');
    % %         % pause; clf;
    % %
    % %     end
    
    
    %%
    
    %     axes(axs(3)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;
    %
    %     xbin = linspace(-1,1,15);
    %     nL = hist(slopesL,xbin);
    %     nA = hist(slopesA,xbin);
    %
    %     bar(xbin,nL/sum(nL),'linestyle','-','edgecolor','k','facecolor','k'); alpha(0.25);
    %     bar(xbin,nA/sum(nA),'linestyle','--','edgecolor','k','facecolor','k'); alpha(0.25);
    %
    %     prcL = quantile(slopesL(slopesL<0),[0.25 0.5 0.75]); % -0.33
    %     prcA = quantile(slopesA(slopesA<0),[0.25 0.5 0.75]); % -0.36
    %
    %     plot(prcL(2),0.20,'r.','markersize',20);
    %     plot([prcL(1) prcL(3)],[0.20 0.20],'k-');
    %
    %     plot(prcA(2),0.18,'rd','markersize',6,'markerfacecolor','r','markeredgecolor','none');
    %     plot([prcA(1) prcA(3)],[0.18 0.18],'k-');
    
    %     pL = signtest(slopesL);
    %     pA = signtest(slopesA);
    pAL = ranksum(slopesL(gsi),slopesA(gsiA))
    
    %     ylabel({'proportion of LC neurons';'or ACC pairs'},'fontname','arial','fontsize',10);
    %     xlabel('slope','fontname','arial','fontsize',10);
    %
    %     axis([-1.2 1.2 0 0.2]);
    %     set(gca,'fontname','arial');
    
    %     intA
    
end


%%

figureNumber = 61; num = 61; wid = 15; hts = [8]; cols = {2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 3, 1, [12], 'Joshi and Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
movegui(fig_,[1,1]); % Move figure so it is visible on screen.

%%

axes(axs(1)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;

all_dat_LC_ZS_Pre = all_dat_LC_ZS(:,1:11,:);
% all_dat_LC_ZS_Pre = all_dat_LC_ZS;

meanP = squeeze(nanmedian(all_dat_LC_ZS_Pre,2));
plot(pax, nanmedian(meanP,2),'k-','linewidth',2);
[ci_P ci_P_m] = bootci(2000,{@median,meanP'},'alpha',0.05,'type','norm'); % mean_ze = nanmean(ci_ze_m);

for kk = 1:length(pax)
    plot([pax(kk) pax(kk)],[ci_P(1,kk) ci_P(2,kk)],'k-');
end
% axis([-10 360 -0.18 0.25]);


axes(axs(2)); cla reset; ax = gca; disableDefaultInteractivity(ax); hold on;
all_rsc_ACC_ZS_Pre = all_rsc_ACC_ZS(:,1:11,:);
% all_rsc_ACC_ZS_Pre = all_rsc_ACC_ZS;

meanP = squeeze(nanmedian(all_rsc_ACC_ZS_Pre,2));
gi = find(isfinite(mean(meanP,1)))
plot(pax, nanmedian(meanP(:,gi),2),'k-','linewidth',2);
[ci_P ci_P_m] = bootci(2000,{@median,meanP(:,gi)'},'alpha',0.05,'type','norm'); % mean_ze = nanmean(ci_ze_m);

for kk = 1:length(pax)
    plot([pax(kk) pax(kk)],[ci_P(1,kk) ci_P(2,kk)],'k-');
end
% axis([-10 360 -0.18 0.25]);












