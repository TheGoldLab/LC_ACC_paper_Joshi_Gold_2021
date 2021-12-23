% Supp_4_3_ACC_LC_Dual_Recording_NoBeep_rsc
%
% INFO: Do stats, summarize and plot the rsc data that was organized in "Supp_4_3_ACC_LC_Dual_Recording_NoBeep_rsc.m".
%
% Added statistical testing for within mnky data.

clear; clear all;

cd ('C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2020\LC_ACC'); % Base directory for brain area.

%% Population average binned rsc.

xbs = linspace(200,1000,5); % 032220

nBins = length(xbs);
load('ACC_noBeep_rscData_051621');

all_mnk_dat = [];

for mm = 1:2
    
    scatter_allA1 = []; scatter_zeA1 = []; scatter_hiA1 = [];
    scatter_allA2 = []; scatter_zeA2 = []; scatter_hiA2 = [];
    scatter_allA3 = []; scatter_zeA3 = []; scatter_hiA3 = [];
    scatter_allA4 = []; scatter_zeA4 = []; scatter_hiA4 = [];
    all_datA_b = [];
    all_datZ_b = [];
    
    for kk = 1:nBins
        
        all_datZ_temp  = unit_all_lo{mm}(:,kk);
        all_datA_temp  = unit_all_all{mm}(:,kk);
        all_datNZ_temp  = unit_all_hi{mm}(:,kk);
        
        % All LC >0, for terciles.
        giA = find(isfinite(all_datA_temp) & abs(all_datA_temp)<1); all_datA_b{kk} = all_datA_temp(giA); giA = [];
        
        % All LC > 0, for LC zero terciles.
        giA = find(isfinite(all_datA_temp) & abs(all_datA_temp)<1 & isfinite(all_datZ_temp) & abs(all_datZ_temp)<1);
        all_datZ_b{kk} = all_datZ_temp(giA); all_datAZ_b{kk} = all_datA_temp(giA); giA = [];
        
        % All LC > 0, for LC zero/nonzero comparison.
        %         giA = find(isfinite(all_datA_temp) & abs(all_datA_temp)<1 & isfinite(all_datZ_temp) & abs(all_datZ_temp)<1 & isfinite(all_datNZ_temp) & abs(all_datNZ_temp)<1);
        %         all_datZ_b2{kk} = all_datZ_temp(giA); all_datAZ_b2{kk} = all_datA_temp(giA); all_datNZ_b2{kk} = all_datNZ_temp(giA); giA = [];
        all_datZ_b2{kk} = all_datZ_temp; all_datAZ_b2{kk} = all_datA_temp; all_datNZ_b2{kk} = all_datNZ_temp;
        
    end
    all_mnk_dat{mm} = {all_datA_b all_datZ_b all_datAZ_b all_datZ_b2 all_datAZ_b2 all_datNZ_b2}; % 13-18.
end

% pause

%% For paper text stats:

% all_100 = [all_mnk_dat{1}{13}{1}; all_mnk_dat{2}{13}{1}];
% all_1000 = [all_mnk_dat{1}{13}{10}; all_mnk_dat{2}{13}{10}];

all_200 = [all_mnk_dat{1}{5}{1}; all_mnk_dat{2}{5}{1}];
all_1000 = [all_mnk_dat{1}{5}{5}; all_mnk_dat{2}{5}{5}];

prc_200 = prctile(all_200,[25 50 75])
prc_1000 = prctile(all_1000,[25 50 75])


%% Plot rsc differences for each tercile.
plotAll = 0;
plotAll = 1;

if plotAll
    
    %% Organize data for plot.
    
    outmean_loNZD = []; outmean_midNZD = []; outmean_hiNZD = [];
    outse_lonz = []; outse_midnz = []; outse_hinz = [];
    outse_loZ = []; outse_midZ = []; outse_hiZ = [];
    outse_loNZD = []; outse_midNZD = []; outse_hiNZD = [];
    sigVal_lonz = []; sigVal_midnz = []; sigVal_hinz = [];
    sigVal_lonz1 = []; sigVal_midnz1 = []; sigVal_hinz1 = [];
    sigVal_lonz2 = []; sigVal_midnz2 = []; sigVal_hinz2 = [];
    
    sigValbnzD = [];
    sigValbnz = []; sigValb1 = []; sigValb2 = []; sigValb3 = []; sigValb4 = [];
    sigValDnz = []; sigValD1 = []; sigValD2 = []; sigValD3 = []; sigValD4 = [];
    outmean_lozb = []; outmean_midzb = []; outmean_hizb = [];
    
    outmean_lonz = []; outmean_midnz = []; outmean_hinz = [];
    outmean_lonz1 = []; outmean_midnz1 = []; outmean_hinz1 = [];
    outmean_lonz2 = []; outmean_midnz2 = []; outmean_hinz2 = [];
    outse_lo1D = []; outse_mid1D = []; outse_hi1D = [];
    outse_lo2D = []; outse_mid2D = []; outse_hi2D = [];
    outmean_loz = []; outmean_midz = []; outmean_hiz = [];
    
    sa_dat_nz = [];
    limVal = 99; % 95 until 042320.
    
    for jj = 1:nBins
        LC_NZ_loD = []; LC_NZ_midD = []; LC_NZ_hiD = [];
        LC_NZ_loD2 = []; LC_NZ_midD2 = []; LC_NZ_hiD2 = [];
        LC_Z_lo = []; LC_Z_mid = []; LC_Z_hi = []; LC_NZ_lo = []; LC_NZ_mid = []; LC_NZ_hi = [];
        
        for mm = 1:2
            %             aDat = all_mnk_dat{mm}{1}{jj}; gti = find(abs(aDat)>0); terc_vals = quantile(aDat(gti),2);
            %             aDat = all_mnk_dat{mm}{2}{jj}; aDatA = all_mnk_dat{mm}{3}{jj};
            %             t1 = find(aDatA<=terc_vals(1)); t2 = find(aDatA>terc_vals(1) & aDatA<terc_vals(2)); t3 = find(aDatA>=terc_vals(2));
            %             LC_Z_lo = [LC_Z_lo; aDat(t1)]; LC_Z_mid = [LC_Z_mid;  aDat(t2)]; LC_Z_hi = [LC_Z_hi;  aDat(t3)];
            
            % ********************
            % LC NZ.
            % Central measures for data.
            D1 = []; D2 = []; D3 = [];
            aDat = all_mnk_dat{mm}{5}{jj}; mdat = all_mnk_dat{mm}{6}{jj}; zdat = all_mnk_dat{mm}{4}{jj};
            terc_vals = quantile(aDat,2);
            t1 = find(aDat<=terc_vals(1)); t2 = find(aDat>terc_vals(1) & aDat<terc_vals(2)); t3 = find(aDat>=terc_vals(2));
            
            LC_Z_lo = [LC_Z_lo; zdat(t1)]; LC_Z_mid = [LC_Z_mid; zdat(t2)]; LC_Z_hi = [LC_Z_hi; zdat(t3)];
            LC_NZ_lo = [LC_NZ_lo; mdat(t1)]; LC_NZ_mid = [LC_NZ_mid; mdat(t2)]; LC_NZ_hi = [LC_NZ_hi; mdat(t3)];
            
            A1=aDat(t1); A2 = aDat(t2); A3 = aDat(t3);
            D1=mdat(t1)-zdat(t1); D2 = mdat(t2)-zdat(t2); D3 = mdat(t3)-zdat(t3);
            % D1 = D1(abs(D1)<prctile(abs(D1),limVal)); D2 = D2(abs(D2)<prctile(abs(D2),limVal)); D3 = D3(abs(D3)<prctile(D3,limVal));
            LC_NZ_loD = [LC_NZ_loD; D1]; LC_NZ_midD = [LC_NZ_midD; D2]; LC_NZ_hiD = [LC_NZ_hiD; D3]; % pause
            p1 = ranksum(mdat(t1),zdat(t1)); p2= ranksum(mdat(t2),zdat(t2)); p3 = ranksum(mdat(t3),zdat(t3));
            sigVal_lonz(mm) = p1; sigVal_midnz(mm) = p2; sigVal_hinz(mm) = p3;
            sa_dat_nz_m{mm} = [mdat(t3) zdat(t3)];
            p1 = signrank(D1(isfinite(D1))); p2= signrank(D2(isfinite(D2))); p3 = signrank(D3(isfinite(D3)));
            sigVal_loZD(mm) = p1; sigVal_midZD(mm) =  p2; sigVal_hiZD(mm) = p3;
            datD_lo{mm} = D1;
            datD_mid{mm} = D2;
            datD_hi{mm} = D3;
            datD_all{mm} = [D1; D2; D3];
            datA_all{mm} = [A1; A2; A3];
            
        end % Mnks for this bin.
        
        binDiffValsL{jj} = {datD_lo datD_mid datD_hi};

        
        % Next: Bootstrapped estimates for confidence intervals (for error bars).
        % LC zero.
        
        LC_Z_lo = LC_Z_lo(isfinite(LC_Z_lo)); LC_Z_mid = LC_Z_mid(isfinite(LC_Z_mid)); LC_Z_hi = LC_Z_hi(isfinite(LC_Z_hi));
        outmean_loz = [outmean_loz nanmedian(LC_Z_lo)]; outmean_midz = [outmean_midz nanmedian(LC_Z_mid)]; outmean_hiz = [outmean_hiz nanmedian(LC_Z_hi)];
        [ci_lo1 ci_lo1_m]= bootci(2000,{@median,LC_Z_lo},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_mid1 ci_mid1_m]= bootci(2000,{@median,LC_Z_mid},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_hi1 ci_hi1_m]= bootci(2000,{@median,LC_Z_hi},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        outse_loZ = [outse_loZ ci_lo1]; outse_midZ = [outse_midZ ci_mid1]; outse_hiZ = [outse_hiZ ci_hi1];
        
        % ********************
        % LC NZ.
        % Central measures for data.
        LC_NZ_lo = LC_NZ_lo(isfinite(LC_NZ_lo)); LC_NZ_mid = LC_NZ_mid(isfinite(LC_NZ_mid)); LC_NZ_hi = LC_NZ_hi(isfinite(LC_NZ_hi));
        outmean_lonz = [outmean_lonz nanmedian(LC_NZ_lo)]; outmean_midnz = [outmean_midnz nanmedian(LC_NZ_mid)]; outmean_hinz = [outmean_hinz nanmedian(LC_NZ_hi)];
        [ci_lo1 ci_lo1_m]= bootci(2000,{@median,LC_NZ_lo},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_mid1 ci_mid1_m]= bootci(2000,{@median,LC_NZ_mid},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_hi1 ci_hi1_m]= bootci(2000,{@median,LC_NZ_hi},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        outse_lonz = [outse_lonz ci_lo1]; outse_midnz = [outse_midnz ci_mid1]; outse_hinz = [outse_hinz ci_hi1];
        
        % Central measures for difference data.
        LC_NZ_loD = LC_NZ_loD(isfinite(LC_NZ_loD)); LC_NZ_midD = LC_NZ_midD(isfinite(LC_NZ_midD)); LC_NZ_hiD = LC_NZ_hiD(isfinite(LC_NZ_hiD));
        gi1=find(isfinite((LC_NZ_loD))); gi2=find(isfinite((LC_NZ_midD))); gi3=find(isfinite((LC_NZ_hiD)));
        outmean_loNZD = [outmean_loNZD nanmedian(LC_NZ_loD(gi1))];
        outmean_midNZD = [outmean_midNZD nanmedian(LC_NZ_midD(gi2))];
        outmean_hiNZD = [outmean_hiNZD nanmedian(LC_NZ_hiD(gi3))];
        [ci_lo1D ci_lo1_m]= bootci(2000,{@median,LC_NZ_loD(gi1)},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_mid1D ci_mid1_m]= bootci(2000,{@median,LC_NZ_midD(gi2)},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_hi1D ci_hi1_m]= bootci(2000,{@median,LC_NZ_hiD(gi3)},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        outse_loNZD = [outse_loNZD ci_lo1D]; outse_midNZD = [outse_midNZD ci_mid1D]; outse_hiNZD = [outse_hiNZD ci_hi1D];
        
        dL1 = datD_lo{1};    dL1 = dL1(isfinite(dL1));  outmean_lonz1 = [outmean_lonz1 nanmedian(dL1)];
        dM1 = datD_mid{1}; dM1 = dM1(isfinite(dM1)); outmean_midnz1 = [outmean_midnz1 nanmedian(dM1)];
        dH1 = datD_hi{1};   dH1 = dH1(isfinite(dH1));  outmean_hinz1 = [outmean_hinz1 nanmedian(dH1)];
        [ci_lo1D ci_lo1_m]      = bootci(2000,{@median,dL1},'alpha',0.05,'type','per'); outse_lo1D = [outse_lo1D ci_lo1D];
        [ci_mid1D ci_mid1_m] = bootci(2000,{@median,dM1},'alpha',0.05,'type','per'); outse_mid1D = [outse_mid1D ci_mid1D];
        [ci_hi1D ci_hi1_m]      = bootci(2000,{@median,dH1},'alpha',0.05,'type','per'); outse_hi1D = [outse_hi1D ci_hi1D];
        
        dL2 = datD_lo{2};    dL2 = dL2(isfinite(dL2));  outmean_lonz2 = [outmean_lonz2 nanmedian(dL2)];
        dM2 = datD_mid{2}; dM2 = dM2(isfinite(dM2)); outmean_midnz2 = [outmean_midnz2 nanmedian(dM2)];
        dH2 = datD_hi{2};   dH2 = dH2(isfinite(dH2));  outmean_hinz2 = [outmean_hinz2 nanmedian(dH2)];
        [ci_lo2D ci_lo1_m]      = bootci(2000,{@median,dL2},'alpha',0.05,'type','per'); outse_lo2D = [outse_lo2D ci_lo2D];
        [ci_mid2D ci_mid1_m] = bootci(2000,{@median,dM2},'alpha',0.05,'type','per'); outse_mid2D = [outse_mid2D ci_mid2D];
        [ci_hi2D ci_hi1_m]      = bootci(2000,{@median,dH2},'alpha',0.05,'type','per'); outse_hi2D = [outse_hi2D ci_hi2D];
        
        % Tests
        pv_loD(jj) = signrank(LC_NZ_loD);
        pv_midD(jj) = signrank(LC_NZ_midD);
        pv_hiD(jj) = signrank(LC_NZ_hiD);
        
        % **********************************
        % Collect all sig vals.
        sigValbnz = [sigValbnz [sigVal_lonz'; sigVal_midnz'; sigVal_hinz']];
        sigValDnz = [sigValDnz [sigVal_loZD'; sigVal_midZD'; sigVal_hiZD']];
        
        datD_A{jj} = datD_all;
        datA_A{jj} = datA_all;
        
    end % Bins.
    
    %     1234
    
    %     pause
    
    %% SUPPLEMENT for Figure 4: Raw data shown as scatter.
    
    %     % SUPPLEMENT Fig 4: Scatter plots with raw data for largest binsize, bins 1-5.
    %
    %     figureNumber = 41; num = 41; wid = 15; hts = [4]; cols = {2 2 2 2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,1,1, [12], 'Joshi & Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    %     movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    %
    %     axlim = 0.7; msiz = 8;
    %     bsiz=1:nBins;
    %     for jj = 1:5
    %
    %         axes(axs(2*jj -1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %         plot(sa_dat_nz{jj}{1}(:,2),sa_dat_nz{jj}{1}(:,1),'k.'); plot([-1 1],[-1 1],'g--','linewidth',.5); plot(nanmedian(sa_dat_nz{jj}{1}(:,2)),nanmedian(sa_dat_nz{jj}{1}(:,1)),'r.','markersize',msiz);
    %         % p = ranksum(sa_dat_nz{jj}{1}(:,2),sa_dat_nz{jj}{1}(:,1));
    %         p = sigValDnz(5,jj);
    %         text(-0.18,-0.18,strcat('p=',num2str(p),', binsize=',num2str(xbs(jj))),'fontsize',8);
    %         axis([-axlim axlim -axlim axlim]); axis square;
    %
    %         axes(axs(2*jj)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %         plot(sa_dat_nz{jj}{2}(:,2),sa_dat_nz{jj}{2}(:,1),'k.'); plot([-1 1],[-1 1],'g--','linewidth',.5); plot(nanmedian(sa_dat_nz{jj}{2}(:,2)),nanmedian(sa_dat_nz{jj}{2}(:,1)),'r.','markersize',msiz);
    %         % p = ranksum(sa_dat_nz{jj}{2}(:,2),sa_dat_nz{jj}{2}(:,1));
    %         p = sigValDnz(6,jj);
    %         text(-0.18,-0.18,strcat('p=',num2str(p),', binsize=',num2str(xbs(jj))),'fontsize',8);
    %         axis([-axlim axlim -axlim axlim]); axis square;
    %
    %     end
    %
    %     % SUPPLEMENT Fig 4: Scatter plots with raw data for largest binsize, bins 6-10.
    %
    %     figureNumber = 42; num = 42; wid = 17.6; hts = [4]; cols = {2 2 2 2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,1,1, [12], 'Joshi & Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    %     movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    %
    %     for jj = 6:nBins
    %
    %         kk = jj-5;
    %
    %         axes(axs(2*kk -1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %         plot(sa_dat_nz{jj}{1}(:,2),sa_dat_nz{jj}{1}(:,1),'k.'); plot([-1 1],[-1 1],'g--','linewidth',.5); plot(nanmedian(sa_dat_nz{jj}{1}(:,2)),nanmedian(sa_dat_nz{jj}{1}(:,1)),'r.','markersize',msiz);
    %         % p = ranksum(sa_dat_nz{jj}{1}(:,2),sa_dat_nz{jj}{1}(:,1));
    %         p = sigValDnz(5,jj);
    %         text(-0.18,-0.18,strcat('p=',num2str(p),', binsize=',num2str(xbs(jj))),'fontsize',8,'fontname','arial');
    %         axis([-axlim axlim -axlim axlim]); axis square;
    %         p = sigValDnz(5,jj);
    %
    %         axes(axs(2*kk)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %         plot(sa_dat_nz{jj}{2}(:,2),sa_dat_nz{jj}{2}(:,1),'k.'); plot([-1 1],[-1 1],'g--','linewidth',.5); plot(nanmedian(sa_dat_nz{jj}{2}(:,2)),nanmedian(sa_dat_nz{jj}{2}(:,1)),'r.','markersize',msiz);
    %         % p = ranksum(sa_dat_nz{jj}{2}(:,2),sa_dat_nz{jj}{2}(:,1));
    %         p = sigValDnz(6,jj);
    %         text(-0.18,-0.18,strcat('p=',num2str(p),', binsize=',num2str(xbs(jj))),'fontsize',8,'fontname','arial');
    %         axis([-axlim axlim -axlim axlim]); axis square;
    %
    %     end
    
    %         pause
    
    %% Supp Figure 4-3: ACC linked LC rsc.
    
    % Setup figure.
    
    figureNumber = 43; num = 43; wid = 17.6; hts = [6]; cols = {2 2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,2,1, [12], 'Joshi and Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    %% Now plot.
    
    xbs = linspace(200,1000,5); % 032220
    
    outmean_all=[]; sigVal_all1 = []; sigVal_all2 = []; outse_all = []; outmean_allD=[];
    
    outmean_all{1} = [outmean_loz; outmean_lonz];
    outse_all{1} = [outse_loZ; outse_lonz];
    outmean_all{2} = [outmean_midz; outmean_midnz];
    outse_all{2} = [outse_midZ; outse_midnz];
    outmean_all{3} = [outmean_hiz; outmean_hinz];
    outse_all{3} = [outse_hiZ; outse_hinz];
    
    % Create binsizes spaced equally in log10 space:
    %     xbs = round(logspace(log10(100),log10(1000),10),1);
    ms = [8 5 20 9 11];
    ms = 10;
    colors = {0.5.*ones(1,3), zeros(1,3)};
    xVals = 100*[-0.15 0.15];
    yMax = 0.35; yMin = -0.13;
    
    for ti = 1:3 %Terciles.
        axes(axs(2*ti -1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
        for uj = 1:nBins % Bins.
            for ui = 1:2 % LC spike=0,1,2,3,>3.
                
                %                 plot(xbs(uj)-15,mean_lo(uj),'ko','markersize',5,'linewidth',1);
                %                 plot(xbs(uj)+15,mean_hi(uj),'kd','markersize',5,'markerfacecolor','k');
                
                if ui==1
                    plot(xbs(uj)+xVals(ui),outmean_all{ti}(ui,uj),'ko','markersize',5,'linewidth',1);
                end
                if ui==2
                    plot(xbs(uj)+xVals(ui),outmean_all{ti}(ui,uj),'kd','markersize',5,'markerfacecolor','k');
                end
                
                plot([xbs(uj)+xVals(ui) xbs(uj)+xVals(ui)],[outse_all{ti}(2*ui - 1,uj) outse_all{ti}(2*ui,uj)],'k-');
                
                %                 plot(xbs(uj)+xVals(ui),outmean_all{ti}(ui,uj),'.','color',colors{ui},'markersize',ms(3));
                %                 plot([xbs(uj)+xVals(ui) xbs(uj)+xVals(ui)],[outse_all{ti}(2*ui - 1,uj) outse_all{ti}(2*ui,uj)],'k-');
                
                
                %                 if ui>1
                %                     if (sigVal_all1{ti}(ui-1,uj)<0.05 & sigVal_all2{ti}(ui-1,uj)<0.05)
                %                         plot(xbs(uj)+xVals(ui),0.21,'k*');
                %                     end
                %                 end
            end
        end
        title(['tercile',num2str(ti)]);
        ylabel('ACC rsc','fontname','arial','fontsize',10);
        xlabel('binsize (msec)','fontname','arial','fontsize',10);
        axis([100 1100 yMin yMax]); plot([100 1100], [0 0],'k--');% ylabel('ACC rsc diff');
    end
    
    set(gca,'fontname','arial');
    
    
    % save forShufPlot outmean_all outse_all;
    
    % **********************************************************************************
    %% PLOT BARS.
    
    %% PLOT BARS.
    colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    
    outmean_allD{1} = outmean_loNZD; outmean_allD{2} = outmean_midNZD; outmean_allD{3} = outmean_hiNZD;
    outmean_allD1{1} = outmean_lonz1; outmean_allD1{2} = outmean_midnz1; outmean_allD1{3} = outmean_hinz1;
    outmean_allD2{1} = outmean_lonz2; outmean_allD2{2} = outmean_midnz2; outmean_allD2{3} = outmean_hinz2;
    
    outse_allD{1} = outse_loNZD;
    outse_allD{2} = outse_midNZD;
    outse_allD{3} = outse_hiNZD;
    
    outse_allD1{1} = outse_lo1D;
    outse_allD1{2} = outse_mid1D;
    outse_allD1{3} = outse_hi1D;
    
    outse_allD2{1} = outse_lo2D;
    outse_allD2{2} = outse_mid2D;
    outse_allD2{3} = outse_hi2D;
    
    sigValDA1{1} = sigValDnz(1,:);
    sigValDA2{1} = sigValDnz(2,:);
    sigValDA1{2} = sigValDnz(3,:);
    sigValDA2{2} = sigValDnz(4,:);
    sigValDA1{3} = sigValDnz(5,:);
    sigValDA2{3} = sigValDnz(6,:);
    
    pv_bothMnks{1}=pv_loD;
    pv_bothMnks{2}=pv_midD;
    pv_bothMnks{3}=pv_hiD;
    
    xjit = 30;
    yMax = 55; yMin = -55;
    for ti = 1:3 % Terciles.
        axes(axs(2*ti)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
        for uj = 1:nBins % Bins.
            ui=1;
            
            % Plot bars.
            if (sigValDA1{ti}(1,uj)<0.05)
                bar(xbs(uj)-xjit,outmean_allD1{ti}(ui,uj),2*(xjit-2),'facecolor',colors{2},'edgecolor',colors{2},'linewidth',1);
            end
            if (sigValDA1{ti}(1,uj)>=0.05)
                bar(xbs(uj)-xjit,outmean_allD1{ti}(ui,uj),2*(xjit-2),'edgecolor',colors{2},'facecolor','none','linewidth',1);
            end
            if (sigValDA2{ti}(1,uj)<0.05)
                bar(xbs(uj)+xjit,outmean_allD2{ti}(ui,uj),2*(xjit-2),'facecolor',colors{3},'edgecolor',colors{3},'linewidth',1);
            end
            if (sigValDA2{ti}(1,uj)>=0.05)
                bar(xbs(uj)+xjit,outmean_allD2{ti}(ui,uj),2*(xjit-2),'edgecolor',colors{3},'facecolor','none','linewidth',1);
            end
            
            % Plot CI's.
            plot([xbs(uj)-xjit xbs(uj)-xjit],[outse_allD1{ti}(1,uj) outse_allD1{ti}(2,uj)],'k-');
            plot([xbs(uj)+xjit xbs(uj)+xjit],[outse_allD2{ti}(1,uj) outse_allD2{ti}(2,uj)],'k-');
            
            %             if (sigValDA1{ti}(1,uj)<0.05 & sigValDA2{ti}(1,uj)<0.05)
            %                 plot(xbs(uj),0.04,'k+','markersize',6);
            %             end
            if (pv_bothMnks{ti}(1,uj)<0.05)
                plot(xbs(uj),0.17,'*','color',colors{4},'markersize',6);
            end
            
        end % Bins loop.
        
        axis([100 1100 -0.14 0.17]);
        ylabel('\DeltaLC rsc re:LCzero','fontname','arial','fontsize',10);
        xlabel('binsize (msec)','fontname','arial','fontsize',10);
    end
    
    %     %%
    %     colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    %
    %     xjit = 30;
    %
    %     outmean_allD{1} = outmean_loNZD*100;
    %     outse_allD{1} = outse_loNZD;
    %     outmean_allD{2} = outmean_midNZD*100;
    %     outse_allD{2} = outse_midNZD;
    %     outmean_allD{3} = outmean_hiNZD*100;
    %     outse_allD{3} = outse_hiNZD;
    %
    %     %     sigVal_all1{1} = [sigValD1(1,:);sigValD2(1,:);sigValD3(1,:);sigValD4(1,:)];
    %     %     sigVal_all2{1} = [sigValD1(2,:);sigValD2(2,:);sigValD3(2,:);sigValD4(2,:)];
    %     %     sigVal_all1{2} = [sigValD1(3,:);sigValD2(3,:);sigValD3(3,:);sigValD4(3,:)];
    %     %     sigVal_all2{2} = [sigValD1(4,:);sigValD2(4,:);sigValD3(4,:);sigValD4(4,:)];
    %     %     sigVal_all1{3} = [sigValD1(5,:);sigValD2(5,:);sigValD3(5,:);sigValD4(5,:)];
    %     %     sigVal_all2{3} = [sigValD1(6,:);sigValD2(6,:);sigValD3(6,:);sigValD4(6,:)];
    %
    %     sigValDA1{1} = sigValDnz(1,:);
    %     sigValDA2{1} = sigValDnz(2,:);
    %     sigValDA1{2} = sigValDnz(3,:);
    %     sigValDA2{2} = sigValDnz(4,:);
    %     sigValDA1{3} = sigValDnz(5,:);
    %     sigValDA2{3} = sigValDnz(6,:);
    %
    %     sigValD{1} = pv_loD;
    %     sigValD{2} = pv_midD;
    %     sigValD{3} = pv_hiD;
    %
    %
    %     %     ylimVals = [-150 100;-110 100; -65 60];
    %     yMax = 100; yMin = -100;
    %     for ti = 1:3 % Terciles.
    %         axes(axs(2*ti)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %         for uj = 1:nBins % Bins.
    %             % for ui = 1:4 % LC spike=1,2,3,>3 relative to LC = 0.
    %             ui=1;
    %             bar(xbs(uj)+xjit,outmean_allD{ti}(ui,uj),15,'facecolor',colors{ui+1},'edgecolor','none');
    %             plot([xbs(uj)+xjit xbs(uj)+xVals(ui+1)],100*[outse_allD{ti}(2*ui-1,uj) outse_allD{ti}(2*ui,uj)],'k-');
    %
    %             % Plot sig symbols if sig for BOTH mnks.
    %             if (sigValDA1{ti}(ui,uj)<0.05 & sigValDA2{ti}(ui,uj)<0.05)
    %                 plot(xbs(uj)+xVals(ui+1),-75,'k+');
    %             end
    %
    %             % end % LC spike conditions loop.
    %
    %             % Plot sig symbols for LC zero vs nonzero if sig for BOTH mnks.
    %             % if (sigValDA1{ti}(1,uj)<0.05 & sigValDA2{ti}(1,uj)<0.05)
    %             if (sigValD{ti}(1,uj)<0.05)
    %                 plot(xbs(uj),-50,'k*');
    %             end
    %
    %         end % Bins loop.
    %         % axis([50 1050 yMin yMax]); plot([50 1050], [0 0],'k--');% ylabel('ACC rsc diff');
    %         %         axis([50 1050+0.5 ylimVals(ti,1) ylimVals(ti,2)]); plot([50 1050], [0 0],'k--');% ylabel('ACC rsc diff');
    %
    %         axis([250 1050+0.5 yMin yMax]); plot([50 1050], [0 0],'k--');% ylabel('ACC rsc diff');
    %         ylabel('\DeltaACC rsc re:LCzero','fontname','arial','fontsize',10);
    %         xlabel('binsize (msec)','fontname','arial','fontsize',10);
    %     end
    %
    %     set(gca,'fontname','arial');
    
%             pause
    
    %%
    
    all_bdat_LC = cell(1,length(xbs));
    
    for jj = 1:length(xbs)
        
        dat_all = [datA_A{jj}{1}; datA_A{jj}{2}];
        diff_A1 = [datD_A{jj}{1}; datD_A{jj}{2}];
        
        gi = find(isfinite(dat_all) & isfinite(diff_A1));
        diff_A1 = diff_A1(gi); dat_all = dat_all(gi);
        
        terc_vals = quantile(dat_all,2);
        t1 = find(dat_all<=terc_vals(1)); t2 = find(dat_all>terc_vals(1) & dat_all<terc_vals(2)); t3 = find(dat_all>=terc_vals(2));
        
        %
        
        %         bin_yd{1} = diff_A1(t1);
        %         bin_yd2{2} = diff_A1(t2);
        %         bin_yd3{3}= diff_A1(t3);
        
        all_bdat_LC{jj} = {diff_A1(t1) diff_A1(t2) diff_A1(t3)};
    end
    
    cd C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC;
    save Fig4_ACC_LC_Compare all_bdat_LC binDiffValsL;
    
    %% SUPPLEMENT Fig 4: Relationship between unconditioned ACC rsc and ACC rsc change between LC zero and nonzero conditions.
    
    figureNumber = 44; num = 44; wid = 17.6; hts = [5]; cols = {5}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,1,1, [12], 'Joshi and Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    %% Plot.
    hi_slope = zeros(5,10);
    lo_slope = zeros(5,10);
    hi_slope_stat = zeros(5,10);
    lo_slope_stat = zeros(5,10);
    pout = [];
    
    yM = .5;
    xM = .5;
    % ax_ind = [0 5 10 15 20];
    
    xb = linspace(-0.5,0.5,15);
    yb = linspace(-0.5,0.5,15);
    
    %     xbs = linspace(100,1000,10); % 032220
    
    ax_ind = 0;
    
    % for jj = 1:10
    for jj = 1:nBins
        %     for jj = 6:10
        
        ax_ind = ax_ind + 1;
        
        axes(axs(ax_ind)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
        
        %         dat_all = [all_mnk_dat{1}{5}{jj};all_mnk_dat{2}{5}{jj}];
        %         dat_Z = [all_mnk_dat{1}{4}{jj};all_mnk_dat{2}{4}{jj}];
        %         dat_hi = [all_mnk_dat{1}{6}{jj};all_mnk_dat{2}{6}{jj}];
        %         diff_A1 = dat_hi-dat_Z;
        %
        %         gi = find(isfinite(diff_A1) & isfinite(dat_all));
        %         %         gi = find(abs(dat_all)<prctile(abs(dat_all),99) & abs(diff_A1)<prctile(abs(diff_A1),99));
        %         diff_A1 = diff_A1(gi);
        %         dat_all = dat_all(gi);
        
        dat_all = [datA_A{jj}{1}; datA_A{jj}{2}];
        diff_A1 = [datD_A{jj}{1}; datD_A{jj}{2}];
        %         size(diff_A1)
        gi = find(isfinite(dat_all) & isfinite(diff_A1));
        diff_A1 = diff_A1(gi);
        dat_all = dat_all(gi);
        
        plot(dat_all,diff_A1,'.','markersize',4);
        % histogram2(dat_all,diff_A1,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');
        % colormap(gray);
        
        [b i] = sort(dat_all); mm1 = movmean(diff_A1(i),100);
        plot(dat_all(i),mm1,'-','linewidth',2,'color',[0.5 0.5 0.5]);
        % xd = [ones(size(dat_all)) dat_all];
        yd = diff_A1;
        
        %         [b,bint,r,rint,stats]=regress(yd,xd);
        
        xd = dat_all;
        mdl1 = fitlm(xd,yd);
        coeffs = mdl1.Coefficients.Estimate;
        pvals = mdl1.Coefficients.pValue;
        slope = coeffs(2); interc = coeffs(1);
        b = [interc slope];
        stats(3) = pvals(2);
        slope
        interc
        pout = [pout pvals(2)];
        
        x_plot = linspace(min(dat_all),max(dat_all),10);
        if stats(3) < 0.05
            plot(x_plot,b(1)+b(2)*x_plot,'k-','linewidth',1);
        else
            plot(x_plot,b(1)+b(2)*x_plot,'k--','linewidth',1);
        end
        hi_slope(kk,jj) = b(2); hi_slope_stat(kk,jj) = stats(3);
        plot([-0.5 0.5],[0 0],'k--'); plot([0 0],[-0.5 0.5],'k--');
        % text(-xM,yM,strcat([num2str(length(dat_all)),', ',num2str(b(2)),'bsiz=',num2str(round(xbs(jj)))]),'fontsize',7,'color','k','fontname','arial');
        text(-xM,yM,strcat([num2str(b(2)),'bsiz=',num2str(round(xbs(jj))),'p=',num2str(stats(3))]),'fontsize',7,'color','k','fontname','arial');
        axis([-xM xM -yM yM]);
        
        if ax_ind == 3 xlabel('unconditioned LC_r_s_c'); end
        if ax_ind == 1 ylabel({'ACC conditioned LC r_s_c';'LC r_s_c (ACC_h_i_g_h) -';'LC r_s_c (ACC_l_o_w)'}); end
        
    end
    
    set(gca,'fontname','arial');
    
    cd C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC;
    save ACC_LC_scatter_042421 dat_all diff_A1;
    
    %% SUPPLEMENT:
    
    %     figureNumber = 45; num = 45; wid = 14; hts = [10]; cols = {1}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,1,2, [12], 'Joshi & Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    %     movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    %
    %     ms = [5 7 9 11 7];
    %
    %     % Now Plot.
    %     axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %
    %     for kk = 1:10
    %
    %         sv1 = hi_slope(1,kk); sv1_s = hi_slope_stat(1,kk);
    %         if sv1_s<0.05 plot(xbs(kk),sv1,'ko','markersize',ms(1),'linewidth',2); end
    %         if sv1_s>0.05 plot(xbs(kk),sv1,'ko','markersize',ms(1),'linewidth',.5); end
    %
    %         sv2 = hi_slope(2,kk); sv2_s = hi_slope_stat(2,kk);
    %         if sv2_s<0.05 plot(xbs(kk),sv2,'ko','markersize',ms(2),'linewidth',2); end
    %         if sv2_s>0.05 plot(xbs(kk),sv2,'ko','markersize',ms(2),'linewidth',.5); end
    %
    %         sv3 = hi_slope(3,kk); sv3_s = hi_slope_stat(3,kk);
    %         if sv3_s<0.05 plot(xbs(kk),sv3,'ko','markersize',ms(3),'linewidth',2); end
    %         if sv3_s>0.05 plot(xbs(kk),sv3,'ko','markersize',ms(3),'linewidth',.5); end
    %
    %         sv4 = hi_slope(4,kk); sv4_s = hi_slope_stat(4,kk);
    %         if sv4_s<0.05 plot(xbs(kk),sv4,'ko','markersize',ms(4),'linewidth',2); end
    %         if sv4_s>0.05 plot(xbs(kk),sv4,'ko','markersize',ms(4),'linewidth',.5); end
    %
    %     end
    %
    %     plot([50 1050],[0 0],'k--');
    %     axis([50 1050 -0.25 0.4]);
    %     xlabel('binsize (ms)');
    %     ylabel('slope (\DeltaACC rsc vs unconditioned ACC rsc)');
    %     % axis square;
    %
    %     % pause
    
    
end

