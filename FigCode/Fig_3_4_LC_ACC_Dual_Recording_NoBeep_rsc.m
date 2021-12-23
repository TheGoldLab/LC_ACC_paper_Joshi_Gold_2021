% Fig_3_4_LC_ACC_Dual_Recording_NoBeep_rsc
%
% INFO: Do stats, summarize and plot the rsc data that was organized in "Fig_2_LC_ACC_Basics.m".
%
% Added statistical testing for within mnky data.

clear; clear all;

%% Example binned rsc: Scatter plots.

load('LC_noBeep_rscData_031119'); % Example data only!
rsc_ze = [unit_all_ze{1}; unit_all_ze{2}];
rsc_hi = [unit_all_hi{1}; unit_all_hi{2}];

xbs = linspace(100,1000,10); % 032320.
nBins = length(xbs);

figureNumber = 3; num = 3; wid = 14; hts = [6 6]; cols = {2 1}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,3,4, [12], 'Joshi & Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
movegui(fig_,[1,1]); % Move figure so it is visible on screen.

% Plot example scatter.

axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

c1='k';
% scatter(eg_ze{1},eg_ze{3},6,'markerfacecolor','m','markeredgecolor','m');
% scatter(eg_hi{1},eg_hi{3},6,'markerfacecolor','none','markeredgecolor','k');
plot(eg_ze{1},eg_ze{3},'ms','markerfacecolor','m','markeredgecolor','none','markersize',5);
plot(eg_hi{1},eg_hi{3},'kd','markersize',5);

legend('LC zero','LC nonzero','autoupdate','off','fontsize',7,'location','northwest','fontname','arial'); % legend('boxoff');

% pause;

br = robustfit(eg_ze{1},eg_ze{3});
xv = linspace(min(eg_ze{1}),max(eg_ze{1}),10);
fitY = br(1)+br(2)*xv;
plot(xv,fitY,'k-','linewidth',0.5);
% text(5,11,num2str(br(2)),'color','k','fontsize',7);
arrr = corrcoef(eg_ze{1},eg_ze{3});
text(5,10,strcat('rsc ze=',num2str(arrr(1,2))),'color','k','fontsize',7,'fontname','arial');

% scatter(eg_hi{1},eg_hi{3},5,'markerfacecolor',c1,'markeredgecolor','none');
br = robustfit(eg_hi{1},eg_hi{3});
xv = linspace(min(eg_hi{1}),max(eg_hi{1}),10);
fitY = br(1)+br(2)*xv;
% text(5,9,num2str(br(2)),'color','r','fontsize',7);
arrr = corrcoef(eg_hi{1},eg_hi{3});
text(5,8,strcat('rsc nz=',num2str(arrr(1,2))),'color','k','fontsize',8,'fontname','arial');
plot(xv,fitY,'k-','linewidth',2);
axis([0 11 0 11]);
xlabel('spike count, unit 1','fontsize',10,'fontname','arial');
ylabel('spike count, unit 2','fontsize',10,'fontname','arial');
% title(['ACC; binsize=',num2str(278),'ms'],'fontsize',7);
% title('ACC','fontsize',10);

axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

% scatter(eg_ze{2},eg_ze{4},12,'markerfacecolor','none','markeredgecolor',c1);
% scatter(eg_hi{2},eg_hi{4},5,'markerfacecolor',c1,'markeredgecolor','none');
plot(eg_ze{2},eg_ze{4},'ms','markerfacecolor','m','markeredgecolor','none','markersize',5);
plot(eg_hi{2},eg_hi{4},'kd','markersize',5);

br = robustfit(eg_ze{2},eg_ze{4});
xv = linspace(min(eg_ze{2}),max(eg_ze{2}),10);
fitY = br(1)+br(2)*xv;
plot(xv,fitY,'k--','linewidth',0.5);
% text(10,20,num2str(br(2)),'color','k','fontsize',7);
arrr = corrcoef(eg_ze{2},eg_ze{4});
text(10,20,strcat('rsc ze=',num2str(arrr(1,2))),'color','k','fontsize',8,'fontname','arial');

% scatter(eg_hi{2},eg_hi{4},5,'markerfacecolor',c1,'markeredgecolor','none');
br = robustfit(eg_hi{2},eg_hi{4});
xv = linspace(min(eg_hi{2}),max(eg_hi{2}),10);
fitY = br(1)+br(2)*xv;
% text(10,16,num2str(br(2)),'color','r','fontsize',7);
arrr = corrcoef(eg_hi{2},eg_hi{4});
text(10,16,strcat('rsc nz=',num2str(arrr(1,2))),'color','k','fontsize',8,'fontname','arial');
plot(xv,fitY,'k-','linewidth',2);
axis([0 20 0 20]);
xlabel('spike count, unit 1','fontsize',10,'fontname','arial');
ylabel('spike count, unit 2','fontsize',10,'fontname','arial');
% title(['binsize=',num2str(1000),'ms'],'fontsize',7);

% pause

% Example binned rsc.

axes(axs(3)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax); hold on;
kk = 1112; example_ze = rsc_ze(kk,:); example_hi = rsc_hi(kk,:);
plot(xbs,example_ze,'ms','linewidth',1,'markerfacecolor','m','markeredgecolor','none','markersize',5);
plot(xbs,example_hi,'kd','linewidth',1,'markersize',5);
plot([50 1050],[0 0],'k--','linewidth',0.5);
legend({'LC zero','LC nonzero'},'autoupdate','off','fontsize',8,'location','northwest','fontname','arial'); legend('boxoff');
xlabel('binsize','fontsize',10,'fontname','arial');
ylabel('rsc','fontsize',10,'fontname','arial');
xlim([50 1050]); ylim([-0.1 0.52]);


%% Population average binned rsc.

% bin = linspace(100,700,10); % 021620
% bin = linspace(100,700,5); % 021620
% bin = linspace(100,700,7); % 032220

% xbs = linspace(100,1000,10); % 032220
xbs = linspace(200,1000,5); % 032220b.

nBins = length(xbs);

% cd ('C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2020\LC_ACC'); % Base directory for brain area.
% load('LC_noBeep_rscData_032320');
% load('LC_noBeep_rscData_042020');
% load('LC_noBeep_rscData_050120');
% load('LC_noBeep_rscData_071820b');
% load('LC_noBeep_rscData_042020b');

% load('LC_noBeep_rscData_071820a');
% load('LC_noBeep_rscData_100120');

cd ('C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC'); % Base directory for brain area.

% load('LC_noBeep_rscData_042721b');
% load('LC_noBeep_rscData_042821d');
% load('LC_noBeep_rscData_042921d');
% load('LC_noBeep_rscData_043021b');
% load('LC_noBeep_rscData_052521');
load('LC_noBeep_rscData_100421');

all_mnk_dat = [];
all_rate_diff = [];
all_FF_diff = [];
all_znz_diff = [];
all_znz_diff2 = [];

for mm = 1:2
    
    giA = []; giA2 = [];
    scatter_allA1 = []; scatter_zeA1 = []; scatter_hiA1 = [];
    scatter_allA2 = []; scatter_zeA2 = []; scatter_hiA2 = [];
    scatter_allA3 = []; scatter_zeA3 = []; scatter_hiA3 = [];
    scatter_allA4 = []; scatter_zeA4 = []; scatter_hiA4 = [];
    all_datA_b = []; all_datZ_b = [];
    all_datAZ_b = []; all_datZ_b = [];
    
    for kk = 1:nBins
        
        all_datZ_temp  = unit_all_ze{mm}(:,kk);
        all_datA_temp  = unit_all_all{mm}(:,kk);
        all_datNZ_temp  = unit_all_hi{mm}(:,kk);
        
        if kk == 4
            all_rate_ze = nanmean(unit_dat_rate_ze{mm},2);
            all_rate_nz = nanmean(unit_dat_rate_hi{mm},2);
            all_FF_ze = nanmean(unit_dat_var_ze{mm},2);
            all_FF_nz = nanmean(unit_dat_var_hi{mm},2);
            all_rate_diffT =  (all_rate_nz - all_rate_ze); % ./abs(all_rate_ze(giA));
            all_FF_diffT =  (all_FF_nz - all_FF_ze); % ./abs(all_rate_ze(giA));
            all_znz_diffT = (all_datNZ_temp - all_datZ_temp); % ./(abs(all_datZ_temp(giA)));
            all_znz_diffT2 = (all_datNZ_temp - all_datZ_temp); % ./(abs(all_datZ_temp(giA)));
            giA = find(isfinite(all_rate_diffT) & isfinite(all_znz_diffT));
            giA2 = find(isfinite(all_FF_diffT) & isfinite(all_znz_diffT));
            
            all_rate_diff{mm} = all_rate_diffT(giA);
            all_FF_diff{mm} = all_FF_diffT(giA);
            all_znz_diff{mm} = all_znz_diffT(giA);
            all_znz_diff2{mm} = all_znz_diffT2(giA2);
            
            % giA = [];
        end
        
        % All LC >0, for terciles.
        % giA = find(isfinite(all_datA_temp) & abs(all_datA_temp)<1); all_datA_b{kk} = all_datA_temp(giA); giA = [];
        
        % All LC > 0, for LC zero terciles.
        % giA = find(isfinite(all_datA_temp) & abs(all_datA_temp)<1 & isfinite(all_datZ_temp) & abs(all_datZ_temp)<1);
        % all_datZ_b{kk} = all_datZ_temp(giA); all_datAZ_b{kk} = all_datA_temp(giA); giA = [];
        
        % All LC > 0, for LC zero/nonzero comparison.
        %         giA = find(isfinite(all_datA_temp) & abs(all_datA_temp)<1 & isfinite(all_datZ_temp) & abs(all_datZ_temp)<1 & isfinite(all_datNZ_temp) & abs(all_datNZ_temp)<1);
        %         all_datZ_b2{kk} = all_datZ_temp(giA); all_datAZ_b2{kk} = all_datA_temp(giA); all_datNZ_b2{kk} = all_datNZ_temp(giA);
        all_datZ_b2{kk} = all_datZ_temp; all_datAZ_b2{kk} = all_datA_temp; all_datNZ_b2{kk} = all_datNZ_temp;
        
        % LC = 1.
        all_datA11T = unit_all_all{mm}(:,kk); zedatA11T = unit_all_ze{mm}(:,kk); hidatA11T = unit_all_1{mm}(:,kk);
        % giA1 = find(isfinite(zedatA11T) & isfinite(hidatA11T) & isfinite(all_datA11T) & abs(zedatA11T)<1 & abs(hidatA11T)<1  & abs(all_datA11T)<1);
        % zedatA11 = zedatA11T(giA1); hidatA11 = hidatA11T(giA1); all_datA11 = all_datA11T(giA1); giA1 = [];
        zedatA11 = zedatA11T; hidatA11 = hidatA11T; all_datA11 = all_datA11T;
        
        % LC = 2.
        all_datA21T = unit_all_all{mm}(:,kk); zedatA21T = unit_all_ze{mm}(:,kk); hidatA21T = unit_all_2{mm}(:,kk);
        % giA2 = find(isfinite(zedatA21T) & isfinite(hidatA21T) & isfinite(all_datA21T) & abs(zedatA21T)<1 & abs(hidatA21T)<1  & abs(all_datA21T)<1);
        % zedatA21 = zedatA21T(giA2); hidatA21 = hidatA21T(giA2); all_datA21 = all_datA21T(giA2); giA2 = [];
        zedatA21 = zedatA21T; hidatA21 = hidatA21T; all_datA21 = all_datA21T;
        
        % LC = 3.
        all_datA31T = unit_all_all{mm}(:,kk); zedatA31T = unit_all_ze{mm}(:,kk); hidatA31T = unit_all_3{mm}(:,kk);
        % giA3 = find(isfinite(zedatA31T) & isfinite(hidatA31T) & isfinite(all_datA31T) & abs(zedatA31T)<1 & abs(hidatA31T)<1  & abs(all_datA31T)<1);
        % zedatA31 = zedatA31T(giA3); hidatA31 = hidatA31T(giA3); all_datA31 = all_datA31T(giA3); giA3 = [];
        zedatA31 = zedatA31T; hidatA31 = hidatA31T; all_datA31 = all_datA31T;
        
        % LC = 4.
        all_datA41T = unit_all_all{mm}(:,kk); zedatA41T = unit_all_ze{mm}(:,kk); hidatA41T = unit_all_4{mm}(:,kk);
        % giA4 = find(isfinite(zedatA41T) & isfinite(hidatA41T) & isfinite(all_datA41T) & abs(zedatA41T)<1 & abs(hidatA41T)<1  & abs(all_datA41T)<1);
        % zedatA41 = zedatA41T(giA4); hidatA41 = hidatA41T(giA4); all_datA41 = all_datA41T(giA4); giA4 = [];
        zedatA41 = zedatA41T; hidatA41 = hidatA41T; all_datA41 = all_datA41T;
        
        scatter_allA1{kk} = all_datA11; scatter_zeA1{kk} = zedatA11; scatter_hiA1{kk} = hidatA11;
        scatter_allA2{kk} = all_datA21; scatter_zeA2{kk} = zedatA21; scatter_hiA2{kk} = hidatA21;
        scatter_allA3{kk} = all_datA31; scatter_zeA3{kk} = zedatA31; scatter_hiA3{kk} = hidatA31;
        scatter_allA4{kk} = all_datA41; scatter_zeA4{kk} = zedatA41; scatter_hiA4{kk} = hidatA41;
        
    end
    
    all_mnk_dat{mm} = {scatter_allA1 scatter_zeA1 scatter_hiA1 scatter_allA2 scatter_zeA2 scatter_hiA2 ...
        scatter_allA3 scatter_zeA3 scatter_hiA3 scatter_allA4 scatter_zeA4 scatter_hiA4 ...
        all_datA_b all_datZ_b all_datAZ_b all_datZ_b2 all_datAZ_b2 all_datNZ_b2}; % 13-18.
    
end

% pause

%% For paper text stats:

% all_100 = [all_mnk_dat{1}{13}{1}; all_mnk_dat{2}{13}{1}];
% all_1000 = [all_mnk_dat{1}{13}{10}; all_mnk_dat{2}{13}{10}];

all_200 = [all_mnk_dat{1}{17}{1}; all_mnk_dat{2}{17}{1}];
all_1000 = [all_mnk_dat{1}{17}{5}; all_mnk_dat{2}{17}{5}];

prc_200 = prctile(all_200,[25 50 75])
prc_1000 = prctile(all_1000,[25 50 75])

% aDat1 = all_mnk_dat{1}{13}{1}; % gti = find(abs(aDat1)>0); aDat1 = aDat1(gti);
% aDat2 = all_mnk_dat{2}{13}{1}; % gti = find(abs(aDat2)>0); aDat2 = aDat2(gti);
% prc_all1 = prctile([aDat1; aDat2],[25 50 75]);
%
% aDat1 = all_mnk_dat{1}{13}{10}; % gti = find(abs(aDat1)>0); aDat1 = aDat1(gti);
% aDat2 = all_mnk_dat{2}{13}{10}; % gti = find(abs(aDat2)>0); aDat2 = aDat2(gti);
% prc_all10 = prctile([aDat1; aDat2],[25 50 75]);


%% Plot rsc differences for each tercile.
plotAll = 0;
plotAll = 1;

if plotAll
    
    %% Organize data for plot.
    
    outmean_loNZD = []; outmean_midNZD = []; outmean_hiNZD = [];
    outmean_lo1D = []; outmean_mid1D = []; outmean_hi1D = [];
    outmean_lo2D = []; outmean_mid2D = []; outmean_hi2D = [];
    outmean_lo3D = []; outmean_mid3D = []; outmean_hi3D = [];
    outmean_lo4D = []; outmean_mid4D = []; outmean_hi4D = [];
    
    outse_lonz = []; outse_midnz = []; outse_hinz = [];
    outse_loZ = []; outse_midZ = []; outse_hiZ = [];
    outse_lo1 = []; outse_mid1 = []; outse_hi1 = [];
    outse_lo2 = []; outse_mid2 = []; outse_hi2 = [];
    outse_lo3 = []; outse_mid3 = []; outse_hi3 = [];
    outse_lo4 = []; outse_mid4 = []; outse_hi4 = [];
    
    outse_loNZD = []; outse_midNZD = []; outse_hiNZD = [];
    outse_lo1D = []; outse_mid1D = []; outse_hi1D = [];
    outse_lo2D = []; outse_mid2D = []; outse_hi2D = [];
    outse_lo3D = []; outse_mid3D = []; outse_hi3D = [];
    outse_lo4D = []; outse_mid4D = []; outse_hi4D = [];
    
    sigVal_lo12 = []; sigVal_mid12 = []; sigVal_hi12 = [];
    sigVal_lo22 = []; sigVal_mid22 = []; sigVal_hi22 = [];
    sigVal_lo32 = []; sigVal_mid32 = []; sigVal_hi32 = [];
    sigVal_lo42 = []; sigVal_mid42 = []; sigVal_hi42 = [];
    
    sigVal_lonz = []; sigVal_midnz = []; sigVal_hinz = [];
    sigVal_lonz1 = []; sigVal_midnz1 = []; sigVal_hinz1 = [];
    sigVal_lonz2 = []; sigVal_midnz2 = []; sigVal_hinz2 = [];
    
    % pause
    % all_mnk_dat{mm} = {scatter_allA1 scatter_zeA1 scatter_hiA1 scatter_allA2 scatter_zeA2 scatter_hiA2 scatter_allA3 scatter_zeA3 scatter_hiA3 scatter_allA4 scatter_zeA4 scatter_hiA4 all_datA_b};
    
    sigValbnzD = [];
    sigValbnz = []; sigValb1 = []; sigValb2 = []; sigValb3 = []; sigValb4 = [];
    sigValDnz = []; sigValD1 = []; sigValD2 = []; sigValD3 = []; sigValD4 = [];
    outmean_lozb = []; outmean_midzb = []; outmean_hizb = [];
    
    outmean_lonz = []; outmean_midnz = []; outmean_hinz = [];
    outmean_lonz1 = []; outmean_midnz1 = []; outmean_hinz1 = [];
    outmean_lonz2 = []; outmean_midnz2 = []; outmean_hinz2 = [];
    outmean_loz = []; outmean_midz = []; outmean_hiz = [];
    outmean_lo1 = []; outmean_mid1 = []; outmean_hi1 = [];
    outmean_lo2 = []; outmean_mid2 = []; outmean_hi2 = [];
    outmean_lo3 = []; outmean_mid3 = []; outmean_hi3 = [];
    outmean_lo4 = []; outmean_mid4 = []; outmean_hi4 = [];
    
    pv_loD_bothMnks = [];
    pv_midD_bothMnks = [];
    pv_hiD_bothMnks = [];
    
    sa_dat_nz = [];
    limVal = 99; % 95 until 042320.
    limVal2 = 1000;
    
    pvalsBM = [];
    pValsB_BothMnks = [];
    
    for jj = 1:nBins
        LC_1_lo = []; LC_1_mid = []; LC_1_hi = []; LC_2_lo = []; LC_2_mid = []; LC_2_hi = []; LC_3_lo = []; LC_3_mid = []; LC_3_hi = []; LC_4_lo = []; LC_4_mid = []; LC_4_hi = [];
        LC_1_loD = []; LC_1_midD = []; LC_1_hiD = []; LC_2_loD = []; LC_2_midD = []; LC_2_hiD = []; LC_3_loD = []; LC_3_midD = []; LC_3_hiD = []; LC_4_loD = []; LC_4_midD = []; LC_4_hiD = [];
        LC_NZ_loD = []; LC_NZ_midD = []; LC_NZ_hiD = [];
        LC_Z_lo = []; LC_Z_mid = []; LC_Z_hi = []; LC_NZ_lo = []; LC_NZ_mid = []; LC_NZ_hi = [];
        sa_dat_nz_m = [];
        
        % Per monkey ANOVA for testing LC spiking level contribution to rate comparison with LC zero.
        pValsB = [];
        
        dat1_LCspInd_BothMnks = [];
        datD_BothMnks = [];
        
        for mm = 1:2
            
            datD = []; p1A=[]; p2A=[];
            %
            %             aDat = all_mnk_dat{mm}{13}{jj}; % gti = find(abs(aDat)>0); terc_vals = quantile(aDat(gti),2);
            %             aDat = all_mnk_dat{mm}{14}{jj}; aDatA = all_mnk_dat{mm}{15}{jj};
            %             t1 = find(aDatA<=terc_vals(1)); t2 = find(aDatA>terc_vals(1) & aDatA<terc_vals(2)); t3 = find(aDatA>=terc_vals(2));
            %             LC_Z_lo = [LC_Z_lo; aDat(t1)]; LC_Z_mid = [LC_Z_mid;  aDat(t2)]; LC_Z_hi = [LC_Z_hi;  aDat(t3)];
            
            % ********************
            % LC NZ.
            % Central measures for data.
            D1 = []; D2 = []; D3 = [];
            aDat = all_mnk_dat{mm}{17}{jj}; mdat = all_mnk_dat{mm}{18}{jj}; zdat = all_mnk_dat{mm}{16}{jj};
            % gti = find(isfinite(aDat) & isfinite(mdat) & isfinite(zdat));
            terc_vals = quantile(aDat,2);
            t1 = find(aDat<=terc_vals(1)); t2 = find(aDat>terc_vals(1) & aDat<terc_vals(2)); t3 = find(aDat>=terc_vals(2));
            
            LC_Z_lo = [LC_Z_lo; zdat(t1)]; LC_Z_mid = [LC_Z_mid; zdat(t2)]; LC_Z_hi = [LC_Z_hi; zdat(t3)];
            LC_NZ_lo = [LC_NZ_lo; mdat(t1)]; LC_NZ_mid = [LC_NZ_mid; mdat(t2)]; LC_NZ_hi = [LC_NZ_hi; mdat(t3)];
            
            A1=aDat(t1); A2 = aDat(t2); A3 = aDat(t3);
            D1=mdat(t1)-zdat(t1); D2 = mdat(t2)-zdat(t2); D3 = mdat(t3)-zdat(t3);
            LC_NZ_loD = [LC_NZ_loD; D1]; LC_NZ_midD = [LC_NZ_midD; D2]; LC_NZ_hiD = [LC_NZ_hiD; D3]; % pause
            p1 = ranksum(mdat(t1),zdat(t1)); p2= ranksum(mdat(t2),zdat(t2)); p3 = ranksum(mdat(t3),zdat(t3));
            sigVal_lonz(mm) = p1; sigVal_midnz(mm) = p2; sigVal_hinz(mm) = p3;
            sa_dat_nz_m{mm} = [mdat(t3) zdat(t3)];
            p1 = signrank(D1(isfinite(D1))); p2= signrank(D2(isfinite(D2))); p3 = signrank(D3(isfinite(D3)));
            sigVal_loZD(mm) = p1; sigVal_midZD(mm) =  p2; sigVal_hiZD(mm) = p3;
            LC_sp_ind_ze = 1*ones(length(D1),1);
            datD_lo{mm} = D1;
            datD_mid{mm} = D2;
            datD_hi{mm} = D3;
            datD_all{mm} = [D1; D2; D3];
            datA_all{mm} = [A1; A2; A3];
            
            datD_all2{mm} = mdat-zdat;
            datA_all2{mm} = aDat;
            
            % ********************
            
            % ********************
            % LC = 1.
            % Central measures for data.
            % aDat = all_mnk_dat{mm}{1}{jj}; mdat = all_mnk_dat{mm}{3}{jj}; zdat = all_mnk_dat{mm}{2}{jj};
            aDat = all_mnk_dat{mm}{17}{jj}; mdat = all_mnk_dat{mm}{3}{jj}; zdat = all_mnk_dat{mm}{2}{jj};
            terc_vals = quantile(aDat,2);
            t1 = find(aDat<=terc_vals(1)); t2 = find(aDat>terc_vals(1) & aDat<terc_vals(2)); t3 = find(aDat>=terc_vals(2));
            LC_1_lo = [LC_1_lo; mdat(t1)]; LC_1_mid = [LC_1_mid; mdat(t2)]; LC_1_hi = [LC_1_hi; mdat(t3)];
            p1 = ranksum(mdat(t1),zdat(t1)); p2= ranksum(mdat(t2),zdat(t2)); p3 = ranksum(mdat(t3),zdat(t3));
            sigVal_lo1(mm) = p1; sigVal_mid1(mm) = p2; sigVal_hi1(mm) = p3;
            D1 = []; D2 = []; D3 = [];
            D1=mdat(t1)-zdat(t1); D2 = mdat(t2)-zdat(t2); D3 = mdat(t3)-zdat(t3);
            D1 = D1(isfinite(D1)); D2 = D2(isfinite(D2)); D3 = D3(isfinite(D3));
            LC_sp_ind_1 = 1*ones(length(D3),1);
            datD = [datD; D3];
            
            % ********************
            % LC = 2.
            % Central measures for data.
            % aDat = all_mnk_dat{mm}{4}{jj}; mdat = all_mnk_dat{mm}{6}{jj}; zdat = all_mnk_dat{mm}{5}{jj};
            aDat = all_mnk_dat{mm}{17}{jj}; mdat = all_mnk_dat{mm}{6}{jj}; zdat = all_mnk_dat{mm}{5}{jj};
            terc_vals = quantile(aDat,2);
            t1 = find(aDat<=terc_vals(1)); t2 = find(aDat>terc_vals(1) & aDat<terc_vals(2)); t3 = find(aDat>=terc_vals(2));
            LC_2_lo = [LC_2_lo; mdat(t1)]; LC_2_mid = [LC_2_mid; mdat(t2)]; LC_2_hi = [LC_2_hi; mdat(t3)];
            p1 = ranksum(mdat(t1),zdat(t1)); p2= ranksum(mdat(t2),zdat(t2)); p3 = ranksum(mdat(t3),zdat(t3));
            sigVal_lo2(mm) = p1; sigVal_mid2(mm) = p2; sigVal_hi2(mm) = p3;
            D1 = []; D2 = []; D3 = [];
            D1=mdat(t1)-zdat(t1); D2 = mdat(t2)-zdat(t2); D3 = mdat(t3)-zdat(t3);
            D1 = D1(isfinite(D1)); D2 = D2(isfinite(D2)); D3 = D3(isfinite(D3));
            LC_sp_ind_2 = 1*ones(length(D3),1);
            datD = [datD; D3];
            
            % ********************
            % LC = 3.
            % Central measures for data.
            % aDat = all_mnk_dat{mm}{7}{jj}; mdat = all_mnk_dat{mm}{9}{jj}; zdat = all_mnk_dat{mm}{8}{jj};
            aDat = all_mnk_dat{mm}{17}{jj}; mdat = all_mnk_dat{mm}{9}{jj}; zdat = all_mnk_dat{mm}{8}{jj};
            terc_vals = quantile(aDat,2);
            t1 = find(aDat<=terc_vals(1)); t2 = find(aDat>terc_vals(1) & aDat<terc_vals(2)); t3 = find(aDat>=terc_vals(2));
            LC_3_lo = [LC_3_lo; mdat(t1)]; LC_3_mid = [LC_3_mid; mdat(t2)]; LC_3_hi = [LC_3_hi; mdat(t3)];
            p1 = ranksum(mdat(t1),zdat(t1)); p2= ranksum(mdat(t2),zdat(t2)); p3 = ranksum(mdat(t3),zdat(t3));
            sigVal_lo3(mm) = p1; sigVal_mid3(mm) = p2; sigVal_hi3(mm) = p3;
            D1 = []; D2 = []; D3 = [];
            D1=mdat(t1)-zdat(t1); D2 = mdat(t2)-zdat(t2); D3 = mdat(t3)-zdat(t3);
            D1 = D1(isfinite(D1)); D2 = D2(isfinite(D2)); D3 = D3(isfinite(D3));
            LC_sp_ind_3 = 1*ones(length(D3),1);
            datD = [datD; D3];
            
            % ********************
            % LC = 4.
            % Central measures for data.
            % aDat = all_mnk_dat{mm}{10}{jj}; mdat = all_mnk_dat{mm}{12}{jj}; zdat = all_mnk_dat{mm}{11}{jj};
            aDat = all_mnk_dat{mm}{17}{jj}; mdat = all_mnk_dat{mm}{12}{jj}; zdat = all_mnk_dat{mm}{11}{jj};
            terc_vals = quantile(aDat,2);
            t1 = find(aDat<=terc_vals(1)); t2 = find(aDat>terc_vals(1) & aDat<terc_vals(2)); t3 = find(aDat>=terc_vals(2));
            LC_4_lo = [LC_4_lo; mdat(t1)]; LC_4_mid = [LC_4_mid; mdat(t2)]; LC_4_hi = [LC_4_hi; mdat(t3)];
            p1 = ranksum(mdat(t1),zdat(t1)); p2= ranksum(mdat(t2),zdat(t2)); p3 = ranksum(mdat(t3),zdat(t3));
            sigVal_lo4(mm) = p1; sigVal_mid4(mm) = p2; sigVal_hi4(mm) = p3;
            D1 = []; D2 = []; D3 = [];
            D1=mdat(t1)-zdat(t1); D2 = mdat(t2)-zdat(t2); D3 = mdat(t3)-zdat(t3);
            D1 = D1(isfinite(D1)); D2 = D2(isfinite(D2)); D3 = D3(isfinite(D3));
            LC_sp_ind_4 = 1*ones(length(D3),1);
            datD = [datD; D3];
            
            % *************************************
            % *************************************
            % Set up ANOVA.
            % Number of data samples for each LC spike condition and each monkey.
            % nze = LC_sp_ind_ze;
            n1 = LC_sp_ind_1;
            n2 = LC_sp_ind_2;
            n3 = LC_sp_ind_3;
            n4 = LC_sp_ind_4;
            
            % LC=0,1,2,3,>3 correspond to indices 1,2,3,4,5.
            % LC_sp_ind_ze = 1*nze;
            LC_sp_ind_1 = 2*n1;
            LC_sp_ind_2 = 3*n2;
            LC_sp_ind_3 = 4*n3;
            LC_sp_ind_4 = 5*n4;
            
            %             pause
            
            % dat1_LCspInd = [LC_sp_ind_ze; LC_sp_ind_1; LC_sp_ind_2; LC_sp_ind_3; LC_sp_ind_4];
            dat1_LCspInd = [LC_sp_ind_1; LC_sp_ind_2; LC_sp_ind_3; LC_sp_ind_4];
            dat1_LCspInd_BothMnks = [dat1_LCspInd_BothMnks; LC_sp_ind_1; LC_sp_ind_2; LC_sp_ind_3; LC_sp_ind_4];
            datD_BothMnks = [datD_BothMnks; datD];
            % 1-way ANOVA with display OFF.
            [pA,tbl,stats] = anova1(datD, dat1_LCspInd,'off');
            pValsB = [pValsB; pA];
            
            % *************************************
            % *************************************
            
        end % Mnks for this bin.
        
        binDiffVals{jj} = {datD_lo datD_mid datD_hi};
        
        % 1-way ANOVA with display OFF.
        [pA,tbl,stats] = anova1(datD_BothMnks, dat1_LCspInd_BothMnks,'off');
        pValsB_BothMnks = [pValsB_BothMnks; pA];
        
        pvalsBM = [pvalsBM pValsB];
        
        sa_dat_nz{jj} = sa_dat_nz_m;
        
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
        
        %         outmean_lonz2 = [outmean_lonz2 nanmedian(datD_lo{2})];
        %         outmean_midnz2 = [outmean_midnz2 nanmedian(datD_mid{2})];
        %         outmean_hinz2 = [outmean_hinz2 nanmedian(datD_hi{2})];
        
        
        %         [h p] = ttest(ci_lo1_m); pv_lo(jj) = p;
        %         [h p] = ttest(ci_mid1_m); pv_mid(jj) = p;
        %         [h p] = ttest(ci_hi1_m); pv_hi(jj) = p;
        
        
        %         pv_lo(jj) = signrank(LC_NZ_loD2(LC_NZ_loD2<prctile(LC_NZ_loD2,limVal)));
        %         pv_mid(jj) = signrank(LC_NZ_midD2(LC_NZ_midD2<prctile(LC_NZ_midD2,limVal)));
        %         pv_hi(jj) = signrank(LC_NZ_hiD2(LC_NZ_hiD2<prctile(LC_NZ_hiD2,limVal)));
        
        % Tests
        pv_loD(jj) = signrank(LC_NZ_loD);
        pv_midD(jj) = signrank(LC_NZ_midD);
        pv_hiD(jj) = signrank(LC_NZ_hiD);
        
        %         pv_loD(jj) = signrank([LC_1_loD; LC_2_loD; LC_3_loD; LC_4_loD]);
        %         pv_midD(jj) = signrank([LC_1_midD; LC_2_midD; LC_3_midD; LC_4_midD]);
        %         pv_hiD(jj) = signrank([LC_1_hiD; LC_2_hiD; LC_3_hiD; LC_4_hiD]);
        %
        %         pv_loD_bothMnks(jj,:) = [signrank(LC_1_loD); signrank(LC_2_loD); signrank(LC_3_loD); signrank(LC_4_loD)];
        %         pv_midD_bothMnks(jj,:) = [signrank(LC_1_midD); signrank(LC_2_midD); signrank(LC_3_midD); signrank(LC_4_midD)];
        %         pv_hiD_bothMnks(jj,:) = [signrank(LC_1_hiD); signrank(LC_2_hiD); signrank(LC_3_hiD); signrank(LC_4_hiD)];
        
        %         %                 pv_lo(jj) = signrank([LC_1_loD; LC_2_loD; LC_3_loD; LC_4_loD]);
        %         %                 pv_mid(jj) = signrank([LC_1_midD; LC_2_midD; LC_3_midD; LC_4_midD]);
        %         %                 pv_hi(jj) = signrank([LC_1_hiD; LC_2_hiD; LC_3_hiD; LC_4_hiD]);
        
        
        
        pv_lo(jj) = ranksum(LC_Z_lo,LC_NZ_lo);
        pv_mid(jj) = ranksum(LC_Z_mid,LC_NZ_mid);
        pv_hi(jj) = ranksum(LC_Z_hi,LC_NZ_hi);
        
        % ********************
        
        % ********************
        % LC =1.
        % Central measures for data.
        LC_1_lo = LC_1_lo(isfinite(LC_1_lo)); LC_1_mid = LC_1_mid(isfinite(LC_1_mid)); LC_1_hi = LC_1_hi(isfinite(LC_1_hi));
        outmean_lo1 = [outmean_lo1 nanmean(LC_1_lo)]; outmean_mid1 = [outmean_mid1 nanmean(LC_1_mid)]; outmean_hi1 = [outmean_hi1 nanmean(LC_1_hi)];
        [ci_lo1 ci_lo1_m]= bootci(2000,{@mean,LC_1_lo},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_mid1 ci_mid1_m]= bootci(2000,{@mean,LC_1_mid},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_hi1 ci_hi1_m]= bootci(2000,{@mean,LC_1_hi},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        outse_lo1 = [outse_lo1 ci_lo1]; outse_mid1 = [outse_mid1 ci_mid1]; outse_hi1 = [outse_hi1 ci_hi1];
        
        % LC =2.
        % Central measures for data.
        LC_2_lo = LC_2_lo(isfinite(LC_2_lo)); LC_2_mid = LC_2_mid(isfinite(LC_2_mid)); LC_2_hi = LC_2_hi(isfinite(LC_2_hi));
        outmean_lo2 = [outmean_lo2 nanmean(LC_2_lo)]; outmean_mid2 = [outmean_mid2 nanmean(LC_2_mid)]; outmean_hi2 = [outmean_hi2 nanmean(LC_2_hi)];
        [ci_lo2 ci_lo2_m]= bootci(2000,{@mean,LC_2_lo},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_mid2 ci_mid2_m]= bootci(2000,{@mean,LC_2_mid},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_hi2 ci_hi2_m]= bootci(2000,{@mean,LC_2_hi},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        outse_lo2 = [outse_lo2 ci_lo2]; outse_mid2 = [outse_mid2 ci_mid2]; outse_hi2 = [outse_hi2 ci_hi2];
        
        % ********************
        % LC =3.
        % Central measures for data.
        LC_3_lo = LC_3_lo(isfinite(LC_3_lo)); LC_3_mid = LC_3_mid(isfinite(LC_3_mid)); LC_3_hi = LC_3_hi(isfinite(LC_3_hi));
        outmean_lo3 = [outmean_lo3 nanmean(LC_3_lo)]; outmean_mid3 = [outmean_mid3 nanmean(LC_3_mid)]; outmean_hi3 = [outmean_hi3 nanmean(LC_3_hi)];
        [ci_lo3 ci_lo3_m]= bootci(2000,{@mean,LC_3_lo},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_mid3 ci_mid3_m]= bootci(2000,{@mean,LC_3_mid},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_hi3 ci_hi3_m]= bootci(2000,{@mean,LC_3_hi},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        outse_lo3 = [outse_lo3 ci_lo3]; outse_mid3 = [outse_mid3 ci_mid3]; outse_hi3 = [outse_hi3 ci_hi3];
        
        % ********************
        % LC >= 4.
        % Central measures for data.
        LC_4_lo = LC_4_lo(isfinite(LC_4_lo)); LC_4_mid = LC_4_mid(isfinite(LC_4_mid)); LC_4_hi = LC_4_hi(isfinite(LC_4_hi));
        outmean_lo4 = [outmean_lo4 nanmean(LC_4_lo)]; outmean_mid4 = [outmean_mid4 nanmean(LC_4_mid)]; outmean_hi4 = [outmean_hi4 nanmean(LC_4_hi)];
        [ci_lo4 ci_lo4_m]= bootci(2000,{@mean,LC_4_lo},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_mid4 ci_mid4_m]= bootci(2000,{@mean,LC_4_mid},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        [ci_hi4 ci_hi4_m]= bootci(2000,{@mean,LC_4_hi},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
        outse_lo4 = [outse_lo4 ci_lo4]; outse_mid4 = [outse_mid4 ci_mid4]; outse_hi4 = [outse_hi4 ci_hi4];
        
        % **********************************
        % Collect all sig vals.
        sigValbnz = [sigValbnz [sigVal_lonz'; sigVal_midnz'; sigVal_hinz']];
        sigValDnz = [sigValDnz [sigVal_loZD'; sigVal_midZD'; sigVal_hiZD']];
        
        sigValb1 = [sigValb1 [sigVal_lo1'; sigVal_mid1'; sigVal_hi1']];
        sigValb2 = [sigValb2 [sigVal_lo2'; sigVal_mid2'; sigVal_hi2']];
        sigValb3 = [sigValb3 [sigVal_lo3'; sigVal_mid3'; sigVal_hi3']];
        sigValb4 = [sigValb4 [sigVal_lo4'; sigVal_mid4'; sigVal_hi4']];
        
        % datD_t3{jj} = datD_hi;
        datD_A{jj} = datD_all;
        datA_A{jj} = datA_all;
        datD_A2{jj} = datD_all2;
        datA_A2{jj} = datA_all2;
        
    end % Bins.
    
    % ANOVA for effects of binsize.
    % datDD = nans(10000,1); NN  = nans(10000,1);
    datDD = []; NN  = [];
    ci = 1;
    for kk = 1:nBins
        datM1 = datD_A{kk}{1}; datM1 = datM1(isfinite(datM1));
        datM2 = datD_A{kk}{2}; datM2 = datM2(isfinite(datM2));
        NN = [NN; kk*ones(length(datM1),1); kk*ones(length(datM2),1)];
        datDD = [datDD; datM1; datM2];
    end
    %     NN(ci+1:end,:) = [];
    %     datDD(ci+1:end,:) = [];
    
    % 1-way ANOVA with display OFF.
    [pZNZ,tbl,stats] = anova1(datDD, NN,'off');
    pValsZNZ_Anova = pZNZ;
    
    %     pause
    
    %% Main Figure 4: LC linked ACC rsc.
    
    % Setup figure.
    
    figureNumber = 4; num = 4; wid = 17.6; hts = [6]; cols = {2 2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,2,2, [12], 'Joshi and Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    % Now plot.
    
    %%
    
    outmean_all=[]; sigVal_all1 = []; sigVal_all2 = []; outse_all = []; outmean_allD=[];
    
    outmean_all{1} = [outmean_loz; outmean_lo1; outmean_lo2; outmean_lo3; outmean_lo4];
    outse_all{1} = [outse_loZ; outse_lo1; outse_lo2; outse_lo3; outse_lo4];
    outmean_all{2} = [outmean_midz; outmean_mid1; outmean_mid2; outmean_mid3; outmean_mid4];
    outse_all{2} = [outse_midZ; outse_mid1; outse_mid2; outse_mid3; outse_mid4];
    outmean_all{3} = [outmean_hiz; outmean_hi1; outmean_hi2; outmean_hi3; outmean_hi4];
    outse_all{3} = [outse_hiZ; outse_hi1; outse_hi2; outse_hi3; outse_hi4];
    
    %     sigVal_all1{1} = [sigValb1(1,:);sigValb2(1,:);sigValb3(1,:);sigValb4(1,:)];
    %     sigVal_all2{1} = [sigValb1(2,:);sigValb2(2,:);sigValb3(2,:);sigValb4(2,:)];
    %     sigVal_all1{2} = [sigValb1(3,:);sigValb2(3,:);sigValb3(3,:);sigValb4(3,:)];
    %     sigVal_all2{2} = [sigValb1(4,:);sigValb2(4,:);sigValb3(4,:);sigValb4(4,:)];
    %     sigVal_all1{3} = [sigValb1(5,:);sigValb2(5,:);sigValb3(5,:);sigValb4(5,:)];
    %     sigVal_all2{3} = [sigValb1(6,:);sigValb2(6,:);sigValb3(6,:);sigValb4(6,:)];
    %
    %     sigValA1{1} = sigValbnz(1,:);
    %     sigValA2{1} = sigValbnz(2,:);
    %     sigValA1{2} = sigValbnz(3,:);
    %     sigValA2{2} = sigValbnz(4,:);
    %     sigValA1{3} = sigValbnz(5,:);
    %     sigValA2{3} = sigValbnz(6,:);
    %
    %     sigVal{1} = pv_lo;
    %     sigVal{2} = pv_mid;
    %     sigVal{3} = pv_hi;
    
    % Create binsizes spaced equally in log10 space:
    %     xbs = round(logspace(log10(100),log10(1000),10),1);
    ms = [8 5 20 9 11];
    colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    % xVals = 100*[-0.3 -0.15 0 0.15 0.3];
    xVals = 100*[-0.4 -0.2 0 0.2 0.4];
    yMax = 0.25; yMin = -0.13;
    
    for ti = 1:3 %Terciles.
        axes(axs(2*ti -1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
        for uj = 1:nBins % Bins.
            for ui = 1:5 % LC spike=0,1,2,3,>3.
                plot(xbs(uj)+xVals(ui),outmean_all{ti}(ui,uj),'.','color',colors{ui},'markersize',25);
                plot([xbs(uj)+xVals(ui) xbs(uj)+xVals(ui)],[outse_all{ti}(2*ui - 1,uj) outse_all{ti}(2*ui,uj)],'k-');
                
                %                 if sigVal{ti}(1,uj)<0.05
                %                     plot(xbs(uj)+xVals(ui),0.21,'k*','markersize',4);
                %                 end
                
            end
            
            %             if (sigValA1{ti}(1,uj)<0.05 & sigValA2{ti}(1,uj)<0.05)
            %                 % if (sigVal{ti}(1,uj)<0.05)
            %                 plot(xbs(uj),0.23,'k+','markersize',4);
            %             end
            
        end
        title(['tercile',num2str(ti)]);
        ylabel('ACC rsc','fontname','arial','fontsize',10);
        xlabel('binsize (msec)','fontname','arial','fontsize',10);
        axis([100 1100 -0.13 0.26]);
        plot([50 1050], [0 0],'k--');% ylabel('ACC rsc diff');
    end
    
    save forShufPlot outmean_all outse_all;
    
    %         pause
    
    % **********************************************************************************
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
                plot(xbs(uj),0.03,'*','color',colors{4},'markersize',6);
            end
            
        end % Bins loop.
        
        axis([100 1100 -0.115 0.05]);
        ylabel('\DeltaACC rsc re:LCzero','fontname','arial','fontsize',10);
        xlabel('binsize (msec)','fontname','arial','fontsize',10);
    end
    
    %% Calculate proportion significant for comparison with shuffled analysis.
    
    svM1 = [sigValDA1{1};sigValDA1{2};sigValDA1{3}];
    svM2 = [sigValDA2{1};sigValDA2{2};sigValDA2{3}];
    sigiM1 = nans(size(svM1));
    sigiM2 = nans(size(svM2));
    
    for ti = 1:3
        sigiM1(ti,:) = svM1(ti,:) < 0.05;
        sigiM2(ti,:) = svM2(ti,:) < 0.05;
    end
    mnkSig{1} = sigiM1;
    mnkSig{2} = sigiM2;
    
    sig_propA = 0; sig_propB = 0;
    sig_propA = sig_propA + sum(sum(mnkSig{1}));
    sig_propB = sig_propB + sum(sum(mnkSig{2}));
    pA1 = (1/15)*(sig_propA);
    pB1 = (1/15)*(sig_propB);
    
    sig_propA = 0; sig_propB = 0;
    for ti = 1:3
        sig_propA = sig_propA + (mnkSig{1}(ti,1) * mnkSig{1}(ti,2)) + (mnkSig{1}(ti,2) * mnkSig{1}(ti,3)) + (mnkSig{1}(ti,3) * mnkSig{1}(ti,4)) + (mnkSig{1}(ti,4) * mnkSig{1}(ti,5));
        sig_propB = sig_propB + (mnkSig{2}(ti,1) * mnkSig{2}(ti,2)) + (mnkSig{2}(ti,2) * mnkSig{2}(ti,3)) + (mnkSig{2}(ti,3) * mnkSig{2}(ti,4)) + (mnkSig{2}(ti,4) * mnkSig{2}(ti,5));
    end
    pA2 = (1/15)*(sig_propA);
    pB2 = (1/15)*(sig_propB);
    
    sig_propA = 0; sig_propB = 0;
    for ti = 1:3
        sig_propA = sig_propA + (mnkSig{1}(ti,1) * mnkSig{1}(ti,2) * mnkSig{1}(ti,3)) + (mnkSig{1}(ti,2) * mnkSig{1}(ti,3) * mnkSig{1}(ti,4)) + (mnkSig{1}(ti,3) * mnkSig{1}(ti,4) * mnkSig{1}(ti,5));
        sig_propB = sig_propB + (mnkSig{2}(ti,1) * mnkSig{2}(ti,2) * mnkSig{2}(ti,3)) + (mnkSig{2}(ti,2) * mnkSig{2}(ti,3) * mnkSig{2}(ti,4)) + (mnkSig{2}(ti,3) * mnkSig{2}(ti,4) * mnkSig{2}(ti,5));
    end
    pA3 = (1/15)*(sig_propA);
    pB3 = (1/15)*(sig_propB);
    
    sig_propA = 0; sig_propB = 0;
    for ti = 1:3
        sig_propA = sig_propA + (mnkSig{1}(ti,1) * mnkSig{1}(ti,2) * mnkSig{1}(ti,3) * mnkSig{1}(ti,4)) + (mnkSig{1}(ti,2) * mnkSig{1}(ti,3) * mnkSig{1}(ti,4) * mnkSig{1}(ti,5));
        sig_propB = sig_propB + (mnkSig{2}(ti,1) * mnkSig{2}(ti,2) * mnkSig{2}(ti,3) * mnkSig{2}(ti,4)) + (mnkSig{2}(ti,2) * mnkSig{2}(ti,3) * mnkSig{2}(ti,4) * mnkSig{2}(ti,5));
    end
    pA4 = (1/15)*(sig_propA);
    pB4 = (1/15)*(sig_propB);
    
    sig_propA = 0; sig_propB = 0;
    for ti = 1:3
        sig_propA = sig_propA + (mnkSig{1}(ti,1) * mnkSig{1}(ti,2) * mnkSig{1}(ti,3) * mnkSig{1}(ti,4) * mnkSig{1}(ti,5));
        sig_propB = sig_propB + (mnkSig{2}(ti,1) * mnkSig{2}(ti,2) * mnkSig{2}(ti,3) * mnkSig{2}(ti,4) * mnkSig{2}(ti,5));
    end
    pA5 = (1/15)*(sig_propA);
    pB5 = (1/15)*(sig_propB);
    
    pa_All = [[pA1 pB1]; [pA2 pB2]; [pA3 pB3]; [pA4 pB4]; [pA5 pB5]];
    
    save propSig_ACC_051521 pa_All; % Used by: Fig_2_LC_ACC_Basics_ShuffOnly.m
    
    %% Zero-Nonzero ONLY.
    
    %     % Setup figure.
    %
    %     figureNumber = 405; num = 404; wid = 15; hts = [6]; cols = {1 1 1}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,2,2, [12], 'Joshi & Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    %     % figureNumber = 405; num = 404; wid = 15; hts = [5]; cols = {2 2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,2,2, [12], 'Joshi & Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    %     movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    %% Plot.
    
    outmean_all=[]; sigVal_all1 = []; sigVal_all2 = []; outse_all = []; outmean_allD=[];
    
    outmean_all{1} = [outmean_loz; outmean_lonz]; outse_all{1} = [outse_loZ; outse_lonz];
    sigVal_all1{1} = sigValbnz(1,:); sigVal_all2{1} = sigValbnz(2,:);
    outmean_all{2} = [outmean_midz; outmean_midnz]; outse_all{2} = [outse_midZ; outse_midnz];
    sigVal_all1{2} = sigValbnz(3,:); sigVal_all2{2} = sigValbnz(4,:);
    outmean_all{3} = [outmean_hiz; outmean_hinz]; outse_all{3} = [outse_hiZ; outse_hinz];
    sigVal_all1{3} = sigValbnz(5,:); sigVal_all2{3} = sigValbnz(6,:);
    
    outmean_allD{1} = 100*outmean_loNZD; outse_allD{1} = 100*outse_loNZD;
    outmean_allD{2} = 100*outmean_midNZD; outse_allD{2} = 100*outse_mid1D;
    outmean_allD{3} = 100*outmean_hiNZD; outse_allD{3} = 100*outse_hiNZD;
    
    % Create binsizes spaced equally in log10 space:
    %     xbs = round(logspace(log10(100),log10(1000),10),1);
    %         ms = [8 5 20 9 11];
    %         colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    %         xVals = [-0.3 -0.15 0 0.15 0.3];
    %         yMax = 0.22; yMin = -0.1;
    
    %         for ti = 1:3 %Terciles.
    %             axes(axs(ti)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %             % axes(axs(2*ti -1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %             for uj = 1:nBins % Bins.
    %                 for ui = 1:2 % LC spike=0,>0.
    %                     plot(uj+xVals(ui),outmean_all{ti}(ui,uj),'.','color',colors{ui},'markersize',ms(3));
    %                     plot([uj+xVals(ui) uj+xVals(ui)],[outse_all{ti}(2*ui - 1,uj) outse_all{ti}(2*ui,uj)],'k-');
    %                     if ui>1
    %                         if (sigVal_all1{ti}(ui-1,uj)<0.05 & sigVal_all2{ti}(ui-1,uj)<0.05)
    %                             plot(uj+xVals(ui),0.21,'k*');
    %                         end
    %                     end
    %                 end
    %             end
    %             axis([0.5 nBins+0.5 yMin yMax]); plot([0.5 nBins+0.5], [0 0],'k--');% ylabel('ACC rsc diff');
    %         end
    
    % save plotShuf_ACCrsc outmean_all outse_all;
    
    % save plotShuf_ACCrsc_081520 outmean_all outse_all;
    save plotShuf_ACCrsc_042921 outmean_all outse_all;
    
    %     pause
    
    % **********************************************************************************
    %% PLOT BARS.
    
    %     ui=1;
    %     yMax = 55; yMin = -55;
    %     for ti = 1:3 % Terciles.
    %         axes(axs(2*ti)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %         for uj = 1:nBins % Bins.
    %             % for ui = 1:4 % LC spike=1,2,3,>3 relative to LC = 0.
    %             bar(uj+xVals(ui+1),outmean_allD{ti}(ui,uj),0.15,'facecolor',colors{ui+1},'edgecolor','none');
    %             plot([uj+xVals(ui+1) uj+xVals(ui+1)],[outse_allD{ti}(2*ui-1,uj) outse_allD{ti}(2*ui,uj)],'k-');
    %             % end
    %         end
    %         axis([0.5 10.5 yMin yMax]); plot([0.5 10.5], [0 0],'k--');% ylabel('ACC rsc diff');
    %     end
    
    %     pause
    
    
    %%
    
    all_bdat = cell(1,length(xbs));
    
    for jj = 1:length(xbs)
        
        dat_all = [datA_A2{jj}{1}; datA_A2{jj}{2}];
        diff_A1 = [datD_A2{jj}{1}; datD_A2{jj}{2}];
        
        gi = find(isfinite(dat_all) & isfinite(diff_A1));
        diff_A1 = diff_A1(gi); dat_all = dat_all(gi);
        
        terc_vals = quantile(dat_all,2);
        t1 = find(dat_all<=terc_vals(1)); t2 = find(dat_all>terc_vals(1) & dat_all<terc_vals(2)); t3 = find(dat_all>=terc_vals(2));
        
        %
        
        %         bin_yd{1} = diff_A1(t1);
        %         bin_yd2{2} = diff_A1(t2);
        %         bin_yd3{3}= diff_A1(t3);
        
        all_bdat{jj} = {diff_A1(t1) diff_A1(t2) diff_A1(t3)};
        
    end
    
    cd C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC;
    save Fig4_LC_ACC_Compare all_bdat binDiffVals;
    
    
    %% NEW Fig 5 (MAYBE): Relationship between unconditioned ACC rsc and ACC rsc change between LC zero and nonzero conditions.
    figureNumber = 5; num = 5; wid = 17.6; hts = [20]; cols = {5}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,1,1, [12], 'Joshi and Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    %% Oct 2021.
    
    %         mdl1 = fitlm(xd,yd);
    %         coeffs = mdl1.Coefficients.Estimate;
    %         pvals = mdl1.Coefficients.pValue;
    %         slope = coeffs(2); interc = coeffs(1);
    %         b = [interc slope];
    %         statVal = pvals(2);
    
    all_dat_P = zeros(5,4); all_dat_R = zeros(5,4);
    slope_val = zeros(5,4); slope_stat = zeros(5,4);
    yM = 0.5; xM = 0.5;
    % yM = .5; xM = .5;
    
    
    for jj = 1:length(xbs)
        
        dat_all = [datA_A2{jj}{1}; datA_A2{jj}{2}];
        diff_A1 = [datD_A2{jj}{1}; datD_A2{jj}{2}];
        
        gi = find(isfinite(dat_all) & isfinite(diff_A1));
        diff_A1 = diff_A1(gi); dat_all = dat_all(gi);
        
        terc_vals = quantile(dat_all,2);
        t1 = find(dat_all<=terc_vals(1)); t2 = find(dat_all>terc_vals(1) & dat_all<terc_vals(2)); t3 = find(dat_all>=terc_vals(2));
        
        %
        
        %         bin_yd{1} = diff_A1(t1);
        %         bin_yd2{2} = diff_A1(t2);
        %         bin_yd3{3}= diff_A1(t3);
        
        all_bdat{jj} = {diff_A1(t1) diff_A1(t2) diff_A1(t3)};
        
        %
        
        
        axes(axs(jj)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
        plot(dat_all,diff_A1,'.','markersize',2);
        [b i] = sort(dat_all); mm1 = movmean(diff_A1(i),500);
        plot(dat_all(i),mm1,'-','linewidth',2,'color',[1 0 0]);
        
        % Fit all data.
        xd = dat_all; yd = diff_A1;
        xdat = [ones(size(xd)) xd];
        [b,bint,r,rint,stats]=regress(yd,xdat);
        statVal = stats(3);
        x_plot = linspace(min(xd),max(xd),10);
        if statVal < 0.05
            plot(x_plot,b(1)+b(2)*x_plot,'k-','linewidth',1);
        else
            plot(x_plot,b(1)+b(2)*x_plot,'k--','linewidth',1);
        end
        plot([-0.5 0.5],[0 0],'k--'); plot([0 0],[-0.5 0.5],'k--');
        text(-0.45,0.4,strcat(['bsiz=',num2str(round(xbs(jj)))]),'fontsize',7,'color','k','fontname','arial');
        text(-0.45,1,strcat(['slope all:',num2str(b(2)),'; p=',num2str(statVal)]),'fontsize',7,'color','k','fontname','arial');
        axis([-xM xM -yM yM]); axis square;
        slope_val(jj,1) = b(2); slope_stat(jj,1) = statVal;
        [r p] = corr(diff_A1,dat_all,'type','Spearman'); all_dat_R(jj,1) = r; all_dat_P(jj,1) = p;
        
        % Fit data for unconditioned rsc <= 0.
        xd = dat_all(dat_all<=0); yd = diff_A1(dat_all<=0);
        %         % Fit t1.
        %         xd = dat_all(t1); yd = diff_A1(t1);
        xdat = [ones(size(xd)) xd];
        [b,bint,r,rint,stats]=regress(yd,xdat);
        statVal = stats(3);
        x_plot = linspace(min(xd),max(xd),10);
        %         if statVal < 0.05
        %             plot(x_plot,b(1)+b(2)*x_plot,'g-','linewidth',1);
        %         else
        %             plot(x_plot,b(1)+b(2)*x_plot,'g--','linewidth',1);
        %         end
        plot([-0.5 0.5],[0 0],'k--'); plot([0 0],[-0.5 0.5],'k--');
        % text(-0.45,0.9,strcat(['slope t1:',num2str(b(2)),'p=',num2str(statVal)]),'fontsize',7,'color','k','fontname','arial');
        axis([-xM xM -yM yM]); axis square;
        slope_val(jj,2) = b(2); slope_stat(jj,2) = statVal;
        [r p] = corr(diff_A1,dat_all,'type','Spearman');
        all_dat_R(jj,2) = r; all_dat_P(jj,2) = p;
        
        % Fit data for unconditioned rsc > 0.
        xd = dat_all(dat_all>0); yd = diff_A1(dat_all>0);
        %         % Fit t2.
        %         xd = dat_all(t2); yd = diff_A1(t2);
        xdat = [ones(size(xd)) xd];
        [b,bint,r,rint,stats]=regress(yd,xdat);
        statVal = stats(3);
        x_plot = linspace(min(xd),max(xd),10);
        %         if statVal < 0.05
        %             plot(x_plot,b(1)+b(2)*x_plot,'g-','linewidth',1);
        %         else
        %             plot(x_plot,b(1)+b(2)*x_plot,'g--','linewidth',1);
        %         end
        plot([-0.5 0.6],[0 0],'k--'); plot([0 0],[-0.5 0.5],'k--');
        % text(-0.45,0.8,strcat(['slope t2:',num2str(b(2)),'p=',num2str(statVal)]),'fontsize',7,'color','k','fontname','arial');
        axis([-xM xM -yM yM]); axis square;
        slope_val(jj,3) = b(2); slope_stat(jj,3) = statVal;
        [r p] = corr(diff_A1,dat_all,'type','Spearman');
        all_dat_R(jj,3) = r; all_dat_P(jj,3) = p;
        
        %         % Fit data for unconditioned rsc > 0.
        %         % xd = dat_all(dat_all>0); yd = diff_A1(dat_all>0);
        %         % Fit t3.
        %         xd = dat_all(t3); yd = diff_A1(t3);
        %         mdl1 = fitlm(xd,yd);
        %         coeffs = mdl1.Coefficients.Estimate;
        %         pvals = mdl1.Coefficients.pValue;
        %         slope = coeffs(2); interc = coeffs(1);
        %         b = [interc slope];
        %         statVal = pvals(2);
        %         x_plot = linspace(min(dat_all(t3)),max(dat_all(t3)),10);
        %         if statVal < 0.05
        %             plot(x_plot,b(1)+b(2)*x_plot,'g-','linewidth',1);
        %         else
        %             plot(x_plot,b(1)+b(2)*x_plot,'g--','linewidth',1);
        %         end
        %         slope_val(1,jj) = b(2); slope_stat(1,jj) = statVal;
        %         plot([-0.5 0.5],[0 0],'k--'); plot([0 0],[-0.5 0.5],'k--');
        %         text(-0.45,1.4,strcat(['slope t3:',num2str(b(2)),'p=',num2str(statVal)]),'fontsize',7,'color','k','fontname','arial');
        %         %         axis([-xM xM -yM yM]);
        %         %         axis square;
        %         slope_val(jj,4) = b(2); slope_stat(jj,4) = statVal;
        %         [r p] = corr(diff_A1,dat_all,'type','Spearman');
        %         all_dat_R(jj,4) = r; all_dat_P(jj,4) = p;
        
        
        if jj == 3 xlabel('unconditioned ACC r_s_c'); end
        if jj == 1 ylabel({'LC conditioned ACC r_s_c';'ACC r_s_c (LC_n_o_n_z_e_r_o) -';'ACC r_s_c (LC_z_e_r_o)'}); end
        
        [r p] = corr(diff_A1,dat_all,'type','Spearman'); all_dat_R(jj) = r; all_dat_P(jj) = p;
        
        [r p] = corrcoef(xd(xd<0),yd(xd<0));
        rrN(jj) = r(1,2); ppN(jj) = p(1,2);
        
        [r p] = corrcoef(xd(xd>0),yd(xd>0));
        rrP(jj) = r(1,2); ppP(jj) = p(1,2);
        
    end    
    
    %%
    
    ax_ind = 0;
    
    for jj = 1:length(xbs)
        
        ax_ind = ax_ind + 1;
        kk = 5;
        % indd = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 16 17 18];
        indd = [17 2 3; 17 5 6; 17 8 9; 17 11 12; 17 16 18];
        
        % ***************************************
        
        dat_1 = all_mnk_dat{1}{indd(kk,1)}{jj};
        dat_Z1 = all_mnk_dat{1}{indd(kk,2)}{jj};
        dat_hi1 = all_mnk_dat{1}{indd(kk,3)}{jj};
        diff_A11 = dat_hi1-dat_Z1;
        xd1 = [ones(size(dat_1)) dat_1];
        yd1 = diff_A11; [b1,bint1,r,rint1,stats1]=regress(yd1,xd1);
        pv1(ax_ind) = stats1(3);
        
        dat_2 = all_mnk_dat{2}{indd(kk,1)}{jj};
        dat_Z2 = all_mnk_dat{2}{indd(kk,2)}{jj};
        dat_hi2 = all_mnk_dat{2}{indd(kk,3)}{jj};
        diff_A12 = dat_hi2-dat_Z2;
        xd2 = [ones(size(dat_2)) dat_2];
        yd2 = diff_A12; [b2,bint2,r,rint2,stats2]=regress(yd2,xd2);
        pv2(ax_ind) = stats2(3);
        
        % ***************************************
    end
    %     pause
    
    hi_slope = zeros(1,5);
    hi_slope_stat = zeros(1,5);
    allR = nans(3,5);
    allP = nans(3,5);
    %     indd = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 16 17 18];
    yM = .5; xM = .5;
    xb = linspace(-0.5,0.5,15); yb = linspace(-0.5,0.5,15);
    ax_ind = 0;
    % for jj = 2:2:10
    p_out = [];
    all_dat_R = [];
    all_dat_P = [];
    for jj = 1:length(xbs)
        b = []; bint = []; r = []; rint = []; stats = [];
        ax_ind = ax_ind + 1; kk = 5;
        axes(axs(ax_ind)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
        
        % ****************************
        %         dat_all = [all_mnk_dat{1}{17}{jj};all_mnk_dat{2}{17}{jj}];
        %         dat_Z = [all_mnk_dat{1}{16}{jj};all_mnk_dat{2}{16}{jj}];
        %         dat_hi = [all_mnk_dat{1}{18}{jj};all_mnk_dat{2}{18}{jj}];
        
        %         gi = find(isfinite(dat_all) & isfinite(dat_Z) & isfinite(dat_hi));
        %         diff_A1 = dat_hi(gi)-dat_Z(gi); dat_all = dat_all(gi);
        % ****************************
        
        % **************************** 043021
        dat_all = [datA_A{jj}{1}; datA_A{jj}{2}];
        diff_A1 = [datD_A{jj}{1}; datD_A{jj}{2}];
        
        gi = find(isfinite(dat_all) & isfinite(diff_A1));
        diff_A1 = diff_A1(gi);
        dat_all = dat_all(gi);
        % **************************** 043021
        
        terc_vals = quantile(dat_all,2);
        t1 = find(dat_all<=terc_vals(1)); t2 = find(dat_all>terc_vals(1) & dat_all<terc_vals(2)); t3 = find(dat_all>=terc_vals(2));
        [r1 p1] = corr(diff_A1(t1),dat_all(t1),'type','Spearman'); allR(1,jj) = r1; allP(1,jj) = p1;
        [r2 p2] = corr(diff_A1(t2),dat_all(t2),'type','Spearman'); allR(2,jj) = r2; allP(2,jj) = p2;
        [r3 p3] = corr(diff_A1(t3),dat_all(t3),'type','Spearman'); allR(3,jj) = r3; allP(3,jj) = p3;
        
        %         prc_val1 = prctile(abs(dat_all),99);
        %         prc_val2 = prctile(abs(diff_A1),99);
        %         prci = find(abs(dat_all)<prc_val & abs(diff_A1)<prc_val2);
        %         diff_A1 = diff_A1(prci); dat_all = dat_all(prci);
        
        plot(dat_all,diff_A1,'.','markersize',2);
        [b i] = sort(dat_all); mm1 = movmean(diff_A1(i),500);
        % plot(dat_all(i),mm1,'-','linewidth',2,'color',[0.5 0.5 0.5]);
        plot(dat_all(i),mm1,'-','linewidth',2,'color',[1 0 0]);
        
        %         xd = [ones(size(dat_all)) dat_all];
        %         yd = diff_A1;
        %         [b,bint,r,rint,stats]=regress(yd,xd);
        
        xd = dat_all;
        yd = diff_A1;
        mdl1 = fitlm(xd,yd);
        coeffs = mdl1.Coefficients.Estimate;
        pvals = mdl1.Coefficients.pValue;
        slope = coeffs(2); interc = coeffs(1);
        b = [interc slope];
        stats(3) = pvals(2);
        p_out = [p_out pvals(2)];
        
        x_plot = linspace(min(dat_all),max(dat_all),10);
        if stats(3) < 0.05
            plot(x_plot,b(1)+b(2)*x_plot,'k-','linewidth',1);
        else
            plot(x_plot,b(1)+b(2)*x_plot,'k--','linewidth',1);
        end
        hi_slope(ax_ind) = b(2); hi_slope_stat(ax_ind) = stats(3);
        plot([-0.5 0.5],[0 0],'k--'); plot([0 0],[-0.5 0.5],'k--');
        % text(-xM,yM,strcat([num2str(length(dat_all)),', ',num2str(b(2)),'bsiz=',num2str(round(xbs(jj)))]),'fontsize',7,'color','k','fontname','arial');
        text(-xM,yM,strcat([num2str(b(2)),'bsiz=',num2str(round(xbs(jj))),'p=',num2str(stats(3))]),'fontsize',7,'color','k','fontname','arial');
        axis([-xM xM -yM yM]);
        axis square;
        
        if ax_ind == 3 xlabel('unconditioned ACC r_s_c'); end
        if ax_ind == 1 ylabel({'LC conditioned ACC r_s_c';'ACC r_s_c (LC_n_o_n_z_e_r_o) -';'ACC r_s_c (LC_z_e_r_o)'}); end
        
        [r p] = corr(diff_A1,dat_all,'type','Spearman'); all_dat_R(jj) = r; all_dat_P(jj) = p;
        
        [r p] = corrcoef(xd(xd<0),yd(xd<0));
        rrN(jj) = r(1,2); ppN(jj) = p(1,2);
        
        [r p] = corrcoef(xd(xd>0),yd(xd>0));
        rrP(jj) = r(1,2); ppP(jj) = p(1,2);
        
    end
    
    %         pause
    
    %% TESTING Fig 5 (MAYBE): Relationship between unconditioned ACC rsc and ACC rsc change between LC zero and nonzero conditions.
    
    % 070420
    
    figureNumber = 55; num = 55; wid = 17.6; hts = [10]; cols = {5 5}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,1,1, [12], 'Joshi & Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    hi_slope = zeros(1,5);
    hi_slope_stat = zeros(1,5);
    
    %     indd = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 16 17 18];
    yM = .5;
    xM = .5;
    
    xb = linspace(-0.5,0.5,15);
    yb = linspace(-0.5,0.5,15);
    
    ax_ind = 0;
    stats_mnk = nans(2,5);
    stats_mnkP = nans(2,5);
    stats_mnkN = nans(2,5);
    stats_mnkPp = nans(2,5);
    stats_mnkNp = nans(2,5);
    
    for mm = 1:2
        % for jj = 2:2:10
        for jj = 1:length(xbs)
            
            ax_ind = ax_ind + 1;
            kk = 5;
            
            axes(axs(ax_ind)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
            dat_all = all_mnk_dat{mm}{17}{jj};
            dat_Z = all_mnk_dat{mm}{16}{jj};
            dat_hi = all_mnk_dat{mm}{18}{jj};
            
            gi = find(isfinite(dat_all) & isfinite(dat_Z) & isfinite(dat_hi));
            diff_A1 = dat_hi(gi)-dat_Z(gi);
            dat_all = dat_all(gi);
            
            plot(dat_all,diff_A1,'.','markersize',2);
            % histogram2(dat_all,diff_A1,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');
            % colormap(gray);
            
            [b i] = sort(dat_all); mm1 = movmean(diff_A1(i),1000);
            plot(dat_all(i),mm1,'-','linewidth',2,'color',[0.5 0.5 0.5]);
            xd = [ones(size(dat_all)) dat_all];
            yd = diff_A1; [b,bint,r,rint,stats]=regress(yd,xd);
            x_plot = linspace(min(dat_all),max(dat_all),10);
            if stats(3) < 0.05
                plot(x_plot,b(1)+b(2)*x_plot,'k-','linewidth',1);
            else
                plot(x_plot,b(1)+b(2)*x_plot,'k--','linewidth',1);
            end
            %             hi_slope(ax_ind) = b(2); hi_slope_stat(ax_ind) = stats(3);
            hi_slope(jj) = b(2); hi_slope_stat(jj) = stats(3);
            plot([-0.5 0.5],[0 0],'k--'); plot([0 0],[-0.5 0.5],'k--');
            [arrr peee] = corr(dat_all,diff_A1,'type','pearson');
            % text(-xM,yM,strcat([num2str(length(dat_all)),', ',num2str(b(2)),'bsiz=',num2str(round(xbs(jj)))]),'fontsize',7,'color','k','fontname','arial');
            % text(-xM,yM,strcat([num2str(b(2)),'bsiz=',num2str(round(xbs(jj))),'p=',num2str(stats(3)),'r=',num2str(arrr),'p=',num2str(peee)]),'fontsize',7,'color','k','fontname','arial');
            text(-xM,yM,strcat([num2str(b(2)),'p=',num2str(stats(3)),'r=',num2str(arrr)]),'fontsize',7,'color','k','fontname','arial');
            axis([-xM xM -yM yM]);
            axis square;
            
            if ax_ind == 3 xlabel('unconditioned ACC rsc'); end
            if ax_ind == 1 ylabel('\Delta ?ACC rsc re: unconditioned ACC rsc'); end
            
            stats
            
            [r p] = corrcoef(dat_all(dat_all<0),yd(dat_all<0));
            rrN(jj) = r(1,2); ppN(jj) = p(1,2);
            
            [r p] = corrcoef(dat_all(dat_all>0),yd(dat_all>0));
            rrP(jj) = r(1,2); ppP(jj) = p(1,2);
            
        end
        stats_mnk(mm,:) = hi_slope_stat;
        stats_mnkP(mm,:) = rrP;
        stats_mnkN(mm,:) = rrN;
        stats_mnkPp(mm,:) = ppP;
        stats_mnkNp(mm,:) = ppN;
    end
    
    %% SUPPLEMENT Fig 4: Relationship between unconditioned ACC rsc and ACC rsc change between LC zero and nonzero conditions.
    
    figureNumber = 41; num = 41; wid = 17.6; hts = [3]; cols = {5 5 5 5 5}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,1,2, [12], 'Joshi & Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    hi_slope = zeros(5,5);
    lo_slope = zeros(5,5);
    hi_slope_stat = zeros(5,5);
    lo_slope_stat = zeros(5,5);
    
    % indd = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 17 16 18];
    indd = [17 2 3; 17 5 6; 17 8 9; 17 11 12; 17 16 18];
    yM = .5;
    xM = .5;
    ax_ind = [0 5 10 15 20];
    
    xb = linspace(-0.5,0.5,15);
    yb = linspace(-0.5,0.5,15);
    
    for jj = 1:5 % Bin size.
        
        ax_ind = ax_ind+1;
        
        for kk = 1:5 % LC spike condition.
            
            axes(axs(ax_ind(kk))); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
            dat_all = [all_mnk_dat{1}{indd(kk,1)}{jj}; all_mnk_dat{2}{indd(kk,1)}{jj}];
            dat_Z = [all_mnk_dat{1}{indd(kk,2)}{jj}; all_mnk_dat{2}{indd(kk,2)}{jj}];
            dat_hi = [all_mnk_dat{1}{indd(kk,3)}{jj}; all_mnk_dat{2}{indd(kk,3)}{jj}];
            
            gi = find(isfinite(dat_all) & isfinite(dat_Z) & isfinite(dat_hi));
            diff_A1 = dat_hi(gi)-dat_Z(gi);
            dat_all = dat_all(gi);
            terc_vals = quantile(dat_all,2);
            t1 = find(dat_all<=terc_vals(1)); t2 = find(dat_all>terc_vals(1) & dat_all<terc_vals(2)); t3 = find(dat_all>=terc_vals(2));
            
            plot(dat_all,diff_A1,'.','markersize',2);
            % histogram2(dat_all,diff_A1,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');
            % colormap(gray);
            
            [b i] = sort(dat_all); mm1 = movmean(diff_A1(i),1000);
            plot(dat_all(i),mm1,'-','linewidth',2,'color',[0.5 0.5 0.5]);
            
            % Plot tercile 1.
            dat_all1 = dat_all; diff_A11 = diff_A1;
            
            %             xd = [ones(size(dat_all1)) dat_all1];
            %             yd = diff_A11;
            %             [b,bint,r,rint,stats]=regress(yd,xd);
            
            xd = dat_all;
            yd = diff_A1;
            mdl1 = fitlm(xd,yd);
            coeffs = mdl1.Coefficients.Estimate;
            pvals = mdl1.Coefficients.pValue;
            slope = coeffs(2); interc = coeffs(1);
            b = [interc slope];
            stats(3) = pvals(2);
            
            
            x_plot = linspace(min(dat_all1),max(dat_all1),10);
            if stats(3) < 0.05
                plot(x_plot,b(1)+b(2)*x_plot,'k-','linewidth',1);
            else
                plot(x_plot,b(1)+b(2)*x_plot,'k--','linewidth',1);
            end
            hi_slope(kk,jj) = b(2); hi_slope_stat(kk,jj) = stats(3);
            
            %             % Plot tercile 3
            %             dat_all3 = dat_all(t3); diff_A13 = diff_A1(t3);
            %             xd = [ones(size(dat_all3)) dat_all3]; yd = diff_A13; [b,bint,r,rint,stats]=regress(yd,xd);
            %             x_plot = linspace(min(dat_all3),max(dat_all3),10);
            %             if stats(3) < 0.05
            %                 plot(x_plot,b(1)+b(2)*x_plot,'k-','linewidth',1);
            %             else
            %                 plot(x_plot,b(1)+b(2)*x_plot,'k--','linewidth',1);
            %             end
            %             hi_slope(kk,jj) = b(2); hi_slope_stat(kk,jj) = stats(3);
            %
            plot([-0.5 0.5],[0 0],'k--'); plot([0 0],[-0.5 0.5],'k--');
            text(-xM,yM,strcat([num2str(b(2)),'bsiz=',num2str(round(xbs(jj))),'p=',num2str(stats(3))]),'fontsize',7,'color','k','fontname','arial');
            axis([-xM xM -yM yM]);
            
        end
    end
    
    %     figureNumber = 42; num = 42; wid = 17.6; hts = [3]; cols = {5 5 5 5 5}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,1,2, [12], 'Joshi & Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    %     movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    %
    %     ax_ind = [0 5 10 15 20];
    %
    %     for jj = 6:10
    %
    %         ax_ind = ax_ind+1;
    %
    %         for kk = 1:5
    %
    %             axes(axs(ax_ind(kk))); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %             dat_all = [all_mnk_dat{1}{indd(kk,1)}{jj}; all_mnk_dat{2}{indd(kk,1)}{jj}];
    %             dat_Z = [all_mnk_dat{1}{indd(kk,2)}{jj}; all_mnk_dat{2}{indd(kk,2)}{jj}];
    %             dat_hi = [all_mnk_dat{1}{indd(kk,3)}{jj}; all_mnk_dat{2}{indd(kk,3)}{jj}];
    %
    %             gi = find(isfinite(dat_all) & isfinite(dat_Z) & isfinite(dat_hi));
    %             diff_A1 = dat_hi(gi)-dat_Z(gi);
    %             dat_all = dat_all(gi);
    %             terc_vals = quantile(dat_all,2);
    %             t1 = find(dat_all<=terc_vals(1)); t2 = find(dat_all>terc_vals(1) & dat_all<terc_vals(2)); t3 = find(dat_all>=terc_vals(2));
    %
    %             plot(dat_all,diff_A1,'.','markersize',2);
    %             % histogram2(dat_all,diff_A1,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');
    %             % colormap(gray);
    %
    %             [b i] = sort(dat_all); mm1 = movmean(diff_A1(i),1000);
    %             plot(dat_all(i),mm1,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    %
    %             % Plot tercile 1.
    %             dat_all1 = dat_all(t1); diff_A11 = diff_A1(t1);
    %             xd = [ones(size(dat_all1)) dat_all1]; yd = diff_A11; [b,bint,r,rint,stats]=regress(yd,xd);
    %             x_plot = linspace(min(dat_all1),max(dat_all1),10);
    %             if stats(3) < 0.05
    %                 plot(x_plot,b(1)+b(2)*x_plot,'k-','linewidth',1);
    %             else
    %                 plot(x_plot,b(1)+b(2)*x_plot,'k--','linewidth',1);
    %             end
    %             hi_slope(kk,jj) = b(2); hi_slope_stat(kk,jj) = stats(3);
    %
    %             % Plot tercile 3
    %             dat_all3 = dat_all(t3); diff_A13 = diff_A1(t3);
    %             xd = [ones(size(dat_all3)) dat_all3]; yd = diff_A13; [b,bint,r,rint,stats]=regress(yd,xd);
    %             x_plot = linspace(min(dat_all3),max(dat_all3),10);
    %             if stats(3) < 0.05
    %                 plot(x_plot,b(1)+b(2)*x_plot,'k-','linewidth',1);
    %             else
    %                 plot(x_plot,b(1)+b(2)*x_plot,'k--','linewidth',1);
    %             end
    %             hi_slope(kk,jj) = b(2); hi_slope_stat(kk,jj) = stats(3);
    %
    %             plot([-0.6 0.6],[0 0],'k--'); plot([0 0],[-0.5 0.5],'k--');
    %             % text(-xM,yM,strcat([num2str(length(dat_all)),', ',num2str(b(2)),'bsiz=',num2str(round(xbs(jj)))]),'fontsize',7,'color','k','fontname','arial');
    %             text(-xM,yM,strcat([num2str(b(2)),'bsiz=',num2str(round(xbs(jj))),'p=',num2str(stats(3))]),'fontsize',7,'color','k','fontname','arial');
    %             axis([-xM xM -yM yM]);
    %
    %         end
    %     end
    
    %     pause
        
    %% Supplement: rsc and pair firing rate/pair Fano factor.
    
    % figureNumber = 223; num = 223; wid = 17.6; hts = [8]; cols = {2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,3,1, [12], 'Joshi & Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
    figureNumber = 45; num = 45; wid = 17.6; hts = [12]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,2,2, [12], 'Joshi & Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    %%
    
    nb = 21;
    yb = linspace(-.5,.5,nb);
    xb = linspace(-2,2,nb);
    % giii = isfinite(all_rate_diff) & isfinite(all_znz_diff);
    %     xDat = all_rate_diff(giii);
    %     yDat = all_znz_diff(giii);
    
    for mm = 1:2
        
        xDat = all_rate_diff{mm};
        yDat = all_znz_diff{mm};
        [r, p] = corr(xDat,yDat,'type','Spearman');
        
        r_out(mm) = r;
        p_out(mm) = p;
        
    end
    
    xDat = [all_rate_diff{1};all_rate_diff{2}];
    yDat = [all_znz_diff{1}; all_znz_diff{2}];
    
    %     nx = hist(xDat,xb); % mx = nanmedian(xDat);
    %     ny = hist(yDat,yb); % my = nanmedian(yDat);
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    % figure; hold on;
    
    histogram2(xDat,yDat,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgealpha',0,'edgecolor','none');
    %     mean_datS = imgaussfilt(mean_datS,sm_sigma);
    %     imagesc(tax,pax,mean_datS); colorbar('box','off');
    
    % histogram2(xDat,yDat,'DisplayStyle','tile','ShowEmptyBins','on','edgealpha',0,'edgecolor','none');
    colormap(gray);
    
    
    % plot(xDat,yDat,'k.','markersize',10,'color',[.5 .5 .5]);
    % plot(mF,mR,'r.','markersize',10);
    % plot([-1 1],[-1 1],'g:');
    plot([-2 2],[0 0],'w:');
    plot([0 0],[-.5 .5],'w:');
    % plot(xx,yy,'r-');
    title({'ACC: Relation between changes in';'pair firing rate and rsc'},'fontsize',10);
    xlabel('change in pair firing rate','fontsize',10,'color','k');
    ylabel('change in rsc','fontsize',10,'color','k');
    cbh = colorbar; set(cbh,'box','off');
    ylabel(cbh,'number of ACC neuron pairs','fontsize',10,'color','k');
    axis([-2 2 -.5 .5]);
    axis square;
    set(gca,'fontname','arial');
    
    %%
    
    nb = 21;
    yb = linspace(-.5,.5,nb);
    xb = linspace(-2,2,nb);
    
    %     giii = isfinite(all_FF_diff) & isfinite(all_znz_diff2);
    %     xDat = all_FF_diff(giii);
    %     yDat = all_znz_diff2(giii);
    
    for mm = 1:2
        
        xDat = all_FF_diff{mm};
        yDat = all_znz_diff{mm};
        [r, p] = corr(xDat,yDat,'type','Spearman');
        
        r_out2(mm) = r;
        p_out2(mm) = p;
        
        
    end
    
    xDat = [all_FF_diff{1};all_FF_diff{2}];
    yDat = [all_znz_diff2{1}; all_znz_diff2{2}];
    
    %     nx = hist(xDat,xb); % mx = nanmedian(xDat);
    %     ny = hist(yDat,yb); % my = nanmedian(yDat);
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    % figure; hold on;
    
    histogram2(xDat,yDat,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgealpha',0,'edgecolor','none');
    % histogram2(xDat,yDat,'DisplayStyle','tile','ShowEmptyBins','on','edgealpha',0,'edgecolor','none');
    colormap(gray);
    
    [r2, p2] = corr(xDat,yDat,'type','Spearman');
    
    % plot(xDat,yDat,'k.','markersize',10,'color',[.5 .5 .5]);
    % plot(mF,mR,'r.','markersize',10);
    % plot([-1 1],[-1 1],'g:');
    plot([-2 2],[0 0],'w:');
    plot([0 0],[-.5 .5],'w:');
    % plot(xx,yy,'r-');
    title({'ACC: Relation between changes in';'pair Fano factor and rsc'},'fontsize',10);
    xlabel('change in pair Fano factor','fontsize',10,'color','k');
    ylabel('change in rsc','fontsize',10,'color','k');
    cbh = colorbar; set(cbh,'box','off');
    ylabel(cbh,'number of ACC neuron pairs','fontsize',10,'color','k');
    %     axis([-2 2 -.5 .5]);
    axis square;
    set(gca,'fontname','arial');
    
    
    
    %% SUPPLEMENT:
    
    figureNumber = 46; num = 46; wid = 14; hts = [10]; cols = {1}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,1,2, [12], 'Joshi & Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    ms = [5 7 9 11 7];
    
    % Now Plot.
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    for kk = 1:length(xbs)
        
        sv1 = hi_slope(1,kk); sv1_s = hi_slope_stat(1,kk);
        if sv1_s<0.05 plot(xbs(kk)-20,sv1,'ko','markersize',ms(1),'linewidth',2); end
        if sv1_s>0.05 plot(xbs(kk)-20,sv1,'ko','markersize',ms(1),'linewidth',.5); end
        
        sv2 = hi_slope(2,kk); sv2_s = hi_slope_stat(2,kk);
        if sv2_s<0.05 plot(xbs(kk)-10,sv2,'ko','markersize',ms(2),'linewidth',2); end
        if sv2_s>0.05 plot(xbs(kk)-10,sv2,'ko','markersize',ms(2),'linewidth',.5); end
        
        sv3 = hi_slope(3,kk); sv3_s = hi_slope_stat(3,kk);
        if sv3_s<0.05 plot(xbs(kk)+10,sv3,'ko','markersize',ms(3),'linewidth',2); end
        if sv3_s>0.05 plot(xbs(kk)+10,sv3,'ko','markersize',ms(3),'linewidth',.5); end
        
        sv4 = hi_slope(4,kk); sv4_s = hi_slope_stat(4,kk);
        if sv4_s<0.05 plot(xbs(kk)+20,sv4,'ko','markersize',ms(4),'linewidth',2); end
        if sv4_s>0.05 plot(xbs(kk)+20,sv4,'ko','markersize',ms(4),'linewidth',.5); end
        
        sv5 = hi_slope(5,kk); sv5_s = hi_slope_stat(5,kk);
        if sv5_s<0.05 plot(xbs(kk)+20,sv5,'kd','markersize',ms(4),'linewidth',2); end
        if sv5_s>0.05 plot(xbs(kk)+20,sv5,'kd','markersize',ms(4),'linewidth',.5); end
        
        if kk == 1
            legend('LC=1','LC=2','LC=3','LC=4','LC>0','autoupdate','off','fontsize',10,'location','southwest','fontname','arial'); legend('boxoff');
        end
        
    end
    
    plot([50 1050],[0 0],'k--');
    axis([50 1050 -0.41 0.05]);
    xlabel('binsize (ms)');
    ylabel('slope (\DeltaACC rsc vs unconditioned ACC rsc)');
    
    %%
    
    figureNumber = 49; num = 49; wid = 14; hts = [5]; cols = {3 3}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,3,1, [12], 'Joshi & Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    %%
    
    clear; clear all;
    
    cd C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC;
    ACC_rsc = load('Fig4_LC_ACC_Compare'); ACC_dat = ACC_rsc.all_bdat; binDVals = ACC_rsc.binDiffVals;
    LC_rsc = load('Fig4_ACC_LC_Compare'); LC_dat = LC_rsc.all_bdat_LC; binDValsL = LC_rsc.binDiffValsL;
    
    nb = length(ACC_dat);
    pvals = nans(nb,3); % num bins x 3 terciles.
    hvals = nans(nb,3); % num bins x 3 terciles.
    
    mpVals = [];
    for mm = 1:2
        pvals = nans(nb,3); % num bins x 3 terciles.
        hvals = nans(nb,3); % num bins x 3 terciles.
        for jj = 1:length(ACC_dat)
            [h p] = kstest2(binDVals{jj}{1}{mm},binDValsL{jj}{1}{mm}); pvals(jj,1) = p; hvals(jj,1) = h;
            [h p] = kstest2(binDVals{jj}{2}{mm},binDValsL{jj}{2}{mm}); pvals(jj,2) = p; hvals(jj,2) = h;
            [h p] = kstest2(binDVals{jj}{3}{mm},binDValsL{jj}{3}{mm}); pvals(jj,3) = p; hvals(jj,3) = h;
        end
        mpVals{mm} = {pvals hvals};
    end
    
    for jj = 1:length(ACC_dat)
        
        [h p] = kstest2(ACC_dat{jj}{1},LC_dat{jj}{1}); pvals(jj,1) = p; hvals(jj,1) = h;
        [h p] = kstest2(ACC_dat{jj}{2},LC_dat{jj}{2}); pvals(jj,2) = p; hvals(jj,2) = h;
        [h p] = kstest2(ACC_dat{jj}{3},LC_dat{jj}{3}); pvals(jj,3) = p; hvals(jj,3) = h;
        %         pvals(jj,1) = ranksum(ACC_dat{jj}{1},LC_dat{jj}{1});
        %         pvals(jj,2) = ranksum(ACC_dat{jj}{2},LC_dat{jj}{2});
        %         pvals(jj,3) = ranksum(ACC_dat{jj}{3},LC_dat{jj}{3});
        
    end
    
    figure;
    ddatA = ACC_dat{5}{3}; subplot(2,1,1); hold on; hist(ddatA); plot([nanmedian(ddatA) nanmedian(ddatA)],[0 300],'r-'); % xlim([-0.5 0.5]);
    ddatL = LC_dat{5}{3}; subplot(2,1,2); hold on; hist(ddatL); plot([nanmedian(ddatL) nanmedian(ddatL)],[0 10],'r-'); % xlim([-0.5 0.5]);
    ranksum(ddatA,ddatL);
    signrank(ddatA)
    signrank(ddatL)
    
    for jj = 1:5
        aM(jj) = nanmedian(ACC_dat{jj}{3});
        aL(jj) = nanmedian(LC_dat{jj}{3});
    end
    
    % Rank-sum test:
    %     Bin     Tercile1   Tercile2  Tercile3
    %
    % 200     0.1611    0.6134    0.9841
    % 400     0.0273    0.2535    0.7932
    % 600     0.0237    0.5962    0.3955
    % 800     0.3427    0.9061    0.4502
    % 1000   0.2800    0.2888    0.3786
    
    % KS test:
    % hvals =
    %
    %      1     0     0
    %      1     0     0
    %      1     0     0
    %      1     0     0
    %      1     0     0
    %
    % pvals =
    %
    %     0.0023    0.2681    0.9186
    %     0.0001    0.1704    0.5931
    %     0.0084    0.4081    0.1925
    %     0.0339    0.1965    0.4007
    %     0.0321    0.2267    0.2064
    
    bini = 5;
    d1A = ACC_dat{bini}{1}; d1L = LC_dat{bini}{1};
    d2A = ACC_dat{bini}{2}; d2L = LC_dat{bini}{2};
    d3A = ACC_dat{bini}{3}; d3L = LC_dat{bini}{3};
    
    minVal = min([d1A;d2A;d3A;d1L;d2L;d3L]);
    maxVal = max([d1A;d2A;d3A;d1L;d2L;d3L]);
    minVal = -0.8;
    maxVal = 1;
    
    nBin=21; % # of bins between X_min and X_max
    dB=(maxVal-minVal)/nBin; % bin size
    BE=(minVal-dB):dB:(maxVal+dB); % bin edges; first and last edges are at X_min-dB and X_max+dB
    
    % Now Plot.
    
    % ACC.
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    H1A=histcounts(d1A,BE); % histogram of X
    MDF=H1A/length(d1A); % mass density function
    Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
    CMF=cumsum(MDF); % Cumulative mass function
    h=bar(Bc,CMF);
    set(h,'FaceColor',0.75*[1 1 1],'EdgeColor','none');
    axis([minVal-dB maxVal+dB 0 1.05]);
    xlabel('r_s_c difference');
    ylabel('Cumulative mass function');
    plot([0 0],[0 1.05],'k--','linewidth',0.5);
    text(-2,0.5,'ACC');
    plot([minVal-dB maxVal+dB],[0.5 0.5],'k--','linewidth',0.5);
    text(-0.25,1.2,'Tercile 1');
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    H2A=histcounts(d2A,BE); % histogram of X
    MDF=H2A/length(d2A); % mass density function
    Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
    CMF=cumsum(MDF); % Cumulative mass function
    h=bar(Bc,CMF);
    set(h,'FaceColor',0.75*[1 1 1],'EdgeColor','none');
    axis([minVal-dB maxVal+dB 0 1.05]);
    xlabel('r_s_c difference');
    %     ylabel('Cumulative mass function');
    plot([0 0],[0 1.05],'k--','linewidth',0.5);
    plot([minVal-dB maxVal+dB],[0.5 0.5],'k--','linewidth',0.5);
    text(-0.25,1.2,'Tercile 2');
    text(-0.5,1.4,'Bin size = 1000 ms');
    
    axes(axs(3)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    H3A=histcounts(d3A,BE); % histogram of X
    MDF=H3A/length(d3A); % mass density function
    Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
    CMF=cumsum(MDF); % Cumulative mass function
    h=bar(Bc,CMF);
    set(h,'FaceColor',0.75*[1 1 1],'EdgeColor','none');
    axis([minVal-dB maxVal+dB 0 1.05]);
    xlabel('r_s_c difference');
    %     ylabel('Cumulative mass function');
    plot([0 0],[0 1.05],'k--','linewidth',0.5);
    plot([minVal-dB maxVal+dB],[0.5 0.5],'k--','linewidth',0.5);
    text(-0.25,1.2,'Tercile 3');
    
    % LC.
    axes(axs(4)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    H1A=histcounts(d1L,BE); % histogram of X
    MDF=H1A/length(d1L); % mass density function
    Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
    CMF=cumsum(MDF); % Cumulative mass function
    h=bar(Bc,CMF);
    set(h,'FaceColor',0.75*[1 1 1],'EdgeColor','none');
    axis([minVal-dB maxVal+dB 0 1.05]);
    xlabel('r_s_c difference');
    ylabel('Cumulative mass function');
    plot([0 0],[0 1.05],'k--','linewidth',0.5);
    text(-2,0.5,'LC');
    plot([minVal-dB maxVal+dB],[0.5 0.5],'k--','linewidth',0.5);
    
    axes(axs(5)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    H2A=histcounts(d2L,BE); % histogram of X
    MDF=H2A/length(d2L); % mass density function
    Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
    CMF=cumsum(MDF); % Cumulative mass function
    h=bar(Bc,CMF);
    set(h,'FaceColor',0.75*[1 1 1],'EdgeColor','none');
    axis([minVal-dB maxVal+dB 0 1.05]);
    xlabel('r_s_c difference');
    %     ylabel('Cumulative mass function');
    plot([0 0],[0 1.05],'k--','linewidth',0.5);
    plot([minVal-dB maxVal+dB],[0.5 0.5],'k--','linewidth',0.5);
    
    axes(axs(6)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    H3A=histcounts(d3L,BE); % histogram of X
    MDF=H3A/length(d3L); % mass density function
    Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
    CMF=cumsum(MDF); % Cumulative mass function
    h=bar(Bc,CMF);
    set(h,'FaceColor',0.75*[1 1 1],'EdgeColor','none');
    axis([minVal-dB maxVal+dB 0 1.05]);
    xlabel('r_s_c difference');
    %     ylabel('Cumulative mass function');
    plot([0 0],[0 1.05],'k--','linewidth',0.5);
    plot([minVal-dB maxVal+dB],[0.5 0.5],'k--','linewidth',0.5);
    
    
end

