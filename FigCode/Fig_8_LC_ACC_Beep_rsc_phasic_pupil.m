% Fig_8_LC_ACC_Beep_rsc_phasic_pupil.m

% INFO: Calculate binned rsc for ACC pairs conditioned on phasic or no phasic LC responses to beep.
% 
% OLD: Also does calculation in tertiles depending on unconditioned ACC rsc (ie, rsc calculated for all trials, regardless of LC response).
%
% 022019: Sidd
% Added Fano and spike rate calculation for LC phasic/no phasic responses.
%
% Modified from: Fig8_LC_Dual_Recording_rsc.m
% Mod 122818 - Sidd. From: "LC_Dual_Recording_xcor.m".
% Mod: 080118 - Sidd.
% Mod: 092718: Added rCCG calculation as done for LC, IC and ACC pairs.
%
% Fixation expts ONLY.
%
% Calls "getLC_DualArea_Beep_binned_pupil_rsc3" to do the heavy lifting.

%% Setup stuff:

% Clear everything?
% clear; clear all;

% For cross correlograms, need to bin at 1msec resolution.
% backTime = 1500; fwdTime = 1500;
% tmin   = -backTime; tmax   = fwdTime; tsize  = 100; tstep  = 20;
% binsX  = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
% binsX = -backTime:fwdTime;
% tax    = mean(binsX,2);
% nb = size(binsX,1);

% *********************************

%% Calculate phasic/no-phasic rsc.

% Don't reanalyze data?
reanalyze = false;
% % Reanalyze data?
reanalyze = true;

if reanalyze
    clear; clear all;
    NTL=5;
    %         NTL=10;
    monks = {'Sprout','Cicero'}; % Add mnks as needed.
    sites = {'LC_ACC_Fixation'}; % Only LC for now.
    base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.
    nMonks = length(monks); nSites = length(sites);
    mnkRes = []; binSiz = 1;
    
    for mm = 1:nMonks % Loop over monkeys.
        
        sitRes = cell(1,nSites);
        
        for ss = 1:nSites  % Loop over brainstem sites.
            
            inDir= strcat([base_dir,'\',monks{mm},'\',sites{ss},'\clean']); % Create dir name for input (clean) files.
            cd(inDir);  dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
            fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
            
            if ~isempty(fnames)
                nf = length(fnames);
                xcor_dat = cell(1,nf);
                for ff = 1:nf
                    % for ff = 1:10
                    % ff = 3
                    load(fnames{ff});
                    % pause % For Testing.
                    % rsc_dat{ff} = getLC_DualArea_Beep_binned_pupil_rsc_BACUP(siteData,NTL,ff); % Pairwise spike count correlations.
                    
                    % rsc_dat{ff} = getLC_DualArea_Beep_binned_pupil_rsc4(siteData,NTL,ff); % Pairwise spike count correlations; grouped by LC spiking condition.
                    rsc_dat{ff} = getLC_DualArea_Beep_binned_pupil_rsc3(siteData,NTL,ff); % Pairwise spike count correlations; grouped by pupil size.
                    
                    disp(sprintf('File %d%sof%s%d', ff,' ',' ',nf));
                end
            end
            sitRes{ss} = rsc_dat;
            clear rsc_dat;
        end % Sites.
        mnkRes{mm} = sitRes;
        
        % Save analysis:
        savedir = strcat([base_dir,'\Results\Results_2021\LC_ACC']);
        cd(savedir);
        clear sitRes;
    end % Mnks.
    
    % save LC_ACC_Results_beep_rsc_100120 mnkRes; % With new way of dividing trials into zero, mid, high.
    % save LC_ACC_Results_beep_rsc_031721 mnkRes; % Pupil Slope.
    % save LC_ACC_Results_beep_pupil_rsc_031721b mnkRes; %
    % save LC_ACC_Results_beep_pupil_rsc_043021 mnkRes; %
    % save LC_ACC_Results_beep_pupil_rsc_050521b mnkRes; % With baseline and evoked pupil. [25 50 75].
    % save LC_ACC_Results_beep_pupil_rsc_050721 mnkRes; % With baseline and evoked pupil. Median split.
    % save LC_ACC_Results_beep_pupil_rsc_050721b mnkRes; % With baseline and evoked pupil. [25 50 75].
    % save LC_ACC_Results_beep_pupil_rsc_050721c mnkRes; % With baseline and evoked pupil. [40 50 60].
    % save LC_ACC_Results_beep_pupil_rsc_051021 mnkRes; % Median.
    % save LC_ACC_Results_beep_pupil_rsc_051021b mnkRes; % Quartile.
    %     save LC_ACC_Results_beep_pupil_rsc_051521 mnkRes; % With baseline and evoked pupil. Median split. 50ms baseline.
    %     save LC_ACC_Results_beep_pupil_rsc_051521 mnkRes; % With baseline and evoked pupil. Median split. 50ms baseline.
    % save LC_ACC_Results_beep_pupil_rsc_052121 mnkRes; % With baseline and evoked pupil. Median split. 50ms baseline.
    % save LC_ACC_Results_beep_pupil_rsc_100321 mnkRes; % With baseline and evoked pupil. Median split. 50ms baseline.
    save LC_ACC_Results_beep_pupil_rsc_112421 mnkRes; % 11/24/21: after re-submission. Now also saved baseline pupil for partial correlation.
    
    clear mnkRes;
end % If reanalyze.


%%

clear; clear all;
base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.
savedir = strcat([base_dir,'\Results\Results_2021\LC_ACC\']); %  Point to the current results directory.
cd(savedir);

monks = {'Sprout','Cicero'}; % Add monks as needed.
sites = {'LC_ACC_Fixation'}; % Add sites as needed.
nmonks = length(monks);
mnkFigDat = cell(1,2);

%% Organize rsc plot data for beep trials.

organizeYesno = 0;
organizeYesno = 1;

bs = linspace(200,1000,5); % 043021.
nbs = length(bs);

if organizeYesno
    
    %     saveFileName = 'LC_ACC_Results_beep_pupil_rsc_051021b';
    %     saveFileName = 'LC_ACC_Results_beep_pupil_rsc_051521';
    % saveFileName = 'LC_ACC_Results_beep_pupil_rsc_052121';
    % saveFileName = 'LC_ACC_Results_beep_pupil_rsc_100321';
    saveFileName = 'LC_ACC_Results_beep_pupil_rsc_112421'; % With baseline pupil also saved.
    
    load(saveFileName);
    
    for mm = 1:nmonks
        
        unit_index = 0; incindex = 0;
        ns = length(mnkRes{mm}{1}); % Number of sessions for this mnky.
        yyi = 0;
        all_loBpp = []; all_loApp = []; all_hiBpp = []; all_hiApp = [];
        all_loBpp3 = []; all_loApp3 = []; all_hiBpp3 = []; all_hiApp3 = [];
        all_loBppB = []; all_loAppB = []; all_hiBppB = []; all_hiAppB = [];
        all_loBpp3B = []; all_loApp3B = []; all_hiBpp3B = []; all_hiApp3B = [];
        si = 0; si3 = 0;
        sib = 0; sib3 = 0;
        ph_resp_pupil1 = [];
        ph_resp_pupil = [];
        al_pp_lo_p = []; al_pp_hi_p = []; al_pp_lo_np = []; al_pp_hi_np = [];
        % al_pp_lo_pM = []; al_pp_hi_pM = []; al_pp_lo_npM = []; al_pp_hi_npM = [];
        al_pp_lo_pB = []; al_pp_hi_pB = []; al_pp_lo_npB = []; al_pp_hi_npB = [];
        
        both_phasic = []; both_no_phasic = []; only_one_phasic = []; any_one_phasic = [];
        
        % pause
        
        for ii = 1:ns
            nu = length(mnkRes{mm}{1}{ii}{1}); % Get site data (one cell for each LC unit).
            pair_phasic = [];
            pair_no_phasic = [];
            for jj = 1:nu
                
                unit_index = unit_index+1;
                
                % pause
                rpp1 = corrcoef(mnkRes{mm}{1}{ii}{5}{jj}', mnkRes{mm}{1}{ii}{10}');
                rpp = partialcorr(mnkRes{mm}{1}{ii}{5}{jj}', mnkRes{mm}{1}{ii}{10}',mnkRes{mm}{1}{ii}{11}','type','Spearman');
                ph_resp_pupil1 = [ph_resp_pupil1 rpp1(1,2)];
                ph_resp_pupil = [ph_resp_pupil rpp];
                
                if nu > 1
                    pair_phasic = [pair_phasic; mnkRes{mm}{1}{ii}{2}{jj}];
                    pair_no_phasic = [pair_no_phasic; mnkRes{mm}{1}{ii}{3}{jj}];
                end
                
                bDat = mnkRes{mm}{1}{ii}{1}{jj}; % pause
                bDat3 = mnkRes{mm}{1}{ii}{4}{jj};
                
                bDatB = mnkRes{mm}{1}{ii}{8}{jj}; % pause
                bDat3B = mnkRes{mm}{1}{ii}{9}{jj};
                
                if ~isempty(bDat) & ~isempty(bDat3)
                    
                    al_pp_lo_p = [al_pp_lo_p; bDat{nbs}{5}];
                    al_pp_hi_p = [al_pp_hi_p; bDat{nbs}{6}];
                    %                     al_pp_lo_pM = [al_pp_lo_pM nanmedian(bDat{nbs}{5})];
                    %                     al_pp_hi_pM = [al_pp_hi_pM nanmedian(bDat{nbs}{6})];
                    
                    al_pp_lo_np = [al_pp_lo_np; bDat3{nbs}{5}(:)];
                    al_pp_hi_np = [al_pp_hi_np; bDat3{nbs}{6}(:)];
                    %                     al_pp_lo_npM = [al_pp_lo_npM nanmedian(bDat3{nbs}{5})];
                    %                     al_pp_hi_npM= [al_pp_hi_npM nanmedian(bDat3{nbs}{6})];
                end
                
                if ~isempty(bDatB) & ~isempty(bDat3B)
                    
                    al_pp_lo_pB = [al_pp_lo_p; bDatB{nbs}{5}];
                    al_pp_hi_pB = [al_pp_hi_p; bDatB{nbs}{6}];
                    %                     al_pp_lo_pM = [al_pp_lo_pM nanmedian(bDat{nbs}{5})];
                    %                     al_pp_hi_pM = [al_pp_hi_pM nanmedian(bDat{nbs}{6})];
                    
                    al_pp_lo_npB = [al_pp_lo_np; bDat3B{nbs}{5}(:)];
                    al_pp_hi_npB = [al_pp_hi_np; bDat3B{nbs}{6}(:)];
                    %                     al_pp_lo_npM = [al_pp_lo_npM nanmedian(bDat3{nbs}{5})];
                    %                     al_pp_hi_npM= [al_pp_hi_npM nanmedian(bDat3{nbs}{6})];
                end
                
                % Evoked pupil.
                if ~isempty(bDat) & ~isempty(bDat3)
                    nbins = length(bDat);
                    for kk = 1:nbins
                        ei = length(bDat{kk}{1}(:,1));
                        
                        % Pupil linked rsc.
                        all_loBpp(1+si:si+ei,kk) = bDat{kk}{1}(:,1);
                        all_loApp(1+si:si+ei,kk) = bDat{kk}{2}(:,1);
                        all_hiBpp(1+si:si+ei,kk) = bDat{kk}{3}(:,1);
                        all_hiApp(1+si:si+ei,kk) = bDat{kk}{4}(:,1);
                        
                    end % Bins.
                    
                    si = si+ei;
                    
                    nbins = length(bDat3);
                    for kk = 1:nbins
                        ei = length(bDat3{kk}{1}(:,1));
                        
                        % Pupil linked rsc.
                        all_loBpp3(1+si3:si3+ei,kk) = bDat3{kk}{1}(:,1);
                        all_loApp3(1+si3:si3+ei,kk) = bDat3{kk}{2}(:,1);
                        all_hiBpp3(1+si3:si3+ei,kk) = bDat3{kk}{3}(:,1);
                        all_hiApp3(1+si3:si3+ei,kk) = bDat3{kk}{4}(:,1);
                        
                    end % Bins.
                    
                    si3 = si3+ei;
                end % Check for empty data struct.
                
                % Baseline pupil.
                if ~isempty(bDatB) & ~isempty(bDat3B)
                    nbins = length(bDatB);
                    for kk = 1:nbins
                        ei = length(bDatB{kk}{1}(:,1));
                        
                        % Pupil linked rsc.
                        all_loBppB(1+sib:sib+ei,kk) = bDatB{kk}{1}(:,1);
                        all_loAppB(1+sib:sib+ei,kk) = bDatB{kk}{2}(:,1);
                        all_hiBppB(1+sib:sib+ei,kk) = bDatB{kk}{3}(:,1);
                        all_hiAppB(1+sib:sib+ei,kk) = bDatB{kk}{4}(:,1);
                        
                    end % Bins.
                    
                    sib = sib+ei;
                    
                    nbins = length(bDat3B);
                    for kk = 1:nbins
                        ei = length(bDat3B{kk}{1}(:,1));
                        
                        % Pupil linked rsc.
                        all_loBpp3B(1+sib3:sib3+ei,kk) = bDat3B{kk}{1}(:,1);
                        all_loApp3B(1+sib3:sib3+ei,kk) = bDat3B{kk}{2}(:,1);
                        all_hiBpp3B(1+sib3:sib3+ei,kk) = bDat3B{kk}{3}(:,1);
                        all_hiApp3B(1+sib3:sib3+ei,kk) = bDat3B{kk}{4}(:,1);
                        
                    end % Bins.
                    
                    sib3 = sib3+ei;
                end % Check for empty data struct.
                
                
            end % Units for this session.
            
            % pause
            if nu>1
                % pause
                nn = size(pair_phasic,1);
                pind = nchoosek(1:nn,2);
                
                for ppi = 1:size(pind,1)
                    both_phasic = [both_phasic (1/2)*nansum(pair_phasic(pind(ppi,1),:)+pair_phasic(pind(ppi,2),:))/size(pair_phasic,2)];
                    both_no_phasic = [both_no_phasic (1/2)*nansum((pair_no_phasic(pind(ppi,1),:)+pair_no_phasic(pind(ppi,2),:)))/size(pair_no_phasic,2)];
                    
                    n1 = pair_phasic(pind(ppi,1),:); n1_p = (n1==1);
                    n2 = pair_phasic(pind(ppi,2),:); n2_p = (n2==1);
                    
                    spp = n1_p+n2_p;
                    
                    only_one_phasic = [only_one_phasic sum(abs(n1_p-n2_p))/size(pair_phasic,2)];
                    any_one_phasic = [any_one_phasic sum(spp>0)/size(pair_phasic,2)];
                    
                end
            end
            
        end % Sessions for this monkey.
        
        mnkFigDat{mm} = {
            all_loBpp all_loApp all_hiBpp all_hiApp ... % 4/47
            all_loBpp3 all_loApp3 all_hiBpp3 all_hiApp3 ... % 4/51
            all_loBppB all_loAppB all_hiBppB all_hiAppB ... % 4/47
            all_loBpp3B all_loApp3B all_hiBpp3B all_hiApp3B ... % 4/51
            ph_resp_pupil ... % 1/52
            al_pp_lo_p al_pp_hi_p al_pp_lo_np al_pp_hi_np ... % 4/56
            both_phasic both_no_phasic only_one_phasic any_one_phasic ... % 4/60
            al_pp_lo_pB al_pp_hi_pB al_pp_lo_npB al_pp_hi_npB ... % 4/56
            ph_resp_pupil1
            };
        
        %         mnkFigDat{mm} = {
        %             all_loBnp all_loAnp all_loBp all_loAp all_file all_lci all_acpi ... % 7
        %             all_loBnp2 all_loAnp2 all_loBp2 all_loAp2 all_file2 all_lci2 all_acpi2 ... % 7/14
        %             all_loBnp3 all_loAnp3 all_loBp3 all_loAp3 all_file3 all_lci3 all_acpi3 ... % 7/21
        %             all_loBnp4 all_loAnp4 all_loBp4 all_loAp4 all_file4 all_lci4 all_acpi4 ... % 7/28
        %             pnp_dat pnp_dat2 pnp_dat3 pnp_dat4 ... % 4/32
        %             sp_dat all_mean_B all_mean_A all_var_B all_var_A all_mean_B3 all_mean_A3 all_var_B3 all_var_A3 phasicUnit phasicUnit2 ... % 11/43
        %             all_loBpp all_loApp all_hiBpp all_hiApp ... 4/47
        %             all_loBpp3 all_loApp3 all_hiBpp3 all_hiApp3 ... 4/51
        %             ph_resp_pupil ... 1/52
        %             al_pp_lo_p al_pp_hi_p al_pp_lo_np al_pp_hi_np ...
        %             both_phasic both_no_phasic only_one_phasic any_one_phasic ...
        %             al_pp_lo_pM al_pp_hi_pM al_pp_lo_npM al_pp_hi_npM
        %             };
    end % Monkeys.
    
    clear mnkRes; % Free up memory by only retaining the plot data in "mnkFigDat".
    
end

% pause

%% Plot stuff:

% Don't plot results?
plotYesNo = false;
% % Plot results?
plotYesNo = true;

if plotYesNo
    %%
    %     bs = linspace(100,1000,10); % 021720.
    xAx = bs;
    %     limVal=0;
    %     limVal2=99;
    %     limVal3 = 1e10; % Use as absolute criterion because sometimes the % difference blows up due to a infinitesimal denominator.
    %     diff_p11_lo = [];    diff_p12_lo = [];    diff_p31_lo = [];    diff_p32_lo = [];
    pval_both_np = nans(2,5);
    pval_both_npB = nans(2,5);
    pval1p = nans(2,5);
    pval2p = nans(2,5);
    pval1pB = nans(2,5);
    pval2pB = nans(2,5);
    % pval_both_np = nans(6,5);
    
    % pause
    
    %% ANOVA for baseline pupil.
    
    % SMALL EVOKED PUPIL, Phasic response.
    all_Bp1_lo = mnkFigDat{1}{9}; all_Bp2_lo =mnkFigDat{2}{9}; all_Ap1_lo = mnkFigDat{1}{10}; all_Ap2_lo = mnkFigDat{2}{10};
    % LARGE EVOKED PUPIL, Phasic response.
    all_Bp1_hi = mnkFigDat{1}{11}; all_Bp2_hi =mnkFigDat{2}{11}; all_Ap1_hi = mnkFigDat{1}{12}; all_Ap2_hi = mnkFigDat{2}{12};
    % SMALL EVOKED PUPIL, No phasic response.
    all_Bnp1_lo = mnkFigDat{1}{13}; all_Bnp2_lo =mnkFigDat{2}{13}; all_Anp1_lo = mnkFigDat{1}{14}; all_Anp2_lo = mnkFigDat{2}{14};
    % LARGE EVOKED PUPIL, No phasic response.
    all_Bnp1_hi = mnkFigDat{1}{15}; all_Bnp2_hi =mnkFigDat{2}{15}; all_Anp1_hi = mnkFigDat{1}{16}; all_Anp2_hi = mnkFigDat{2}{16};
    
    nbins = size(all_Anp2_hi,2); all_rsc_datB = []; main_groupB = []; bin_indexB = [];
    
    for kk = 1:nbins
        
        drsc1 = [all_Ap1_lo(:,kk) - all_Bp1_lo(:,kk); all_Ap2_lo(:,kk) - all_Bp2_lo(:,kk)];
        drsc2 = [all_Ap1_hi(:,kk) - all_Bp1_hi(:,kk); all_Ap2_hi(:,kk) - all_Bp2_hi(:,kk)];
        drsc3 = [all_Anp1_lo(:,kk) - all_Bnp1_lo(:,kk); all_Anp2_lo(:,kk) - all_Bnp2_lo(:,kk)];
        drsc4 = [all_Anp1_hi(:,kk) - all_Bnp1_hi(:,kk); all_Anp2_hi(:,kk) - all_Bnp2_hi(:,kk)];
        
        DDD1 = drsc2-drsc3;
        DDD2 = drsc1-drsc4;
        
        all_rsc_datB = [all_rsc_datB; DDD1; DDD2];
        main_groupB = [main_groupB; zeros(length(DDD1),1); ones(length(DDD2),1)];
        bin_indexB = [bin_indexB; kk*[ones(length(DDD1),1); ones(length(DDD2),1)]];
        
    end
    
    gi = find(isfinite(all_rsc_datB));
    all_rsc_datB = all_rsc_datB(gi);
    main_groupB = main_groupB(gi);
    bin_indexB = bin_indexB(gi);
    
    % ANOVA with display OFF.
    [pA,tbl,stats] = anovan(all_rsc_datB, {main_groupB, bin_indexB},'model','interaction','varnames',{'mainGroup','binIndex'});
    
    % ANOVA results for baseline pupil analysis. Quartiles.
    %     Source                                SumSq   d.f.  MeanSq      F          p
    %     pupil                                    1.545      1   1.54539     6.86   0.0089
    %     LCPhasic                              3.543      1   3.5429     15.73   0.0001
    %     binIndex                              1.457      4   0.36436     1.62   0.1675
    %     baselinePupil*LCPhasic         0.063      1   0.06329     0.28   0.5962
    %     baselinePupil*binIndex         1.614      4   0.40343     1.79   0.1283
    %     LCPhasic*binIndex               2.053      4   0.51319     2.28   0.059
    
    %% ANOVA for evoked pupil.
    
    % SMALL EVOKED PUPIL, Phasic response.
    all_Bp1_lo = mnkFigDat{1}{1}; all_Bp2_lo =mnkFigDat{2}{1}; all_Ap1_lo = mnkFigDat{1}{2}; all_Ap2_lo = mnkFigDat{2}{2};
    % LARGE EVOKED PUPIL, Phasic response.
    all_Bp1_hi = mnkFigDat{1}{3}; all_Bp2_hi =mnkFigDat{2}{3}; all_Ap1_hi = mnkFigDat{1}{4}; all_Ap2_hi = mnkFigDat{2}{4};
    % SMALL EVOKED PUPIL, No phasic response.
    all_Bnp1_lo = mnkFigDat{1}{5}; all_Bnp2_lo =mnkFigDat{2}{5}; all_Anp1_lo = mnkFigDat{1}{6}; all_Anp2_lo = mnkFigDat{2}{6};
    % LARGE EVOKED PUPIL, No phasic response.
    all_Bnp1_hi = mnkFigDat{1}{7}; all_Bnp2_hi =mnkFigDat{2}{7}; all_Anp1_hi = mnkFigDat{1}{8}; all_Anp2_hi = mnkFigDat{2}{8};
    
    nbins = size(all_Anp2_hi,2); all_rsc_dat = []; main_group = []; bin_index = [];
    
    for kk = 1:nbins
        
        drsc1 = [all_Ap1_lo(:,kk) - all_Bp1_lo(:,kk); all_Ap2_lo(:,kk) - all_Bp2_lo(:,kk)];
        drsc2 = [all_Ap1_hi(:,kk) - all_Bp1_hi(:,kk); all_Ap2_hi(:,kk) - all_Bp2_hi(:,kk)];
        drsc3 = [all_Anp1_lo(:,kk) - all_Bnp1_lo(:,kk); all_Anp2_lo(:,kk) - all_Bnp2_lo(:,kk)];
        drsc4 = [all_Anp1_hi(:,kk) - all_Bnp1_hi(:,kk); all_Anp2_hi(:,kk) - all_Bnp2_hi(:,kk)];
        
        DDD1 = drsc2-drsc3;
        DDD2 = drsc1-drsc4;
        
        all_rsc_dat = [all_rsc_dat; DDD1; DDD2];
        main_group = [main_group; zeros(length(DDD1),1); ones(length(DDD2),1)];
        bin_index = [bin_index; kk*[ones(length(DDD1),1); ones(length(DDD2),1)]];
        
    end
    
    gi = find(isfinite(all_rsc_dat));
    all_rsc_dat = all_rsc_dat(gi);
    main_group = main_group(gi);
    bin_index = bin_index(gi);
    
    % ANOVA with display OFF.
    [pA,tbl,stats] = anovan(all_rsc_dat, {main_group, bin_index},'model','interaction','varnames',{'mainGroup','binIndex'});
    
    % ANOVA results for evoked pupil analysis. Quartiles.
    %     Source                                SumSq   d.f.  MeanSq      F          p
    %     pupil                                    0.003      1   0.00298    0.01   0.9112
    %     LCPhasic                              1.649      1   1.64906    6.87   0.0089
    %     binIndex                              0.818      4   0.20453    0.85   0.4921
    %     evokedPupil*LCPhasic           0.596      1   0.59626    2.48   0.1152
    %     evokedPpupil*binIndex         0.152      4   0.03789    0.16   0.9594
    %     LCPhasic*binIndex                3.25       4   0.81256    3.39   0.0091
    
    %% ANOVA for combined baseline and evoked pupil.
    
    bp_index = [zeros(length(all_rsc_datB),1); ones(length(all_rsc_dat),1)];
    [pA,tbl,stats] = anovan([all_rsc_datB; all_rsc_dat], {[main_groupB; main_group], [bin_indexB; bin_index],bp_index},'model','interaction','varnames',{'mainGroup','binIndex','b_or_p'});
    %     [pA,tbl,stats] = anovan([all_rsc_datB; all_rsc_dat], {[pupil_indexB; pupil_index+2], [LCphasic_indexB; LCphasic_index+2], [bin_indexB; bin_index]},'model','interaction','varnames',{'sizePupil','LCPhasic','binIndex'});
    
    %   sizePupil                                   0.848      1   0.84783     3.65   0.0562
    %   LCPhasic                                   5.038      1   5.03753     21.69     0
    %   binIndex                                    1.96       4   0.48994     2.11   0.0771
    %   baseline_or_evoked                   0.405      1   0.40475     1.74   0.1869
    %   sizePupil*LCPhasic                     0.132      1   0.13229     0.57   0.4505
    %   sizePupil*binIndex                     1.294      4   0.32339     1.39   0.234
    %   sizePupil*baseline_or_evoked     0.595      1   0.59544     2.56   0.1095
    %   LCPhasic*binIndex                     5.099      4   1.2748       5.49   0.0002
    %   LCPhasic*baseline_or_evoked     0.218      1   0.21828     0.94   0.3324
    %   binIndex*baseline_or_evoked     0.248      4   0.06198     0.27   0.8994
    
    % pause
    
    %% Setup figure.
    
    figureNumber = 9; num = 9; wid = 17.6; hts = [10]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 2, 2, [12], 'Joshi&Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    % figureNumber = 9; num = 9; wid = 17.6; hts = [7]; cols = {6 6}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 2, 1, [12], 'Joshi&Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    %%
    
    % SMALL EVOKED PUPIL, Phasic response.
    all_Bp1_lo = mnkFigDat{1}{1}; all_Bp2_lo =mnkFigDat{2}{1}; all_Ap1_lo = mnkFigDat{1}{2}; all_Ap2_lo = mnkFigDat{2}{2};
    % LARGE EVOKED PUPIL, Phasic response.
    all_Bp1_hi = mnkFigDat{1}{3}; all_Bp2_hi =mnkFigDat{2}{3}; all_Ap1_hi = mnkFigDat{1}{4}; all_Ap2_hi = mnkFigDat{2}{4};
    % SMALL EVOKED PUPIL, No phasic response.
    all_Bnp1_lo = mnkFigDat{1}{5}; all_Bnp2_lo =mnkFigDat{2}{5}; all_Anp1_lo = mnkFigDat{1}{6}; all_Anp2_lo = mnkFigDat{2}{6};
    % LARGE EVOKED PUPIL, No phasic response.
    all_Bnp1_hi = mnkFigDat{1}{7}; all_Bnp2_hi =mnkFigDat{2}{7}; all_Anp1_hi = mnkFigDat{1}{8}; all_Anp2_hi = mnkFigDat{2}{8};
    
    %% Calculate plot things.
    
    % 3. Phasic response, large pupil - no-phasic response, small pupil.
    
    for kk = 1:length(bs)
        D1a = all_Anp1_lo(:,kk)-all_Bnp1_lo(:,kk); D1b = all_Anp2_lo(:,kk)-all_Bnp2_lo(:,kk);
        D1c = all_Ap1_hi(:,kk)-all_Bp1_hi(:,kk); D1d = all_Ap2_hi(:,kk)-all_Bp2_hi(:,kk);
        DD1 = D1c-D1a; DD1 = DD1(isfinite(DD1));
        DD2 = D1d-D1b; DD2 = DD2(isfinite(DD2));
        % pval1p(1,kk) = signrank(DD1); pval2p(1,kk) = signrank(DD2);
        pval_both_np(1,kk) = signrank([DD1; DD2]);
        pval_both_npRS(1,kk) = ranksum([D1c; D1d],[D1a; D1b]);
        diff_p5{kk} = DD1; diff_p6{kk} = DD2;
        
        
        % 4. No phasic response, large pupil - Phasic response, small pupil.
        
        D1aw = all_Ap1_lo(:,kk)-all_Bp1_lo(:,kk); D1bw = all_Ap2_lo(:,kk)-all_Bp2_lo(:,kk);
        D1cw = all_Anp1_hi(:,kk)-all_Bnp1_hi(:,kk); D1dw = all_Anp2_hi(:,kk)-all_Bnp2_hi(:,kk);
        DD1w = D1aw-D1cw; DD1w = DD1w(isfinite(DD1w));
        DD2w = D1bw-D1dw; DD2w = DD2w(isfinite(DD2w));
        % pval1p(2,kk) = signrank(DD1); pval2p(2,kk) = signrank(DD2);
        pval_both_np(2,kk) = signrank([DD1w; DD2w]);
        pval_both_npRS(2,kk) = ranksum([D1cw; D1dw],[D1aw; D1bw]);
        diff_p7{kk} = DD1w; diff_p8{kk} = DD2w;
        
        p_across_conds(kk) = ranksum([DD1; DD2],[DD1w; DD2w]);
        
    end
    
    % pause
    
    %% Plot.
    
    xjit = [-10 10];
    allD = {diff_p5 diff_p6 diff_p7 diff_p8};
    title_text = {'Evoked Pupil'};
    colors = {[0 0 0] [0.5 0.5 0.5]};
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    %     allPvals_befAft = [pvalphi;pvalnplo];
    
    for pind = 1:2
        for kk = 1:length(bs)
            
            ii = 2*pind;
            yd = [allD{ii-1}{kk}; allD{ii}{kk}]; mean_yd = nanmedian(yd);
            
            if pval_both_np(pind,kk) < 0.05
                plot(xjit(pind)+xAx(kk),mean_yd,'ko','markerfacecolor',colors{pind},'markeredgecolor',colors{pind},'markersize',8);
            end
            if pval_both_np(pind,kk) >= 0.05
                plot(xjit(pind)+xAx(kk),mean_yd,'ko','markerfacecolor','none','markeredgecolor',colors{pind},'markersize',8);
            end
            
            [ciV ciV_m] = bootci(1000,{@median,yd},'alpha',0.05,'type','norm'); % mean_ze = nanmean(ci_ze_m);
            plot([xjit(pind)+xAx(kk) xjit(pind)+xAx(kk)],[ciV(1) ciV(2)],'-','color',colors{pind}); % 95% CI, bootstrapped.
            
        end
    end
    
    for kk = 1:length(bs)
        if p_across_conds(kk)<0.05 plot(xjit(pind)+xAx(kk),0.07,'*','color',colors{pind}); end
    end
    
    plot([xAx(1) xAx(end)],[0 0],'k--');
    xlim([100 1100]); ylim([-0.03 0.08]);
    % ylim([-0.3 0.3]);
    set(gca,'fontname','arial','fontsize',10);
    title(title_text{1},'fontsize',12,'fontname','arial','fontweight','normal');
    ylabel({'ACC r_s_c';'After beep - before beep'},'fontsize',10,'fontname','arial');
    xlabel('bin size (ms)','fontsize',10,'fontname','arial');
    
    % *****************************************************************
    % pause
    
    %% Baseline.
    
    % SMALL EVOKED PUPIL, Phasic response.
    all_Bp1_lo = mnkFigDat{1}{9}; all_Bp2_lo =mnkFigDat{2}{9}; all_Ap1_lo = mnkFigDat{1}{10}; all_Ap2_lo = mnkFigDat{2}{10};
    % LARGE EVOKED PUPIL, Phasic response.
    all_Bp1_hi = mnkFigDat{1}{11}; all_Bp2_hi =mnkFigDat{2}{11}; all_Ap1_hi = mnkFigDat{1}{12}; all_Ap2_hi = mnkFigDat{2}{12};
    % SMALL EVOKED PUPIL, No phasic response.
    all_Bnp1_lo = mnkFigDat{1}{13}; all_Bnp2_lo =mnkFigDat{2}{13}; all_Anp1_lo = mnkFigDat{1}{14}; all_Anp2_lo = mnkFigDat{2}{14};
    % LARGE EVOKED PUPIL, No phasic response.
    all_Bnp1_hi = mnkFigDat{1}{15}; all_Bnp2_hi =mnkFigDat{2}{15}; all_Anp1_hi = mnkFigDat{1}{16}; all_Anp2_hi = mnkFigDat{2}{16};
    
    %% 3. Phasic response, large pupil - no-phasic response, small pupil.
    
    %% Calculate plot things.
    
    % 3. Phasic response, large pupil - no-phasic response, small pupil.
    
    for kk = 1:length(bs)
        D1a = all_Anp1_lo(:,kk)-all_Bnp1_lo(:,kk); D1b = all_Anp2_lo(:,kk)-all_Bnp2_lo(:,kk);
        D1c = all_Ap1_hi(:,kk)-all_Bp1_hi(:,kk); D1d = all_Ap2_hi(:,kk)-all_Bp2_hi(:,kk);
        DD1 = D1c-D1a; DD1 = DD1(isfinite(DD1));
        DD2 = D1d-D1b; DD2 = DD2(isfinite(DD2));
        % pval1p(1,kk) = signrank(DD1); pval2p(1,kk) = signrank(DD2);
        pval_both_np(1,kk) = signrank([DD1; DD2]);
        pval_both_npRS(1,kk) = ranksum([D1c; D1d],[D1a; D1b]);
        diff_p5{kk} = DD1; diff_p6{kk} = DD2;
        
        
        % 4. No phasic response, large pupil - Phasic response, small pupil.
        
        D1aw = all_Ap1_lo(:,kk)-all_Bp1_lo(:,kk); D1bw = all_Ap2_lo(:,kk)-all_Bp2_lo(:,kk);
        D1cw = all_Anp1_hi(:,kk)-all_Bnp1_hi(:,kk); D1dw = all_Anp2_hi(:,kk)-all_Bnp2_hi(:,kk);
        DD1w = D1aw-D1cw; DD1w = DD1w(isfinite(DD1w));
        DD2w = D1bw-D1dw; DD2w = DD2w(isfinite(DD2w));
        % pval1p(2,kk) = signrank(DD1); pval2p(2,kk) = signrank(DD2);
        pval_both_np(2,kk) = signrank([DD1w; DD2w]);
        pval_both_npRS(2,kk) = ranksum([D1cw; D1dw],[D1aw; D1bw]);
        diff_p7{kk} = DD1w; diff_p8{kk} = DD2w;
        
        p_across_conds(kk) = ranksum([DD1; DD2],[DD1w; DD2w]);
        
    end
    
    % pause
    
    %
    %     for kk = 1:length(bs)
    %         D1a = all_Anp1_lo(:,kk)-all_Bnp1_lo(:,kk); D1b = all_Anp2_lo(:,kk)-all_Bnp2_lo(:,kk);
    %         D1c = all_Ap1_hi(:,kk)-all_Bp1_hi(:,kk); D1d = all_Ap2_hi(:,kk)-all_Bp2_hi(:,kk);
    %         DD1 = D1c-D1a; DD1 = DD1(isfinite(DD1));
    %         DD2 = D1d-D1b; DD2 = DD2(isfinite(DD2));
    %
    %         % pval1pB(1,kk) = signrank(DD1); pval2pB(1,kk) = signrank(DD2);
    %
    %         pval_both_npB(1,kk) = signrank([DD1; DD2]);
    %         diff_p5{kk} = DD1; diff_p6{kk} = DD2;
    %
    %
    %
    %     end
    %
    %     %% 4. No phasic response, large pupil - Phasic response, small pupil.
    %
    %     for kk = 1:length(bs)
    %         D1a = all_Ap1_lo(:,kk)-all_Bp1_lo(:,kk); D1b = all_Ap2_lo(:,kk)-all_Bp2_lo(:,kk);
    %         D1c = all_Anp1_hi(:,kk)-all_Bnp1_hi(:,kk); D1d = all_Anp2_hi(:,kk)-all_Bnp2_hi(:,kk);
    %         DD1 = D1a-D1c; DD1 = DD1(isfinite(DD1));
    %         DD2 = D1b-D1d; DD2 = DD2(isfinite(DD2));
    %
    %         % pval1pB(2,kk) = signrank(DD1); pval2pB(2,kk) = signrank(DD2);
    %
    %         pval_both_npB(2,kk) = signrank([DD1; DD2]);
    %         diff_p7{kk} = DD1; diff_p8{kk} = DD2;
    %
    %         %         pvalploB(kk) = ranksum([all_Ap1_lo(:,kk);all_Ap2_lo(:,kk)],[all_Bp1_lo(:,kk);all_Bp2_lo(:,kk)]);
    %         %         pvalnphiB(kk) = ranksum([all_Anp1_hi(:,kk);all_Anp2_hi(:,kk)],[all_Bnp1_hi(:,kk);all_Bnp2_hi(:,kk)]);
    %
    %     end
    
    %% Plot.
    
    xjit = [-10 10];
    allD = {diff_p5 diff_p6 diff_p7 diff_p8};
    title_text = {'Baseline Pupil'};
    colors = {[0 0 0] [0.5 0.5 0.5]};
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    %     allPvals_befAft = [pvalphi;pvalnplo];
    
    for pind = 1:2
        for kk = 1:length(bs)
            
            ii = 2*pind;
            yd = [allD{ii-1}{kk}; allD{ii}{kk}]; mean_yd = nanmedian(yd);
            
            if pval_both_npRS(pind,kk) < 0.05
                plot(xjit(pind)+xAx(kk),mean_yd,'ko','markerfacecolor',colors{pind},'markeredgecolor',colors{pind},'markersize',8);
            end
            if pval_both_npRS(pind,kk) >= 0.05
                plot(xjit(pind)+xAx(kk),mean_yd,'ko','markerfacecolor','none','markeredgecolor',colors{pind},'markersize',8);
            end
            
            [ciV ciV_m] = bootci(1000,{@median,yd},'alpha',0.05,'type','norm'); % mean_ze = nanmean(ci_ze_m);
            plot([xjit(pind)+xAx(kk) xjit(pind)+xAx(kk)],[ciV(1) ciV(2)],'-','color',colors{pind}); % 95% CI, bootstrapped.
            
        end
    end
    
    for kk = 1:length(bs)
        if p_across_conds(kk)<0.05 plot(xjit(pind)+xAx(kk),0.07,'*','color',colors{pind}); end
    end
    
    plot([xAx(1) xAx(end)],[0 0],'k--');
    xlim([100 1100]); ylim([-0.03 0.08]);
    % xlim([100 1100]); ylim([-0.05 0.15]);
    % ylim([-0.3 0.3]);
    set(gca,'fontname','arial','fontsize',10);
    title(title_text{1},'fontsize',12,'fontname','arial','fontweight','normal');
    ylabel({'ACC r_s_c';'After beep - before beep'},'fontsize',10,'fontname','arial');
    xlabel('bin size (ms)','fontsize',10,'fontname','arial');
    
    %     xjit = [-10 10];
    %     allD = {diff_p5 diff_p6 diff_p7 diff_p8};
    %     % title_text = {'LP, LC P - SP, LC NP', 'LP, LC P - LP, LC NP', 'LP, LC P - SP, LC P', 'SP, LC P - LP, LC NP', 'LP, LC NP - SP, LC NP', 'SP, LC P - SP, LC NP'};
    %     title_text = {'Baseline Pupil'};
    %     colors = {[0 0 0] [0.5 0.5 0.5]};
    %     %     allPvals_befAftB = [pvalphiB;pvalnploB];
    %
    %     axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %
    %     for pind = 1:2
    %         for kk = 1:length(bs)
    %
    %             ii = 2*pind;
    %             yd = [allD{ii-1}{kk}; allD{ii}{kk}]; mean_yd = nanmedian(yd);
    %
    %             %             if allPvals_befAftB(pind,kk) < 0.05
    %                 plot(xjit(pind)+xAx(kk),mean_yd,'ko','markerfacecolor',colors{pind},'markeredgecolor',colors{pind},'markersize',8);
    %                 %             end
    %                 %             if allPvals_befAftB(pind,kk) >= 0.05
    %                 %                 plot(xjit(pind)+xAx(kk),mean_yd,'ko','markerfacecolor','none','markeredgecolor',colors{pind},'markersize',8);
    %                 %             end
    %
    %             [ciV ciV_m] = bootci(1000,{@median,yd},'alpha',0.05,'type','norm'); % mean_ze = nanmean(ci_ze_m);
    %             plot([xjit(pind)+xAx(kk) xjit(pind)+xAx(kk)],[ciV(1) ciV(2)],'-','color',colors{pind}); % 95% CI, bootstrapped.
    %
    %
    %             if pind == 1
    %                 if pval_both_npB(pind,kk)<0.05 plot(xjit(pind)+xAx(kk),0.15,'*','color',colors{pind}); end
    %             end
    %
    %             if pind == 2
    %                 if pval_both_npB(pind,kk)<0.05 plot(xjit(pind)+xAx(kk),0.14,'*','color',colors{pind}); end
    %             end
    %         end
    %     end
    %
    %     plot([xAx(1) xAx(end)],[0 0],'k--');
    %     xlim([100 1100]); ylim([-0.05 0.15]);
    %     set(gca,'fontname','arial','fontsize',10);
    %     title(title_text{1},'fontsize',12,'fontname','arial','fontweight','normal');
    %     % ylabel({'change in beep linked r_s_c';'phasic re:no phasic LC response'},'fontsize',10,'fontname','arial');
    %     xlabel('bin size (ms)','fontsize',10,'fontname','arial');
    
    plot(135,0.07,'.','markersize',30,'color',colors{1});
    plot(135,0.06,'.','markersize',30,'color',colors{2});
    text(180,0.07,'Large pupil, LC phasic - small pupil, LC no phasic','fontname','arial','fontsize',8);
    text(180,0.06,' Small pupil, LC phasic - large pupil, LC no phasic','fontname','arial','fontsize',8);
    
    %% Plot.
    %
    %     xjit = 10;
    %
    %     allD = {diff_p5 diff_p6 diff_p1 diff_p2 diff_p3 diff_p4 diff_p7 diff_p8 diff_p9 diff_p10 diff_p11 diff_p12};
    %
    %     title_text = {'LB, LC P - SB, LC NP', 'LB, LC P - LB, LC NP', 'LB, LC P - SB, LC P', 'SB, LC P - LB, LC NP', 'LB, LC NP - SB, LC NP', 'SB, LC P - SB, LC NP'};
    %
    %     for pind = 1:6
    %         axes(axs(pind+6)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %         for kk = 1:length(bs)
    %
    %             ii = 2*pind;
    %             yd = [allD{ii-1}{kk}; allD{ii}{kk}]; mean_yd = nanmedian(yd);
    %             plot(-xjit+xAx(kk),mean_yd,'k.','markersize',20,'color','k');
    %
    %             [ciV ciV_m] = bootci(1000,{@median,yd},'alpha',0.05,'type','norm'); % mean_ze = nanmean(ci_ze_m);
    %             plot([-xjit+xAx(kk) -xjit+xAx(kk)],[ciV(1) ciV(2)],'-','markersize',20,'color','k'); % 95% CI, bootstrapped.
    %
    %             if pval_both_np(pind,kk)<0.05 plot(-xjit+xAx(kk),0.14,'k*'); end
    %
    %         end
    %
    %         title(title_text{pind},'fontsize',8,'fontname','arial','fontweight','normal');
    %
    %         %         if pind==1
    %         %             ylabel({'change in beep linked r_s_c';'phasic re:no phasic LC response'},'fontsize',8,'fontname','arial');
    %         %         end
    %         plot([xAx(1) xAx(end)],[0 0],'k--');
    %         xlim([100 1100]);
    %         % ylim([-0.5 0.7]); % Pre 051321.
    %         ylim([-0.1 0.15]); % 051321.
    %         set(gca,'fontname','arial','fontsize',8);
    %
    %         xlabel('bin size (ms)','fontsize',10);
    %
    %     end
    %
    %     % *****************************************************************
    
end

%%

Plot_Supp_YesNo = 1;
Plot_Supp_YesNo = 0;

if Plot_Supp_YesNo
    
    %% Figure 9. Supplementary Figure: correlation between evoked pupil and evoked spiking.
    
    figureNumber = 91; num = 91; wid = 17.6; hts = [7]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 2, 2, [12], 'Joshi and Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    d1 = mnkFigDat{1}{17}; d2 = mnkFigDat{2}{17};
    %     d1 = mnkFigDat{1}{30}; d2 = mnkFigDat{2}{30};
    
    minB = min([min(d1) min(d2)]);
    maxB = max([max(d1) max(d2)]);
    % binsX =linspace(minB,maxB,9);
    binsX =linspace(-0.5,0.5,11);
    n1 = hist(d1,binsX);
    n2 = hist(d2,binsX);
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    % plot(binsX,n1/sum(n1),'k-');
    bar(binsX,n1/sum(n1),'edgecolor','none','facecolor','k');
    plot([nanmedian(d1) nanmedian(d1)],[0 .3],':','color',[.6 .6 .6],'linewidth',2); % 0.1306
    plot([0 0],[0 .35],'k--');
    xlabel('correlation between evoked response and evoked pupil','fontsize',10,'fontname','arial');
    ylabel('proportion of LC neurons','fontsize',10,'fontname','arial');
    title('mnky Sp','fontsize',12,'fontname','arial');
    axis([-0.55 0.55 0 0.3]);
    signrank(d1); % 1.0217e-06
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    % plot(binsX,n2/sum(n2),'k-');
    bar(binsX,n2/sum(n2),'edgecolor','none','facecolor','k');
    plot([nanmedian(d2) nanmedian(d2)],[0 .3],':','color',[.6 .6 .6],'linewidth',2); % 0.0562
    plot([0 0],[0 .35],'k--');
    xlabel('correlation between evoked response and evoked pupil','fontsize',10,'fontname','arial');
    ylabel('proportion of LC neurons','fontsize',10,'fontname','arial');
    title('mnky Ci','fontsize',12,'fontname','arial');
    set(gca,'fontname','arial');
    axis([-0.55 0.55 0 0.3]);
    signrank(d2); % 8.7788e-04
    
    
    %% Figure 9. SUPPLEMENTARY. Pupil size for each condition:
    
    figureNumber = 92; num = 92; wid = 17.6; hts = [15]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1, 1, [12], 'Joshi&Gold, 2019', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    % ***************************
    % Baseline.
    
    pupil_low_phasic_1 = mnkFigDat{1}{18}; n1 = 2*ones(1,length(pupil_low_phasic_1));
    pupil_high_phasic_1 = mnkFigDat{1}{19}; n2 = ones(1,length(pupil_high_phasic_1));
    pupil_low_no_phasic_1 = mnkFigDat{1}{20}; n3 = 4*ones(1,length(pupil_low_no_phasic_1));
    pupil_high_no_phasic_1 = mnkFigDat{1}{21}; n4 = 3*ones(1,length(pupil_high_no_phasic_1));
    all_pdat1 = [pupil_high_phasic_1; pupil_low_phasic_1; pupil_high_no_phasic_1; pupil_low_no_phasic_1];
    NN1 = [n2 n1 n4 n3];
    med1 = nanmedian(pupil_high_phasic_1);
    
    pupil_low_phasic_2 = mnkFigDat{2}{18}; n1 = 2*ones(1,length(pupil_low_phasic_2));
    pupil_high_phasic_2 = mnkFigDat{2}{19}; n2 = ones(1,length(pupil_high_phasic_2));
    pupil_low_no_phasic_2 = mnkFigDat{2}{20}; n3 = 4*ones(1,length(pupil_low_no_phasic_2));
    pupil_high_no_phasic_2 = mnkFigDat{2}{21}; n4 = 3*ones(1,length(pupil_high_no_phasic_2));
    all_pdat2 = [pupil_high_phasic_2; pupil_low_phasic_2; pupil_high_no_phasic_2; pupil_low_no_phasic_2];
    NN2 = [n2 n1 n4 n3];
    med2 = nanmedian(pupil_high_phasic_2);
    
    all_pdat = [all_pdat1; all_pdat2];
    all_NN = [NN1 NN2];
    all_med = nanmedian([pupil_high_phasic_1; pupil_high_phasic_2]);
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    boxplot(all_pdat,all_NN);
    plot([0 4.5],[0 0],'k:');
    % plot([0 4.5],[all_med all_med],'k:');
    title('Beep trials, peak pupil');
    xticklabels({'LC phasic, large evoked pupil','LC phasic, small evoked pupil','LC no phasic, large evoked pupil','LC no phasic, small evoked pupil'});
    xtickangle(45);
    ylim([-4 4]);
    
    
    % ***************************
    % Baseline.
    pupil_low_phasic_1 = mnkFigDat{1}{26}; n1 = 2*ones(1,length(pupil_low_phasic_1));
    pupil_high_phasic_1 = mnkFigDat{1}{27}; n2 = ones(1,length(pupil_high_phasic_1));
    pupil_low_no_phasic_1 = mnkFigDat{1}{28}; n3 = 4*ones(1,length(pupil_low_no_phasic_1));
    pupil_high_no_phasic_1 = mnkFigDat{1}{29}; n4 = 3*ones(1,length(pupil_high_no_phasic_1));
    all_pdat1 = [pupil_high_phasic_1; pupil_low_phasic_1; pupil_high_no_phasic_1; pupil_low_no_phasic_1];
    NN1 = [n2 n1 n4 n3];
    med1 = nanmedian(pupil_high_phasic_1);
    
    pupil_low_phasic_2 = mnkFigDat{2}{26}; n1 = 2*ones(1,length(pupil_low_phasic_2));
    pupil_high_phasic_2 = mnkFigDat{2}{27}; n2 = ones(1,length(pupil_high_phasic_2));
    pupil_low_no_phasic_2 = mnkFigDat{2}{28}; n3 = 4*ones(1,length(pupil_low_no_phasic_2));
    pupil_high_no_phasic_2 = mnkFigDat{2}{29}; n4 = 3*ones(1,length(pupil_high_no_phasic_2));
    all_pdat2 = [pupil_high_phasic_2; pupil_low_phasic_2; pupil_high_no_phasic_2; pupil_low_no_phasic_2];
    NN2 = [n2 n1 n4 n3];
    med2 = nanmedian(pupil_high_phasic_2);
    
    all_pdat = [all_pdat1; all_pdat2];
    all_NN = [NN1 NN2];
    all_med = nanmedian([pupil_high_phasic_1; pupil_high_phasic_2]);
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    boxplot(all_pdat,all_NN);
    plot([0 4.5],[0 0],'k:');
    % plot([0 4.5],[all_med all_med],'k:');
    title('Beep trials, baseline pupil');
    xticklabels({'LC phasic, large baseline pupil','LC phasic, small baseline pupil','LC no phasic, large baseline pupil','LC no phasic, small baseline pupil'});
    xtickangle(45);
    ylim([-4 4]);
    
    %%
    
    %     p1 = signrank(([pupil_high_phasic_1 pupil_high_phasic_2] - [pupil_low_phasic_1 pupil_low_phasic_2])./[pupil_low_phasic_1 pupil_low_phasic_2]);
    %     p2 = signrank(([pupil_high_phasic_1 pupil_high_phasic_2] - [pupil_high_no_phasic_1 pupil_high_no_phasic_2])./[pupil_high_no_phasic_1 pupil_high_no_phasic_2]);
    %     p3 = signrank(([pupil_high_phasic_1 pupil_high_phasic_2] - [pupil_low_no_phasic_1 pupil_low_no_phasic_2])./[pupil_low_no_phasic_1 pupil_low_no_phasic_2]);
    
    prop_both_phasic = [mnkFigDat{1}{22} mnkFigDat{2}{22}];
    prop_both_no_phasic = [mnkFigDat{1}{23} mnkFigDat{2}{23}];
    prop_only_one_phasic = [mnkFigDat{1}{24} mnkFigDat{2}{24}];
    prop_any_one_phasic = [mnkFigDat{1}{25} mnkFigDat{2}{25}];
    
    % pause
    
    % 16+-1
    m_bp = 100*nanmean(prop_both_phasic);
    sem_bp = 100*nanse(prop_both_phasic,2);
    
    % 17+-1
    m_bnp = 100*nanmean(prop_both_no_phasic);
    sem_bnp = 100*nanse(prop_both_no_phasic,2);
    
    % 47+-1
    m_oop = 100*nanmean(prop_only_one_phasic);
    sem_oop = 100*nanse(prop_only_one_phasic,2);
    
    % 63+-1
    m_aop = 100*nanmean(prop_any_one_phasic);
    sem_aop = 100*nanse(prop_any_one_phasic,2);
    
    %%
    
end % plotYesNo.





