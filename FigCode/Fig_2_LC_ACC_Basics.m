% Fig2_LC_ACC_Basics.m
%
% *****************************
% INFO: Calculate ACC rsc conditioned on LC spiking.
% NO BEEP trials only.
% *****************************
% "rsc" is the term used by Bair, Kohn et al for pairwise spike-count (sc) correlations (r).
% Dual area recording - Analysis script; FIXATION task.
% Calls "getLC_DualArea_rsc" and "getACC_DualArea_rsc" to do the heavy lifting.
%
% Calculates pairwise correlations (pwc) between ACC neurons.
% Plots distributions of spike rates in ACC, conditioned on zero, low or high LC counts.
% Saves pwc conditioned on LC/ACC low vs high spiking regimes.
% Plots variability of trial spike counts and trial ISI's using Fano Factor (CV dropped for now).
% *****************************
% HISTORY.
% *****************************
% Origin: 062116 - Sidd.
%
% Mod from: Fig4_LC_ACC_Dual_Recording_Firing_Rate.m
% 012820: Changed marker types and added bootstrap analysis for errors.
% 120219: Modified to add plot of variance of ACC  spike counts; also, create LC conditioned ACC figure.
% 080719: Modified to use "getLC_DualArea_rsc_2" - to address reviewer concern
% Instead of LC zero/nonzero split, we will try dividing trials based on spike counts of 1,2,3,4.... see if that is possible.
%
% 011219: Modified - now split trials based on spike/no spike rather than zero/median split.
% 011819: Finalizing figure layouts. Ha ha. Finalizing....
% 121218. Modified from: LC_Dual_Recording_rsc.m and ACC_Dual_Recording_rsc.m

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

min_trials = 5; % a
% min_trials = 10; % b

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
                    % ff = 28% For troubleshooting.
                    load(fnames{ff});
                    %                     pause % For troubleshooting.
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
    % save LC_ACC_Results_rsc_032521 mnkRes; % With new way of dividing trials into zero, mid, high.
    % save LC_ACC_Results_rsc_042921d mnkRes; % With new way of dividing trials into zero, mid, high.
    %     save LC_ACC_Results_rsc_043021b mnkRes; % With new way of dividing trials into zero, mid, high.
    % save LC_ACC_Results_rsc_052521 mnkRes; % With new way of dividing trials into zero, mid, high.
    save LC_ACC_Results_rsc_100421 mnkRes; % With new way of dividing trials into zero, mid, high.
    
    clear mnkRes;
end % If reanalyze LC loop.

% pause

%% ACC linked LC measurements.

min_trials = 5;

if reanalyzeACC
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
                    % ff = 35
                    load(fnames{ff});
                    %                     pause % For troubleshooting.
                    rsc_dat{ff} = getACC_DualArea_rsc2(siteData,min_trials,ff); % Added analyses for window size effects.
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
    savedir = strcat([base_dir,'\Results\Results_2020\LC_ACC']); cd(savedir);
    
    % save ACC_LC_Results_rsc_100120 mnkRes;
    save ACC_LC_Results_rsc_051621 mnkRes;
    
    clear mnkRes;
end % If reanalyze ACC loop.


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
    
    % base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2020\LC_ACC'; % Base directory for brain area.
    base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC'; % Base directory for brain area.
    savedir = base_dir;
    
    % inFile = 'LC_ACC_Results_rsc_021720f'; % New way of dividing ACC trials into low, mid, high.
    % inFile = 'LC_ACC_Results_rsc_032220b'; % New way of dividing ACC trials into low, mid, high.
    % inFile = 'LC_ACC_Results_rsc_032320'; % New way of dividing ACC trials into low, mid, high.
    % inFile = 'LC_ACC_Results_rsc_042020'; % New way of dividing ACC trials into low, mid, high.
    % inFile = 'LC_ACC_Results_rsc_042320'; % New way of dividing ACC trials into low, mid, high.
    % inFile = 'LC_ACC_Results_rsc_071820a';
    % inFile = 'LC_ACC_Results_rsc_100120';
    % inFile = 'LC_ACC_Results_rsc_042821d';
    % inFile = 'LC_ACC_Results_rsc_042921d';
    % inFile = 'LC_ACC_Results_rsc_043021b';
    % inFile = 'LC_ACC_Results_rsc_052521';
    inFile = 'LC_ACC_Results_rsc_100421';
    
    cd(savedir); load(inFile);
    
    monks = {'Sprout', 'Cicero'}; % Add mnks as needed...
    sites = {'LC_ACC_Fixation'}; % For now LC only.
    nMonks = length(monks); nSites = length(sites);
    
    % **********************************************
    oki=0; pairs_mnk = []; unit_all_lo = []; unit_all_hi = []; unit_all_all = [];
    % *****************************************
    
    %     nMonks = 1;
    for mm = 1:nMonks % Loop over mnks....
        
        all_ze_c_meanT = []; all_hi_c_meanT = []; all_all_c_meanT = []; % Spike count matrix.
        var_zeT = []; var_hiT = []; var_ALLT = []; % Variance of spike counts.
        
        unit_fii = []; unit_lci = []; unit_api = [];
        
        unit_all_1T = []; unit_all_2T = [];
        unit_all_3T = []; unit_all_4T = [];
        
        unit_all_zeT = []; unit_all_hiT = []; unit_all_allT = []; % Pairwise spike count correlation, rsc.
        unit_all_ze_pvT = []; unit_all_hi_pvT = []; unit_all_all_pvT = []; % p-value from crosscorr calculation of rsc.
        % unit_all_ze_pT = []; unit_all_hi_pT = []; unit_all_all_pT = []; % rsc values from shuffled trials.
        unit_all_ze_pT = []; unit_all_hi_pT = []; unit_all_all_pT = []; % rsc values from shuffled trials.
        unit_all_rate_hi = [];
        unit_allZi = []; unit_allHi = [];
        acc_bc_zeT = []; acc_bc_hiT = []; acc_bc_allT = []; % Mean binned spike count for each ACC unit.
        unit_dat_rate_zeT = []; unit_dat_rate_hiT = []; unit_dat_rate_alT = []; % ACC neuron spike rates.
        unit_dat_var_zeT = []; unit_dat_var_hiT = []; unit_dat_var_alT = []; % ACC neuron spike rates.
        
        FF_b_zeT = []; FF_b_allT = []; FF_b_hiT = []; FF_b_1T = []; FF_b_2T = []; FF_b_3T = []; FF_b_4T = [];
        var_b_zeT = []; var_b_allT = []; var_b_hiT = []; var_b_1T = []; var_b_2T = []; var_b_3T = []; var_b_4T = [];
        mean_b_zeT = []; mean_b_allT = []; mean_b_hiT = [];  mean_b_1T = [];  mean_b_2T = [];  mean_b_3T = [];  mean_b_4T = [];
        all_ze_c_medianT = []; all_lo_c_medianT = []; all_hi_c_medianT = []; all_all_c_medianT = [];
        std_zeT = []; std_loT = []; std_hiT = []; std_allT = []; %  std_mid{2} T = [];
        all_ze_count1T =[]; all_lo_count1T =[]; all_hi_count1T =[]; all_ze_count2T =[]; all_lo_count2T =[]; all_hi_count2T =[];
        
        mnkDat = mnkRes{mm};
        
        lci = 0;
        lc_u_lo = [];
        lc_u_hi = [];
        pair_ind = [];
        
        pairs = [];
        
        for ss = 1:nSites % Loop over brainstem sites.... LC and/or IC.
            
            sitDat = mnkDat{ss};
            nnuu = 0; nnAA = 0; % Housekeeping - how many LC and ACC units did we have?
            
            if ~isempty(sitDat)
                for ap = 1:length(sitDat) % Loop over sessions.
                    dat = sitDat{ap}; % Get session data.
                    if ~isempty(dat)
                        allZ = []; allH = []; allA = [];
                        % acc_count = dat{3};
                        % nnuu = nnuu + length(dat{1});
                        nnAA = nnAA + length(dat{2});
                        
                        % *********************************************
                        % Get pair indices.
                        pair_ind = dat{6};
                        
                        for sp = 1:length(dat{1}) % Loop over LC units in this session (can be 1 or more).
                            
                            LC_lohi = dat{4}{sp};
                            LC_nz = dat{7}{sp}; % Added: 112319.
                            % pause
                            zei = find(LC_lohi==0);
                            onei = find(LC_lohi==1);
                            twoi = find(LC_lohi==2);
                            threei = find(LC_lohi==3);
                            fouri = find(LC_lohi==4);
                            % pause
                            hii = LC_nz;
                            
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
                                
                                all_ze = []; all_hi = []; all_all = [];
                                all_1 = []; all_2 = []; all_3 = []; all_4 = []; all_nz = [];
                                all_ze_pv = []; all_hi_pv = []; all_all_pv = [];
                                all_ze_p = []; all_hi_p = []; all_all_p = [];
                                dat_rate_ze = []; dat_rate_hi = []; dat_rate_al = [];
                                dat_var_ze = []; dat_var_hi = []; dat_var_al = [];
                                all_hi_rate = [];
                                
                                dat0 = dat{1}{sp};
                                
                                dat1_LC_spikes = dat{5}{sp}; % LC trial spike counts.
                                % LC_sp_hi = (dat1_LC_spikes(LC_lohi==1));
                                LC_sp_hi = dat1_LC_spikes; % 100219; get all spike counts for distribution.
                                
                                % all_hi_rate = [all_hi_rate LC_sp_hi];
                                
                                
                                if ~isempty(dat0)
                                    fii = dat0{1}{1}(:,7);
                                    lci = dat0{1}{1}(:,8);
                                    api = dat0{1}{1}(:,9);
                                end
                                if isempty(dat0)
                                    fii = [];
                                    lci = [];
                                    api = [];
                                end
                                
                                for bhi= 1:length(dat0) % Iterate over bins.
                                    dat1 = dat0{bhi}{1}; % Pairwise spike-count correlation results.
                                    dat1_pze = dat0{bhi}{3}; % Zero LC. Pairwise spike-count correlation results - shuffle tests.
                                    dat1_phi = dat0{bhi}{4}; % High LC. Pairwise spike-count correlation results - shuffle tests.
                                    dat1_pall = dat0{bhi}{5}; % All trials. Pairwise spike-count correlation results - shuffle tests.
                                    
                                    % pause
                                    
                                    % For binsize 800, get the ACC pair spike rate.
                                    if bhi ==4
                                        dat_rate = dat0{bhi}{2}; % Get mean spike rate of each neuron in each pair.
                                        dat_var = dat0{bhi}{9}; % Get mean spike rate of each neuron in each pair.
                                        for pp = 1:length(dat_rate) % Iterate over ACC pairs.
                                            s_rate = dat_rate{pp};
                                            s_var = dat_var{pp};
                                            dat_rate_ze = [dat_rate_ze; s_rate(1) s_rate(2)];
                                            dat_rate_hi = [dat_rate_hi; s_rate(3) s_rate(4)];
                                            dat_rate_al = [dat_rate_al; s_rate(5) s_rate(6)];
                                            dat_var_ze = [dat_var_ze; s_var(1) s_var(2)];
                                            dat_var_hi = [dat_var_hi; s_var(3) s_var(4)];
                                        end
                                    end
                                    
                                    % Organize the pairwise corr data:
                                    
                                    % all_p = find(dat1(:,6)>0.05);
                                    % dat1(all_p,1) = nan; dat1(all_p,2) = nan;  dat1(all_p,3) = nan;
                                    % dat1(all_p,10) = nan; dat1(all_p,11) = nan; dat1(all_p,12) = nan; dat1(all_p,13) = nan;
                                    
                                    all_ze = [all_ze dat1(:,1)]; all_ze_pv = [all_ze_pv dat1(:,4)];
                                    all_hi = [all_hi dat1(:,2)]; all_hi_pv = [all_hi_pv dat1(:,5)];
                                    all_all = [all_all dat1(:,3)]; all_all_pv = [all_all_pv dat1(:,6)];
                                    
                                    % Added 080819.
                                    % all_1p = find(dat1(:,14)<0.05); all_2p = find(dat1(:,15)<0.05); all_3 = find(dat1(:,16)<0.05); all_4 = find(dat1(:,17)<0.05);
                                    all_1 = [all_1 dat1(:,10)]; all_2 = [all_2 dat1(:,11)]; all_3 = [all_3 dat1(:,12)]; all_4 = [all_4 dat1(:,13)];
                                    
                                    
                                    % % Shuffled trials.
                                    % unit_all_ze_pT{bhi} = [unit_all_ze_pT{bhi} dat1_pze'];
                                    % unit_all_hi_pT{bhi} = [unit_all_hi_pT{bhi} dat1_phi'];
                                    % unit_all_all_pT{bhi} = [unit_all_all_pT{bhi} dat1_pall'];
                                    % Shuffled trials.
                                    all_ze_p = [all_ze_p dat1_pze'];
                                    all_hi_p = [all_hi_p dat1_phi'];
                                    all_all_p = [all_all_p dat1_pall'];
                                    
                                end % Loop over number of bins for rsc data.
                                
                                
                                % Collect the ACC pair spike rates for this LC unit.
                                % if bhi == 6
                                % bhi
                                % pause
                                
                                % if bhi == 4
                                unit_dat_rate_zeT = [unit_dat_rate_zeT; dat_rate_ze];
                                unit_dat_rate_hiT = [unit_dat_rate_hiT; dat_rate_hi];
                                unit_dat_rate_alT = [unit_dat_rate_alT; dat_rate_al];
                                unit_dat_var_zeT = [unit_dat_var_zeT; dat_var_ze];
                                unit_dat_var_hiT = [unit_dat_var_hiT; dat_var_hi];
                                % end
                                
                                % Collect rsc for this binsize, this LC unit.
                                unit_all_zeT = [unit_all_zeT; all_ze]; unit_all_ze_pvT = [unit_all_ze_pvT; all_ze_pv];
                                unit_all_hiT = [unit_all_hiT; all_hi];   unit_all_hi_pvT = [unit_all_hi_pvT; all_hi_pv];
                                unit_all_allT = [unit_all_allT; all_all]; unit_all_all_pvT = [unit_all_all_pvT; all_all_pv];
                                
                                unit_all_1T = [unit_all_1T; all_1];
                                unit_all_2T = [unit_all_2T; all_2];
                                unit_all_3T = [unit_all_3T; all_3];
                                unit_all_4T = [unit_all_4T; all_4];
                                
                                
                                % Shuffled trials.
                                unit_all_ze_pT = [unit_all_ze_pT; all_ze_p];
                                unit_all_hi_pT = [unit_all_hi_pT; all_hi_p];
                                unit_all_all_pT = [unit_all_all_pT; all_all_p];
                                
                                unit_fii = [unit_fii; fii]; unit_lci = [unit_lci; lci]; unit_api = [unit_api; api];
                                
                                unit_all_rate_hi = [unit_all_rate_hi LC_sp_hi];
                                
                                unit_allZi = [unit_allZi zei];
                                unit_all1i = [unit_allZi onei];
                                unit_all2i = [unit_allZi twoi];
                                unit_all3i = [unit_allZi threei];
                                unit_all4i = [unit_allZi fouri];
                                % unit_allHi = [unit_allHi hii];
                                
                                % if size(unit_all_zeT,1)>178
                                % if size(unit_all_zeT,1)>63
                                %                                 if size(unit_all_zeT,1)>1112
                                %                                     oki=oki+1;
                                %                                     if oki==1
                                %                                         % eg_i = size(all_ze,1) - (size(unit_all_zeT,1)-178);
                                %                                         % eg_i = size(all_ze,1) - (size(unit_all_zeT,1)-63);
                                %                                         eg_i = size(all_ze,1) - (size(unit_all_zeT,1)-1112);
                                %                                         pair_i = pair_ind(eg_i,:);
                                %                                     end
                                %                                 end
                                
                                % pause
                                
                                % *********************************************
                                % Now get binned ACC spike count data to calculate FF.
                                
                                nMinTrials = 10;
                                % if (length(zei)>nMinTrials && length(onei)>nMinTrials && length(threei)>nMinTrials && length(fouri)>nMinTrials && length(hii)>nMinTrials)
                                
                                ACC_binned = dat{3}; % Binned ACC unit spike counts.
                                % Loop over ACC units for spike counts and calculate FF.
                                for ac = 1:size(ACC_binned,1) % Loop over ACC units to get binned spike counts conditioned on this LC unit.
                                    
                                    acu_bd = ACC_binned(ac,:);
                                    
                                    FF_b_hi_temp = []; FF_b_ze_temp = []; FF_b_all_temp = [];
                                    FF_b_1_temp = []; FF_b_2_temp = []; FF_b_3_temp = []; FF_b_4_temp = [];
                                    
                                    mean_b_hi_temp = []; mean_b_ze_temp = []; mean_b_all_temp = [];
                                    mean_b_1_temp = []; mean_b_2_temp = []; mean_b_3_temp = []; mean_b_4_temp = [];
                                    
                                    var_b_hi_temp = []; var_b_ze_temp = []; var_b_all_temp = [];
                                    var_b_1_temp = []; var_b_2_temp = []; var_b_3_temp = []; var_b_4_temp = [];
                                    
                                    acc_uc_ze = []; acc_uc_hi = []; acc_uc_all = [];
                                    acc_uc_1 = []; acc_uc_2 = []; acc_uc_3 = []; acc_uc_4 = [];
                                    
                                    
                                    for acb = 1:length(acu_bd) % Loop over binsizes.
                                        
                                        % pause
                                        
                                        % ************************************
                                        % ALL:
                                        acu_bd_all = acu_bd{acb}; % pause
                                        mbc_v = reshape(acu_bd_all,size(acu_bd_all,1)*size(acu_bd_all,2),1);
                                        mbc_mean = nanmean(mbc_v);
                                        mbc_var = nanvar(mbc_v);
                                        mbc = nanmean(acu_bd_all,2);
                                        
                                        acc_uc_all = [acc_uc_all mbc_mean];
                                        
                                        % if min(mbc)>=0
                                        if mbc_mean>0 FF_b_all_temp = [FF_b_all_temp mbc_var/mbc_mean]; end
                                        if mbc_mean==0 FF_b_all_temp = [FF_b_all_temp nan]; end
                                        mean_b_all_temp = [mean_b_all_temp mbc_mean];
                                        var_b_all_temp = [var_b_all_temp mbc_var];
                                        % end
                                        
                                        % ***********************************
                                        
                                        % ************************************
                                        % ZERO:
                                        acu_bd_ze = acu_bd{acb}(zei,:); % pause
                                        % mbc_v = reshape(acu_bd_ze,size(acu_bd_ze,1)*size(acu_bd_ze,2),1);
                                        mbc_v = acu_bd_ze(:);
                                        mbc_mean = nanmean(mbc_v);
                                        mbc_var = nanvar(mbc_v);
                                        mbc = nanmean(acu_bd_ze,2);
                                        
                                        acc_uc_ze = [acc_uc_ze mbc_mean];
                                        
                                        if length(zei)>=nMinTrials
                                            % if min(mbc)>=0
                                            if mbc_mean>0 FF_b_ze_temp = [FF_b_ze_temp mbc_var/mbc_mean]; end
                                            if mbc_mean==0 FF_b_ze_temp = [FF_b_ze_temp nan]; end
                                            mean_b_ze_temp = [mean_b_ze_temp mbc_mean];
                                            var_b_ze_temp = [var_b_ze_temp mbc_var];
                                            % end
                                        end
                                        if length(zei)<nMinTrials
                                            mean_b_ze_temp = [mean_b_ze_temp nan];
                                            var_b_ze_temp = [var_b_ze_temp nan];
                                            FF_b_ze_temp = [FF_b_ze_temp nan];
                                        end
                                        
                                        % ***********************************
                                        
                                        % One LC spike.
                                        acu_bd_1 = acu_bd{acb}(onei,:); % pause.
                                        % mbc_v = reshape(acu_bd_1,size(acu_bd_1,1)*size(acu_bd_1,2),1);
                                        mbc_v = acu_bd_1(:);
                                        mbc_mean = nanmean(mbc_v);
                                        mbc_var = nanvar(mbc_v);
                                        mbc = nanmean(acu_bd_1,2);
                                        
                                        acc_uc_1 = [acc_uc_1 mbc_mean];
                                        
                                        if length(onei)>=nMinTrials
                                            % if min(mbc)>=0
                                            if mbc_mean>0 FF_b_1_temp = [FF_b_1_temp mbc_var/mbc_mean]; end
                                            if mbc_mean==0 FF_b_1_temp = [FF_b_1_temp nan]; end
                                            mean_b_1_temp = [mean_b_1_temp mbc_mean];
                                            var_b_1_temp = [var_b_1_temp mbc_var];
                                            % end
                                        end
                                        if length(onei)<nMinTrials
                                            mean_b_1_temp = [mean_b_1_temp nan];
                                            var_b_1_temp = [var_b_1_temp nan];
                                            FF_b_1_temp = [FF_b_1_temp nan];
                                        end
                                        
                                        % ************************************
                                        
                                        % Two LC spikes.
                                        acu_bd_2 = acu_bd{acb}(twoi,:); % pause.
                                        % mbc_v = reshape(acu_bd_2,size(acu_bd_2,1)*size(acu_bd_2,2),1);
                                        mbc_v = acu_bd_2(:);
                                        mbc_mean = nanmean(mbc_v);
                                        mbc_var = nanvar(mbc_v);
                                        mbc = nanmean(acu_bd_2,2);
                                        
                                        acc_uc_2 = [acc_uc_2 mbc_mean];
                                        
                                        if length(twoi)>=nMinTrials
                                            % if min(mbc)>=0
                                            if mbc_mean>0 FF_b_2_temp = [FF_b_2_temp mbc_var/mbc_mean]; end
                                            if mbc_mean==0 FF_b_2_temp = [FF_b_2_temp nan]; end
                                            mean_b_2_temp = [mean_b_2_temp mbc_mean];
                                            var_b_2_temp = [var_b_2_temp mbc_var];
                                            % end
                                        end
                                        if length(twoi)<nMinTrials
                                            mean_b_2_temp = [mean_b_2_temp nan];
                                            var_b_2_temp = [var_b_2_temp nan];
                                            FF_b_2_temp = [FF_b_2_temp nan];
                                        end
                                        
                                        % ************************************
                                        
                                        % Three LC spikes.
                                        acu_bd_3 = acu_bd{acb}(threei,:); % pause.
                                        % mbc_v = reshape(acu_bd_3,size(acu_bd_3,1)*size(acu_bd_3,2),1);
                                        mbc_v = acu_bd_3(:);
                                        mbc_mean = nanmean(mbc_v);
                                        mbc_var = nanvar(mbc_v);
                                        mbc = nanmean(acu_bd_3,2);
                                        
                                        acc_uc_3 = [acc_uc_3 mbc_mean];
                                        
                                        if length(threei)>=nMinTrials
                                            % if min(mbc)>=0
                                            if mbc_mean>0 FF_b_3_temp = [FF_b_3_temp mbc_var/mbc_mean]; end
                                            if mbc_mean==0 FF_b_3_temp = [FF_b_3_temp nan]; end
                                            mean_b_3_temp = [mean_b_3_temp mbc_mean];
                                            var_b_3_temp = [var_b_3_temp mbc_var];
                                            % end
                                        end
                                        if length(threei)<nMinTrials
                                            mean_b_3_temp = [mean_b_3_temp nan];
                                            var_b_3_temp = [var_b_3_temp nan];
                                            FF_b_3_temp = [FF_b_3_temp nan];
                                        end
                                        
                                        % ************************************
                                        
                                        % 4 or more LC spikes.
                                        acu_bd_4 = acu_bd{acb}(fouri,:); % pause.
                                        % mbc_v = reshape(acu_bd_4,size(acu_bd_4,1)*size(acu_bd_4,2),1);
                                        mbc_v = acu_bd_4(:);
                                        mbc_mean = nanmean(mbc_v);
                                        mbc_var = nanvar(mbc_v);
                                        mbc = nanmean(acu_bd_4,2);
                                        
                                        acc_uc_4 = [acc_uc_4 mbc_mean];
                                        
                                        if length(fouri)>=nMinTrials
                                            % if min(mbc)>=0
                                            if mbc_mean>0 FF_b_4_temp = [FF_b_4_temp mbc_var/mbc_mean]; end
                                            if mbc_mean==0 FF_b_4_temp = [FF_b_4_temp nan]; end
                                            mean_b_4_temp = [mean_b_4_temp mbc_mean];
                                            var_b_4_temp = [var_b_4_temp mbc_var];
                                            % end
                                        end
                                        if length(fouri)<nMinTrials
                                            mean_b_4_temp = [mean_b_4_temp nan];
                                            var_b_4_temp = [var_b_4_temp nan];
                                            FF_b_4_temp = [FF_b_4_temp nan];
                                        end
                                        
                                        % ************************************
                                        
                                        % Nonzero LC spikes:
                                        acu_bd_hi = acu_bd{acb}(hii,:); % pause.
                                        % mbc_v = reshape(acu_bd_hi,size(acu_bd_hi,1)*size(acu_bd_hi,2),1);
                                        mbc_v = acu_bd_hi(:);
                                        mbc_mean = nanmean(mbc_v);
                                        mbc_var = nanvar(mbc_v);
                                        mbc = nanmean(acu_bd_hi,2);
                                        
                                        acc_uc_hi = [acc_uc_hi mbc_mean];
                                        
                                        if length(hii)>=nMinTrials
                                            % if min(mbc)>=0
                                            if mbc_mean>0 FF_b_hi_temp = [FF_b_hi_temp mbc_var/mbc_mean]; end
                                            if mbc_mean==0 FF_b_hi_temp = [FF_b_hi_temp nan]; end
                                            mean_b_hi_temp = [mean_b_hi_temp mbc_mean];
                                            var_b_hi_temp = [var_b_hi_temp mbc_var];
                                            % end
                                        end
                                        if length(hii)<nMinTrials
                                            mean_b_hi_temp = [mean_b_hi_temp nan];
                                            var_b_hi_temp = [var_b_hi_temp nan];
                                            FF_b_hi_temp = [FF_b_hi_temp nan];
                                        end
                                        
                                        % ************************************
                                        
                                    end % Binsizes.
                                    
                                    % pause
                                    acc_bc_zeT = [acc_bc_zeT; acc_uc_ze];
                                    acc_bc_hiT = [acc_bc_hiT; acc_uc_hi];
                                    acc_bc_allT = [acc_bc_allT; acc_uc_all];
                                    
                                    FF_b_zeT = [FF_b_zeT; FF_b_ze_temp];
                                    FF_b_1T = [FF_b_1T; FF_b_1_temp];
                                    FF_b_2T = [FF_b_2T; FF_b_2_temp];
                                    FF_b_3T = [FF_b_3T; FF_b_3_temp];
                                    FF_b_4T = [FF_b_4T; FF_b_4_temp];
                                    FF_b_hiT = [FF_b_hiT; FF_b_hi_temp];
                                    FF_b_allT = [FF_b_allT; FF_b_all_temp];
                                    
                                    mean_b_zeT = [mean_b_zeT; mean_b_ze_temp];
                                    mean_b_1T = [mean_b_1T; mean_b_1_temp];
                                    mean_b_2T = [mean_b_2T; mean_b_2_temp];
                                    mean_b_3T = [mean_b_3T; mean_b_3_temp];
                                    mean_b_4T = [mean_b_4T; mean_b_4_temp];
                                    mean_b_hiT = [mean_b_hiT; mean_b_hi_temp];
                                    mean_b_allT = [mean_b_allT; mean_b_all_temp];
                                    
                                    var_b_zeT = [var_b_zeT; var_b_ze_temp];
                                    var_b_1T = [var_b_1T; var_b_1_temp];
                                    var_b_2T = [var_b_2T; var_b_2_temp];
                                    var_b_3T = [var_b_3T; var_b_3_temp];
                                    var_b_4T = [var_b_4T; var_b_4_temp];
                                    var_b_hiT = [var_b_hiT; var_b_hi_temp];
                                    var_b_allT = [var_b_allT; var_b_all_temp];
                                    
                                end % ACC units is did.
                                
                                dat2 = dat{2}; % ACC unit spike counts - over whole window.
                                
                                for ac = 1:length(dat2) % Loop over ACC units to get counts for this brainstem unit.
                                    
                                    % Get ACC counts corresponding to this brainstem unit.
                                    dat3 = dat2{ac};
                                    allZ = [allZ dat3(zei)]; allH = [allH dat3(hii)]; allA = [allA dat3];
                                    
                                    % Calculate the mean spike count over trials for each ACC unit, for this LC unit.
                                    all_ze_c_meanT = [all_ze_c_meanT nanmean(dat3(zei))];
                                    all_hi_c_meanT = [all_hi_c_meanT nanmean(dat3(hii))];
                                    all_all_c_meanT = [all_all_c_meanT nanmean(dat3)];
                                    
                                    all_ze_c_medianT = [all_ze_c_medianT nanmedian(dat3(zei))];
                                    all_hi_c_medianT = [all_hi_c_medianT nanmedian(dat3(hii))];
                                    all_all_c_medianT = [all_all_c_medianT nanmedian(dat3)];
                                    
                                    std_zeT = [std_zeT nanstd(dat3(zei))];
                                    var_zeT = [var_zeT nanvar(dat3(zei))];
                                    std_hiT = [std_hiT nanstd(dat3(hii))];
                                    var_hiT = [var_hiT nanvar(dat3(hii))];
                                    std_allT = [std_allT nanstd(dat3)];
                                    var_ALLT = [var_ALLT nanvar(dat3)];
                                    
                                end % All ACC units for this brainstem unit.
                            end
                        end % All brainstem units for this session.
                    end % Check for empty dat structure.
                end % All sessions.
            end % Check for empty sitDat structure.
        end % Sites.
        
        %         all_pairs{mm} = pair_ind;
        
        all_ze_c_mean{mm} = all_ze_c_meanT; all_hi_c_mean{mm} = all_hi_c_meanT; all_all_c_mean{mm} = all_all_c_meanT; % Spike count matrix.
        var_ze{mm} = var_zeT; var_hi{mm} = var_hiT; var_ALL{mm} = var_ALLT; % Variance of spike counts.
        
        unit_all_fii{mm} = unit_fii; unit_all_lci{mm} = unit_lci; unit_all_api{mm} = unit_api;
        
        unit_all_1{mm} = unit_all_1T; unit_all_2{mm} = unit_all_2T; unit_all_3{mm} = unit_all_3T; unit_all_4{mm} = unit_all_4T;
        
        unit_all_ze{mm} = unit_all_zeT; unit_all_hi{mm} = unit_all_hiT; unit_all_all{mm} = unit_all_allT; % Pairwise spike count correlation, rsc.
        unit_all_ze_pv{mm} = unit_all_ze_pvT; unit_all_hi_pv{mm} = unit_all_hi_pvT; unit_all_all_pv{mm} = unit_all_all_pvT; % p-value from crosscorr calculation of rsc.
        unit_all_ze_p{mm} = unit_all_ze_pT; unit_all_hi_p{mm} = unit_all_hi_pT; unit_all_all_p{mm} = unit_all_all_pT; % rsc values from shuffled trials.
        
        acc_bc_ze{mm} = acc_bc_zeT; acc_bc_hi{mm} = acc_bc_hiT; % Mean binned spike count for each ACC unit.
        unit_dat_rate_ze{mm} = unit_dat_rate_zeT; unit_dat_rate_hi{mm} = unit_dat_rate_hiT; unit_dat_rate_al{mm} = unit_dat_rate_alT; % ACC neuron spike rates.
        unit_dat_var_ze{mm} = unit_dat_var_zeT; unit_dat_var_hi{mm} = unit_dat_var_hiT; % unit_dat_rate_al{mm} = unit_dat_rate_alT; % ACC neuron FF.
        
        FF_b_ze{mm} = FF_b_zeT; FF_b_nz{mm} = FF_b_hiT; FF_b_all{mm} = FF_b_allT;
        FF_b_1{mm} = FF_b_1T; FF_b_2{mm} = FF_b_2T; FF_b_3{mm} = FF_b_3T; FF_b_4{mm} = FF_b_4T;
        
        var_b_ze{mm} = var_b_zeT; var_b_hi{mm} = var_b_hiT; var_b_all{mm} = var_b_allT;
        var_b_1{mm} = var_b_1T; var_b_2{mm} = var_b_2T; var_b_3{mm} = var_b_3T; var_b_4{mm} = var_b_4T;
        
        mean_b_ze{mm} = mean_b_zeT; mean_b_nz{mm} = mean_b_hiT; mean_b_all{mm} = mean_b_allT;
        mean_b_1{mm} = mean_b_1T; mean_b_2{mm} = mean_b_2T; mean_b_3{mm} = mean_b_3T; mean_b_4{mm} = mean_b_4T;
        
        all_ze_c_median{mm} = all_ze_c_medianT; all_hi_c_median{mm} = all_hi_c_medianT; all_all_c_median{mm} = all_all_c_medianT;
        std_ze{mm} = std_zeT; std_hi{mm} = std_hiT; std_all{mm} = std_allT;
        all_ze_count1{mm} = all_ze_count1T; all_hi_count1{mm} = all_hi_count1T;
        all_ze_count2{mm} =all_ze_count2T; all_hi_count2{mm} = all_hi_count2T;
        
        
        lc_sp_hi{mm} = unit_all_rate_hi;
        
        lc_AZ{mm} = unit_allZi;
        lc_AH{mm} = unit_allHi;
        
        nOut(mm) = nnuu;
        nOutA(mm) = nnAA;
        
        pairs_mnk{mm} = pairs;
        
        %         eg_ze = {ze_pair_small_1 ze_pair_large_1 ze_pair_small_2 ze_pair_large_2};
        %         eg_hi = {hi_pair_small_1 hi_pair_large_1 hi_pair_small_2 hi_pair_large_2};
        
    end % Mnks.
    
    eg_ze = [];
    eg_hi = [];
    
    save('LC_noBeep_rscData_100421','unit_all_ze','unit_all_hi','unit_all_all','eg_ze','eg_hi','unit_all_ze_p','unit_all_hi_p','unit_all_all_p','unit_all_1','unit_all_2','unit_all_3','unit_all_4','unit_dat_rate_ze','unit_dat_rate_hi','unit_dat_rate_al','unit_dat_var_ze','unit_dat_var_hi');
    % save('LC_noBeep_rscData_052521','unit_all_ze','unit_all_hi','unit_all_all','eg_ze','eg_hi','unit_all_ze_p','unit_all_hi_p','unit_all_all_p','unit_all_1','unit_all_2','unit_all_3','unit_all_4','unit_dat_rate_ze','unit_dat_rate_hi','unit_dat_rate_al','unit_dat_var_ze','unit_dat_var_hi');
    % save('LC_noBeep_rscData_043021b','unit_all_ze','unit_all_hi','unit_all_all','eg_ze','eg_hi','unit_all_ze_p','unit_all_hi_p','unit_all_all_p','unit_all_1','unit_all_2','unit_all_3','unit_all_4','unit_dat_rate_ze','unit_dat_rate_hi','unit_dat_rate_al');
    % save('LC_noBeep_indices_043021b','unit_all_all','unit_all_fii','unit_all_lci','unit_all_api');
    
end

% pause

%% Now PLOT!

% Don't plot results?
plotYesNo = false;
% Plot results?
plotYesNo = true;

if plotYesNo
    %% Main Figure 2: LC low/high; ACC spike counts, rates, variance and Fano factor.
    
    % ********************************************
    % Supplementary Figure 2 starts here.
    % ********************************************
    
    % Setup figure:
    figureNumber = 2; num = 2; wid = 17.6; hts = [6 4 4 4]; cols = {1 2 2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1, 2, [12], 'Joshi and Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    % bs = round(logspace(log10(100),log10(1000),10),1);
    xbs = linspace(200,1000,5); % 032220b.
    %     xbs = linspace(100,1000,10); % 032320.
    nbs = length(xbs);
    
    %% low/high spike rate distributions.
    
    %     pause
    
    nZ = 0; hi_sp = [];
    
    for jjm = 1:nMonks
        
        all_sp = lc_sp_hi{jjm};
        all_sp_Z = all_sp(all_sp==0);
        all_sp_NZ = all_sp(all_sp>0);
        nZ = nZ + length(all_sp_Z);
        % hi_sp = [hi_sp all_sp_NZ/1.1]; % hi_sp = hi_sp(~isnan(hi_rsc));
        hi_sp = [hi_sp all_sp_NZ/1.1]; % hi_sp = hi_sp(~isnan(hi_rsc));
        
    end
    
    % xMax = length(xbs); xMin = 1;
    % hist(hi_sp);
    % Look at hist of nonzero spike rates... find a sensible range to calculate distribution.
    xMax = 10; xMin = 1;
    xBins = linspace(xMin,xMax,11);
    
    prcHi = prctile(hi_sp,[25 50 75]);
    
    nhi = hist(hi_sp,xBins);
    NNN = sum(nhi)+ nZ;
    
    c2='k';
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    % plot(0,nZ/NNN,'o','markersize',8,'markerfacecolor','none','markeredgecolor',c2,'linewidth',2);
    % plot(0,nZ/NNN,'ms','markersize',8,'linewidth',2,'markerfacecolor','none');
    bar(0,nZ/NNN,'edgecolor','none','facecolor','m');
    
    % plot(xBins,nhi/NNN,'k-','linewidth',2);
    bar(xBins,nhi/NNN,'edgecolor','none','facecolor','k');
    
    legend({'zero','> 0'},'autoupdate','off','fontsize',8); legend('boxoff');
    % plot([prcHi(1) prcHi(3)],[0.3 0.3],'-','linewidth',2,'color',c2);
    plot([prcHi(1) prcHi(3)],[0.3 0.3],'k-','linewidth',2);
    plot([prcHi(2)],[0.3],'kd','markersize',10,'linewidth',2);
    xlabel('LC firing rate','fontsize',10);
    ylabel('proportion of trials','fontsize',10);
    % title('LC','fontsize',12);
    % if jjm == 1 title('LC spike counts'); end
    n_neurons = nOut(1)+nOut(2);
    text(7,0.05,strcat('# neurons=',num2str(n_neurons)),'fontsize',8);
    text(7,0.15,strcat('# trials =',num2str(NNN)),'fontsize',8);
    axis([-1 xMax+1 -0.01 0.36]);
    title('LC','fontsize',10);
    set(gca,'fontname','arial');
    
    %     pause
    
    %% ACC mean spike count for each LC spiking condition (0,1,2,3,>3).
    
    
    ms = [8 5 20 9 11];
    colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    barwidth=15;
    
    % Plot data.
    all_nz = [mean_b_nz{1};  mean_b_nz{2}];  ginz=find(isfinite(sum(all_nz,2))); mean_nz = nanmedian(all_nz(ginz,:));
    all_ze = [mean_b_ze{1};  mean_b_ze{2}];  gi0=find(isfinite(sum(all_ze,2))); mean_ze = nanmedian(all_ze(gi0,:));
    all_1 = [mean_b_1{1}; mean_b_1{2}]; gi1=find(isfinite(sum(all_1,2))); mean_1 = nanmedian(all_1(gi1,:));
    all_2 = [mean_b_2{1}; mean_b_2{2}]; gi2=find(isfinite(sum(all_2,2))); mean_2 = nanmedian(all_2(gi2,:));
    all_3 = [mean_b_3{1}; mean_b_3{2}]; gi3=find(isfinite(sum(all_3,2))); mean_3 = nanmedian(all_3(gi3,:));
    all_4 = [mean_b_4{1}; mean_b_4{2}]; gi4=find(isfinite(sum(all_4,2))); mean_4 = nanmedian(all_4(gi4,:));
    
    [ci_ze ci_ze_m] = bootci(2000,{@median,all_ze(gi0,:)},'alpha',0.05,'type','norm'); % mean_ze = nanmean(ci_ze_m);
    [ci_1 ci_1_m]= bootci(2000,{@median,all_1(gi1,:)},'alpha',0.05,'type','norm'); % mean_1 = nanmean(ci_1_m);
    [ci_2 ci_2_m]= bootci(2000,{@median,all_2(gi2,:)},'alpha',0.05,'type','norm'); % mean_2 = nanmean(ci_2_m);
    [ci_3 ci_3_m]= bootci(2000,{@median,all_3(gi3,:)},'alpha',0.05,'type','norm'); % mean_3 = nanmean(ci_3_m);
    [ci_4 ci_4_m]= bootci(2000,{@median,all_4(gi4,:)},'alpha',0.05,'type','norm'); % mean_4 = nanmean(ci_4_m);
    
    
    % *************
    % Set up ANOVA.
    % Number of data samples for each LC spike condition and each monkey.
    nze1 = length(mean_b_ze{1}); nze2 = length(mean_b_ze{2});
    n11 = length(mean_b_1{1}); n12 = length(mean_b_1{2});
    n21 = length(mean_b_2{1}); n22 = length(mean_b_2{2});
    n31 = length(mean_b_3{1}); n32 = length(mean_b_3{2});
    n41 = length(mean_b_4{1}); n42 = length(mean_b_4{2});
    
    % LC=0,1,2,3,>3 correspond to indices 1,2,3,4,5.
    LC_sp_ind_ze1 = 1*ones(nze1,1); LC_sp_ind_ze2 = 1*ones(nze2,1);
    LC_sp_ind_11 = 2*ones(n11,1); LC_sp_ind_12 = 2*ones(n12,1);
    LC_sp_ind_21 = 3*ones(n21,1); LC_sp_ind_22 = 3*ones(n22,1);
    LC_sp_ind_31 = 4*ones(n31,1); LC_sp_ind_32 = 4*ones(n32,1);
    LC_sp_ind_41 = 5*ones(n41,1); LC_sp_ind_42 = 5*ones(n42,1);
    
    dat1_LCspInd = [LC_sp_ind_ze1; LC_sp_ind_11; LC_sp_ind_21; LC_sp_ind_31; LC_sp_ind_41];
    dat2_LCspInd = [LC_sp_ind_ze2; LC_sp_ind_12; LC_sp_ind_22; LC_sp_ind_32; LC_sp_ind_42];
    
    % Per monkey ANOVA for testing LC spiking level contribution to rate comparison with LC zero.
    pvals = nans(length(xbs),2);
    p1A=[]; p2A=[]; pm1=[]; pm2=[];
    p1=[]; p2=[]; pnz=[];
    pvalsAll = nans(4,10);
    pvalsAllznz = nans(1,10);
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    for uj = 1:length(xbs)
        dat1 = [mean_b_ze{1}(:,uj); mean_b_1{1}(:,uj); mean_b_2{1}(:,uj); mean_b_3{1}(:,uj); mean_b_4{1}(:,uj)];
        dat2 = [mean_b_ze{2}(:,uj); mean_b_1{2}(:,uj); mean_b_2{2}(:,uj); mean_b_3{2}(:,uj); mean_b_4{2}(:,uj)];
        
        % 1-way ANOVA with display OFF.
        [p1A,tbl,stats] = anova1(dat1, dat1_LCspInd,'off');
        [p2A,tbl,stats] = anova1(dat2, dat2_LCspInd,'off');
        pvals(uj,:) = [p1A p2A];
        
        % pm1(uj) = ranksum(all_nz1(:,uj),mean_b_ze{1}(:,uj));
        % pm2(uj) = ranksum(all_nz2(:,uj),mean_b_ze{2}(:,uj));
        
        % Calculate sig level for comparison between LC zero and non zero cases.
        p1 = ranksum([mean_b_1{1}(:,uj);mean_b_1{2}(:,uj)],[mean_b_ze{1}(:,uj);mean_b_ze{2}(:,uj)]);
        p2 = ranksum([mean_b_2{1}(:,uj);mean_b_2{2}(:,uj)],[mean_b_ze{1}(:,uj);mean_b_ze{2}(:,uj)]);
        p3 = ranksum([mean_b_3{1}(:,uj);mean_b_3{2}(:,uj)],[mean_b_ze{1}(:,uj);mean_b_ze{2}(:,uj)]);
        p4 = ranksum([mean_b_4{1}(:,uj);mean_b_4{2}(:,uj)],[mean_b_ze{1}(:,uj);mean_b_ze{2}(:,uj)]);
        
        p1znz = ranksum(mean_b_ze{1}(:,uj),mean_b_nz{1}(:,uj));
        p2znz = ranksum(mean_b_ze{2}(:,uj),mean_b_nz{2}(:,uj));
        pznz = ranksum([mean_b_ze{1}(:,uj);mean_b_ze{2}(:,uj)],[mean_b_nz{1}(:,uj);mean_b_nz{2}(:,uj)]);
        
        % Plot.
        plot([xbs(uj)-20 xbs(uj)-20],ci_1(:,uj),'k-');
        plot([xbs(uj) xbs(uj)],ci_2(:,uj),'k-');
        plot([xbs(uj)+20 xbs(uj)+20],ci_3(:,uj),'k-');
        plot([xbs(uj)+20 xbs(uj)+40],ci_4(:,uj),'k-');
        
        plot(xbs(uj)-20,mean_1(uj),'.','color',colors{2},'markersize',25);
        plot(xbs(uj),mean_2(uj),'.','color',colors{3},'markersize',25);
        plot(xbs(uj)+20,mean_3(uj),'.','color',colors{4},'markersize',25);
        plot(xbs(uj)+40,mean_4(uj),'.','color',colors{5},'markersize',25);
        
        plot(xbs(uj)-40,mean_ze(uj),'.','color',colors{1},'markersize',25);
        plot([xbs(uj)-40 xbs(uj)-40],ci_ze(:,uj),'m-');
        
    end
    
    % Plots done, ANOVA done.
    % *************
    
    % Legend.
    yval = [0.9 1.2 1.5 1.8 2.1]; textVals = {'0' '1' '2' '3' '>3'};
    for ui = 1:5
        plot(150,yval(ui),'.','color',colors{ui},'markersize',ms(3));
        text(175,yval(ui),textVals(ui),'fontsize',8);
    end
    text(150,2.4,'LC spike','fontsize',8,'fontweight','bold');
    
    % Make plot sensible.
    axis([100 1100 -0.1 2.4]);
    % axis([100 1100 -0.1 3]);
    xlabel('binsize (ms)','fontsize',10);
    ylabel('ACC spike count','fontsize',10);
    set(gca,'fontname','arial');
    
    %     pause
    
    % *************
    
    %% Plot bars for change in spike count relative to LC zero.
    
    xjit = 30;
    
    all_dnz1 = mean_b_nz{1}-mean_b_ze{1};
    all_dnz2 = mean_b_nz{2}-mean_b_ze{2};
    all_dnz = [all_dnz1; all_dnz2];
    
    ci_1 = nans(10,2);  ci_2 = nans(10,2);  ci_3 = nans(10,2);  ci_4 = nans(10,2);
    mean_d1 = nans(1,10); mean_d2 = nans(1,10); mean_d3 = nans(1,10); mean_d4 = nans(1,10);
    
    for uj = 1:length(xbs)
        
        datT = all_dnz(:,uj); datT= datT(isfinite(datT));
        datT1 = all_dnz1(:,uj); datT1 = datT1(isfinite(datT1));
        datT2 = all_dnz2(:,uj); datT2 = datT2(isfinite(datT2));
        pvalznz(3,uj) = signrank(datT);
        pvalznz(1,uj) = signrank(all_dnz1(:,uj));
        pvalznz(2,uj) = signrank(all_dnz2(:,uj));
        
        [ci_1T ci_1_mT]= bootci(2000,{@median,datT},'alpha',0.05,'type','per'); mean_d1(uj) = nanmedian(datT);
        ci_NZ(uj,:) = ci_1T';
        
        [ci_1T1 ci_1_mT]= bootci(2000,{@median,datT2},'alpha',0.05,'type','per'); mean_d11(uj) = nanmedian(datT1);
        ci_NZ1(uj,:) = ci_1T1';
        
        [ci_1T2 ci_1_mT]= bootci(2000,{@median,datT2},'alpha',0.05,'type','per'); mean_d12(uj) = nanmedian(datT2);
        ci_NZ2(uj,:) = ci_1T2';
        
    end
    
    axes(axs(3)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    for uj = 1:length(xbs)
        
        %         plot([xbs(uj) xbs(uj)],ci_NZ(uj,:),'k-');
        %         bar(xbs(uj),mean_d1(uj),60,'facecolor',colors{2},'edgecolor','none');
        
        plot([xbs(uj)-xjit xbs(uj)-xjit],ci_NZ1(uj,:),'k-');
        if pvalznz(1,uj)<0.05
            bar(xbs(uj)-xjit,mean_d11(uj),2*(xjit-2),'facecolor',colors{2},'edgecolor','none','linewidth',1);
        end
        if pvalznz(1,uj)>=0.05
            bar(xbs(uj)-xjit,mean_d11(uj),2*(xjit-2),'edgecolor',colors{2},'facecolor','none','linewidth',1);
        end
        
        plot([xbs(uj)+xjit xbs(uj)+xjit],ci_NZ2(uj,:),'k-');
        if pvalznz(2,uj)<0.05
            bar(xbs(uj)+xjit,mean_d12(uj),2*(xjit-2),'facecolor',colors{3},'edgecolor','none','linewidth',1);
        end
        if pvalznz(2,uj)>=0.05
            bar(xbs(uj)+xjit,mean_d12(uj),2*(xjit-2),'edgecolor',colors{3},'facecolor','none','linewidth',1);
        end
        
        if pvalznz(3,uj)<0.05
            %         if pvalznz(1,uj)<0.05 & pvalznz(2,uj)<0.05
            plot(xbs(uj),0.05,'k*');
        end
        
        %         if pvalznz(1,uj)<0.05
        %             plot(xbs(uj)-xjit,0.035,'k*');
        %         end
        %
        %         if pvalznz(2,uj)<0.05
        %             plot(xbs(uj)+xjit,0.035,'k*');
        %         end
        
    end
    
    axis([100 1100 -0.052 0.052]);
    xlabel('binsize (ms)','fontsize',10);
    ylabel({'ACC spike count';'difference re: LC=0'},'fontsize',10);
    set(gca,'fontname','arial');
    
    %% ACC variance of spike count for each LC spiking condition (0,1,2,3,>3).
    
    % Plot data.
    all_ze = [var_b_ze{1};  var_b_ze{2}]; gi0=find(isfinite(sum(all_ze,2))); mean_ze = nanmedian(all_ze);
    all_1 = [var_b_1{1}; var_b_1{2}]; gi1=find(isfinite(sum(all_1,2))); mean_1 = nanmedian(all_1);
    all_2 = [var_b_2{1}; var_b_2{2}]; gi2=find(isfinite(sum(all_2,2))); mean_2 = nanmedian(all_2);
    all_3 = [var_b_3{1}; var_b_3{2}]; gi3=find(isfinite(sum(all_3,2))); mean_3 = nanmedian(all_3);
    all_4 = [var_b_4{1}; var_b_4{2}]; gi4=find(isfinite(sum(all_4,2))); mean_4 = nanmedian(all_4);
    
    [ci_ze ci_ze_m] = bootci(2000,{@median,all_ze(gi0,:)},'alpha',0.05,'type','per');
    [ci_1 ci_1_m]= bootci(2000,{@median,all_1(gi1,:)},'alpha',0.05,'type','per');
    [ci_2 ci_2_m]= bootci(2000,{@median,all_2(gi2,:)},'alpha',0.05,'type','per');
    [ci_3 ci_3_m]= bootci(2000,{@median,all_3(gi3,:)},'alpha',0.05,'type','per');
    [ci_4 ci_4_m]= bootci(2000,{@median,all_4(gi4,:)},'alpha',0.05,'type','per');
    
    % *************
    % Set up ANOVA.
    % Number of data samples for each LC spike condition and each monkey.
    nze1 = length(mean_b_ze{1}); nze2 = length(mean_b_ze{2});
    n11 = length(mean_b_1{1}); n12 = length(mean_b_1{2});
    n21 = length(mean_b_2{1}); n22 = length(mean_b_2{2});
    n31 = length(mean_b_3{1}); n32 = length(mean_b_3{2});
    n41 = length(mean_b_4{1}); n42 = length(mean_b_4{2});
    
    % LC=0,1,2,3,>3 correspond to indices 1,2,3,4,5.
    LC_sp_ind_ze1 = 1*ones(nze1,1); LC_sp_ind_ze2 = 1*ones(nze2,1);
    LC_sp_ind_11 = 2*ones(n11,1); LC_sp_ind_12 = 2*ones(n12,1);
    LC_sp_ind_21 = 3*ones(n21,1); LC_sp_ind_22 = 3*ones(n22,1);
    LC_sp_ind_31 = 4*ones(n31,1); LC_sp_ind_32 = 4*ones(n32,1);
    LC_sp_ind_41 = 5*ones(n41,1); LC_sp_ind_42 = 5*ones(n42,1);
    
    dat1_LCspInd = [LC_sp_ind_ze1; LC_sp_ind_11; LC_sp_ind_21; LC_sp_ind_31; LC_sp_ind_41];
    dat2_LCspInd = [LC_sp_ind_ze2; LC_sp_ind_12; LC_sp_ind_22; LC_sp_ind_32; LC_sp_ind_42];
    
    % Per monkey ANOVA for testing LC spiking level contribution to variance comparison with LC zero.
    pvals = nans(10,2);
    p1A=[]; p2A=[]; pm1=[]; pm2=[];
    p1=[]; p2=[]; pnz = [];
    pvalsAllV = nans(4,10);
    pvalsAllznz = nans(1,10);
    
    axes(axs(4)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    for uj = 1:length(xbs)
        dat1 = [var_b_ze{1}(:,uj); var_b_1{1}(:,uj); var_b_2{1}(:,uj); var_b_3{1}(:,uj); var_b_4{1}(:,uj)];
        dat2 = [var_b_ze{2}(:,uj); var_b_1{2}(:,uj); var_b_2{2}(:,uj); var_b_3{2}(:,uj); var_b_4{2}(:,uj)];
        
        % 1-way ANOVA with display OFF.
        [p1A,tbl,stats] = anova1(dat1, dat1_LCspInd,'off');
        [p2A,tbl,stats] = anova1(dat2, dat2_LCspInd,'off');
        pvals(uj,:) = [p1A p2A];
        
        %         pm1(uj) = ranksum(all_nz1(:,uj),var_b_ze{1}(:,uj));
        %         pm2(uj) = ranksum(all_nz2(:,uj),var_b_ze{2}(:,uj));
        
        % Calculate sig level for comparison between LC zero and non zero cases.
        p11 = ranksum(var_b_1{1}(:,uj),var_b_ze{1}(:,uj));
        p12 = ranksum(var_b_1{2}(:,uj),var_b_ze{2}(:,uj));
        p21 = ranksum(var_b_2{1}(:,uj),var_b_ze{1}(:,uj));
        p22 = ranksum(var_b_2{2}(:,uj),var_b_ze{2}(:,uj));
        p31 = ranksum(var_b_3{1}(:,uj),var_b_ze{1}(:,uj));
        p32 = ranksum(var_b_3{2}(:,uj),var_b_ze{2}(:,uj));
        p41 = ranksum(var_b_4{1}(:,uj),var_b_ze{1}(:,uj));
        p42 = ranksum(var_b_4{2}(:,uj),var_b_ze{2}(:,uj));
        
        p1znz = ranksum(var_b_ze{1}(:,uj),var_b_hi{1}(:,uj));
        p2znz = ranksum(var_b_ze{2}(:,uj),var_b_hi{2}(:,uj));
        
        %         if p1znz<0.05 & p2znz<0.05 plot(xbs(uj),10,'k+'); end
        %
        %         % Plot * if sig.
        %         if p11<0.05 & p12<0.05 plot(xbs(uj)-15,8,'k*'); end
        %         if p21<0.05 & p22<0.05 plot(xbs(uj),8,'k*'); end
        %         if p31<0.05 & p32<0.05 plot(xbs(uj)+15,8,'k*'); end
        %         if p41<0.05 & p42<0.05 plot(xbs(uj)+30,8,'k*'); end
        
        % Plot.
        plot([xbs(uj)-20 xbs(uj)-20],ci_1(:,uj),'k-');
        plot([xbs(uj) xbs(uj)],ci_2(:,uj),'k-');
        plot([xbs(uj)+20 xbs(uj)+20],ci_3(:,uj),'k-');
        plot([xbs(uj)+40 xbs(uj)+40],ci_4(:,uj),'k-');
        
        plot(xbs(uj)-20,mean_1(uj),'.','color',colors{2},'markersize',25);
        plot(xbs(uj),mean_2(uj),'.','color',colors{3},'markersize',25);
        plot(xbs(uj)+20,mean_3(uj),'.','color',colors{4},'markersize',25);
        plot(xbs(uj)+40,mean_4(uj),'.','color',colors{5},'markersize',25);
        
        plot(xbs(uj)-40,mean_ze(uj),'.','color',colors{1},'markersize',25);
        plot([xbs(uj)-40 xbs(uj)-40],ci_ze(:,uj),'m-');
        
    end
    
    % Plots done, ANOVA done.
    % *************
    
    % Make plot sensible.
    axis([100 1100 -0.2 4.5]);
    xlabel('binsize(ms)','fontsize',10);
    ylabel('ACC variance of spike count','fontsize',10);
    set(gca,'fontname','arial');
    
    %% Plot bars for change in spike count variance relative to LC zero.
    
    all_dnz1 = var_b_hi{1}-var_b_ze{1};
    all_dnz2 = var_b_hi{2}-var_b_ze{2};
    all_dnz = [all_dnz1; all_dnz2];
    
    ci_NZ = nans(10,2);
    ci_1 = nans(10,2);  ci_2 = nans(10,2);  ci_3 = nans(10,2);  ci_4 = nans(10,2);
    mean_d1 = nans(1,10); mean_d2 = nans(1,10); mean_d3 = nans(1,10); mean_d4 = nans(1,10);
    pval = [];
    
    
    for uj = 1:length(xbs)
        
        datT = all_dnz(:,uj); datT= datT(isfinite(datT));
        datT1 = all_dnz1(:,uj); datT1 = datT1(isfinite(datT1));
        datT2 = all_dnz2(:,uj); datT2 = datT2(isfinite(datT2));
        pvalznz(3,uj) = signrank(datT);
        pvalznz(1,uj) = signrank(datT1);
        pvalznz(2,uj) = signrank(datT2);
        
        [ci_1T ci_1_mT]= bootci(2000,{@median,datT},'alpha',0.05,'type','per'); mean_d1(uj) = nanmedian(datT);
        ci_NZ(uj,:) = ci_1T';
        
        [ci_1T1 ci_1_mT]= bootci(2000,{@median,datT2},'alpha',0.05,'type','per'); mean_d11(uj) = nanmedian(datT1);
        ci_NZ1(uj,:) = ci_1T1';
        
        [ci_1T2 ci_1_mT]= bootci(2000,{@median,datT2},'alpha',0.05,'type','per'); mean_d12(uj) = nanmedian(datT2);
        ci_NZ2(uj,:) = ci_1T2';
        
    end
    
    axes(axs(5)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    for uj = 1:length(xbs)
        
        %         plot([xbs(uj) xbs(uj)],ci_NZ(uj,:),'k-');
        %         bar(xbs(uj),mean_d1(uj),60,'facecolor',colors{2},'edgecolor','none');
        
        plot([xbs(uj)-xjit xbs(uj)-xjit],ci_NZ1(uj,:),'k-');
        if pvalznz(1,uj)<0.05
            bar(xbs(uj)-xjit,mean_d11(uj),2*(xjit-2),'facecolor',colors{2},'edgecolor','none','linewidth',1);
        end
        if pvalznz(1,uj)>=0.05
            bar(xbs(uj)-xjit,mean_d11(uj),2*(xjit-2),'edgecolor',colors{2},'facecolor','none','linewidth',1);
        end
        
        plot([xbs(uj)+xjit xbs(uj)+xjit],ci_NZ2(uj,:),'k-');
        if pvalznz(2,uj)<0.05
            bar(xbs(uj)+xjit,mean_d12(uj),2*(xjit-2),'facecolor',colors{3},'edgecolor','none','linewidth',1);
        end
        if pvalznz(2,uj)>=0.05
            bar(xbs(uj)+xjit,mean_d12(uj),2*(xjit-2),'edgecolor',colors{3},'facecolor','none','linewidth',1);
        end
        
        if pvalznz(3,uj)<0.05
            %         if pvalznz(1,uj)<0.05 & pvalznz(2,uj)<0.05
            plot(xbs(uj),0.1,'k*');
        end
        
        %         if pvalznz(1,uj)<0.05
        %             plot(xbs(uj)-xjit,0.011,'k*');
        %         end
        %
        %         if pvalznz(2,uj)<0.05
        %             plot(xbs(uj)+xjit,0.011,'k*');
        %         end
    end
    
    axis([100 1100 -0.12 0.12]);
    xlabel('binsize (ms)','fontsize',10);
    ylabel({'ACC variance of spike count';'difference re: LC=0'},'fontsize',10);
    
    
    %% ACC Fano factor for each LC spiking condition (0,1,2,3,>3).
    
    %     pause
    
    % Plot data.
    all_ze = [FF_b_ze{1};  FF_b_ze{2}]; mean_ze = nanmedian(all_ze); gi0=find(isfinite(sum(all_ze,2)));
    all_1 = [FF_b_1{1}; FF_b_1{2}]; mean_1 = nanmedian(all_1); gi1=find(isfinite(sum(all_1,2)));
    all_2 = [FF_b_2{1}; FF_b_2{2}]; mean_2 = nanmedian(all_2); gi2=find(isfinite(sum(all_2,2)));
    all_3 = [FF_b_3{1}; FF_b_3{2}]; mean_3 = nanmedian(all_3); gi3=find(isfinite(sum(all_3,2)));
    all_4 = [FF_b_4{1}; FF_b_4{2}]; mean_4 = nanmedian(all_4); gi4=find(isfinite(sum(all_4,2)));
    
    all_all_FF = [all_ze;all_1;all_2;all_3;all_4];
    all_all_FF = all_all_FF(:);
    
    mean_ze_ALL = [mean_ze;mean_1;mean_2;mean_3;mean_4];
    mean_ze_ALL = mean_ze_ALL(:);
    
    [ci_ze ci_ze_m] = bootci(2000,{@median,all_ze(gi0,:)},'alpha',0.05,'type','per');
    [ci_1 ci_1_m]= bootci(2000,{@median,all_1(gi1,:)},'alpha',0.05,'type','per');
    [ci_2 ci_2_m]= bootci(2000,{@median,all_2(gi2,:)},'alpha',0.05,'type','per');
    [ci_3 ci_3_m]= bootci(2000,{@median,all_3(gi3,:)},'alpha',0.05,'type','per');
    [ci_4 ci_4_m]= bootci(2000,{@median,all_4(gi4,:)},'alpha',0.05,'type','per');
    
    % *************
    % Set up ANOVA.
    % Number of data samples for each LC spike condition and each monkey.
    nze1 = length(FF_b_ze{1}); nze2 = length(FF_b_ze{2});
    n11 = length(FF_b_1{1}); n12 = length(FF_b_1{2});
    n21 = length(FF_b_2{1}); n22 = length(FF_b_2{2});
    n31 = length(FF_b_3{1}); n32 = length(FF_b_3{2});
    n41 = length(FF_b_4{1}); n42 = length(FF_b_4{2});
    
    % LC=0,1,2,3,>3 correspond to indices 1,2,3,4,5.
    LC_sp_ind_ze1 = 1*ones(nze1,1); LC_sp_ind_ze2 = 1*ones(nze2,1);
    LC_sp_ind_11 = 2*ones(n11,1); LC_sp_ind_12 = 2*ones(n12,1);
    LC_sp_ind_21 = 3*ones(n21,1); LC_sp_ind_22 = 3*ones(n22,1);
    LC_sp_ind_31 = 4*ones(n31,1); LC_sp_ind_32 = 4*ones(n32,1);
    LC_sp_ind_41 = 5*ones(n41,1); LC_sp_ind_42 = 5*ones(n42,1);
    
    dat1_LCspInd = [LC_sp_ind_ze1; LC_sp_ind_11; LC_sp_ind_21; LC_sp_ind_31; LC_sp_ind_41];
    dat2_LCspInd = [LC_sp_ind_ze2; LC_sp_ind_12; LC_sp_ind_22; LC_sp_ind_32; LC_sp_ind_42];
    
    % Per monkey ANOVA for testing LC spiking level contribution to variance comparison with LC zero.
    pvalsFF = nans(10,2);
    p1A=[]; p2A=[]; pm1=[]; pm2=[];
    p1=[]; p2=[]; pnz=[];
    pvalsAllF = nans(4,10);
    pvalsAllznz = nans(1,10);
    
    axes(axs(6)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    for uj = 1:length(xbs)
        dat1 = [FF_b_ze{1}(:,uj); FF_b_1{1}(:,uj); FF_b_2{1}(:,uj); FF_b_3{1}(:,uj); FF_b_4{1}(:,uj)];
        dat2 = [FF_b_ze{2}(:,uj); FF_b_1{2}(:,uj); FF_b_2{2}(:,uj); FF_b_3{2}(:,uj); FF_b_4{2}(:,uj)];
        
        % 1-way ANOVA with display OFF.
        [p1A,tbl,stats] = anova1(dat1, dat1_LCspInd,'off');
        [p2A,tbl,stats] = anova1(dat2, dat2_LCspInd,'off');
        pvalsFF(uj,:) = [p1A p2A];
        
        %         pm1(uj) = ranksum(all_nz1(:,uj),var_b_ze{1}(:,uj));
        %         pm2(uj) = ranksum(all_nz2(:,uj),var_b_ze{2}(:,uj));
        
        % Calculate sig level for comparison between LC zero and non zero cases.
        p11 = ranksum(FF_b_1{1}(:,uj),FF_b_ze{1}(:,uj));
        p12 = ranksum(FF_b_1{2}(:,uj),FF_b_ze{2}(:,uj));
        p21 = ranksum(FF_b_2{1}(:,uj),FF_b_ze{1}(:,uj));
        p22 = ranksum(FF_b_2{2}(:,uj),FF_b_ze{2}(:,uj));
        p31 = ranksum(FF_b_3{1}(:,uj),FF_b_ze{1}(:,uj));
        p32 = ranksum(FF_b_3{2}(:,uj),FF_b_ze{2}(:,uj));
        p41 = ranksum(FF_b_4{1}(:,uj),FF_b_ze{1}(:,uj));
        p42 = ranksum(FF_b_4{2}(:,uj),FF_b_ze{2}(:,uj));
        
        p1znz = ranksum(FF_b_ze{1}(:,uj),FF_b_nz{1}(:,uj));
        p2znz = ranksum(FF_b_ze{2}(:,uj),FF_b_nz{2}(:,uj));
        
        %         if p1znz<0.05 & p2znz<0.05 pvalsAllznz(1,uj) = 1; end
        %         if p1znz<0.05 & p2znz<0.05 plot(xbs(uj),2.5,'k+'); end
        
        %         if p11<0.05 & p12<0.05 pvalsAllF(1,uj) = 1; end
        %         if p21<0.05 & p22<0.05 pvalsAllF(2,uj) = 1; end
        %         if p31<0.05 & p32<0.05 pvalsAllF(3,uj) = 1; end
        %         if p41<0.05 & p42<0.05 pvalsAllF(4,uj) = 1; end
        
        % Plot * if sig.
        %                 if p11<0.05 & p12<0.05 plot(xbs(uj)-15,2.3,'k*'); end
        %                 if p21<0.05 & p22<0.05 plot(xbs(uj),2.3,'k*'); end
        %                 if p31<0.05 & p32<0.05 plot(xbs(uj)+15,2.3,'k*'); end
        %                 if p41<0.05 & p42<0.05 plot(xbs(uj)+30,2.3,'k*'); end
        
        % Plot.
        plot([xbs(uj)-20 xbs(uj)-20],ci_1(:,uj),'k-');
        plot([xbs(uj) xbs(uj)],ci_2(:,uj),'k-');
        plot([xbs(uj)+20 xbs(uj)+20],ci_3(:,uj),'k-');
        plot([xbs(uj)+40 xbs(uj)+40],ci_4(:,uj),'k-');
        
        plot(xbs(uj)-20,mean_1(uj),'.','color',colors{2},'markersize',25);
        plot(xbs(uj),mean_2(uj),'.','color',colors{3},'markersize',25);
        plot(xbs(uj)+20,mean_3(uj),'.','color',colors{4},'markersize',25);
        plot(xbs(uj)+40,mean_4(uj),'.','color',colors{5},'markersize',25);
        
        plot(xbs(uj)-40,mean_ze(uj),'.','color',colors{1},'markersize',25);
        plot([xbs(uj)-40 xbs(uj)-40],ci_ze(:,uj),'m-');
        
    end
    
    % Plots done, ANOVA done.
    % *************
    
    % Make plot sensible.
    axis([100 1100 0.9 2]);
    xlabel('binsize(ms)','fontsize',10);
    ylabel('ACC Fano factor','fontsize',10);
    set(gca,'fontname','arial');
    
    % pause
    
    %% Plot bars for change in FF relative to LC zero.
    
    all_dnz1 = FF_b_nz{1}-FF_b_ze{1};
    all_dnz2 = FF_b_nz{2}-FF_b_ze{2};
    all_dnz = [all_dnz1; all_dnz2];
    
    
    ci_NZ = nans(10,2);
    ci_1 = nans(10,2);  ci_2 = nans(10,2);  ci_3 = nans(10,2);  ci_4 = nans(10,2);
    mean_d1 = nans(1,10); mean_d2 = nans(1,10); mean_d3 = nans(1,10); mean_d4 = nans(1,10);
    pval = [];
    
    %     pause
    
    for uj = 1:length(xbs)
        
        datT = all_dnz(:,uj); datT= datT(isfinite(datT));
        datT1 = all_dnz1(:,uj); datT1 = datT1(isfinite(datT1));
        datT2 = all_dnz2(:,uj); datT2 = datT2(isfinite(datT2));
        pvalznz(3,uj) = signrank(datT);
        pvalznz(1,uj) = signrank(datT1);
        pvalznz(2,uj) = signrank(datT2);
        
        [ci_1T ci_1_mT]= bootci(2000,{@median,datT},'alpha',0.05,'type','per'); mean_d1(uj) = nanmedian(datT);
        ci_NZ(uj,:) = ci_1T';
        
        [ci_1T1 ci_1_mT]= bootci(2000,{@median,datT2},'alpha',0.05,'type','per'); mean_d11(uj) = nanmedian(datT1);
        ci_NZ1(uj,:) = ci_1T1';
        
        [ci_1T2 ci_1_mT]= bootci(2000,{@median,datT2},'alpha',0.05,'type','per'); mean_d12(uj) = nanmedian(datT2);
        ci_NZ2(uj,:) = ci_1T2';
        
    end
    
    axes(axs(7)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    for uj = 1:length(xbs)
        
        %         plot([xbs(uj) xbs(uj)],ci_NZ(uj,:),'k-');
        %         bar(xbs(uj),mean_d1(uj),60,'facecolor',colors{2},'edgecolor','none');
        
        plot([xbs(uj)-xjit xbs(uj)-xjit],ci_NZ1(uj,:),'k-');
        if pvalznz(1,uj)<0.05
            bar(xbs(uj)-xjit,mean_d11(uj),2*(xjit-2),'facecolor',colors{2},'edgecolor','none','linewidth',1);
        end
        if pvalznz(1,uj)>=0.05
            bar(xbs(uj)-xjit,mean_d11(uj),2*(xjit-2),'facecolor',colors{2},'edgecolor','none','linewidth',1);
        end
        
        plot([xbs(uj)+xjit xbs(uj)+xjit],ci_NZ2(uj,:),'k-');
        if pvalznz(2,uj)<0.05
            bar(xbs(uj)+xjit,mean_d12(uj),2*xjit,'facecolor',colors{3},'edgecolor','none','linewidth',1);
        end
        if pvalznz(2,uj)>=0.05
            bar(xbs(uj)+xjit,mean_d12(uj),2*xjit,'edgecolor',colors{3},'facecolor','none','linewidth',1);
        end
        
        if pvalznz(3,uj)<0.05
            %         if pvalznz(1,uj)<0.05 & pvalznz(2,uj)<0.05
            plot(xbs(uj),0.1,'k*');
        end
        
        %         if pvalznz(1,uj)<0.05
        %             plot(xbs(uj)-xjit,0.011,'k*');
        %         end
        %
        %         if pvalznz(2,uj)<0.05
        %             plot(xbs(uj)+xjit,0.011,'k*');
        %         end
        
    end
    
    axis([100 1110 -0.11 0.11]);
    xlabel('binsize(ms)','fontsize',10);
    ylabel({'ACC FF';'difference re: LC=0'},'fontsize',10);
    set(gca,'fontname','arial');
    
    % ********************************************
    % Figure 2 ends here.
    % ********************************************
    
    %         pause
    
    %% Supplement: Fano and firing rate.
    
    % figureNumber = 223; num = 223; wid = 17.6; hts = [8]; cols = {2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,3,1, [12], 'Joshi & Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
    figureNumber = 21; num = 21; wid = 17.6; hts = [12]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,2,2, [12], 'Joshi & Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    
    %%
    
    
    all_ze = [mean_b_ze{1}(:,nbs); mean_b_ze{2}(:,nbs)];
    all_nz = [mean_b_nz{1}(:,nbs); mean_b_nz{2}(:,nbs)];
    
    all_zeF = [FF_b_ze{1}(:,nbs); FF_b_ze{2}(:,nbs)];
    all_nzF = [FF_b_nz{1}(:,nbs); FF_b_nz{2}(:,nbs)];
    
    all_diff_Rate = (all_nz-all_ze)./all_ze;
    all_diff_FF = (all_nzF-all_zeF)./all_zeF;
    
    %     gi = find(isfinite(all_diff_FF) & isfinite(all_diff_Rate) & abs(all_diff_FF)<2 & abs(all_diff_Rate)<2);
    gi = find(isfinite(all_diff_FF) & isfinite(all_diff_Rate));
    
    mF = nanmedian(all_diff_FF(gi));
    mR = nanmedian(all_diff_Rate(gi));
    
    %     [b,bint,r,rint,stats] = regress(all_diff_FF(gi),[ones(length(all_diff_Rate(gi)),1) all_diff_Rate(gi)]);
    %     xx = linspace(min(all_diff_Rate(gi)),max(all_diff_Rate(gi)),10); yy = b(1)+b(2)*xx;
    
    xb = linspace(-1,1,30); yb = xb;
    xDat = all_diff_Rate(gi);
    yDat = all_diff_FF(gi);
    nx = hist(xDat,xb); mx = nanmedian(xDat);
    ny = hist(yDat,yb); my = nanmedian(yDat);
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    % figure; hold on;
    
    % histogram2(xDat,yDat,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');
    histogram2(xDat,yDat,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgealpha',0,'edgecolor','none');
    colormap(gray);
    
    [r, p] = corr(xDat,yDat,'type','Spearman'); % 0.27; p=2.6e-17
    
    %     pause
    
    % plot(xDat,yDat,'k.','markersize',10,'color',[.5 .5 .5]);
    % plot(mF,mR,'r.','markersize',10);
    % plot([-1 1],[-1 1],'g:');
    plot([-1 1],[0 0],'w:');
    plot([0 0],[-1 1],'w:');
    % plot(xx,yy,'r-');
    title({'Relation between changes in';'ACC neuron firing rate and Fano factor'},'fontsize',10);
    xlabel('\Delta spike count','fontsize',10);
    ylabel('\Delta Fano factor','fontsize',10);
    cbh = colorbar; set(cbh,'box','off');
    ylabel(cbh,'number of neurons','fontsize',10);
    axis([-1 1 -1 1]);
    axis square;
    set(gca,'fontname','arial');
    
    %             pause
    
    %% Supplement: Fano and variance of spike count.
    
    all_ze = [var_b_ze{1}(:,nbs); var_b_ze{2}(:,nbs)];
    all_nz = [var_b_hi{1}(:,nbs); var_b_hi{2}(:,nbs)];
    
    all_zeF = [FF_b_ze{1}(:,nbs); FF_b_ze{2}(:,nbs)];
    all_nzF = [FF_b_nz{1}(:,nbs); FF_b_nz{2}(:,nbs)];
    
    all_diff_Rate = (all_nz-all_ze)./all_ze;
    all_diff_FF = (all_nzF-all_zeF)./all_zeF;
    
    %     gi = find(isfinite(all_diff_FF) & isfinite(all_diff_Rate) & abs(all_diff_FF)<2 & abs(all_diff_Rate)<2);
    gi = find(isfinite(all_diff_FF) & isfinite(all_diff_Rate));
    
    mF = nanmedian(all_diff_FF(gi));
    mR = nanmedian(all_diff_Rate(gi));
    
    %     [b,bint,r,rint,stats] = regress(all_diff_FF(gi),[ones(length(all_diff_Rate(gi)),1) all_diff_Rate(gi)]);
    %     xx = linspace(min(all_diff_Rate(gi)),max(all_diff_Rate(gi)),10); yy = b(1)+b(2)*xx;
    
    xb = linspace(-1,1,30); yb = xb;
    xDat = all_diff_Rate(gi);
    yDat = all_diff_FF(gi);
    nx = hist(xDat,xb); mx = nanmedian(xDat);
    ny = hist(yDat,yb); my = nanmedian(yDat);
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    % figure; hold on;
    
    % histogram2(xDat,yDat,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');
    histogram2(xDat,yDat,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgealpha',0,'edgecolor','auto');
    colormap(gray);
    
    [r, p] = corr(xDat,yDat,'type','Spearman'); % r=0.85; p=0.
    
    % plot(xDat,yDat,'k.','markersize',10,'color',[.5 .5 .5]);
    % plot(mF,mR,'r.','markersize',10);
    % plot([-1 1],[-1 1],'g:');
    plot([-1 1],[0 0],'w:');
    plot([0 0],[-1 1],'w:');
    % plot(xx,yy,'r-');
    title({'Relation between changes in';'ACC neuron spike count variance and Fano factor'},'fontsize',10);
    xlabel('\Delta variance of spike count','fontsize',10);
    ylabel('\Delta Fano factor','fontsize',10);
    cbh = colorbar; set(cbh,'box','off');
    ylabel(cbh,'number of neurons','fontsize',10);
    axis([-1 1 -1 1]);
    axis square;
    set(gca,'fontname','arial');
    
end % plotYesNo

% pause

%% ACC linked LC results.

% Don't organize results?
orgYesNo = false;
% Plot results?
orgYesNo = true;

if orgYesNo
    
    clear; clear all;
    
    base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2020\LC_ACC'; % Base directory for brain area.
    savedir = base_dir;
    % inFile = 'ACC_LC_Results_rsc_120219'; % New way of dividing ACC trials into low, mid, high.
    % inFile = 'ACC_LC_Results_rsc_040320'; % New way of dividing ACC trials into low, mid, high.
    % inFile = 'ACC_LC_Results_rsc_042420'; % New way of dividing ACC trials into low, mid, high.
    
    % inFile = 'ACC_LC_Results_rsc_081620'; % New way of dividing ACC trials into low, mid, high.
    
    % inFile = 'ACC_LC_Results_rsc_100120'; % New way of dividing ACC trials into low, mid, high.
    inFile = 'ACC_LC_Results_rsc_051621'; % New way of dividing ACC trials into low, mid, high.
    
    cd(savedir); load(inFile);
    
    monks = {'Sprout', 'Cicero'}; % Add mnks as needed...
    sites = {'LC_ACC_Fixation'}; % For now LC only.
    nMonks = length(monks); nSites = length(sites);
    
    % **********************************************
    unit_all_lo = []; unit_all_hi = []; unit_all_all = [];
    
    ze_pair_small_1 = []; ze_pair_large_1 = []; ze_pair_small_2 = []; ze_pair_large_2 = [];
    hi_pair_small_1 = []; hi_pair_large_1 = []; hi_pair_small_2 = []; hi_pair_large_2= [];
    
    nMinTrials = 5;
    for mm = 1:nMonks % Loop over mnks....
        
        all_lo_c_meanT = []; all_mid_c_meanT = []; all_hi_c_meanT = []; all_all_c_meanT = []; % Spike count matrix.
        var_loT = []; var_midT = []; var_hiT = []; var_ALLT = []; % Variance of spike counts.
        
        unit_all_loT = []; unit_all_midT = []; unit_all_hiT = []; unit_all_allT = []; % Pairwise spike count correlation, rsc.
        unit_all_lo_pvT = []; unit_all_mid_pvT = []; unit_all_hi_pvT = []; unit_all_all_pvT = []; % p-value from crosscorr calculation of rsc.
        %         unit_all_lo_pT = []; unit_all_hi_pT = []; unit_all_all_pT = []; % rsc values from shuffled trials.
        unit_all_lo_pT = cell(1,10); unit_all_hi_pT = cell(1,10); unit_all_all_pT = cell(1,10); % rsc values from shuffled trials.
        unit_all_rate_lo = []; unit_all_rate_mid = []; unit_all_rate_hi = [];
        unit_allZi = []; unit_allHi = [];
        acc_bc_loT = []; acc_bc_midT = []; acc_bc_hiT = []; % Mean binned spike count for each ACC unit.
        unit_dat_rate_loT = []; unit_dat_rate_hiT = []; unit_dat_rate_alT = []; % ACC neuron spike rates.
        
        FF_b_loT = []; FF_b_midT = []; FF_b_hiT = [];
        var_b_loT = []; var_b_midT = []; var_b_hiT = [];
        mean_b_loT = []; mean_b_midT = []; mean_b_hiT = [];
        %     all_lo_c_medianT = []; all_lo_c_medianT = []; all_hi_c_medianT = []; all_all_c_medianT = [];
        %     std_loT = []; std_loT = []; std_hiT = []; std_allT = []; %  std_mid{2} T = [];
        all_lo_count1T =[]; all_lo_count1T =[]; all_hi_count1T =[]; all_lo_count2T =[]; all_lo_count2T =[]; all_hi_count2T =[];
        
        mnkDat = mnkRes{mm};
        
        lci = 0;
        lc_u_lo = [];
        lc_u_hi = [];
        
        for ss = 1:nSites % Loop over brainstem sites.... LC and/or IC.
            
            sitDat = mnkDat{ss};
            nnuu = 0; nnAA = 0; % Housekeeping - how many LC and ACC units did we have?
            
            if ~isempty(sitDat)
                for ap = 1:length(sitDat) % Loop over sessions.
                    dat = sitDat{ap}; % Get session data.
                    if ~isempty(dat)
                        allZ = []; allM = []; allH = []; allA = [];
                        % acc_count = dat{3};
                        % nnuu = nnuu + length(dat{1});
                        nnAA = nnAA + length(dat{2});
                        
                        % *********************************************
                        % Get pair indices.
                        pair_ind = dat{6};
                        % pause
                        % for sp = 1:length(dat{1}) % Loop over LC units in this session (can be 1 or more).
                        
                        LC_lohi = dat{4};
                        
                        loi = find(LC_lohi==0);
                        %                         midi = find(LC_lohi==1);
                        %                         hii = find(LC_lohi==2);
                        hii = find(LC_lohi==1);
                        
                        if ~isempty(loi) && ~isempty(hii)
                            
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
                            
                            all_lo = []; all_mid = []; all_hi = []; all_all = [];
                            all_lo_pv = []; all_mid_pv = []; all_hi_pv = []; all_all_pv = [];
                            all_lo_p = []; all_mid_p = []; all_hi_p = []; all_all_p = [];
                            dat_rate_lo = []; dat_rate_mid = []; dat_rate_hi = []; dat_rate_al = [];
                            % all_lo_rate = []; all_mid_rate = [];  all_hi_rate = [];
                            
                            % dat0 = dat{1}{sp};
                            dat0 = dat{1}{1};
                            
                            dat1_LC_spikes = dat{5}; % LC trial spike counts.
                            LC_sp_lo = (dat1_LC_spikes(LC_lohi==0));
                            %                             LC_sp_mid = (dat1_LC_spikes(LC_lohi==1));
                            % LC_sp_hi = (dat1_LC_spikes(LC_lohi==2));
                            LC_sp_hi = (dat1_LC_spikes(LC_lohi==1));
                            
                            for bhi= 1:length(dat0) % Iterate over bins.
                                dat1 = dat0{bhi}{1}; % Pairwise spike-count correlation results.
                                dat1_plo = dat0{bhi}{3}; % Zero LC. Pairwise spike-count correlation results - shuffle tests.
                                dat1_phi = dat0{bhi}{4}; % High LC. Pairwise spike-count correlation results - shuffle tests.
                                dat1_pall = dat0{bhi}{5}; % All trials. Pairwise spike-count correlation results - shuffle tests.
                                
                                % For largest binsize, get the ACC pair spike rate.
                                if bhi == 5
                                    dat_rate = dat0{bhi}{2}; % Get mean spike rate of each neuron in each pair.
                                    for pp = 1:length(dat_rate) % Iterate over ACC pairs.
                                        s_rate = dat_rate{pp};
                                        dat_rate_lo = [dat_rate_lo; s_rate(1) s_rate(2)];
                                        dat_rate_hi = [dat_rate_hi; s_rate(3) s_rate(4)];
                                        dat_rate_al = [dat_rate_al; s_rate(5) s_rate(6)];
                                    end
                                end
                                
                                % Organize the pairwise corr data:
                                
                                all_lo = [all_lo dat1(:,1)];
                                % all_lo_p = [all_lo_p dat1_plo'];
                                all_lo_pv = [all_lo_pv dat1(:,4)];
                                
                                %                                 all_mid = [all_mid dat1(:,2)];
                                %                                 % all_mid_p = [all_mid_p dat1_pmid'];
                                %                                 all_mid_pv = [all_mid_pv dat1(:,6)];
                                
                                all_hi = [all_hi dat1(:,2)];
                                % all_hi_p = [all_hi_p dat1_phi'];
                                all_hi_pv = [all_hi_pv dat1(:,5)];
                                
                                all_all = [all_all dat1(:,3)];
                                % all_all_p = [all_all_p dat1_pall'];
                                all_all_pv = [all_all_pv dat1(:,6)];
                                
                                % Shuffled trials.
                                unit_all_lo_pT{bhi} = [unit_all_lo_pT{bhi}; dat1_plo];
                                unit_all_hi_pT{bhi} = [unit_all_hi_pT{bhi}; dat1_phi];
                                unit_all_all_pT{bhi} = [unit_all_all_pT{bhi}; dat1_pall];
                                
                            end % Loop over number of bins for rsc data.
                            
                            
                            % Collect the ACC pair spike rates for this LC unit.
                            if bhi == 10
                                unit_dat_rate_loT = [unit_dat_rate_loT; dat_rate_lo];
                                unit_dat_rate_hiT = [unit_dat_rate_hiT; dat_rate_hi];
                                unit_dat_rate_alT = [unit_dat_rate_alT; dat_rate_al];
                            end
                            
                            % prev_pnum = size(unit_all_loT,1);
                            % Collect rsc for this binsize, this LC unit.
                            unit_all_loT = [unit_all_loT; all_lo];   unit_all_lo_pvT = [unit_all_lo_pvT; all_lo_pv];
                            %                             unit_all_midT = [unit_all_midT; all_mid];   unit_all_mid_pvT = [unit_all_mid_pvT; all_mid_pv];
                            unit_all_hiT = [unit_all_hiT; all_hi];   unit_all_hi_pvT = [unit_all_hi_pvT; all_hi_pv];
                            unit_all_allT = [unit_all_allT; all_all]; unit_all_all_pvT = [unit_all_all_pvT; all_all_pv];
                            
                            % this_pnum = size(unit_all_loT,1);
                            
                            
                            unit_all_rate_lo = [unit_all_rate_lo LC_sp_lo];
                            %                             unit_all_rate_mid = [unit_all_rate_mid LC_sp_mid];
                            unit_all_rate_hi = [unit_all_rate_hi LC_sp_hi];
                            
                            unit_allZi = [unit_allZi loi];
                            unit_allHi = [unit_allHi hii];
                            
                            % *********************************************
                            
                            % Now get binned LC spike count data to calculate FF.
                            
                            LC_binned = dat{3}; % Binned LC unit spike counts.
                            % Loop over ACC units for spike counts and calculate FF.
                            for ac = 1:size(LC_binned,1) % Loop over ACC units to get binned spike counts conditioned on this LC unit.
                                
                                acu_bd = LC_binned(ac,:);
                                FF_b_hi_temp = []; FF_b_mid_temp = []; FF_b_lo_temp = [];
                                mean_b_hi_temp = []; mean_b_mid_temp = []; mean_b_lo_temp = [];
                                var_b_hi_temp = []; var_b_mid_temp = []; var_b_lo_temp = [];
                                
                                acc_uc_lo = []; acc_uc_mid = []; acc_uc_hi = [];
                                
                                for acb = 1:length(acu_bd) % Loop over binsizes.
                                    
                                    % ************************************
                                    % Low ACC spiking case.
                                    acu_bd_lo = acu_bd{acb}(loi,:); % pause
                                    % mbc_v = reshape(acu_bd_lo,size(acu_bd_lo,1)*size(acu_bd_lo,2),1);
                                    mbc_v = acu_bd_lo(:);
                                    mbc_mean = nanmean(mbc_v);
                                    mbc_var = nanvar(mbc_v);
                                    mbc = nanmean(acu_bd_lo,2);
                                    
                                    acc_uc_lo = [acc_uc_lo mbc_mean];
                                    
                                    if length(loi)>=nMinTrials
                                        % if min(mbc)>=0
                                        if mbc_mean>0 FF_b_lo_temp = [FF_b_lo_temp mbc_var/mbc_mean]; end
                                        if mbc_mean==0 FF_b_lo_temp = [FF_b_lo_temp nan]; end
                                        mean_b_lo_temp = [mean_b_lo_temp mbc_mean];
                                        var_b_lo_temp = [var_b_lo_temp mbc_var];
                                    end
                                    if length(loi)<nMinTrials
                                        mean_b_lo_temp = [mean_b_lo_temp nan];
                                        var_b_lo_temp = [var_b_lo_temp nan];
                                        FF_b_lo_temp = [FF_b_lo_temp nan];
                                    end
                                    
                                    % HIGH:
                                    acu_bd_hi = acu_bd{acb}(hii,:); % pause.
                                    mbc_v = reshape(acu_bd_hi,size(acu_bd_hi,1)*size(acu_bd_hi,2),1);
                                    mbc_mean = nanmean(mbc_v);
                                    mbc_var = nanvar(mbc_v);
                                    mbc = nanmean(acu_bd_hi,2);
                                    
                                    acc_uc_hi = [acc_uc_hi mbc_mean];
                                    
                                    if length(hii)>=nMinTrials
                                        % if min(mbc)>=0
                                        if mbc_mean>0 FF_b_hi_temp = [FF_b_hi_temp mbc_var/mbc_mean]; end
                                        if mbc_mean==0 FF_b_hi_temp = [FF_b_hi_temp nan]; end
                                        mean_b_hi_temp = [mean_b_hi_temp mbc_mean];
                                        var_b_hi_temp = [var_b_hi_temp mbc_var];
                                    end
                                    if length(hii)<nMinTrials
                                        mean_b_hi_temp = [mean_b_hi_temp nan];
                                        var_b_hi_temp = [var_b_hi_temp nan];
                                        FF_b_hi_temp = [FF_b_hi_temp nan];
                                    end
                                    
                                    % ************************************
                                    
                                end % Binsizes.
                                
                                % pause
                                acc_bc_loT = [acc_bc_loT; acc_uc_lo];
                                %                                 acc_bc_midT = [acc_bc_midT; acc_uc_mid];
                                acc_bc_hiT = [acc_bc_hiT; acc_uc_hi];
                                
                                FF_b_loT = [FF_b_loT; FF_b_lo_temp];
                                %                                 FF_b_midT = [FF_b_midT; FF_b_mid_temp];
                                FF_b_hiT = [FF_b_hiT; FF_b_hi_temp];
                                
                                mean_b_loT = [mean_b_loT; mean_b_lo_temp];
                                %                                 mean_b_midT = [mean_b_midT; mean_b_mid_temp];
                                mean_b_hiT = [mean_b_hiT; mean_b_hi_temp];
                                
                                var_b_loT = [var_b_loT; var_b_lo_temp];
                                %                                 var_b_midT = [var_b_midT; var_b_mid_temp];
                                var_b_hiT = [var_b_hiT; var_b_hi_temp];
                                
                            end % LC units is did.
                            
                            dat2 = dat{2}; % ACC unit spike counts - over whole window.
                            
                            for ac = 1:length(dat2) % Loop over LC units to get counts for this brainstem unit.
                                
                                % Get LC counts corresponding to this ACC session.
                                dat3 = dat2{ac};
                                allZ = [allZ dat3(loi)]; allH = [allH dat3(hii)]; allA = [allA dat3];
                                
                                % Calculate the mean spike count over trials for each ACC unit, for this LC unit.
                                all_lo_c_meanT = [all_lo_c_meanT nanmean(dat3(loi))];
                                all_hi_c_meanT = [all_hi_c_meanT nanmean(dat3(hii))];
                                all_all_c_meanT = [all_all_c_meanT nanmean(dat3)];
                                
                            end % All LC units for this ACC session.
                        end
                    end % Check for empty dat structure.
                end % All sessions.
            end % Check for empty sitDat structure.
        end % Sites.
        
        %         all_pairs{mm} = pair_ind;
        
        all_lo_c_mean{mm} = all_lo_c_meanT; all_hi_c_mean{mm} = all_hi_c_meanT; all_all_c_mean{mm} = all_all_c_meanT; % Spike count matrix.
        var_lo{mm} = var_loT; var_hi{mm} = var_hiT; var_ALL{mm} = var_ALLT; % Variance of spike counts.
        
        unit_all_lo{mm} = unit_all_loT; unit_all_hi{mm} = unit_all_hiT; unit_all_all{mm} = unit_all_allT; % Pairwise spike count correlation, rsc.
        unit_all_lo_pv{mm} = unit_all_lo_pvT; unit_all_hi_pv{mm} = unit_all_hi_pvT; unit_all_all_pv{mm} = unit_all_all_pvT; % p-value from crosscorr calculation of rsc.
        unit_all_lo_p{mm} = unit_all_lo_pT; unit_all_hi_p{mm} = unit_all_hi_pT; unit_all_all_p{mm} = unit_all_all_pT; % rsc values from shuffled trials.
        
        acc_bc_lo{mm} = acc_bc_loT; acc_bc_hi{mm} = acc_bc_hiT; % Mean binned spike count for each ACC unit.
        unit_dat_rate_lo{mm} = unit_dat_rate_loT; unit_dat_rate_hi{mm} = unit_dat_rate_hiT; unit_dat_rate_al{mm} = unit_dat_rate_alT; % ACC neuron spike rates.
        
        FF_b_lo{mm} = FF_b_loT; FF_b_hi{mm} = FF_b_hiT;
        var_b_lo{mm} = var_b_loT; var_b_hi{mm} = var_b_hiT;
        mean_b_lo{mm} = mean_b_loT; mean_b_hi{mm} = mean_b_hiT;
        
        lc_sp_lo{mm} = unit_all_rate_lo;
        %         lc_sp_mid{mm} = unit_all_rate_mid;
        lc_sp_hi{mm} = unit_all_rate_hi;
        
        lc_AZ{mm} = unit_allZi;
        lc_AH{mm} = unit_allHi;
        
        nOut(mm) = nnuu;
        nOutA(mm) = nnAA;
        
    end % Mnks.
    
    %     eg_ze = {ze_pair_small_1 ze_pair_large_1 ze_pair_small_2 ze_pair_large_2};
    %     eg_hi = {hi_pair_small_1 hi_pair_large_1 hi_pair_small_2 hi_pair_large_2};
    
    eg_ze = [];
    eg_hi = [];
    
    % save('ACC_noBeep_rscData_100120','unit_all_lo','unit_all_hi','unit_all_all');
    save('ACC_noBeep_rscData_051621','unit_all_lo','unit_all_hi','unit_all_all');
    
end

%%

%     pause

% *************************************************
% % Now PLOT!
% *************************************************

% Don't plot results?
plotYesNo = false;
% Plot results?
plotYesNo = true;

% if plotYesNo

%% Figure 2: LC low/high; ACC spike counts, rates, variance and Fano factor.

% Setup figure:
% function [axs_,fig_,cap_] = getPLOT_axes(num, wid, hts, cols, psh, psw, fs, al, cap,centerTitle, fNum)
% figureNumber = 2; num = 2; wid = 17.6; hts = [5]; cols = {2 4 4}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 2, 2, [12], 'Joshi & Gold, 2018', true,1,figureNumber); set(axs,'Units','normalized');
figureNumber = 22; num = 22; wid = 17.6; hts = [6 4 4 4]; cols = {1 2 2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1, 2, [12], 'Joshi and Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
movegui(fig_,[1,1]); % Move figure so it is visible on screen.

colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};

%% low/high spike rate distributions.

%     bs = round(logspace(log10(100),log10(1000),10),1);
% xbs = linspace(100,1000,5); % 032320.
xbs = linspace(200,1000,5); % April 2021.
ms = [10 15 20];

hi_sp = [];

for jjm = 1:nMonks
    all_sp = [lc_sp_lo{jjm} lc_sp_hi{jjm}];
    hi_sp = [hi_sp all_sp/1.1];
end

xMax = 10; xMin = 0;
xBins = linspace(xMin,xMax,11);

all_lo = [lc_sp_lo{1}/1.1 lc_sp_lo{2}/1.1];
all_hi = [lc_sp_hi{1}/1.1 lc_sp_hi{2}/1.1];

prcLo = prctile(all_lo,[25 50 75]);
prcHi = prctile(all_hi,[25 50 75]);

n_all = hist(hi_sp,xBins);
NNN = sum(n_all);

axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
plot([prcLo(2)],[0.29],'ko','markersize',7,'linewidth',2);
% plot([prcLo(2)],[0.29],'cs','markersize',7,'linewidth',1);
plot([prcHi(2)],[0.35],'kd','markersize',7,'linewidth',2);
legend({'low','high'},'autoupdate','off','fontsize',8,'box','off'); % legend('boxoff');

plot([prcLo(1) prcLo(3)],[0.29 0.29],'k-','linewidth',2);
plot([prcHi(1) prcHi(3)],[0.35 0.35],'k-','linewidth',2);
% plot(xBins,n_all/NNN,'k-','linewidth',2);
bar(xBins,n_all/NNN,'edgecolor','none','facecolor','k');


xlabel('ACC firing rate (sp/s)','fontsize',10);
ylabel('proportion of trials','fontsize',10);
% title('LC','fontsize',12);
% if jjm == 1 title('LC spike counts'); end
% n_neurons = nOut(1)+nOut(2);
n_sessions = nOut(1)+nOut(2);
text(7,0.1,strcat('# num sessions=',num2str(n_sessions)),'fontsize',10);
text(7,0.15,strcat('# trials =',num2str(NNN)),'fontsize',10);
axis([-1 xMax+1 -0.01 0.36]);
text(5,0.45,'ACC','fontsize',12,'fontweight','bold');
set(gca,'fontname','arial');

% pause

%% LC mean spike count for each ACC spiking condition (Low and High).

axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

all_lo = [mean_b_lo{1};  mean_b_lo{2}]; gi1=find(isfinite(sum(all_lo,2))); mean_lo = nanmedian(all_lo(gi1,:)); % se_lo = nanse(all_lo(gi1,:),1);
all_hi = [mean_b_hi{1}; mean_b_hi{2}]; gi3=find(isfinite(sum(all_hi,2))); mean_hi = nanmedian(all_hi(gi3,:)); % se_hi = nanse(all_hi(gi3,:),1);

[ci_1 ci_1_m]= bootci(2000,{@median,all_lo(gi1,:)},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
[ci_3 ci_3_m]= bootci(2000,{@median,all_hi(gi3,:)},'alpha',0.05,'type','per'); % mean_3 = nanmean(ci_3_m);

% *************
% Set up ANOVA.
% Number of data samples for each LC spike condition and each monkey.
nlo1 = length(mean_b_lo{1}); nlo2 = length(mean_b_lo{2});
nhi1 = length(mean_b_hi{1}); nhi2 = length(mean_b_hi{2});

% LC=0,1,2,3,>3 correspond to indices 1,2,3,4,5.
LC_sp_ind_lo1 = 1*ones(nlo1,1); LC_sp_ind_lo2 = 1*ones(nlo2,1);
LC_sp_ind_hi1 = 3*ones(nhi1,1); LC_sp_ind_hi2 = 3*ones(nhi2,1);

dat1_LCspInd = [LC_sp_ind_lo1; LC_sp_ind_hi1];
dat2_LCspInd = [LC_sp_ind_lo2; LC_sp_ind_hi2];

pvals = nans(length(xbs),2);
pvAll = nans(1,length(xbs));

for uj = 1:length(xbs)
    
    dat1 = [mean_b_lo{1}(:,uj); mean_b_hi{1}(:,uj)];
    dat2 = [mean_b_lo{2}(:,uj); mean_b_hi{2}(:,uj)];
    
    % 1-way ANOVA with display OFF.
    [p1,tbl,stats] = anova1(dat1, dat1_LCspInd,'off');
    [p2,tbl,stats] = anova1(dat2, dat2_LCspInd,'off');
    pvals(uj,:) = [p1 p2];
    
    p21 = ranksum(mean_b_hi{1}(:,uj),mean_b_lo{1}(:,uj));
    p22 = ranksum(mean_b_hi{2}(:,uj),mean_b_lo{2}(:,uj));
    
    plot(xbs(uj)-15,mean_lo(uj),'ko','markersize',5,'linewidth',1);
    plot(xbs(uj)+15,mean_hi(uj),'kd','markersize',5,'linewidth',1);
    
    plot([xbs(uj)-15 xbs(uj)-15],ci_1(:,uj),'k-','markersize',5,'linewidth',0.5);
    plot([xbs(uj)+15 xbs(uj)+15],ci_3(:,uj),'k-','markersize',5,'linewidth',0.5);
    
    if p21<0.05 & p22<0.05 pvAll(uj) = 1; end
    
    if uj == 1
        legend({'ACC low','ACC high'},'autoupdate','off','fontsize',8,'box','off','location','northwest'); % legend('boxoff');
    end
end

% ANOVA is done.
% *************
% axis([50 1050 -0.2 3]);
axis([100 1100 0 2.2]);
xlabel('bin size (ms)','fontsize',10); ylabel('LC spike count','fontsize',10);
set(gca,'fontname','arial');

%% Plot bars for change in LC spike count relative to ACC low.

D1_hi1 = mean_b_hi{1}-mean_b_lo{1};
D1_hi2 = mean_b_hi{2}-mean_b_lo{2};

D1 = [D1_hi1; D1_hi2]; % ginz=find(isfinite(sum(all_dnz,2))); mean_dnz = nanmean(all_dnz(ginz,:));

mean_dhi1 = nans(5,1); ci_hi1 = nans(5,2); ppv1 = nans(5,1);
mean_dhi2 = nans(5,1); ci_hi2 = nans(5,2); ppv2 = nans(5,1);
ppv_all = nans(5,1);
xjit = 30;

for ti = 1:length(xbs)
    
    datTemp1 = D1_hi1(:,ti); gihi1 = find(isfinite(datTemp1)); datTemp1 = datTemp1(gihi1); % datTemp = datTemp(datTemp<prctile(datTemp,limitVal));
    [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@median,datTemp1},'alpha',0.05,'type','per');
    mean_dhi1(ti) = nanmedian(datTemp1);
    ci_hi1(ti,:) = ci_hi_temp1';
    ppv1(ti) = signrank(datTemp1);
    
    datTemp2 = D1_hi2(:,ti); gihi2 = find(isfinite(datTemp2)); datTemp2 = datTemp2(gihi2); % datTemp = datTemp(datTemp<prctile(datTemp,limitVal));
    [ci_hi_temp2 ci_hi_m_temp2]= bootci(2000,{@median,datTemp2},'alpha',0.05,'type','per');
    mean_dhi2(ti) = nanmedian(datTemp2);
    ci_hi2(ti,:) = ci_hi_temp2';
    ppv2(ti) = signrank(datTemp2);
    
    ppv_all(ti) = signrank(D1(:,ti));
    
end

axes(axs(3)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

xjit = 30;
for uj = 1:length(xbs)
    
    plot([xbs(uj)-xjit xbs(uj)-xjit],ci_hi1(uj,:),'k-');
    if ppv1(uj) < 0.05
        bar(xbs(uj)-xjit,mean_dhi1(uj),2*(xjit-2),'facecolor',colors{2},'edgecolor','none','linewidth',1);
    end
    if ppv1(uj) >= 0.05
        bar(xbs(uj)-xjit,mean_dhi1(uj),2*(xjit-2),'facecolor','none','edgecolor',colors{2},'linewidth',1);
    end
    
    plot([xbs(uj)+xjit xbs(uj)+xjit],ci_hi2(uj,:),'k-');
    if ppv2(uj) < 0.05
        bar(xbs(uj)+xjit,mean_dhi2(uj),2*(xjit-2),'facecolor',colors{3},'edgecolor','none','linewidth',1);
    end
    if ppv2(uj) >= 0.05
        bar(xbs(uj)+xjit,mean_dhi2(uj),2*(xjit-2),'facecolor','none','edgecolor',colors{3},'linewidth',1);
    end
    
    if ppv_all(uj)<0.05
        plot(xbs(uj),0.1,'k*');
    end
    
end

axis([100 1100 -0.1 0.1]);
xlabel('bin size (ms)','fontsize',10); ylabel({'LC spike count';'difference re: ACC low'},'fontsize',10);
set(gca,'fontname','arial');

%% LC variance of spike count for each ACC spiking condition (Low, mid, high).

axes(axs(4)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

all_lo = [var_b_lo{1}; var_b_lo{2}]; gi1=find(isfinite(sum(all_lo,2))); mean_lo = nanmedian(all_lo(gi1,:));
all_hi = [var_b_hi{1}; var_b_hi{2}]; gi3=find(isfinite(sum(all_hi,2))); mean_hi = nanmedian(all_hi(gi3,:));

[ci_1 ci_1_m]= bootci(2000,{@median,all_lo(gi1,:)},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
[ci_3 ci_3_m]= bootci(2000,{@median,all_hi(gi3,:)},'alpha',0.05,'type','per'); % mean_3 = nanmean(ci_3_m);

% *************
% Set up ANOVA.
% Number of data samples for each LC spike condition and each monkey.
nlo1 = length(var_b_lo{1}); nlo2 = length(var_b_lo{2});
nhi1 = length(var_b_hi{1}); nhi2 = length(var_b_hi{2});

% LC=0,1,2,3,>3 correspond to indices 1,2,3,4,5.
LC_sp_ind_lo1 = 1*ones(nlo1,1); LC_sp_ind_lo2 = 1*ones(nlo2,1);
LC_sp_ind_hi1 = 3*ones(nhi1,1); LC_sp_ind_hi2 = 3*ones(nhi2,1);

dat1_LCspInd = [LC_sp_ind_lo1; LC_sp_ind_hi1];
dat2_LCspInd = [LC_sp_ind_lo2; LC_sp_ind_hi2];

pvals = nans(length(xbs),2);
pvAllV = nans(1,length(xbs));

for uj = 1:length(xbs)
    dat1 = [var_b_lo{1}(:,uj); var_b_hi{1}(:,uj)];
    dat2 = [var_b_lo{2}(:,uj); var_b_hi{2}(:,uj)];
    % 1-way ANOVA with display OFF.
    [p1,tbl,stats] = anova1(dat1, dat1_LCspInd,'off');
    [p2,tbl,stats] = anova1(dat2, dat2_LCspInd,'off');
    pvals(uj,:) = [p1 p2];
    
    p21 = ranksum(var_b_hi{1}(:,uj),var_b_lo{1}(:,uj));
    p22 = ranksum(var_b_hi{2}(:,uj),var_b_lo{2}(:,uj));
    
    plot(xbs(uj)-15,mean_lo(uj),'ko','markersize',5,'linewidth',1);
    plot(xbs(uj)+15,mean_hi(uj),'kd','markersize',5,'linewidth',1);
    
    plot([xbs(uj)-15 xbs(uj)-15],ci_1(:,uj),'k-','markersize',5,'linewidth',1);
    plot([xbs(uj)+15 xbs(uj)+15],ci_3(:,uj),'k-','markersize',5,'linewidth',1);
    
    if p21<0.05 & p22<0.05 plot(xbs(uj),4.4,'k*'); end
    if p21<0.05 & p22<0.05 pvAllV(uj) =1; end
    
    %     if uj == 1
    %         legend({'ACC low','ACC high'},'autoupdate','off','fontsize',8,'box','off','location','best'); % legend('boxoff');
    %     end
end

% ANOVA is done.
% *************
axis([100 1100 -0.5 2.5]);
xlabel('bin size (ms)','fontsize',10); ylabel('Variance of LC spike count','fontsize',10);
set(gca,'fontname','arial');

%% Plot bars for change in variance of spike count relative to ACC low.

D1_hi1 = var_b_hi{1}-var_b_lo{1};
D1_hi2 = var_b_hi{2}-var_b_lo{2};

D1 = [D1_hi1; D1_hi2]; % ginz=find(isfinite(sum(all_dnz,2))); mean_dnz = nanmean(all_dnz(ginz,:));

mean_dhi1 = nans(5,1); ci_hi1 = nans(5,2); ppv1 = nans(5,1);
mean_dhi2 = nans(5,1); ci_hi2 = nans(5,2); ppv2 = nans(5,1);
ppv_all = nans(5,1);

for ti = 1:length(xbs)
    
    datTemp1 = D1_hi1(:,ti); gihi1 = find(isfinite(datTemp1)); datTemp1 = datTemp1(gihi1); % datTemp = datTemp(datTemp<prctile(datTemp,limitVal));
    [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@median,datTemp1},'alpha',0.05,'type','per');
    mean_dhi1(ti) = nanmedian(datTemp1);
    ci_hi1(ti,:) = ci_hi_temp1';
    ppv1(ti) = signrank(datTemp1);
    
    datTemp2 = D1_hi2(:,ti); gihi2 = find(isfinite(datTemp2)); datTemp2 = datTemp2(gihi2); % datTemp = datTemp(datTemp<prctile(datTemp,limitVal));
    [ci_hi_temp2 ci_hi_m_temp2]= bootci(2000,{@median,datTemp2},'alpha',0.05,'type','per');
    mean_dhi2(ti) = nanmedian(datTemp2);
    ci_hi2(ti,:) = ci_hi_temp2';
    ppv2(ti) = signrank(datTemp2);
    
    ppv_all(ti) = signrank(D1(:,ti));
    
end

axes(axs(5)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

for uj = 1:length(xbs)
    
    plot([xbs(uj)-xjit xbs(uj)-xjit],ci_hi1(uj,:),'k-');
    if ppv1(uj) < 0.05
        bar(xbs(uj)-xjit,mean_dhi1(uj),2*(xjit-2),'facecolor',colors{2},'edgecolor','none','linewidth',1);
    end
    if ppv1(uj) >= 0.05
        bar(xbs(uj)-xjit,mean_dhi1(uj),2*(xjit-2),'facecolor','none','edgecolor',colors{2},'linewidth',1);
    end
    
    plot([xbs(uj)+xjit xbs(uj)+xjit],ci_hi2(uj,:),'k-');
    if ppv2(uj) < 0.05
        bar(xbs(uj)+xjit,mean_dhi2(uj),2*(xjit-2),'facecolor',colors{3},'edgecolor','none','linewidth',1);
    end
    if ppv2(uj) >= 0.05
        bar(xbs(uj)+xjit,mean_dhi2(uj),2*(xjit-2),'facecolor','none','edgecolor',colors{3},'linewidth',1);
    end
    
    if ppv_all(uj)<0.05
        plot(xbs(uj),0.25,'k*');
    end
    
end

axis([100 1100 -0.25 0.25]);
xlabel('bin size (ms)','fontsize',10); ylabel({'Variance of';'LC spike count';'difference re: ACC low'},'fontsize',10);
set(gca,'fontname','arial');

%% LC Fano factor for each ACC spiking condition (Low, mid, high).

axes(axs(6)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

all_lo = [FF_b_lo{1}; FF_b_lo{2}]; gi1=find(isfinite(sum(all_lo,2))); mean_lo = nanmedian(all_lo(gi1,:));
all_hi = [FF_b_hi{1}; FF_b_hi{2}]; gi3=find(isfinite(sum(all_hi,2))); mean_hi = nanmedian(all_hi(gi3,:));

[ci_1 ci_1_m]= bootci(2000,{@median,all_lo(gi1,:)},'alpha',0.05,'type','per'); % mean_1 = nanmean(ci_1_m);
[ci_3 ci_3_m]= bootci(2000,{@median,all_hi(gi3,:)},'alpha',0.05,'type','per'); % mean_3 = nanmean(ci_3_m);

% *************
% Set up ANOVA.
% Number of data samples for each LC spike condition and each monkey.
nlo1 = length(FF_b_lo{1}); nlo2 = length(FF_b_lo{2});
nhi1 = length(FF_b_hi{1}); nhi2 = length(FF_b_hi{2});

% LC=0,1,2,3,>3 correspond to indices 1,2,3,4,5.
LC_sp_ind_lo1 = 1*ones(nlo1,1); LC_sp_ind_lo2 = 1*ones(nlo2,1);
LC_sp_ind_hi1 = 3*ones(nhi1,1); LC_sp_ind_hi2 = 3*ones(nhi2,1);

dat1_LCspInd = [LC_sp_ind_lo1; LC_sp_ind_hi1];
dat2_LCspInd = [LC_sp_ind_lo2; LC_sp_ind_hi2];

pvals = nans(length(xbs),2);
pvAllF = nans(1,length(xbs));

for uj = 1:length(xbs)
    dat1 = [FF_b_lo{1}(:,uj); FF_b_hi{1}(:,uj)];
    dat2 = [FF_b_lo{2}(:,uj); FF_b_hi{2}(:,uj)];
    % 1-way ANOVA with display OFF.
    [p1,tbl,stats] = anova1(dat1, dat1_LCspInd,'off');
    [p2,tbl,stats] = anova1(dat2, dat2_LCspInd,'off');
    pvals(uj,:) = [p1 p2];
    
    p21 = ranksum(FF_b_hi{1}(:,uj),FF_b_lo{1}(:,uj));
    p22 = ranksum(FF_b_hi{2}(:,uj),FF_b_lo{2}(:,uj));
    
    plot(xbs(uj)-15,mean_lo(uj),'ko','markersize',5,'linewidth',1);
    plot(xbs(uj)+15,mean_hi(uj),'kd','markersize',5,'linewidth',1);
    
    plot([xbs(uj)-15 xbs(uj)-15],ci_1(:,uj),'k-','markersize',5,'linewidth',1);
    plot([xbs(uj)+15 xbs(uj)+15],ci_3(:,uj),'k-','markersize',5,'linewidth',1);
    
    % if p21<0.05 & p22<0.05 plot(xbs(uj),0.25+mean_mid(uj),'k*'); end
    if p21<0.05 & p22<0.05 pvAllF(uj) =1; end
    
    %     if uj == 1
    %         legend({'ACC low','ACC high'},'autoupdate','off','fontsize',8,'box','off','location','best'); % legend('boxoff');
    %     end
end

% ANOVA is done.
% *************
axis([100 1100 0.9 1.25]);
xlabel('bin size (ms)','fontsize',10); ylabel('LC Fano factor','fontsize',10);
set(gca,'fontname','arial');

%% Plot bars for %change in Fano factor relative to ACC low.

D1_hi1 = FF_b_hi{1}-FF_b_lo{1};
D1_hi2 = FF_b_hi{2}-FF_b_lo{2};

D1 = [D1_hi1; D1_hi2]; % ginz=find(isfinite(sum(all_dnz,2))); mean_dnz = nanmean(all_dnz(ginz,:));

mean_dhi1 = nans(5,1); ci_hi1 = nans(5,2); ppv1 = nans(5,1);
mean_dhi2 = nans(5,1); ci_hi2 = nans(5,2); ppv2 = nans(5,1);
ppv_all = nans(5,1);

for ti = 1:length(xbs)
    
    datTemp1 = D1_hi1(:,ti); gihi = find(isfinite(datTemp1)); datTemp1 = datTemp1(gihi); % datTemp = datTemp(datTemp<prctile(datTemp,limitVal));
    [ci_hi_temp1 ci_hi_m_temp1]= bootci(2000,{@median,datTemp1},'alpha',0.05,'type','per');
    mean_dhi1(ti) = nanmedian(datTemp1);
    ci_hi1(ti,:) = ci_hi_temp1';
    ppv1(ti) = signrank(datTemp1);
    
    datTemp2 = D1_hi2(:,ti); gihi = find(isfinite(datTemp2)); datTemp2 = datTemp2(gihi); % datTemp = datTemp(datTemp<prctile(datTemp,limitVal));
    [ci_hi_temp2 ci_hi_m_temp2]= bootci(2000,{@median,datTemp2},'alpha',0.05,'type','per');
    mean_dhi2(ti) = nanmedian(datTemp2);
    ci_hi2(ti,:) = ci_hi_temp2';
    ppv2(ti) = signrank(datTemp2);
    
    ppv_all(ti) = signrank(D1(:,ti));
    
end

axes(axs(7)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

for uj = 1:length(xbs)
    
    plot([xbs(uj)-xjit xbs(uj)-xjit],ci_hi1(uj,:),'k-');
    if ppv1(uj) < 0.05
        bar(xbs(uj)-xjit,mean_dhi1(uj),2*(xjit-2),'facecolor',color{2},'edgecolor','none','linewidth',1);
    end
    if ppv1(uj) >= 0.05
        bar(xbs(uj)-xjit,mean_dhi1(uj),2*(xjit-2),'facecolor','none','edgecolor',colors{2},'linewidth',1);
    end
    
    plot([xbs(uj)+xjit xbs(uj)+xjit],ci_hi2(uj,:),'k-');
    if ppv2(uj) < 0.05
        bar(xbs(uj)+xjit,mean_dhi2(uj),2*(xjit-2),'facecolor',color{3},'edgecolor','none','linewidth',1);
    end
    if ppv2(uj) >= 0.05
        bar(xbs(uj)+xjit,mean_dhi2(uj),2*(xjit-2),'facecolor','none','edgecolor',colors{3},'linewidth',1);
    end
    
    if ppv_all(uj)<0.05
        plot(xbs(uj),0.12,'k*');
    end
    
end

axis([100 1100 -0.12 0.12]);
xlabel('bin size (ms)','fontsize',10); ylabel({'LC Fano factor';'difference re: ACC low'},'fontsize',10);
set(gca,'fontname','arial');


%% Supplement: Fano and firing rate.

% figureNumber = 223; num = 223; wid = 17.6; hts = [8]; cols = {2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,3,1, [12], 'Joshi & Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
figureNumber = 23; num = 23; wid = 17.6; hts = [12]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,2,2, [12], 'Joshi and Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
movegui(fig_,[1,1]); % Move figure so it is visible on screen.

all_ze = [mean_b_lo{1}(:,10); mean_b_lo{2}(:,10)];
all_nz = [mean_b_hi{1}(:,10); mean_b_hi{2}(:,10)];

all_zeF = [FF_b_lo{1}(:,10); FF_b_lo{2}(:,10)];
all_nzF = [FF_b_hi{1}(:,10); FF_b_hi{2}(:,10)];

all_diff_Rate = (all_nz-all_ze)./all_ze;
all_diff_FF = (all_nzF-all_zeF)./all_zeF;

%     gi = find(isfinite(all_diff_FF) & isfinite(all_diff_Rate) & abs(all_diff_FF)<2 & abs(all_diff_Rate)<2);
gi = find(isfinite(all_diff_FF) & isfinite(all_diff_Rate));

mF = nanmedian(all_diff_FF(gi));
mR = nanmedian(all_diff_Rate(gi));

%     [b,bint,r,rint,stats] = regress(all_diff_FF(gi),[ones(length(all_diff_Rate(gi)),1) all_diff_Rate(gi)]);
%     xx = linspace(min(all_diff_Rate(gi)),max(all_diff_Rate(gi)),10); yy = b(1)+b(2)*xx;

xb = linspace(-1,1,30); yb = xb;
xDat = all_diff_Rate(gi);
yDat = all_diff_FF(gi);
nx = hist(xDat,xb); mx = nanmedian(xDat);
ny = hist(yDat,yb); my = nanmedian(yDat);

axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

% figure; hold on;

% histogram2(xDat,yDat,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');
histogram2(xDat,yDat,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgealpha',0);
colormap(gray);

% plot(xDat,yDat,'k.','markersize',10,'color',[.5 .5 .5]);
% plot(mF,mR,'r.','markersize',10);
% plot([-1 1],[-1 1],'g:');
plot([-1 1],[0 0],'w:');
plot([0 0],[-1 1],'w:');
% plot(xx,yy,'r-');
title({'Relation between changes in';'LC neuron firing rate and Fano factor'},'fontsize',10);
xlabel('\Delta spike count','fontsize',10);
ylabel('\Delta Fano factor','fontsize',10);
cbh = colorbar; set(cbh,'box','off');
ylabel(cbh,'number of neurons','fontsize',10);
axis([-1 1 -1 1]);
axis square;
set(gca,'fontname','arial');

%% Supplement: Fano and variance of spike count.

%     % figureNumber = 223; num = 223; wid = 17.6; hts = [8]; cols = {2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,3,1, [12], 'Joshi & Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
%     figureNumber = 21; num = 21; wid = 17.6; hts = [12]; cols = {1}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols,1,1, [12], 'Joshi & Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
%     movegui(fig_,[1,1]); % Move figure so it is visible on screen.

all_ze = [var_b_lo{1}(:,10); var_b_lo{2}(:,10)];
all_nz = [var_b_hi{1}(:,10); var_b_hi{2}(:,10)];

all_zeF = [FF_b_lo{1}(:,10); FF_b_lo{2}(:,10)];
all_nzF = [FF_b_hi{1}(:,10); FF_b_hi{2}(:,10)];

all_diff_Rate = (all_nz-all_ze)./all_ze;
all_diff_FF = (all_nzF-all_zeF)./all_zeF;

%     gi = find(isfinite(all_diff_FF) & isfinite(all_diff_Rate) & abs(all_diff_FF)<2 & abs(all_diff_Rate)<2);
gi = find(isfinite(all_diff_FF) & isfinite(all_diff_Rate));

mF = nanmedian(all_diff_FF(gi));
mR = nanmedian(all_diff_Rate(gi));

%     [b,bint,r,rint,stats] = regress(all_diff_FF(gi),[ones(length(all_diff_Rate(gi)),1) all_diff_Rate(gi)]);
%     xx = linspace(min(all_diff_Rate(gi)),max(all_diff_Rate(gi)),10); yy = b(1)+b(2)*xx;

xb = linspace(-1,1,30); yb = xb;
xDat = all_diff_Rate(gi);
yDat = all_diff_FF(gi);
nx = hist(xDat,xb); mx = nanmedian(xDat);
ny = hist(yDat,yb); my = nanmedian(yDat);

axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);

% figure; hold on;

% histogram2(xDat,yDat,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgecolor','none');
histogram2(xDat,yDat,xb,yb,'DisplayStyle','tile','ShowEmptyBins','on','edgealpha',0);
colormap(gray);

% plot(xDat,yDat,'k.','markersize',10,'color',[.5 .5 .5]);
% plot(mF,mR,'r.','markersize',10);
% plot([-1 1],[-1 1],'g:');
plot([-1 1],[0 0],'w:');
plot([0 0],[-1 1],'w:');
% plot(xx,yy,'r-');
title({'Relation between changes in';'LC neuron spike count variance and Fano factor'},'fontsize',10);
xlabel('\Delta variance of spike count','fontsize',10);
ylabel('\Delta Fano factor','fontsize',10);
cbh = colorbar; set(cbh,'box','off');
ylabel(cbh,'number of neurons','fontsize',10);
axis([-1 1 -1 1]);
axis square;
set(gca,'fontname','arial');

% Finis.


