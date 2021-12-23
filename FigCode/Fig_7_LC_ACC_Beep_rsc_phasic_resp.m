% Fig7_LC_ACC_Beep_rsc_phasic_resp.m

% Origin: 022019, Sidd.
% Calculate binned rsc for ACC pairs conditioned on phasic or no phasic LC responses to beep.
% 
% Added Fano and spike rate calculation for LC phasic/no phasic responses.
%
% Modified from: Fig8_LC_Dual_Recording_rsc.m
% Mod 122818 - Sidd. From: "LC_Dual_Recording_xcor.m".
% Mod: 080118 - Sidd.
% Mod: 092718: Added rCCG calculation as done for LC, IC and ACC pairs.
%
% Fixation expts ONLY.
% BEEP trials only.
% Calls "getLC_DualArea_Beep_binned_rsc.m" and "getLC_DualArea_FakeBeep_binned_rsc.m" to do the heavy lifting.

%% Calculate phasic/no-phasic rsc.

% Clear everything?
clear; clear all;

% Don't reanalyze data?
reanalyze = false;
% % Reanalyze data?
% reanalyze = true;

if reanalyze
    clear; clear all;
    %     NTL=10; % Before 042921.
    NTL=5;
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
                    load(fnames{ff});
                    % pause
                    rsc_dat{ff} = getLC_DualArea_Beep_binned_rsc(siteData,NTL,ff); % Pairwise spike count correlations.
                    % rsc_dat{ff} = getLC_DualArea_Beep_rsc(siteData,NTL); % Pairwise spike count correlations.
                    disp(sprintf('File %d%sof%s%d', ff,' ',' ',nf));
                end
            end
            sitRes{ss} = rsc_dat;
            clear rsc_dat;
        end % Sites.
        mnkRes{mm} = sitRes;
        
        % Save analysis:
        clear sitRes;
        
    end % Mnks.
    
    % savedir = strcat([base_dir,'\Results\Results_2020\LC_ACC']);
    savedir = strcat([base_dir,'\Results\Results_2021\LC_ACC']);
    cd(savedir);
    
    % saveFileName = 'LC_ACC_Results_Beep_rsc_031721';
    % saveFileName = 'LC_ACC_Results_Beep_rsc_042121';
    % saveFileName = 'LC_ACC_Results_Beep_rsc_042921c';
    % saveFileName = 'LC_ACC_Results_Beep_rsc_050521';
    % saveFileName = 'LC_ACC_Results_Beep_rsc_052121b';
    % saveFileName = 'LC_ACC_Results_Beep_rsc_052521';
    saveFileName = 'LC_ACC_Results_Beep_rsc_100321';
    
    save(saveFileName,'mnkRes'); % CCG results are LARGE... need to use -v7.3
    clear mnkRes;
end

%% Calculate phasic/no-phasic rsc for FAKE "phasic" events.

% Don't reanalyze data?
reanalyze = false;
% % Reanalyze data?
% reanalyze = true;

if reanalyze
    clear; clear all;
    %     NTL=10;
    NTL=5;
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
                    load(fnames{ff});
                    % pause
                    rsc_dat{ff} = getLC_DualArea_FakeBeep_binned_rsc(siteData,NTL,ff); % Pairwise spike count correlations.
                    % rsc_dat{ff} = getLC_DualArea_Beep_rsc(siteData,NTL); % Pairwise spike count correlations.
                    disp(sprintf('File %d%sof%s%d', ff,' ',' ',nf));
                end
            end
            sitRes{ss} = rsc_dat;
            clear rsc_dat;
        end % Sites.
        mnkRes{mm} = sitRes;
        
        % Save analysis:
        clear sitRes;
        
    end % Mnks.
    
    % savedir = strcat([base_dir,'\Results\Results_2020\LC_ACC']);
    savedir = strcat([base_dir,'\Results\Results_2021\LC_ACC']);
    cd(savedir);
    
    % saveFileName = 'LC_ACC_Results_FakeBeep_rsc_031721';
    % saveFileName = 'LC_ACC_Results_FakeBeep_rsc_042121';
    % saveFileName = 'LC_ACC_Results_FakeBeep_rsc_042921c';
    % saveFileName = 'LC_ACC_Results_FakeBeep_rsc_052121b';
    saveFileName = 'LC_ACC_Results_FakeBeep_rsc_100321';
    
    save(saveFileName,'mnkRes'); % CCG results are LARGE... need to use -v7.3
    clear mnkRes;
end

%%

clear; clear all;

base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.
% savedir = strcat([base_dir,'\Results\Results_2020\LC_ACC\']); %  Point to the current results directory.
savedir = strcat([base_dir,'\Results\Results_2021\LC_ACC\']); %  Point to the current results directory.
cd(savedir);

monks = {'Sprout','Cicero'}; % Add monks as needed.
sites = {'LC_ACC_Fixation'}; % Add sites as needed.
nmonks = length(monks);
mnkFigDat = cell(1,2);

%% Organize rsc plot data for beep trials.

organizeYesno = 0;
organizeYesno = 1;

if organizeYesno
    
    binSiz = 1; pairs_mnk2 = [];
    
    % saveFileName = 'LC_ACC_Results_Beep_rsc_042121';
    % saveFileName = 'LC_ACC_Results_Beep_rsc_031721';
    % saveFileName = 'LC_ACC_Results_Beep_rsc_042921c';
    % saveFileName = 'LC_ACC_Results_Beep_rsc_050521';
    % saveFileName = 'LC_ACC_Results_Beep_rsc_052121b';
    % saveFileName = 'LC_ACC_Results_Beep_rsc_052521';
    saveFileName = 'LC_ACC_Results_Beep_rsc_100321';
    
    load(saveFileName);
    % bs = linspace(100,1000,10); % 021720.
    bs = linspace(200,1000,5); % 042921.
    
    for mm = 1:nmonks
        
        pairs = [];
        unit_index = 0; incindex = 0;
        ns = length(mnkRes{mm}{1}); % Number of sessions for this mnky.
        yyi = 0;
        all_loBnp = []; all_loAnp = []; all_loBp = [];  all_loAp = [];
        all_loBnp2 = []; all_loAnp2 = []; all_loBp2 = [];  all_loAp2 = [];
        all_loBnp3 = []; all_loAnp3 = []; all_loBp3 = [];  all_loAp3 = [];
        all_loBnp4 = []; all_loAnp4 = []; all_loBp4 = [];  all_loAp4 = [];
        all_lohi_count = [];
        si = 0; si2 = 0; si3 = 0; si4 = 0;
        si11 = 0; si12 = 0; si13 = 0; si14 = 0;
        pnp_dat = []; pnp_dat2 = []; pnp_dat3 = []; pnp_dat4 = [];
        phasicUnit = []; phasicUnit2 = [];
        sp_dat = []; mean_sp_dat = [];
        all_file = []; all_lci = []; all_acpi = [];
        all_file2 = []; all_lci2 = []; all_acpi2 = [];
        all_file3 = []; all_lci3 = []; all_acpi3 = [];
        all_file4 = []; all_lci4 = []; all_acpi4 = [];
        
        all_loApE=[]; all_loApL=[]; all_loAnpE=[]; all_loAnpL=[];
        
        all_mean_B = []; all_mean_A = []; all_var_B = []; all_var_A = [];
        all_mean_B3 = []; all_mean_A3 = []; all_var_B3 = []; all_var_A3 = [];
        
        % pause
        
        for ii = 1:ns
            nu = length(mnkRes{mm}{1}{ii}{1}); % Get site data (one cell for each LC unit).
            for jj = 1:nu
                
                unit_index = unit_index+1;
                
                sp_dat = [sp_dat; zscore(mnkRes{mm}{1}{ii}{2}{jj}')'];
                pnp_dat = [pnp_dat mnkRes{mm}{1}{ii}{3}{jj}];
                phasicYesNo = mnkRes{mm}{1}{ii}{3}{jj};
                phasicYesNoProportion = sum(isfinite(phasicYesNo))/length(phasicYesNo);
                phasicUnit{unit_index} = phasicYesNoProportion;
                % phasicUnit{unit_index} = mnkRes{mm}{1}{ii}{3}{jj}; % pause
                % phasicUnit2{unit_index} = mnkRes{mm}{1}{ii}{4}{jj};
                bDat = mnkRes{mm}{1}{ii}{1}{jj}; % pause
                pnp_dat3 = [pnp_dat3 mnkRes{mm}{1}{ii}{5}{jj}];
                bDat3 = mnkRes{mm}{1}{ii}{8}{jj};
                
                if ~isempty(bDat) & ~isempty(bDat3)
                    nbins = length(bDat);
                    for kk = 1:nbins
                        
                        ei = length(bDat{kk}{1}(:,1));
                        all_loBnp(1+si:si+ei,kk) = bDat{kk}{3}(:,1);
                        all_loAnp(1+si:si+ei,kk) = bDat{kk}{4}(:,1);
                        all_loBp(1+si:si+ei,kk) = bDat{kk}{1}(:,1);
                        all_loAp(1+si:si+ei,kk) = bDat{kk}{2}(:,1);
                        
                        all_loApE(1+si:si+ei,kk) = bDat{kk}{9}(:,1);
                        all_loApL(1+si:si+ei,kk) = bDat{kk}{10}(:,1);
                        
                        ei11 = length(bDat{kk}{5});
                        all_mean_B(1+si11:si11+ei11,kk) = bDat{kk}{5};
                        all_var_B(1+si11:si11+ei11,kk) = bDat{kk}{6};
                        all_mean_A(1+si11:si11+ei11,kk) = bDat{kk}{7};
                        all_var_A(1+si11:si11+ei11,kk) = bDat{kk}{8};
                        
                    end % Bins.
                    
                    all_file(1+si:si+ei,1) = bDat{1}{1}(:,3);
                    all_lci(1+si:si+ei,1) = bDat{1}{1}(:,4);
                    all_acpi(1+si:si+ei,1) = bDat{1}{1}(:,5);
                    
                    si = si+ei;
                    si11 = si11+ei11;
                    
                    %                 end % Check for empty data struct.
                    
                    nbins = length(bDat3);
                    for kk = 1:nbins
                        ei = length(bDat3{kk}{1}(:,1));
                        all_loBnp3(1+si3:si3+ei,kk) = bDat3{kk}{3}(:,1);
                        all_loAnp3(1+si3:si3+ei,kk) = bDat3{kk}{4}(:,1);
                        all_loBp3(1+si3:si3+ei,kk) = bDat3{kk}{1}(:,1);
                        all_loAp3(1+si3:si3+ei,kk) = bDat3{kk}{2}(:,1);
                        
                        all_loAnpE(1+si3:si3+ei,kk) = bDat3{kk}{9}(:,1);
                        all_loAnpL(1+si3:si3+ei,kk) = bDat3{kk}{10}(:,1);
                        
                        ei13 = length(bDat3{kk}{5});
                        all_mean_B3(1+si13:si13+ei13,kk) = bDat3{kk}{5};
                        all_var_B3(1+si13:si13+ei13,kk) = bDat3{kk}{6};
                        all_mean_A3(1+si13:si13+ei13,kk) = bDat3{kk}{7};
                        all_var_A3(1+si13:si13+ei13,kk) = bDat3{kk}{8};
                        
                    end % Bins.
                    
                    all_file3(1+si3:si3+ei,1) = bDat3{1}{1}(:,3);
                    all_lci3(1+si3:si3+ei,1) = bDat3{1}{1}(:,4);
                    all_acpi3(1+si3:si3+ei,1) = bDat3{1}{1}(:,5);
                    
                    si3 = si3+ei;
                    si13 = si13+ei13;
                    
                end % Check for empty data struct.
                
            end % Units for this session.
        end % Sessions for this monkey.
        mnkFigDat{mm} = {all_loBnp all_loAnp all_loBp all_loAp all_file all_lci all_acpi ...
            all_loBnp2 all_loAnp2 all_loBp2 all_loAp2 all_file2 all_lci2 all_acpi2 ...
            all_loBnp3 all_loAnp3 all_loBp3 all_loAp3 all_file3 all_lci3 all_acpi3 ...
            all_loBnp4 all_loAnp4 all_loBp4 all_loAp4 all_file4 all_lci4 all_acpi4 ...
            pnp_dat pnp_dat2 pnp_dat3 pnp_dat4 ...
            sp_dat all_mean_B all_mean_A all_var_B all_var_A all_mean_B3 all_mean_A3 all_var_B3 all_var_A3 phasicUnit...
            all_loApE(:,1) all_loApL(:,1) all_loAnpE(:,1) all_loAnpL(:,1)}; % phasicUnit2};
        
        pairs_mnk2{mm} = pairs;
        
    end % Monkeys.
    
    %     clear mnkRes; % Free up memory by only retaining the plot data in "mnkFigDat".
    
    %     pause
    
    %% Organize rsc plot data for "Fake beep" trials.
    
    binSiz = 1; pairs_mnk2 = [];
    
    % saveFileName = 'LC_ACC_Results_FakeBeep_rsc_042121';
    % saveFileName = 'LC_ACC_Results_FakeBeep_rsc_031721';
    % saveFileName = 'LC_ACC_Results_FakeBeep_rsc_042921c';
    % saveFileName = 'LC_ACC_Results_FakeBeep_rsc_052121b';
    saveFileName = 'LC_ACC_Results_FakeBeep_rsc_100321';
    
    load(saveFileName);
    
    for mm = 1:nmonks
        
        pairs = [];
        unit_index = 0; incindex = 0;
        ns = length(mnkRes{mm}{1}); % Number of sessions for this mnky.
        yyi = 0;
        all_loBnp = []; all_loAnp = []; all_loBp = [];  all_loAp = [];
        all_loBnp2 = []; all_loAnp2 = []; all_loBp2 = [];  all_loAp2 = [];
        all_loBnp3 = []; all_loAnp3 = []; all_loBp3 = [];  all_loAp3 = [];
        all_loBnp4 = []; all_loAnp4 = []; all_loBp4 = [];  all_loAp4 = [];
        all_lohi_count = [];
        si = 0; si2 = 0; si3 = 0; si4 = 0;
        si11 = 0; si12 = 0; si13 = 0; si14 = 0;
        pnp_dat = []; pnp_dat2 = []; pnp_dat3 = []; pnp_dat4 = [];
        phasicUnit = []; phasicUnit2 = [];
        sp_dat = []; mean_sp_dat = [];
        all_file = []; all_lci = []; all_acpi = [];
        all_file2 = []; all_lci2 = []; all_acpi2 = [];
        all_file3 = []; all_lci3 = []; all_acpi3 = [];
        all_file4 = []; all_lci4 = []; all_acpi4 = [];
        
        all_mean_B = []; all_mean_A = []; all_var_B = []; all_var_A = [];
        all_mean_B3 = []; all_mean_A3 = []; all_var_B3 = []; all_var_A3 = [];
        
        % pause
        
        for ii = 1:ns
            nu = length(mnkRes{mm}{1}{ii}{1}); % Get site data (one cell for each LC unit).
            for jj = 1:nu
                
                unit_index = unit_index+1;
                
                sp_dat = [sp_dat; zscore(mnkRes{mm}{1}{ii}{2}{jj}')'];
                pnp_dat = [pnp_dat mnkRes{mm}{1}{ii}{3}{jj}];
                phasicUnit{unit_index} = mnkRes{mm}{1}{ii}{3}{jj};
                phasicUnit2{unit_index} = mnkRes{mm}{1}{ii}{4}{jj};
                bDat = mnkRes{mm}{1}{ii}{1}{jj}; % pause
                pnp_dat3 = [pnp_dat3 mnkRes{mm}{1}{ii}{5}{jj}];
                bDat3 = mnkRes{mm}{1}{ii}{8}{jj};
                if ~isempty(bDat) & ~isempty(bDat3)
                    
                    nbins = length(bDat);
                    for kk = 1:nbins
                        
                        ei = length(bDat{kk}{1}(:,1));
                        all_loBnp(1+si:si+ei,kk) = bDat{kk}{3}(:,1);
                        all_loAnp(1+si:si+ei,kk) = bDat{kk}{4}(:,1);
                        all_loBp(1+si:si+ei,kk) = bDat{kk}{1}(:,1);
                        all_loAp(1+si:si+ei,kk) = bDat{kk}{2}(:,1);
                        
                        ei11 = length(bDat{kk}{5});
                        all_mean_B(1+si11:si11+ei11,kk) = bDat{kk}{5};
                        all_var_B(1+si11:si11+ei11,kk) = bDat{kk}{6};
                        all_mean_A(1+si11:si11+ei11,kk) = bDat{kk}{7};
                        all_var_A(1+si11:si11+ei11,kk) = bDat{kk}{8};
                        
                    end % Bins.
                    
                    all_file(1+si:si+ei,1) = bDat{1}{1}(:,3);
                    all_lci(1+si:si+ei,1) = bDat{1}{1}(:,4);
                    all_acpi(1+si:si+ei,1) = bDat{1}{1}(:,5);
                    
                    si = si+ei;
                    si11 = si11+ei11;
                    
                    %                 end % Check for empty data struct.
                    
                    nbins = length(bDat3);
                    for kk = 1:nbins
                        ei = length(bDat3{kk}{1}(:,1));
                        all_loBnp3(1+si3:si3+ei,kk) = bDat3{kk}{3}(:,1);
                        all_loAnp3(1+si3:si3+ei,kk) = bDat3{kk}{4}(:,1);
                        all_loBp3(1+si3:si3+ei,kk) = bDat3{kk}{1}(:,1);
                        all_loAp3(1+si3:si3+ei,kk) = bDat3{kk}{2}(:,1);
                        
                        ei13 = length(bDat3{kk}{5});
                        all_mean_B3(1+si13:si13+ei13,kk) = bDat3{kk}{5};
                        all_var_B3(1+si13:si13+ei13,kk) = bDat3{kk}{6};
                        all_mean_A3(1+si13:si13+ei13,kk) = bDat3{kk}{7};
                        all_var_A3(1+si13:si13+ei13,kk) = bDat3{kk}{8};
                        
                    end % Bins.
                    
                    all_file3(1+si3:si3+ei,1) = bDat3{1}{1}(:,3);
                    all_lci3(1+si3:si3+ei,1) = bDat3{1}{1}(:,4);
                    all_acpi3(1+si3:si3+ei,1) = bDat3{1}{1}(:,5);
                    
                    si3 = si3+ei;
                    si13 = si13+ei13;
                    
                end % Check for empty data struct.
                
            end % Units for this session.
        end % Sessions for this monkey.
        mnkFigDat2{mm} = {all_loBnp all_loAnp all_loBp all_loAp all_file all_lci all_acpi ...
            all_loBnp2 all_loAnp2 all_loBp2 all_loAp2 all_file2 all_lci2 all_acpi2 ...
            all_loBnp3 all_loAnp3 all_loBp3 all_loAp3 all_file3 all_lci3 all_acpi3 ...
            all_loBnp4 all_loAnp4 all_loBp4 all_loAp4 all_file4 all_lci4 all_acpi4 ...
            pnp_dat pnp_dat2 pnp_dat3 pnp_dat4 ...
            sp_dat all_mean_B all_mean_A all_var_B all_var_A all_mean_B3 all_mean_A3 all_var_B3 all_var_A3 phasicUnit phasicUnit2};
        
        pairs_mnk2{mm} = pairs;
        
    end % Monkeys.
    
    clear mnkRes; % Free up memory by only retaining the plot data in "mnkFigDat".
    
end

%% Plot stuff:

% Don't plot results?
plotYesNo = false;
% % Plot results?
plotYesNo = true;

if plotYesNo
    
    xAx = bs;
    %     limVal = 99; % Use as absolute criterion because sometimes the % difference blows up due to a infinitesimal denominator.
    
    %% BEEP: Phasic/no-phasic response.
    
    xAx = bs;
    mean_D = nans(4,length(bs));
    sem_D = nans(4,length(bs));
    
    % Phasic response.
    all_Bp1 = mnkFigDat{1}{3};   all_Bp2 = mnkFigDat{2}{3};   all_Ap1 = mnkFigDat{1}{4};   all_Ap2 = mnkFigDat{2}{4};
    % No phasic response.
    all_Bnp1 = mnkFigDat{1}{17};   all_Bnp2 = mnkFigDat{2}{17};   all_Anp1 = mnkFigDat{1}{18};   all_Anp2 = mnkFigDat{2}{18};
    
    
    all_Ap1E = mnkFigDat{1}{43};   all_Ap2E = mnkFigDat{2}{43};
    all_Ap1L = mnkFigDat{1}{44};   all_Ap2L = mnkFigDat{2}{44};
    all_Anp1E = mnkFigDat{1}{45};   all_Anp2E = mnkFigDat{2}{45};
    all_Anp1L = mnkFigDat{1}{46};   all_Anp2L = mnkFigDat{2}{46};
    
    % all_loApE(:,1) all_loApL(:,1) all_loAnpE(:,1) all_loAnpL(:,1)}; % phasicUnit2};
    %       43                  44                  45                  46
    for kk = 1:length(bs)
        
        if kk ==5
            
            gi1 = find(isfinite(all_Ap1(:,kk)) & isfinite(all_Bp1(:,kk)) & isfinite(all_Anp1(:,kk)) & isfinite(all_Bnp1(:,kk)));
            gi2 = find(isfinite(all_Ap2(:,kk)) & isfinite(all_Bp2(:,kk)) & isfinite(all_Anp2(:,kk)) & isfinite(all_Bnp2(:,kk)));
            
            D11E = all_Ap1E(gi1)-all_Bp1(gi1,kk);
            D12E = all_Ap2E(gi2)-all_Bp2(gi2,kk);
            D11nE = all_Anp1E(gi1)-all_Bnp1(gi1,kk);
            D12nE = all_Anp2E(gi2)-all_Bnp2(gi2,kk);
            DD1E = D11E-D11nE; % DD1 = DD1(abs(DD1)<prctile(abs(DD1),limVal));
            DD2E = D12E-D12nE; % DD2 = DD2(abs(DD2)<prctile(abs(DD2),limVal));
            
            gi1 = find(abs(all_Ap1L)>0 & abs(all_Ap1L)<1 & isfinite(all_Ap1L) & abs(all_Bp1(:,kk))<1 & isfinite(all_Bp1(:,kk)) & abs(all_Bp1(:,kk))>0 ...
                & abs(all_Anp1L)>0 & abs(all_Anp1L)<1 & isfinite(all_Anp1L) & abs(all_Bnp1(:,kk))<1 & isfinite(all_Bnp1(:,kk)) & abs(all_Bnp1(:,kk))>0);
            
            gi2 = find(abs(all_Ap2L)>0 & abs(all_Ap2L)<1 & isfinite(all_Ap2L) & abs(all_Bp2(:,kk)<1) & isfinite(all_Bp2(:,kk)) & abs(all_Bp2(:,kk))>0 ...
                & abs(all_Anp2L)>0 & abs(all_Anp2L)<1 & isfinite(all_Anp2L) & abs(all_Bnp2(:,kk))<1 & isfinite(all_Bnp2(:,kk)) & abs(all_Bnp2(:,kk))>0);
            
            D11L = all_Ap1L(gi1)-all_Bp1(gi1,kk);
            D12L = all_Ap2L(gi2)-all_Bp2(gi2,kk);
            D11nL = all_Anp1L(gi1)-all_Bnp1(gi1,kk);
            D12nL = all_Anp2L(gi2)-all_Bnp2(gi2,kk);
            DD1L = D11L-D11nL; % DD1 = DD1(abs(DD1)<prctile(abs(DD1),limVal));
            DD2L = D12L-D12nL; % DD2 = DD2(abs(DD2)<prctile(abs(DD2),limVal));
            
        end
        
        gi1 = find(isfinite(all_Ap1(:,kk)) & isfinite(all_Bp1(:,kk)) & isfinite(all_Anp1(:,kk)) & isfinite(all_Bnp1(:,kk)));
        gi2 = find(isfinite(all_Ap2(:,kk)) & isfinite(all_Bp2(:,kk)) & isfinite(all_Anp2(:,kk)) & isfinite(all_Bnp2(:,kk)));
        
        D11 = all_Ap1(gi1,kk)-all_Bp1(gi1,kk);
        D12 = all_Ap2(gi2,kk)-all_Bp2(gi2,kk);
        D11n = all_Anp1(gi1,kk)-all_Bnp1(gi1,kk);
        D12n = all_Anp2(gi2,kk)-all_Bnp2(gi2,kk);
        DD1 = D11-D11n; % DD1 = DD1(abs(DD1)<prctile(abs(DD1),limVal));
        DD2 = D12-D12n; % DD2 = DD2(abs(DD2)<prctile(abs(DD2),limVal));
        
        pval1p(kk) = signrank(DD1);
        pval2p(kk) = signrank(DD2);
        pval_both_np(kk) = signrank([DD1; DD2]);
        pval_both_npRS(kk) = ranksum([D11;D12], [D11n;D12n]);
        
        
        diff_p1{kk} = DD1;
        diff_p2{kk} = DD2;
        
        bef_p1{kk} = all_Bp1(gi1,kk);
        bef_p2{kk} = all_Bp2(gi2,kk);
        aft_p1{kk} = all_Ap1(gi1,kk);
        aft_p2{kk} = all_Ap2(gi2,kk);
        
        bef_np1{kk} = all_Bnp1(gi1,kk);
        bef_np2{kk} = all_Bnp2(gi2,kk);
        aft_np1{kk} = all_Anp1(gi1,kk);
        aft_np2{kk} = all_Anp2(gi2,kk);
        
    end
    
    %         pause
    
    %% Fake Beep.
    
    xAx = bs;
    mean_D = nans(4,length(bs));
    sem_D = nans(4,length(bs));
    
    % Phasic response.
    all_Bp1 = mnkFigDat2{1}{3};   all_Bp2 = mnkFigDat2{2}{3};   all_Ap1 = mnkFigDat2{1}{4};   all_Ap2 = mnkFigDat2{2}{4};
    % No phasic response.
    all_Bnp1 = mnkFigDat2{1}{17};   all_Bnp2 = mnkFigDat2{2}{17};   all_Anp1 = mnkFigDat2{1}{18};   all_Anp2 = mnkFigDat2{2}{18};
    
    for kk = 1:length(bs)
        
        gi1 = find(isfinite(all_Ap1(:,kk)) & isfinite(all_Bp1(:,kk)) & isfinite(all_Anp1(:,kk)) & isfinite(all_Bnp1(:,kk)));
        gi2 = find(isfinite(all_Ap2(:,kk)) & isfinite(all_Bp2(:,kk)) & isfinite(all_Anp2(:,kk)) & isfinite(all_Bnp2(:,kk)));
        
        D11 = all_Ap1(gi1,kk)-all_Bp1(gi1,kk);
        D12 = all_Ap2(gi2,kk)-all_Bp2(gi2,kk);
        D11n = all_Anp1(gi1,kk)-all_Bnp1(gi1,kk);
        D12n = all_Anp2(gi2,kk)-all_Bnp2(gi2,kk);
        
        %         D11 = D11t; % (abs(D11t)<prctile(abs(D11t),limVal) & abs(D11nt)<prctile(abs(D11nt),limVal));
        %         D11n = D11nt; % (abs(D11t)<prctile(abs(D11t),limVal) & abs(D11nt)<prctile(abs(D11nt),limVal));
        %         D12 = D12t; % (abs(D12t)<prctile(abs(D12t),limVal) & abs(D12nt)<prctile(abs(D12nt),limVal));
        %         D12n = D12nt; % (abs(D12t)<prctile(abs(D12t),limVal) & abs(D12nt)<prctile(abs(D12nt),limVal));
        
        DD1 = D11-D11n; % DD1 = DD1(abs(DD1)<prctile(abs(DD1),limVal));
        DD2 = D12-D12n; % DD2 = DD2(abs(DD2)<prctile(abs(DD2),limVal));
        
        %         pval1pF(kk) = signrank(DD1);
        %         pval2pF(kk) = signrank(DD2);
        pval_both_npF(kk) = signrank([DD1; DD2]);
        pval_both_npFRS(kk) = ranksum([D11;D12], [D11n;D12n]);
        
        
        %         diff_p1F{kk} = DD1;
        %         diff_p2F{kk} = DD2;
        
        bef_p1F{kk} = all_Bp1(gi1,kk);
        bef_p2F{kk} = all_Bp2(gi2,kk);
        aft_p1F{kk} = all_Ap1(gi1,kk);
        aft_p2F{kk} = all_Ap2(gi2,kk);
        
        bef_np1F{kk} = all_Bnp1(gi1,kk);
        bef_np2F{kk} = all_Bnp2(gi2,kk);
        aft_np1F{kk} = all_Anp1(gi1,kk);
        aft_np2F{kk} = all_Anp2(gi2,kk);
        
    end
    %     pause
    
    %%
    figureNumber = 8; num = 8; wid = 17.6; hts = [10]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1, 2, [12], 'Joshi&Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    xjit = 10;
    
    %% Beep
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    for kk = 1:5
        
        allBefP = [bef_p1{kk};bef_p2{kk}]; allAftP = [aft_p1{kk};aft_p2{kk}];
        diffP = allAftP - allBefP;
        allBefNP = [bef_np1{kk};bef_np2{kk}]; allAftNP = [aft_np1{kk};aft_np2{kk}];
        diffNP = allAftNP - allBefNP;
        
        pp_all(kk) = signrank(diffP);
        pnp_all(kk) = signrank(diffNP);
        %         pp1(kk) = signrank(aft_p1{kk}-bef_p1{kk});
        %         pp2(kk) = signrank(aft_p2{kk}-bef_p2{kk});
        %         pnp1(kk) = signrank(aft_np1{kk}-bef_np1{kk});
        %         pnp2(kk) = signrank(aft_np2{kk}-bef_np2{kk});
        
        mean_P = nanmedian(diffP)
        if pp_all(kk) < 0.05
            plot(1*xjit+xAx(kk),mean_P,'ko','markerfacecolor','r','markeredgecolor','r','markersize',8);
        end
        if pp_all(kk) >= 0.05
            plot(1*xjit+xAx(kk),mean_P,'ko','markerfacecolor','none','markeredgecolor','r','markersize',8);
        end
        
        mean_NP = nanmedian(diffNP)
        if pnp_all(kk) < 0.05
            plot(-1*xjit+xAx(kk),mean_NP,'ko','markerfacecolor','k','markeredgecolor','k','markersize',8);
        end
        if pnp_all(kk) >= 0.05
            plot(-1*xjit+xAx(kk),mean_NP,'ko','markerfacecolor','none','markeredgecolor','k','markersize',8);
        end
        
        if kk == 1 legend({'LC phasic','LC no phasic'},'autoupdate','off','fontsize',8,'location','northwest'); legend('boxoff'); end
        
        [ciV1 ciV_m] = bootci(1000,{@median,diffP},'alpha',0.05,'type','norm');
        plot([1*xjit+xAx(kk) 1*xjit+xAx(kk)],[ciV1(1) ciV1(2)],'-','color','r'); % 95% CI, bootstrapped.
        
        [ciV2 ciV_m] = bootci(1000,{@median,diffNP},'alpha',0.05,'type','norm');
        plot([-1*xjit+xAx(kk) -1*xjit+xAx(kk)],[ciV2(1) ciV2(2)],'-','color','k'); % 95% CI, bootstrapped.
        
        %         if pp_all(kk) < 0.05 plot(xAx(kk),0.040,'r*'); end
        %         if pnp_all(kk) < 0.05 plot(xAx(kk),0.039,'k*'); end
        %
        %         if pp1(kk) < 0.05 & pp2(kk) < 0.05 plot(xAx(kk),0.043,'r+'); end
        %         if pnp1(kk) < 0.05 & pnp2(kk) < 0.05 plot(xAx(kk),0.041,'k+'); end
        if pval_both_npRS(kk)<0.05 % & mval1 > 0 & mval2 > 0
            plot(xAx(kk),0.043,'k*');
        end
        
    end
    
    plot([100 1100],[0 0],'k--');
    axis([100 1100 -0.02 0.045]);
    set(gca,'fontname','arial','fontsize',10);
    xlabel('bin size (ms)','fontname','arial','fontsize',12);
    ylabel('ACC r_s_c','fontname','arial','fontsize',12);
    title('Beep','fontname','arial','fontsize',12,'fontweight','normal');
    
    %% Fake Beep.
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    for kk = 1:5
        
        allBefP = [bef_p1F{kk};bef_p2F{kk}]; allAftP = [aft_p1F{kk};aft_p2F{kk}];
        diffP = allAftP - allBefP;
        allBefNP = [bef_np1F{kk};bef_np2F{kk}]; allAftNP = [aft_np1F{kk};aft_np2F{kk}];
        diffNP = allAftNP - allBefNP;
        
        pp_all(kk) = signrank(diffP);
        pnp_all(kk) = signrank(diffNP);
        %         pp1(kk) = signrank(aft_p1F{kk}-bef_p1F{kk});
        %         pp2(kk) = signrank(aft_p2F{kk}-bef_p2F{kk});
        %         pnp1(kk) = signrank(aft_np1F{kk}-bef_np1F{kk});
        %         pnp2(kk) = signrank(aft_np2F{kk}-bef_np2F{kk});
        
        mean_P = nanmedian(diffP)
        if pp_all(kk) < 0.05
            plot(1*xjit+xAx(kk),mean_P,'ko','markerfacecolor','r','markeredgecolor','r','markersize',8);
        end
        if pp_all(kk) >= 0.05
            plot(1*xjit+xAx(kk),mean_P,'ko','markerfacecolor','none','markeredgecolor','r','markersize',8);
        end
        
        mean_NP = nanmedian(diffNP)
        if pnp_all(kk) < 0.05
            plot(-1*xjit+xAx(kk),mean_NP,'ko','markerfacecolor','k','markeredgecolor','k','markersize',8);
        end
        if pnp_all(kk) >= 0.05
            plot(-1*xjit+xAx(kk),mean_NP,'ko','markerfacecolor','none','markeredgecolor','k','markersize',8);
        end
        
        %         if kk == 1 legend({'before beep, LC phasic','before beep, LC no phasic','after beep, LC phasic','after beep, LC no phasic'},'autoupdate','off','fontsize',8,'location','northwest'); legend('boxoff'); end
        
        [ciV1 ciV_m] = bootci(1000,{@median,diffP},'alpha',0.05,'type','norm');
        plot([1*xjit+xAx(kk) 1*xjit+xAx(kk)],[ciV1(1) ciV1(2)],'-','color','r'); % 95% CI, bootstrapped.
        
        [ciV2 ciV_m] = bootci(1000,{@median,diffNP},'alpha',0.05,'type','norm');
        plot([-1*xjit+xAx(kk) -1*xjit+xAx(kk)],[ciV2(1) ciV2(2)],'-','color','k'); % 95% CI, bootstrapped.
        
        %         if pp_all(kk) < 0.05 plot(xAx(kk),0.040,'r*'); end
        %         if pnp_all(kk) < 0.05 plot(xAx(kk),0.039,'k*'); end
        
        %         if pp1(kk) < 0.05 & pp2(kk) < 0.05 plot(xAx(kk),0.043,'r+'); end
        %         if pnp1(kk) < 0.05 & pnp2(kk) < 0.05 plot(xAx(kk),0.041,'k+'); end
        
        if pval_both_npFRS(kk)<0.05 % & mval1 > 0 & mval2 > 0
            plot(xAx(kk),0.043,'k*');
        end
    end
    
    plot([100 1100],[0 0],'k--');
    axis([100 1100 -0.02 0.045]);
    set(gca,'fontname','arial','fontsize',10);
    xlabel('bin size (ms)','fontname','arial','fontsize',12);
    ylabel('ACC r_s_c','fontname','arial','fontsize',12);
    title('Fake Beep','fontname','arial','fontsize',12,'fontweight','normal');
    
    %%
    
    %     pause
    
    figureNumber = 82; num = 82; wid = 17.6; hts = [12]; cols = {1}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1, 1, [12], 'Joshi&Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    xjit = 10;
    
    %%
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    for kk = 1:5
        
        allBefP = [bef_p1{kk};bef_p2{kk}]; allAftP = [aft_p1{kk};aft_p2{kk}];
        allBefNP = [bef_np1{kk};bef_np2{kk}]; allAftNP = [aft_np1{kk};aft_np2{kk}];
        
        %         pp_all(kk) = ranksum(allBefP,allBefNP);
        %         ppB1(kk) = ranksum(bef_p1{kk},bef_np1{kk});
        %         ppB2(kk) = ranksum(bef_p2{kk},bef_np2{kk});
        
        pB_all(kk) = ranksum(allBefP,allBefNP);
        pA_all(kk) = ranksum(allAftP,allAftNP);
        pB1(kk) = ranksum(bef_p1{kk},bef_np1{kk});
        pB2(kk) = ranksum(bef_p2{kk},bef_np2{kk});
        pA1(kk) = ranksum(aft_p1{kk},aft_np1{kk});
        pA2(kk) = ranksum(aft_p2{kk},aft_np2{kk});
        
        mean_BefP = nanmedian(allBefP);
        plot(-3*xjit+xAx(kk),mean_BefP,'ro','markerfacecolor','none','markeredgecolor','r','markersize',7);
        
        mean_BefNP = nanmedian(allBefNP);
        plot(-1*xjit+xAx(kk),mean_BefNP,'ko','markerfacecolor','none','markeredgecolor','k','markersize',7);
        
        mean_AftP = nanmedian(allAftP);
        plot(1*xjit+xAx(kk),mean_AftP,'ro','markerfacecolor','r','markeredgecolor','none','markersize',7);
        
        mean_AftNP = nanmedian(allAftNP);
        plot(3*xjit+xAx(kk),mean_AftNP,'ko','markerfacecolor','k','markeredgecolor','none','markersize',7);
        
        if kk == 1 legend({'before beep, LC phasic','before beep, LC no phasic','after beep, LC phasic','after beep, LC no phasic'},'autoupdate','off','fontsize',8,'location','northwest'); legend('boxoff'); end
        
        [ciV1 ciV_m] = bootci(1000,{@median,allBefP},'alpha',0.05,'type','norm');
        plot([-3*xjit+xAx(kk) -3*xjit+xAx(kk)],[ciV1(1) ciV1(2)],'-','color','r'); % 95% CI, bootstrapped.
        
        [ciV2 ciV_m] = bootci(1000,{@median,allBefNP},'alpha',0.05,'type','norm');
        plot([-1*xjit+xAx(kk) -1*xjit+xAx(kk)],[ciV2(1) ciV2(2)],'-','color','k'); % 95% CI, bootstrapped.
        
        [ciV3 ciV_m] = bootci(1000,{@median,allAftP},'alpha',0.05,'type','norm');
        plot([1*xjit+xAx(kk) 1*xjit+xAx(kk)],[ciV3(1) ciV3(2)],'-','color','r'); % 95% CI, bootstrapped.
        
        [ciV4 ciV_m] = bootci(1000,{@median,allAftNP},'alpha',0.05,'type','norm');
        plot([3*xjit+xAx(kk) 3*xjit+xAx(kk)],[ciV4(1) ciV4(2)],'-','color','k'); % 95% CI, bootstrapped.
        
        if pB_all(kk) < 0.05 plot(xAx(kk),0.057,'k^'); end
        if pA_all(kk) < 0.05 plot(xAx(kk),0.055,'kd'); end
        
        if pB1(kk) < 0.05 & pB2(kk) < 0.05 plot(xAx(kk),0.061,'k+'); end
        if pA1(kk) < 0.05 & pA2(kk) < 0.05 plot(xAx(kk),0.059,'m+'); end
        
    end
    
    axis([100 1100 0.01 0.061]);
    set(gca,'fontname','arial','fontsize',10);
    xlabel('bin size (ms)','fontname','arial','fontsize',12);
    ylabel('ACC r_s_c','fontname','arial','fontsize',12);
    title('Beep','fontname','arial','fontsize',12,'fontweight','normal');
    
    % pause
    
    %% Setup SUPP figure.
    
    figureNumber = 83; num = 83; wid = 17.6; hts = [10]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 3, 1, [12], 'Joshi&Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    
    
    %% Plot Beep.
    
    xjit = 10;
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    ydE1p = [D11E; D12E]; ydE1p = ydE1p(isfinite(ydE1p)); mean_ydE1p = nanmedian(ydE1p);
    ydE1np = [D11nE; D12nE]; ydE1np = ydE1np(isfinite(ydE1np)); mean_ydE1np = nanmedian(ydE1np);
    plot(1,mean_ydE1p,'ro','markerfacecolor','r','markeredgecolor','none','markersize',10);
    plot(2,mean_ydE1np,'ko','markerfacecolor','k','markeredgecolor','none','markersize',10);
    [ciVp ciV_m] = bootci(1000,{@median,ydE1p},'alpha',0.05,'type','norm');
    plot([1 1],[ciVp(1) ciVp(2)],'-','color','r'); % 95% CI, bootstrapped.
    [ciVnp ciV_m] = bootci(1000,{@median,ydE1np},'alpha',0.05,'type','norm');
    plot([2 2],[ciVnp(1) ciVnp(2)],'-','color','k'); % 95% CI, bootstrapped.
    
    ydL1p = [D11L; D12L]; ydL1p = ydL1p(isfinite(ydL1p)); mean_ydL1p = nanmedian(ydL1p);
    ydL1np = [D11nL; D12nL]; ydL1np = ydL1np(isfinite(ydL1np)); mean_ydL1np = nanmedian(ydL1np);
    plot(5,mean_ydL1p,'ro','markerfacecolor','r','markeredgecolor','none','markersize',10);
    plot(6,mean_ydL1np,'ko','markerfacecolor','k','markeredgecolor','none','markersize',10);
    [ciVp ciV_m] = bootci(1000,{@median,ydL1p},'alpha',0.05,'type','norm');
    plot([5 5],[ciVp(1) ciVp(2)],'-','color','r'); % 95% CI, bootstrapped.
    [ciVnp ciV_m] = bootci(1000,{@median,ydL1np},'alpha',0.05,'type','norm');
    plot([6 6],[ciVnp(1) ciVnp(2)],'-','color','k'); % 95% CI, bootstrapped.
    
    p1 = ranksum(ydE1p,ydL1p);
    p2 = ranksum(ydE1np,ydL1np);
    
    p11 = ranksum(D11E,D11L);
    p12 = ranksum(D12E,D12L);
    
    p11n = ranksum(D11nE,D11nL);
    p12n = ranksum(D12nE,D12nL);
    
    p1pnp = ranksum(ydE1p,ydE1np);
    p2pnp = ranksum(ydL1p,ydL1np);
    
    p11pnp = ranksum(D11L,D11nL);
    p12pnp = ranksum(D12L,D12nL);
    p21pnp = ranksum(D11E,D11nE);
    p22pnp = ranksum(D12E,D12nE);
    
    %     p2 = ranksum(ydE2,ydL2);
    %     p1E = signrank(DD1E);
    %     p1L = signrank(DD1L);
    %     p2E = signrank(DD2E);
    %     p2L = signrank(DD2L);
    
    %     pE = signrank([ydE1]);%;ydE2]);
    %     pL = signrank([ydL1]);%;ydE2]);
    
    if p1 < 0.05 plot(5,0.044,'r*','markersize',8); end
    if p2 < 0.05 plot(6,0.044,'k*','markersize',8); end
    if p11 < 0.05 & p12 < 0.05 plot(5,0.046,'r+','markersize',8); end
    if p11n < 0.05 & p12n < 0.05 plot(6,0.046,'k+','markersize',8); end
    
    ylabel({'change in beep linked r_s_c';'after re:before beep'},'fontsize',8,'fontname','arial');
    title('Beep trials','fontsize',10,'fontname','arial');
    plot([0 7],[0 0],'k--');
    axis([0 7 -0.046 0.046]);
    set(gca,'fontname','arial','fontsize',10);
    set(gca,'xtick',[1.5 5.5]);
    % set(gca,'xtick',[1.5 5.5]);
    set(gca,'xticklabel',{'early','late'});
    xlabel('bin size (ms)');
    
    %%
    figureNumber = 84; num = 84; wid = 17.6; hts = [10]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1, 2, [12], 'Joshi&Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    xjit = 10;
    
    %%
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    for kk = 1:5
        
        allBefP = [bef_p1{kk};bef_p2{kk}]; allAftP = [aft_p1{kk};aft_p2{kk}];
        allBefNP = [bef_np1{kk};bef_np2{kk}]; allAftNP = [aft_np1{kk};aft_np2{kk}];
        diffP = allAftP - allAftNP;
        diffNP = allBefP - allBefNP;
        
        pp_all(kk) = signrank(diffP);
        pnp_all(kk) = signrank(diffNP);
        pp1(kk) = signrank(aft_p1{kk}-aft_np1{kk});
        pp2(kk) = signrank(aft_p2{kk}-aft_np2{kk});
        pnp1(kk) = signrank(bef_p1{kk}-bef_np1{kk});
        pnp2(kk) = signrank(bef_p2{kk}-bef_np2{kk});
        
        mean_P = nanmedian(diffP)
        plot(1*xjit+xAx(kk),mean_P,'ko','markerfacecolor','r','markeredgecolor','none','markersize',8);
        
        mean_NP = nanmedian(diffNP)
        plot(-1*xjit+xAx(kk),mean_NP,'ko','markerfacecolor','k','markeredgecolor','none','markersize',8);
        
        %         if kk == 1 legend({'before beep, LC phasic','before beep, LC no phasic','after beep, LC phasic','after beep, LC no phasic'},'autoupdate','off','fontsize',8,'location','northwest'); legend('boxoff'); end
        
        [ciV1 ciV_m] = bootci(1000,{@median,diffP},'alpha',0.05,'type','norm');
        plot([1*xjit+xAx(kk) 1*xjit+xAx(kk)],[ciV1(1) ciV1(2)],'-','color','r'); % 95% CI, bootstrapped.
        
        [ciV2 ciV_m] = bootci(1000,{@median,diffNP},'alpha',0.05,'type','norm');
        plot([-1*xjit+xAx(kk) -1*xjit+xAx(kk)],[ciV2(1) ciV2(2)],'-','color','k'); % 95% CI, bootstrapped.
        
        if pp_all(kk) < 0.05 plot(1*xjit+xAx(kk),0.021,'r*'); end
        if pnp_all(kk) < 0.05 plot(-1*xjit+xAx(kk),0.020,'k*'); end
        
        if pp1(kk) < 0.05 & pp2(kk) < 0.05 plot(1*xjit+xAx(kk),0.023,'r+'); end
        if pnp1(kk) < 0.05 & pnp2(kk) < 0.05 plot(-1*xjit+xAx(kk),0.022,'k+'); end
        
    end
    
    axis([100 1100 -0.029 0.025]);
    set(gca,'fontname','arial','fontsize',10);
    xlabel('bin size (ms)','fontname','arial','fontsize',12);
    ylabel('ACC r_s_c','fontname','arial','fontsize',12);
    title('Beep','fontname','arial','fontsize',12,'fontweight','normal');
    
    %%
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    for kk = 1:5
        
        allBefP = [bef_p1F{kk};bef_p2F{kk}]; allAftP = [aft_p1F{kk};aft_p2F{kk}];
        diffP = allAftP - allBefP;
        allBefNP = [bef_np1F{kk};bef_np2F{kk}]; allAftNP = [aft_np1F{kk};aft_np2F{kk}];
        diffNP = allAftNP - allBefNP;
        
        pp_all(kk) = signrank(diffP);
        pnp_all(kk) = signrank(diffNP);
        pp1(kk) = signrank(aft_p1F{kk}-bef_p1F{kk});
        pp2(kk) = signrank(aft_p2F{kk}-bef_p2F{kk});
        pnp1(kk) = signrank(aft_np1F{kk}-bef_np1F{kk});
        pnp2(kk) = signrank(aft_np2F{kk}-bef_np2F{kk});
        
        mean_P = nanmedian(diffP)
        plot(1*xjit+xAx(kk),mean_P,'ko','markerfacecolor','r','markeredgecolor','none','markersize',8);
        
        mean_NP = nanmedian(diffNP)
        plot(-1*xjit+xAx(kk),mean_NP,'ko','markerfacecolor','k','markeredgecolor','none','markersize',8);
        
        %         if kk == 1 legend({'before beep, LC phasic','before beep, LC no phasic','after beep, LC phasic','after beep, LC no phasic'},'autoupdate','off','fontsize',8,'location','northwest'); legend('boxoff'); end
        
        [ciV1 ciV_m] = bootci(1000,{@median,diffP},'alpha',0.05,'type','norm');
        plot([1*xjit+xAx(kk) 1*xjit+xAx(kk)],[ciV1(1) ciV1(2)],'-','color','r'); % 95% CI, bootstrapped.
        
        [ciV2 ciV_m] = bootci(1000,{@median,diffNP},'alpha',0.05,'type','norm');
        plot([-1*xjit+xAx(kk) -1*xjit+xAx(kk)],[ciV2(1) ciV2(2)],'-','color','k'); % 95% CI, bootstrapped.
        
        if pp_all(kk) < 0.05 plot(xAx(kk),0.040,'r*'); end
        if pnp_all(kk) < 0.05 plot(xAx(kk),0.039,'k*'); end
        
        if pp1(kk) < 0.05 & pp2(kk) < 0.05 plot(xAx(kk),0.043,'r+'); end
        if pnp1(kk) < 0.05 & pnp2(kk) < 0.05 plot(xAx(kk),0.041,'k+'); end
        
    end
    
    axis([100 1100 -0.02 0.045]);
    set(gca,'fontname','arial','fontsize',10);
    xlabel('bin size (ms)','fontname','arial','fontsize',12);
    ylabel('ACC r_s_c','fontname','arial','fontsize',12);
    title('Fake Beep','fontname','arial','fontsize',12,'fontweight','normal');
    
    %% Setup OLD Main figure.
    
    %     figureNumber = 888; num = 888; wid = 17.6; hts = [7]; cols = {2 2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 3, 1, [12], 'Joshi&Gold, 2020', true,1,figureNumber); set(axs,'Units','normalized');
    %     movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    %     ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    %
    %     colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    %
    %
    %     %% Plot Beep.
    %
    %     xjit = 25;
    %     axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %
    %     for kk = 1:length(bs)
    %
    %         yd = [diff_p1{kk};diff_p2{kk}]; yd = yd(isfinite(yd)); mean_yd = nanmedian(yd);
    %         [ciV ciV_m] = bootci(1000,{@median,yd},'alpha',0.05,'type','norm');
    %         plot(xAx(kk),mean_yd,'ko','markerfacecolor','k','markeredgecolor','none');
    %         plot([xAx(kk) xAx(kk)],[ciV(1) ciV(2)],'-','markersize',20,'color','k'); % 95% CI, bootstrapped.
    %
    %         %         yd1 = diff_p1{kk}; yd1 = yd1(isfinite(yd1)); mean_yd1 = nanmedian(yd1);
    %         %         [ciV1 ciV_m] = bootci(1000,{@median,yd1},'alpha',0.05,'type','norm'); % mean_ze = nanmean(ci_ze_m);
    %         %         bar(xAx(kk)-xjit,mean_yd1,2*xjit,'facecolor',colors{3},'edgecolor','none');
    %         %         plot([xAx(kk)-xjit xAx(kk)-xjit],[ciV1(1) ciV1(2)],'-','markersize',20,'color','k'); % 95% CI, bootstrapped.
    %         %
    %         %         yd2 = diff_p2{kk}; yd2 = yd2(isfinite(yd2)); mean_yd2 = nanmedian(yd2);
    %         %         [ciV2 ciV_m] = bootci(1000,{@median,yd2},'alpha',0.05,'type','norm'); % mean_ze = nanmean(ci_ze_m);
    %         %         bar(xAx(kk)+xjit,mean_yd2,2*xjit,'facecolor',colors{4},'edgecolor','none');
    %         %         plot([xAx(kk)+xjit xAx(kk)+xjit],[ciV2(1) ciV2(2)],'-','markersize',20,'color','k'); % 95% CI, bootstrapped.
    %
    %         if pval1p(kk)<0.05 & pval2p(kk)<0.05 % & mval1 > 0 & mval2 > 0
    %             plot(xAx(kk),0.06,'k+');
    %         end
    %         if pval_both_np(kk)<0.05 % & mval1 > 0 & mval2 > 0
    %             plot(xAx(kk),0.055,'k*');
    %         end
    %
    %     end
    %     ylabel({'change in beep linked r_s_c';'phasic re:no phasic LC response'},'fontsize',8,'fontname','arial');
    %     title('Beep trials','fontsize',10,'fontname','arial');
    %     plot([xAx(1) xAx(end)],[0 0],'k--');
    %     axis([100 1100 -0.005 0.06]);
    %     set(gca,'fontname','arial','fontsize',10);
    %     xlabel('bin size (ms');
    %
    %
    %     %% Plot FAKE Beep.
    %
    %     axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    %
    %     for kk = 1:length(bs)
    %
    %         yd = [diff_p1F{kk};diff_p2F{kk}]; yd = yd(isfinite(yd)); mean_yd = nanmedian(yd);
    %         [ciV ciV_m] = bootci(1000,{@median,yd},'alpha',0.05,'type','norm'); % mean_ze = nanmean(ci_ze_m);
    %         % bar(xAx(kk),mean_yd1,2*xjit,'facecolor',colors{3},'edgecolor','none');
    %         plot(xAx(kk),mean_yd,'ko','markerfacecolor','k','markeredgecolor','none');
    %         plot([xAx(kk) xAx(kk)],[ciV(1) ciV(2)],'-','markersize',20,'color','k'); % 95% CI, bootstrapped.
    %
    %         if pval1pF(kk)<0.05 & pval2pF(kk)<0.05 % & mval1 > 0 & mval2 > 0
    %             plot(xjit+xAx(kk),0.06,'k+');
    %         end
    %         if pval_both_npF(kk) < 0.05 % & mval1 > 0 & mval2 > 0
    %             plot(xjit+xAx(kk),0.055,'k*');
    %         end
    %
    %     end
    %     % ylabel({'change in beep linked r_s_c';'phasic re:no phasic LC response'},'fontsize',10,'fontname','arial');
    %     title('FAKE beep','fontsize',10,'fontname','arial');
    %     plot([xAx(1) xAx(end)],[0 0],'k--');
    %     axis([100 1100 -0.005 0.06]);
    %     set(gca,'fontname','arial','fontsize',10);
    %     xlabel('bin size (ms');
    %
    %     pause
    
    
    %%
    figureNumber = 85; num = 85; wid = 17.6; hts = [10]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1, 2, [12], 'Joshi&Gold, 2021', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    colors = {'m', 0.75.*ones(1,3), 0.5.*ones(1,3), 0.25.*ones(1,3), zeros(1,3)};
    xjit = 10;
    
    %% Fano.
    
    Fano_b1 = []; Fano_b2 = []; Fano_a1 = []; Fano_a2 = [];
    Fano_b1n = []; Fano_b2n = []; Fano_a1n = []; Fano_a2n = [];
    Fano_diff1 =  []; Fano_diff2 =  [];
    Fano_diff1n =  []; Fano_diff2n =  [];
    
    mb1 = mnkFigDat{1}{34}; ma1 = mnkFigDat{1}{35}; vb1 = mnkFigDat{1}{36}; va1 = mnkFigDat{1}{37};
    mb1n = mnkFigDat{1}{38}; ma1n = mnkFigDat{1}{39}; vb1n = mnkFigDat{1}{40}; va1n = mnkFigDat{1}{41};
    mb2 = mnkFigDat{2}{34}; ma2 = mnkFigDat{2}{35}; vb2 = mnkFigDat{2}{36}; va2 = mnkFigDat{2}{37};
    mb2n = mnkFigDat{2}{38}; ma2n = mnkFigDat{2}{39}; vb2n = mnkFigDat{2}{40}; va2n = mnkFigDat{2}{41};
    Fano_se = [];
    Fano_seN = [];
    
    for kk = 1:length(bs)
        
        FFb1 = vb1(:,kk)./mb1(:,kk); FFa1 = va1(:,kk)./ma1(:,kk);
        FFb2 = vb2(:,kk)./mb2(:,kk); FFa2 = va2(:,kk)./ma2(:,kk);
        
        FFb1n = vb1n(:,kk)./mb1n(:,kk); FFa1n = va1n(:,kk)./ma1n(:,kk);
        FFb2n = vb2n(:,kk)./mb2n(:,kk); FFa2n = va2n(:,kk)./ma2n(:,kk);
        
        %         gi1 = find(isfinite(FFb1) & isfinite(FFa1) & FFb1>0 & FFa1>0 & FFb1<10 & FFa1<20 ...
        %             & isfinite(FFb1n) & isfinite(FFa1n) & FFb1n>0 & FFa1n>0 & FFb1n<10 & FFa1n<20);
        %         gi2 = find(isfinite(FFb2) & isfinite(FFa2) & FFb2>0 & FFa2>0 & FFb2<10 & FFa2<20 ...
        %             & isfinite(FFb2n) & isfinite(FFa2n) & FFb2n>0 & FFa2n>0 & FFb2n<10 & FFa2n<20);
        
        %         gi1 = find(isfinite(FFb1n) & isfinite(FFa1n) & FFb1n>0 & FFa1n>0 ...
        %             & isfinite(FFb1) & isfinite(FFa1) & FFb1>0 & FFa1>0);
        %         gi2 = find(isfinite(FFb2n) & isfinite(FFa2n) & FFb2n>0 & FFa2n>0 ...
        %             & isfinite(FFb2) & isfinite(FFa2) & FFb2>0 & FFa2>0);
        
        gi1 = find(isfinite(FFb1n) & isfinite(FFa1n) & isfinite(FFb1) & isfinite(FFa1));
        gi2 = find(isfinite(FFb2n) & isfinite(FFa2n) & isfinite(FFb2) & isfinite(FFa2));
        
        D11 = FFa1(gi1)-FFb1(gi1);
        D11n = FFa1n(gi1)-FFb1n(gi1);
        D12 = FFa2(gi2)-FFb2(gi2);
        D12n = FFa2n(gi2)-FFb2n(gi2);
        
        DD1 = D11-D11n;
        DD2 = D12-D12n;
        %         DD1 = DD1(abs(DD1)<prctile(abs(DD1),limVal));
        %         DD2 = DD2(abs(DD2)<prctile(abs(DD2),limVal));
        
        %         pval1p_Fano(kk) = signrank(DD1);
        %         pval2p_Fano(kk) = signrank(DD2);
        pval_both_np_Fano(kk) = signrank([DD1; DD2]);
        pval_both_npRS(kk) = ranksum([D11;D12], [D11n;D12n]);
        
        %         diff_p1_Fano{kk} = DD1;
        %         diff_p2_Fano{kk} = DD2;
        
        bef_p1F{kk} = FFb1(gi1);
        bef_p2F{kk} = FFb2(gi2);
        aft_p1F{kk} = FFa1(gi1);
        aft_p2F{kk} = FFa2(gi2);
        
        bef_p1Fn{kk} = FFb1n(gi1);
        bef_p2Fn{kk} = FFb2n(gi2);
        aft_p1Fn{kk} = FFa1n(gi1);
        aft_p2Fn{kk} = FFa2n(gi2);
        
    end
    
    %% Plot Beep Fano.
    
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    for kk = 1:5
        
        allBefP = [bef_p1F{kk};bef_p2F{kk}]; allAftP = [aft_p1F{kk};aft_p2F{kk}];
        diffP = allAftP - allBefP;
        allBefNP = [bef_p1Fn{kk};bef_p2Fn{kk}]; allAftNP = [aft_p1Fn{kk};aft_p2Fn{kk}];
        diffNP = allAftNP - allBefNP;
        
        pp_all(kk) = signrank(diffP);
        pnp_all(kk) = signrank(diffNP);
        %         pp1(kk) = signrank(aft_p1{kk}-bef_p1{kk});
        %         pp2(kk) = signrank(aft_p2{kk}-bef_p2{kk});
        %         pnp1(kk) = signrank(aft_np1{kk}-bef_np1{kk});
        %         pnp2(kk) = signrank(aft_np2{kk}-bef_np2{kk});
        
        mean_P = nanmedian(diffP)
        if pp_all(kk) < 0.05
            plot(1*xjit+xAx(kk),mean_P,'ko','markerfacecolor','r','markeredgecolor','r','markersize',8);
        end
        if pp_all(kk) >= 0.05
            plot(1*xjit+xAx(kk),mean_P,'ko','markerfacecolor','none','markeredgecolor','r','markersize',8);
        end
        
        mean_NP = nanmedian(diffNP)
        if pnp_all(kk) < 0.05
            plot(-1*xjit+xAx(kk),mean_NP,'ko','markerfacecolor','k','markeredgecolor','r','markersize',8);
        end
        if pnp_all(kk) >= 0.05
            plot(-1*xjit+xAx(kk),mean_NP,'ko','markerfacecolor','none','markeredgecolor','k','markersize',8);
        end
        
        if kk == 1 legend({'LC phasic','LC no phasic'},'autoupdate','off','fontsize',8,'location','southwest'); legend('boxoff'); end
        
        [ciV1 ciV_m] = bootci(1000,{@median,diffP},'alpha',0.05,'type','norm');
        plot([1*xjit+xAx(kk) 1*xjit+xAx(kk)],[ciV1(1) ciV1(2)],'-','color','r'); % 95% CI, bootstrapped.
        
        [ciV2 ciV_m] = bootci(1000,{@median,diffNP},'alpha',0.05,'type','norm');
        plot([-1*xjit+xAx(kk) -1*xjit+xAx(kk)],[ciV2(1) ciV2(2)],'-','color','k'); % 95% CI, bootstrapped.
        
        %         if pp_all(kk) < 0.05 plot(xAx(kk),0.040,'r*'); end
        %         if pnp_all(kk) < 0.05 plot(xAx(kk),0.039,'k*'); end
        %
        %         if pp1(kk) < 0.05 & pp2(kk) < 0.05 plot(xAx(kk),0.043,'r+'); end
        %         if pnp1(kk) < 0.05 & pnp2(kk) < 0.05 plot(xAx(kk),0.041,'k+'); end
        if pval_both_npRS(kk)<0.05 % & mval1 > 0 & mval2 > 0
            plot(xAx(kk),0.02,'k*');
        end
        
    end
    
    plot([100 1100],[0 0],'k--');
    axis([100 1100 -0.14 0.029]);
    set(gca,'fontname','arial','fontsize',10);
    xlabel('bin size (ms)','fontname','arial','fontsize',12);
    ylabel('ACC Fano factor','fontname','arial','fontsize',12);
    title('Beep','fontname','arial','fontsize',12,'fontweight','normal');
    
    %% Fake beep Fano.
    
    Fano_b1 = []; Fano_b2 = []; Fano_a1 = []; Fano_a2 = [];
    Fano_b1n = []; Fano_b2n = []; Fano_a1n = []; Fano_a2n = [];
    Fano_diff1 =  []; Fano_diff2 =  [];
    Fano_diff1n =  []; Fano_diff2n =  [];
    
    mb1 = mnkFigDat2{1}{34}; ma1 = mnkFigDat2{1}{35}; vb1 = mnkFigDat2{1}{36}; va1 = mnkFigDat2{1}{37};
    mb1n = mnkFigDat2{1}{38}; ma1n = mnkFigDat2{1}{39}; vb1n = mnkFigDat2{1}{40}; va1n = mnkFigDat2{1}{41};
    mb2 = mnkFigDat2{2}{34}; ma2 = mnkFigDat2{2}{35}; vb2 = mnkFigDat2{2}{36}; va2 = mnkFigDat2{2}{37};
    mb2n = mnkFigDat2{2}{38}; ma2n = mnkFigDat2{2}{39}; vb2n = mnkFigDat2{2}{40}; va2n = mnkFigDat2{2}{41};
    Fano_se = [];
    Fano_seN = [];
    
    for kk = 1:length(bs)
        
        FFb1 = vb1(:,kk)./mb1(:,kk); FFa1 = va1(:,kk)./ma1(:,kk);
        FFb2 = vb2(:,kk)./mb2(:,kk); FFa2 = va2(:,kk)./ma2(:,kk);
        
        FFb1n = vb1n(:,kk)./mb1n(:,kk); FFa1n = va1n(:,kk)./ma1n(:,kk);
        FFb2n = vb2n(:,kk)./mb2n(:,kk); FFa2n = va2n(:,kk)./ma2n(:,kk);
        
        %         gi1 = find(isfinite(FFb1) & isfinite(FFa1) & FFb1>0 & FFa1>0 & FFb1<10 & FFa1<20 ...
        %             & isfinite(FFb1n) & isfinite(FFa1n) & FFb1n>0 & FFa1n>0 & FFb1n<10 & FFa1n<20);
        %         gi2 = find(isfinite(FFb2) & isfinite(FFa2) & FFb2>0 & FFa2>0 & FFb2<10 & FFa2<20 ...
        %             & isfinite(FFb2n) & isfinite(FFa2n) & FFb2n>0 & FFa2n>0 & FFb2n<10 & FFa2n<20);
        
        %         gi1 = find(isfinite(FFb1n) & isfinite(FFa1n) & FFb1n>0 & FFa1n>0 ...
        %             & isfinite(FFb1) & isfinite(FFa1) & FFb1>0 & FFa1>0);
        %         gi2 = find(isfinite(FFb2n) & isfinite(FFa2n) & FFb2n>0 & FFa2n>0 ...
        %             & isfinite(FFb2) & isfinite(FFa2) & FFb2>0 & FFa2>0);
        
        gi1 = find(isfinite(FFb1n) & isfinite(FFa1n) & isfinite(FFb1) & isfinite(FFa1));
        gi2 = find(isfinite(FFb2n) & isfinite(FFa2n) & isfinite(FFb2) & isfinite(FFa2));
        
        D11 = FFa1(gi1)-FFb1(gi1);
        D11n = FFa1n(gi1)-FFb1n(gi1);
        D12 = FFa2(gi2)-FFb2(gi2);
        D12n = FFa2n(gi2)-FFb2n(gi2);
        
        DD1 = D11-D11n;
        DD2 = D12-D12n;
        %         DD1 = DD1(abs(DD1)<prctile(abs(DD1),limVal));
        %         DD2 = DD2(abs(DD2)<prctile(abs(DD2),limVal));
        
        %         pval1p_Fano(kk) = signrank(DD1);
        %         pval2p_Fano(kk) = signrank(DD2);
        pval_both_np_Fano(kk) = signrank([DD1; DD2]);
        pval_both_npRS(kk) = ranksum([D11;D12], [D11n;D12n]);
        
        %         diff_p1_Fano{kk} = DD1;
        %         diff_p2_Fano{kk} = DD2;
        
        bef_p1F{kk} = FFb1(gi1);
        bef_p2F{kk} = FFb2(gi2);
        aft_p1F{kk} = FFa1(gi1);
        aft_p2F{kk} = FFa2(gi2);
        
        bef_p1Fn{kk} = FFb1n(gi1);
        bef_p2Fn{kk} = FFb2n(gi2);
        aft_p1Fn{kk} = FFa1n(gi1);
        aft_p2Fn{kk} = FFa2n(gi2);
        
    end
    
    %% Plot Fake beep Fano.
    
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    
    for kk = 1:5
        
        allBefP = [bef_p1F{kk};bef_p2F{kk}]; allAftP = [aft_p1F{kk};aft_p2F{kk}];
        diffP = allAftP - allBefP;
        allBefNP = [bef_p1Fn{kk};bef_p2Fn{kk}]; allAftNP = [aft_p1Fn{kk};aft_p2Fn{kk}];
        diffNP = allAftNP - allBefNP;
        
        pp_all(kk) = signrank(diffP);
        pnp_all(kk) = signrank(diffNP);
        %         pp1(kk) = signrank(aft_p1{kk}-bef_p1{kk});
        %         pp2(kk) = signrank(aft_p2{kk}-bef_p2{kk});
        %         pnp1(kk) = signrank(aft_np1{kk}-bef_np1{kk});
        %         pnp2(kk) = signrank(aft_np2{kk}-bef_np2{kk});
        
        mean_P = nanmedian(diffP)
        if pp_all(kk) < 0.05
            plot(1*xjit+xAx(kk),mean_P,'ko','markerfacecolor','r','markeredgecolor','r','markersize',8);
        end
        if pp_all(kk) >= 0.05
            plot(1*xjit+xAx(kk),mean_P,'ko','markerfacecolor','none','markeredgecolor','r','markersize',8);
        end
        
        mean_NP = nanmedian(diffNP)
        if pnp_all(kk) < 0.05
            plot(-1*xjit+xAx(kk),mean_NP,'ko','markerfacecolor','k','markeredgecolor','k','markersize',8);
        end
        if pnp_all(kk) >= 0.05
            plot(-1*xjit+xAx(kk),mean_NP,'ko','markerfacecolor','none','markeredgecolor','k','markersize',8);
        end
        
        if kk == 1 legend({'LC phasic','LC no phasic'},'autoupdate','off','fontsize',8,'location','southwest'); legend('boxoff'); end
        
        [ciV1 ciV_m] = bootci(1000,{@median,diffP},'alpha',0.05,'type','norm');
        plot([1*xjit+xAx(kk) 1*xjit+xAx(kk)],[ciV1(1) ciV1(2)],'-','color','r'); % 95% CI, bootstrapped.
        
        [ciV2 ciV_m] = bootci(1000,{@median,diffNP},'alpha',0.05,'type','norm');
        plot([-1*xjit+xAx(kk) -1*xjit+xAx(kk)],[ciV2(1) ciV2(2)],'-','color','k'); % 95% CI, bootstrapped.
        
        %         if pp_all(kk) < 0.05 plot(xAx(kk),0.040,'r*'); end
        %         if pnp_all(kk) < 0.05 plot(xAx(kk),0.039,'k*'); end
        %
        %         if pp1(kk) < 0.05 & pp2(kk) < 0.05 plot(xAx(kk),0.043,'r+'); end
        %         if pnp1(kk) < 0.05 & pnp2(kk) < 0.05 plot(xAx(kk),0.041,'k+'); end
        if pval_both_npFRS(kk)<0.05 % & mval1 > 0 & mval2 > 0
            plot(xAx(kk),0.02,'k*');
        end
        
    end
    
    plot([100 1100],[0 0],'k--');
    axis([100 1100 -0.14 0.029]);
    set(gca,'fontname','arial','fontsize',10);
    xlabel('bin size (ms)','fontname','arial','fontsize',12);
    ylabel('ACC Fano factor','fontname','arial','fontsize',12);
    title('FakeBeep','fontname','arial','fontsize',12,'fontweight','normal');
    
    
    %% LC phasic proportions.
    
    figureNumber = 81; num = 81; wid = 17.6; hts = [8]; cols = {2}; [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1, 2, [12], 'Joshi&Gold, 2018', true,1,figureNumber); set(axs,'Units','normalized');
    movegui(fig_,[1,1]); % Move figure so it is visible on screen.
    ax = gca; removeToolbarExplorationButtons(ax); disableDefaultInteractivity(ax);
    
    %%
    
    pud1 = mnkFigDat{1}{42};
    for kk = 1:length(pud1)
        pnp1 = pud1{kk};
        prop_pnp1(kk) = pnp1;
        % prop_pnp1(kk) = sum(isfinite(pnp))/length(pnp);
    end
    
    pud2 = mnkFigDat{2}{42};
    for kk = 1:length(pud2)
        pnp2 = pud2{kk};
        % prop_pnp2(kk) = sum(isfinite(pnp2))/length(pnp2);
        prop_pnp2(kk) = pnp2;
    end
    
    all_pnp = [prop_pnp1 prop_pnp2];
    xBins = linspace(0,1,11);
    % xBins = linspace(min(all_pnp),max(all_pnp),15);
    n1 = hist(prop_pnp1,xBins);
    n2 = hist(prop_pnp2,xBins);
    
    yPlotLim = 0.31;
    
    med1 = nanmedian(prop_pnp1); % 0.5413
    axes(axs(1)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    % plot(xBins,n1/sum(n1),'k','linewidth',1);
    bar(xBins,n1/sum(n1),'k','edgecolor','none','facecolor','k');
    plot([med1 med1],[0 yPlotLim],'r--','linewidth',2);
    
    axis([-0.1 1.1 -0.005 yPlotLim]);
    title('Monkey Sp');
    xlabel({'proportion of trials with phasic response';'followed by suppression'});
    ylabel('proportion of LC neurons');
    text(.1,yPlotLim,['n=',num2str(sum(n1))]);
    
    med2 = nanmedian(prop_pnp2); % 0.4381
    axes(axs(2)); cla reset; hold on; ax = gca; disableDefaultInteractivity(ax);
    % plot(xBins,n2/sum(n2),'k','linewidth',1);
    bar(xBins,n2/sum(n2),'k','edgecolor','none','facecolor','k');
    plot([med2 med2],[0 yPlotLim],'r--','linewidth',2);
    
    axis([-0.1 1.1 -0.005 yPlotLim]);
    title('Monkey Ci');
    xlabel({'proportion of trials with phasic response';'followed by suppression'});
    ylabel('proportion of LC neurons');
    text(.1,yPlotLim,['n=',num2str(sum(n2))]);
    
    
end % plotYesNo.





