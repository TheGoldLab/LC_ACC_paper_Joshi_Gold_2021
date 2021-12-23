
% LC_Dual_Recording_CP_Collect.m
% 063016: LC-ACC dual area recording, preliminary analysis - Change Point Reward Task
%
% Modified from:
% Sidd_Dual_Recording_CP_Collect.m
%
% This script runs through the raw data files (.nex) and generates FIRA structures if needed.
%
% It then generates "clean" data structures as we did for the LC-pupil paper.


%% LC Data Conversion:

% First - get FIRA structure from nex:

clear; clear all;

monks = {'sprout'}; % Add monks as needed....
% sites = ???? % Add sites as needed... for now, LC-ACC only.
% base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.
base_dir = 'E:\NexHannah'

% CHANGE THIS TO REBUILD ALL FIRA FILES FROM NEX:
REMAKE_FIRA = true;
% REMAKE_FIRA = false;

global FIRA

%%

% Now go through the raw data files and create FIRA and clean structures:
for mm = 1:length(monks)
    
    % outDir = strcat([base_dir,'\',monks{mm},'\LC_ACC_ChangePoint\mat']); % Create dir name for output (mat) files.
    outDir = strcat([base_dir,'\',monks{mm},'\LC_ACC_ChangePoint\mat']); % Create dir name for output (mat) FIRA files.
    cleanDir = strcat([base_dir,'\',monks{mm},'\LC_ACC_ChangePoint\clean']); % Create dir name for output (mat) clean files.
    
    if REMAKE_FIRA
        
        inDir= strcat([base_dir,'\',monks{mm},'\LC_ACC_ChangePoint\nex']); % Create dir name for input (nex) files.
        cd(inDir);  dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
        fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
        
        cd(strcat([base_dir,'\',monks{mm}])); % Git to the mnky directory.
        
        nn= length(fnames);
        for uu = 1:nn % Go through the nex files and bNex them (ie convert to FIRA).

            tfName = fnames{uu};
            bNex(fullfile(inDir, tfName), 'spmSiddCP', fullfile(outDir, strcat([tfName(1:end-4), '_OUT'])), 'all', 'all', 0, 1, [], []);
            %             pause
            disp(sprintf('File %d%sof%s%d%s%s', uu,' ',' ',nn,', ',tfName));
            
        end
    end
    
    % FIRA's are created... mat files shd be in the mat subdirectory.
    
    % clear; clear all;
    cd(outDir);   dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
    fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
    nn = length(fnames);
    
    for dc = 1:nn
        
        tfName = fnames{dc};
        aa = load(tfName);
        
        if isfield(aa, 'FIRA')
            FIRA = aa.FIRA;
        elseif isfield(aa, 'data')
            FIRA = aa.data;
        end
        
        clear aa;
        
        % Now call cleanLC_DualArea_FIRA to do the heavy lifting...
        % This will create the clean data structures similar to what we had for
        % the LC-Pupil analysis:
        basename = tfName(1:end-8);
        %         pause
        cleanLC_DualArea_CP_FIRA(fullfile(cleanDir, basename));
        
        disp(sprintf('File %d%sof%s%d%s%s', dc,' ',' ',nn,', ',basename));
        
    end
    
end

% *******************

% %% IC Data Conversion:
%
% % First - get FIRA structure from nex:
%
% clear; clear all;
%
% monks = {'sprout'};
% base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.
%
% % CHANGE THIS TO REBUILD ALL FIRA FILES FROM NEX:
% REMAKE_FIRA = true;
% % REMAKE_FIRA = false;
% global FIRA
%
% % Now go through the data files and create
% for mm = 1:length(monks)
%
%     outDir = strcat([base_dir,'\',monks{mm},'\IC_ACC_Fixation\mat']); % Create dir name for output (mat) files.
%     cleanDir = strcat([base_dir,'\',monks{mm},'\IC_ACC_Fixation\clean']); % Create dir name for output (mat) files.
%
%     if REMAKE_FIRA
%
%         inDir= strcat([base_dir,'\',monks{mm},'\IC_ACC_Fixation\nex']); % Create dir name for input (nex) files.
%         cd(inDir);  dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
%         fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
%
%         cd(strcat([base_dir,'\',monks{mm}])); % Git to the mnky directory.
%
%         nn = length(fnames);
%         for uu = 1:nn % Go through the nex files and bNex them.
%
%             tfName = fnames{uu};
%             bNex(fullfile(inDir, tfName), 'spmRKLCPupil', fullfile(outDir, strcat([tfName(1:end-4), '_OUT'])), 'all', 'all', 0, 1, [], []);
%
%             disp(sprintf('File %d%sof%s%d%s%s', uu,' ',' ',nn,', ',tfName));
%
%         end
%
%     end
%
%     % FIRA's are created... mat files shd be in the mat subdirectory.
%
%     cd(outDir);   dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
%     fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
%     nn = length(fnames);
%
%     for dc = 1:nn
%
%         tfName = fnames{dc};
%         aa = load(tfName);
%
%         if isfield(aa, 'FIRA')
%             FIRA = aa.FIRA;
%         elseif isfield(aa, 'data')
%             FIRA = aa.data;
%         end
%
%         clear aa;
%
%         % Now call cleanLC_DualArea_FIRA to do the heavy lifting...
%         % This will create the clean data structures similar to what we had for
%         % the LC-Pupil analysis:
%         basename = tfName(1:end-8);
%         cleanLC_DualArea_FIRA(fullfile(cleanDir, basename));
%
%     end
%
%     disp(sprintf('File %d%sof%s%d%s%s', dc,' ',' ',nn,', ',basename));
%
% end
%
%
%

