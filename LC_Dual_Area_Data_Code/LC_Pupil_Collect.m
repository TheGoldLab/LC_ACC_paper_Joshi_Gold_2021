
% LC_Dual_Recording_Collect.m
%
% Sidd - 062116: LC-ACC dual area recording, preliminary analysis - fixation data ONLY!
%
% This script runs through the raw (.nex) data files and generates FIRA structures if needed.
%
% It then generates "clean" data structures as we did for the LC-pupil paper.
%

%% Note - LC and IC data are FIRA'ed separately for now... thus the two sections below...

%% First - get FIRA structure from nex: for LC data

clear; clear all; % If needed....

% monks = {'Sprout'}; % Add more monks as needed...
% monks = {'Cicero'}; % Add more monks as needed...
monks = {'Cicero'}; % Add more monks as needed...
% monks = {'OZ'}; % Add more monks as needed...

% base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data'; % Base directory for brain area.

% ***************************************************************************
% For temp testing of daily data:
% base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Sprout\IC_ACC_Fixation';
% base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Sprout\LC_ACC_Fixation';
% base_dir = 'C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Sprout\LC_ACC_uStim';

base_dir = 'C:\Sidd\PostDoc2\Data\ALL_PUPIL_DATA_ALIGNED\LC';

% ***************************************************************************

% CHANGE THIS TO REBUILD ALL FIRA FILES FROM NEX:

REMAKE_FIRA = false;
% REMAKE_FIRA = true;
global FIRA

% Now go through the data files and create
for mm = 1:length(monks)
    
    % outDir = strcat([base_dir,'\',monks{mm},'\LC_ACC_Fixation\mat']); % Create dir name for output (mat) files.
    % cleanDir = strcat([base_dir,'\',monks{mm},'\LC_ACC_Fixation\clean']); % Create dir name for output (mat) files.
    outDir = strcat([base_dir,'\',monks{mm},'\','mat']); % Create dir name for output (mat) files.
    cleanDir = strcat([base_dir,'\',monks{mm},'\','clean']); % Create dir name for output (mat) files.
    
    if REMAKE_FIRA
        
        % inDir= strcat([base_dir,'\',monks{mm},'\LC_ACC_Fixation\nex']); % Create dir name for input (nex) files.
        inDir= strcat([base_dir,'\',monks{mm},'\','nex']); % Create dir name for input (nex) files.
        cd(inDir);  dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
        fnames = {dirData(~dirIndex).name}'; % Get nex data filenames.
        % pause
        % cd(strcat([base_dir,'\',monks{mm}])); % Git to the mnky directory.
        
        nn= length(fnames);
        for uu = 1:nn % Go through the nex files and bNex them.
            % pause
            tfName = fnames{uu};
            bNex(fullfile(inDir, tfName), 'spmRKLCPupil', fullfile(outDir, strcat([tfName(1:end-4), '_OUT'])), 'all', 'all', 0, 1, [], []);
            
            disp(sprintf('File %d%sof%s%d%s%s', uu,' ',' ',nn,', ',tfName));
            
        end
        
    end
    
    % FIRA's are created... mat files shd be in the mat subdirectory.
    
    cd(outDir);   dirData = dir(pwd); dirIndex = [dirData.isdir]; % Get dir listing.
    fnames = {dirData(~dirIndex).name}'; % Get mat data filenames.
    nn = length(fnames);
    
    for dc = 1:nn
        %         dc=1
        tfName = fnames{dc};
        aa = load(tfName);
        
        if isfield(aa, 'FIRA')
            FIRA = aa.FIRA;
        elseif isfield(aa, 'data')
            FIRA = aa.data;
        end
        
        clear aa;
        
        % Now call "cleanLC_DualArea_FIRA" to do the heavy lifting...
        
        % This will create the clean data structures similar to what we had for
        % the LC-Pupil analysis:
        
        basename = tfName(1:end-8);
        cleanLC_DualArea_FIRA(fullfile(cleanDir, basename));
        disp(sprintf('File %d%sof%s%d%s%s', dc,' ',' ',nn,', ',basename));
        
    end
    
end

% % FIRA and/or "clean" data structures are created and ready for consumption! Bon appetit...

%% Next - get FIRA structure from nex: FOR IC data

% % clear; clear all;
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
%         nn= length(fnames);
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
%         % Now call "cleanLC_DualArea_FIRA" to do the heavy lifting...
%         
%         % This will create the clean data structures similar to what we had for
%         % the LC-Pupil analysis:
% 
%         basename = tfName(1:end-8);
%         cleanLC_DualArea_FIRA(fullfile(cleanDir, basename));
%         disp(sprintf('File %d%sof%s%d%s%s', dc,' ',' ',nn,', ',basename));
%         
%     end
%     
% end
% 
% % FIRA and/or "clean" data structures are created and ready for consumption! Bon appetit...


