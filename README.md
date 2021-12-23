# LC_ACC_paper_Joshi_Gold_2021
Code and data for generating figures for this paper.

NOTES for code and data accompanying “Context-Dependent Relationships between Locus Coeruleus Firing Patterns and Coordinated Neural Activity in the Anterior Cingulate Cortex”, S. Joshi and JI Gold.

Contacts: thesidjoshi@gmail.com and jigold@pennmedicine.upenn.edu

1. Use ONLY the “DataFiles” and “FigCode” folders. The “DataCode” scripts are used ONLY for creating data structures from experiment raw data files.

2. Set the data directory appropriately so it points to wherever you save the data structures in the “DataFiles” folder. The “Fig” MATLAB scripts in the “FigCode” will generate all figures in the paper and in some cases, a few supplements that were not in the final manuscript.

3. Run figure scripts with the first section commented as in the uploaded files (reproduced below):

% Don't reanalyze data?

reanalyzeLC = false;

reanalyzeACC = false;

% Reanalyze data?

% reanalyzeLC = true;

% reanalyzeACC = true;
