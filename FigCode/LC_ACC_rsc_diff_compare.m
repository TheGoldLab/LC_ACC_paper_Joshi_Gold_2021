% Sidd.
%
% Nov. 2020
% Compare LC conditioned ACC rsc differences (and ACC conditioned LC rsc differences) per tercile.

clear; clear all;

cd C:\Sidd\PostDoc2\Data\LC_Dual_Area_Data\Results\Results_2021\LC_ACC;
ACC_rsc = load('Fig4_LC_ACC_Compare'); ACC_dat = ACC_rsc.all_bdat; % binDVals = ACC_rsc.binDiffVals;
LC_rsc = load('Fig4_ACC_LC_Compare'); LC_dat = LC_rsc.all_bdat_LC; % binDValsL = LC_rsc.binDiffValsL;

nb = length(ACC_dat);
pvals = nans(nb,3); % num bins x 3 terciles.
hvals = nans(nb,3); % num bins x 3 terciles.


% Monkeys combined.
for jj = 1:length(ACC_dat)

    % kstest2.
    [h p] = kstest2(ACC_dat{jj}{1},LC_dat{jj}{1}); pvals(jj,1) = p; hvals(jj,1) = h;
    [h p] = kstest2(ACC_dat{jj}{2},LC_dat{jj}{2}); pvals(jj,2) = p; hvals(jj,2) = h;
    [h p] = kstest2(ACC_dat{jj}{3},LC_dat{jj}{3}); pvals(jj,3) = p; hvals(jj,3) = h;

    % ranksum.
    %         pvals(jj,1) = ranksum(ACC_dat{jj}{1},LC_dat{jj}{1});
    %         pvals(jj,2) = ranksum(ACC_dat{jj}{2},LC_dat{jj}{2});
    %         pvals(jj,3) = ranksum(ACC_dat{jj}{3},LC_dat{jj}{3});
    
end

% % Per-monkey?
% mpVals = [];
% for mm = 1:2
%     pvals = nans(nb,3); % num bins x 3 terciles.
%     hvals = nans(nb,3); % num bins x 3 terciles.
%     for jj = 1:length(ACC_dat)
%         [h p] = kstest2(binDVals{jj}{1}{mm},binDValsL{jj}{1}{mm}); pvals(jj,1) = p; hvals(jj,1) = h;
%         [h p] = kstest2(binDVals{jj}{2}{mm},binDValsL{jj}{2}{mm}); pvals(jj,2) = p; hvals(jj,2) = h;
%         [h p] = kstest2(binDVals{jj}{3}{mm},binDValsL{jj}{3}{mm}); pvals(jj,3) = p; hvals(jj,3) = h;
%     end
%     mpVals{mm} = {pvals hvals};
% end

