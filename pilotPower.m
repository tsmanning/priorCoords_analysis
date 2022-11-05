
% Find number of subjects needed for each plot in poster
load('/home/tyler/Documents/MATLAB/cooperLab/1-Experiments_stimuli/VSS2022archive/data/xSub/groupStats.mat')

alpha = 0.05;
tails = 2;

pairs = [1 3;...    % within, slow, t<r
         2 4;...    % within, fast, t<r
         5 7;...    % within, slow, t>r
         6 8;...    % within, fast, t>r
         9 11;...   % between, fast
         2 13;...   % size cntl, t<r 
         6 14];     % size cntl, t>r

numGroups = size(pairs,1);

for ii = 1:numGroups

    thisx = stats.PSEs(pairs(ii,1),:);
    thisy = stats.PSEs(pairs(ii,2),:);

    [nCrit(ii),theseStats(ii)] = doPowerAnalysis(thisx,thisy,alpha,tails);
end