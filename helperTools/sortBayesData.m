function [sortedPars,dataStruct] = sortBayesData(saveDir,subjID,incTrTypes)

%% Organizes data and necessary parameters for Bayesian ideal observer model fitting

% Grab response data
load([saveDir,subjID,'.mat']);

a1 = arrayfun(@(x) x.respMat,stimStruct,'UniformOutput',false); % Make array of respmatss
a2 = reshape(a1,[numel(a1),1]);                                 % Linearize array
a3 = cellfun(@(x) ~isempty(x),a2);                              % Cull empty inds
a2 = a2(a3);

% Get trial types & dists
a4 = stimStruct(:);
mt = arrayfun(@(x) isempty(x.testDist),a4);
a4 = a4(~mt);
withinBlocks  = arrayfun(@(x) strcmp(x.blockType,'within'),a4);
betweenBlocks = arrayfun(@(x) strcmp(x.blockType,'between'),a4);
trueDist      = arrayfun(@(x) (x.testDist == x.trueTestDist) & (x.refDist == x.trueRefDist),a4);

trueTestD     = arrayfun(@(x) x.trueTestDist,a4);
trueRefD      = arrayfun(@(x) x.trueRefDist,a4);

% Decide which trials we want to feed to model
switch incTrTypes
    case 'Within'
        incInds = withinBlocks;
    case 'WithinBetweenPers'
        incInds = withinBlocks | (betweenBlocks & trueDist);
    case 'WithinBetweenAll'
        incInds = withinBlocks | betweenBlocks;
    case 'BetweenAll'
        incInds = betweenBlocks;
end

a2 = a2(incInds);
dataStruct = a4(incInds);

% ... and sort into stimV/stimC/stimD/r [ref, test]
fullRespMat = cell2mat(a2);                                     % Concatenate arrays

stimV = fullRespMat(:,[5 3]);
stimC = fullRespMat(:,[4 2]);
stimD = fullRespMat(:,[19 20]);
r     = fullRespMat(:,21);

% Define x/y positions of stimuli
trueTestD = trueTestD(incInds);
trueRefD  = trueRefD(incInds);
sz = cellfun(@(x) size(x,1),a2);
trd = repelem(trueRefD,sz);
ttd = repelem(trueTestD,sz);

farxy  = [tand(15) tand(6)];

%         refXY  = farxy.*trd;
%         testXY = farxy.*ttd;

% Use a virtual position based on viewing angle and simulated
% distance
refXY  = farxy.*stimD(:,1);
testXY = farxy.*stimD(:,2);

% Reduce respMat down to pars of interest
% (test cont, test vel, ref vel, ref dist, test dist, true ref dist, true test dist)
redRespMat = [fullRespMat(:,[2 3 5 19 20]) trd ttd];
% Find unique stimulus pairs (including test velocity)
[b,bi]  = unique(redRespMat,'rows');
% Find unique trial conditions (i.e. without test velocity)
[c,~,cii]  = unique(b(:,[1 3 4 5 6 7]),'rows');

% Sort respMat
[~,sortInds1] = sortrows(c,'ascend');

% Sort dataStruct
dsOrder = arrayfun(@(x) [x.testCont x.refVel x.refDist x.testDist x.trueRefDist x.trueTestDist],dataStruct,'uniformOutput',0);
dsOrder = cell2mat(dsOrder);
[~,sortInds2] = sortrows(dsOrder,'ascend');

sortedPars.stimV  = stimV;
sortedPars.stimC  = stimC;
sortedPars.stimD  = stimD;
sortedPars.r      = r;
sortedPars.testXY = testXY;
sortedPars.refXY  = refXY;
sortedPars.inds   = {bi,cii,sortInds1,sortInds2};


end