function [prs0,constrStr] = designOptimization(geoOn,zeroMean,numPC,numConts,numDists,numStartPts)

% Make constraints/starting points for solver based on desired options

%% Initialize parameters

nuLim = 2;

if zeroMean
    pNu0 = [];
else
    pNu0 = zeros(numPC,1);               % initial values of component centers
end

pGam0    = log(2*ones(numPC,1));         % initial values of component widths
pW0      = 0.5*ones(numPC,1);
pW0      = pW0./sum(pW0);                % initial value of MoG weights

pSigNse0 = log(0.5*ones(numConts,1));    % initial estimate of noise stdev

if geoOn
    pSigDist0 = [];
else
    pSigDist0 = log(ones(numDists,1));   % initial estimate of dist-dependent noise
end

prs0 = [pNu0;pGam0;pW0;pSigNse0;pSigDist0];

if numStartPts > 1

    nuScF  = nuLim/3;
    gamScF = (-10+log(10))/4;
    wScF   = 0.3;
    conScF = (-10+log(10))/4;
    disScF = (-10+log(10))/4;

    nuOff  = 0;
    gamOff = 0;
    wOff   = 0;
    conOff = 0;
    disOff = 0;

    randOffs = [nuScF*rand(numel(pNu0),numStartPts) + nuOff;
                gamScF*rand(numPC,numStartPts) + gamOff;
                wScF*rand(numPC,numStartPts) + wOff;
                conScF*rand(numConts,numStartPts) + conOff;
                disScF*rand(numel(pSigDist0),numStartPts) + disOff];

    prs0 = prs0 + randOffs;

    prs0(numel(pNu0)+numel(pGam0)+1:numel(pNu0)+numel(pGam0)+numel(pW0),:) = ...
        prs0(numel(pNu0)+numel(pGam0)+1:numel(pNu0)+numel(pGam0)+numel(pW0),:) ./ ...
        repmat(sum(prs0(numel(pNu0)+numel(pGam0)+1:numel(pNu0)+numel(pGam0)+numel(pW0),:)),[numel(pGam0) 1]);

end


%% Setup constraints
% Constraint matrix to make likelihood widths decrease in size as a function of increasing contrast
% and components of prior larger than likelihoods
contConst = zeros(numConts-1,numConts) + [-eye(numConts-1) zeros(numConts-1,1)] + [zeros(numConts-1,1) eye(numConts-1)];
contConst = [zeros(numConts-1,numPC*3) contConst];

likeConst = [zeros(numPC*numConts,numPC) repmat(-eye(numPC),[numConts 1]) zeros(numPC*numConts,numPC) repelem(eye(numConts),numPC,1)];

combConst = [contConst;likeConst];

if geoOn == 1
    % [priorMeans priorWidths priorWeights stimUncertainties]
    LB  = [ zeros(numPC,1)-eps;     -10*ones(numPC,1);    zeros(numPC,1);     -10*ones(numConts,1)]; % lower bound
    UB  = [nuLim*ones(numPC,1); log(10)*ones(numPC,1); ones(numPC,1)+eps; log(10)*ones(numConts,1)]; % upper bound

    Aeq = [zeros(1,2*numPC)                             ones(1,numPC)         zeros(1,numConts)]; % prior weights sum to 1
    beq = 1;                                           

    A   = combConst; % sum of prior widths must be greater than sum of likelihoods
    b   = zeros(size(A,1),1);
else
    % [priorMeans priorWidths priorWeights stimUncertainties distUncertainties]
    % constraint on distuncertainties ensures reasonable distance values
    % are extracted
    LB  = [ zeros(numPC,1)-eps;     -10*ones(numPC,1);    zeros(numPC,1);     -10*ones(numConts,1); -1.3246*ones(numDists,1)]; % lower bound
    UB  = [nuLim*ones(numPC,1); log(10)*ones(numPC,1); ones(numPC,1)+eps; log(10)*ones(numConts,1);  0.0294*ones(numDists,1)]; % upper bound

    Aeq = [zeros(1,2*numPC)                             ones(1,numPC)                        zeros(1,numConts + numDists)]; % prior weights sum to 1
    beq = 1;

    A   = [combConst zeros(size(combConst,1),numDists)]; % sum of prior widths must be greater than sum of likelihoods
    b   = zeros(size(A,1),1);
end

% If constraining prior components to be zero mean, just lop off the front
% of the constraint arrays
if zeroMean
    LB = LB(numPC+1:end,:);
    UB = UB(numPC+1:end,:);

    Aeq = Aeq(numPC+1:end);

    A = A(:,numPC+1:end);
end


%% Package up constraint structure
constrStr.LB  = LB;
constrStr.UB  = UB;
constrStr.Aeq = Aeq;
constrStr.beq = beq;
constrStr.A   = A;
constrStr.b   = b;


end