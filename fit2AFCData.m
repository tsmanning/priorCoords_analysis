 
% Numerically estimate the best fitting Bayesian parameters for the speed
% estimation x distance task. 

function [outMat] = fit2AFCData(sortedPars,ipd,numPC,xgrid,geoOn,numStartPts,solverOpts)

% Fit Bayesian ideal observer model with fixed Gaussian noise from observer
% 2AFC data, consisting of triplets of stimuli [x1,x2] and observer choices
% 'x2>x1: yes/no'. Numerically calculates p("yes"|x1,x2) by summing
% probability of joint likelihood above decision boundary.
%
% Inputs:
% -------
%          stim [N x 2] - array of presented stimuli
%             r [N x 1] - column vector of observer choices
%                    nB - scalar defining number of Gaussians in mixture
%        xgrid [Nx x 1] - grid of x evenly-spaced values on which prior is defined
%       prs0 [Nb+1 x 1] - initial parameters: [signse0; bwts0] (OPTIONAL)
%
% Output:
% -------
%    pNuHat          - estimated prior on x grid
%    pGamHat         - estimated prior on x grid
%    pWHat           - estimated prior on x grid
%    pSigNse(1/2)Hat - estimates of stdev of encoding noise
%    priHat          - estimated prior on x grid
%    nll             - negative log-likelihood of data, given parameters (OPTIONAL)

%%
addpath(genpath('/home/tyler/Documents/MATLAB/lib/bads'))
addpath(genpath('/home/tyler/Documents/MATLAB/cooperLab/4-Papers/bayesIdealObserverMoG/'));

% ----------------------------
% Extract sizes & initialize 
% ----------------------------

stimV = sortedPars.stimV;
stimC = sortedPars.stimC;
stimD = sortedPars.stimD;
refXY = sortedPars.refXY;
testXY = sortedPars.testXY;
r      = sortedPars.r;

[conts,~,contInds]    = unique(stimC(:));
contInds = reshape(contInds,size(stimC));
numConts = numel(conts);

[dists,~,distInds]    = unique(stimD(:));
distInds = reshape(distInds,size(stimC));
numDists = numel(dists);

% Put response in term of TEST selected y/n
r = ~r;


% ----------------------------
% Determine whether to fit distance parameter or use geometric noise prop
% model
% ----------------------------

if geoOn == 1
    % Ref
    pos1     = [refXY stimD(:,1)];
    % Test
    pos2     = [testXY stimD(:,2)];
    eyePos   = ipd/2*[-1 1];
    stimD_SF = [getVarY(pos1,eyePos,1) getVarY(pos2,eyePos,1)];
else
    stimD_SF = [];
end


% ----------------------------
% Numerically optimize negative log-likelihood
% ----------------------------

zeroMean    = 1;

[prs0,constrStr] = designOptimization(geoOn,zeroMean,numPC,numConts,numDists,numStartPts);

LB  = constrStr.LB;
UB  = constrStr.UB;
Aeq = constrStr.Aeq;
beq = constrStr.beq;
A   = constrStr.A;
b   = constrStr.b;

% No consts
A = [];
b = [];
Aeq = [];
beq = [];

% No prior vs likelihood consts
% A = A(1:numConts-1,:);
% b = b(1:numConts-1);

% No contrast-likelihood inverse consts
% A = A(numConts-1:end,:);
% b = b(numConts-1:end);


% ------------------------
% Define loss function
%-------------------------

switch solverOpts.fxnForm
    case 'numerical'
        % Numerical
        lossfun = @(prs)(neglogli_BaysObsModel_numerical(prs,stimV,contInds,distInds,r,numPC,stimD_SF,geoOn,zeroMean));
    case 'analytical'
        % Analytical
        lossfun = @(prs)(neglogli_BaysObsModel_analytical(prs,stimV,contInds,distInds,r,numPC,stimD_SF,geoOn,zeroMean));
end

% optimization options
switch solverOpts.solver
    case 'fmincon'
        switch solverOpts.fxnForm
            case 'numerical'
                opts = optimset('display','iter','algorithm','interior-point','UseParallel',1);
            case 'analytical'
                opts = optimset('display','off','algorithm','interior-point');
        end
    case 'bads'
        opts = bads('defaults');
        opts.Display = 'off';
end


% ------------------------
% Perform optimization 
%-------------------------

for ii = 1:numStartPts
    disp(['Optimization run ',num2str(ii),'/',num2str(numStartPts),'...']);

    switch solverOpts.solver
        case 'fmincon'
            % Use fmincon
            [prshatSet(:,ii),nllSet(ii)] = fmincon(lossfun,prs0(:,ii),A,b,Aeq,beq,LB,UB,[],opts);
        case 'bads'
            % Use BADS
            [prshatSet(:,ii),nllSet(ii)] = bads(lossfun,prs0(:,ii)',LB',UB',LB',UB',[],opts);
    end
end

[nll,minInd] = min(nllSet);

prshat = prshatSet(:,minInd);

[~,p] = testerNLL(prshat,stimV,contInds,distInds,r,numPC,stimD_SF,geoOn,zeroMean);

% ----------------------------
% Extract fitted parameters 
% ----------------------------

[pNuHat,pGamHat,pWHat,pSigNseHat,pSigDistHat] = extractParameters(prshat,geoOn,zeroMean,numPC,numConts,numDists);

% optimal fit of prior
priHat = buildMoGPrior(pGamHat,pNuHat,pWHat,xgrid);

% output best fit pars, nll, estimate of prior
outMat.prs0        = prs0;
outMat.bestInd     = minInd;
outMat.pNuHat      = pNuHat;
outMat.pGamHat     = pGamHat;
outMat.pWHat       = pWHat;
outMat.pSigNseHat  = pSigNseHat;
outMat.pSigDistHat = pSigDistHat;
outMat.stimD_SF    = unique(stimD_SF);
outMat.stimD_SFall = stimD_SF;
outMat.nll         = nll;
outMat.priHat      = priHat;
outMat.xgrid       = xgrid;
outMat.contVals    = conts;
outMat.distVals    = dists;
outMat.contInds    = contInds;
outMat.distInds    = distInds;
outMat.constrStr   = constrStr;
outMat.p           = p';

end


% ===================================================================
% LOSS FUNCTION: negative log-likelihood (numerical)
% ===================================================================
function [nll] = neglogli_BaysObsModel_numerical(prs,stimV,contInds,distInds,r,nPC,stimD_SF,geoOn,zeroMean)

% Computes negative log-likelihood of data (for optimization)

numConts = numel(unique(contInds));
numDists = numel(unique(distInds));

% Extract parameters
[pNu,pGam,pW,sigC,sigD] = extractParameters(prs,geoOn,zeroMean,nPC,numConts,numDists);

mu1  = stimV(:,1);
mu2  = stimV(:,2);

%%%%%%%%%%%%% really should change this so we're giving distance as an
%%%%%%%%%%%%% argument
% Convert scale factor into distance and scale factor
if ~geoOn
    % estDis is an Nx2 matrix
    estDis = sqrt(2)*sigD(distInds).^2;
else
    estDis = sqrt(2*stimD_SF);
end
scFac  = sqrt(stimD_SF);

% Transform retinal speeds into world speeds
% (technically this overestimates by ~3.2%)
mu1 = mu1.*estDis(:,1);
mu2 = mu2.*estDis(:,2);

% Transform retinal uncertainties into world uncertainties
if ~geoOn
    sigs = sigC(contInds) .* sigD(distInds);
else
    sigs = sigC(contInds) .* sqrt(scFac);
end

sig1 = sigs(:,1);
sig2 = sigs(:,2);

% Log transform data
pNu  = getLogXform(pNu,0.3);
pGam = getLogXform(pGam,0.3);
mu1  = getLogXform(mu1,0.3);
mu2  = getLogXform(mu2,0.3);
sig1 = getLogXform(sig1,0.3);
sig2 = getLogXform(sig2,0.3);

% Generate supports
dx      = 100;

nTrials  = numel(mu1);

% Calculate psychometric function values
for ii = 1:nTrials
    
    suppLB1 = mu1(ii) - 4*sig1(ii);
    suppUB1 = mu1(ii) + 4*sig1(ii);
    suppLB2 = mu2(ii) - 4*sig2(ii);
    suppUB2 = mu2(ii) + 4*sig2(ii);
    
    sup1 = linspace(suppLB1,suppUB1,dx);
    sup2 = linspace(suppLB2,suppUB2,dx);
    
    p(ii) = calcMoGPFxn_Numeric(sup1,sup2,pNu,pGam,pW,mu1(ii),sig1(ii),mu2(ii),sig2(ii),0);

    if p(ii) == 0
        p(ii) = eps;
    end
    if p(ii) == 1
        p(ii) = 1-eps;
    end
end

% Compute negative log-likelihood of data
nll = -r'*log(p') - (1 - r')*log(1-p');

end

% ===================================================================
% LOSS FUNCTION: negative log-likelihood (analytical)
% ===================================================================
function [nll] = neglogli_BaysObsModel_analytical(prs,stimV,contInds,distInds,r,nPC,stimD_SF,geoOn,zeroMean)

% Computes negative log-likelihood of data (for optimization)

numConts = numel(unique(contInds));
numDists = numel(unique(distInds));

% Extract parameters
[pNu,pGam,pW,sigC,sigD] = extractParameters(prs,geoOn,zeroMean,nPC,numConts,numDists);

mu1  = stimV(:,1);
mu2  = stimV(:,2);

%%%%%%%%%%%%% really should change this so we're giving distance as an
%%%%%%%%%%%%% argument
% Convert scale factor into distance and scale factor
if ~geoOn
    % estDis is an Nx2 matrix
    estDis = sqrt(2)*sigD(distInds).^2;
else
    estDis = sqrt(2*stimD_SF);
end
scFac  = sqrt(stimD_SF);

% Transform retinal speeds into world speeds
% (technically this overestimates by ~3.2%)
mu1 = mu1.*estDis(:,1);
mu2 = mu2.*estDis(:,2);

% Transform retinal uncertainties into world uncertainties
if ~geoOn
    sigs = sigC(contInds) .* sigD(distInds);
else
    sigs = sigC(contInds) .* sqrt(scFac);
end

sig1 = sigs(:,1);
sig2 = sigs(:,2);

% Log transform data
pNu  = getLogXform(pNu,0.3);
pGam = getLogXform(pGam,0.3);
mu1  = getLogXform(mu1,0.3);
mu2  = getLogXform(mu2,0.3);
sig1 = getLogXform(sig1,0.3);
sig2 = getLogXform(sig2,0.3);

nTrials  = numel(mu1);

% Calculate psychometric function values
for ii = 1:nTrials
    
    p(ii) = calcMoGPFxn_Analytic(pNu,pGam,pW,mu1(ii),sig1(ii),mu2(ii),sig2(ii));

    if p(ii) == 0
        p(ii) = eps;
    end
    if p(ii) == 1
        p(ii) = 1-eps;
    end
    if isnan(p(ii))
        p(ii) = eps;
    end
end

% Compute negative log-likelihood of data
nll = -r'*log(p') - (1 - r')*log(1-p');

end