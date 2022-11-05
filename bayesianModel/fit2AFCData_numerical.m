
% Numerically estimate the best fitting Bayesian parameters for the speed
% estimation x distance task. 

function [outMat] = fit2AFCData_numerical(stimV,stimC,stimD,r,nB,xgrid,geoOn,prs0)

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


% ----------------------------
% Extract sizes & initialize 
% ----------------------------

[conts,~,contInds]    = unique(stimC(:));
contInds = reshape(contInds,size(stimC));
numConts = numel(conts);

[dists,~,distInds]    = unique(stimD(:));
distInds = reshape(distInds,size(stimC));
numDists = numel(dists);

% Initialize parameters, if necessary
if nargin < 8
    pNu0      = zeros(nB,1);                  % initial values of component centers
    pGam0     = 2*ones(nB,1);                 % initial values of component widths
    pW0       = normpdf(1:nB,(nB+1)/2,nB/4)';
    pW0       = pW0./sum(pW0);                % initial value of MoG weights

    pSigNse0  = log(ones(numConts,1));        % initial estimate of noise stdev
    
    if geoOn ~= 1
        pSigDist0 = log(ones(numDists,1));   % initial estimate of dist-dependent noise
        prs0      = [pNu0;pGam0;pW0;pSigNse0;pSigDist0];
    else
        prs0      = [pNu0;pGam0;pW0;pSigNse0];
    end
end


% ----------------------------
% Determine whether to fit distance parameter or use geometric noise prop
% model
% ----------------------------

if geoOn == 1
    stimD_SF = distScaleFactor(stimD);
end


% ----------------------------
% Numerically optimize negative log-likelihood
% ----------------------------

zeroCent = 1;

if zeroCent == 1
    nuLim = eps;
else
    nuLim = 10;
end

% bounds for parameters
if geoOn ~= 1
    % [priorMeans priorWidths priorWeights stimUncertainties]
    LB  = [-nuLim*ones(nB,1); -10*zeros(nB,1); zeros(nB,1);   -10*ones(numConts,1)]; % lower bound
    UB  = [nuLim*ones(nB,1);  10*ones(nB,1);  ones(nB,1)+eps; 10*ones(numConts,1)]; % upper bound
    Aeq = [zeros(1,2*nB) ones(1,nB) zeros(1,numConts + numDists)]; % equality constraint (prior weights sum to 1)
    beq = 1;                              % equality constraint
else
    % [priorMeans priorWidths priorWeights stimUncertainties distUncertainties]
    LB  = [-nuLim*ones(nB,1); -10*zeros(nB,1); zeros(nB,1);   -10*ones(numConts,1); -10*ones(numDists,1)]; % lower bound
    UB  = [nuLim*ones(nB,1);  10*ones(nB,1);  ones(nB,1)+eps; 10*ones(numConts,1); 10*ones(numDists,1)]; % upper bound
    Aeq = [zeros(1,2*nB) ones(1,nB) zeros(1,numConts + numDists)]; % equality constraint (prior weights sum to 1)
    beq = 1;                              % equality constraint
end

% loss function
%-------------------------

% Numerical
lossfun = @(prs)(neglogli_BaysObsModel_numerical(prs,stimV,contInds,distInds,r,nB));

% Analytical
% lossfun = @(prs)(neglogli_BaysObsModel_analytical(prs,stimV,contInds,distInds,r,nB));

% optimization options
opts = optimset('display', 'iter','algorithm','interior-point'); 
% opts = optimset('display', 'iter','algorithm','sqp'); 
% opts = bads('defaults');

% Perform optimization 
%-------------------------
% Use fmincon
prshat = fmincon(lossfun,prs0,[],[],Aeq,beq,LB,UB,[],opts);

% % Use BADS
% prshat = bads(lossfun,prs0,LB,UB,[],[],[],opts);

% ----------------------------
% Extract fitted parameters 
% ----------------------------

pNuHat  = prshat(1:nB);
pGamHat = exp(prshat(nB + 1:2*nB));
pWHat   = prshat(2*nB + 1:3*nB);

pSigNseHat  = exp(prshat(3*nB + 1:3*nB + numConts));
pSigDistHat = exp(prshat(3*nB + numConts + 1:3*nB + numConts + numDists));

% log-likelihood at optimum
nll = lossfun(prshat);

% optimal fit of prior
priHat = buildMoGPrior(pGamHat,pNuHat,pWHat,xgrid);

% output best fit pars, nll, estimate of prior
outMat.pNuHat      = pNuHat;
outMat.pGamHat     = pGamHat;
outMat.pWHat       = pWHAt;
outMat.pSigNseHat  = pSigNseHat;
outMat.pSigDistHat = pSigDistHat;
outMat.nll         = nll;
outMat.priHat      = priHat;

end

% ===================================================================
% LOSS FUNCTION: negative log-likelihood (numerical)
% ===================================================================
function [nll] = neglogli_BaysObsModel_numerical(prs,stimV,contInds,distInds,r,nB)

% Computes negative log-likelihood of data (for optimization)

% Extract parameters
pNu  = prs(1:nB);
pGam = exp(prs(nB + 1:2*nB));
pW   = prs(2*nB + 1:3*nB);

sigC = exp(prs(3*nB + 1:3*nB + numConts));
sigD = exp(prs(3*nB + numConts + 1:3*nB + numConts + numDists));

mu1  = stimV(:,1);
mu2  = stimV(:,2);

% Assign parameters to individual trials
sigs = sigC(contInds) + sigD(distInds);
sig1 = sigs(:,1);
sig2 = sigs(:,2);

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
end

% Compute negative log-likelihood of data
nll = -r'*log(p') - (1 - r')*log(1-p');

end

% ===================================================================
% LOSS FUNCTION: negative log-likelihood (analytical)
% ===================================================================
function [nll] = neglogli_BaysObsModel_analytical(prs,stimV,contInds,distInds,r,nB)

% Computes negative log-likelihood of data (for optimization)

% Extract parameters
pNu  = prs(1:nB);
pGam = exp(prs(nB + 1:2*nB));
pW   = prs(2*nB + 1:3*nB);

sigC = exp(prs(3*nB + 1:3*nB + numConts));
sigD = exp(prs(3*nB + numConts + 1:3*nB + numConts + numDists));

mu1  = stimV(:,1);
mu2  = stimV(:,2);

% Assign parameters to individual trials
sigs = sigC(contInds) + sigD(distInds);
sig1 = sigs(:,1);
sig2 = sigs(:,2);

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
    
    p(ii) = calcMoGPFxn_Analytic(sup1,sup2,pNu,pGam,pW,mu1(ii),sig1(ii),mu2(ii),sig2(ii),0);
end

% Compute negative log-likelihood of data
nll = -r'*log(p') - (1 - r')*log(1-p');

end