function [nll,p] = testerNLL(prs,stimV,contInds,distInds,r,nPC,stimD_SF,geoOn,zeroMean)

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