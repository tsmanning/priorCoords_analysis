function [pNuHat,pGamHat,pWHat,pSigNseHat,pSigDistHat] = extractParameters(prshat,geoOn,zeroMean,nPC,numConts,numDists)

% Extract optimized parameters from solver output vecotr depending upon
% options

% NOTE: ASSUMES GAUSSIAN COMPONENT SIGMAS ARE PRE-LOG-XFORMED
if zeroMean
    pNuHat      = zeros(nPC,1);
    pGamHat     = exp(prshat(1:nPC));
    pWHat       = prshat(nPC + 1:2*nPC);
    pSigNseHat  = exp(prshat(2*nPC + 1:2*nPC + numConts));

    if geoOn
        pSigDistHat = [];
    else
        pSigDistHat = exp(prshat(2*nPC + numConts + 1:2*nPC + numConts + numDists));
    end
else
    pNuHat  = prshat(1:nPC);
    pGamHat = exp(prshat(nPC + 1:2*nPC));
    pWHat   = prshat(2*nPC + 1:3*nPC);
    pSigNseHat  = exp(prshat(3*nPC + 1:3*nPC + numConts));

    if geoOn
        pSigDistHat = [];
    else
        pSigDistHat = exp(prshat(3*nPC + numConts + 1:3*nPC + numConts + numDists));
    end
end


end