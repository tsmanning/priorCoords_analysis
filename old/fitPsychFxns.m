function [stimStruct,f] = fitPsychFxns(stimV,stimC,stimD,r,fitMatSG,fitMatMoG,expType)

% Takes data from psychophysics exp, fits with cumulative Gaussian and also
% plots Bayesian fit for comparison

%% Parse Data

% Number of psychometric functions
[uniqTrialTypes,uniqInd] = unique([(stimV(:,1)) stimC stimD],'rows');
numPsychFxns             = size(uniqTrialTypes,1);

fitOn = 1;
plotOn = 1;


%% Loop over psychometric functions (unique trial types)

if fitOn
for ii = 1:numPsychFxns

    % Find stimulus values for this psych function
    thisRefVel   = uniqTrialTypes(ii,1);
    thisRefCon   = uniqTrialTypes(ii,2);
    thisRefDis   = uniqTrialTypes(ii,4);
    thisTestCont = uniqTrialTypes(ii,3);
    thisTestDist = uniqTrialTypes(ii,5);

    theseInds = (stimV(:,1) == thisRefVel) & ...
                (stimC(:,1) == thisRefCon) & ...
                (stimD(:,1) == thisRefDis) & ...
                (stimC(:,2) == thisTestCont) & ...
                (stimD(:,2) == thisTestDist);

    theseTestVel = stimV(theseInds,2);
    theseResps   = ~r(theseInds);

    % Find parameter indices for this psych function
    refTestCInd = fitMatSG.contInds(uniqInd(ii),:);
    refTestDInd = fitMatSG.distInds(uniqInd(ii),:);

    % Calculate proportion of "test faster"
    [uniqueTestVel,~,uniInds] = unique(theseTestVel);

    respP  = nan(numel(uniqueTestVel),1);
    numInp = nan(numel(uniqueTestVel),1);

    for ll = 1:numel(uniqueTestVel)
        
        % Find all responses to this test vel
        inputs     = theseResps(uniInds == ll);
        respP(ll)  = sum(inputs)/numel(inputs);
        numInp(ll) = numel(inputs);

    end

    % Fit cum. Gauss with psignifit & get CIs on fit parameters
    respData{ii} = [getLogXform(uniqueTestVel,0.3) respP.*numInp numInp];

    options                = struct;
    options.confP          = 0.95;
    options.estimateType   = 'MAP';
    options.expType        = 'YesNo';
    options.sigmoidName    = 'norm';
    options.stepN          = [80,80,1,1,1];
    options.threshPC       = 0.5;
    options.fixedPars      = [nan;nan;0;0;0];
    options.nblocks        = size(respData{ii},1);
    options.maxBorderValue = exp(-20);

    result  = psignifit(respData{ii},options);
    pFxn    = result.psiHandle(getLogXform(uniqueTestVel,0.3));

    % Calculate NLL for cumulative Gaussian fit & coin flip
    nllCG = nan(numel(uniqueTestVel),1);
    nllCF = nan(numel(uniqueTestVel),1);

    for ll = 1:numel(uniqueTestVel)
        inputs = theseResps(uniInds == ll);
        pCGfit = pFxn(ll);
        pCF    = 0.5;

        nllCG(ll) = -sum(inputs*log(pCGfit) + (1 - inputs)*log(1-pCGfit));
        nllCF(ll) = -sum(inputs*log(pCF) + (1 - inputs)*log(1-pCF));
    end

    nllCG = sum(nllCG);
    nllCF = sum(nllCF);

    % Plot estimates of PSE and slope with CIs
    PSE   = result.Fit(1);
    PSECI = result.conf_Intervals(1,:);

    % width = SD for cumulative Gaussian fit
    width   = result.Fit(2);
    widthCI = result.conf_Intervals(2,:);

    % Get slope and slope CIs
    alpha               = 0.05;
    C                   = norminv(1-alpha,0,1) - norminv(alpha,0,1);
    stimLevel           = PSE;
    normalizedStimLevel = (stimLevel-PSE)/width.*C;
    slopeNormalized     = normpdf(normalizedStimLevel,0,1);

    slope               = slopeNormalized *C./width;
    slopeCI             = slopeNormalized *C./widthCI;

    % Calculate psychometric function using Bayesian ideal observer
    % model
    priorMusSG  = fitMatSG.pNuHat;
    priorSigsSG = fitMatSG.pGamHat;
    priorWsSG   = fitMatSG.pWHat;

    priorMus  = fitMatMoG.pNuHat;
    priorSigs = fitMatMoG.pGamHat;
    priorWs   = fitMatMoG.pWHat;

    priorMus  = getLogXform(priorMus,0.3);
    priorSigs = getLogXform(priorSigs,0.3);
    priorMusSG  = getLogXform(priorMusSG,0.3);
    priorSigsSG = getLogXform(priorSigsSG,0.3);

    if ~isempty(fitMatMoG.pSigDistHat)
        % Ind 1: Ref, ind 2: test
        likeSigs   = [fitMatMoG.pSigNseHat(refTestCInd(1)) * fitMatMoG.pSigDistHat(refTestDInd(1));
                      fitMatMoG.pSigNseHat(refTestCInd(2)) * fitMatMoG.pSigDistHat(refTestDInd(2))];
        likeSigsSG = [fitMatSG.pSigNseHat(refTestCInd(1)) * fitMatMoG.pSigDistHat(refTestDInd(1));
                      fitMatSG.pSigNseHat(refTestCInd(2)) * fitMatMoG.pSigDistHat(refTestDInd(2))];

        estDist    = sqrt(2)*[fitMatMoG.pSigDistHat(refTestDInd(1)) fitMatMoG.pSigDistHat(refTestDInd(2))];
    else
        likeSigs   = [fitMatMoG.pSigNseHat(refTestCInd(1)) * sqrt(fitMatMoG.stimD_SF(refTestDInd(1)));
                      fitMatMoG.pSigNseHat(refTestCInd(2)) * sqrt(fitMatMoG.stimD_SF(refTestDInd(2)))];
        likeSigsSG = [fitMatSG.pSigNseHat(refTestCInd(1)) * sqrt(fitMatMoG.stimD_SF(refTestDInd(1)));
                      fitMatSG.pSigNseHat(refTestCInd(2)) * sqrt(fitMatMoG.stimD_SF(refTestDInd(2)))];
        estDist    = sqrt(2*[fitMatMoG.stimD_SF(refTestDInd(1)) fitMatMoG.stimD_SF(refTestDInd(2))]);
    end

    pFxnBayes = nan(numel(uniqueTestVel),1);
    pFxnBayesSG = nan(numel(uniqueTestVel),1);

    for ll = 1:numel(uniqueTestVel)

        % Ind 1: Ref, ind 2: test
        mu1 = thisRefVel*estDist(1);
        mu2 = uniqueTestVel(ll)*estDist(2);

        dV = 50;

        sig1 = likeSigs(1);
        sig2 = likeSigs(2);
        sig1SG = likeSigsSG(1);
        sig2SG = likeSigsSG(2);

        mu1  = getLogXform(mu1,0.3);
        mu2  = getLogXform(mu2,0.3);
        sig1 = getLogXform(sig1,0.3);
        sig2 = getLogXform(sig2,0.3);
        sig1SG = getLogXform(sig1SG,0.3);
        sig2SG = getLogXform(sig2SG,0.3);

        support1 = linspace(mu1-3*sig1,mu1+3*sig1,dV);
        support2 = linspace(mu2-3*sig2,mu2+3*sig2,dV);
        support1SG = linspace(mu1-3*sig1SG,mu1+3*sig1SG,dV);
        support2SG = linspace(mu2-3*sig2SG,mu2+3*sig2SG,dV);

        pFxnBayesSG(ll) = calcMoGPFxn_Numeric(support1SG,support2SG,priorMusSG,priorSigsSG,priorWsSG,mu1,sig1SG,mu2,sig2SG,0);
        pFxnBayes(ll)   = calcMoGPFxn_Numeric(support1,support2,priorMus,priorSigs,priorWs,mu1,sig1,mu2,sig2,0);
%         pFxnBayesSG(ll) = calcMoGPFxn_Analytic(priorMusSG,priorSigsSG,priorWsSG,mu1,sig1SG,mu2,sig2SG);
%         pFxnBayes(ll) = calcMoGPFxn_Analytic(priorMus,priorSigs,priorWs,mu1,sig1,mu2,sig2);
    end

    % Output all fits and data as a structure
    stimStruct(ii).PSE          = PSE;
    stimStruct(ii).PSECI        = PSECI;
    stimStruct(ii).slope        = slope;
    stimStruct(ii).slopeCI      = slopeCI;
    stimStruct(ii).respData     = respData{ii};
    stimStruct(ii).pFxn         = pFxn;
    stimStruct(ii).pFxnBayesSG  = pFxnBayesSG;
    stimStruct(ii).pFxnBayes    = pFxnBayes;
    stimStruct(ii).nllCG        = nllCG;
    stimStruct(ii).nllCF        = nllCF;
    stimStruct(ii).trialType    = uniqTrialTypes(ii,:);

end
end


%% Plot

% Plotting pars
sTickLin = [0.1 0.5 1 2 4 8 12 16 20];
speedTick = getLogXform(sTickLin,0.3);
for ii = 1:numel(sTickLin)
    sTickLab{ii} = sTickLin(ii);
end
speedLims = getLogXform([0.1 16],0.3);

[uConts,~,contInds] = unique(uniqTrialTypes(:,2));
[uDists,~,distInds] = unique(uniqTrialTypes(:,4));
numRefC = numel(uConts);
numRefD = numel(uDists);

switch expType
    case 'within'
        % Within screen: # figs = #ref speeds
        [~,~,figInds]       = unique(uniqTrialTypes(:,1));
        numRows = numRefC;
        numColumns = numRefD;
        spInds = reshape(1:numRows*numColumns,[numRows numColumns]);

    case 'between'
        % Between screens: # figs = #far distances
        [~,~,figInds]       = unique(uniqTrialTypes(:,1));
        numColumns = numRefD;
        numRows    = (size(uniqTrialTypes,1)+1)/numColumns;
        spInds     = reshape(1:numRows*numColumns,[numColumns numRows])';

    case 'sizeControl'
        % Within screen: # figs = #ref speeds
        [~,~,figInds]       = unique(uniqTrialTypes(:,1));
        numRows    = numRefC;
        numColumns = numRefD;
        spInds = reshape(1:numRows*numColumns,[numRows numColumns]);
end

% Linearize indexing
spInds = spInds(:);

% Skip second index for retinally identical case
matchedInd = find(uniqTrialTypes(:,4) == uniqTrialTypes(:,5))+1;
keepInds   = ones(size(spInds),'logical');
keepInds(matchedInd)   = false;
spInds     = spInds(keepInds);

fixorder   = [1 2 3 4 6 5 7];
spInds     = spInds(fixorder);

position   = [100 100 450*numColumns+200 450*numRows];

if plotOn
for ii = 1:numPsychFxns

    fInd    = figInds(ii);

    f{fInd} = figure(fInd);
    f{fInd}.Position = position + [500 0 0 0]*(fInd-1);

    theseStim = stimStruct(ii).respData(:,1);
    respP  = stimStruct(ii).respData(:,2)./stimStruct(ii).respData(:,3);
    numInp = stimStruct(ii).respData(:,3);

    % Cull bins with no trials
    respP     = respP(numInp ~= 0);
    theseStim = theseStim(numInp ~= 0);
    numInp    = numInp(numInp ~= 0);

    subplot(numRows,numColumns,spInds(ii));
    hold on;

    plot(theseStim,stimStruct(ii).pFxn,'k','linewidth',2);
    plot(theseStim,stimStruct(ii).pFxnBayesSG,'g','linewidth',2);
    plot(theseStim,stimStruct(ii).pFxnBayes,'b','linewidth',2);
    scatter(theseStim,respP,30*numInp,'k');
    plot(getLogXform(stimStruct(ii).trialType(1),0.3)*[1 1],[0 1],'--k');

    set(gca,'fontsize',15,'plotboxaspectratio',[1 1 1],'xtick',speedTick,'xticklabel',sTickLab,'xlim',speedLims,'ylim',[0 1]);
    ylabel('p("test faster")');
    xlabel('Test speed (deg/s)');
    title(['C_{test}=',num2str(stimStruct(ii).trialType(2)),', D_{ref}=',num2str(stimStruct(ii).trialType(4)),...
           ', D_{test}=',num2str(stimStruct(ii).trialType(5))]);
end

for ii = 1:5
     % Save figures
%     saveas(f{ii},['/home/tyler/Documents/MATLAB/cooperLab/1-Experiments_stimuli/dataAnalysisCode/RDS/figs/psychFxn',num2str(ii),'.svg']);
end

end


end