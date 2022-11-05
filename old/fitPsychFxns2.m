function [dataStruct,f1,f2] = fitPsychFxns2(subjID,geoOn,dataStruct,fitMatSG,fitMatMoG,saveDir,plotOn)

% Takes data from psychophysics exp, fits with cumulative Gaussian and also
% plots Bayesian fit for comparison

%% Parse Data
% 
% % Number of psychometric functions
% numPsychFxns = numel(dataStruct);
% 

%% Loop over psychometric functions (unique trial types)
% for ii = 1:numPsychFxns
% 
%     % Find stimulus values for this psych function
%     thisRefVel   = dataStruct(ii).refVel;
%     theseTestVel = dataStruct(ii).respData(:,1);
% 
%     % Find parameter indices for this psych function (Ind 1: Ref, ind 2: test)
%     switch dataStruct(ii).blockType
%         case 'within'
%             %%%% Should really just make this a field in analyzeExpRDS.m
%             dataStruct(ii).refCont = fitMatSG.contVals(dataStruct(ii).testCont ~= fitMatSG.contVals);
%         case 'between'
%             dataStruct(ii).refCont = dataStruct(ii).testCont;
%     end
%     refTestCInd = [find(dataStruct(ii).refCont==fitMatSG.contVals) find(dataStruct(ii).testCont==fitMatSG.contVals)];
%     refTestDInd = [find(dataStruct(ii).refDist==fitMatSG.distVals) find(dataStruct(ii).testDist==fitMatSG.distVals)];
% 
%     % Calculate NLL for cumulative Gaussian fit & coin flip
%     nllCG = nan(numel(theseTestVel),1);
%     nllCF = nan(numel(theseTestVel),1);
% 
%     for ll = 1:numel(theseTestVel)
%         inputs = zeros(dataStruct(ii).respData(ll,3),1);
%         inputs(1:dataStruct(ii).respData(ll,2),1) = 1;
%         pCGfit = min([dataStruct(ii).pFxn(ll),1-eps]);  % avoids log(0) below if p=1
%         pCF    = 0.5;
% 
%         nllCG(ll) = -sum(inputs*log(pCGfit) + (1 - inputs)*log(1-pCGfit));
%         nllCF(ll) = -sum(inputs*log(pCF) + (1 - inputs)*log(1-pCF));
%     end
% 
%     nllCG = sum(nllCG);
%     nllCF = sum(nllCF);
% 
%     if ~isempty(fitMatMoG.pSigDistHat)
%         % Geo = 0
%         estDistSG    = sqrt(2)*[fitMatSG.pSigDistHat(refTestDInd(2)) fitMatSG.pSigDistHat(refTestDInd(1))].^2;
%         estDistMoG   = sqrt(2)*[fitMatMoG.pSigDistHat(refTestDInd(2)) fitMatMoG.pSigDistHat(refTestDInd(1))].^2;
%     else
%         % Geo = 1Near 
%         estDistSG    = sqrt(2*[fitMatSG.stimD_SF(refTestDInd(2)) fitMatSG.stimD_SF(refTestDInd(1))]);
%         estDistMoG   = sqrt(2*[fitMatMoG.stimD_SF(refTestDInd(2)) fitMatMoG.stimD_SF(refTestDInd(1))]);
%     end
% 
%     % Output all fits and data as a structure
%     dataStruct(ii).nllCG        = nllCG;
%     dataStruct(ii).nllCF        = nllCF;
%     dataStruct(ii).estDistSG    = estDistSG;
%     dataStruct(ii).estDistMoG   = estDistMoG;
% 
% end

%% Plot psychometric functions (cum Gauss, SG, MoG)

if plotOn

    [types,b]    = sort(arrayfun(@(x) x.blockType,dataStruct,'UniformOutput',false));
    betweenInds  = find(strcmp('between',types)==1);
    withinInds   = find(strcmp('within',types)==1);

    dataStruct = dataStruct(b);

    vels = [0.5 1 2 4 8 16 32];
    for ii = 1:numel(vels)
        velLabs{ii} = num2str(vels(ii));
    end

    if sum(betweenInds)>0
    % Plot Between trials
    f1 = figure;
    f1.Position = [100 100 1075 450];
    n = 1;
    for ii = betweenInds'

        subplot(2,4,n);
        n = n+1;
        hold on;

        linVels = getLinXform(dataStruct(ii).respData(:,1),0.3);
        linRef = dataStruct(ii).refVel;

        scatter(linVels,dataStruct(ii).respData(:,2)./dataStruct(ii).respData(:,3),10*dataStruct(ii).respData(:,3),'k');
        % Cum. Gauss
        plot(linVels,dataStruct(ii).pFxn,'k','linewidth',2);
        % Single Gauss prior
        plot(linVels,dataStruct(ii).pFxnBayesSG,'g','linewidth',2);
        % MoG prior
        plot(linVels,dataStruct(ii).pFxnBayes,'b','linewidth',2);

        plot(linRef*[1 1],[0 1],'--k');
        set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0.1 30],'ylim',[0 1],'xscale','log','xtick',vels,'xticklabel',velLabs);
        title([dataStruct(ii).blockType,', D_{ref}=',num2str(dataStruct(ii).testDist),', D_{test}=',num2str(dataStruct(ii).testDist)]);

    end

    % Save figures
    saveas(f1,[saveDir,'figures/',subjID,'_geo',num2str(geoOn),'_betweenTrialsFits.svg']);
    end

    if sum(withinInds)>0
    % Plot Within trials
    f2 = figure;
    f2.Position = [100 100 1400 930];
    n = 1;
    for ii = withinInds'

        subplot(4,5,n);
        n = n+1;
        hold on;

        linVels = getLinXform(dataStruct(ii).respData(:,1),0.3);
        linRef = dataStruct(ii).refVel;

        scatter(linVels,dataStruct(ii).respData(:,2)./dataStruct(ii).respData(:,3),10*dataStruct(ii).respData(:,3),'k');
        % Cum. Gauss
        plot(linVels,dataStruct(ii).pFxn,'k','linewidth',2);
        % Single Gauss prior
        plot(linVels,dataStruct(ii).pFxnBayesSG,'g','linewidth',2);
        % MoG prior
        plot(linVels,dataStruct(ii).pFxnBayes,'b','linewidth',2);

        plot(linRef*[1 1],[0 1],'--k');
        set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0.1 30],'ylim',[0 1],'xscale','log','xtick',vels,'xticklabel',velLabs);
        title([dataStruct(ii).blockType,', C_{test}=',num2str(dataStruct(ii).testCont),', D=',num2str(dataStruct(ii).testDist)]);

    end

    % Save figures
    saveas(f2,[saveDir,'figures/',subjID,'geo',num2str(geoOn),'_withinTrialsFits.svg']);
    end

end

end