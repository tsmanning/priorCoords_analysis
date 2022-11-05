% Run analysis pipeline for Late Summer/Fall 2022 iteration of experiment

clear all
close all

addpath(genpath('/home/tyler/Documents/MATLAB/cooperLab/1-Experiments_stimuli/dataAnalysisCode/RDS'))
rmpath(genpath('/home/tyler/Documents/MATLAB/cooperLab/1-Experiments_stimuli/dataAnalysisCode/RDS/old'));

saveDir = '/home/tyler/Documents/MATLAB/cooperLab/3-Data/priorCoordsRDS/xSub/';

singleSub    = 0;
plotSubPFxns = 0;
xSub         = 0;
ssBayes      = 0;
ssBayesPlot  = 1;
xSubBayes    = 0;


%% Define subject pool
% Clearly good
% subjIDs = {'PCJ','DJB','FHP','MLZ','JTD','BCG','JSG','GHT','WIE','UCK',...
%            'FTB','MUB','KBU'};

% Clearly good + questionable
% subjIDs = {'PCJ','DJB','FHP','RAY','MLZ','JTD','YMY','ZDG','BCG','JSG',...
%            'AUI','GHT','GDZ','WIE','UCK','FTI','FTB','MUB','KBU'};

% Complete datasets
% subjIDs = {'PCJ','FHP','MLZ','JSG','AUI','GDZ','PSZ','WIE','UCK','YKB','FTB','FTI','KBU','IXK'};

% Left screen nears
% subjIDs = {'PCJ','DJB','FHP','KTP','RAY','GDZ','WIE','UCK','YKB','FTI','FTB'};

% Right screen near
% subjIDs = {'MLZ','JTD','YMY','BCG','JSG','AUI'};

% All
% subjIDs = {'PCJ','LGL','DJB','FHP','XAZ','KTP','RAY','MLZ','JTD','YMY',...
%            'ZDG','BCG','JSG','AUI','GHT','GDZ','PSZ','WIE','UCK','YKB',...
%            'FTI','FTB','MUB','HCL','IXK','KBU'};

% Third or more bad p Fxns
% subjIDs = {'LGL','XAZ','KTP','PSZ','YKB','HCL','IXK'};

% Testing
subjIDs = {'BIV'};

% NOTES:
% DJB has weird between screen results
numSubs = numel(subjIDs);


%% Run data extraction and single subject analysis
if singleSub
    
    saveOn = 1;

    for ii = 1:numSubs
        
        [data] = analyzeExpRDS(subjIDs{ii},saveOn);

    end
end

%% Plot a subject's psychometric curves
if plotSubPFxns

    saveOn = 0;

    for ii = 1:numSubs

        subjID = subjIDs{ii};

        plotAllPsychFxns(saveDir,subjID,saveOn);

    end
end

%% Cross subject frequentist analysis
if xSub

    saveOn = 0;

%     stats = analyzeXsubRDS2(saveOn,subjIDs,saveDir);
    stats = analyzeXsubRDS2(saveOn,subjIDs);

end


%% Run Bayesian analysis on single subjects
if ssBayes

    % Get IPDs
    [ipd] = loadInIPD(subjIDs,saveDir);

    % Define velocity support
    xgrid = getLogXform(linspace(0,16,100),0.3);

    % Define number of starting points in obj function space
    numStartPts = 5;

    % Choose numerical/analytical forms of MoG model
    % solverOpts.fxnForm = 'numerical';
    solverOpts.fxnForm = 'analytical';

    % Choose optimization solver
    solverOpts.solver = 'fmincon';
    % solverOpts.solver = 'bads';

    % Define which trial types we want to fit ideal observer with
    incTrTypes = 'Within';
%     incTrTypes = 'WithinBetweenPers';
%     incTrTypes = 'WithinBetweenAll';
%     incTrTypes = 'BetweenAll';

    saveOn = 1;

    for ii = 1:numSubs

        % Grab sorted data and parameters
        [sortedPars,dataStruct] = sortBayesData(saveDir,subjIDs{ii},incTrTypes);
        bi        = sortedPars.inds{1};
        cii       = sortedPars.inds{2};
        sortInds1 = sortedPars.inds{3};
        sortInds2 = sortedPars.inds{4};

        % Loop over model types (geometric/nonparametric)
        for jj = 1:2

            if jj == 1
                % Fit ideal observer with noise scaling parameters
                geoOn = 0;
            else
                % Fit ideal observer with geometric noise propagation model
                geoOn = 1;
            end

            % Run for single Gaussian
            nB = 1;
            [outMatSG] = fit2AFCData(sortedPars,ipd(ii),nB,xgrid,geoOn,numStartPts,solverOpts);

            % Run for MoG
            nB = 2;
            [outMatMoG] = fit2AFCData(sortedPars,ipd(ii),nB,xgrid,geoOn,numStartPts,solverOpts);

            % Collect unique psychometric function values
            fitPM = outMatMoG.p(bi)';
            fitPS = outMatSG.p(bi)';

            for kk = 1:numel(dataStruct)

                % Collect individual psychometric function values into complete
                % fxns and place into dataStruct
                theseInds = cii == sortInds1(kk);

                dsInd = sortInds2(kk);

                dataStruct(dsInd).pFxnBayes   = fitPM(theseInds);
                dataStruct(dsInd).pFxnBayesSG = fitPS(theseInds);

                % Extract estimated distances from distance component of
                % likelihood
                refTestDInd = [find(dataStruct(kk).refDist==outMatSG.distVals),...
                               find(dataStruct(kk).testDist==outMatSG.distVals)];

                if ~isempty(outMatMoG.pSigDistHat)
                    % Geo = 0
                    estDistSG    = sqrt(2)*[outMatSG.pSigDistHat(refTestDInd(2)) outMatSG.pSigDistHat(refTestDInd(1))].^2;
                    estDistMoG   = sqrt(2)*[outMatMoG.pSigDistHat(refTestDInd(2)) outMatMoG.pSigDistHat(refTestDInd(1))].^2;
                else
                    % Geo = 1
                    %%%%%%% note this is a simplification and we really
                    %%%%%%% should do this properly (assumes the stimulus
                    %%%%%%% is along midsaggital plane, when this assumption wasn't used to fit)
                    estDistSG    = sqrt(2*[outMatSG.stimD_SF(refTestDInd(2)) outMatSG.stimD_SF(refTestDInd(1))]);
                    estDistMoG   = sqrt(2*[outMatMoG.stimD_SF(refTestDInd(2)) outMatMoG.stimD_SF(refTestDInd(1))]);
                end

                dataStruct(ii).estDistSG    = estDistSG;
                dataStruct(ii).estDistMoG   = estDistMoG;
            end

            % Extract estimated distances

            if saveOn
                % Save fits
                save([saveDir,subjIDs{ii},'geoOn',num2str(geoOn),'_',incTrTypes,'.mat'],...
                    'dataStruct','sortedPars','outMatMoG','outMatSG');
            end

        end

    end

end

%% Plot psychometric functions with Bayesian and Gaussian fits, Bayesian parameters
if ssBayesPlot

    saveOn = 1;
    
    % Define which trial types we want to plot fits for
    incTrTypes = 'Within';
%     incTrTypes = 'WithinBetweenPers';
%     incTrTypes = 'WithinBetweenAll';
%     incTrTypes = 'BetweenAll';

    for ii = 1:numSubs
        
        close all

        for jj = 1:2

            if jj == 1
                geoOn = 0;
            else
                geoOn = 1;
            end

            % Grab response data and sort into stimV/C/D/r
            load([saveDir,subjIDs{ii},'geoOn',num2str(geoOn),'_',incTrTypes,'.mat']);

            % Best fit cum. Gauss and Bayesian psychometric fxns
            plotBayesPsychFxns(subjIDs{ii},geoOn,dataStruct,saveDir,saveOn);

            % Best fit Bayesian parameters
            plotBayesFit(subjIDs{ii},outMatSG,outMatMoG,dataStruct,saveDir,geoOn,saveOn);

        end
    end

end


%% Cross subject Bayesian analysis
if xSubBayes

    saveOn = 0;

    % Define which trial types we want to plot fits for
    incTrTypes = 'Within';
%     incTrTypes = 'WithinBetweenPers';
%     incTrTypes = 'WithinBetweenAll';
%     incTrTypes = 'BetweenAll';

    analyzeXsubBayes(saveOn,subjIDs,saveDir);

end

