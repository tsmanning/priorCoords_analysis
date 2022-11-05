function [stimStruct] = analyzeExpRDS_noPlot(subjID,saveData)

% Function to analyze data from priorCoord RDS exp

close all

%% Load in and sort data
% Set up directories
splPath = regexp(which('analyzeExpRDS'),filesep,'split');

dataDir = [filesep,fullfile(splPath{1:numel(splPath)-4}),filesep,'3-Data/priorCoordsRDS/'];
subjDir = [dataDir,subjID,filesep];

% Find datsets in subject directory
files     = dir(subjDir);
fileNames = arrayfun(@(x)x.name,files,'uniformoutput',false);
datFiles  = fileNames(arrayfun(@(x)x.bytes>100000,files));

% Load in data and assign block names
for ii = 1:numel(datFiles)
    
    blockID{ii} = datFiles{ii}(end-4);
    
    switch blockID{ii}
        case '1'
            blockType{ii} = 'within';
        case '2'
            blockType{ii} = 'between';
        case '3'
            blockType{ii} = 'sizeCntl';
    end
    
    data(ii)       = load([subjDir,datFiles{ii}],'pa','ds');
    staircases(ii) = load([subjDir,datFiles{ii}],'gSCell');
    
end

% Initial data extraction/sorting
for ii = 1:numel(data)
    
    % Select one dataset
    thisRespMat  = data(ii).pa.response;
    
    % Tack on another column for ref order (1: ref first, 2: ref second)
    if size(thisRespMat,2) == 15
        thisRespMat = [thisRespMat [~data(ii).pa.stimOrder + 1]'];
    end

    % Sort data into stimulus conditions
    numStair(ii) = data(ii).pa.numStaircases;
    
    [contCombs,contInd]     = unique(thisRespMat(:,2));             % test contrast (reference is other)
    [testScr]               = unique(thisRespMat(:,6));             % presentation screen (test distance)
    [refVels]               = unique(thisRespMat(:,5));             % reference velocity
    if size(thisRespMat,2) >= 15
    [refSzs]                = unique(thisRespMat(:,15));             % reference sizes
    end
    [stairTypes,~,stairIDs] = unique(thisRespMat(:,11:12),'row');   % staircase target %
    
    % Find number of distinct staircase types (e.g. 1u/1d, 2u/1d)
    numStairTypes(ii)       = size(stairTypes,1);
    
    % Group staircase triplets together
    stimInds     = makeCombos([numel(contInd) numel(testScr) numel(refVels)]);

    if size(thisRespMat,2) >= 15
    stimInds     = makeCombos([numel(contInd) numel(testScr) numel(refVels) numel(refSzs)]);
    end
    numUniqStim  = size(stimInds,1);
    
    for jj = 1:numUniqStim
        
        if size(thisRespMat,2) < 15
            theseInds = (thisRespMat(:,2) == contCombs(stimInds(jj,1))) & ...
                        (thisRespMat(:,6) == testScr(stimInds(jj,2))) & ...
                        (thisRespMat(:,5) == refVels(stimInds(jj,3)));


            respMat{ii,jj} = [thisRespMat(theseInds,:) stairIDs(theseInds)];
            stimStruct(ii,jj).blockType = blockType{ii};
            stimStruct(ii,jj).testCont  = contCombs(stimInds(jj,1));
            stimStruct(ii,jj).testDist  = data(ii).pa.distances(testScr(stimInds(jj,2)));
            stimStruct(ii,jj).refVel    = refVels(stimInds(jj,3));
            stimStruct(ii,jj).respMat   = respMat{ii,jj};
        else
            theseInds = (thisRespMat(:,2) == contCombs(stimInds(jj,1))) & ...
                        (thisRespMat(:,6) == testScr(stimInds(jj,2))) & ...
                        (thisRespMat(:,5) == refVels(stimInds(jj,3)));
                        (thisRespMat(:,15) == refSzs(stimInds(jj,4)));

            respMat{ii,jj} = [thisRespMat(theseInds,:) stairIDs(theseInds)];
            stimStruct(ii,jj).blockType = blockType{ii};
            stimStruct(ii,jj).testCont  = contCombs(stimInds(jj,1));
            stimStruct(ii,jj).testDist  = data(ii).pa.distances(testScr(stimInds(jj,2)));
            stimStruct(ii,jj).refVel    = refVels(stimInds(jj,3));
            stimStruct(ii,jj).refSz    = refSzs(stimInds(jj,4));
            stimStruct(ii,jj).respMat   = respMat{ii,jj};
        end

    end
    
end


%% Plot staircases, psychometric functions, and parameter estimates for different conditions

% Loop over task blocks (within/between)
for ii = 1:numel(data)
    
    thisNumStim = sum(cellfun(@(x) ~isempty(x),respMat(ii,:)));
    
    % Loop over staircase triplets
    for jj = 1:thisNumStim
        
        % Grab response matrix
        theseResps    = respMat{ii,jj};
        
        % Grab staircase IDs
        theseStairIDs = respMat{ii,jj}(:,end);
        
        % Define reference velocity for this set of responses
        vRef          = theseResps(1,5);
        
        % Get test speeds for each staircase trial
        %-------------------------------------------------------------------------------------%
        stairVels     = cell(numStairTypes(ii),1);
        reversalInds  = cell(numStairTypes(ii),1);
        thresh        = nan(numStairTypes(ii),1);
        numTrials     = nan(numStairTypes(ii),1);
        
        for sid = 1:numStairTypes(ii)
            
            % Get indices of all trials run for this staircase
            temp1           = find(theseStairIDs == sid);
            thisStairInd    = theseResps(temp1(1),10);
            
            if ~isstruct(staircases(ii).gSCell{thisStairInd})
                % Get all test velocities presented in this staircase
                stairVels{sid}  = get(staircases(ii).gSCell{thisStairInd},'values');
                
                % Find trial indices where a reversal occurred
                reversalInds{sid} = logical([get(staircases(ii).gSCell{thisStairInd},'reversals')]);
            else
                stairVels{sid}  = staircases(ii).gSCell{thisStairInd}.values;
                reversalInds{sid} = logical(staircases(ii).gSCell{thisStairInd}.reversals');
            end

            % Calculate threshold as mean of speeds after 1st reversal
            temp2       = find(reversalInds{sid} == 1);
            
            % Sometimes there aren't many reversals before the max trials are
            % hit - let's just use last third of the reversals
            meaningInds = (length(temp2) - round(length(temp2)/3)):length(temp2);
            
            if sum(meaningInds) ~= 0
                thresh(sid) = mean(stairVels{sid}(temp2(meaningInds)));
            else
                % One subject just had no reversals??
                thresh(sid) = nan;
            end
            
            % Get number of trials run for this staircase
            numTrials(sid) = numel(stairVels{sid});
            
        end
        
        % Recode responses in terms of whether test stimulus was seen faster
        [uniqueVel,~,uniInds] = unique(theseResps(:,3));
        respDir   = theseResps(:,9);
        
        switch blockType{ii}
            case 'within'
                % Check reference position within screen
                refDir = theseResps(:,8);
            case 'sizeCntl'
                if size(theseResps,2) <= 15
                    % Check reference position within screen
                    refDir = theseResps(:,8);
                else
                    % Check reference order
                    refDir = theseResps(:,16);
                end
            case 'between'
                % Check reference screen
                refDir = ~(theseResps(:,6)-1) + 1;
        end
        
        respStim  = respDir ~= refDir;
        
        % Calculate proportion of "test faster" responses for each unique velocity
        respP  = nan(numel(uniqueVel),1);
        numInp = nan(numel(uniqueVel),1);
        
        for ll = 1:numel(uniqueVel)
            
            inputs     = respStim(uniInds == ll);
            respP(ll)  = sum(inputs)/numel(inputs);
            numInp(ll) = numel(inputs);
            
        end
        
        % Fit data with Bayesian ideal observer model (w/ prior made of 1+
        % Gaussians)
        

        % Fit a cumulative Gaussian to the data
        [pFxnA,PSEA,widthA] = fitCumGauss(uniqueVel,respP,numInp/sum(numInp));
        slopeA              = (1/widthA)*(1/sqrt(2*pi));
        
        % Fit cum. Gauss with psignifit & get CIs on fit parameters
        respData{ii,jj} = [getLogXform(uniqueVel,0.3) respP.*numInp numInp];
        
        options                = struct;
        options.confP          = 0.95;
        options.estimateType   = 'MAP';
        options.expType        = 'YesNo';
        options.sigmoidName    = 'norm';
        options.stepN          = [80,80,1,1,1];
        options.threshPC       = 0.5;
        options.fixedPars      = [nan;nan;0;0;0];
        options.nblocks        = size(respData{ii,jj},1);
        options.maxBorderValue = exp(-20);
        
        result  = psignifit(respData{ii,jj},options);
        pFxn    = result.psiHandle(getLogXform(uniqueVel,0.3));
        
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
        
        % Output all fits and data as a structure
        stimStruct(ii,jj).PSE          = PSE;
        stimStruct(ii,jj).PSECI        = PSECI;
        stimStruct(ii,jj).slope        = slope;
        stimStruct(ii,jj).slopeCI      = slopeCI;
        stimStruct(ii,jj).respData     = respData{ii,jj};
        stimStruct(ii,jj).pFxn         = pFxn;
        stimStruct(ii,jj).stairVels    = stairVels;
        stimStruct(ii,jj).reversalInds = reversalInds;
        stimStruct(ii,jj).thresh       = thresh;
        stimStruct(ii,jj).pFxnA        = pFxnA;
        stimStruct(ii,jj).PSEA         = PSEA;
        stimStruct(ii,jj).slopeA       = slopeA;
        
    end
    
end

if saveData
    
    if ~exist([dataDir,'xSub',filesep])
        mkdir([dataDir,'xSub',filesep]);
    end
    
    save([dataDir,'xSub',filesep,subjID],'stimStruct');
    
end

end