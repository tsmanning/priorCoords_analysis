function [fit,respData,stimStruct] = analyzeExpRDS(subjID,saveFigs,saveData)

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
datFiles  = fileNames(arrayfun(@(x)x.bytes>250000,files));

% Load in data and assign block names
for ii = 1:numel(datFiles)
    
    blockID{ii} = datFiles{ii}(end-4);
    
    switch blockID{ii}
        case '1'
            blockType{ii} = 'within';
        case '2'
            blockType{ii} = 'between';
    end
    
    data(ii)       = load([subjDir,datFiles{ii}],'pa','ds');
    staircases(ii) = load([subjDir,datFiles{ii}],'gSCell');
    
end

% Initial data extraction/sorting
for ii = 1:numel(data)
    
    % Select one dataset
    thisRespMat      = data(ii).pa.response;
    
    % Sort data into stimulus conditions
    numStair(ii) = data(ii).pa.numStaircases;
    
    [contCombs,contInd]     = unique(thisRespMat(:,2));             % test contrast (reference is other)
    [testScr]               = unique(thisRespMat(:,6));             % presentation screen (test distance)
    [refVels]               = unique(thisRespMat(:,5));             % reference velocity
    [stairTypes,~,stairIDs] = unique(thisRespMat(:,11:12),'row');   % staircase target %
    
    % Find number of distinct staircase types (e.g. 1u/1d, 2u/1d)
    numStairTypes(ii)       = size(stairTypes,1);
    
    % Group staircase triplets together
    stimInds     = makeCombos([numel(contInd) numel(testScr) numel(refVels)]);
    numUniqStim  = size(stimInds,1);
    
    for jj = 1:numUniqStim
        
        theseInds = (thisRespMat(:,2) == contCombs(stimInds(jj,1))) & ...
                    (thisRespMat(:,6) == testScr(stimInds(jj,2))) & ...
                    (thisRespMat(:,5) == refVels(stimInds(jj,3)));
        
        if sum(theseInds) < 25        
                keyboard
        end
        
        respMat{ii,jj} = [thisRespMat(theseInds,:) stairIDs(theseInds)];
        stimStruct(ii,jj).blockType = blockType{ii};
        stimStruct(ii,jj).testCont  = contCombs(stimInds(jj,1));
        stimStruct(ii,jj).testDist  = data(ii).pa.distances(testScr(stimInds(jj,2)));
        stimStruct(ii,jj).refVel    = refVels(stimInds(jj,3));
        q(ii,jj).respMat   = respMat{ii,jj};
    end
    
end


%% Plot staircases, psychometric functions, and parameter estimates for different conditions

% Stimulus pairings
stimPairs1  = [1 3;...  % Low Cont; low speed
              2 4;...  % High Cont; low speed
              5 7;...  % Low Cont; high speed
              6 8];    % High Cont; high speed

stimPairs2  = [1 2;...  % High Cont; low speed
              3 4];    % High Cont; High speed   
     
     
% Plotting pars
speedTick = [0.1 0.5 1 2 4 8 12 16 20];
colorMat  = colororder;
xticklabs = {'0.5m','1m','0.5m','1m','0.5m','1m','0.5m','1m'};
leglabs   = {'C_{low}, V_{low}','C_{high}, V_{low}','C_{low}, V_{high}','C_{high}, V_{high}'};
stimLabs  = {'CLo_SLo','CHi_SLo','CLo_SHi','CHi_SHi'};

% Loop over task blocks (within/between)
for ii = 1:numel(data)
    
    if sum(cellfun(@(x) isempty(x),respMat(ii,:))) > 1
        stimPairs = stimPairs2;
    else
        stimPairs = stimPairs1;
    end
    
    numStimPairs = size(stimPairs,1);
    
    % Initialize PSE plots
    f2{ii} = figure;
    f2{ii}.Position = [1050 (300+ii*100) 810 400];
    hold on;
    
    % Initialize slope plots
    f4{ii} = figure;
    f4{ii}.Position = [1250 (300+ii*100) 810 400];
    hold on;
    
    % Loop over stimulus pairings (contrast level/reference speed)
    for jj = 1:numStimPairs
        
        % Intialize staircase/psychometric fxn plot
        f1{ii,jj} = figure;
        f1{ii,jj}.Position = [(50 + jj*100) (50 + ii*100) 1100 1075];
        hold on;
        
        figLab{ii,jj} = ['_',blockType{ii},'_',stimLabs{jj}];
        
        thisStim = stimPairs(jj,:);
        
        % Loop over elements of pair (i.e. near/far)
        for kk = 1:numel(thisStim)
            
            % Grab response matrix
            theseResps    = respMat{ii,thisStim(kk)};
            
            % Grab staircase IDs
            theseStairIDs = respMat{ii,thisStim(kk)}(:,end);
            
            % Define reference velocity for this set of responses
            vRef          = theseResps(1,5);
            
            % Record indexing
            stimStruct(ii,thisStim(kk)).cellInd = [ii,jj,kk];
            
            switch vRef
                % Define speed limits for plotting
                case 1
                    speedLims    = [0.1 8];
                case 8
                    speedLims    = [4 18];
            end
            
            % Get test speeds for each staircase trial
            %-------------------------------------------------------------------------------------%
            stairVels     = cell(numStairTypes(ii),1);
            reversalInds  = cell(numStairTypes(ii),1);
            thresh        = nan(numStair(ii),1);
            numTrials     = nan(numStair(ii),1);
            p             = nan(numStair(ii),1);
            
            for sid = 1:numStairTypes(ii)
                
                % Get indices of all trials run for this staircase
                temp1           = find(theseStairIDs == sid);
                thisStairInd    = theseResps(temp1(1),10);
                
                % Get all test velocities presented in this staircase
                stairVels{sid}  = get(staircases(ii).gSCell{thisStairInd},'values');
                
                % Find trial indices where a reversal occurred
                reversalInds{sid} = logical([get(staircases(ii).gSCell{thisStairInd},'reversals')]);
                
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
            
            % Plot test speeds for all trials run in this staircase
            %-------------------------------------------------------------------------------------%
            set(0,'CurrentFigure',f1{ii,jj});
            hold on;
            
            subplot(2,2,2*(kk-1)+1);
            hold on;
            
            cind = 1;
            
            for sid = 1:numStairTypes
                % Plot all trials
                p(cind)     = scatter(1:numel(stairVels{sid}),stairVels{sid},30,colorMat(sid,:));
                
                % Plot reversal trials in filled dots
                trialInds = 1:numel(stairVels{sid});
                scatter(trialInds(reversalInds{sid}),stairVels{sid}(reversalInds{sid}),30,...
                        colorMat(sid,:),'filled');
                plot([0 max(numTrials)],thresh(sid)*[1 1],'color',colorMat(cind,:));
                cind = cind + 1;
            end
            
            % Plot reference velocity
            p(4) = plot([0 max(numTrials)],theseResps(1,5)*[1 1],'--k');
            
            switch theseResps(1,5)
                case 8
                    locale = 'southeast';
                case 1
                    locale = 'northeast';
            end
            
%             legend([p(1) p(2) p(3) p(4)],{'1u/1d','1u/2d','2u/1d','Ref. vel.'},'location',locale);
            xlabel('Trial');
            ylabel('Test speed (deg/s)');
            title(['V_{ref}=',num2str(theseResps(1,5)),'\circ/s; C_{ref}=',num2str(theseResps(1,4)),...
                   '; D_{test}=',num2str(data(ii).ds.screenDistance(theseResps(1,6))),'m']);
            set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'yscale','log',...
                'ytick',speedTick,'ylim',speedLims,'xlim',[0 max(numTrials)]);
            
            % Plot responses and psychometric functions
            %-------------------------------------------------------------------------------------%
            
            % Recode responses in terms of whether test stimulus was seen faster
            [uniqueVel,~,uniInds] = unique(theseResps(:,3));
            respDir   = theseResps(:,9);
            
            switch blockType{ii}
                case 'within'
                    % Check reference position within screen
                    refDir = theseResps(:,8);
                    xlab   = 'Viewing distance';
                case 'between'
                    % Check reference screen
                    refDir = ~(theseResps(:,6)-1) + 1;
                    xlab   = 'Reference distance';
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
            
            % Fit a cumulative Gaussian to the data
%             [pFxnA,PSEA,sensA] = fitCumGauss(uniqueVel,respP,numInp/sum(numInp));
            
            % Get CIs on fit parameters
            respData{ii,jj,kk} = [uniqueVel respP.*numInp numInp];
            
            options                = struct;
            options.confP          = 0.95;
            options.estimateType   = 'MAP';
            options.expType        = 'YesNo';
            options.sigmoidName    = 'norm';
            options.stepN          = [80,80,1,1,1];
            options.threshPC       = 0.5;
            options.fixedPars      = [nan;nan;0;0;0];
            options.nblocks        = size(respData{ii,jj,kk},1);
            options.maxBorderValue = exp(-20);
            
            result  = psignifit(respData{ii,jj,kk},options);
            pFxn    = result.psiHandle(uniqueVel);
            
            fit{ii,jj,kk} = result;
            fit{ii,jj,kk}.condition = [];
            
            % Plot fit and resampled fits
            subplot(2,2,2*(kk-1)+2);
            hold on;
            
            plot(uniqueVel,pFxn,'k','linewidth',2);
            
            % Calculate where each staircase would fall on psychometric function
            s1 = 0.5;
            s2 = sqrt(0.5);
            s3 = 1 - sqrt(0.5);
            
            % Plot all response proportions,thresholds fit to staircases, and psychometric functions
            scatter(uniqueVel,respP,30*numInp,'k');
            scatter(thresh(1),[s1],100,colorMat(1,:),'filled');
            scatter(thresh(2),[s2],100,colorMat(2,:),'filled');
            scatter(thresh(3),[s3],100,colorMat(3,:),'filled');
            plot(theseResps(1,5)*[1 1],[0 1],'--k');
            set(gca,'xscale','log','fontsize',20,'xtick',speedTick,'xlim',speedLims);
            ylabel('p("test faster")');
            xlabel('Test speed (deg/s)');
            title(blockType{ii});
            set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
            
            % Plot estimates of PSE and slope with CIs
            PSE   = result.Fit(1);
            PSECI = result.conf_Intervals(1,:);
            
            % width = SD for cumulative Gaussian fit
            width   = result.Fit(2);
            widthCI = result.conf_Intervals(2,:);
            
            slope   = (1/width)*(1/sqrt(2*pi));
            slopeCI = (1./widthCI)*(1/sqrt(2*pi));
            
            fit{ii,jj,kk}.slope   = slope;
            fit{ii,jj,kk}.slopeCI = slopeCI;
            
            %%% Should just output single struct like this with everything
            %%% together as fields
            stimStruct(ii,thisStim(kk)).PSE      = PSE;
            stimStruct(ii,thisStim(kk)).PSECI    = PSECI;
            stimStruct(ii,thisStim(kk)).slope    = slope;
            stimStruct(ii,thisStim(kk)).slopeCI  = slopeCI;
            stimStruct(ii,thisStim(kk)).respData = respData{ii,jj,kk};
            stimStruct(ii,thisStim(kk)).pFxn     = pFxn;
            
            % Define current index for plot (stimulus pairing + element of pair)
            thisPlotInd = 2*(jj-1) + kk;
            
            % PSE plot
            set(0,'CurrentFigure',f2{ii});
            hold on;
            ps{jj,kk} = errorbar(thisPlotInd,PSE,PSECI(1)-PSE,PSECI(2)-PSE,...
                's','MarkerFaceColor',colorMat(jj,:),'MarkerEdgeColor',colorMat(jj,:),'linewidth',3,'color',colorMat(jj,:));
            
            % Slope plot
            set(0,'CurrentFigure',f4{ii});
            hold on;
            sl{jj,kk} = errorbar(thisPlotInd,slope,slope-slopeCI(1),slope-slopeCI(2),...
                's','MarkerFaceColor',colorMat(jj,:),'MarkerEdgeColor',colorMat(jj,:),'linewidth',3,'color',colorMat(jj,:));
            
        end
        
    end
    
    set(0,'CurrentFigure',f2{ii});
    hold on;
    plot([0.5 8.5],1*[1 1],'--k','linewidth',2);
    plot([0.5 8.5],8*[1 1],'--k','linewidth',2);
    title(['PSE (',blockType{ii},')']);
    set(gca,'xlim',[0 numStimPairs*2]+0.5,'xtick',1:numStimPairs*2,'xticklabel',xticklabs,'ytick',speedTick,'ylim',speedLims,'yscale','log','fontsize',20,'plotboxaspectratio',[4 2 1]);
    xlabel(xlab);
    ylabel('Test speed (deg/s)');
%     legend([ps{1,1} ps{2,1} ps{3,1} ps{4,1}],leglabs,'location','southeast');
    
    set(0,'CurrentFigure',f4{ii});
    hold on;
    title(['Slope (',blockType{ii},')']);
    set(gca,'xlim',[0 numStimPairs*2]+0.5,'xtick',1:numStimPairs*2,'xticklabel',xticklabs,'ylim',[0 1],'fontsize',20,'plotboxaspectratio',[4 2 1]);
    xlabel(xlab);
    ylabel('\Delta p/\Delta V');
%     legend([sl{1,1} sl{2,1} sl{3,1} sl{4,1}],leglabs,'location','northeast');
    
end

if saveFigs
    
%     figLab = figLab(:);
%     
%     if ~exist([subjDir,'figures',filesep])
%         mkdir([subjDir,'figures',filesep]);
%     end
%     
%     for ii = 1:numel(f1)
%         try
%             saveas(f1{ii},[subjDir,'figures/f1_',figLab{ii},'.svg']);
%         catch
%             continue
%         end
%     end
%     
%     for ii = 1:numel(f2)
%         saveas(f2{ii},[subjDir,'figures/f2_',num2str(ii),'.svg']);
%     end
%     
%     for ii = 1:numel(f4)
%         saveas(f4{ii},[subjDir,'figures/f4_',num2str(ii),'.svg']);
%     end
%     
end

if saveData
    
    if ~exist([dataDir,'xSub',filesep])
        mkdir([dataDir,'xSub',filesep]);
    end
    
    save([dataDir,'xSub',filesep,subjID],'respData','fit','stimStruct');
    
end

end