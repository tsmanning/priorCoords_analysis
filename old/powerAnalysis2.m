% close all
clear all

%% Load in data

datDir = '/media/tyler/Data/MATLAB/cooperLab/3-Data/priorCoordsRDS/';
% datStr = 'tyler/tyler-20220311T170416-1.mat';
% datStr = 'tyler/tyler-20220314T165526-1.mat';
% datStr = 'CTW/CTW-20220316T144044-1.mat';
% datStr = 'CTW/CTW-20220316T151417-2.mat';
% datStr = 'tyler/tyler-20220318T182022-1.mat';
% datStr = 'tyler/tyler-20220321T130518-1.mat';
% datStr = 'tyler/tyler-20220321T153604-2.mat';
% datStr = 'tyler/tyler-20220321T164350-1.mat';
datStr = 'MJX/MJX-20220323T164126-1.mat';

load([datDir,datStr]);

blockID = datStr(length(datStr)-4);

% Response Mat
respMat       = pa.response;

% Sort data into stimulus conditions
numStair      = pa.numStaircases;

[contCombs,contInd,~]   = unique(respMat(:,2));             % contrast combinations
[testScr,~,~]           = unique(respMat(:,6));             % presentation screen
[refVels,refVInd,~]     = unique(respMat(:,5));             % reference velocity
[stairTypes,~,stairIDs] = unique(respMat(:,11:12),'row');   % staircase target %

numStairTypes = size(stairTypes,1);

% Group staircase triplets together
stimInds     = makeCombos([numel(contInd) numel(testScr) numel(refVels)]);
numUniqStim  = size(stimInds,1);

for ii = 1:numUniqStim
    
    theseInds = (respMat(:,2) == contCombs(stimInds(ii,1))) & ...
                (respMat(:,6) == testScr(stimInds(ii,2))) & ...
                (respMat(:,5) == refVels(stimInds(ii,3)));
    
    resps{ii} = [respMat(theseInds,:) stairIDs(theseInds)];
    
end

% Which stimulus to plot?
% stim         = [1 3];     % Low Cont; low speed
% stim         = [2 4];     % High Cont; low speed
stim         = [5 7];     % Low Cont; high speed
% stim         = [6 8];     % High Cont; high speed

f1 = figure;
f1.Position = [100 100 1100 1075];
hold on;

f2 = figure;
f2.Position = [1050 300 810 400];
hold on;

for conds = 1:numel(stim)
    
    theseResps    = resps{stim(conds)};
    theseStairIDs = resps{stim(conds)}(:,end);
    vRef          = theseResps(1,5);
    
    %% Plot test speeds for each staircase trial
    stairVels    = cell(numStairTypes,1);
    reversalInds = cell(numStairTypes,1);
    thresh       = nan(numStair,1);
    numTrials    = nan(numStair,1);
    p            = nan(numStair,1);
    
    for sid = 1:numStairTypes
        
        % Get indices of all trials run for this staircase
        temp1           = find(theseStairIDs == sid);
        thisStairInd    = theseResps(temp1(1),10);
        
        % Get all test velocities presented in this staircase
        stairVels{sid}  = get(gSCell{thisStairInd},'values');
        
        % Find trial indices where a reversal occurred
        reversalInds{sid} = logical([get(gSCell{thisStairInd},'reversals')]);
        
        % Calculate threshold as mean of speeds after 1st reversal
        temp2       = find(reversalInds{sid} == 1);
        
        % Sometimes there aren't many reversals before the max trials are
        % hit - determine where we're meaning based on percentage
        meaningInds = (length(temp2) - round(length(temp2)/3)):length(temp2);
        
        thresh(sid) = mean(stairVels{sid}(temp2(meaningInds)));
        
        % Get number of trials run for this staircase
        numTrials(sid) = numel(stairVels{sid});
        
    end
    
    % Plotting pars
    speedTick    = [0.1 0.5 1 2 4 8 12 16 20];
%     speedLims    = [min([theseResps(:,3);theseResps(1,5)]) 15];
    switch vRef
        case 1
            speedLims    = [0.25 8];
        case 8
            speedLims    = [4 18];
    end
    colorMat     = colororder;
    
    % Plot test speeds for all trials run in this staircase
    set(0,'CurrentFigure',f1);
    hold on;
    
    subplot(2,2,2*(conds-1)+1);
    hold on;
    
    cind = 1;
    
    for sid = 1:numStairTypes
        % Plot all trials
        p(cind)     = scatter(1:numel(stairVels{sid}),stairVels{sid},30,colorMat(sid,:));
        
        % Plot reversal trials in filled dots
        trialInds = 1:numel(stairVels{sid});
        scatter(trialInds(reversalInds{sid}),stairVels{sid}(reversalInds{sid}),30,colorMat(sid,:),'filled');
        plot([0 max(numTrials)],thresh(sid)*[1 1],'color',colorMat(cind,:));
        cind = cind + 1;
    end
    
    p(4) = plot([0 max(numTrials)],theseResps(1,5)*[1 1],'--k');
    
    switch theseResps(1,5)
        case 8
            locale = 'southeast';
        case 1
            locale = 'northeast';
    end
    
    legend([p(1) p(2) p(3) p(4)],{'1u/1d','1u/2d','2u/1d','Ref. vel.'},'location',locale);
    xlabel('Trial');
    ylabel('Test speed (deg/s)');
    title(['V_{ref}=',num2str(theseResps(1,5)),'\circ/s; C_{ref}=',num2str(theseResps(1,2)),'; Dist.=',num2str(ds.screenDistance(theseResps(1,6))),'m']);
    set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'yscale','log',...
        'ytick',speedTick,'ylim',speedLims,'xlim',[0 max(numTrials)]);
    
    
    %% Plot responses and psychometric functions
    
    % Recode responses in terms of whether test stimulus was seen faster
    [uniqueVel,~,uniInds] = unique(theseResps(:,3));
    respDir   = theseResps(:,9);
    switch blockID
        case '1'
            % Check reference position within screen
            refDir = theseResps(:,8);
            xlab   = 'Viewing distance';
        case '2'
            % Check reference screen
            refDir = ~(theseResps(:,6)-1) + 1;
            xlab   = 'Reference distance';
    end
    
    respStim  = respDir ~= refDir;

    % Calculate proportion of "test faster" responses for each unique velocity
    respP  = nan(numel(uniqueVel),1);
    numInp = nan(numel(uniqueVel),1);
    
    for ii = 1:numel(uniqueVel)
        
        inputs     = respStim(uniInds == ii);
        respP(ii)  = sum(inputs)/numel(inputs);
        numInp(ii) = numel(inputs);
        
    end
    
    % Fit a cumulative Gaussian to the data
    [pFxnA,PSEA,sensA] = fitCumGauss(uniqueVel,respP,numInp/sum(numInp));
    
    % Get CIs on fit parameters
    data = [uniqueVel respP.*numInp numInp];
    
    options                = struct;
    options.confP          = 0.95;
    options.estimateType   = 'MAP';
    options.expType        = 'YesNo';
    options.sigmoidName    = 'norm';
    options.stepN          = [80,80,1,1,1];
    options.threshPC       = 0.5;
    options.fixedPars      = [nan;nan;0;0;0];
    options.nblocks        = size(data,1);
    
    result   = psignifit(data,options); 
    pFxn     = result.psiHandle(uniqueVel);
    
    fit{conds} = result;
    
    % Plot fit and resampled fits
    subplot(2,2,2*(conds-1)+2);
    hold on;
    
    plot(uniqueVel,pFxn,'k','linewidth',2);
%     plot(uniqueVel,pFxnA,'--k','linewidth',2);
    
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
    set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
    
    
    %% Plot estimates of PSE and slope with CIs
    set(0,'CurrentFigure',f2);
    hold on;
    
    PSE   = result.Fit(1);
    PSECI = result.conf_Intervals(1,:);
    
    % width = SD for cumulative Gaussian fit
    width   = result.Fit(2);
    widthCI = result.conf_Intervals(2,:);
    
    slope   = (1/width)*(1/sqrt(2*pi));
    slopeCI = (1./widthCI)*(1/sqrt(2*pi));
    
    subplot(1,2,1);
    hold on;
    errorbar(conds,PSE,PSECI(1)-PSE,PSECI(2)-PSE,...
        's','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'linewidth',3,'color',[0 0 0]);
%     scatter(conds+0.1,PSEA,'k','filled');
    plot([0.5 2.5],vRef*[1 1],'--k','linewidth',2);
    title('PSE');
    set(gca,'xlim',[0.5 2.5],'xtick',[1 2],'xticklabel',{'0.5m','1m'},'ytick',speedTick,'ylim',theseResps(1,5)*[0.4 2.2],'yscale','log','fontsize',20,'plotboxaspectratio',[1 1 1]);
    xlabel(xlab);
    ylabel('Test speed (deg/s)');
    
    subplot(1,2,2);
    hold on;
    errorbar(conds,slope,slope-slopeCI(1),slope-slopeCI(2),...
        's','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'linewidth',3,'color',[0 0 0]);
%     scatter(conds+0.1,(1/sensA)*(1/sqrt(2*pi)),'k','filled');
    title('Slope');
    set(gca,'xlim',[0.5 2.5],'xtick',[1 2],'xticklabel',{'0.5m','1m'},'ylim',[0 0.5],'fontsize',20,'plotboxaspectratio',[1 1 1]);
    xlabel(xlab);
    ylabel('\Delta p/\Delta V');
    
end


