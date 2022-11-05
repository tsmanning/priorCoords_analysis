function [stats] = analyzeXsubRDS(saveOn,varargin)

% Compare prior coords experiment across subjects (RDS variant)


%% Load in and sort data
% Set up directories
splPath = regexp(which('analyzeXsubRDS'),filesep,'split');

dataDir = [filesep,fullfile(splPath{1:numel(splPath)-4}),filesep,'3-Data/priorCoordsRDS/xSub/'];

if nargin < 2
    % Find datsets in subject directory
    files     = dir(dataDir);
    fileNames = arrayfun(@(x)x.name,files,'uniformoutput',false);
    datFiles  = fileNames( cellfun(@(x) contains(x,'mat'),fileNames) );
else
    % Specify which subjects to include in analysis
    datFiles  = varargin{1};
end

% Load in data and assign block names
for ii = 1:numel(datFiles)
    
    % make array of subjIDs
    if numel(datFiles{ii})>3
        subjID{ii} = datFiles{ii}(1:3);
    else
        subjID{ii} = datFiles{ii};
    end
    
    data(ii)            = load([dataDir,datFiles{ii}],'stimStruct');
    data(ii).stimStruct = data(ii).stimStruct(:);
    
    % Cull empty fields
    incInds = arrayfun(@(x) ~isempty(x.testCont),data(ii).stimStruct);
    data(ii).stimStruct = data(ii).stimStruct(incInds);
    
    % Reorganize fields into standard order
    blockType = arrayfun(@(x) x.blockType,data(ii).stimStruct(:),'uniformoutput',false);
    testCont  = arrayfun(@(x) x.testCont,data(ii).stimStruct(:));
    testDist  = arrayfun(@(x) x.testDist,data(ii).stimStruct(:));
    refVel    = arrayfun(@(x) x.refVel,data(ii).stimStruct(:));
    
    blockID = nan(numel(blockType),1);
    
    for jj = 1:numel(blockType)
        
        if strcmp(blockType{jj},'within')
            blockID(jj) = 1;
        elseif strcmp(blockType{jj},'between')
            blockID(jj) = 2;
        elseif strcmp(blockType{jj},'sizeCntl')
            blockID(jj) = 3;
        end
        
    end
    
    theseStimInds = [blockID testCont testDist refVel];
    [~,sortInds]  = sortrows(theseStimInds);

    data(ii).stimStruct = data(ii).stimStruct(sortInds);
    
end

% Sort array of stim parameters
testCont = testCont(sortInds);
testDist = testDist(sortInds);
refVel   = refVel(sortInds);
blockID  = blockID(sortInds);

numSubs = numel(data);

% Make pairings
highVel = max(refVel);
lowVel  = min(refVel);

exp1_lowContLoVel = (refVel == lowVel)  & (testCont == 0.1) & (blockID == 1);
exp1_lowContHiVel = (refVel == highVel) & (testCont == 0.1) & (blockID == 1);
exp1_HiContLoVel  = (refVel == lowVel)  & (testCont == 0.5) & (blockID == 1);
exp1_HiContHiVel  = (refVel == highVel) & (testCont == 0.5) & (blockID == 1);

exp2_HiContLoVel  = (refVel == lowVel)  & (testCont == 0.5) & (blockID == 2);
exp2_HiContHiVel  = (refVel == highVel) & (testCont == 0.5) & (blockID == 2);

% exp1_lowContLoVel = (refVel ~= highVel)  & (testCont == 0.1) & (blockID == 1);
% exp1_lowContHiVel = (refVel == highVel) & (testCont == 0.1) & (blockID == 1);
% exp1_HiContLoVel  = (refVel ~= highVel)  & (testCont == 0.5) & (blockID == 1);
% exp1_HiContHiVel  = (refVel == highVel) & (testCont == 0.5) & (blockID == 1);
% 
% exp2_HiContLoVel  = (refVel ~= highVel)  & (testCont == 0.5) & (blockID == 2);
% exp2_HiContHiVel  = (refVel == highVel) & (testCont == 0.5) & (blockID == 2);

exp3_lowContHiVel = (refVel == highVel) & (testCont == 0.1) & (blockID == 3 | blockID == 1) & (testDist == 0.5);
exp3_HiContHiVel  = (refVel == highVel) & (testCont == 0.5) & (blockID == 3 | blockID == 1) & (testDist == 0.5);

% Initial data extraction/sorting
for ii = 1:numSubs

    % Select dataset
    thisDataset = data(ii).stimStruct(:);
    
    % Grab PSEs
    PSEs(:,ii) = arrayfun(@(x) x.PSE,thisDataset);
    PSECILB(:,ii) = arrayfun(@(x) x.PSECI(1),thisDataset);
    PSECIUB(:,ii) = arrayfun(@(x) x.PSECI(2),thisDataset);
    
    % Grab slopes
    slopes(:,ii) = arrayfun(@(x) x.slope,thisDataset);
    slopeCILB(:,ii) = arrayfun(@(x) x.slopeCI(1),thisDataset);
    slopeCIUB(:,ii) = arrayfun(@(x) x.slopeCI(2),thisDataset);
    
end 


%% Run intersubject statistics



%% Plot

spread   = linspace(-0.1,0.1,numSubs);
colorMat = colororder;
ytixLin  = [0.001 0.1 0.5 1 2 4 8 16 32];
ytix     = getLogXform(ytixLin,0.3);
for ii = 1:numel(ytixLin)
    ytixlabel{ii} = num2str(ytixLin(ii));
end

% yLims = [getLogXform([0.1 24],0.3);...
%          getLogXform([0.5 36],0.3)];
yLims = [getLogXform([0.1 24],0.3);...
         getLogXform([2 36],0.3)];
     
% Comparison inds
if numel(unique(blockID)) == 1
    if blockID(1) == 1
        compInds   = [exp1_lowContLoVel exp1_lowContHiVel exp1_HiContLoVel exp1_HiContHiVel];
        compFiles  = {'exp1_HiContLoVel','exp1_HiContHiVel','exp1_LoContLoVel','exp1_LoContHiVel'};
        compLabels = {'Within: C_{Ref}=0.5, V_{Ref}=1\circ/s','Within: C_{Ref}=0.5, V_{Ref}=8\circ/s',...
                      'Within: C_{Ref}=0.1, V_{Ref}=1\circ/s','Within: C_{Ref}=0.1, V_{Ref}=8\circ/s'};
    else
        compInds   = [exp2_HiContLoVel exp2_HiContHiVel];
        compFiles  = {'exp2_HiContLoVel','exp2_HiContHiVel'};
        compLabels = {'Between: C_{Ref}=0.5, V_{Ref}=1\circ/s','Between: C_{Ref}=0.5, V_{Ref}=8\circ/s'};
    end
elseif numel(unique(blockID)) == 2
    compInds   = [exp1_lowContLoVel exp1_lowContHiVel exp1_HiContLoVel exp1_HiContHiVel ...
                  exp2_HiContLoVel exp2_HiContHiVel];
    compFiles  = {'exp1_HiContLoVel','exp1_HiContHiVel','exp1_LoContLoVel','exp1_LoContHiVel',...
                  'exp2_HiContLoVel','exp2_HiContHiVel'};
    compLabels = {'Within: C_{Ref}=0.5, V_{Ref}=1\circ/s','Within: C_{Ref}=0.5, V_{Ref}=8\circ/s',...
                  'Within: C_{Ref}=0.1, V_{Ref}=1\circ/s','Within: C_{Ref}=0.1, V_{Ref}=8\circ/s',...
                  'Between: C_{Ref}=0.5, V_{Ref}=1\circ/s','Between: C_{Ref}=0.5, V_{Ref}=8\circ/s'};
              
elseif numel(unique(blockID)) == 3
    compInds   = [exp1_lowContLoVel exp1_lowContHiVel exp1_HiContLoVel exp1_HiContHiVel ...
                  exp2_HiContLoVel exp2_HiContHiVel exp3_lowContHiVel exp3_HiContHiVel];
    compFiles  = {'exp1_HiContLoVel','exp1_HiContHiVel','exp1_LoContLoVel','exp1_LoContHiVel',...
                  'exp2_HiContLoVel','exp2_HiContHiVel','exp3_HiContHiVel','exp3_LoContHiVel'};
    compLabels = {'Within: C_{Ref}=0.5, V_{Ref}=1\circ/s','Within: C_{Ref}=0.5, V_{Ref}=8\circ/s',...
                  'Within: C_{Ref}=0.1, V_{Ref}=1\circ/s','Within: C_{Ref}=0.1, V_{Ref}=8\circ/s',...
                  'Between: C_{Ref}=0.5, V_{Ref}=1\circ/s','Between: C_{Ref}=0.5, V_{Ref}=8\circ/s',...
                  'Within (SC): C_{Ref}=0.5, V_{Ref}=8\circ/s','Within (SC): C_{Ref}=0.1, V_{Ref}=8\circ/s'};    
end

% Cull empty comparisons
compInds   = compInds(:,sum(compInds) ~= 0);
compFiles  = compFiles(sum(compInds) ~= 0);
compLabels = compLabels(sum(compInds) ~= 0);

legEnts = [subjID,'Group mean'];

for hh = 1:size(compInds,2)
    
    theseDists  = testDist(compInds(:,hh));
    theseBlocks = blockID(compInds(:,hh));
    thisRefVel  = refVel(compInds(:,hh));
    theseConts  = testCont(compInds(:,hh));
    
    % Flip to put in terms of reference distance for Exp 2
    xlab = 'Viewing Distance (m)';
    
    if theseBlocks(1) == 2
        theseDists = flipud(theseDists);
        xlab = 'Reference Distance (m)';
    end
    
    if theseBlocks(2) ~= 3
    % PSE xsubs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f1(hh) = figure;
    f1(hh).Position = [100 100 650 650];
    hold on;
    
    plot([0 3],getLogXform(thisRefVel(1),0.3)*[1 1],'--k')
    
    for ii = 1:numSubs
        
%         thisColor  = colorMat(ii,:);
        thisColor  = 0.75*[1 1 1];
        
        thesePSEs   = PSEs(compInds(:,hh),ii);
        thesePSELBs = PSECILB(compInds(:,hh),ii);
        thesePSEUBs = PSECIUB(compInds(:,hh),ii);
        
        errorbar(theseDists + spread(ii),thesePSEs,thesePSEs-thesePSELBs,thesePSEUBs-thesePSEs,'color',thisColor,'linewidth',2);
      
        if (theseConts(1) == 0.5) && (theseBlocks(1) == 1)
            p1(ii) = scatter(theseDists + spread(ii),thesePSEs,150,thisColor,'linewidth',2);
        else
            p1(ii) = scatter(theseDists + spread(ii),thesePSEs,150,thisColor,'filled');
        end
        
    end
    
    % Plot group mean and STD
    thisPSEData   = PSEs(compInds(:,hh),:);
    
    meanPSEs               = mean(thisPSEData,2);
    stdPSEs                = std(thisPSEData,[],2);
    [~,~,groupPSECI(1,:)]  = ttest(thisPSEData(1,:));
    [~,~,groupPSECI(2,:)]  = ttest(thisPSEData(2,:));
    
    errorbar(theseDists,meanPSEs,meanPSEs-groupPSECI(:,1),groupPSECI(:,2)-meanPSEs,'color',[0 0 0],'linewidth',2);
    if (theseConts(1) == 0.5) && (theseBlocks(1) == 1)
        p1b = scatter(theseDists,meanPSEs,200,[0 0 0],'linewidth',2);
    else
        p1b = scatter(theseDists,meanPSEs,200,[0 0 0],'filled');
    end
    
    pSRPSE(hh)        = signrank(thisPSEData(1,:),thisPSEData(2,:));
    [~,pNormPSE(hh)]  = ttest(thisPSEData(1,:),thisPSEData(2,:));
    
    % Set y-axis limits
    if thisRefVel(1) < 6
        theseyLims = yLims(1,:);
    else
        theseyLims = yLims(2,:);
    end
    
    set(gca,'fontsize',20,'xlim',[0 1.5],'xtick',[0.5 1],'ytick',ytix,'yticklabel',ytixlabel,'ylim',theseyLims);
    xlabel(xlab);
    ylabel('PSE (\circ/s)');
    title(compLabels{hh});
    text(0.25,0.2*abs(diff(theseyLims))+theseyLims(1),['p = ',num2str(round(pNormPSE(hh),5))],'fontsize',20);
    
%     if hh == 1
%        legend([p1(1) p1b],legEnts,'location','northeast'); 
%     end
    
    % Slope xsubs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f2(hh) = figure;
    f2(hh).Position = [800 100 650 650];
    hold on;
    
    for ii = 1:numSubs
        
%         thisColor  = colorMat(ii,:);
        thisColor  = 0.75*[1 1 1];
        
        theseSlopes   = slopes(compInds(:,hh),ii);
        theseSlopeLBs = slopeCILB(compInds(:,hh),ii);
        theseSlopeUBs = slopeCIUB(compInds(:,hh),ii);
        
        errorbar(theseDists + spread(ii),theseSlopes,theseSlopes-theseSlopeLBs,theseSlopeUBs-theseSlopes,'color',thisColor,'linewidth',2);
        
        if (testCont(hh) == 0.5) && (blockID(hh) == 1)
            p2(ii) = scatter(theseDists + spread(ii),theseSlopes,150,thisColor);
        else
            p2(ii) = scatter(theseDists + spread(ii),theseSlopes,150,thisColor,'filled');
        end
        
    end
    
    % Plot group mean and STD
    thisSlopeData = slopes(compInds(:,hh),:);
    
    meanSlopes    = mean(thisSlopeData,2);
    stdSlopes     = std(thisSlopeData,[],2);
    [~,~,groupSlopesCI(1,:)]  = ttest(thisSlopeData(1,:));
    [~,~,groupSlopesCI(2,:)]  = ttest(thisSlopeData(2,:));

    slopedata{hh} = thisSlopeData;
    slopedataLB{hh} = slopeCILB(compInds(:,hh),:);
    slopedataUB{hh} = slopeCIUB(compInds(:,hh),:);

%     groupSlopesCI = std(slopes(compInds(:,hh),:),0,2)*0.95/2/sqrt(numSubs);
    errorbar(theseDists,meanSlopes,meanSlopes-groupSlopesCI(:,1),groupSlopesCI(:,2)-meanSlopes,'color',[0 0 0],'linewidth',2);
    p2b = scatter(theseDists,meanSlopes,200,[0 0 0],'filled');
    
    pSRSlope(hh)       = signrank(thisSlopeData(1,:),thisSlopeData(2,:));
    [~,pNormSlope(hh)] = ttest(thisSlopeData(1,:),thisSlopeData(2,:));
    
    set(gca,'fontsize',20,'xlim',[0 1.5],'xtick',[0.5 1]);
    xlabel('Test Distance (m)');
    ylabel('Slope (\Delta p/\Delta log(V))');
    title(compLabels{hh});
    text(0.25,2,['p = ',num2str(round(pNormSlope(hh),5))],'fontsize',20);
    
%     if hh == 1
%        legend([p2 p2b],legEnts,'location','northeast'); 
%     end
    else
        f1(hh) = figure;
        f1(hh).Position = [100 100 650 650];
        hold on;
        
        plot([0 3],getLogXform(thisRefVel(1),0.3)*[1 1],'--k')
        theseDists = [0.5 1];
        
        for ii = 1:numSubs
            
            thisColor  = 0.75*[1 1 1];
            
            thesePSEs   = PSEs(compInds(:,hh),ii);
            thesePSELBs = PSECILB(compInds(:,hh),ii);
            thesePSEUBs = PSECIUB(compInds(:,hh),ii);
            
            errorbar(theseDists + spread(ii),thesePSEs,thesePSEs-thesePSELBs,thesePSEUBs-thesePSEs,'color',thisColor,'linewidth',2);
            
            if (theseConts(1) == 0.5) && (theseBlocks(1) == 1)
                p1(ii) = scatter(theseDists + spread(ii),thesePSEs,150,thisColor,'linewidth',2);
            else
                p1(ii) = scatter(theseDists + spread(ii),thesePSEs,150,thisColor,'filled');
            end
            
        end
        
        % Plot group mean and STD
        thisPSEData   = PSEs(compInds(:,hh),:);
        
        meanPSEs               = mean(thisPSEData,2);
        stdPSEs                = std(thisPSEData,[],2);
        [~,~,groupPSECI(1,:)]  = ttest(thisPSEData(1,:));
        [~,~,groupPSECI(2,:)]  = ttest(thisPSEData(2,:));
        
        errorbar(theseDists,meanPSEs,meanPSEs-groupPSECI(:,1),groupPSECI(:,2)-meanPSEs,'color',[0 0 0],'linewidth',2);
        if (theseConts(1) == 0.5) && (theseBlocks(1) == 1)
            p1b = scatter(theseDists,meanPSEs,200,[0 0 0],'linewidth',2);
        else
            p1b = scatter(theseDists,meanPSEs,200,[0 0 0],'filled');
        end
        
        pSRPSE(hh)        = signrank(thisPSEData(1,:),thisPSEData(2,:));
        [~,pNormPSE(hh)]  = ttest(thisPSEData(1,:),thisPSEData(2,:));
        
        % Set y-axis limits
        if thisRefVel(1) < 6
            theseyLims = yLims(1,:);
        else
            theseyLims = yLims(2,:);
        end
        
        set(gca,'fontsize',20,'xlim',[0 1.5],'xtick',[0.5 1],'ytick',ytix,'yticklabel',ytixlabel,...
            'ylim',theseyLims,'xticklabel',{'large','small'});
        xlabel('Stimulus size');
        ylabel('PSE (\circ/s)');
        title(compLabels{hh});
        text(0.25,0.2*abs(diff(theseyLims))+theseyLims(1),['p = ',num2str(round(pNormPSE(hh),5))],'fontsize',20);
        
    end

    meanPSEsAll{hh}   = meanPSEs;
    stdPSEsAll{hh}   = stdPSEs;
    groupPSECIAll{hh} = groupPSECI;

    if theseBlocks(2) ~= 3
        meanSlopesAll{hh} = meanSlopes;
        stdSlopesAll{hh} = stdSlopes;
        groupSlopesCIAll{hh} = groupSlopesCI;
    end
    
end

% Since we were lazy in assigning cell array inds, cull the empty ones
inc1 = cellfun(@(x) ~isempty(x),meanPSEsAll);
inc2 = cellfun(@(x) ~isempty(x),groupPSECIAll);
inc3 = cellfun(@(x) ~isempty(x),meanSlopesAll);
inc4 = cellfun(@(x) ~isempty(x),groupSlopesCIAll);
inc5 = cellfun(@(x) ~isempty(x),stdPSEsAll);
inc6 = cellfun(@(x) ~isempty(x),stdSlopesAll);

% Output stats
stats.pSRSlope      = pSRSlope;
stats.pNormSlope    = pNormSlope;
stats.pSRPSE        = pSRPSE;
stats.pNormPSE      = pNormPSE;
stats.slopeData     = slopedata;
stats.slopeDataLB   = slopedataLB;
stats.slopeDataUB   = slopedataUB;
stats.PSEs          = PSEs;
stats.PSECILB       = PSECILB;
stats.PSECIUB       = PSECIUB;
stats.meanPSEs      = meanPSEsAll(inc1);
stats.groupPSECI    = groupPSECIAll(inc2);
stats.meanSlopes    = meanSlopesAll(inc3);
stats.groupSlopesCI = groupSlopesCIAll(inc4);
stats.stdPSEs       = stdPSEsAll(inc5);
stats.stdSlopes     = stdSlopesAll(inc6);


%% Save plots

if saveOn
   
    if nargin>2
        subdir = [varargin{2},'/'];
    else
        subdir = '';
    end
    
    if ~exist([dataDir,'figures/',subdir])
        mkdir([dataDir,'figures/',subdir]);
    end
    
    for ii = 1:numel(p1)   
        saveas(f1(ii),[dataDir,'figures/',subdir,compFiles{ii},'_PSE.svg']);
%         saveas(f1(ii),[dataDir,'figures/',compFiles{ii},'_PSE.png']);
    end
    
    for ii = 1:numel(p2)
        saveas(f2(ii),[dataDir,'figures/',subdir,compFiles{ii},'_slope.svg']);
%         saveas(f2(ii),[dataDir,'figures/',compFiles{ii},'_slope.png']);
    end
    
end

end