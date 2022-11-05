
% Define what to do

singleSub = 0;
xSub      = 1;
plotStair = 0;


%% Single Subject
if singleSub

    % List of subjectIDs to include
    perspective = 1;
    
    if perspective
        % Perspective cue cohort
        subjIDs = {'PXD','KGJ','GGK','EPI','ICK','OYT'};    %%% CULLED PSE OUT OF RANGE
        
        %     subjIDs = {'PXD','KGJ','GGK','EPI','ICK','JQN','VXM','ZQC','DRH','HOZ'};  %% w/SOME SPEED NOT BLOCKED
    else
        % Matched retinal size cohort
        %     subjIDs = {'NYY','XLR','GDZ','TRQ','ZCJ','HEV','BCG','KKB','WVU','JSO'};
        
        subjIDs = {'TRQ','ZCJ','HEV','BCG','KKB','WVU'};
    end

    numSubs = numel(subjIDs);
    saveOn  = 1;
    
    for ii = 1:numSubs
        
        [data] = analyzeExpRDS_noPlot(subjIDs{ii},saveOn);
        
    end
end


%% Cross subject
if xSub
    
%     close all
    saveOn  = 0;
    
    % List of subjectIDs to include
    perspective = 1;
    
    if perspective
        % Perspective cue cohort
        subjIDs = {'PXD','KGJ','GGK','EPI','ICK','OYT'};    %%% CULLED PSE OUT OF RANGE
        
%         subjIDs = {'PXD','KGJ','GGK','EPI','ICK','JQN','VXM','ZQC','DRH','HOZ'};  %% w/SOME SPEED NOT BLOCKED
        saveDir = 'perspectiveCue';
    else
        % Matched retinal size cohort
        %     subjIDs = {'NYY','XLR','GDZ','TRQ','ZCJ','HEV','BCG','KKB','WVU','JSO'};
        
        subjIDs = {'TRQ','ZCJ','HEV','BCG','KKB','WVU'};

%         subjIDs = {'TSM'};
        saveDir = 'matchedRetinalSize';
    end
    
%     analyzeXsubRDS(saveOn);
    stats = analyzeXsubRDS(saveOn,subjIDs,saveDir);
    
end


%% Plot a staircase triplet/psychometric curve from a subject

if plotStair

%     close all
    
    sub   = 'TSM';
%     sub   = 'PXD';
    exp   = 'within';
%     exp   = 'between';
%     exp   = 'sizeCntl';
    cTest = [0.5];
    vRef  = [8];
%     vRef  = [2;8];
    dTest = [0.5];
    
    stimInds = makeCombos([numel(cTest) numel(dTest) numel(vRef)]);
    stims    = [cTest(stimInds(:,1)) dTest(stimInds(:,2)) vRef(stimInds(:,3))];
    
    for ii = 1:size(stims,1)
        [f1] = plotThisData(sub,exp,stims(ii,:));
    end
end