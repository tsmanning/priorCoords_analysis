function [f1,f2,f3] = plotAllPsychFxns(dataDir,subjID,saveOn)

% Plot all data from experimental sessions and best fitting cumulative
% Gaussian functions

load([dataDir,subjID,'.mat']);

a = stimStruct(:);
inds = arrayfun(@(x)~isempty(x.blockType),a);
a = a(inds);
[types,b] = sort(arrayfun(@(x) x.blockType,a,'UniformOutput',false));

betweenInds  = find(strcmp('between',types)==1);
withinInds   = find(strcmp('within',types)==1);
sizeCntlInds = find(strcmp('sizeCntl',types)==1);

a = a(b);

vels = [0.5 1 2 4 8 16 32];
for ii = 1:numel(vels)
    velLabs{ii} = num2str(vels(ii));
end

% ID bad fits/psychometric functions
badFits = arrayfun(@(x) x.badFit,a);


%% Plot
% Plot Between trials
if ~isempty(betweenInds)
    f1 = figure;
    f1.Position = [100 100 1075 450];
    n = 1;
    for ii = betweenInds'

        subplot(2,4,n);
        n = n+1;
        hold on;

        linVels = getLinXform(a(ii).respData(:,1),0.3);
        linRef = a(ii).refVel;

        if badFits(ii)
            plot(linVels,a(ii).pFxn,'r','linewidth',2);
            scatter(linVels,a(ii).respData(:,2)./a(ii).respData(:,3),10*a(ii).respData(:,3),'r');
        else
            plot(linVels,a(ii).pFxn,'k','linewidth',2);
            scatter(linVels,a(ii).respData(:,2)./a(ii).respData(:,3),10*a(ii).respData(:,3),'k');
        end
            
        plot(linRef*[1 1],[0 1],'--k');
        set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0.1 30],'ylim',[0 1],'xscale','log','xtick',vels,'xticklabel',velLabs);
        title([a(ii).blockType,', D_{ref}=',num2str(a(ii).testDist),', D_{test}=',num2str(a(ii).testDist)]);

    end
else
    f1 = [];
end

% Plot Size control trials
if ~isempty(sizeCntlInds)
    f2 = figure;
    f2.Position = [100 100 300 300];
    for ii = sizeCntlInds'

        hold on;

        linVels = getLinXform(a(ii).respData(:,1),0.3);
        linRef = a(ii).refVel;

        if badFits(ii)
            plot(linVels,a(ii).pFxn,'r','linewidth',2);
            scatter(linVels,a(ii).respData(:,2)./a(ii).respData(:,3),10*a(ii).respData(:,3),'r');
        else
            plot(linVels,a(ii).pFxn,'k','linewidth',2);
            scatter(linVels,a(ii).respData(:,2)./a(ii).respData(:,3),10*a(ii).respData(:,3),'k');
        end

        plot(linRef*[1 1],[0 1],'--k');
        set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0.1 30],'ylim',[0 1],'xscale','log','xtick',vels,'xticklabel',velLabs);
        title('Size control (Diam. Ref > test)');

    end
else
    f2 = [];
end

% Plot Within trials
if ~isempty(withinInds)
    f3 = figure;
    f3.Position = [100 100 1400 930];
    n = 1;
    for ii = withinInds'

        subplot(4,5,n);
        n = n+1;
        hold on;

        linVels = getLinXform(a(ii).respData(:,1),0.3);
        linRef = a(ii).refVel;

        if badFits(ii)
            plot(linVels,a(ii).pFxn,'r','linewidth',2);
            scatter(linVels,a(ii).respData(:,2)./a(ii).respData(:,3),10*a(ii).respData(:,3),'r');
        else
            plot(linVels,a(ii).pFxn,'k','linewidth',2);
            scatter(linVels,a(ii).respData(:,2)./a(ii).respData(:,3),10*a(ii).respData(:,3),'k');
        end

        plot(linRef*[1 1],[0 1],'--k');
        set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0.1 30],'ylim',[0 1],'xscale','log','xtick',vels,'xticklabel',velLabs);
        title([a(ii).blockType,', C_{test}=',num2str(a(ii).testCont),', D=',num2str(a(ii).testDist)]);

    end
else
    f3 = [];
end

%% Save figures
if saveOn
    
    subjDir = [dataDir,'/figs/',subjID,'/'];

    if ~exist(subjDir)
        mkdir(subjDir);
    end

    if ~isempty(f1)
        saveas(f1,[subjDir,'betweenTrials.svg']);
    end
    if ~isempty(f2)
        saveas(f2,[subjDir,'sizeContTrials.svg']);
    end
    if ~isempty(f3)
        saveas(f3,[subjDir,'withinTrials.svg']);
    end
    
end

end