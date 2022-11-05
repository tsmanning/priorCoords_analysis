function [dataStruct,f1,f2] = plotBayesPsychFxns(subjID,geoOn,dataStruct,saveDir,saveOn)

% Takes data from psychophysics exp, fits with cumulative Gaussian and also
% plots Bayesian fit for comparison

%% Plot psychometric functions (cum Gauss, SG, MoG)

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

    if saveOn
        % Save figures
        saveas(f1,[saveDir,'figures/',subjID,'_geo',num2str(geoOn),'_betweenTrialsFits.svg']);
    end
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

    if saveOn
        % Save figures
        saveas(f2,[saveDir,'figures/',subjID,'geo',num2str(geoOn),'_withinTrialsFits.svg']);
    end
end


end