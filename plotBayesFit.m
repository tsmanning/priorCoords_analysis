function [SGnllScaled,MoGnllScaled,allSigsSG,allSigs,priHatSG,priHat,f1,f2,f3,f4] = plotBayesFit(subjID,fitMatSG,fitMat,stimStruct,saveDir,geoOn,saveOn)

% Fix this later on when doing multiple subs
numSubs = 1;

%% Extract parameters
xgrid = fitMat.xgrid;
xtix = getLogXform([5*1e-3 0.5 1 2 4 8 16],0.3);

xgridLin = getLinXform(xtix,0.3);
for ii = 1:numel(xtix)
    xLinLab{ii} = num2str(xgridLin(ii));
end

priHatSG = fitMatSG.priHat;
priHat   = fitMat.priHat;

contSigsSG = fitMatSG.pSigNseHat';
contSigs = fitMat.pSigNseHat';

if geoOn
    distSigsSG = unique(sqrt(fitMatSG.stimD_SF));
    distSigs = unique(sqrt(fitMat.stimD_SF));
else
    distSigsSG = unique(fitMatSG.pSigDistHat);
    distSigs = unique(fitMat.pSigDistHat);
end

% Conts: columns, Dists: rows
allSigsSG  = repmat(contSigsSG,[numel(distSigsSG) 1]) .* repmat(distSigsSG,[1 numel(contSigsSG)]);
allSigs  = repmat(contSigs,[numel(distSigs) 1]) .* repmat(distSigs,[1 numel(contSigs)]);

if iscolumn(fitMat.contVals)
    fitMat.contVals = fitMat.contVals';
end

contVals = repmat(fitMat.contVals,[numel(distSigsSG) 1]);
distVals = repmat(fitMat.distVals,[1 numel(contSigsSG)]);

contVals = contVals(:);
distVals = distVals(:);

if numel(distVals) < numel(distSigs)
    distVals = [0.25 0.5 0.5 0.75 1]';
end

allSigsSG  = allSigsSG(:);
allSigs  = allSigs(:);

% Calculate scaled goodness of fit
pCumGauss = -sum(arrayfun(@(x) x.nllCG,stimStruct));
pCoinFlip = -sum(arrayfun(@(x) x.nllCF,stimStruct));

SGnllScaled  = (-fitMatSG.nll - pCoinFlip) / (pCumGauss - pCoinFlip);
MoGnllScaled = (-fitMat.nll - pCoinFlip) / (pCumGauss - pCoinFlip);

% Define title based on whether or not geometric noise propagation model
% was used
if geoOn
    thisTitle = 'Geometric noise model';
else
    thisTitle = 'Nonparametric';
end


%% Plot priors
f1 = figure;
f1.Position = [100 100 650 600];
hold on;

plot(xgrid,priHat,'b','linewidth',2);
plot(xgrid,priHatSG,'g','linewidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',xtix,'xticklabel',xLinLab,'yscale','log',...
        'xlim',[xtix(1) xtix(end)],'fontsize',20,'ylim',[1e-20 1]);
xlabel('velocity (deg/s)');
ylabel('probability');
title(['Speed prior(',thisTitle,')']);
legend('MoG','SG');

f2 = figure;
f2.Position = [600 100 650 600];
hold on;

plot(xgrid,priHat,'b','linewidth',2);
plot(xgrid,priHatSG,'g','linewidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',xtix,'xticklabel',xLinLab,...
        'xlim',[xtix(1) xtix(end)],'fontsize',20,'ylim',[0 1]);
xlabel('velocity (deg/s)');
ylabel('probability');
title(['Speed prior(',thisTitle,')'])
legend('MoG','SG');


%% Plot likelihood widths
f3 = figure;
f3.Position = [100 600 650 600];
hold on;

for ii = 1:numel(distVals)
    if ii ~= 1
        if distVals(ii) == distVals(ii-1)
            suff = '_alt';
        else
            suff = '';
        end
    else
        suff = '';
    end
    stimLabs{ii} = ['C=',num2str(contVals(ii)),', D=',num2str(distVals(ii)),suff];
end

X = categorical(stimLabs);
X = reordercats(X,stimLabs);

Y = [allSigsSG';allSigs'];

b = bar(X,Y,'FaceColor','flat');
colors = [0 1 0; 0 0 1];

for ii = 1:2
    b(ii).CData = colors(ii,:);
end
ylabel('Likelihood \sigma');
xtickangle(45);
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',15);
legend('SG','MoG');
title(thisTitle);


%% Plot model comparison
f4 = figure;
f4.Position = [600 600 800 600];
hold on;

X = 1;
Y = [SGnllScaled;MoGnllScaled];

b2 = bar(X,Y,'FaceColor','flat');
for ii = 1:2
    b2(ii).CData = colors(ii,:);
end

plot([0 numSubs],pCumGauss*[1 1],'--k');
plot([0 numSubs],pCoinFlip*[1 1],'--k');

ylabel('Log-probability of data');
xlabel('Subject');
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',15,'xtick',1,'ylim',[0 1],'ytick',[0 1],'yticklabel',{'Coin flip','Cum. Gauss.'});
legend('SG','MoG');
title(thisTitle);

%% Plot estimated distance comparison
f5 = figure;
f5.Position = [200 200 800 600];
hold on;

distsSG = unique(cell2mat(arrayfun(@(x) x.estDistSG,stimStruct,'UniformOutput',false)));
distsMoG = unique(cell2mat(arrayfun(@(x) x.estDistMoG,stimStruct,'UniformOutput',false)));

scatter(unique(distVals),distsSG,100,'g','filled');
scatter(unique(distVals),distsMoG,100,'b','filled');

plot([0 2],[0 2],'--k');

set(gca,'plotboxaspectratio',[1 1 1],'fontsize',15,'xlim',[0 2],'ylim',[0 2]);
xlabel('True distance');
ylabel('Perceived distance');
legend('SG','MoG');
title(thisTitle);


%% Save figures
if saveOn
    saveas(f1,[saveDir,'/figures/',subjID,'geo',num2str(geoOn),'_priorFitsLog.svg']);
    saveas(f2,[saveDir,'/figures/',subjID,'geo',num2str(geoOn),'_priorFitsLin.svg']);
    saveas(f3,[saveDir,'/figures/',subjID,'geo',num2str(geoOn),'_likelihoodFits.svg']);
    saveas(f4,[saveDir,'/figures/',subjID,'geo',num2str(geoOn),'_modelComparison.svg']);
    saveas(f5,[saveDir,'/figures/',subjID,'geo',num2str(geoOn),'_distanceEst.svg']);
end

end