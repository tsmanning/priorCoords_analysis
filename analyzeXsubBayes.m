function [] = analyzeXsubBayes(saveOn,subjIDs,saveDir)

% Aggregate Bayesian ideal observer model fits from all subjects and plot
% in single plots

cullBF = 1;

%% Load in all subjects' data/fits

nonpar_mogFits = [];
nonpar_sgFits  = [];
nonpar_datastr = [];
allSortedPars  = [];
geo_mogFits    = [];
geo_sgFits     = [];
geo_datastr    = [];

numSubs = numel(subjIDs);

for ii = 1:numSubs

    load([saveDir,subjIDs{ii},'geoOn0']);

    nonpar_mogFits = [nonpar_mogFits; outMatMoG];
    nonpar_sgFits  = [nonpar_sgFits; outMatSG];
    nonpar_datastr{ii,1} = dataStruct;
    allSortedPars  = [allSortedPars; sortedPars];

    load([saveDir,subjIDs{ii},'geoOn1']);

    geo_mogFits = [geo_mogFits; outMatMoG];
    geo_sgFits  = [geo_sgFits; outMatSG];
    geo_datastr{ii,1} = dataStruct;

end

xgrid = outMatSG.xgrid;



%% Aggregate fits of interest

% Priors
for jj = 1:numSubs

    thisPrior(1,jj,:) = geo_mogFits(jj).priHat;
    thisPrior(2,jj,:) = nonpar_mogFits(jj).priHat;

end


% Estimated Distances
dists = [0.25 0.5 0.75 1];
estDistsGeo = nan(numSubs,numel(dists));
estDistsNP  = nan(numSubs,numel(dists));

for ii = 1:numSubs

    theseDists{ii}= unique(geo_mogFits(ii).distVals);
    theseDistsGeo = unique(cell2mat(arrayfun(@(x) x.estDistMoG,geo_datastr{ii},'UniformOutput',false)));
    theseDistsNP  = unique(cell2mat(arrayfun(@(x) x.estDistMoG,nonpar_datastr{ii},'UniformOutput',false)));

    for jj = 1:numel(theseDists{ii})
    
        % Not everyone has a complete distance dataset, gotta do this goofy
        % matching thing
        thisDistInd = find(theseDists{ii}(jj) == dists);

        estDistsGeo(ii,thisDistInd) = theseDistsGeo(jj);
        estDistsNP(ii,thisDistInd)  = theseDistsNP(jj);

    end

end

numDists = numel(dists);

% Likelihood widths
for ii = 1:numSubs

    % Get contrast components
    geo_contSigsSG     = geo_sgFits(ii).pSigNseHat';
    geo_contSigsMoG    = geo_mogFits(ii).pSigNseHat';
    nonPar_contSigsSG  = nonpar_sgFits(ii).pSigNseHat';
    nonPar_contSigsMoG = nonpar_mogFits(ii).pSigNseHat';

    numConts = numel(geo_contSigsSG);

    % Get distance components
    geo_distSigsSG     = unique(sqrt(geo_sgFits(ii).stimD_SF));
    geo_distSigsMoG    = unique(sqrt(geo_mogFits(ii).stimD_SF));
    nonPar_distSigsSG  = nonpar_sgFits(ii).pSigDistHat;
    nonPar_distSigsMoG = nonpar_mogFits(ii).pSigDistHat;

    theseDists{ii}     = unique(geo_mogFits(ii).distVals);
    sg_geoDists        = nan(1,numDists);
    mog_geoDists       = nan(1,numDists);
    sg_npDists         = nan(1,numDists);
    mog_npDists        = nan(1,numDists);

    for jj = 1:numel(theseDists{ii})
    
        % Not everyone has a complete distance dataset, gotta do this goofy
        % matching thing
        thisDistInd = find(theseDists{ii}(jj) == dists);

        sg_geoDists(thisDistInd)  = geo_distSigsSG(jj);
        mog_geoDists(thisDistInd) = geo_distSigsMoG(jj);
        sg_npDists(thisDistInd)   = nonPar_distSigsSG(jj);
        mog_npDists(thisDistInd)  = nonPar_distSigsMoG(jj);

    end 

    % Expand out everything to 1x(numContsxnumDists) vectors
    geo_contSigsSG     = repelem(geo_contSigsSG,numDists);
    geo_distSigsSG     = repmat(sg_geoDists,[1 numConts]);
    geo_contSigsMoG    = repelem(geo_contSigsMoG,numDists);
    geo_distSigsMoG    = repmat(mog_geoDists,[1 numConts]);
    nonPar_contSigsSG  = repelem(nonPar_contSigsSG,numDists);
    nonPar_distSigsSG  = repmat(sg_npDists,[1 numConts]); 
    nonPar_contSigsMoG = repelem(nonPar_contSigsMoG,numDists);
    nonPar_distSigsMoG = repmat(mog_npDists,[1 numConts]);

    % Make matrix of possible likelihood widths (2xConts, 4xDists, SG v. MoG, Geo v. NonPar)
    geoSGSigs(ii,:)     = geo_contSigsSG.*geo_distSigsSG;
    geoMoGSigs(ii,:)    = geo_contSigsMoG.*geo_distSigsMoG;
    nonParSGSigs(ii,:)  = nonPar_contSigsSG.*nonPar_distSigsSG;
    nonParMoGSigs(ii,:) = nonPar_contSigsMoG.*nonPar_distSigsMoG;

end


% Fit quality
for ii = 1:numSubs

    %%%% should really re-run fitPsychFxns2 and get rid of this
    theseCGNlls = arrayfun(@(x) x.nllCG,nonpar_datastr{ii,1});
    theseCFNlls = arrayfun(@(x) x.nllCF,nonpar_datastr{ii,1});

    incNlls = ~isnan(theseCGNlls) & ~isnan(theseCFNlls);

    pCumGauss = -sum(theseCGNlls(incNlls));
    pCoinFlip = -sum(theseCFNlls(incNlls));
    %%%%

    geo_SGnllScaled(ii)  = (-geo_sgFits(ii).nll - pCoinFlip) / (pCumGauss - pCoinFlip);
    geo_MoGnllScaled(ii) = (-geo_mogFits(ii).nll - pCoinFlip) / (pCumGauss - pCoinFlip);
    nonPar_SGnllScaled(ii)  = (-nonpar_sgFits(ii).nll - pCoinFlip) / (pCumGauss - pCoinFlip);
    nonPar_MoGnllScaled(ii) = (-nonpar_mogFits(ii).nll - pCoinFlip) / (pCumGauss - pCoinFlip);

    numPars(ii,:) = [size(geo_sgFits(ii).prs0,1);...
                     size(geo_mogFits(ii).prs0,1);...
                     size(nonpar_sgFits(ii).prs0,1);...
                     size(nonpar_mogFits(ii).prs0,1)];
    numDpts = size(geo_sgFits(ii).stimD_SFall,1);

    CFbic(ii)         = -2*pCoinFlip;
    CGbic(ii)         = -2*pCumGauss + sum(incNlls)*log(numDpts);

    geo_SGbic(ii)     = 2*geo_sgFits(ii).nll + numPars(ii,1)*log(numDpts);
    geo_MoGbic(ii)    = 2*geo_mogFits(ii).nll + numPars(ii,2)*log(numDpts);
    nonPar_SGbic(ii)  = 2*nonpar_sgFits(ii).nll + numPars(ii,3)*log(numDpts);
    nonPar_MoGbic(ii) = 2*nonpar_mogFits(ii).nll + numPars(ii,4)*log(numDpts);

%     geo_SGbic(ii)     = (geo_SGbic(ii) - CFbic(ii)) / (CGbic(ii) - CFbic(ii));
%     geo_MoGbic(ii)    = (geo_MoGbic(ii) - CFbic(ii)) / (CGbic(ii) - CFbic(ii));
%     nonPar_SGbic(ii)  = (nonPar_SGbic(ii) - CFbic(ii)) / (CGbic(ii) - CFbic(ii));
%     nonPar_MoGbic(ii) = (nonPar_MoGbic(ii) - CFbic(ii)) / (CGbic(ii) - CFbic(ii));

end



%% Optionally cull bad fits
if cullBF

    % ID bad fits based on non-parametric MoG NLL
    cutoff = 0.3;
    goodFits = nonPar_MoGnllScaled >= cutoff;

    % Select only subjects with full datasets
    fullData = cellfun(@(x) numel(x) == 28,nonpar_datastr)';
    
    incInds = goodFits & fullData;
    numSubs = sum(incInds);

    % Output subs that didn't make the cut
    for ii = 1:sum(~incInds)
        badSubs = subjIDs(~incInds);
        badSubsNll = nonPar_MoGnllScaled(~incInds);

        disp(['Excluded: ',badSubs{ii},': ',num2str(round(badSubsNll(ii),3))]);
    end
    % Output subs that didn't make the cut
    for ii = 1:sum(incInds)
        goodSubs = subjIDs(incInds);

        disp(['Included: ',goodSubs{ii}]);
    end

    % Cull priors
    thisPrior = thisPrior(:,incInds,:);

    % Cull likelihoods
    geoSGSigs     = geoSGSigs(incInds,:);
    geoMoGSigs    = geoMoGSigs(incInds,:);
    nonParSGSigs  = nonParSGSigs(incInds,:);
    nonParMoGSigs = nonParMoGSigs(incInds,:);

    % Cull distances
    theseDists  = theseDists(incInds);
    estDistsGeo = estDistsGeo(incInds,:);
    estDistsNP  = estDistsNP(incInds,:);

    % Cull fits
    geo_SGnllScaled     = geo_SGnllScaled(incInds);
    geo_MoGnllScaled    = geo_MoGnllScaled(incInds);
    nonPar_SGnllScaled  = nonPar_SGnllScaled(incInds);
    nonPar_MoGnllScaled = nonPar_MoGnllScaled(incInds);
    numPars             = numPars(incInds);
    geo_SGbic           = geo_SGbic(incInds);
    geo_MoGbic          = geo_MoGbic(incInds);
    nonPar_SGbic        = nonPar_SGbic(incInds);
    nonPar_MoGbic       = nonPar_MoGbic(incInds);

end

%% Medians and IQRs for scatterplots across participants

% Distances
medDistGeo = median(estDistsGeo,'omitnan');
medDistNP  = median(estDistsNP,'omitnan');

iqrDisNP   = [prctile(estDistsNP,25);prctile(estDistsNP,75)];

% Fit qualities
medSGnll_geo      = median(geo_SGnllScaled,'omitnan');
medMoGnll_geo     = median(geo_MoGnllScaled,'omitnan');
medSGnll_nonPar   = median(nonPar_SGnllScaled,'omitnan');
medMoGnll_nonPar  = median(nonPar_MoGnllScaled,'omitnan');

iqrSGnll_geo      = [prctile(geo_SGnllScaled,25);prctile(geo_SGnllScaled,75)];
iqrMoGnll_geo     = [prctile(geo_MoGnllScaled,25);prctile(geo_MoGnllScaled,75)];
iqrSGnll_nonPar   = [prctile(nonPar_SGnllScaled,25);prctile(nonPar_SGnllScaled,75)];
iqrMoGnll_nonPar  = [prctile(nonPar_MoGnllScaled,25);prctile(nonPar_MoGnllScaled,75)];

medSGbic_geo     = median(geo_SGbic,'omitnan');
medMoGbic_geo    = median(geo_MoGbic,'omitnan');
medSGbic_nonPar  = median(nonPar_SGbic,'omitnan');
medMoGbic_nonPar = median(nonPar_MoGbic,'omitnan');

iqrSGbic_geo      = [prctile(geo_SGbic,25);prctile(geo_SGbic,75)];
iqrMoGbic_geo     = [prctile(geo_MoGbic,25);prctile(geo_MoGbic,75)];
iqrSGbic_nonPar   = [prctile(nonPar_SGbic,25);prctile(nonPar_SGbic,75)];
iqrMoGbic_nonPar  = [prctile(nonPar_MoGbic,25);prctile(nonPar_MoGbic,75)];

% Likelihoods
medSGSigs_geo     = median(geoSGSigs,'omitnan');
medMoGSigs_geo    = median(geoMoGSigs,'omitnan');
medSGSigs_nonPar  = median(nonParSGSigs,'omitnan');
medMoGSigs_nonPar = median(nonParMoGSigs,'omitnan');

iqrSGSigs_geo     = [prctile(geoSGSigs,25);prctile(geoSGSigs,75)];
iqrMoGSigs_geo    = [prctile(geoMoGSigs,25);prctile(geoMoGSigs,75)];
iqrSGSigs_nonPar  = [prctile(nonParSGSigs,25);prctile(nonParSGSigs,75)];
iqrMoGSigs_nonPar = [prctile(nonParMoGSigs,25);prctile(nonParMoGSigs,75)];


%% Plot all MoG Priors

pTitles = {'Geometric','Nonparametric'};
xtix = getLogXform([5*1e-3 0.5 1 2 4 8 16],0.3);
colors = [0 1 0; 0 0 1];

xgridLin = getLinXform(xtix,0.3);
for ii = 1:numel(xtix)
    xLinLab{ii} = num2str(xgridLin(ii));
end

f1 = figure;
f1.Position = [100 100 1200 600];

for ii = 1:2
    subplot(1,2,ii);
    hold on;

    for jj = 1:numSubs
        plot(xgrid,squeeze(thisPrior(ii,jj,:)),'b','linewidth',2);
    end

    set(gca,'plotboxaspectratio',[1 1 1],'xtick',xtix,'xticklabel',xLinLab,'yscale','log',...
        'xlim',[xtix(1) xtix(end)],'fontsize',20,'ylim',[1e-10 1.5]);
    xlabel('velocity (deg/s)');
    ylabel('probability');
    title(pTitles{ii});
end


%% Plot best fit likelihood widths
offset = [-0.375 -0.125 0.125 0.375];
xVals  = linspace(1,15,8);
xLabs  = {'C=0.1,D=0.25','C=0.1,D=0.5','C=0.1,D=0.75','C=0.1,D=1.0',...
          'C=0.5,D=0.25','C=0.5,D=0.5','C=0.5,D=0.75','C=0.5,D=1.0'};

f2 = figure;
f2.Position = [700 100 1200 700];
hold on;

swarmchart(repelem(xVals+offset(1),numSubs),geoSGSigs(:),50,[0.5 1 0.5],XJitterWidth=0.1);
swarmchart(repelem(xVals+offset(2),numSubs),nonParSGSigs(:),50,[0.5 1 0.5],'filled',XJitterWidth=0.1);
swarmchart(repelem(xVals+offset(3),numSubs),geoMoGSigs(:),50,[0.5 0.5 1],XJitterWidth=0.1);
swarmchart(repelem(xVals+offset(4),numSubs),nonParMoGSigs(:),50,[0.5 0.5 1],'filled',XJitterWidth=0.1);

p(1) = scatter(xVals+offset(1),medSGSigs_geo,100,[0 0.8 0],'linewidth',2);
p(2) = scatter(xVals+offset(2),medMoGSigs_geo,100,[0 0.8 0],'filled');
p(3) = scatter(xVals+offset(3),medSGSigs_nonPar,100,[0 0 1],'linewidth',2);
p(4) = scatter(xVals+offset(4),medMoGSigs_nonPar,100,[0 0 1],'filled');
errorbar(xVals+offset(1),medSGSigs_geo,medSGSigs_geo-iqrSGSigs_geo(1,:),...
    iqrSGSigs_geo(2,:)-medSGSigs_geo,'color',[0 0.8 0],'linewidth',2,'LineStyle','none');
errorbar(xVals+offset(2),medMoGSigs_geo,medMoGSigs_geo-iqrMoGSigs_geo(1,:),...
    iqrMoGSigs_geo(2,:)-medMoGSigs_geo,'color',[0 0.8 0],'linewidth',2,'LineStyle','none');
errorbar(xVals+offset(3),medSGSigs_nonPar,medSGSigs_nonPar-iqrSGSigs_nonPar(1,:),...
    iqrSGSigs_nonPar(2,:)-medSGSigs_nonPar,'color',[0 0 1],'linewidth',2,'LineStyle','none');
errorbar(xVals+offset(4),medMoGSigs_nonPar,medMoGSigs_nonPar-iqrMoGSigs_nonPar(1,:),...
    iqrMoGSigs_nonPar(2,:)-medMoGSigs_nonPar,'color',[0 0 1],'linewidth',2,'LineStyle','none');

ylabel('Likelihood \sigma');
xtickangle(45);
set(gca,'plotboxaspectratio',[2 1 1],'fontsize',15,'xtick',xVals,'xticklabel',xLabs,...
        'xlim',[xVals(1)-0.5 xVals(end)+0.5],'yscale','log','ylim',[1e-2 1e1]);
legend(p,{'SG (geo)','SG (NonPar)','MoG (geo)','MoG (NonPar)'},'location','southeast');


%% Plot estimated distances

maxDist = 1.5;

f3 = figure;
f3.Position = [100 700 650 600];
hold on;

for ii = 1:numSubs

    theseinds = ~isnan(estDistsGeo(ii,:));

    if ii == 1
        p2(1) = scatter(theseDists{ii},estDistsGeo(ii,theseinds),100,[0.5 0.5 1],'filled');
        p2(2) = scatter(theseDists{ii},estDistsNP(ii,theseinds),100,[0.5 0.5 1],'linewidth',2);
    else
        scatter(theseDists{ii},estDistsGeo(ii,theseinds),100,[0.5 0.5 1],'filled');
        scatter(theseDists{ii},estDistsNP(ii,theseinds),100,[0.5 0.5 1],'linewidth',2);
    end
end

p2(3) = scatter(dists,medDistNP,100,'b','filled');
errorbar(dists,medDistNP,medDistNP-iqrDisNP(1,:),iqrDisNP(2,:)-medDistNP,'b','linewidth',2);

plot([0 maxDist],[0 maxDist],'--k');
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'xlim',[0 maxDist],'ylim',[0 maxDist]);
xlabel('True distance');
ylabel('Perceived distance estimate');
legend(p2,{'Geo','NonPar','Median+IQR'},'location','northwest');


%% Plot Fit quality

f4 = figure;
f4.Position = [700 700 1200 650];

subplot(1,2,1);
hold on;

xSG  = [ones(1,numSubs) 2*ones(1,numSubs)];
xMoG = [3*ones(1,numSubs) 4*ones(1,numSubs)];
y1 = geo_SGnllScaled;
y2 = nonPar_SGnllScaled;
y3 = geo_MoGnllScaled;
y4 = nonPar_MoGnllScaled;
ySG = [y1 y2];
yMoG = [y3 y4];
swarmchart(xSG,ySG,50,[0.5 1 0.5],'filled',XJitterWidth=0.2);
swarmchart(xMoG,yMoG,50,[0.5 0.5 1],'filled',XJitterWidth=0.2);

scatter([1 2],[medSGnll_geo medSGnll_nonPar],90,[0 0.75 0],'filled');
scatter([3 4],[medMoGnll_geo medMoGnll_nonPar],90,[0 0 1],'filled');

errorbar([1 2],[medSGnll_geo medSGnll_nonPar],...
               [medSGnll_geo medSGnll_nonPar] - [iqrSGnll_geo(1) iqrSGnll_nonPar(1)],...
               [iqrSGnll_geo(2) iqrSGnll_nonPar(2)] - [medSGnll_geo medSGnll_nonPar],'color',[0 0.75 0],'linewidth',2,'LineStyle','none');
errorbar([3 4],[medMoGnll_geo medMoGnll_nonPar],...
               [medMoGnll_geo medMoGnll_nonPar] - [iqrMoGnll_geo(1) iqrMoGnll_nonPar(1)],...
               [iqrMoGnll_geo(2) iqrMoGnll_nonPar(2)] - [medMoGnll_geo medMoGnll_nonPar],'color',[0 0 1],'linewidth',2,'LineStyle','none');

ylabel('Log-probability of data');
xlabel('Model');
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',15,'xtick',[1 2 3 4],'ylim',[0 1],'ytick',[0 1],...
    'yticklabel',{'Coin flip','Cum. Gauss.'},'xticklabel',{'Geo. SG.','Nonpar. SG.','Geo. MoG.','Nonpar. MoG.'});

subplot(1,2,2);
hold on;

y1 = geo_SGbic;
y2 = nonPar_SGbic;
y3 = geo_MoGbic;
y4 = nonPar_MoGbic;
ySG = [y1 y2];
yMoG = [y3 y4];
swarmchart(xSG,ySG,50,[0.5 1 0.5],'filled',XJitterWidth=0.2);
swarmchart(xMoG,yMoG,50,[0.5 0.5 1],'filled',XJitterWidth=0.2);

scatter([1 2],[medSGbic_geo medSGbic_nonPar],90,[0 0.75 0],'filled');
scatter([3 4],[medMoGbic_geo medMoGbic_nonPar],90,[0 0 1],'filled');

errorbar([1 2],[medSGbic_geo medSGbic_nonPar],...
               [medSGbic_geo medSGbic_nonPar] - [iqrSGbic_geo(1) iqrSGbic_nonPar(1)],...
               [iqrSGbic_geo(2) iqrSGbic_nonPar(2)] - [medSGbic_geo medSGbic_nonPar],'color',[0 0.75 0],'linewidth',2,'LineStyle','none');
errorbar([3 4],[medMoGbic_geo medMoGbic_nonPar],...
               [medMoGbic_geo medMoGbic_nonPar] - [iqrMoGbic_geo(1) iqrMoGbic_nonPar(1)],...
               [iqrMoGbic_geo(2) iqrMoGbic_nonPar(2)] - [medMoGbic_geo medMoGbic_nonPar],'color',[0 0 1],'linewidth',2,'LineStyle','none');

ylabel('BIC');
xlabel('Model');
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',15,'xtick',[1 2 3 4],...
    'xticklabel',{'Geo. SG.','Nonpar. SG.','Geo. MoG.','Nonpar. MoG.'});


medSGnll_geo      = median(geo_SGnllScaled,'omitnan');
medMoGnll_geo     = median(geo_MoGnllScaled,'omitnan');
medSGnll_nonPar   = median(nonPar_SGnllScaled,'omitnan');
medMogGnll_nonPar = median(nonPar_MoGnllScaled,'omitnan');

iqrSGnll_geo      = [prctile(geo_SGnllScaled,25);prctile(geo_SGnllScaled,75)];
iqrMoGnll_geo     = [prctile(geo_MoGnllScaled,25);prctile(geo_MoGnllScaled,75)];
iqrSGnll_nonPar   = [prctile(nonPar_SGnllScaled,25);prctile(nonPar_SGnllScaled,75)];
iqrMogGnll_nonPar = [prctile(nonPar_MoGnllScaled,25);prctile(nonPar_MoGnllScaled,75)];


%% Save plots

% plots
if saveOn

    saveas(f1,[saveDir,'/figures/','xSubBayes_MoGPriors.svg']);
    saveas(f2,[saveDir,'/figures/','xSubBayes_LikelihoodWidths.svg']);
    saveas(f3,[saveDir,'/figures/','xSubBayes_estimatedDistances.svg']);
    saveas(f4,[saveDir,'/figures/','xSubBayes_fitQuality.svg']);

end

% Fits
if saveOn
    
    xSubBayes.subjIDs        = goodSubs;
    xSubBayes.nonpar_mogFits = nonpar_mogFits;
    xSubBayes.nonpar_sgFits  = nonpar_sgFits;
    xSubBayes.nonpar_datastr = nonpar_datastr;
    xSubBayes.allSortedPars  = allSortedPars;
    xSubBayes.geo_mogFits    = geo_mogFits;
    xSubBayes.geo_sgFits     = geo_sgFits;
    xSubBayes.geo_datastr    = geo_datastr;
    xSubBayes.thisPrior      = thisPrior;


    save([saveDir,'xSubBayes'],'xSubBayes');

end


end