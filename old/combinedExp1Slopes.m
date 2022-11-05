
close all

combinedData = [stats.slopeData{2} stats.slopeData{4}]';
combinedDataLB = [stats.slopeDataLB{2} stats.slopeDataLB{4}]';
combinedDataUB = [stats.slopeDataUB{2} stats.slopeDataUB{4}]';

%%%%%%%%%

numSubs = size(combinedData,1);
theseDists = [0.5 1];
spread   = linspace(-0.1,0.1,numSubs);

fig = figure;
fig.Position = [800 100 650 650];
hold on;

for ii = 1:numSubs

    thisColor  = 0.75*[1 1 1];

    theseSlopes   = combinedData(ii,:);
    theseSlopeLBs = combinedDataLB(ii,:);
    theseSlopeUBs = combinedDataUB(ii,:);

    errorbar(theseDists + spread(ii),theseSlopes,theseSlopes-theseSlopeLBs,theseSlopeUBs-theseSlopes,'color',thisColor,'linewidth',2);

    if ii > 6
        % Low Cont Ref
        p2(ii) = scatter(theseDists + spread(ii),theseSlopes,150,thisColor);
    else
        p2(ii) = scatter(theseDists + spread(ii),theseSlopes,150,thisColor,'filled');
    end

end

% Plot group mean and STD
thisSlopeData = combinedData';

meanSlopes    = mean(thisSlopeData,2);
[~,~,groupSlopesCI(1,:)]  = ttest(thisSlopeData(1,:));
[~,~,groupSlopesCI(2,:)]  = ttest(thisSlopeData(2,:));

errorbar(theseDists,meanSlopes,meanSlopes-groupSlopesCI(:,1),groupSlopesCI(:,2)-meanSlopes,'color',[0 0 0],'linewidth',2);
p2b = scatter(theseDists,meanSlopes,200,[0 0 0],'filled');

pSRSlope       = signrank(thisSlopeData(1,:),thisSlopeData(2,:));
[~,pNormSlope] = ttest(thisSlopeData(1,:),thisSlopeData(2,:));

set(gca,'fontsize',20,'xlim',[0 1.5],'xtick',[0.5 1]);
xlabel('Test Distance (m)');
ylabel('Slope (\Delta p/\Delta log(V))');
text(0.25,2,['p = ',num2str(round(pNormSlope,5))],'fontsize',20);