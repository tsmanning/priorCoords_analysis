function [ipd] = loadInIPD(subjIDs,saveDir)

%% Loads in CSV file of IPDs and converts to meters

% Load in csv as a table and convert to array
ipdTab  = readtable([saveDir,'ipds.csv']);
ipdSubs = table2cell(ipdTab(:,1));
ipdArr  = table2array(ipdTab(:,2));

% Match subject IDs to those in array and extract IPDs
numSubs = numel(subjIDs);

for ii = 1:numSubs
    ipdInd  = strcmp(subjIDs{ii},ipdSubs);
    ipd(ii) = ipdArr(ipdInd);
end

% Convert from mm to m
ipd = ipd*0.001;

end