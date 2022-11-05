function [pFxn,PSE,sens] = fitCumGauss(x,p,varargin)

% Fit weibull function to psychometric data

if numel(varargin) == 0
    wts = ones(size(p,1),size(p,2));
else
    wts = varargin{1};
end

% Make sure all inputs are column vectors
if ~iscolumn(x)
    x = x';
end

if ~iscolumn(p)
    p = p';
end

if ~iscolumn(wts)
    wts = wts';
end

lb = [0 0];
ub = [100 100];

opts = optimset('Algorithm', 'interior-point', 'Display', 'off', ...
    'MaxFunEvals', 5000, 'MaxIter', 500, 'GradObj', 'off');

% [mu sig]
parVec0 = [median(x) median(x)-min(x)];

lossFxn = @(prs) sum( wts.*(normcdf(x,prs(1),prs(2)) - p).^2 ,'omitnan');

pars  = fmincon(lossFxn,parVec0,[],[],[],[],lb,ub,[],opts);

if numel(varargin) == 2
    pFxn = normcdf(varargin{2},pars(1),pars(2));
else
    pFxn = normcdf(x,pars(1),pars(2));
end

PSE = pars(1);

sens = pars(2);

end