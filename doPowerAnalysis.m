function [nCrit,stats] = doPowerAnalysis(x,y,alpha,tails)

% Find number of subjects to reach critical value given paired means/STDs
if ~iscolumn(x)
    x = x';
end
if ~iscolumn(y)
    y = y';
end

m1 = mean(x);
m2 = mean(y);
s1 = std(x);
s2 = std(y);
rho = corr(x,y);

tFxn = @(m1,m2,s1,s2,rho,n) (abs(m2-m1)/sqrt(s1^2 + s2^2 - 2*rho*s1*s2))*sqrt(n-1);
% pFxn = @(df,t) gamma((df+1)/2)/(sqrt(df*pi)*gamma(df/2))*(1 + (t^2)/df)^-((df+1)/2);

n = 1:500;

for ii = 1:numel(n)

    thisdf = n(ii) - 1;

    lowTail = icdf('T',alpha/tails,thisdf-1);

    x = lowTail + tFxn(m1,m2,s1,s2,rho,n(ii)+1);

    p(1,ii) = cdf('T',x,thisdf-1);

    if p(1,ii) >= (1-alpha) %(1-alpha/tails)
        nCrit = n(ii);
        break
    end

    if ii == 500
        nCrit = nan;
    end

end

stats.means = [m1,m2];
stats.stds  = [s1,s2];
stats.corr  = rho;
stats.p     = p(ii);
stats.effSz = tFxn(m1,m2,s1,s2,rho,n(ii))/sqrt(n(ii)-1);
stats.nvp   = [n(1:ii);p];

end