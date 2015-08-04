function [beta,kappa]=betamim(X,y,h)
% BETAMIM        Kolmogorov Smirnov Moment Independent Measure
%
% B=BETAMIM(X,Y) computes the mean distance between Y and Y given X
%                        using the two-sided KS test statistics
% [B,K]=BETAMIM(X,Y) also returns the mean Kuiper distance K
%
% ...=BETAMIM(X,Y,H) additionally specifies the partition size H
%
% written by Emanuele Borgonovo, emanuele.borgonovo@unibocconi.it
%        and Elmar Plischke , elmar.plischke@tu-clausthal.de
[n,k]=size(X);
if(nargin<3)
h=min(ceil(n^(1/(2+(n<1500)))),50);
end
W=round(n/h);

[~,Indx]=sort(X);MM=y(Indx);

s=zeros(h,k);
t=zeros(h,k);
% loop on the partitions
for j=1:h
    fX=MM((j-1)*W+1:min(j*W,n),:);
    % loop on the factors
    for i=1:k
     [s(j,i),t(j,i)]=ksstat(y,fX(:,i));
    end
end
beta=mean(s);kappa=mean(t);
end

function [KS,Kui]=ksstat(x1,x2)
% Calculate F1(x) and F2(x), the empirical (i.e., sample) CDFs.
% (from kstest2.m)

binEdges    =  [-inf ; sort([x1;x2]) ; inf];

binCounts1  =  histc (x1 , binEdges, 1);
binCounts2  =  histc (x2 , binEdges, 1);

sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);

deltaCDF    = sampleCDF1 - sampleCDF2

KS  =  max(abs(deltaCDF));
Kui =  max(deltaCDF) - min(deltaCDF);
end

