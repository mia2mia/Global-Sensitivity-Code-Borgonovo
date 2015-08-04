function [b,d,t,e,w,bm,dm,tm,em,wm]=betaKS2(x,y,M,gfx)
% BETAKS Kolmogorov-Smirnov Sensitivity.
%[B,D,K,E,W]=BETAKS2(X,Y) returns the following sensitivity indicators
%      B Kolmogorov-Smirnov Distance Beta
%      D Borgonovo Delta
%      K Kuiper Discrepancy Kappa
%      E Pearson Correlation Ratio Eta2
%      W Wald Wolfowitz number of runs (normalized)

%      Data X and Y must contain unique values per column

%      BM, DM, KM, EM, WM Conditional Shift/Separation Measurements

%      
% written by elmar.plischke@tu-clausthal.de
[n,k]=size(x);
if(size(y,1)~=n), error('BETAKS: size mismatch.'); end
gfxrows=2;

if(nargin<3), M=32;   end
if(M==0)
    % test debiasing: predict value at partition size M=0
    [b1,d1,t1,e1]=betaKS2(x,y,32); 
    [b2,d2,t2,e2]=betaKS2(x,y,64); % double partition
    b=2*(b1-.5*b2);
    d=2*(d1-.5*d2);
    t=2*(t1-.5*t2);
    e=2*(e1-.5*e2);
    return
end
% cumulative distribution: with ties
[ys,indx,rindx]=unique(y);
nn=size(ys,1);
if(nn==n)
 cdf=(1:n)/n;
else
 rle=zeros(1,nn); 
 % run length encoder
 for i=rindx'; rle(i)=rle(i)+1; end; 
% cdf_=(cumsum([0,rle(1:end-1)])+(rle+1)/2)/n;
 cdf_=cumsum(rle)/n;
 cdf= cdf_(sort(rindx));
 [ys,indx]=sort(y);  
end
xs=zeros(n,k);
for i=1:k
 xs(:,i)=x(indx,i);
end

VY=var(y,1);

if(nargin<4), gfx=''; else cols=jet(M);end
m=linspace(0,1,M+1);
bm=zeros(k,M);
dm=zeros(k,M);
tm=zeros(k,M); % test
wm=zeros(k,M); % more test
em=ones(k,M)*VY/n; % and even more
nm=zeros(k,M);
%% keep only values from the partition
[xr,indxx]=sort(xs);
%% test for constant input
factorlist=find(xr(1,:)~=xr(end,:));
if(any(xr(1,:)==xr(end,:)))
    disp('BETAKS: Constant input factors detected.');
end
%% 
for i=factorlist % 1:k
%   xr(indxx(:,i),i)=1:n; % ranks (no ties)
   xr(indxx(:,i),i)=cumsum([1;diff(xr(:,i))~=0]); % cheap tied ranks
end
xr=xr./(ones(n,1)*max(xr)); % scale to 1
%%
for j=1:M
   indx= (m(j)<xr) & (xr <= m(j+1));
   %xtreme=xor(indx(2:end,:),indx(1:end-1,:));
   %xtreme(2:end,:)=xtreme(2:end,:)|xtreme(1:end-1,:);
   xtreme=indx; % test
   cdfc=cumsum(indx); 
   nm(:,j)=cdfc(end,:)'; %sum(indx); % indx(end,:); % if no ties: always same nr. of realizations
   %% wenn y(1) nicht drin, dann y(1) mit 0 codieren, ebenso y(end) mit 1
   xtreme(1,:)=1;xtreme(n,:)=1;
   %%
   for i=factorlist %1:k
   % sum of conditionals 
   scc=nm(i,j);
   if ~scc, continue, end
   %difference of conditionals
   dcdf=cdf(xtreme(:,i))-cdfc(xtreme(:,i),i)'/scc;
   if(~isempty(gfx))
    subplot(gfxrows,ceil(k/gfxrows),i)
    %plot(ys(xtreme(:,i)),dcdf,'Color',cols(j,:));hold on
    if j==1
	 stairs(ys,cdf,'k','LineWidth',4);hold on;
    end
	[ss,tt]=stairs(ys(xtreme(:,i)),cdfc(xtreme(:,i),i)'/scc);
    plot(ss,tt,'Color',cols(j,:));
   end
   mx=max(dcdf);mn=min(dcdf);
   % Kolmogorov Smirnov distance 
   bm(i,j)=max([0,mx,-mn]);
   % Kuiper dicrepancy 
   tm(i,j)=max([0,mx-mn]);    
   %%  find max in each positive run, min in each negative
   delta=0.0;
   sgn=0; exval=0.0;
   for val=dcdf
    sg=sign(val);
    if sg==sgn
     exval=max([exval,sgn*val]);
    else 
	 if(exval>0.1*sqrt(1/n+1/scc))
      delta=delta+exval;
	 end
	 sgn  =sg;
	 exval=sg*val;
	end %what if zero?
   end
   delta=delta+exval;
   dm(i,j)=delta;
   %% and more indicators
   % Wald-Wolfowitz number of runs
   mju=2*scc*(n-scc)/n-1;
   wm(i,j)=(sum(xtreme(:,i))-mju)/sqrt(mju*(mju+1)/(n-1));
   % Correlation ratio
   em(i,j)=var(ys(indx(:,i)),1);
   end % for
end
if(~isempty(gfx))
 for i=factorlist %1:k;
     subplot(gfxrows,ceil(k/gfxrows),i);
     %title(gfx);xlabel('y');ylabel('\Delta cdf');
     title(gfx);xlabel('y');ylabel('cdf');
     hold off;
 end
end
% b=mean(bm,2)';
% d=mean(dm,2)';
% t=mean(tm,2)';
% w=mean(wm,2)';
% f=1-mean(em,2)'./var(y);
 b=sum(bm.*nm,2)'/n;
 d=sum(dm.*nm,2)'/n;
 t=sum(tm.*nm,2)'/n;
 w=sum(wm.*nm,2)'/n;
 f=1-sum(em.*nm,2)'./n/VY;
 e=zeros(1,k);
 e(factorlist)=f(factorlist);
% e=1-(n-M)/(n-1)*mean(em,2)'/var(y); % bias correction
end

function testbeta
%%
ishigami;
n=2^14-1;
x=trafo(sobolseq(n,k));y=model(x);
%%
bs=[];ds=[];ts=[];
Ds=[];Ts=[];
Ps=2:200;
h=waitbar(0,'Computation in progress');
for i=Ps;
    waitbar(sqrt(i*(i+1))/Ps(end),h);
    %if(mod(i,50)==0), fprintf('.'); drawnow('update'); end
    [b,d,t]=betaKS2(x,y,i);bs=[bs;b];ds=[ds;d];ts=[ts;t];
%    Ds=[Ds; deltabork(x,y,struct('PartitionSize',i,...
%        'KSLevel',0,'KDWidth','clone','QuadraturePoints',500))];
    Ts=[Ts; deltafast2(x,y,i)];
end
close(h);
%%
clf
plot(Ps,ts);hold on; 
plot(Ps,ds,'--');plot(Ps,Ts,'-.');hold off
%%

end
 