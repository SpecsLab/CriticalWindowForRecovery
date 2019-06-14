% Performs difference test on fixed intervals. Change
% the test in the file to use different tests

% Author Armin Duff armin.duff@gmail.com 2018
%-------------------------------------------------------

function [p,windowLimits,stats, pT1, pT2, pT3] = windowDiffTest(rate,rgs,time,n)
nV = numel(rate) ;

if(numel(n)==1)
    indexLimits=floor(linspace(1,nV+1,n+1));
else
    indexLimits(1)=1;
    for ii=1:numel(n)
        indexLimits(ii+1)=find(time>=n(ii),1);
    end
    indexLimits(end+1)=numel(time);
end

for i=1:numel(indexLimits)-1
    windowLimits(i,:)=[time(indexLimits(i)),time(indexLimits(i+1)-1)];
    
    rate_=rate(indexLimits(i):indexLimits(i+1)-1);
    rgs_=rgs(indexLimits(i):indexLimits(i+1)-1);
    
    stats(i,:)=[mean(rate_(rgs_==0)),mean(rate_(rgs_==1)),...
        std(rate_(rgs_==0))/sqrt(numel(rate_(rgs_==0))),std(rate_(rgs_==1))/sqrt(numel(rate_(rgs_==1))),...
        numel(rate_(rgs_==0)),numel(rate_(rgs_==1))];
    
    [p(i),h] = ranksum(rate_(rgs_==0),rate_(rgs_==1));
   
end

for i=1:numel(indexLimits)-1    
    rate_=rate(indexLimits(i):indexLimits(i+1)-1);
    rgs_=rgs(indexLimits(i):indexLimits(i+1)-1);
    aux_rate= rate;
    aux_rate(rgs==0)=nan;
    aux_rate(rgs==2)=nan;
    
    [pT1(i),h1] = ranksum(aux_rate(indexLimits(1):indexLimits(2)-1),rate_(rgs_==1));
    [pT2(i),h2] = ranksum(aux_rate(indexLimits(2):indexLimits(3)-1),rate_(rgs_==1));
    [pT3(i),h3] = ranksum(aux_rate(indexLimits(3):indexLimits(4)-1),rate_(rgs_==1));
   
end



    
