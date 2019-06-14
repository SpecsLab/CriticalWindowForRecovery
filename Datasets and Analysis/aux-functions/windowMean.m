% Calculates Mean over fixed intervals
% Author Armin Duff armin.duff@gmail.com 2018
% ---------------------------------------------------------------- %

function [means,windowLimits] = slidettest2(rate,rgs,time,n)
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
    
    means(i,:)=[mean(rate_(rgs_==0)), mean(rate_(rgs_==1))];
end
    
