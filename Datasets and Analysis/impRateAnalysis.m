%function [cleanData] = impRateAnalysis(data, plotting)
%CLEANDATA Reads, cleans, and structures data from:
%
%   Workbook: Clinical_Scales_all.csv
%
%   For analysis in:
%
%   For the recovery rates analysis in:
%
%   "The impact of VR-based rehabilitation on 
%   post-stroke functional recovery: a retrospective meta-analysis"
% 
%   This script prepares data and plots Figure 2 if var plotting set to
%   true
% 
%   Author Armin Duff armin.duff@gmail.com 
%   and Belen Rubio Ballester belen.rubio.ballester@gmail.com 2018
% ---------------------------------------------------------------- %

function cleanData=impRateAnalysis(data, plotting)

    if(nargin<2)
        plotting=false;
    end

    %% wirte DATA to readable variables
    scoreAll = cell2mat(data(:,6));
    rgsAll=cell2mat(data(:,7));
    rateAll=cell2mat(data(:,8));
    timeAll=cell2mat(data(:,9));
    patientAll=cell2mat(data(:,10));
    studyAll=cell2mat(data(:,12));
    

    %% Names and Max for scales
    scaleName = {'Fugl-Meyer','Cahai', 'Barthel'};
    scaleMax={66,91, 100};

    %% prepare DATA
    cleanData=cell(2,12);
    
    %% Load corrected statistics
    if(exist('statsCorrected.mat', 'file')==2&&exist('statsCorrected.mat', 'file')==2)
        load pCorrected.mat
        load statsCorrected.mat
    else
        plotting=false;
    end

    for ii=1:2 %  for all scales Fuglmeyer %Cahai

        rate=rateAll(:,ii);
        score=scoreAll(:,ii);
        time=timeAll;
        patient=patientAll;
        study=studyAll;

        cleanIndex=~isnan(rate)&~isnan(time);

        %rgs = rgsAll(cleanIndex)==1; %RGS yes no
        rgs = rgsAll(cleanIndex); % CHANGE BY BR OCTOBER
        %rgs(rgs==0)=0;
        %rgs(rgs==2)=2;
        %rgs(rgs==0)=3;
        %rgs(rgs==2)=0;
        controlOnTreatment = rgsAll(cleanIndex)==2; %RGS yes no
        rate=rate(cleanIndex);
        time=time(cleanIndex);
        patient=patient(cleanIndex);
        study=study(cleanIndex);
        
        cleanIndex_=~isnan(score);
        score=score(cleanIndex);

        [time,sortIndex]=sort(time); %sort following time
        rate=rate(sortIndex,:);
        rgs=rgs(sortIndex);
        patient=patient(sortIndex);
        study=study(sortIndex);
        score=score(sortIndex,:);

        cleanData{ii,1}=time;
        cleanData{ii,2}=rate;
        cleanData{ii,3}=rgs;
        cleanData{ii,13}=patient;
        cleanData{ii,14}=study;
        cleanData{ii,15}=score;

        limits=[180,1.5*365];
        M=3;

        [cleanData{ii,4},cleanData{ii,5},cleanData{ii,6} pT1, pT2, pT3]=windowDiffTest(cleanData{ii,2},rgs,time,limits); % p-value, WindowLimits, means

        W=100;
        WMode='forward';
        rateNORGS=rate;
        rgs = rgsAll(cleanIndex);
        rateNORGS(rgs==1)=nan;
        rateRGS=rate;
        rateRGS(rgs==0)=nan;
        rateRGS(rgs==2)=nan;

        %medians for RGS and no RGS
        cleanData{ii,7}=slidefun(@nanmean,W,rateNORGS, WMode);
        cleanData{ii,8}=slidefun(@(x) prctile(x,25),W,rateNORGS, WMode);
        cleanData{ii,9}=slidefun(@(x) prctile(x,75),W,rateNORGS, WMode);

        cleanData{ii,10}=slidefun(@nanmean,W,rateRGS, WMode);
        cleanData{ii,11}=slidefun(@(x) prctile(x,25),W,rateRGS,WMode);
        cleanData{ii,12}=slidefun(@(x) prctile(x,75),W,rateRGS,WMode);

    end

    for ii=1:2

        % Get DATA for plotting    
        time=cleanData{ii,1};
        rate=cleanData{ii,2};
        rgs=cleanData{ii,3};
        p=cleanData{ii,4};
        pCor=pCorrected(ii,:);
        windowLimits=cleanData{ii,5};

        stats = cleanData{ii,6};

        meansNORGS=cleanData{ii,7};
        lowNORGS=cleanData{ii,8};
        upNORGS=cleanData{ii,9};

        meansRGS=cleanData{ii,10};
        lowRGS=cleanData{ii,11};
        upRGS=cleanData{ii,12};
        
    end

    %% PREPARE ESTIMATES FOR BAR CHART
    meansBarCorrectedT1=statsCorrected{1}(:,1:2);
    meansBarT1=(stats(:,1:2)-meansBarCorrectedT1+repmat(meansBarCorrectedT1(:,1),1,2))*7/scaleMax{1}*100;
    stdsBarT1=stats(:,3:4)*7/scaleMax{1}*100;

    meansBarCorrectedFU=statsCorrected{2}(:,1:2);
    meansBarFU=(stats(:,1:2)-meansBarCorrectedFU+repmat(meansBarCorrectedFU(:,1),1,2))*7/scaleMax{2}*100;
    stdsBarFU=stats(:,3:4)*7/scaleMax{2}*100;

    str = {'Subacute';'Chronic';'Late Chronic'};
    str_y = {'UE-FM Norm. Improvement';'CAHAI Norm. Improvement'};
    
%% PLOTTING FIGURE 2:
    
    if(plotting)
   
        clear g
        f = figure ('position', [100 100 1000 400]);
        g(1,1)=gramm('x',repmat(str, 1, 2),'y',reshape(meansBarT1, 1, 6),...
    'ymin',reshape(meansBarT1, 1, 6)-(reshape(stdsBarT1, 1, 6).*1.16),'ymax',reshape(meansBarT1, 1, 6)+(reshape(stdsBarT1, 1, 6).*1.16),'color', [{'FU'} {'FU'} {'FU'} {'RGS'} {'RGS'} {'RGS'}]);
        g(1,1).set_names('x',' ','y',str_y(1), 'color','# Group');
        g(1,1).geom_bar('dodge',0.8,'width',0.6);
        g(1,1).geom_interval('geom','black_errorbar','dodge',0.8,'width',1);
        g(1,1).axe_property('XTickLabelRotation',60); %Should work for recent Matlab versions
        g(1,1).set_order_options('x',0);
        g(1,1).set_layout_options('Position',[0.07 0.6 0.2 0.37],'legend_pos',[0.25 0.75 0.2 0.2],'margin_height',[0.05 0.02],...
    'margin_width',[0.02 0.02],...
    'redraw',false);
        g(1,1).set_color_options('map','d3_10');
        
        g(2,1)=gramm('x',repmat(str, 1, 2),'y',reshape(meansBarFU, 1, 6),...
    'ymin',reshape(meansBarFU, 1, 6)-(reshape(stdsBarFU, 1, 6).*1.16),'ymax',reshape(meansBarFU, 1, 6)+(reshape(stdsBarFU, 1, 6).*1.16),'color', [{'FU'} {'FU'} {'FU'} {'RGS'} {'RGS'} {'RGS'}]);
        g(2,1).set_names('x',' ','y',str_y(2),'color','# Group');
        g(2,1).geom_bar('dodge',0.8,'width',0.6);
        g(2,1).geom_interval('geom','black_errorbar','dodge',0.8,'width',1);
        g(2,1).axe_property('XTickLabelRotation',60); %Should work for recent Matlab versions
        g(2,1).set_layout_options('Position',[0.07 0.05 0.2 0.37], 'legend_pos',[0.25 0.25 0.2 0.2],'margin_height',[0.2 0.02],...
    'margin_width',[0.02 0.02],...
    'redraw',false);
        g(2,1).set_color_options('map','d3_10');
        g(2,1).set_order_options('x',0);


    end

    for ii=1:2
        % compute stats
        p=cleanData{ii,4};
        pCor=pCorrected(ii,:);
        windowLimits=cleanData{ii,5};

        meansNORGS=cleanData{ii,7}*7/scaleMax{ii}*100;
        lowNORGS=cleanData{ii,8}*7/scaleMax{ii}*100;
        upNORGS=cleanData{ii,9}*7/scaleMax{ii}*100;

        meansRGS=cleanData{ii,10}*7/scaleMax{ii}*100;
        lowRGS=cleanData{ii,11}*7/scaleMax{ii}*100;
        upRGS=cleanData{ii,12}*7/scaleMax{ii}*100;
    end         
   
    % FM
    timeFM=cleanData{1,1};
    rateFM=cleanData{1,2}*7/scaleMax{1}*100;
    rgs=cleanData{1,3};
    cats=repmat({'RGS'}, length(rgs),1);
    cats(find(rgs==0)) = {'FU'};
    cats(find(rgs==2)) = {'OT'};

        
    if(plotting)
            
        g(1,2) = gramm('x', timeFM,'y', rateFM, 'color', [cats]);
        g(1,2).geom_point();
        g(1,2).stat_smooth('geom', {'line'});
        g(1,2).set_names('x','Days Post-Stroke','y',str_y(1), 'color', '# Group');
        g(1,2).axe_property('XScale','log');
        g(1,2).set_layout_options('Position',[0.4 0.6 0.3 0.37], 'legend_pos',[0.65 0.75 0.2 0.2],'margin_height',[0.05 0.02],...
    'margin_width',[0.02 0.02],...
    'redraw',false);
        g(1,2).set_color_options('map','brewer2');
        g(1,2).set_order_options('x',0, 'color',0);
    end
    
    
    %CAHAI
    timeCAHAI=cleanData{2,1};
    rateCAHAI=cleanData{2,2}*7/scaleMax{2}*100;
    rgs=cleanData{2,3};
    cats=repmat({'RGS'}, length(rgs),1);
    cats(find(rgs==0)) = {'FU'};
    cats(find(rgs==2)) = {'OT'};
        
     if(plotting)
        g(2,2) = copy(g(1));
        g(2,2) = gramm('x', timeCAHAI,'y', rateCAHAI, 'color', [cats]);
        g(2,2).geom_point();
        g(2,2).stat_smooth('geom', {'line'});
        g(2,2).set_names('x','Days Post-Stroke','y',str_y(2), 'color', '# Group');
        g(2,2).set_layout_options('Position',[0.4 0.05 0.3 0.37], 'legend_pos',[0.65 0.25 0.2 0.2],'margin_height',[0.2 0.02],...
    'margin_width',[0.02 0.02],...
    'redraw',false);
        g(2,2).set_color_options('map','brewer2');
        g(2,2).axe_property('XScale','log');
        g(2,2).set_order_options('x',0, 'color',0);
        g.draw();

        % Indicate Stages
        text(22,22,{'Acute/Subacute'},'Parent',g(1,2).facet_axes_handles,'FontName','Helvetica');
        line(repmat(limits(1),1,2),[-6 22],'Color','k','LineStyle','--','Parent',g(1,2).facet_axes_handles);
        text(limits(1)+5,22,{'Chronic'},'Parent',g(1,2).facet_axes_handles,'FontName','Helvetica');
        line(repmat(limits(2),1,2),[-6 22],'Color','k','LineStyle','--','Parent',g(1,2).facet_axes_handles);
        text(limits(2)+5,22,{'Late' 'Chronic'},'Parent',g(1,2).facet_axes_handles,'FontName','Helvetica');
        text(22,28,{'Acute/Subacute'},'Parent',g(2,2).facet_axes_handles,'FontName','Helvetica');
        line(repmat(limits(1),1,2),[-10 28],'Color','k','LineStyle','--','Parent',g(2,2).facet_axes_handles);
        text(limits(1)+5,28,{'Chronic'},'Parent',g(2,2).facet_axes_handles,'FontName','Helvetica');
        line(repmat(limits(2),1,2),[-10 28],'Color','k','LineStyle','--','Parent',g(2,2).facet_axes_handles);
        text(limits(2)+5,28,{'Late' 'Chronic'},'Parent',g(2,2).facet_axes_handles,'FontName','Helvetica');

        saveas(gcf,'figures/Figure2.eps','epsc')

     end
         
end
