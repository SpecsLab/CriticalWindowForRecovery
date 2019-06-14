
%   Main Script for the full analysis of:
%
%   "A critical time window for recovery extends beyond one-year post-stroke"
% 
%   This script imports data, calculates the corrected p-values 
%   and expected means, generates demographics table, performs 
%   statistical tests within conditions and plots
%   all Figures and Tables included in the manuscript and Supporting Information.
% 
%   Author Armin Duff armin.duff@gmail.com 
%   and Belen Rubio Ballester belen.rubio.ballester@gmail.com 2018
% ---------------------------------------------------------------- %

% Dependencies
addpath(genpath('dependencies'))
addpath(genpath('aux-functions'))

% Clinical dataset
filename = 'data/Clinical_Scales_all.csv'
[data] = importAll(filename);

% Parse clinical data to structure
Headings = {'patientID' 'Chronicity' 'Treatment' 'ChronicityDays' 'EvalTime' 'Scale' 'TreatmentFlag' 'ImprovementRate' 'ImprovementTime' 'PatientIDs' 'StudyId' 'StudyIDs' 'AllImprovement' 'TreatmentImprovement' 'FollowUpImprovement', 'OrderingIndex', 'NamesConditions'};
dataStructure = cell2struct(data, Headings, 2);
dataStructure = nestedSortStruct(dataStructure, 'OrderingIndex');
[studiesNames, ~, ib_group] = unique({dataStructure.NamesConditions}', 'stable');
conditions=ib_group;

N=10000; %repetitions of analysis
limits=[180,1.5*365]; %predefined limits for analysis half year (subacute to chronic), more than 1.5 years (late chronic)
M=3; %number of intervals

W=25; %window length for data resampling
WMode='central'; %window type for data resampling

out=impRateAnalysis(data, true);

%% Compute Corrected P-Values
pValues=zeros(2,M,N);
pValuesKS=zeros(2,N);
for ii=1:2 %for all three scales 
    
    rgs=out{ii,3};
    rgs(rgs~=1)=0;
    time=out{ii,1};
    rate=out{ii,2};
    rateNoRGS=rate;
    rateNoRGS(rgs==1)=nan;
    rateNoRGS(rgs==2)=nan;
    
    empRandRate=slidefun(@(x,n)emprand(x,n,1),W,rateNoRGS,WMode,N); %generate empirical random data from noRGS data
    
    for kk=1:N  %for all repetitions             
        [h,pValuesKS(ii,kk)]=kstest2(rate(rgs==0),empRandRate(:,kk)); %Single sample Kolmogorov-Smirnov to compare generated data and original data
        pValues(ii,:,kk)=windowDiffTest(empRandRate(:,kk),rgs,time,limits); %calcualte significance of generated data
    end
end

x=zeros(2,M);
fx=zeros(2,M);
p0=0.05;

%Get median p-values
for ii=1:2
    for jj=1:M
        pVal=squeeze(pValues(ii,jj,:));
        pValSort=sort(pVal);
        cutOff=floor(N*p0);
        x(ii,jj)=(pValSort(cutOff)+pValSort(cutOff-1))/2; %calculate corrected p Values :median p-val
        fx(ii,jj)=sum(pVal<x(ii,jj))/N; %test corrected p Values
    end
end

x(x>0.05)=0.05; %crop p Values at 0.05 to avoid higher than 0.05 p values

pKS=sum(sum(pValuesKS<.05))/3/N;  %Check how often the generated data was significantly different from the original data. This number should be low <0.05

pCorrected=x;
save pCorrected.mat pCorrected

%% Compute Expected Means

means=zeros(2,M,2,N);
pValuesKS=zeros(2,N);

for ii=1:2
    time=out{ii,1};
    rate=out{ii,2};
    rgs=out{ii,3};
    rgs(rgs~=1)=0;
    
    rateNoRGS=rate;
    rateNoRGS(rgs==1)=nan;
    
    empRandRate=slidefun(@(x,n)emprand(x,n,1),W,rateNoRGS,WMode,N);
    
    for kk=1:N
        [h,pValuesKS(ii,kk)]=kstest2(rate(rgs==0),empRandRate(:,kk)); 
        means(ii,:,:,kk)=windowMean(empRandRate(:,kk),rgs,time,limits);
    end
end

%plot means to check corrections
for ii =1:2
    meansBar=mean(squeeze(means(ii,:,:,:)),3);
    stdsBar=std(squeeze(means(ii,:,:,:)),0,3);
    statsCorrected{ii}=[meansBar,stdsBar];
    
end

pKS=sum(sum(pValuesKS<.05))/3/N;  %Check how often the generated data was significantly different from the original data. This number should be low <0.05
save statsCorrected.mat statsCorrected;

%% Report on Demographics by Study
filenameDemog = 'patients_Demographics_and_ClinicalScreening.xlsx';
[summaryTable, dataDemographics] = importDemog(filenameDemog, [dataStructure.patientID], conditions);
totalPatients= nansum(summaryTable.GroupCount);

%% Create tables:

% TABLE 1: Demographics Summary

for i=1:length(studiesNames)
    RGSStudies(i) = ~isempty(cell2mat(regexp(studiesNames(i),'RGS|rgs|AM')));
end

ControlStudies= find(RGSStudies==0);
RGSStudies = find(RGSStudies);

%Demographics RGS interventions:
summaryTable(find(RGSStudies),:);
%Demographics Control interventions:
summaryTable(find(ControlStudies),:);

save tables/SummaryTable % store demographics summary table for all conditions


% TABLE 2: Impact on improvement at T1 and FU

for(studyId=1:max(conditions))
    ind= find(conditions==studyId);
    datapoints= reshape([dataStructure(ind).TreatmentImprovement], 3, length([dataStructure(ind).TreatmentImprovement])/3)';
    medianImprPerStudy(studyId, 1:3) = nanmedian(datapoints);
    varImprPerStudy(studyId, 1:3) = mad(datapoints);
    
    [hNorm(studyId, 1) p] = kstest(datapoints(:,1))
    [hNorm(studyId, 2) p] = kstest(datapoints(:,2))
    
    if (hNorm(studyId, 1))
        [pval(studyId, 1), h]= signrank(datapoints(:,1)); %FM
    else
        [h, pval(studyId, 1)]= ttest(datapoints(:,1)); %FM
    end
    if (hNorm(studyId, 2))
        [pval(studyId, 2), h]= signrank(datapoints(:,2)); %CAHAI
    else
        [h, pval(studyId, 2)]= ttest(datapoints(:,2)); %CAHAI
    end
    
    
    datapointsFU= reshape([dataStructure(ind).FollowUpImprovement], 3, length([dataStructure(ind).FollowUpImprovement])/3)';
    medianImprPerStudyFU(studyId, 1:3) = nanmedian(datapointsFU);
    varImprPerStudyFU(studyId, 1:3) = mad(datapointsFU);
    
    [hNormFU(studyId, 1) p] = kstest(datapointsFU(:,1))
    [hNormFU(studyId, 2) p] = kstest(datapointsFU(:,2))
    
    if (hNormFU(studyId, 1))
        [pvalFU(studyId, 1), h]= signrank(datapointsFU(:,1)); %FM
    else
        [h, pvalFU(studyId, 1)]= ttest(datapointsFU(:,1)); %FM
    end
    if (hNormFU(studyId, 2))
        [pvalFU(studyId, 2), h]= signrank(datapointsFU(:,2)); %CAHAI
    else
        [h, pvalFU(studyId, 2)]= ttest(datapointsFU(:,2)); %CAHAI
    end
end

%Create a summary table of the improvement metrics
ImprovementTable = table;
ImprovementTable.conditionId=[1:17]';
ImprovementTable.medianFM_T1=medianImprPerStudy(:,1);
ImprovementTable.madFM_T1=varImprPerStudy(:,1);
ImprovementTable.medianCAHAI_T1=medianImprPerStudy(:,2);
ImprovementTable.madCAHAI_T1=varImprPerStudy(:,2);
ImprovementTable.pFM_T1=pval(:,1);
ImprovementTable.pCAHAI_T1=pval(:,2);
ImprovementTable.medianFM_FU=medianImprPerStudyFU(:,1);
ImprovementTable.madFM_FU=varImprPerStudyFU(:,1);
ImprovementTable.medianCAHAI_FU=medianImprPerStudyFU(:,2);
ImprovementTable.madCAHAI_FU=varImprPerStudyFU(:,2);
ImprovementTable.pFM_FU=pvalFU(:,1);
ImprovementTable.pCAHAI_FU=pvalFU(:,2);

save tables/ImprovementTable

%% Plot violins and boxplots

%RGS T1
stringConditions= {'2' '3' '7' '9' '10' '12' '13' '14' '15' '16' '17'};
ChronicityIndex = ordinal([dataStructure.OrderingIndex],{'Acute', 'Subacute','Chronic','Late Chronic'},...
                       [],[0 30 180 2.5*365 nanmax([dataStructure.OrderingIndex])]);
figure ('position', [100 100 1000 370]);
ScalesImprovements= [dataStructure.TreatmentImprovement];
clear g
g(1,1) = gramm('x', {dataStructure.NamesConditions}','y', [ScalesImprovements(1:3:length(ScalesImprovements))], 'color', ChronicityIndex, 'subset',ismember(conditions,RGSStudies));
g(1,1).set_names('x','Condition ID','y','d Fugl-Meyer','color','# Chronicity');
g(1,1).set_title('Recovery from baseline to end of treatment (T1)');
g(1,1).stat_violin('half',false,'normalization','width','dodge',0,'fill','transparent');
g(1,1).stat_boxplot('width',1);
g(1,1).set_color_options('map','brewer2');
g(1,1).axe_property('xticklabels',stringConditions); %Should work for recent Matlab versions
g(1,1).set_order_options('x',0, 'color',0);
ScalesImprovements= [dataStructure.TreatmentImprovement];
g(2,1) = gramm('x', {dataStructure.NamesConditions}','y', [ScalesImprovements(2:3:length(ScalesImprovements))], 'color', ChronicityIndex, 'subset',ismember(conditions,RGSStudies));
g(2,1).set_names('x','Condition ID','y','d Cahai','color','# Chronicity');
g(2,1).stat_violin('half',false,'normalization','width','dodge',0,'fill','transparent');
g(2,1).stat_boxplot('width',1);
g(2,1).set_order_options('x',0, 'color',0);
g(2,1).axe_property('xticklabels',stringConditions); %Should work for recent Matlab versions
g(2,1).set_color_options('map','brewer2');
g.draw();
saveas(gcf,'figures/Figure1.eps','epsc')

% RGS Follow-Up
figure ('position', [100 100 1000 370]);
ScalesImprovements= [dataStructure.FollowUpImprovement];
clear g
g(1,1) = gramm('x', {dataStructure.NamesConditions}','y', [ScalesImprovements(1:3:length(ScalesImprovements))], 'color', ChronicityIndex, 'subset',ismember(conditions,RGSStudies));
g(1,1).set_names('x','Condition ID','y','d Fugl-Meyer','color','# Chronicity');
g(1,1).set_title('Recovery from end of treatment (T1) to Follow-up');
g(1,1).stat_violin('half',false,'normalization','width','dodge',0,'fill','transparent');
g(1,1).stat_boxplot('width',1);
g(1,1).axe_property('xticklabels',stringConditions); %Should work for recent Matlab versions
g(1,1).set_order_options('x',0, 'color',0);
g(1,1).set_color_options('map','brewer2');

ScalesImprovements= [dataStructure.FollowUpImprovement];
g(2,1) = gramm('x', {dataStructure.NamesConditions}','y', [ScalesImprovements(2:3:length(ScalesImprovements))], 'color', ChronicityIndex, 'subset',ismember(conditions,RGSStudies));
g(2,1).set_names('x','Condition ID','y','d Cahai','color','# Chronicity');
g(2,1).stat_violin('half',false,'normalization','width','dodge',0,'fill','transparent');
g(2,1).stat_boxplot('width',1);
g(2,1).axe_property('xticklabels',stringConditions); %Should work for recent Matlab versions
g(2,1).set_order_options('x',0, 'color',0);
g(2,1).set_color_options('map','brewer2');
g.draw();
saveas(gcf,'figures/Figure_S2.eps','epsc')

% Controls T1
stringConditionsC= {'1' '4' '5' '6' '8' '11'};

figure ('position', [100 100 1000 370]);
ScalesImprovements= [dataStructure.TreatmentImprovement];
clear g
g(1,1) = gramm('x', {dataStructure.NamesConditions}','y', [ScalesImprovements(1:3:length(ScalesImprovements))], 'color', ChronicityIndex, 'subset',ismember(conditions,ControlStudies));
g(1,1).set_names('x','Condition ID','y','d Fugl-Meyer','color','# Chronicity');
g(1,1).set_title('Recovery from baseline to end of treatment (T1)');
g(1,1).stat_violin('half',false,'normalization','width','dodge',0,'fill','transparent');
g(1,1).stat_boxplot('width',1);
%g(1,1).no_legend();
g(1,1).axe_property('xticklabels',stringConditionsC); %Should work for recent Matlab versions
g(1,1).set_order_options('x',0, 'color',0);
g(1,1).set_color_options('map','brewer2');

ScalesImprovements= [dataStructure.TreatmentImprovement];
g(2,1) = gramm('x', {dataStructure.NamesConditions}','y', [ScalesImprovements(2:3:length(ScalesImprovements))], 'color', ChronicityIndex, 'subset',ismember(conditions,ControlStudies));
g(2,1).set_names('x','Condition ID','y','d Cahai','color','# Chronicity');
g(2,1).stat_violin('half',false,'normalization','width','dodge',0,'fill','transparent');
g(2,1).stat_boxplot('width',1);
g(2,1).axe_property('xticklabels',stringConditionsC); %Should work for recent Matlab versions
g(2,1).set_order_options('x',0, 'color',0);
g(2,1).set_color_options('map','brewer2');

% Controls Follow-Up
ScalesImprovements= [dataStructure.FollowUpImprovement];
g(1,2) = gramm('x', {dataStructure.NamesConditions}','y', [ScalesImprovements(1:3:length(ScalesImprovements))], 'color', ChronicityIndex, 'subset',ismember(conditions,ControlStudies));
g(1,2).set_names('x','Condition ID','y','d Fugl-Meyer','color','# Chronicity');
g(1,2).set_title('Recovery from end of treatment (T1) to Follow-up');
g(1,2).stat_violin('half',false,'normalization','width','dodge',0,'fill','transparent');
g(1,2).stat_boxplot('width',1);
g(1,2).set_order_options('x',0, 'color',0);
g(1,2).set_color_options('map','brewer2');
g(1,2).axe_property('xticklabels',stringConditionsC); %Should work for recent Matlab versions

ScalesImprovements= [dataStructure.FollowUpImprovement];
g(2,2) = gramm('x', {dataStructure.NamesConditions}','y', [ScalesImprovements(2:3:length(ScalesImprovements))], 'color', ChronicityIndex, 'subset',ismember(conditions,ControlStudies));
g(2,2).set_names('x','Condition ID','y','d Cahai','color','# Chronicity');
g(2,2).stat_violin('half',false,'normalization','width','dodge',0,'fill','transparent');
g(2,2).stat_boxplot('width',1);
g(2,2).set_order_options('x',0, 'color',0);
g(2,2).set_color_options('map','brewer2');
g(2,2).axe_property('xticklabels',stringConditionsC); %Should work for recent Matlab versions
g.draw();
saveas(gcf,'figures/Figure_S3.eps', 'epsc')

%% Correlations control
age_aux=[dataDemographics.age];
chron_aux = (reshape([dataStructure(:).ChronicityDays], 1, length([dataStructure(:).ChronicityDays]))');
[r_ageChron p_ageChron] = corr(age_aux(~isnan(age_aux)), chron_aux(~isnan(age_aux)), 'type', 'spearman');

% exctract baselines and correlate with chronicity for control
for ii=1:totalPatients
    scores = cell2mat(data(ii,6));
    baselines(ii,1:3) = [scores(1,1), scores(1,2), scores(1,3)];
end
[r_FMChron p_FMChron] = corr(baselines(~isnan(baselines(:,1)), 1), chron_aux(~isnan(baselines(:,1))), 'type', 'spearman');
[r_cahaiChron p_CahaiChron] = corr(baselines(~isnan(baselines(:,2)), 1), chron_aux(~isnan(baselines(:,2))), 'type', 'spearman');

%% Correlations of FM with CAHAI and Barthel at different chronicity quartiles:
FM_improvement_all = [ScalesImprovements(1:3:length(ScalesImprovements))];
CAHAI_improvement_all = [ScalesImprovements(2:3:length(ScalesImprovements))];
BARTHEL_improvement_all = [ScalesImprovements(3:3:length(ScalesImprovements))];
[A,sortIdx] = sort(chron_aux,'ascend');
FM_improvement_all = FM_improvement_all(sortIdx);
CAHAI_improvement_all = CAHAI_improvement_all(sortIdx);
BARTHEL_improvement_all = BARTHEL_improvement_all(sortIdx);

% Correlations of FM and CAHAI across chronicity 
Allr_FC=[];
Allr_FB=[];
Allp_FC=[];
Allp_FB=[];
window=length(FM_improvement_all)/4;
for i=1:window:length(FM_improvement_all)
    % We cannot assume symmetry in any of the datasets (Kolmogorov)
    [rhoC pvalC] = corr(FM_improvement_all(i:i+window)', CAHAI_improvement_all(i:i+window)', 'type', 'spearman')
    [Allr_FC] = cat (1, Allr_FC, rhoC)
    [Allp_FC] = cat (1, Allp_FC, pvalC)
    [rhoB pvalB] = corr(FM_improvement_all(i:i+window)', BARTHEL_improvement_all(i:i+window)', 'type', 'spearman')
    [Allr_FB] = cat (1, Allr_FB, rhoB)
    [Allp_FB] = cat (1, Allp_FB, pvalB)
end

% plot correlations
figure; 
subplot(2,1,1);
plot(Allr_FC, '-or','LineWidth',1, 'MarkerFaceColor','r');
hold on; plot(Allr_FB, '-^b','LineWidth',1, 'MarkerFaceColor','b');
ylabel('Spearman coefficient (ρ)') 
axis([0.5 4.5 0 0.8])
xticks([1 2 3 4])
xticklabels({'','','',''})
title('Correlation of FM-UE improvements with CAHAI and Barthel improvements')
legend('CAHAI','Barthel')
subplot(2,1,2);
plot(Allp_FC, '-or','LineWidth',1, 'MarkerFaceColor','r');
hold on; plot(Allp_FB, '-^b','LineWidth',1, 'MarkerFaceColor','b');
plot([0.5 4.5], ones(1,2)*0.05, '-.k','LineWidth',1);
plot([0.5 4.5], ones(1,2)*0.01, '-.k', 'LineWidth',1);
set(gca, 'YScale', 'log')
xticks([1 2 3 4])
xticklabels({'Q1','Q2','Q3','Q4'})
xlabel('Quartiles') 
ylabel('p-value') 
axis([0.5 4.5 0 1])
saveas(gcf,'figures/Figure3.eps','epsc')
