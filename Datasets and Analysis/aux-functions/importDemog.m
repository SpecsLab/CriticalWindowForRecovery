% ------------------- Main Script -------------------------------- %
%
% Script for importing data from the following spreadsheet:
%
%    Workbook: patients_Demographics_and_ClinicalScreening.xlsx
%
% Reported in:
%
% "The impact of VR-based rehabilitation on 
% post-stroke functional recovery: a retrospective meta-analysis"
% 
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Author Armin Duff armin.duff@gmail.com 
% and Belen Rubio belen.rubio.ballester@gmail.com 2018
% ---------------------------------------------------------------- %

function [statarray, patientsDemographicsandClinicalScreening] = importDemog(filename, idPat, categoryID)


%% Import the data
[~, ~, raw] = xlsread('data/patients_Demographics_and_ClinicalScreening.xlsx');
raw = raw(2:end,1:11);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[2,3,6,10]);
raw = raw(:,[1,4,5,7,8,9,11]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
patientsDemographicsandClinicalScreening = table;

%% Allocate imported array to column variable names
aux = data(~isnan(data(:,1)),1);
IncludedPatients = find(ismember(aux, setdiff(aux,idPat, 'sorted'))==0);
    

patientsDemographicsandClinicalScreening.patientID = data(IncludedPatients,1);
patientsDemographicsandClinicalScreening.category = cellVectors(IncludedPatients,1);
patientsDemographicsandClinicalScreening.hospital = cellVectors(IncludedPatients,2);
patientsDemographicsandClinicalScreening.gender = data(IncludedPatients,2);
patientsDemographicsandClinicalScreening.age = data(IncludedPatients,3);
patientsDemographicsandClinicalScreening.stroke_type = cellVectors(IncludedPatients,3);
patientsDemographicsandClinicalScreening.oxford_classification = data(IncludedPatients,4);
patientsDemographicsandClinicalScreening.affected_arm = data(IncludedPatients,5);
patientsDemographicsandClinicalScreening.dominant_hand = data(IncludedPatients,6);
patientsDemographicsandClinicalScreening.aphasia = cellVectors(IncludedPatients,4);
patientsDemographicsandClinicalScreening.days_after_stroke = data(IncludedPatients,7);
%[a b] = unique(patientsDemographicsandClinicalScreening.patientID, 'stable');
%patientsDemographicsandClinicalScreening.groupIntervention = categoryID;
for(ii=1: length(idPat))
   aux2(find(ismember(patientsDemographicsandClinicalScreening.patientID,idPat(ii)))) = categoryID(ii);
end
patientsDemographicsandClinicalScreening.groupIntervention = [aux2]';

%% Clear temporary variables
clearvars data raw cellVectors R;

%% Parse Demographics
func = @(x) nansum(x==2)/length(x)*100;
statarray = grpstats(patientsDemographicsandClinicalScreening(:,[5 11 12]),'groupIntervention', {@median, 'range'} );
statarray(:,7) = table(splitapply(func,patientsDemographicsandClinicalScreening(:,[4]),patientsDemographicsandClinicalScreening.groupIntervention));
statarray.Properties.VariableNames{'Var7'} = 'Gender';
statarray(:,8) = table(splitapply(func,patientsDemographicsandClinicalScreening.affected_arm,patientsDemographicsandClinicalScreening.groupIntervention));
statarray.Properties.VariableNames{'Var8'} = 'PareticArm';

%% Parse OxfordClass
Oxf=NaN(17,5);
for(cats=1:17)
    tbl = tabulate([patientsDemographicsandClinicalScreening{(aux2==cats),7}]);
    if(~isempty(tbl))
        Oxf(cats,tbl(:,1)) = tbl(:,2)';
    end
end
    