%function [data] = importAll(filename)
%IMPORTALL Reads and structures data from:
%
%   Workbook: Clinical_Scales_all.csv
%
%   For analysis in:
%
%   "The impact of VR-based rehabilitation on 
%   post-stroke functional recovery: a retrospective meta-analysis"
% 
%   IMPORTALL returns a 219x17 cell containing the following columns:
%   'patientID' 
%   'Chronicity' 
%   'Treatment' 
%   'ChronicityDays' 
%   'EvalTime' 
%   'Scale' 
%   'TreatmentFlag' 
%   'ImprovementRate' 
%   'ImprovementTime' 
%   'PatientIDs' 
%   'StudyId' 
%   'StudyIDs' 
%   'AllImprovement' 
%   'TreatmentImprovement' 
%   'FollowUpImprovement', 
%   'OrderingIndex', 
%   'NamesConditions'
% 
% Author Armin Duff armin.duff@gmail.com 
% and Belen Rubio Ballester belen.rubio.ballester@gmail.com 2018
% ---------------------------------------------------------------- %

function [data] = importAll(filename)

    %% Initialize variables.
    delimiter = ',';

    %% Read columns of data as strings:
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%s%s%s%s%s%s%s%s%[^\n\r]';

    %% Open the text file2
    fileID = fopen(filename,'r');

    %% Read columns of data according to format string.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

    %% Close the text file.
    fclose(fileID);

    %% Convert the contents of columns containing numeric strings to numbers.
    % Replace non-numeric strings with NaN.
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));

    for col=[1,4,6,7,8]
        % Converts strings in the input cell array to numbers. Replaced non-numeric
        % strings with NaN.
        rawData = dataArray{col};
        for row=1:size(rawData, 1);
            % Create a regular expression to detect and remove non-numeric prefixes and
            % suffixes.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData{row}, regexstr, 'names');
                numbers = result.numbers;

                % Detected commas in non-thousand locations.
                invalidThousandsSeparator = false;
                if any(numbers==',');
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(thousandsRegExp, ',', 'once'));
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % Convert numeric strings to numbers.
                if ~invalidThousandsSeparator;
                    numbers = textscan(strrep(numbers, ',', ''), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch me
            end
        end
    end

    %% Split data into numeric and cell columns.
    rawNumericColumns = raw(:, [1,4,6,7,8]);
    rawCellColumns = raw(:, [2,3,5]);


    %% Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
    rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

    %% Allocate imported array to column variable names
    patientID = cell2mat(rawNumericColumns(2:end, 1));
    treatment = rawCellColumns(2:end, 1);
    study = rawCellColumns(2:end, 2);
    daysAS = cell2mat(rawNumericColumns(2:end, 2));
    eval_time = rawCellColumns(2:end, 3);
    eval_time(strcmp(eval_time, 'baseline'))={0};
    eval_time(strcmp(eval_time, 'w01'))={7};
    eval_time(strcmp(eval_time, 'w02'))={14};
    eval_time(strcmp(eval_time, 'w03'))={21};
    eval_time(strcmp(eval_time, 'w05'))={35};
    eval_time(strcmp(eval_time, 'w06'))={45};
    eval_time(strcmp(eval_time, 'w12'))={84};
    eval_time(strcmp(eval_time, 'm06'))={183};
    eval_time(strcmp(eval_time, 'm12'))={365};
    evalT=cell2mat(eval_time);
    clear scale
    scale(:,1) = cell2mat(rawNumericColumns(2:end, 3));
    scale(:,2) = cell2mat(rawNumericColumns(2:end, 4));
    scale(:,3) = cell2mat(rawNumericColumns(2:end, 5));


    index=0;
    pID=0;
    uniPatientID=unique(patientID);
    uniPatientID=uniPatientID(histc(patientID(:),unique(patientID))==1); %patients with just one measurement
    dataLength=length(unique(patientID))-length(uniPatientID); %prepare data structure
    data=cell(dataLength,7);


    for ii = 1:length(patientID)
        if sum(patientID(ii)==uniPatientID)==0 %ignore the cases with just one measurement
            if(patientID(ii)>pID)
                index=index+1;
                pID=patientID(ii); %add patient data
                data{index,1}=patientID(ii);
                data{index,2}=study(ii);
                data{index,3}=treatment(ii);
                data{index,4}=daysAS(ii);
            end
            data{index,5}=[data{index,5};evalT(ii)]; %append evalT
            data{index,6}=[data{index,6};scale(ii,:)]; %append scale

            if(evalT(ii)>3) %default = 0
                rgs=0;
            else
                rgs=[]; %ignore baseline
            end
            
            if ~isempty(cell2mat(regexp(treatment(ii),'c|C')))
                if evalT(ii)>3 && evalT(ii)<40 %for w1 w2 w3 and w5 
                    rgs=2;
                end
                if evalT(ii)==84&& ~isempty(cell2mat(regexp(treatment(ii),'c12w')))
                    rgs=2; %for w12 in rgs12w
                end
                if evalT(ii)==45&& ~isempty(cell2mat(regexp(treatment(ii),'AM')))
                    rgs=2; %for w12 in rgs12w

                end
            else
            
                if ~isempty(cell2mat(regexp(treatment(ii),'RGS|rgs')))
                    if evalT(ii)>3 && evalT(ii)<40 %for w1 w2 w3 and w5 
                        rgs=1;
                    end
                    if evalT(ii)==84&& ~isempty(cell2mat(regexp(treatment(ii),'rgs12w')))
                        rgs=1; %for w12 in rgs12w
                    end
                end
                if evalT(ii)==45&& ~isempty(cell2mat(regexp(treatment(ii),'AM')))
                    rgs=1; %for w12 in rgs12w

                end
            end
            
            data{index,7}=[data{index,7};rgs]; %append rgs treatment flag
        end
    end

    %sorting Data 
    for ii =1:length(data)
        [data{ii,5},sortIndex]=sort(data{ii,5});
        data{ii,6}=data{ii,6}(sortIndex,:);
        data{ii,7}=data{ii,7}(sortIndex(sortIndex>1)-1);
    end


    %% Calc Improvement
    data(:,8)=cellfun(@(x,y) diff(x)./repmat(diff(y),1,3), data(:,6),data(:,5),'UniformOutput', false) ;%improvement rate (diff scale/time)
    data(:,9)=cellfun(@(x,y) x+y(1:end-1)+diff(y)/2,data(:,4),data(:,5),'UniformOutput', false); %time of improvment

    %% Duplicate Patient id and study
    studyLabels = {'acuteOld','acute','subacute','chronicOld','chronic','atHome', 'chronicAM'};


    nData=num2cell(cellfun(@length, data(:,7))); %size
    data(:,10)=cellfun(@(x,y)repmat(x,y,1),data(:,1),nData,'UniformOutput',false); %patient ID
    data(:,11)=num2cell(cellfun(@(x)find(strcmp(x,studyLabels)),data(:,2)));
    data(:,12)=cellfun(@(x,y)repmat(x,y,1),data(:,11),nData,'UniformOutput',false); %study
    
    %Identify group per study (i.e. according to intervention and dosage)
    c=[data{:,2}; data{:,3}]';
    groupsNames = strcat(c(:,1), c(:,2));
    [UniqueGroups, ~, ib_group] = unique(groupsNames);
    [meanChronConditions] = grpstats([data{:,4}],ib_group, @median);
    for(cat=1: length(UniqueGroups))
        aux(find(ismember(ib_group,cat))) = meanChronConditions(cat);
    end
    
    %Calc absolute improvement
    data(:,13)=cellfun(@(x) diff(x), data(:,6),'UniformOutput', false); %improvement 
    
    %Compute absolute improvement during RGS treatment and follow-up
    for i = 1 : numel(data(:,13))
    % Apply mask on each of the matrices and store the result back in the cell.
        [r c] = size(cell2mat(data(i,13)));
        if (r>1)
            data(i,14) = mat2cell(nansum(cell2mat(data(i,13)).* (cell2mat(repmat(data(i,7),1,3))~=0)), 1, 3); %RGS
            data(i,15) = mat2cell(nansum(cell2mat(data(i,13)).* ((cell2mat(repmat(data(i,7),1,3)))==0) ), 1, 3); %RGS follow-up
        else
            data(i,14) = mat2cell((cell2mat(data(i,13)).* (cell2mat(repmat(data(i,7),1,3))~=0)), 1, 3);
  
        end
    end
    
    data(:,16)=num2cell(aux); %Add chornicity by group
    data(:,17)= groupsNames;
    
    %% Clear temporary variables
    clearvars filename delimiter eval_time formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;

end
