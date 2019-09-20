function ProcessData_Figs_4D_4E
% 
% Figures 4D and 4E, Behavioral responses to fast and slow decrease in odor concentration 
% 
% September 2019, Sagi Levy
%  
%%  load data and initialization
%%% load 'File' structures. Store in All_Files
NewPath = 'G:\';
load([NewPath,'Training_SB1_Pulses_Setup2_20170501\Training_SB1_Pulses_Setup2_20170501_DataMatrices_AllDelaysCorrected1.mat'],'File'); 
OldPath      = File.BackgroundFile(1:3);
File         = CorrectPathForFile_inline (File, OldPath, NewPath);
All_Files{1} = File;
load([NewPath,'Training_SB1_Pulses_Setup2_20170504\Training_SB1_Pulses_Setup2_20170504_DataMatrices_AllDelaysCorrected1.mat'],'File'); 
OldPath      = File.BackgroundFile(1:3);
File         = CorrectPathForFile_inline (File, OldPath, NewPath);
All_Files{2} = File;
load([NewPath,'Training_SB1_Pulses_Setup2_20170504\Training_SB1_Pulses_Setup2_20170504_DyePatterns'],'FlowDelay','PatternsOnset_Fast')

%%% Animals placement in microfluidic device arenas, per experiment
ConditionsTitles      = {'Naive,fast' , 'AfterTreatment,fast', 'Naive,slow','AfterTreatment,slow'};
ConditionsPerArena{1} = [2 1 3 4]; % Animals placement in experiment #1
ConditionsPerArena{2} = [1 2 4 3]; % Animals placement in experiment #2

%%% load 'Data' structures. Store in All_Data
for exp_ind = 1:2
    File = All_Files{exp_ind};      
    for arena = 1:File.NumArenas
        ConditionsIndex                   = ConditionsPerArena{exp_ind}(arena);        
        X                                 = load(File.FileNames.DataMatrices_AllDelaysCorrected{arena},'Data'); 
        All_Data{exp_ind,ConditionsIndex} = X.Data;
    end
end

%%% Initialization and assignment to ExpInfo
StimuliSequence    = [2 1; 2 4; 2 3];
StimuliTitles      = {'0','5*10^-^8','3*10^-^7'};                                           % Corresponds to columns 1:4 in the data matrices below 
TimeWindowInSec    = [-30 60];
TimeWindowInFrames = TimeWindowInSec * File.FrameRate;
TimeVecInSec       = (TimeWindowInFrames(1):TimeWindowInFrames(2))/File.FrameRate;          % Corresponds to X-axis in the data matrices below 

ExpInfo.FrameRate          = File.FrameRate; 
ExpInfo.StimuliTitles      = StimuliTitles; 
ExpInfo.ConditionsTitles   = ConditionsTitles; 
ExpInfo.TimeVecInSec       = TimeVecInSec; 
ExpInfo.FlowDelay          = FlowDelay; 
ExpInfo.PatternsOnset_Fast = PatternsOnset_Fast; 

%%  Process data
%%% Compute change in aversive behavior probability
AnalysisMode = 1;    % Probability of any aversive behavior 
[ResponseDynamics, ~, ~, AfterMinusBeforeStimulus] = ComputeStimulusResponses (All_Files, All_Data, TimeWindowInSec, StimuliSequence, AnalysisMode); 
ChangeInAversiveBehaviorProbability                = AfterMinusBeforeStimulus;
ChangeInAversiveBehaviorProbability_Dynamics       = ResponseDynamics;

%%% Compute change in speed
AnalysisMode = 2;   % Animal Speed
[ResponseDynamics, ~, ~, ~, DeltaBOverB] = ComputeStimulusResponses (All_Files, All_Data, TimeWindowInSec, StimuliSequence, AnalysisMode); 
ChangeInSpeed          = DeltaBOverB;
ChangeInSpeed_Dynamics = ResponseDynamics;

%% Save behavioral responses
save('BehavioralResponses_Figs_4D_4E.mat','ExpInfo',...
                                          'ChangeInAversiveBehaviorProbability','ChangeInAversiveBehaviorProbability_Dynamics',...
                                          'ChangeInSpeed','ChangeInSpeed_Dynamics','-v7.3');
                                      
return

function [AverageResponse, Responses] = CalculateResponse(Vec, positions, ElapsedTimeInFrames, StimuliSequence,TimeWindowInFrames)

FramesWindow    = TimeWindowInFrames(1):TimeWindowInFrames(2);
AverageResponse = zeros(size(StimuliSequence,1), length(FramesWindow))*NaN;
Responses = cell(1,size(StimuliSequence,1));
for s_ind = 1:size(StimuliSequence,1)
    previous_position = StimuliSequence(s_ind,1);
    next_position     = StimuliSequence(s_ind,2);
    
    positions_eq_previous = find(positions==previous_position);
    positions_eq_next     = find(positions==next_position);
    relevant_indices      = intersect(positions_eq_next, positions_eq_previous+1);
    EventsTiming          = ElapsedTimeInFrames(relevant_indices-1);
    
    if isempty(EventsTiming)
        continue
    end
    Response = zeros(length(EventsTiming), length(FramesWindow))*NaN;
    for e_ind = 1:length(EventsTiming)
        CurrentFrames = EventsTiming(e_ind)+ FramesWindow;
        if CurrentFrames(end)>length(Vec)
            Response = Response(1:end-1,:);
            break
        end
        Response(e_ind,:) = Vec(CurrentFrames);                
    end    
    AverageResponse(s_ind,:) = mean(Response,1);  
    Responses{s_ind} = Response;
end

return

function BehaviorCode
%  My Codes:  
%    Low Level
%             0 'Out of bounds'
%             1 'Forward'
%             2 'Curve'
%             3 'Pause'
%             4 'Reverse'
%             5 'Omega'
%             6 'SharpTurn'
%
%    High Level
%             0 'Out of bounds'
%             1 'Forward'
%             2 'Curve'
%             3 'Pause'
%             4 'Reversal'
%             5 'Omega-Pause'
%             6 'SharpTurn'
%             7 'Pirouette. Initial reversal'
%             8 'Pirouette. After reversal'



return

function [AllRawValues, BeforeStimulus, AfterStimulus, AfterMinusBeforeStimulus, DeltaBOverB] = ...
                      ComputeStimulusResponses (All_Files, All_Data, TimeWindowInSec, StimuliSequence, AnalysisMode) 
                  
%% Notes:
%  Choose AnalysisMode==1 for probability of aversive behavior
%  Choose AnalysisMode==2 for animal speed

%% Initialization
% Correct delay (between camera and odor valves) 
FrameRate         = All_Files{1}.FrameRate;    % Hz 
ManualTimeDelay   = -4*FrameRate; % Delay between camera and valves initiation 
% Filter
smoothvec         = true;
AverageWindowSize = 5*FrameRate+1; % 5.16 seconds
b     = ones(1,AverageWindowSize)/AverageWindowSize;
a     = 1;
delay = round((AverageWindowSize-1)/2);
% Time
TimeWindowInFrames = TimeWindowInSec * FrameRate;
TimeVecInSec       = (TimeWindowInFrames(1):TimeWindowInFrames(2))/FrameRate;
ZeroIndex          = find(TimeVecInSec==0,1);
% Structure initialization
LEN                             = size(StimuliSequence,1);
AllRawValues.Matrices           = cell(All_Files{1}.NumArenas,LEN);
BeforeStimulus.MeanValues       = cell(All_Files{1}.NumArenas,LEN);
AfterStimulus.MeanValues        = cell(All_Files{1}.NumArenas,LEN);
TotalNumberOfWorms_PerCondition = ones(All_Files{1}.NumArenas,1);

%% Extract responses to change in odor concentration for each condition and concentration level
for exp_ind = 1:2
    File = All_Files{exp_ind};      
    for condition_ind = 1:File.NumArenas
        %% Extract relevant data
        Data = All_Data{exp_ind,condition_ind}; 

        TimeInMinutes          = (1:size(Data.Coordinates_X,2))/File.FrameRate/60;
        positions              = File.AttractionTests.positions_valve1;
        positions(File.AttractionTests.positions_FastValve==1) = File.AttractionTests.positions_valve2(File.AttractionTests.positions_FastValve==1);
        ElapsedTimeInMinutes   = cumsum(File.AttractionTests.StimulusTimes);
        RelevantStimuli        = ElapsedTimeInMinutes<=TimeInMinutes(end);
        positions              = positions(RelevantStimuli);
        ElapsedTimeInMinutes   = ElapsedTimeInMinutes(RelevantStimuli);
        ElapsedTimeInFrames    = ElapsedTimeInMinutes*60*File.FrameRate;        
        NumOfWorms             = size(Data.Speed,1);
        TotalNumberOfWorms_PerCondition(condition_ind) = TotalNumberOfWorms_PerCondition(condition_ind) + NumOfWorms;
        if AnalysisMode==1
            Vec                 = sum(Data.BehaviorProbability_LowLevel(3:end,:),1); Ylimits = [0 1];  % all aversive behavior, without curving
        elseif AnalysisMode==2
            Speed = Data.Speed; Speed(Data.BehaviorProbability_LowLevel>2)=NaN;  Vec = nanmean(Speed,1); Ylimits = [0 350];        
        end
        % Smooth and bias corrections
        if smoothvec
            Vec = filter(b,a,Vec);
            Vec(1:end-delay) = Vec(delay+1:end); 
        end
        if ManualTimeDelay>0
            Vec(1:end-ManualTimeDelay) = Vec(ManualTimeDelay+1:end); 
        elseif ManualTimeDelay<0
            Vec(abs(ManualTimeDelay)+1:end) = Vec(1:end-abs(ManualTimeDelay)); 
        end        
        [~,Responses] = CalculateResponse(Vec, positions, ElapsedTimeInFrames, StimuliSequence,TimeWindowInFrames);
        
        %% Compare data before vs. after stimulus initiation 
        for pulse_ind=1:LEN
            if isempty(Responses{pulse_ind})
                continue
            end
            % Before stimulus
            indices = (ZeroIndex-10*File.FrameRate):ZeroIndex;  % [-10,0] seconds
            MeanValues = mean(Responses{pulse_ind}(:,indices),2);
            BeforeStimulus.MeanValues{condition_ind,pulse_ind} = [BeforeStimulus.MeanValues{condition_ind,pulse_ind}; MeanValues];
            % After stimulus
            indices = ZeroIndex:(ZeroIndex+50*File.FrameRate);  % [0,50] seconds
            MeanValues = mean(Responses{pulse_ind}(:,indices),2);
            AfterStimulus.MeanValues{condition_ind,pulse_ind} = [AfterStimulus.MeanValues{condition_ind,pulse_ind}; MeanValues]; 
            % Raw values
            AllRawValues.Matrices{condition_ind,pulse_ind}    = [AllRawValues.Matrices{condition_ind,pulse_ind}; Responses{pulse_ind}]; 
        end        
    end
end
TotalNumberOfWorms_PerCondition'

%% Compute response change (e.g. after vs. before stimulus) and statistics over multiple repeats (responses to independent pulses of similar concentration change) 
AllRawValues.MeanDynamics = zeros(File.NumArenas,pulse_ind,length(TimeVecInSec),'single')*NaN;
AllRawValues.STDDynamics  = zeros(File.NumArenas,pulse_ind,length(TimeVecInSec),'single')*NaN;
AllRawValues.SEMDynamics  = zeros(File.NumArenas,pulse_ind,length(TimeVecInSec),'single')*NaN;
AllRawValues.NumOfRepeats = zeros(File.NumArenas,pulse_ind,'single')*NaN;
BeforeStimulus.MEANS      = zeros(File.NumArenas,LEN)*NaN;
BeforeStimulus.STDs       = zeros(File.NumArenas,LEN)*NaN;
BeforeStimulus.SEMs       = zeros(File.NumArenas,LEN)*NaN;
AfterStimulus.MEANS       = zeros(File.NumArenas,LEN)*NaN;
AfterStimulus.STDs        = zeros(File.NumArenas,LEN)*NaN;
AfterStimulus.SEMs        = zeros(File.NumArenas,LEN)*NaN;
AfterMinusBeforeStimulus  = BeforeStimulus;
DeltaBOverB               = BeforeStimulus;
for condition_ind = 1:File.NumArenas
    for pulse_ind=1:LEN
        % Before
        MeanValues                                    = BeforeStimulus.MeanValues{condition_ind,pulse_ind};
        BeforeStimulus.MEANS(condition_ind,pulse_ind) = nanmean(MeanValues);
        BeforeStimulus.STDs(condition_ind,pulse_ind)  = nanstd(MeanValues,1);
        BeforeStimulus.SEMs(condition_ind,pulse_ind)  = BeforeStimulus.STDs(condition_ind,pulse_ind) / sqrt(length(MeanValues));
        % After
        MeanValues                                      = AfterStimulus.MeanValues{condition_ind,pulse_ind};
        AfterStimulus.MEANS(condition_ind,pulse_ind)  = nanmean(MeanValues);
        AfterStimulus.STDs(condition_ind,pulse_ind)   = nanstd(MeanValues,1);
        AfterStimulus.SEMs(condition_ind,pulse_ind)   = AfterStimulus.STDs(condition_ind,pulse_ind) / sqrt(length(MeanValues));
        % After-Before
        MeanValues                                      = AfterStimulus.MeanValues{condition_ind,pulse_ind} - BeforeStimulus.MeanValues{condition_ind,pulse_ind};
        AfterMinusBeforeStimulus.MeanValues{condition_ind,pulse_ind} = MeanValues;
        AfterMinusBeforeStimulus.MEANS(condition_ind,pulse_ind)      = nanmean(MeanValues);
        AfterMinusBeforeStimulus.STDs(condition_ind,pulse_ind)       = nanstd(MeanValues,1);
        AfterMinusBeforeStimulus.SEMs(condition_ind,pulse_ind)       = AfterMinusBeforeStimulus.STDs(condition_ind,pulse_ind) / sqrt(length(MeanValues));            
        % DeltaBOverB = (After-Before)/Before
        MeanValues =  (AfterStimulus.MeanValues{condition_ind,pulse_ind} - BeforeStimulus.MeanValues{condition_ind,pulse_ind})./ ...
                       BeforeStimulus.MeanValues{condition_ind,pulse_ind} * 100;
        DeltaBOverB.MeanValues{condition_ind,pulse_ind} = MeanValues;
        DeltaBOverB.MEANS(condition_ind,pulse_ind)      = nanmean(MeanValues);
        DeltaBOverB.STDs(condition_ind,pulse_ind)       = nanstd(MeanValues,1);
        DeltaBOverB.SEMs(condition_ind,pulse_ind)       = DeltaBOverB.STDs(condition_ind,pulse_ind) / sqrt(length(MeanValues));    
        % Raw values
        Mat                                                    = AllRawValues.Matrices{condition_ind,pulse_ind};
        AllRawValues.MeanDynamics(condition_ind,pulse_ind,:) = nanmean(Mat,1);
        AllRawValues.STDDynamics(condition_ind,pulse_ind,:)  = nanstd(Mat,1,1);
        NumOfRepeats                                           = size(Mat,1);
        AllRawValues.SEMDynamics(condition_ind,pulse_ind,:)  = nanstd(Mat,1,1)/ sqrt(NumOfRepeats);
        AllRawValues.NumOfRepeats(condition_ind,pulse_ind)   = NumOfRepeats; % Number of pulses per condition (rows) and concentration (columns)   
    end
end

return

function File = CorrectPathForFile_inline (File, OldPath, NewPath)
% example
%  OldPath = 'K:\'
%  NewPath = 'M:\'

FIELDNAMES = fieldnames(File);
for f_ind = 1:length(FIELDNAMES)
    FieldName = FIELDNAMES{f_ind};    
    File.(FieldName) = ChangePathForField (File.(FieldName), OldPath, NewPath);    
end

return

function FieldAfterChange = ChangePathForField (FieldBeforeChange, OldPath, NewPath)
X                = FieldBeforeChange;
FieldAfterChange = FieldBeforeChange;

LengthOldPath = length(OldPath);

ID_X  = whos('X');
classX= ID_X.class;

if strcmpi(classX,'char')        % Possible existence of the path string in FIRST layer fields
    if strfind(X,OldPath)
        FieldAfterChange = [NewPath, X((LengthOldPath+1):end)];
    end

elseif strcmpi(classX,'cell')  
    X     = FieldAfterChange{1};
    ID_X  = whos('X');
    classX= ID_X.class;
    if strcmpi(classX,'char')   % Possible existence of the path string in a cell array at the FIRST field 
        for cell_ind = 1:length(FieldAfterChange);
            X = FieldAfterChange{cell_ind};
            if strfind(X,OldPath)
                X_New = [NewPath, X((LengthOldPath+1):end)];
                FieldAfterChange{cell_ind} = X_New;
            end
        end                 
    end
    
elseif strcmpi(classX,'struct')&& length(X)==1  % Possible existence of the path string in SECOND layer fields. NOT CHECKING STRUCTURE ARRAYS !!!
    SECOND_FIELDNAMES = fieldnames(FieldBeforeChange);
    for f_ind = 1:length(SECOND_FIELDNAMES);
        SecondFieldName                    = SECOND_FIELDNAMES{f_ind};    
        FieldAfterChange.(SecondFieldName) = ChangePathForField (FieldBeforeChange.(SecondFieldName), OldPath, NewPath);   
    end        
end

return





