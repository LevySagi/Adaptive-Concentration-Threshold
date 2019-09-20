function CorrectFlowDelaysForPulseDevices_02(File, ExperimentInformationFunctionName)

% Examples:     CorrectFlowDelaysForPulseDevices_01(File,'Information_Exp1')
%               CorrectFlowDelaysForPulseDevices_01(File,'Information_Exp2')

%% Get experiment-specific information
eval(['AttractionTests = ',ExperimentInformationFunctionName,';'])
warning('off','images:initSize:adjustingMag');

%% Load FlowDelay data
File.FileNames.DyePatterns = [File.TrackFile(1:end-4),'_DyePatterns.mat'];    
if exist(File.FileNames.DyePatterns,'file')==2
    disp(['loading ''FlowDelay'' from ',File.FileNames.DyePatterns]);
    load (File.FileNames.DyePatterns, 'FlowDelay');
    AttractionTests.FlowDelay = FlowDelay;
else
    disp('''FlowDelay'' was not calculated from dye patterns. Aborting.');    
end
File.AttractionTests = AttractionTests;

%% Correct Flow Delay 
disp([datestr(now), '  Correcting flow delay for raw matrices']);
CorrectFlowDelay(File);

warning('on','images:initSize:adjustingMag');
return

%% Flow rate corrections to raw Data matrices
function CorrectFlowDelay(File)
AttractionTests      = File.AttractionTests;
Interpolation_Factor = 1e-4;   % For all-delay-corrected behavior probability matrices 

NumOfArenas = File.NumArenas;
for ar = 1:NumOfArenas
    % load
    load(File.FileNames.DataMatrices{ar},'Data','TracksStats','background');
    
    % normalize
    [Data_AllDelaysCorrected, Data_ArenaDelayCorrected] = NormalizeFlowDelay (Data , ar, File, AttractionTests, Interpolation_Factor, TracksStats.BehaviorCode) ;
    
    % save
    Data = Data_ArenaDelayCorrected;
    save(File.FileNames.DataMatrices_ArenaDelayCorrected{ar},'Data','TracksStats','File','background','-v7.3');
    Data = Data_AllDelaysCorrected;
    save(File.FileNames.DataMatrices_AllDelaysCorrected{ar},'Data','TracksStats','File','background','-v7.3');
end    
    
return

function [Data_AllDelaysCorrected, Data_ArenaDelayCorrected, DetectionAndSegmentationProbabilityVector]  = NormalizeFlowDelay (Data , ArenaID, File, AttractionTests,  Interpolation_Factor, BehaviorCode)  

%% Inputs:
% BehDataPerArena_in  rows = worms, columns= frames, data = behavior code
% FlowDelay           structure with delay time information (see below).
%                     FlowDelay.DelayBetweenValveAndArenas = a vector corresponding to arenas 1-4 (upper left, upper right, lower left, lower right).  
%                                                            Values represent the of between the valve switching time and the time that the new fluid STARTS to flow into the device. 
%                                                            START time is defined as the time that the dye level BEFORE the arena reaches 50% of its maximal level. 
% File                Structure with the movie parameters, including FrameRate, resolution and image size.   
% AttractionTests     Structure that contains the stimuli timing: 
%                     The stimuli timing will be used to 
% 
%% Output:
% 'BehDataPerArena'  structure is similar to the input structure: BehDataPerArena_in, with the following additional fields:
%   DelayBetweenValveAndArenas_InFrames =   a vector of length (NumberOfArenas). [Frames]. Each value corresponds to a different arena.   
%                                           Values are the per arena delay between the valve switch time and onset of flow into the arena.
%                                           For each arena, the corresponding value will be ADDED to ALL WORMS during ALL THE EXPERIMENT  
%   FrameDelayMatrix = A cell array of length (NumberOfArenas). 
%                      Optional parameter, which will be calculated only of 'AttractionTests' input is given. 
%                      Each cell corresponds to a different arena and contains a matrix of size (behdata)= size(xmat). [Frames]. 
%                      FrameDelayMatrix{ar}(w_ind, fr_ind) = In arena 'ar', the worm 'w_ind', in frame 'f_ind' feels the input stimulus that is given at frame: ( f_ind - FrameDelayMatrix{ar}(w_ind, fr_ind)); 
%
%% Function explanation
% This function normalizes the matrices in 'Data' with respect to the flow delay in each arena.  
% Each column in these matrices is mapped to the corresponding frame number in the movie.   
% However, there is a delay from the time that the valve changes its position to the time the flow reaches the worms. 
% The Column index is therefore not a good measure for the time that the worms experience a stimulus. 
%
% Two types of normalizations are needed: 
%
%   1. DelayBetweenValveAndArenas_InFrames field: Normalizing the frames in each arena relative to the time it takes for the stimulus to reach the BEGINNING of the arena.  
%      This is important to devices and has two delay sources:
%           a. Inter-arena delay: Delay between the time stimulus reaches the arenas near the inlets, to the time it reaches the arenas near the outlets.      
%           b. In devices for which the switching between buffers is upstream to the microfluidic devices, there is a delay due to the liquid flow in the tubing between the valve in device.    
%      Similar normalization is given to all worms in a given arena. Use the following field from 'FlowDelay.DelayBetweenValveAndArenas' vector and the movie frame rate 
% 
%   2. Normalizing the frames for each worm seperately relative to it's instantenous positions in the arena AT THE TIME WHEN THE STIMULI ARRIVE. 
%      i.e. Normalization 1 deals only with the delay until the arena entrance and not with the delay from the arena entrance to the actual worm location within the arena.  
%      I will use:
%           The location of each worm:          'Coordinates_X' (not Coordinates_Y !). Corresponds to ImageSize(2) (not 1!!)    
%           The flow rate in the device using:  'FlowDelay.DelayWithinArena_framesPER_pixel' that is calculated below.
%           Stimuli onset and offset times:     given by 
%      NOTE!!! The delay is calculated from the position of the worm AT THE TIME THAT THE STIMULUS REACH THE ARENA !! The worms displacement during the spread of stimulus through the arena is NEGLECTED.   
%      The main idea is: 
%               At the time of each stimulus --> the frames will EXACTLY fit the time delay in the device. this will be done BOTH for the ON and OFF stimuli changes.  (This can be revisited later !!).
%               Since there is a different delay time at each stimulus time (the worms have different positions when different stimuli reach the arena) then the time scale between stimuli need to be corrected
%               The time scale between stimuli will be either streched or shortened with EVEN SPACING. The frames values will then be ROUNDED in order to allow further averaging.  
%
%% Initialization
FlowDelay                = AttractionTests.FlowDelay;
FieldsToCorrect          = fieldnames(Data); 
NumOfWorms               = size(Data.Coordinates_X,1);
NumOfFrames              = size(Data.Coordinates_X,2); 
if NumOfWorms == 1
    FieldsOnlyForMultipleWorms = {'BehaviorProbability_LowLevel','BehaviorProbability_Smoothed_LowLevel','BehaviorProbability_HighLevel','BehaviorProbability_Smoothed_HighLevel'};
    NON_relevant_indices       = false(1,length(FieldsToCorrect));
    for f=1:length(FieldsOnlyForMultipleWorms)
        NON_relevant_indices = NON_relevant_indices | strcmp(FieldsToCorrect,FieldsOnlyForMultipleWorms{f})';
    end
    FieldsToCorrect = FieldsToCorrect(~NON_relevant_indices);
end
Data_ArenaDelayCorrected = Data;
Data_AllDelaysCorrected  = Data; 

Frame_per_sec                                  = File.FrameRate;      % Frames per second
Pixels_per_mm                                  = File.PixelSize;      % pixels per mm !!!!!
Pixels_per_micrometer                          = Pixels_per_mm /1e3;  % pixels per micrometer
FlowDelay.DelayBetweenValveAndArenas_InFrames  = FlowDelay.DelayBetweenValveAndArenas * Frame_per_sec;

%% Normalization 1: Assigning each arena with its delay time based on arena ID 
ArenaDelay = FlowDelay.DelayBetweenValveAndArenas_InFrames(ArenaID);  % Delay until stimulus reach the arena

for f_num = 1:length(FieldsToCorrect)       % Correct all matrices by shifting it 'ArenaDelay' frames BACKWARDS  
    field_name = FieldsToCorrect{f_num};
    Data_ArenaDelayCorrected.(field_name)(:, 1: end-ArenaDelay)        = Data.(field_name)(:, ArenaDelay+1: end);    
    MAT = Data_ArenaDelayCorrected.(field_name); X = whos('MAT');
    if strcmpi(X.class,'single')|| strcmpi(X.class,'double')
        if NumOfWorms > 1 
            Data_ArenaDelayCorrected.(field_name)(:, (end-ArenaDelay+1) : end) = NaN; 
        else
            Data_ArenaDelayCorrected.(field_name)((end-ArenaDelay+1) : end) = NaN; 
        end
    else
        if NumOfWorms > 1 
            Data_ArenaDelayCorrected.(field_name)(:, (end-ArenaDelay+1) : end) = 0;   
        else
            Data_ArenaDelayCorrected.(field_name)((end-ArenaDelay+1) : end) = 0; 
        end
    end
end
if NumOfWorms > 1 
    Data_ArenaDelayCorrected.NaN(:, (end-ArenaDelay+1) : end) = 1; 
else
    Data_ArenaDelayCorrected.NaN((end-ArenaDelay+1) : end) = 1;     
end

%% Normalization 2: Calculating the worms delay at relative to the stimuli entrance to the arena. Calculate the overall FrameDelayMatrix
% For each arena, generate a timing matrix: 'FrameDelayMatrix'.
% size(FrameDelayMatrix(ar)) = size(behmat) = rows are worms, columns are frames.  
% This includes both that delay of the given arena, and the effects of flow per arena.  
            
if exist('AttractionTests','var')    % Stimulus data must be available for normalization 2      

    %% Calculate flow rate in the device. It may be a bit different between left and right sides of the device.
    % For SB1 devices:
    %      INLETS
    %  -------------
    %  |     |     |
    %  |  1  |  2  |
    %  |     |     |
    %  -------------
    %  |     |     |
    %  |  3  |  4  |
    %  |     |     |
    %  -------------
    %     OUTLETS
    %
    %
    if (trcmpi(FlowDelay.DeviceType,'SB1') && isfield(FlowDelay,'FlowRate_micrometerPERsec_Arenas_1and3')  % flow rate is defined differentially for left and right sides        
        if (ArenaID==1)||(ArenaID==3)
            FlowRate_micrometerPERsec  = FlowDelay.FlowRate_micrometerPERsec_Arenas_1and3;
            InterArenaDelay            = FlowDelay.InterArenasDelay_Arenas_1and3;
        elseif (ArenaID==2)||(ArenaID==4)
            FlowRate_micrometerPERsec  = FlowDelay.FlowRate_micrometerPERsec_Arenas_2and4;
            InterArenaDelay            = FlowDelay.InterArenasDelay_Arenas_2and4;
        end
    else
        disp('not SB1 device. The device delay is assumed to be:  min(FlowDelay.DelayBetweenValveAndArenas_InFrames)')
        FlowRate_micrometerPERsec  = FlowDelay.FlowRate_micrometerPERsec;
        InterArenaDelay            = FlowDelay.InterArenasDelay;
    end   

    FlowDelay.FlowRate_pixelsPERframe              = FlowRate_micrometerPERsec * Pixels_per_micrometer / Frame_per_sec;
    FlowDelay.DelayWithinArena_framesPER_pixel     = 1/FlowDelay.FlowRate_pixelsPERframe ;

    %% Find frames in which the stimulus was modified 
    AttractionTests.EndTimeInMinutes   = cumsum( AttractionTests.StimulusTimes);
    AttractionTests.SwitchFrames       = AttractionTests.EndTimeInMinutes * 60 * Frame_per_sec; 
    AttractionTests.FirstFrame         = 1;  
    AttractionTests.LastFrame          = min([AttractionTests.SwitchFrames(end) File.NumberOfFrames]);   % Some switched may correspond to dye tests where worm tracks were not analyzed
    AttractionTests.SwitchFrames       = AttractionTests.SwitchFrames(AttractionTests.SwitchFrames < AttractionTests.LastFrame);      

    %% Find the time when the stimuli reach the arena entrance
    FrameIndexAtSwitchTimes    = AttractionTests.SwitchFrames + ArenaDelay;         % ADDING THE ArenaDelay
    % Avoid possible switches that reach the arena only after the experiment ends:               
    FrameIndexAtSwitchTimes    = FrameIndexAtSwitchTimes (FrameIndexAtSwitchTimes < AttractionTests.LastFrame) ;         
    FrameIndexAtSwitchTimes    = round(FrameIndexAtSwitchTimes);

    %% Find location of worms coordinate at the time when STIMULUS ENTERS THE ARENA (i.e. at ARENA DELAY)  
    xmat_smoothed                                                   = Data.Coordinates_X_Smoothed;     
    MIN_X                                                           = min(xmat_smoothed(:));
    LocationForNaNs                                                 = MIN_X;
    WormsLocationAtSwitchTimes                                      = xmat_smoothed(:,FrameIndexAtSwitchTimes);
    WormsLocationAtSwitchTimes(isnan(WormsLocationAtSwitchTimes))   = LocationForNaNs;
    WormsLocationBeginning                                          = xmat_smoothed(:,ArenaDelay);
    WormsLocationBeginning(isnan(WormsLocationBeginning))           = LocationForNaNs;
    WormsLocationLastFrame                                          = xmat_smoothed(:,end);
    WormsLocationLastFrame(isnan(WormsLocationLastFrame))           = LocationForNaNs;
    % Take into account only the distance from the arena entrance !!
    WormsLocationAtSwitchTimes                                      = WormsLocationAtSwitchTimes - MIN_X;
    WormsLocationBeginning                                          = WormsLocationBeginning     - MIN_X;
    WormsLocationLastFrame                                          = WormsLocationLastFrame     - MIN_X;
    
    %% Find frame delay due to the worms coordinates at these stimulus switch times   
    FramesDelayAtSwitchTimes                                        = round(WormsLocationAtSwitchTimes *  FlowDelay.DelayWithinArena_framesPER_pixel);  % Delay relative to the stimulus at the arena entrance time.
    FramesDelayAtBeginning                                          = round(WormsLocationBeginning *  FlowDelay.DelayWithinArena_framesPER_pixel);  % Delay relative to the stimulus at the arena entrance time.
    FramesDelayAtLastFrame                                          = round(WormsLocationLastFrame *  FlowDelay.DelayWithinArena_framesPER_pixel);  % Delay relative to the stimulus at the arena entrance time.

    %% Find the new frames indices
    FrameIndicesMat = FindIndicesForDelayedMatrices (ArenaDelay, FrameIndexAtSwitchTimes, FramesDelayAtSwitchTimes, FramesDelayAtBeginning, FramesDelayAtLastFrame, NumOfFrames);
    NaNFrames       = isnan(FrameIndicesMat);
    % Generate delayed matrix based on the delay-time information and the original behavior matrix    
    for f_num = 1:length(FieldsToCorrect)       % Correct all matrices by shifting it 'ArenaDelay' frames forward  
        field_name = FieldsToCorrect{f_num};
        % In 'BehaviorProbability' matrices rows are Behavior codes number rather than worm index. The probability matrices will be re-calculated (see below) from the delayed behavior code matrices.  
        if isempty(strfind(field_name,'BehaviorProbability'))
            Data_AllDelaysCorrected.(field_name) = CreateDelayedMatrix ( Data.(field_name), FrameIndicesMat, NaNFrames);    % syntax:   DelayedMatrix = CreateDelayedMatrix ( Matrix, FrameIndicesMat)  
        end
    end
    Data_AllDelaysCorrected.FrameIndicesMat      = FrameIndicesMat;
    Data_AllDelaysCorrected.NaN(NaNFrames)       = true;
    
    % Quality assurance: Worm Delay cannot be more than inter-arena delay, up to numerical fluctuations. 
    MaximumDelayInSeconds = max( [FramesDelayAtSwitchTimes(:)', FramesDelayAtBeginning(:)', FramesDelayAtLastFrame(:)'])/Frame_per_sec;
    
    DelayQA_CannotBeSignificantlyLargerThan1 = (MaximumDelayInSeconds - InterArenaDelay)/MaximumDelayInSeconds; % 
    if DelayQA_CannotBeSignificantlyLargerThan1 > 1.1   % 10% more than the calculated inter-arena delay
        disp(['QA failure. worm-position-delay calculation in arena ',num2str(ArenaID),' is inconsistent.']);
    end
    
    %% re-calculate probability matrices from the delayed behavior code matrices.
    % The function below is taken from an internal function of: 'SegmentBehavior_v04(File, ArenaID)'  
    % Note that this is necessary only for the AllDelayCorrected normalization. 
    BehaviorCodeNumbers                                                         = BehaviorCode.LowLevel_CodeNumbers;  
    BehaviorCodeMAT                                                             = Data_AllDelaysCorrected.BehaviorCode_LowLevel;
    [BehaviorProbability, BehaviorProbability_Smoothed, DetectionAndSegmentationProbabilityVector] = CalculateBehaviorProbability (BehaviorCodeMAT, BehaviorCodeNumbers, Interpolation_Factor);
    Data_AllDelaysCorrected.BehaviorProbability_LowLevel                        = BehaviorProbability;
    Data_AllDelaysCorrected.BehaviorProbability_Smoothed_LowLevel               = BehaviorProbability_Smoothed;

    BehaviorCodeNumbers                                                          = BehaviorCode.HighLevel_CodeNumbers;  
    BehaviorCodeMAT                                                              = Data_AllDelaysCorrected.BehaviorCode_HighLevel;
    [BehaviorProbability, BehaviorProbability_Smoothed]                          = CalculateBehaviorProbability (BehaviorCodeMAT, BehaviorCodeNumbers, Interpolation_Factor);
    Data_AllDelaysCorrected.BehaviorProbability_HighLevel                        = BehaviorProbability;
    Data_AllDelaysCorrected.BehaviorProbability_Smoothed_HighLevel               = BehaviorProbability_Smoothed;

end


return

function FrameIndicesMat = FindIndicesForDelayedMatrices (DeviceDelay_InFrames, FrameIndexAtSwitchTimes, FramesDelayAtSwitchTimes, FramesDelayAtBeginning, FramesDelayAtLastFrame, LastFrame)

FramesDelayAtSwitchTimes = round(FramesDelayAtSwitchTimes);    % Use round indices for the important switch frames
NumOfWorms               = size(FramesDelayAtSwitchTimes,1);
NumOfSwitches            = length(FrameIndexAtSwitchTimes);
FrameIndicesMat          = zeros(NumOfWorms,LastFrame,'single')*NaN;

FrameIndexAtSwitchTimes_StimulusLines = FrameIndexAtSwitchTimes - DeviceDelay_InFrames; 
NumOfFramesBetweenSwitches_REAL     = diff([1, FrameIndexAtSwitchTimes_StimulusLines+1 , LastFrame]);  % The exact switch time is considered the PREVIOUS stimulus

for w_ind = 1:NumOfWorms
    Current_IndexInMatrix = 1;                    
    Current_Frame         = DeviceDelay_InFrames+FramesDelayAtBeginning(w_ind);         
    PositionCorrected_FrameIndexAtSwitchTimes = FrameIndexAtSwitchTimes + FramesDelayAtSwitchTimes(w_ind,:); 
    PositionCorrected_FrameIndexAtSwitchTimes = PositionCorrected_FrameIndexAtSwitchTimes(PositionCorrected_FrameIndexAtSwitchTimes<LastFrame); % in order to avoid possible non-existing frames... 
    PositionCorrected_FrameIndexAtSwitchTimes = [PositionCorrected_FrameIndexAtSwitchTimes LastFrame];
    
    for s_ind = 1: NumOfSwitches  
        SwitchFrameIndex = PositionCorrected_FrameIndexAtSwitchTimes(s_ind);
        NumOfFrames_Real = NumOfFramesBetweenSwitches_REAL(s_ind);  
%         FrameInterval    = IntervalBetweenFrames(s_ind);               
        FrameInterval    = (SwitchFrameIndex - Current_Frame)/(NumOfFrames_Real-1);               
        
        FrameIndicesMat(w_ind, Current_IndexInMatrix: (NumOfFrames_Real+Current_IndexInMatrix-1)) =   Current_Frame :  FrameInterval : SwitchFrameIndex ; 
        
        Current_IndexInMatrix = Current_IndexInMatrix + NumOfFrames_Real;
        Current_Frame         = SwitchFrameIndex + 1;               
    end
    % Last frame
    % neglect frame delay at last switch
    LastIndicesInMatrix = Current_IndexInMatrix:LastFrame;
    FrameInterval       = (LastFrame+DeviceDelay_InFrames+FramesDelayAtLastFrame(w_ind) - Current_Frame)/(length(LastIndicesInMatrix)-1); 
    LastFrameValues     = Current_Frame :  FrameInterval : (LastFrame+DeviceDelay_InFrames+FramesDelayAtLastFrame(w_ind));
    
    FrameIndicesMat(w_ind, LastIndicesInMatrix) =  LastFrameValues ; 
end

FrameIndicesMat                            = round(FrameIndicesMat);
FrameIndicesMat(FrameIndicesMat>LastFrame) = NaN;

return

function DelayedMatrix = CreateDelayedMatrix ( Matrix, FrameIndicesMat, NaNMat) 

X                  = whos('Matrix');
Matrix_Class       = X.class;

FrameIndicesNoNaN         = FrameIndicesMat;
FrameIndicesNoNaN(NaNMat) = 1; 

DelayedMatrix = zeros(size(Matrix),'single')*NaN;
for w_ind = 1:size(FrameIndicesMat,1)
    frame_indices_vec       = FrameIndicesNoNaN(w_ind,:);    
    DelayedMatrix(w_ind, :) = Matrix(w_ind, frame_indices_vec);       
end
DelayedMatrix(NaNMat) = NaN;

 % NOTE the conversion below assign zeros instead of NaNs for uint8, uint16 and logical matrices 
if strcmpi(Matrix_Class,'uint8')
    DelayedMatrix = uint8(DelayedMatrix);
elseif strcmpi(Matrix_Class,'uint16')
    DelayedMatrix = uint16(DelayedMatrix);
elseif strcmpi(Matrix_Class,'logical')
    DelayedMatrix(NaNMat) = 0;
    DelayedMatrix         = logical(DelayedMatrix);    
end

return

function [BehaviorProb, BehaviorProb_Smoothed, DetectionAndSegmentationProbabilityVector] = CalculateBehaviorProbability (BehaviorCode, BehaviorCodeNumbers, Interpolation_Factor) 
%% Extract probabilities from BehaviorCode matrix
% This will work to any Behavior Code ASSUMING Code '0' is OutOfBound !!!!!!!!!!!!!!!!!!!!
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
              
% Real_Values --> [true/false] segmentation for each worm ID (row) and frame (col). 

%%
BehaviorCode             = single(BehaviorCode);
BehaviorCodeNumbers      = single(BehaviorCodeNumbers);
NumOfWorms               = size(BehaviorCode,1);

behhist                    = hist(BehaviorCode,BehaviorCodeNumbers);    % How many tracks per behavior code (row) and frame (col)  --> THIS STARTS FROM CODE '0' that represent out of bound segmentation !!!!!!
behhist_SegmentedTracks    = behhist(2:end,:);                          % Taking only non-out of bound segmentation. In 'behhist_SegmentedTracks' row 'i' corresponds to BehaviorCode 'i' 

NumberOfSegmentedTracks    = sum(behhist_SegmentedTracks,1);  

DetectionAndSegmentationProbabilityVector = single(NumberOfSegmentedTracks/NumOfWorms);
% Probability of each (not out of bound) behavior relative to the total number of tracks that were segmented for behaviors (detected + not out of bound) 
BehaviorProb          = single(behhist_SegmentedTracks ./ repmat(NumberOfSegmentedTracks,size(behhist_SegmentedTracks,1),1)); 

% Correct for possible frames were no worms were detected
BehaviorProb_Smoothed = zeros(size(BehaviorProb),'single')*NaN;
for row = 1:size(BehaviorProb,1)
    Vec          = BehaviorProb(row,:);
    NaNs         = isnan(Vec);
    VecNoNaNs    = double(Vec(~NaNs));              
    x            = 1:length(VecNoNaNs);
    SmoothVec    = single(csaps(x,VecNoNaNs,Interpolation_Factor, x));
    
    BehaviorProb_Smoothed(row,~NaNs) = SmoothVec;  
end

return






