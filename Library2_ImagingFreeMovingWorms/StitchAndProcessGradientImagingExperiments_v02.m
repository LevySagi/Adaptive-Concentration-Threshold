function Data = StitchAndProcessGradientImagingExperiments_v02

% September 2019, Sagi Levy

%% Inputs
MovieNames     = {'20161215_1','20161215_2','20160404_1_CX14215','20160404_2_CX14215','20160412_2',...
                  '20170221_1','20170223_1','20170223_2','20170224_1','20170224_2'};
DirectoryNames = {'20161215_1','20161215_2','20160404_1_CX14215','20160404_2_CX14215','20160412_2_CX14215',...
                  '20170221_1','20170223_1','20170223_2','20170224_1','20170224_2'};
PathName       = 'F:\';

%% Concatenate data from all gradient imaging experiments
AllFiles          = cell(1,length(MovieNames));
AllData           = cell(1,length(MovieNames));
TopConcentrations = zeros(1,length(MovieNames),'single');
SaveFileName      = [PathName,'GradientImagingMetaAnalysis.mat'];

for m_ind = 1:length(MovieNames)
    disp(['loading file ',num2str(m_ind),' / ',num2str(length(MovieNames))]);
    FileName = [PathName,DirectoryNames{m_ind},'\',MovieNames{m_ind},'_DataMatrices_WithGradient.mat'];
    load(FileName,'File','Data');
    AllData{m_ind}           = Data;
    AllFiles{m_ind}          = File;
    TopConcentrations(m_ind) = File.ExpInfo.ConcentrationAtArenaTop;
end
Data                        = ConcatenateData(AllData, AllFiles);
Data.ExpInfo.DirectoryNames = DirectoryNames;
Data.ExpInfo.MovieNames     = MovieNames;
Data.ExpInfo.PathName       = PathName;

%% Additional data processing 
Data                        = AnalyzeGradientNavigationDecisions(Data);         % e.g. extract downgradient termination decisions    
Data                        = Extract_SingleFrameResolution_BehaviorCode(Data); % For exact Behavior code  
Data.GradientFramesMatrix   = CalculateGradientFramesMatrix(Data);              % Calculate gradient features
Data.BehaviorCode_LowLevel(Data.BehaviorCode_LowLevel==2) = 1;                  % Re-segment all forward curves as forward 

save(SaveFileName,'Data','-v7.3')  
return

%% inline functions
function Data = ConcatenateData(AllData, AllFiles)

NumOfWormsPerMovie  = zeros(1,length(AllData),'single');
NumOfFramesPerMovie = zeros(1,length(AllData),'single');
for m_ind = 1:length(AllData)
    NumOfWormsPerMovie(m_ind)  = size(AllData{m_ind}.NaN,1);
    NumOfFramesPerMovie(m_ind) = size(AllData{m_ind}.NaN,2);
end
NumOfWorms  = sum(NumOfWormsPerMovie);
NumOfFrames = max(NumOfFramesPerMovie);

m_ind = 1;
CurrentData = AllData{m_ind}; 
DataFields  = fieldnames(CurrentData);
CurrentNumberOfTracks = NumOfWormsPerMovie(m_ind);
CurrentNumberOfFrames = NumOfFramesPerMovie(m_ind);
%initialization
for f_ind = 1:length(DataFields)
    CurrentField          = DataFields{f_ind};
    MAT = CurrentData.(CurrentField);
    if (size(MAT,1)==CurrentNumberOfTracks) && (size(MAT,2)==CurrentNumberOfFrames)
        X=whos('MAT'); 
        if     strcmpi(X.class,'logical'),  Data.(CurrentField) = false(NumOfWorms, NumOfFrames);
        elseif strcmpi(X.class,'single'),   Data.(CurrentField) = zeros(NumOfWorms, NumOfFrames,'single')*NaN;
        elseif strcmpi(X.class,'uint8'),    Data.(CurrentField) = zeros(NumOfWorms, NumOfFrames,'uint8');
        elseif strcmpi(X.class,'uint16'),   Data.(CurrentField) = zeros(NumOfWorms, NumOfFrames,'uint16');
        end  
    end    
end
Data.DowngradientStartFrames = false(NumOfWorms, NumOfFrames);
Data.DowngradientEndFrames   = false(NumOfWorms, NumOfFrames);
Data.ActivationStartFrames   = false(NumOfWorms, NumOfFrames);
Data.ActivationPeakFrames    = false(NumOfWorms, NumOfFrames);
Data.ActivationFrameOOB      = false(NumOfWorms, NumOfFrames);
Data.ActivationSharp         = false(NumOfWorms, NumOfFrames);

ExpInfo.ConcentrationAtArenaTop    = zeros(1,NumOfWorms,'single'); 
ExpInfo.ConcentrationAtArenaBottom = zeros(1,NumOfWorms,'single'); 
ExpInfo.MovieNumber                = zeros(1,NumOfWorms,'single'); 
ExpInfo.TrackNumber                = zeros(1,NumOfWorms,'single'); 
ExpInfo.GradientStartFrame         = zeros(1,NumOfWorms,'single'); 
ExpInfo.GradientEndFrame           = zeros(1,NumOfWorms,'single'); 
ExpInfo.GradientRiseTime_InitiationTo90Percent_InSec = zeros(1,NumOfWorms,'single'); 

worm_ind = 0;
for m_ind = 1:length(AllData)
    CurrentData =  AllData{m_ind}; 
    CurrentFile =  AllFiles{m_ind}; 
    CurrentNumberOfTracks = NumOfWormsPerMovie(m_ind);
    CurrentNumberOfFrames = NumOfFramesPerMovie(m_ind);
    
    for tr=1:CurrentNumberOfTracks       
        worm_ind = worm_ind+1;
        ExpInfo.ConcentrationAtArenaTop(worm_ind)   = CurrentFile.ExpInfo.ConcentrationAtArenaTop;
        ExpInfo.ConcentrationAtArenaBottom(worm_ind)= CurrentFile.ExpInfo.ConcentrationAtArenaBottom;    
        ExpInfo.MovieNumber(worm_ind)               = m_ind;
        ExpInfo.TrackNumber(worm_ind)               = tr;
        ExpInfo.GradientStartFrame(worm_ind)        = CurrentFile.GradientStartFrame;
        ExpInfo.GradientEndFrame(worm_ind)          = CurrentFile.GradientEndFrame;
        ExpInfo.GradientRiseTime_InitiationTo90Percent_InSec(worm_ind) = CurrentFile.GradientRiseTime_InitiationTo90Percent_InSec;    
        
        for f_ind = 1:length(DataFields)
            CurrentField          = DataFields{f_ind};
            MAT = CurrentData.(CurrentField);
            if (size(MAT,1)==CurrentNumberOfTracks) && (size(MAT,2)==CurrentNumberOfFrames)                
                Data.(CurrentField)(worm_ind,1:CurrentNumberOfFrames)= MAT(tr,1:CurrentNumberOfFrames);
            end    
        end
        Data.DowngradientStartFrames(worm_ind,CurrentData.DowngradientSegments.StartFrames{tr}) = true;
        Data.DowngradientEndFrames(worm_ind,CurrentData.DowngradientSegments.EndFrames{tr})     = true;
        Data.ActivationStartFrames(worm_ind,CurrentData.ActivationSegments.InitiationFrames{tr})= true;
        Data.ActivationPeakFrames(worm_ind,CurrentData.ActivationSegments.PeakFrames{tr})       = true;
        Data.ActivationFrameOOB(worm_ind,CurrentData.ActivationSegments.InitiationFrames{tr})   = CurrentData.ActivationSegments.OOB{tr};
        Data.ActivationFrameOOB(worm_ind,CurrentData.ActivationSegments.PeakFrames{tr})         = CurrentData.ActivationSegments.OOB{tr};
        Data.ActivationSharp(worm_ind,CurrentData.ActivationSegments.InitiationFrames{tr})      = CurrentData.ActivationSegments.SharpActivation{tr};
        Data.ActivationSharp(worm_ind,CurrentData.ActivationSegments.PeakFrames{tr})            = CurrentData.ActivationSegments.SharpActivation{tr};                
    end
end

Data.ExpInfo = ExpInfo;

return

function GradientFramesMatrix = CalculateGradientFramesMatrix(Data)

%% Initialization
MatrixSize                   = size(Data.NaN);
NumberOfFrames               = MatrixSize(2);
NumberOfWorms                = MatrixSize(1);
OOB                          = Data.BehaviorCode_LowLevel==0;  

%% Remove frames before and after gradient is established.
%% Remove control non-gradient experiments 
%% Remove OOB
GradientFramesMatrix = false(MatrixSize);
for tr=1:NumberOfWorms
    FirstFrameOfOdor     = Data.ExpInfo.GradientStartFrame(tr);   % Gradient STARTS to be established at time 0  
    FirstFrameOfGradient = FirstFrameOfOdor + Data.ExpInfo.GradientRiseTime_InitiationTo90Percent_InSec(tr)*30; % Gradient is established   
    LastFrameOfGradient  = Data.ExpInfo.GradientEndFrame(tr);     % Gradient STARTS to vanish at this frame 
    LastFrameOfGradient  = min([LastFrameOfGradient NumberOfFrames]);
            
    if Data.ExpInfo.ConcentrationAtArenaTop(tr) == Data.ExpInfo.ConcentrationAtArenaBottom(tr) % no gradient
        FirstFrameOfGradient = 1;
        LastFrameOfGradient  = NumberOfFrames;
    end

    GradientFrames       = false(1,NumberOfFrames);
    GradientFrames(FirstFrameOfGradient:LastFrameOfGradient) = true;
    GradientFramesMatrix(tr,:) = GradientFrames;
end
ExperimentsWithOdorGradient = Data.ExpInfo.ConcentrationAtArenaTop ~= Data.ExpInfo.ConcentrationAtArenaBottom;
GradientFramesMatrix(~ExperimentsWithOdorGradient, :) = false;

GradientFramesMatrix(OOB) = false;


return

function Data = AnalyzeGradientNavigationDecisions(Data)

%% Initialization of free Parameters
Parameters.IncludeManualCorrection = true;   % Manual correction of downgradient end frames, using Excel file
Parameters.FrameRate        = 30;
Parameters.MinimalCTXperSec = 0.01; % CTX units per second
Parameters.MaxDwellTime     = 5;    % seconds, concatinate segments that have less than 'MaxDwellTime' of dwelling time between them    
Parameters.MinCTX_for_LongGradientMovement     = 0.2;
Parameters.MinCTX_for_NotShortGradientMovement = 0.05;
Parameters.SearchIntervalForFittingNonSmoothed_InSeconds = 3; 
Parameters.LPS   = designfilt('lowpassfir', ...
                              'PassbandFrequency',0.05,'StopbandFrequency',0.5, ...
                              'PassbandRipple',1,'StopbandAttenuation',80, ...
                              'DesignMethod','equiripple','SampleRate',Parameters.FrameRate);
Parameters.D_LPS = round(mean(grpdelay(Parameters.LPS)));
%  fvtool(Parameters.LPS)  
Parameters.plotme = false;

%% Read manual correction file
if Parameters.IncludeManualCorrection
    load('D:\ImagingInGradientDevices\ManualInspection_GradientLocomotion.mat','num','txt','raw')
    ManualCorrection.EndFramesRows   = find(strcmpi(txt(2:end,2),'End Frames'));
    ManualCorrection.StartFramesRows = find(strcmpi(txt(2:end,2),'Start Frames'));
    ManualCorrection.TrackNumbers    = num(:,1);
    ManualCorrection.FrameNumbers    = num(:,3:end);
    Parameters.ManualCorrection      = ManualCorrection;
end

%% Loop over tracks
NumOfTracks = size(Data.NaN,1);
NumOfFrames = size(Data.NaN,2);
SingleTracksStructures                                   = cell(1,NumOfTracks);
Data.DowngradientNavigation                              = false(NumOfTracks,NumOfFrames);
Data.UpgradientNavigation                                = false(NumOfTracks,NumOfFrames);
Data.Dwell                                               = false(NumOfTracks,NumOfFrames);
Data.DowngradientNavigationEndFrames_All                 = false(NumOfTracks,NumOfFrames);
Data.DowngradientNavigationEndFrames_WithoutOOB          = false(NumOfTracks,NumOfFrames);
Data.DowngradientNavigationEndFrames_FastDirectionChange = false(NumOfTracks,NumOfFrames);
Data.UpgradientNavigationEndFrames_All                   = false(NumOfTracks,NumOfFrames);
Data.UpgradientNavigationEndFrames_WithoutOOB            = false(NumOfTracks,NumOfFrames);
Data.UpgradientNavigationEndFrames_FastDirectionChange   = false(NumOfTracks,NumOfFrames);

for tr = 1:NumOfTracks
    CurrentStructureOut = SingleTrackGradientAnalysis(Data, tr, Parameters);  
    Data.DowngradientNavigation(tr,:)                               = CurrentStructureOut.Downgradient;
    Data.UpgradientNavigation(tr,:)                                 = CurrentStructureOut.Upgradient;
    Data.Dwell(tr,:)                                                = CurrentStructureOut.Dwell;
    Data.DowngradientNavigationEndFrames_All(tr,:)                  = CurrentStructureOut.DowngradientEndFrames_All                 ;
    Data.DowngradientNavigationEndFrames_WithoutOOB(tr,:)           = CurrentStructureOut.DowngradientEndFrames_WithoutOOB          ;
    Data.DowngradientNavigationEndFrames_FastDirectionChange(tr,:)  = CurrentStructureOut.DowngradientEndFrames_FastDirectionChange ;
    Data.UpgradientNavigationEndFrames_All(tr,:)                    = CurrentStructureOut.UpgradientEndFrames_All                   ;
    Data.UpgradientNavigationEndFrames_WithoutOOB(tr,:)             = CurrentStructureOut.UpgradientEndFrames_WithoutOOB            ;
    Data.UpgradientNavigationEndFrames_FastDirectionChange(tr,:)    = CurrentStructureOut.UpgradientEndFrames_FastDirectionChange   ;    
    SingleTracksStructures{tr}                                      = CurrentStructureOut;
end
Data.GradientNavigationStructures = SingleTracksStructures;

return

function Data = Extract_SingleFrameResolution_BehaviorCode(Data)

BehaviorCode_LowLevel       = Data.BehaviorCode_LowLevel;
BehaviorCode_LowLevel_NEW   = Data.BehaviorCode_LowLevel;
NumOfWorms                  = size(Data.BehaviorCode_LowLevel,1);
NumOfFrames                 = size(Data.BehaviorCode_LowLevel,2);
%%  Free Parameters
MaximumPauseTimeAllowedForSwitching   = 2;     % seconds
MinimumTimeForFwdOrRevMovement        = 0.333; % 1/3 seconds
MinimumFramesForOOB                   = 3;     % OOB with 2 frames or less will get the last behavior

MaximumPauseFramesAllowedForSwitching = 30 * MaximumPauseTimeAllowedForSwitching; % frames
MinimumFramesForFwdOrRevMovement      = 30 * MinimumTimeForFwdOrRevMovement;      % frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Get rid of small Out-of-bound segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
OutOfBounds = (Data.BehaviorCode_LowLevel==0);   %  Initialization

for tr = 1:NumOfWorms
    CurrentOOB      = OutOfBounds(tr,:);
    CurrentLongOOB  = FindSegment(CurrentOOB, MinimumFramesForOOB);   
    CurrentShortOOB = CurrentOOB & ~CurrentLongOOB;
    if isempty(find(CurrentShortOOB,1,'first'))
        continue
    end
    [~, ShortOOBStartFrames, ShortOOBEndFrames]= FindSegment(CurrentShortOOB, 0);    
    for start_ind = 1:length(ShortOOBStartFrames)
        CurrentStartFrame = ShortOOBStartFrames(start_ind);
        if CurrentStartFrame==1
            continue
        end
        CurrentEndFrame   = ShortOOBEndFrames(start_ind);
        Indices           = CurrentStartFrame:CurrentEndFrame;
        BehaviorCode_LowLevel_NEW(tr,Indices) = BehaviorCode_LowLevel_NEW(tr,CurrentStartFrame-1);        
    end    
end
% load('BehaviorColormapLowLevelForImaging2.mat','BehaviorColormapLowLevel');
% figure('name','old'); imagesc(BehaviorCode_LowLevel); colormap(BehaviorColormapLowLevel); colorbar;
% figure('name','new'); imagesc(BehaviorCode_LowLevel_NEW); colormap(BehaviorColormapLowLevel); colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Find single-frame resolution direction switching in several steps:
%  Step 1:  Find *short enough* pauses between forward and reverse and vice versa   
%  Step 2:  If there is frame of "sharp turn" within the segment - mark it.     
%           If there is NO frame of "sharp turn" within the segment then estimate and then mark it:     
%               2A. by velocity amplitude and angular velocity   
%               2B. If results are not clear take the midpoint
%           Marking of these Sharp turns: RevFwdSwitch or FwdRevSwitch 
%  Step 3:  Elongate reversal and forward movements until the switch points     
%  Step 4:  Find: reversal initiation, reversal termination    
%                 (regardless of future direction, i.e. including pauses, omegas BOT NOT INCLUDING OOB)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

ForwardOrCurve              = (BehaviorCode_LowLevel_NEW==1)|(BehaviorCode_LowLevel_NEW==2);
Pause                       = (BehaviorCode_LowLevel_NEW==3);
Reverse                     = (BehaviorCode_LowLevel_NEW==4);
SharpTurns                  = (BehaviorCode_LowLevel_NEW==6);
PauseOrSharpTurn            = Pause | SharpTurns;
AngularVelocity             = Data.AngularVelocity;

SwitchingRevToFwd    = false(size(Pause));
SwitchingFwdToRev    = false(size(Pause));
for tr = 1:NumOfWorms
    CurrentReverse          = Reverse(tr,:);
    CurrentFwd              = ForwardOrCurve(tr,:);
    CurrentPauseOrSharpTurn = PauseOrSharpTurn(tr,:);
    CurrentSharpTurn        = SharpTurns(tr,:);
    CurrentAngularVelocity  = AngularVelocity(tr,:);
    
    [CurrentLongRev, ~, ReversalEndFrames] = FindSegment(CurrentReverse, MinimumFramesForFwdOrRevMovement);    
    [CurrentLongFwd, ~,  ForwardEndFrames] = FindSegment(CurrentFwd,     MinimumFramesForFwdOrRevMovement); 
    
    %% Reversal --> Forward 
    EndFrames        = ReversalEndFrames;
    OpposingMovement = CurrentLongFwd;
    
    for end_ind = 1:length(EndFrames)
        CurrentFrame = EndFrames(end_ind);

        %%%%  Step 1:  Find *short enough* pauses between reverse and forward   
        NextOppositeDirectionFrame   = CurrentFrame + find(OpposingMovement(CurrentFrame:end),1,'first')-1;
        if ~isempty(NextOppositeDirectionFrame)
            FramesBetweenOppositeDirections = (CurrentFrame+1):(NextOppositeDirectionFrame-1);
%             if isempty(FramesBetweenOppositeDirections) % this can actually occur after correction of single OOB frames 
%                 disp('debug me...')
%                 keyboard;
%             end
            if length(FramesBetweenOppositeDirections) <= MaximumPauseFramesAllowedForSwitching   % closeby direction change
                if all(CurrentPauseOrSharpTurn(FramesBetweenOppositeDirections))  % all frames within interval or pause or sharpturn
                    %%%% Step 2:  Find or estimate "sharp turn" and mark it.    
                    SharpSwitchFrame = CurrentFrame + find(CurrentSharpTurn(FramesBetweenOppositeDirections),1,'first')-1;    
                    if isempty(SharpSwitchFrame) % Not found, estimate it
                        AngularVelocityWithinPauseInterval = CurrentAngularVelocity(FramesBetweenOppositeDirections);                   
%                         figure; plot(AngularVelocityWithinPauseInterval);
                        [~,MaxInd] = max(abs(AngularVelocityWithinPauseInterval));
                        SharpSwitchFrame = CurrentFrame + MaxInd -1; 
                    end                                                          
                    SwitchingRevToFwd(tr,SharpSwitchFrame) = true;                    
                end
            end                               
        end          
    end
    %%%%  Step 3:  Elongate reversal and forward movements until the switch points     
    CurrentSwitchingFramesRevToFwd = find(SwitchingRevToFwd(tr,:));
    for Frame = CurrentSwitchingFramesRevToFwd
        PreviousReverseFrame = find(CurrentLongRev(1:Frame),1,'last');
        BehaviorCode_LowLevel_NEW(tr,PreviousReverseFrame:(Frame-1)) = 4; % reverse code before switch
        NextForwardFrame     = find(CurrentLongFwd(Frame:end),1,'first') + Frame -1;
        BehaviorCode_LowLevel_NEW(tr,Frame:NextForwardFrame)         = 1; % forward code after switch        
    end
        
    %% Forward --> Reversal
    EndFrames        = ForwardEndFrames;
    OpposingMovement = CurrentLongRev;
    
    for end_ind = 1:length(EndFrames)
        CurrentFrame = EndFrames(end_ind);

        %%%%  Step 1:  Find *short enough* pauses between forward and reverse 
        NextOppositeDirectionFrame   = CurrentFrame + find(OpposingMovement(CurrentFrame:end),1,'first')-1;
        if ~isempty(NextOppositeDirectionFrame)
            FramesBetweenOppositeDirections = (CurrentFrame+1):(NextOppositeDirectionFrame-1);
%             if isempty(FramesBetweenOppositeDirections)
%                 disp('debug me...')
%                 keyboard;
%             end
            if length(FramesBetweenOppositeDirections) <= MaximumPauseFramesAllowedForSwitching   % closeby direction change
                if all(CurrentPauseOrSharpTurn(FramesBetweenOppositeDirections))  % all frames within interval or pause or sharpturn
                    %%%% Step 2:  Find or estimate "sharp turn" and mark it.    
                    SharpSwitchFrame = CurrentFrame + find(CurrentSharpTurn(FramesBetweenOppositeDirections),1,'first')-1;    
                    if isempty(SharpSwitchFrame) % Not found, estimate it
                        AngularVelocityWithinPauseInterval = CurrentAngularVelocity(FramesBetweenOppositeDirections);                   
%                         figure; plot(AngularVelocityWithinPauseInterval);
                        [~,MaxInd] = max(abs(AngularVelocityWithinPauseInterval));
                        SharpSwitchFrame = CurrentFrame + MaxInd -1; 
                    end                                                          
                    SwitchingFwdToRev(tr,SharpSwitchFrame) = true;                    
                end
            end                               
        end          
    end
    %%%%  Step 3:  Elongate reversal and forward movements until the switch points     
    CurrentSwitchingFramesFwdToRev = find(SwitchingFwdToRev(tr,:));
    for Frame = CurrentSwitchingFramesFwdToRev
        PreviousForwardFrame = find(CurrentLongFwd(1:Frame),1,'last');
        BehaviorCode_LowLevel_NEW(tr,PreviousForwardFrame:(Frame-1)) = 1; % forward code before switch
        NextReverseFrame     = find(CurrentLongRev(Frame:end),1,'first') + Frame -1;
        BehaviorCode_LowLevel_NEW(tr,Frame:NextReverseFrame)         = 4; % reverse code after switch        
    end    
end

%%  Step 4:  Find: reversal initiation, reversal termination, Switch points (defined as Rev-->Fwd or Fwd-->Rev)
ReversalInitiationNotOOB = false(size(BehaviorCode_LowLevel_NEW));
ReversalEndingNotOOB     = false(size(BehaviorCode_LowLevel_NEW));
for tr = 1:NumOfWorms
    CurrentReverse = BehaviorCode_LowLevel_NEW(tr,:)==4;
    CurrentOOB     = BehaviorCode_LowLevel_NEW(tr,:)==0;
    [~, CurrentReversalStartFrames, CurrentReversalEndFrames] = FindSegment(CurrentReverse, MinimumFramesForFwdOrRevMovement);    
    % Check start frames
    CurrentReversalStartFrames       = CurrentReversalStartFrames(CurrentReversalStartFrames>1);
    StartAfterOOB                    = CurrentOOB(CurrentReversalStartFrames-1);
    CurrentReversalStartFramesNotOOB = CurrentReversalStartFrames(~StartAfterOOB);
    ReversalInitiationNotOOB(tr,CurrentReversalStartFramesNotOOB) = true;
    % Check end frames
    CurrentReversalEndFrames         = CurrentReversalEndFrames(CurrentReversalEndFrames<NumOfFrames);
    EndBeforeOOB                     = CurrentOOB(CurrentReversalEndFrames+1);
    CurrentReversalEndFramesNotOOB   = CurrentReversalEndFrames(~EndBeforeOOB);
    ReversalEndingNotOOB(tr,CurrentReversalEndFramesNotOOB)       = true;
end

%% Assign to structure
Data.BehaviorCode_LowLevel_OLD = BehaviorCode_LowLevel;
Data.BehaviorCode_LowLevel     = BehaviorCode_LowLevel_NEW;
Data.SwitchingFwdToRev         = SwitchingFwdToRev;
Data.SwitchingRevToFwd         = SwitchingRevToFwd;
Data.ReversalInitiationNotOOB  = ReversalInitiationNotOOB;
Data.ReversalEndingNotOOB      = ReversalEndingNotOOB;

% load('BehaviorColormapLowLevelForImaging2.mat','BehaviorColormapLowLevel');
% figure('name','old'); imagesc(BehaviorCode_LowLevel);     colormap(BehaviorColormapLowLevel); colorbar;
% figure('name','new'); imagesc(BehaviorCode_LowLevel_NEW); colormap(BehaviorColormapLowLevel); colorbar;
% figure; plot(SwitchingFwdToRev(1,:),'k*');
return

function StructureOut = SingleTrackGradientAnalysis(Data, tr, Parameters)
% Example:  tr = 14; % track number

%% initialization
plotme                  = Parameters.plotme;
FrameRate               = Parameters.FrameRate;
MinimalCTXperSec        = Parameters.MinimalCTXperSec; % CTX units per second
MaxDwellTime            = Parameters.MaxDwellTime;     % seconds, concatinate segments that have less than 'MaxDwellTime' of dwelling time between them    
MinCTX_for_LongGradientMovement               = Parameters.MinCTX_for_LongGradientMovement;
MinCTX_for_NotShortGradientMovement           = Parameters.MinCTX_for_NotShortGradientMovement;
SearchIntervalForFittingNonSmoothed_InSeconds = Parameters.SearchIntervalForFittingNonSmoothed_InSeconds; 
LPS   = Parameters.LPS;
D_LPS = Parameters.D_LPS;
% manual correction of downgradient end frames
IncludeManualCorrection = Parameters.IncludeManualCorrection;
if IncludeManualCorrection    
    ManualCorrection     = Parameters.ManualCorrection;    
    CurrrentTrackRows    = find(ManualCorrection.TrackNumbers==tr);
    EndFramesRow         = intersect(CurrrentTrackRows,ManualCorrection.EndFramesRows); 
    StartFramesRow       = intersect(CurrrentTrackRows,ManualCorrection.StartFramesRows); % Only for ADDITIONAL segments, NaN otherwise
    EndFramesRowValues   = ManualCorrection.FrameNumbers(EndFramesRow,:);
    
    if ~isempty(StartFramesRow) % New segments were defined
        StartFramesRowValues                            = ManualCorrection.FrameNumbers(StartFramesRow,:);
        StartFramesColumns                              = find(~isnan(StartFramesRowValues));
        ManualCorrection.AdditionalSegments.StartFrames = StartFramesRowValues(StartFramesColumns);
        ManualCorrection.AdditionalSegments.EndFrames   = EndFramesRowValues(StartFramesColumns);
        EndFramesRowValues(StartFramesColumns)          = NaN;             
    end
    ManualCorrection.ExistingSegments.EndFrames   = EndFramesRowValues(~isnan(EndFramesRowValues));  
end
% per frame units
MinimalCTXperFrame = MinimalCTXperSec/FrameRate; % CTX units per frame
MaxDwellFrames     = MaxDwellTime*FrameRate;
SearchIntervalForFittingNonSmoothed_InFrames = SearchIntervalForFittingNonSmoothed_InSeconds * FrameRate;

%% Filter
Y              = Data.StimulusFiltered(tr,:);
Y_Filtered     = filter(LPS,[Y'; zeros(D_LPS,1)]');    
Y_Filtered     = Y_Filtered(D_LPS+1:end);
Y_Filtered     = Y_Filtered*max(abs(Y))/max(abs(Y_Filtered));
NaNs           = Data.BehaviorCode_LowLevel(tr,:)==0;
MaxNumOfFrames = length(Y);
T              = (1:MaxNumOfFrames)/FrameRate;

%% Define 'Downgradient' and 'Upgradient' by high gradient slope.    
%  Let's take an average of 21.5 micrometer/sec
%  arena: 14*300+100=4300 micrometer = 2 CTX units --> 21.5 micrometer = 0.01 CTX units 
%  MinimalCTXperSec = 0.01;

DiffY        = [diff(Y_Filtered) NaN]; 
Downgradient = (DiffY <= -MinimalCTXperFrame) & ~NaNs ;
Upgradient   = (DiffY >=  MinimalCTXperFrame) & ~NaNs ;
[~, DowngradientStartFrames, DowngradientEndFrames]= FindSegment(Downgradient, 1);    
[~, UpgradientStartFrames, UpgradientEndFrames]    = FindSegment(Upgradient, 1);   

%%%%%%% plots %%%%%%%%
if plotme
    Y_OOB  = Y_Filtered; Y_OOB(~NaNs)=NaN;
    Y_down = Y_Filtered; Y_down(~Downgradient)=NaN;
    Y_up   = Y_Filtered; Y_up(~Upgradient)=NaN;
    DiffY_OOB  = DiffY; DiffY_OOB(~NaNs)=NaN;
    DiffY_down = DiffY; DiffY_down(~Downgradient)=NaN;
    DiffY_up   = DiffY; DiffY_up(~Upgradient)=NaN;
    % figure; plot(T,Y,'.','color',[0.7 0.7 0.7]); hold on; plot(T,Y_Filtered,'r'); plot(T(NaNs),Y_Filtered(NaNs),'k.');
    figure; plot(T,Y_Filtered,'k'); hold on; plot(T,Y_OOB,'color',[0.7 0.7 0.7]); plot(T,Y_down,'r-'); plot(T,Y_up,'b-');
    figure; plot(T,DiffY,'k'); hold on; plot(T,DiffY_OOB,'color',[0.5 0.5 0.5]); plot(T,DiffY_down,'r-'); plot(T,DiffY_up,'b-');
    figure; plot(T,Y_Filtered,'k'); hold on; plot(T,Y_OOB,'color',[0.7 0.7 0.7]); 
            plot(T(DowngradientStartFrames),Y_Filtered(DowngradientStartFrames),'ro');
            plot(T(DowngradientEndFrames),  Y_Filtered(DowngradientEndFrames),'rx');
            plot(T(UpgradientStartFrames),  Y_Filtered(UpgradientStartFrames),'bo');
            plot(T(UpgradientEndFrames),    Y_Filtered(UpgradientEndFrames),'bx');
end

%% Concatinate segments that have less than 'MaxDwellTime' of dwelling time between them  
Dwell = ~Downgradient & ~Upgradient & ~NaNs ;
for f_ind = 1:length(DowngradientStartFrames)-1
    EndFrame       = DowngradientEndFrames(f_ind);
    NextStartFrame = DowngradientStartFrames(f_ind+1);
    if (NextStartFrame-EndFrame-1)<MaxDwellFrames
        Indices =  (EndFrame+1):(NextStartFrame-1);
        if all(Dwell(Indices))|| (length(Indices)< FrameRate/2) % short dwelling or less than 0.5 sec
            Downgradient(Indices) = true;
        end        
    end    
end
Dwell = ~Downgradient & ~Upgradient & ~NaNs ;
for f_ind = 1:length(UpgradientStartFrames)-1
    EndFrame       = UpgradientEndFrames(f_ind);
    NextStartFrame = UpgradientStartFrames(f_ind+1);
    if (NextStartFrame-EndFrame-1)<MaxDwellFrames
        Indices =  (EndFrame+1):(NextStartFrame-1);
        if all(Dwell(Indices))|| (length(Indices)< FrameRate/2)
            Upgradient(Indices) = true;
        end        
    end    
end    
[~, DowngradientStartFrames, DowngradientEndFrames]= FindSegment(Downgradient, 1);    
[~, UpgradientStartFrames, UpgradientEndFrames]    = FindSegment(Upgradient, 1);   

%% Extend segments until turning point when the angle is switched (down->up, or up->down)   
for f_ind = 1:length(DowngradientEndFrames)
    EndFrame = DowngradientEndFrames(f_ind);
    Indices  = (EndFrame+1):(EndFrame+MaxDwellFrames);   
    Indices  = Indices(Indices<=MaxNumOfFrames);
    NextUpgradientFrame = Indices(1)+ find(Upgradient(Indices),1,'first') - 1;
    if ~isempty(NextUpgradientFrame) && ~any(NaNs(Indices)) 
        LastDowngradientIndex      = find(DiffY(Indices)>0,1,'first')-1;
        LastDowngradientIndexFrame = Indices(1)+ LastDowngradientIndex - 1;
        Downgradient(EndFrame:LastDowngradientIndexFrame)              = true;
        Upgradient((LastDowngradientIndexFrame+1):NextUpgradientFrame) = true;
    end    
end
for f_ind = 1:length(UpgradientEndFrames)
    EndFrame = UpgradientEndFrames(f_ind);
    Indices  = (EndFrame+1):(EndFrame+MaxDwellFrames);   
    Indices  = Indices(Indices<=MaxNumOfFrames);
    NextDowngradientFrame = Indices(1)+ find(Downgradient(Indices),1,'first') - 1;
    if ~isempty(NextDowngradientFrame) && ~any(NaNs(Indices)) 
        LastUpgradientIndex      = find(DiffY(Indices)<0,1,'first')-1;
        LastUpgradientIndexFrame = Indices(1)+ LastUpgradientIndex - 1;
        Upgradient(EndFrame:LastUpgradientIndexFrame)                    = true;
        Downgradient((LastUpgradientIndexFrame+1):NextDowngradientFrame) = true;
    end    
end
[~, DowngradientStartFrames, DowngradientEndFrames]= FindSegment(Downgradient, 1);    
[~, UpgradientStartFrames, UpgradientEndFrames]    = FindSegment(Upgradient, 1);   

%%%%%%% plots %%%%%%%%
if plotme
    figure; plot(T,Y_Filtered,'k.'); hold on; plot(T,Y_OOB,'color',[0.7 0.7 0.7]); 
            plot(T(DowngradientStartFrames),Y_Filtered(DowngradientStartFrames),'ro');
            plot(T(DowngradientEndFrames),  Y_Filtered(DowngradientEndFrames),'rx');
            plot(T(UpgradientStartFrames),  Y_Filtered(UpgradientStartFrames),'bo');
            plot(T(UpgradientEndFrames),    Y_Filtered(UpgradientEndFrames),'bx');
end
%% Find Long and not-short gradient locomotion segments
DowngradientAmplitudes          = Y_Filtered(DowngradientStartFrames) - Y_Filtered(DowngradientEndFrames);
UpgradientAmplitudes            = Y_Filtered(UpgradientEndFrames)     - Y_Filtered(UpgradientStartFrames);
LongDowngradientStartFrames     = DowngradientStartFrames(DowngradientAmplitudes > MinCTX_for_LongGradientMovement);
LongDowngradientEndFrames       = DowngradientEndFrames  (DowngradientAmplitudes > MinCTX_for_LongGradientMovement);
LongUpgradientStartFrames       = UpgradientStartFrames  (UpgradientAmplitudes   > MinCTX_for_LongGradientMovement);
LongUpgradientEndFrames         = UpgradientEndFrames    (UpgradientAmplitudes   > MinCTX_for_LongGradientMovement);
NotShortDowngradientStartFrames = DowngradientStartFrames(DowngradientAmplitudes > MinCTX_for_NotShortGradientMovement);
NotShortDowngradientEndFrames   = DowngradientEndFrames  (DowngradientAmplitudes > MinCTX_for_NotShortGradientMovement);
NotShortUpgradientStartFrames   = UpgradientStartFrames  (UpgradientAmplitudes   > MinCTX_for_NotShortGradientMovement);
NotShortUpgradientEndFrames     = UpgradientEndFrames    (UpgradientAmplitudes   > MinCTX_for_NotShortGradientMovement);

LongDowngradient = false(size(Downgradient));
for f_ind = 1:length(LongDowngradientEndFrames)
    StartFrame = LongDowngradientStartFrames(f_ind);
    EndFrame   = LongDowngradientEndFrames(f_ind);
    LongDowngradient(StartFrame:EndFrame) = true;    
end
LongUpgradient   = false(size(Upgradient));
for f_ind = 1:length(LongUpgradientEndFrames)
    StartFrame = LongUpgradientStartFrames(f_ind);
    EndFrame   = LongUpgradientEndFrames(f_ind);
    LongUpgradient(StartFrame:EndFrame) = true;    
end

%% Find EXTENDED segments: that have less than 5*'MaxDwellTime' of dwelling time between them.  
DowngradientExtended = Downgradient;
UpgradientExtended   = Upgradient;
Dwell_ForExtended    = ~DowngradientExtended & ~UpgradientExtended & ~NaNs ;
for f_ind = 1:length(DowngradientStartFrames)-1
    EndFrame       = DowngradientEndFrames(f_ind);
    NextStartFrame = DowngradientStartFrames(f_ind+1);
    if (NextStartFrame-EndFrame-1)<5*MaxDwellFrames
        Indices =  (EndFrame+1):(NextStartFrame-1);
        if all(Dwell_ForExtended(Indices))|| (length(Indices)< FrameRate/2) % short dwelling or less than 0.5 sec
            DowngradientExtended(Indices) = true;
        end        
    end    
end
Dwell_ForExtended    = ~DowngradientExtended & ~UpgradientExtended & ~NaNs ;
for f_ind = 1:length(UpgradientStartFrames)-1
    EndFrame       = UpgradientEndFrames(f_ind);
    NextStartFrame = UpgradientStartFrames(f_ind+1);
    if (NextStartFrame-EndFrame-1)<5*MaxDwellFrames
        Indices =  (EndFrame+1):(NextStartFrame-1);
        if all(Dwell_ForExtended(Indices))|| (length(Indices)< FrameRate/2)
            UpgradientExtended(Indices) = true;
        end        
    end    
end   
[~, ExtendedDowngradientStartFrames, ExtendedDowngradientEndFrames]= FindSegment(DowngradientExtended, 1);    
[~, ExtendedUpgradientStartFrames, ExtendedUpgradientEndFrames]    = FindSegment(UpgradientExtended, 1);  

%% Add Long extended segments  
DowngradientAmplitudes              = Y_Filtered(ExtendedDowngradientStartFrames) - Y_Filtered(ExtendedDowngradientEndFrames);
UpgradientAmplitudes                = Y_Filtered(ExtendedUpgradientEndFrames)     - Y_Filtered(ExtendedUpgradientStartFrames);
LongExtendedDowngradientStartFrames = ExtendedDowngradientStartFrames(DowngradientAmplitudes > MinCTX_for_LongGradientMovement);
LongExtendedDowngradientEndFrames   = ExtendedDowngradientEndFrames  (DowngradientAmplitudes > MinCTX_for_LongGradientMovement);
LongExtendedUpgradientStartFrames   = ExtendedUpgradientStartFrames  (UpgradientAmplitudes   > MinCTX_for_LongGradientMovement);
LongExtendedUpgradientEndFrames     = ExtendedUpgradientEndFrames    (UpgradientAmplitudes   > MinCTX_for_LongGradientMovement);

LongExtendedDowngradient = false(size(Downgradient));
for f_ind = 1:length(LongExtendedDowngradientStartFrames)
    StartFrame = LongExtendedDowngradientStartFrames(f_ind);
    EndFrame   = LongExtendedDowngradientEndFrames(f_ind);
    LongExtendedDowngradient(StartFrame:EndFrame) = true;    
end
LongExtendedUpgradient   = false(size(Upgradient));
for f_ind = 1:length(LongExtendedUpgradientStartFrames)
    StartFrame = LongExtendedUpgradientStartFrames(f_ind);
    EndFrame   = LongExtendedUpgradientEndFrames(f_ind);
    LongExtendedUpgradient(StartFrame:EndFrame) = true;    
end

%% Add intermediate segments that are just a continuation of a previous long segment.   
% Use Long Extended segments to know if an intermediate segment is part of a longer strech.   
FramesToAdd                 = ismember(NotShortDowngradientStartFrames,find(LongExtendedDowngradient));
LongDowngradientStartFrames = unique([LongDowngradientStartFrames; NotShortDowngradientStartFrames(FramesToAdd)]);
LongDowngradientEndFrames   = unique([LongDowngradientEndFrames;   NotShortDowngradientEndFrames(FramesToAdd)]);
FramesToAdd                 = ismember(NotShortUpgradientStartFrames,find(LongExtendedUpgradient));
LongUpgradientStartFrames   = unique([LongUpgradientStartFrames; NotShortUpgradientStartFrames(FramesToAdd)]);
LongUpgradientEndFrames     = unique([LongUpgradientEndFrames;   NotShortUpgradientEndFrames(FramesToAdd)]);

% correct the Long-Segments accordingly
LongDowngradient = false(size(Downgradient));
for f_ind = 1:length(LongDowngradientEndFrames)
    StartFrame = LongDowngradientStartFrames(f_ind);
    EndFrame   = LongDowngradientEndFrames(f_ind);
    LongDowngradient(StartFrame:EndFrame) = true;    
end
LongUpgradient   = false(size(Upgradient));
for f_ind = 1:length(LongUpgradientEndFrames)
    StartFrame = LongUpgradientStartFrames(f_ind);
    EndFrame   = LongUpgradientEndFrames(f_ind);
    LongUpgradient(StartFrame:EndFrame) = true;    
end

%%%%%%% plots %%%%%%%%
if plotme
    figure; plot(T,Y_Filtered,'k'); hold on; plot(T,Y_OOB,'color',[0.7 0.7 0.7]); 
            plot(T(DowngradientStartFrames),    Y_Filtered(DowngradientStartFrames),'r.');
            plot(T(DowngradientEndFrames),      Y_Filtered(DowngradientEndFrames),'r.');
            plot(T(UpgradientStartFrames),      Y_Filtered(UpgradientStartFrames),'b.');
            plot(T(UpgradientEndFrames),        Y_Filtered(UpgradientEndFrames),'b.');
            plot(T(LongDowngradientStartFrames),Y_Filtered(LongDowngradientStartFrames),'ro');
            plot(T(LongDowngradientEndFrames),  Y_Filtered(LongDowngradientEndFrames),'rx');
            plot(T(LongUpgradientStartFrames),  Y_Filtered(LongUpgradientStartFrames),'bo');
            plot(T(LongUpgradientEndFrames),    Y_Filtered(LongUpgradientEndFrames),'bx');
end
%% Correct start and end frames to fit non-smoothed stimulus
for f_ind = 1:length(LongDowngradientEndFrames)  % correct downgradient
    StartFrame = LongDowngradientStartFrames(f_ind);
    StartVal   = Y(StartFrame);
    Indices    = StartFrame + (-SearchIntervalForFittingNonSmoothed_InFrames:SearchIntervalForFittingNonSmoothed_InFrames);
    Indices    = Indices(Indices>0); 
    Indices    = Indices(Indices<=MaxNumOfFrames); 
    [MaxVal, loc] = max(Y(Indices));
    if (MaxVal-0.005 > StartVal) 
        NewIndex = Indices(loc);
        if ~NaNs(NewIndex)
            LongDowngradientStartFrames(f_ind) = NewIndex;
        end
    end    
    EndFrame   = LongDowngradientEndFrames(f_ind);
    EndVal     = Y(EndFrame);
    Indices    = EndFrame+(-SearchIntervalForFittingNonSmoothed_InFrames:SearchIntervalForFittingNonSmoothed_InFrames);
    Indices    = Indices(Indices>0); 
    Indices    = Indices(Indices<=MaxNumOfFrames); 
    [MinVal, loc] = min(Y(Indices));
    if (MinVal+0.005 < EndVal) 
        NewIndex = Indices(loc);
        if ~NaNs(NewIndex)
            LongDowngradientEndFrames(f_ind) = NewIndex;
        end
    end
end
for f_ind = 1:length(LongUpgradientEndFrames)   % correct upgradient
    StartFrame = LongUpgradientStartFrames(f_ind);
    StartVal   = Y(StartFrame);
    Indices    = StartFrame+(-SearchIntervalForFittingNonSmoothed_InFrames:SearchIntervalForFittingNonSmoothed_InFrames);
    Indices    = Indices(Indices>0); 
    Indices    = Indices(Indices<=MaxNumOfFrames); 
    [MinVal, loc] = min(Y(Indices));
    if (MinVal+0.005 < StartVal) 
        NewIndex = Indices(loc);
        if ~NaNs(NewIndex)
            LongUpgradientStartFrames(f_ind) = NewIndex;
        end
    end    
    EndFrame   = LongUpgradientEndFrames(f_ind);
    EndVal     = Y(EndFrame);
    Indices    = EndFrame+(-SearchIntervalForFittingNonSmoothed_InFrames:SearchIntervalForFittingNonSmoothed_InFrames);
    Indices    = Indices(Indices>0); 
    Indices    = Indices(Indices<=MaxNumOfFrames); 
    [MaxVal, loc] = max(Y(Indices));
    if (MaxVal-0.005 > EndVal) 
        NewIndex = Indices(loc);
        if ~NaNs(NewIndex)
            LongUpgradientEndFrames(f_ind) = NewIndex;
        end
    end
end
% Concatinate segments that overlap after correction to non-smoothed vector
SegmentsToConnectedToPreviousOne = find((LongDowngradientStartFrames(2:end) - LongDowngradientEndFrames(1:end-1))<0)+1;
if ~isempty(SegmentsToConnectedToPreviousOne)
    for seg_num = SegmentsToConnectedToPreviousOne' 
        Len = length(LongDowngradientStartFrames);
        LongDowngradientStartFrames = LongDowngradientStartFrames(setdiff(1:Len,seg_num));
        LongDowngradientEndFrames   = LongDowngradientEndFrames(setdiff(1:Len,seg_num-1));
    end
end
SegmentsToConnectedToPreviousOne = find((LongUpgradientStartFrames(2:end) - LongUpgradientEndFrames(1:end-1))<0)+1;
if ~isempty(SegmentsToConnectedToPreviousOne)
    for seg_num = SegmentsToConnectedToPreviousOne' 
        Len = length(LongUpgradientStartFrames);
        LongUpgradientStartFrames = LongUpgradientStartFrames(setdiff(1:Len,seg_num));
        LongUpgradientEndFrames   = LongUpgradientEndFrames(setdiff(1:Len,seg_num-1));
    end
end

%%%%%%% plots %%%%%%%%
if plotme
    Y_OOB_nofilter  = Y; Y_OOB_nofilter(~NaNs)=NaN;
    figure; plot(T,Y,'k'); hold on; plot(T,Y_OOB_nofilter,'color',[0.7 0.7 0.7]); 
            plot(T(LongDowngradientStartFrames),Y(LongDowngradientStartFrames),'ro');
            plot(T(LongDowngradientEndFrames),  Y(LongDowngradientEndFrames),'rx');
            plot(T(LongUpgradientStartFrames),  Y(LongUpgradientStartFrames),'bo');
            plot(T(LongUpgradientEndFrames),    Y(LongUpgradientEndFrames),'bx');
end

%% Manual correction 
if IncludeManualCorrection        
    % If new segments were defined
    if isfield(ManualCorrection,'AdditionalSegments') % New segments were defined
        StartFrames = ManualCorrection.AdditionalSegments.StartFrames';
        EndFrames   = ManualCorrection.AdditionalSegments.EndFrames'  ;        
        LongDowngradientStartFrames = unique([LongDowngradientStartFrames; StartFrames]);
        LongDowngradientEndFrames   = unique([LongDowngradientEndFrames;   EndFrames]);
    end
    % Correction of existing segments end frames
    EndFrames = ManualCorrection.ExistingSegments.EndFrames; 
    for Frame = EndFrames
        [MinVal, index] = min(abs(LongDowngradientEndFrames - Frame)); % find the frame that needs correction
%         if MinVal>5*FrameRate
%             disp(['warning: The correction of Track ',num2str(tr),' frame ',num2str(Frame),' is ',num2str(MinVal/FrameRate),' seconds']) ;
%         end
        LongDowngradientEndFrames(index) = Frame;
    end        
end

%% Final assignment to logical vectors
LongDowngradient = false(size(Downgradient));
for f_ind = 1:length(LongDowngradientEndFrames)
    StartFrame = LongDowngradientStartFrames(f_ind);
    EndFrame   = LongDowngradientEndFrames(f_ind);
    LongDowngradient(StartFrame:EndFrame) = true;    
end
LongUpgradient   = false(size(Upgradient));
for f_ind = 1:length(LongUpgradientEndFrames)
    StartFrame = LongUpgradientStartFrames(f_ind);
    EndFrame   = LongUpgradientEndFrames(f_ind);
    LongUpgradient(StartFrame:EndFrame) = true;    
end
Dwell = ~LongDowngradient & ~LongUpgradient & ~NaNs ;

%%%%%%% plots %%%%%%%%
if plotme
    Y_OOB_nofilter  = Y; Y_OOB_nofilter(~NaNs)=NaN;
    figure; plot(T,Y,'k'); hold on; plot(T,Y_OOB_nofilter,'color',[0.7 0.7 0.7]); 
            plot(T(LongDowngradientStartFrames),Y(LongDowngradientStartFrames),'ro');
            plot(T(LongDowngradientEndFrames),  Y(LongDowngradientEndFrames),'rx');
            plot(T(LongUpgradientStartFrames),  Y(LongUpgradientStartFrames),'bo');
            plot(T(LongUpgradientEndFrames),    Y(LongUpgradientEndFrames),'bx');
end

%% Find Characteristics of each segment
% Gradients that finish with a change in direction (with ot without delay)  
% Gradients that finish with a continuation in the same direction (with ot without delay)     
% Gradients that finish with going out of bounds (with ot without delay)    
DowngradientSegments.FinishWithDirectionChange.StartFrames       = [];
DowngradientSegments.FinishWithDirectionChange.EndFrames         = [];
DowngradientSegments.FinishWithDirectionChange.DelayAfterSegment = [];
DowngradientSegments.FinishWithSameDirection    = DowngradientSegments.FinishWithDirectionChange;
DowngradientSegments.FinishWithOutOfBound       = DowngradientSegments.FinishWithDirectionChange;
UpgradientSegments                              = DowngradientSegments;
LongNaNsNotInGradientSegments                                  = FindSegment(NaNs, FrameRate/2);    
LongNaNsNotInGradientSegments(LongDowngradient|LongUpgradient) = false;

for f_ind = 1:length(LongDowngradientEndFrames)  % Downgradient
    StartFrame = LongDowngradientStartFrames(f_ind);
    EndFrame   = LongDowngradientEndFrames(f_ind);
    if EndFrame == MaxNumOfFrames
        Flag = 3;
        NextFrame = MaxNumOfFrames;
    else   
        Indices    = (EndFrame+1):MaxNumOfFrames;   
        NextUpgradient   = Indices(1)+ find(LongUpgradient(Indices),1,'first') - 1;
        NextDowngradient = Indices(1)+ find(LongDowngradient(Indices),1,'first') - 1;
        NextOOBFrame     = Indices(1)+ find(LongNaNsNotInGradientSegments(Indices),1,'first') - 1;
        [NextFrame, Flag]= min([NextUpgradient, NextDowngradient, NextOOBFrame]); 
        if isempty(NextFrame)
            Flag = 3;
            NextFrame = MaxNumOfFrames;
        end
    end
    switch Flag
        case 1
            CurrentField = 'FinishWithDirectionChange';                        
        case 2
            CurrentField = 'FinishWithSameDirection';                           
        case 3
            CurrentField = 'FinishWithOutOfBound';         
    end
    Delay = NextFrame - EndFrame;    
    DowngradientSegments.(CurrentField).StartFrames       = [DowngradientSegments.(CurrentField).StartFrames,       StartFrame];
    DowngradientSegments.(CurrentField).EndFrames         = [DowngradientSegments.(CurrentField).EndFrames,         EndFrame];
    DowngradientSegments.(CurrentField).DelayAfterSegment = [DowngradientSegments.(CurrentField).DelayAfterSegment, Delay];    
end
for f_ind = 1:length(LongUpgradientEndFrames)   % Upgradient
    StartFrame = LongUpgradientStartFrames(f_ind);
    EndFrame   = LongUpgradientEndFrames(f_ind);
    if EndFrame == MaxNumOfFrames
        Flag = 3;
        NextFrame = MaxNumOfFrames;
    else        
        Indices    = (EndFrame+1):MaxNumOfFrames;   
        NextUpgradient   = Indices(1)+ find(LongUpgradient(Indices),1,'first') - 1;
        NextDowngradient = Indices(1)+ find(LongDowngradient(Indices),1,'first') - 1;
        NextOOBFrame     = Indices(1)+ find(LongNaNsNotInGradientSegments(Indices),1,'first') - 1;
        [NextFrame, Flag]= min([NextUpgradient, NextDowngradient, NextOOBFrame]); 
        if isempty(NextFrame)
            Flag = 3;
            NextFrame = MaxNumOfFrames;
        end
    end
    switch Flag
        case 1
            CurrentField = 'FinishWithSameDirection';                        
        case 2
            CurrentField = 'FinishWithDirectionChange';                           
        case 3
            CurrentField = 'FinishWithOutOfBound';         
    end
    Delay = NextFrame - EndFrame;    
    UpgradientSegments.(CurrentField).StartFrames       = [UpgradientSegments.(CurrentField).StartFrames,       StartFrame];
    UpgradientSegments.(CurrentField).EndFrames         = [UpgradientSegments.(CurrentField).EndFrames,         EndFrame];
    UpgradientSegments.(CurrentField).DelayAfterSegment = [UpgradientSegments.(CurrentField).DelayAfterSegment, Delay];    
end

%%%%%%% plots %%%%%%%%
if plotme
    Y_OOB_nofilter  = Y; Y_OOB_nofilter(~LongNaNsNotInGradientSegments)=NaN;
    figure; 
    plot(T,Y_Filtered,'k'); hold on; plot(T,Y_OOB,'color',[0.7 0.7 0.7]); 
    plot(T(DowngradientSegments.FinishWithOutOfBound.StartFrames),  Y_Filtered(DowngradientSegments.FinishWithOutOfBound.StartFrames),'ro');
    plot(T(DowngradientSegments.FinishWithOutOfBound.EndFrames),    Y_Filtered(DowngradientSegments.FinishWithOutOfBound.EndFrames),'rx');
    plot(T(UpgradientSegments.FinishWithOutOfBound.StartFrames),    Y_Filtered(UpgradientSegments.FinishWithOutOfBound.StartFrames),'bo');
    plot(T(UpgradientSegments.FinishWithOutOfBound.EndFrames),      Y_Filtered(UpgradientSegments.FinishWithOutOfBound.EndFrames),'bx');
    plot(T(DowngradientSegments.FinishWithSameDirection.StartFrames),  Y_Filtered(DowngradientSegments.FinishWithSameDirection.StartFrames),'ro');
    plot(T(DowngradientSegments.FinishWithSameDirection.EndFrames),    Y_Filtered(DowngradientSegments.FinishWithSameDirection.EndFrames),'rv');
    plot(T(UpgradientSegments.FinishWithSameDirection.StartFrames),    Y_Filtered(UpgradientSegments.FinishWithSameDirection.StartFrames),'bo');
    plot(T(UpgradientSegments.FinishWithSameDirection.EndFrames),      Y_Filtered(UpgradientSegments.FinishWithSameDirection.EndFrames),'bv');
    plot(T(DowngradientSegments.FinishWithDirectionChange.StartFrames),Y_Filtered(DowngradientSegments.FinishWithDirectionChange.StartFrames),'ro','markerfacecolor','r');
    plot(T(DowngradientSegments.FinishWithDirectionChange.EndFrames),  Y_Filtered(DowngradientSegments.FinishWithDirectionChange.EndFrames),'rv','markerfacecolor','r');
    plot(T(UpgradientSegments.FinishWithDirectionChange.StartFrames),  Y_Filtered(UpgradientSegments.FinishWithDirectionChange.StartFrames),'bo','markerfacecolor','b');
    plot(T(UpgradientSegments.FinishWithDirectionChange.EndFrames),    Y_Filtered(UpgradientSegments.FinishWithDirectionChange.EndFrames),'bv','markerfacecolor','b');

    figure; 
    plot(T,Y,'k'); hold on; plot(T,Y_OOB_nofilter,'color',[0.7 0.7 0.7]); 
    plot(T(DowngradientSegments.FinishWithOutOfBound.StartFrames),  Y(DowngradientSegments.FinishWithOutOfBound.StartFrames),'ro');
    plot(T(DowngradientSegments.FinishWithOutOfBound.EndFrames),    Y(DowngradientSegments.FinishWithOutOfBound.EndFrames),'rx');
    plot(T(UpgradientSegments.FinishWithOutOfBound.StartFrames),    Y(UpgradientSegments.FinishWithOutOfBound.StartFrames),'bo');
    plot(T(UpgradientSegments.FinishWithOutOfBound.EndFrames),      Y(UpgradientSegments.FinishWithOutOfBound.EndFrames),'bx');
    plot(T(DowngradientSegments.FinishWithSameDirection.StartFrames),  Y(DowngradientSegments.FinishWithSameDirection.StartFrames),'ro');
    plot(T(DowngradientSegments.FinishWithSameDirection.EndFrames),    Y(DowngradientSegments.FinishWithSameDirection.EndFrames),'rv');
    plot(T(UpgradientSegments.FinishWithSameDirection.StartFrames),    Y(UpgradientSegments.FinishWithSameDirection.StartFrames),'bo');
    plot(T(UpgradientSegments.FinishWithSameDirection.EndFrames),      Y(UpgradientSegments.FinishWithSameDirection.EndFrames),'bv');
    plot(T(DowngradientSegments.FinishWithDirectionChange.StartFrames),Y(DowngradientSegments.FinishWithDirectionChange.StartFrames),'ro','markerfacecolor','r');
    plot(T(DowngradientSegments.FinishWithDirectionChange.EndFrames),  Y(DowngradientSegments.FinishWithDirectionChange.EndFrames),'rv','markerfacecolor','r');
    plot(T(UpgradientSegments.FinishWithDirectionChange.StartFrames),  Y(UpgradientSegments.FinishWithDirectionChange.StartFrames),'bo','markerfacecolor','b');
    plot(T(UpgradientSegments.FinishWithDirectionChange.EndFrames),    Y(UpgradientSegments.FinishWithDirectionChange.EndFrames),'bv','markerfacecolor','b');

    figure; 
    plot(T,Y,'k'); hold on; plot(T,Y_OOB_nofilter,'color',[0.7 0.7 0.7]); 
    plot(T(DowngradientSegments.FinishWithOutOfBound.StartFrames),  Y(DowngradientSegments.FinishWithOutOfBound.StartFrames),'ro');
    plot(T(UpgradientSegments.FinishWithOutOfBound.StartFrames),    Y(UpgradientSegments.FinishWithOutOfBound.StartFrames),'bo');
    plot(T(DowngradientSegments.FinishWithSameDirection.StartFrames),  Y(DowngradientSegments.FinishWithSameDirection.StartFrames),'ro');
    plot(T(UpgradientSegments.FinishWithSameDirection.StartFrames),    Y(UpgradientSegments.FinishWithSameDirection.StartFrames),'bo');
    plot(T(DowngradientSegments.FinishWithDirectionChange.StartFrames),Y(DowngradientSegments.FinishWithDirectionChange.StartFrames),'ro','markerfacecolor','r');
    plot(T(UpgradientSegments.FinishWithDirectionChange.StartFrames),  Y(UpgradientSegments.FinishWithDirectionChange.StartFrames),'bo','markerfacecolor','b');
    plot(T(DowngradientSegments.FinishWithOutOfBound.EndFrames),    Y(DowngradientSegments.FinishWithOutOfBound.EndFrames),'rx');
    plot(T(UpgradientSegments.FinishWithOutOfBound.EndFrames),      Y(UpgradientSegments.FinishWithOutOfBound.EndFrames),'bx');
    plot(T(DowngradientSegments.FinishWithSameDirection.EndFrames),    Y(DowngradientSegments.FinishWithSameDirection.EndFrames),'rv');
    plot(T(UpgradientSegments.FinishWithSameDirection.EndFrames),      Y(UpgradientSegments.FinishWithSameDirection.EndFrames),'bv');
    plot(T(DowngradientSegments.FinishWithDirectionChange.EndFrames),  Y(DowngradientSegments.FinishWithDirectionChange.EndFrames),'rv','markerfacecolor','r');
    plot(T(UpgradientSegments.FinishWithDirectionChange.EndFrames),    Y(UpgradientSegments.FinishWithDirectionChange.EndFrames),'bv','markerfacecolor','b');

    figure; 
    Y_LongDowngradient = Y; Y_LongDowngradient(~LongDowngradient)= NaN;
    Y_LongUpgradient   = Y; Y_LongUpgradient(~LongUpgradient)    = NaN;
    plot(T,Y,'k'); hold on; plot(T,Y_OOB_nofilter,'color',[0.7 0.7 0.7]); 
    plot(T,Y_LongDowngradient,'r-');
    plot(T,Y_LongUpgradient,'b-');
    plot(T(DowngradientSegments.FinishWithOutOfBound.StartFrames),  Y(DowngradientSegments.FinishWithOutOfBound.StartFrames),'ro');
    plot(T(UpgradientSegments.FinishWithOutOfBound.StartFrames),    Y(UpgradientSegments.FinishWithOutOfBound.StartFrames),'bo');
    plot(T(DowngradientSegments.FinishWithSameDirection.StartFrames),  Y(DowngradientSegments.FinishWithSameDirection.StartFrames),'ro');
    plot(T(UpgradientSegments.FinishWithSameDirection.StartFrames),    Y(UpgradientSegments.FinishWithSameDirection.StartFrames),'bo');
    plot(T(DowngradientSegments.FinishWithDirectionChange.StartFrames),Y(DowngradientSegments.FinishWithDirectionChange.StartFrames),'ro','markerfacecolor','r');
    plot(T(UpgradientSegments.FinishWithDirectionChange.StartFrames),  Y(UpgradientSegments.FinishWithDirectionChange.StartFrames),'bo','markerfacecolor','b');
    plot(T(DowngradientSegments.FinishWithOutOfBound.EndFrames),    Y(DowngradientSegments.FinishWithOutOfBound.EndFrames),'rx');
    plot(T(UpgradientSegments.FinishWithOutOfBound.EndFrames),      Y(UpgradientSegments.FinishWithOutOfBound.EndFrames),'bx');
    plot(T(DowngradientSegments.FinishWithSameDirection.EndFrames),    Y(DowngradientSegments.FinishWithSameDirection.EndFrames),'rv');
    plot(T(UpgradientSegments.FinishWithSameDirection.EndFrames),      Y(UpgradientSegments.FinishWithSameDirection.EndFrames),'bv');
    plot(T(DowngradientSegments.FinishWithDirectionChange.EndFrames),  Y(DowngradientSegments.FinishWithDirectionChange.EndFrames),'rv','markerfacecolor','r');
    plot(T(UpgradientSegments.FinishWithDirectionChange.EndFrames),    Y(UpgradientSegments.FinishWithDirectionChange.EndFrames),'bv','markerfacecolor','b');

    figure; 
    T = 1:MaxNumOfFrames;
    Y_LongDowngradient = Y; Y_LongDowngradient(~LongDowngradient)= NaN;
    Y_LongUpgradient   = Y; Y_LongUpgradient(~LongUpgradient)    = NaN;
    plot(T,Y,'k'); hold on; plot(T,Y_OOB_nofilter,'color',[0.7 0.7 0.7]); 
    plot(T,Y_LongDowngradient,'r-');
    plot(T,Y_LongUpgradient,'b-');
    plot(T(DowngradientSegments.FinishWithOutOfBound.StartFrames),  Y(DowngradientSegments.FinishWithOutOfBound.StartFrames),'ro');
    plot(T(UpgradientSegments.FinishWithOutOfBound.StartFrames),    Y(UpgradientSegments.FinishWithOutOfBound.StartFrames),'bo');
    plot(T(DowngradientSegments.FinishWithSameDirection.StartFrames),  Y(DowngradientSegments.FinishWithSameDirection.StartFrames),'ro');
    plot(T(UpgradientSegments.FinishWithSameDirection.StartFrames),    Y(UpgradientSegments.FinishWithSameDirection.StartFrames),'bo');
    plot(T(DowngradientSegments.FinishWithDirectionChange.StartFrames),Y(DowngradientSegments.FinishWithDirectionChange.StartFrames),'ro','markerfacecolor','r');
    plot(T(UpgradientSegments.FinishWithDirectionChange.StartFrames),  Y(UpgradientSegments.FinishWithDirectionChange.StartFrames),'bo','markerfacecolor','b');
    plot(T(DowngradientSegments.FinishWithOutOfBound.EndFrames),    Y(DowngradientSegments.FinishWithOutOfBound.EndFrames),'rx');
    plot(T(UpgradientSegments.FinishWithOutOfBound.EndFrames),      Y(UpgradientSegments.FinishWithOutOfBound.EndFrames),'bx');
    plot(T(DowngradientSegments.FinishWithSameDirection.EndFrames),    Y(DowngradientSegments.FinishWithSameDirection.EndFrames),'rv');
    plot(T(UpgradientSegments.FinishWithSameDirection.EndFrames),      Y(UpgradientSegments.FinishWithSameDirection.EndFrames),'bv');
    plot(T(DowngradientSegments.FinishWithDirectionChange.EndFrames),  Y(DowngradientSegments.FinishWithDirectionChange.EndFrames),'rv','markerfacecolor','r');
    plot(T(UpgradientSegments.FinishWithDirectionChange.EndFrames),    Y(UpgradientSegments.FinishWithDirectionChange.EndFrames),'bv','markerfacecolor','b');

end
%% Logical Vectors
Logical_AllDowngradientEndFrames                            = false(1,MaxNumOfFrames); 
Logical_AllDowngradientEndFrames(LongDowngradientEndFrames) = true;

OOB_EndFrames = DowngradientSegments.FinishWithOutOfBound.EndFrames;
OOB_EndFrames = OOB_EndFrames(DowngradientSegments.FinishWithOutOfBound.DelayAfterSegment < FrameRate/2); % 0.5 sec. keep only immediate OOB
Logical_DowngradientEndFrames_WithoutOOB                                                   = false(1,MaxNumOfFrames); 
Logical_DowngradientEndFrames_WithoutOOB(setdiff(LongDowngradientEndFrames,OOB_EndFrames)) = true;

EndFrames = DowngradientSegments.FinishWithDirectionChange.EndFrames;
EndFrames = EndFrames(DowngradientSegments.FinishWithDirectionChange.DelayAfterSegment < FrameRate); % 1 sec. keep only immediate change of direction
Logical_DowngradientEndFrames_WithFastDirectionChange            = false(1,MaxNumOfFrames); 
Logical_DowngradientEndFrames_WithFastDirectionChange(EndFrames) = true;

Logical_AllUpgradientEndFrames                          = false(1,MaxNumOfFrames); 
Logical_AllUpgradientEndFrames(LongUpgradientEndFrames) = true;

OOB_EndFrames = UpgradientSegments.FinishWithOutOfBound.EndFrames;
OOB_EndFrames = OOB_EndFrames(UpgradientSegments.FinishWithOutOfBound.DelayAfterSegment < FrameRate/2); % 0.5 sec. keep only immediate OOB
Logical_UpgradientEndFrames_WithoutOOB                                                 = false(1,MaxNumOfFrames); 
Logical_UpgradientEndFrames_WithoutOOB(setdiff(LongUpgradientEndFrames,OOB_EndFrames)) = true;

EndFrames = UpgradientSegments.FinishWithDirectionChange.EndFrames;
EndFrames = EndFrames(UpgradientSegments.FinishWithDirectionChange.DelayAfterSegment < FrameRate); % 1 sec. keep only immediate change of direction
Logical_UpgradientEndFrames_WithFastDirectionChange            = false(1,MaxNumOfFrames); 
Logical_UpgradientEndFrames_WithFastDirectionChange(EndFrames) = true;

%% Assign To Structure
StructureOut.DowngradientSegments = DowngradientSegments;
StructureOut.UpgradientSegments   = UpgradientSegments;
StructureOut.Downgradient         = LongDowngradient;
StructureOut.Upgradient           = LongUpgradient;
StructureOut.Dwell                = Dwell;
StructureOut.OutOfBounds          = ~Dwell & ~LongDowngradient & ~LongUpgradient;

StructureOut.DowngradientEndFrames_All                 = Logical_AllDowngradientEndFrames;
StructureOut.DowngradientEndFrames_WithoutOOB          = Logical_DowngradientEndFrames_WithoutOOB;
StructureOut.DowngradientEndFrames_FastDirectionChange = Logical_DowngradientEndFrames_WithFastDirectionChange;
StructureOut.UpgradientEndFrames_All                   = Logical_AllUpgradientEndFrames;
StructureOut.UpgradientEndFrames_WithoutOOB            = Logical_UpgradientEndFrames_WithoutOOB;
StructureOut.UpgradientEndFrames_FastDirectionChange   = Logical_UpgradientEndFrames_WithFastDirectionChange;

return

function  [VecTrueIfFrameWithinSegment, StartFrames, EndFrames]= FindSegment(LogicalVec, MinimumSegmentLength)
% Example:
% [FramesWithLongNaNSegments, StartNaNSegments, EndNaNSegments]= FindSegment(NaN_Frames, MaxAllowedMissingFrames);

StartFrames = (find(diff(LogicalVec)==1)+1)'; if LogicalVec(1)==1,   StartFrames = [1; StartFrames]; end
EndFrames   = (find(diff(LogicalVec)==-1))';  if LogicalVec(end)==1, EndFrames   = [EndFrames; length(LogicalVec)]; end
if length(StartFrames) ~= length(EndFrames)
    disp('error in segments calculation. keyboard...');
    keyboard;
end
LongSegments_indices = find((EndFrames - StartFrames+1)>MinimumSegmentLength);
if ~isempty(LongSegments_indices)
    VecTrueIfFrameWithinSegment = false(1,length(LogicalVec));
    for ind = LongSegments_indices'
        VecTrueIfFrameWithinSegment(StartFrames(ind):EndFrames(ind)) = true;    
    end
end        
% Correct to keep only long segments
LogicalVec = VecTrueIfFrameWithinSegment;
StartFrames = (find(diff(LogicalVec)==1)+1)'; if LogicalVec(1)==1,   StartFrames = [1; StartFrames]; end
EndFrames   = (find(diff(LogicalVec)==-1))';  if LogicalVec(end)==1, EndFrames   = [EndFrames; length(LogicalVec)]; end
if length(StartFrames) ~= length(EndFrames)
    disp('error in segments calculation. keyboard...');
    keyboard;
end

return



