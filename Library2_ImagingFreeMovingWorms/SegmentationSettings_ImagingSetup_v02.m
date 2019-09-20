function Settings = SegmentationSettings_ImagingSetup_v02(NoPlots)
% original function:   

%% Arenas position definition
% if true  --> Re-define arenas position, even if they were previously defined. 
% If false --> Re-define arenas only if they were not previously defined.    
Settings.RedefineArenaPosition    = false;      
Settings.DefineIDingArenaPosition = false;      

HeadTailRelativeLocation.ALONG_Midline = [0.07 0.11];  % relative the the worm's length.   0=worm edge. 0.5= Middle of its midline.    Compare 'near edges' areas
HeadTailRelativeLocation.FROM_Midline  = 0.3;      % relative to the worms HALF width. 0=midline, 1= worm edge or a bit beyond it. Compare 'near middle' areas

Settings.HeadTailRelativeLocation = HeadTailRelativeLocation;      

%% Free parameters and smoothing information                         
Settings.SmoothingInfo  = struct(  'InterpolationFactor',           0.1, ...              % Interpolation factor for the path coordinates. NOTE! for 30Hz movies, a value of 1e-4 will be used!
                                   'DeltaFrameSmoothed',            1,   ...              % If less than 1, then the locomotion vectors will be sampled between frames and will be LONGER. 
                                   'MaxAllowedTimeForMissingFrames',1,    ...             % seconds. If more than 1 second of frames are missing: don't allow smoothed vectors to predict the worm path  
                                   'HeadTail_InterpolationFactor',  0.9, ...              % less smoothing in head and tail coordinates than the path coordinates to allow fast movements. 
                                   'HeadTail_MaxAllowedTimeForMissingFrames',  0.3, ...    % 3 frames for 10 Hz. Less smoothing in head and tail coordinates than the path coordinates to allow fast movements. 
                                   'HeadTail_MaxAllowedHeadVelocity',  2 ...              % mm/sec. Above that it will be flagged as possible error and treated accordingly.  
                                  );

%% Thresholds for behavior segmentation (these parameters will be used by the function: 'SegmentSingleTrack_SmallArena_SL_vXX')                              
Thresholds.MinimumTimeForBehaviorDefinition  = 1;                           % sec. Minimum time for properly defining the worm behavior. Below this time behavior will be given a '0' value (as Out of bounds)  

Thresholds.Turns.MinVelocityAngleChange      = 80;                          % deg/frame. It must happen within a time interval of 0.3 seconds (up to deg/3-frames)..
Thresholds.Curve.MinVelocityAngleChange      = 2;                           % deg/frame.
Thresholds.Curve.MinTime_ForCurveDefinition  = 0.5;                         % sec. This will be used to define whether the run starts with a curvature.

Thresholds.Pause.MinTimeForShortPause        = 0.3;                         % sec
Thresholds.Pause.MinTimeForLongPause         = 1;                           % sec. This is also used as the maximum value for short pauses.
Thresholds.Pause.MinTimeForVeryLongPause     = 1.5;                         % sec. This is also used as the maximum value for short pauses.

Thresholds.StraightRun.MinVelocityAmplitude_ForRunDefinition        = 0.06; % mm/sec
Thresholds.StraightRun.MinVelocityAmplitude_ForInermediateVelocity  = 0.03; % mm/sec
Thresholds.StraightRun.MinTimeForShortRun                           = 0.3;  % sec. Used for Out of bounds extension  
Thresholds.StraightRun.AllowedExtensionTime                         = 0.6;  % sec. The start/finish of 'extended straight run' segments are determined by the velocity angle change.  
Thresholds.StraightRun.AllowedPauseTimeInExtendedStraightRun        = 0.5;  % sec. Concatinate runs with short pauses  
Thresholds.StraightRun.AllowedPauseTimeInLongExtendedStraightRun    = 1.5;  % sec. Concatinate runs with pauses. THIS WILL BE USED FOR REVERSAL INITATION DEFINITION  

Thresholds.Omega.MaxNormPerimeterForOmega                     = 0.85;
Thresholds.Omega.MinNormPerimeterOverSizeForOmega             = 0.9;
Thresholds.Omega.MaxEccentricityForOmega                      = 0.9;
Thresholds.Omega.MinTimeForOmegaState_EccentricityCriterion   = 0.8;        % sec. May include also a pause that only part of it is an omega.
Thresholds.Omega.MinTimeForOmegaState_PerimeterCriterion      = 3;          % sec. May include also a pause that only part of it is an omega. 
                                                                            %      Here the worm is going SQUEEZE IN BETWEEN POSTS. This tends to be slower than regular omegas. 
Thresholds.Omega.MinOverlapWithPauseForOmegaExtension         = 50;         % [%]. Minimal percent of overlap between a 'pause-like' state and an 'omega-state' for including the pause as part of the omega 
                                                                        
Thresholds.Pirouette.MaxTimeIntervalBetweenReverseAndOmega           = 2.5;   % sec.                                                                         
Thresholds.Pirouette.MaxTimeIntervalBetweenLastReverse_And_Forward   = 2.5; % sec. For the cases of Reverse--> Omega --> reverse --> short pause --> forward
Thresholds.Pirouette.MaxPauseTimeIntervalWithinPirouetteReversal     = 2;   % sec. For the cases of Reverse--> Omega --> reverse --> short pause --> forward

%% How to normalize the the worm size along the tracks?
%  If it's a long track (above the minimum threshold):  normalize with repect to with median along this single track. 
%  If it's a short track (below the minimum threshold): normalize with repect to with median of all tracks. 
Thresholds.MinTrackLengthForSizeNormalization        = 20; % seconds
Thresholds.MinFractionOfWormPixelsInBound            = 0.9;   % At least 90 perecent of the worm body needs to be In-Bound in order to consider the frame as relevant.
Thresholds.MinFractionOfWormPixelsInBound_For_ID     = 0.999;  % For IDing All the worm body needs to be In-Bound in order to consider the frame as relevant.

Settings.Thresholds = Thresholds; 

%% Plots Information: which plots and movies to display. No parameters here.
if exist('NoPlots','var')
    NoPlotsInfo       = NoPlots;  
else
    NoPlotsInfo       = true;
end    
Settings.PlotInfo = GetPlotInfo(NoPlotsInfo);

%% Old parameters for Dirk code. Not including the parameters that were imbedded within the old functions. 
Settings.DirkOldSetting = GetOldSetting;
                
return

function PlotInfo = GetPlotInfo(NoPlotsInfo)

% PlotInfo.NoPlots                                         = true;  % MASTER FIELD, if true --> all PlotInfo will be set to false and NO PLOTS will be drawn. 
% PlotInfo.NoPlots                                         = false;  % MASTER FIELD, if true --> all PlotInfo will be set to false and NO PLOTS will be drawn. 
PlotInfo.NoPlots                                         = NoPlotsInfo;  % MASTER FIELD, if true --> all PlotInfo will be set to false and NO PLOTS will be drawn. 

%% Plots related to tracks association with arenas
PlotInfo_Default.TracksHistogram                         = true;
PlotInfo_Default.TracksAssociationToArenas               = true;

%% Plots related to locomotion  
PlotInfo_Default.Path                                    = false;
PlotInfo_Default.VelocityCartesian                       = false;
PlotInfo_Default.VelocityPolar                           = false;
PlotInfo_Default.AccelerationCartesian                   = false;
PlotInfo_Default.AccelerationPolar                       = false;
PlotInfo_Default.ChangesedInVelocityPolar                = false;
PlotInfo_Default.Path_Velcoty_Acceleration_AllCartesian  = false;
PlotInfo_Default.Path_Velcoty_Acceleration_Polar         = false;
PlotInfo_Default.Path_Velcoty_ChangesInVelocity_Polar    = false;
PlotInfo_Default.histograms                              = false;
PlotInfo_Default.Velocity_vs_ChangesInVelocity           = false;
PlotInfo_Default.HeadTailPath                            = true;

%% Movies related to locomotion  
% PlotInfo_Default.MOVIE_using_ChangesInVelocity           = true;                 % If true: Create movie of the centroid path with the locomotion analysis   
PlotInfo_Default.MOVIE_using_ChangesInVelocity           = false;                  % If true: Create movie of the centroid path with the locomotion analysis   
% PlotInfo_Default.MOVIE_SingleTrack_WithWorms             = true;                 % If true: Create movie of the centroid path with the locomotion analysis on the worm Movie   
PlotInfo_Default.MOVIE_SingleTrack_WithWorms             = false;                  % If true: Create movie of the centroid path with the locomotion analysis on the worm Movie    

%% Plots related to morphology
PlotInfo_Default.Morphology_All                          = false;
PlotInfo_Default.Eccentricity_Perimeter                  = false;

% Internal choice: When ploting the results in 'SegmentSingleTrack_SmallArena_SL_vXX' , which plot configurations to show?   
%   PlotConfiguration = 0; % No Plots
%   PlotConfiguration = 1; % 6 plots
%   PlotConfiguration = 2;   % 4 plots
PlotInfo_Default.PlotConfiguration                        = 0;                  
                 
%% Set all non-defined fields to the default fields OR to false, depending on MASTER field
PlotFields = fieldnames(PlotInfo_Default);
for f_num  = 1:length(PlotFields)
    field_name = PlotFields{f_num};
    if PlotInfo.NoPlots
        PlotInfo.(field_name) = false;
    else
        PlotInfo.(field_name) = PlotInfo_Default.(field_name);
    end        
end

return




