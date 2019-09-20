function [Tracks, File, TracksStats, Data] = SegmentBehavior_ImagingInGradient_v03(File, Tracks, TimeOfGradientStartValveSwitch, TimeOfGradientEndValveSwitch)

% September 2019, Sagi Levy

%% Inputs
% Input variables can be loaded by: % load('MoviePath\MovieName_AfterNeuronPositionManualCorrection','File','Tracks')  
% TimeOfGradientStartValveSwitch, TimeOfGradientEndValveSwitch >>  elapsed time in which gradient valves were turned on and off, respectively 
 
CorrectBasedOnNeuronCoordinatesAndMask = true;
load(File.BackgroundFile,'Mask');
LEN_Tracks = length(Tracks);

%% Important note about Path definition change!!!
%% The Tracks(tr).Path variables has a switched X-Y coordinate relative to all other parameters  !!!  This is corrected in this function. Please note that when running GUI tracker to observe tracking results! 
disp(['The Tracks(tr).Path variables has a switched X-Y coordinate relative to all other parameters  !!!',char(10),'It is corrected in this function but it needs to be addressed in the tracking function! ', char(10)]);
for tr = 1:LEN_Tracks         
    Tracks(tr).Path = Tracks(tr).Path(:,[2 1]);   % FLIPPING TO FIT THE REAL X-Y COORDINATE!!!    
                                                  % Note that Head, Tail and Neuron coordinates don't need to be flipped to match the same consensus    
end

%% Initialize Settings parameters 
Settings           = File.VariablesInformation.SegmentationSettings;           % Call External Function with segmentation parameters. 
Settings.FrameRate = File.FrameRate; 
Settings.PixelSize = File.PixelSize;
Settings.AllowVeryLongTracks = true;  % force allowing long tracks
NumOfFrames        = File.NumberOfFrames; 
ArenaID            = 1;

%% LOOP OVER ALL TRACKS AND CALCULATE LOCOMOTION AND MORPHOLOGY DYNAMICS
disp(['Movie name: ',File.MovieName]); 
disp([datestr(now),' -- Analyzing locomotion data in ',num2str(LEN_Tracks),' tracks']);
h_waitbar  = waitbar(0,['Analyzing locomotion data in ',num2str(LEN_Tracks),' tracks']);

% General locomotion and morphology properties
Sizes                               = [Tracks.Size];
TracksStats.MedianSize              = nanmedian(single(Sizes));
Midlines                            = [Tracks.Midline];
ReliableSkeletons                   = [Midlines.FlagTrueIfReliable];
SkeletonLengths                     = [Tracks.SkeletonLength];
PerimeterLengths                    = [Tracks.PerimeterLength];
TracksStats.MedianSkeleton          = nanmedian(single(SkeletonLengths(ReliableSkeletons)));
TracksStats.MedianPerimeter         = nanmedian(single(PerimeterLengths(ReliableSkeletons)));
TracksStats.MedianPerimeterOverSize = nanmedian(single(PerimeterLengths(ReliableSkeletons))./single(Sizes(ReliableSkeletons)));

% % Convert 'Tracks' to fit the movie number of frames
Tracks = ExtendTracksVectorsToMovieLength(Tracks, NumOfFrames);

%% Correct Head/Tail and OOB based on Neuron Coordinates, and find the Limit For Gradient axis 
if CorrectBasedOnNeuronCoordinatesAndMask
    plotme = false;
%     plotme = true;
    Tracks = CorrectHeadTailSegmentationBasedOnNeuronCoordinates(Tracks, plotme);
    [Tracks, GradientAxisLimits, FlowAxisLimits] = CorrectOOBBasedOnNeuronCoordinatesAndMask(Tracks, Mask, plotme);
else
    GradientAxisLimits = [];
    FlowAxisLimits     = [];
end

for tr = 1:LEN_Tracks         
    waitbar(tr/LEN_Tracks , h_waitbar);  

    % Calculate locomotion and morphology dynamics using EXTERNAL FUNCTIONS 
    TrackLocomotionData = AnalyzeLocomotionProperties_ImagingSetup_v2(Tracks(tr), Settings);               % Tracks(tr).Path =  Worm Coordinates XY(t), time in units of frames   
    TrackMorphologyData = AnalyzeMorphologyProperties_ImagingSetup_v3(Tracks(tr), Settings);           
   
    %% Assign Vectors to Tracks(tr) 
    % Locomotion
    Tracks(tr).Speed                                 = TrackLocomotionData.Velocity.Amplitude';         % Speed [pixels/frame] == radial velocity == tangential velocity == velocity amplitude
    Tracks(tr).AngularDisplacement                   = TrackLocomotionData.Velocity.Angle';             % Range is [-180 180]. Angular displacement [degrees] == velocity angle == change in radial velocity in direction of its normal.    
    Tracks(tr).Velocity_X                            = TrackLocomotionData.Velocity.X';                 % X-axis velocity amplitude [pixels/frame] 
    Tracks(tr).Velocity_Y                            = TrackLocomotionData.Velocity.Y';                 % Y-axis velocity amplitude [pixels/frame] 
    Tracks(tr).AngularDisplacement_ContinuousAngle   = TrackLocomotionData.Velocity.ContinuousAngle';   % CONTINUOUS range. Angular displacement [degrees] == velocity angle == change in radial velocity in direction of its normal.        
    Tracks(tr).Acceleration                          = TrackLocomotionData.ChangeInVelocity.Amplitude'; % Acceleration [pixels/frame^2] == radial Acceleration == tangential acceleration 
    Tracks(tr).Jerk                                  = TrackLocomotionData.Jerk.Amplitude';             % Jerk [pixels/frame^3] == radial jerk == tangential jerk 
    Tracks(tr).AngularVelocity                       = TrackLocomotionData.ChangeInVelocity.Angle';     % AngularVelocity [degrees/frame] == derivative of the angular displacement
    Tracks(tr).Path_smoothed_X_coordinate            = TrackLocomotionData.Path.X';
    Tracks(tr).Path_smoothed_Y_coordinate            = TrackLocomotionData.Path.Y';  
    Tracks(tr).Curvature                             = Tracks(tr).AngularVelocity ./ Tracks(tr).Speed;  % Curvature [degree/pixel] == angular velocity/ speed.   
    
    % Morphology
    Tracks(tr).PerimeterOverSize                     = TrackMorphologyData.PerimeterOverSize;
    Tracks(tr).NormalizedSize                        = TrackMorphologyData.NormalizedSize;
    Tracks(tr).NormalizedSkeletonLength              = TrackMorphologyData.NormalizedSkeletonLength;
    Tracks(tr).NormalizedPerimeterLength             = TrackMorphologyData.NormalizedPerimeterLength;
    Tracks(tr).NormalizedPerimeterOverSize           = TrackMorphologyData.NormalizedPerimeterOverSize;
    
    if isfield(Tracks(tr),'Midline')
        midline_first_angle            = mod(180- Tracks(tr).Midline.AngleFirstPoint,360); 
        HeadDirection                  = mod(-midline_first_angle+180,360)-180;                     % 'HeadDirection' has the same angles definition as the Angluar Displacement !!
        HeadDirection_ContinuousAngle  = FindContinuousAngle(HeadDirection,'degrees');              % FindContinuousAngle(Angle,'degrees',[],true); % plotfigure = true          
        Tracks(tr).HeadDirection       = HeadDirection;                                             % Angle [degree]
        Tracks(tr).HeadDirectionChange = single([NaN diff(HeadDirection_ContinuousAngle)]);         % Angle derivative [degree/frame]
        
        HeadAngleRelativeToDirection   = Tracks(tr).HeadDirection - Tracks(tr).AngularDisplacement; 
        HeadAngleRelativeToDirection   = asin(sin(HeadAngleRelativeToDirection/180*pi))*180/pi;     % Angle of the head relative to movement direction between [-90 to 90] degrees == [left right]
        Tracks(tr).HeadAngleRelativeToDirection = single(HeadAngleRelativeToDirection);             % Angle [degree]
    end

    %% Calculate Stimulus and FlowAxis location: from (-1) to (1)
    % For Stimulus:  -1 corresponds to bottom if the field of view (High values of Dimension 1)   
    % For Flow:      -1 corresponds to near the outlet barriers    (High values of Dimension 2)   
    
    Tracks = CalculateRelativeCoordinates(Tracks, GradientAxisLimits, FlowAxisLimits);      
end       

%% Segment behavior properties
disp([datestr(now),' -- Segmenting behavior in ',num2str(LEN_Tracks),' tracks']);
waitbar(0,h_waitbar,['Segmenting behavior in ',num2str(LEN_Tracks),' tracks']);
warning off

for tr = 1:LEN_Tracks  
    waitbar(tr/LEN_Tracks , h_waitbar);      
   [BehaviorCode, Segments_HighLevelBehavior, Segments_LowLevelBehavior] = SegmentBehavior_SingleTrack_ImagingSetup_v2(Tracks(tr), tr, Settings);    
%     PlotLocomotion_And_Behavior(Tracks(tr), BehaviorCode);
    
    Tracks(tr).Segments_LowLevelBehavior    = Segments_LowLevelBehavior;
    Tracks(tr).Segments_HighLevelBehavior   = Segments_HighLevelBehavior;
    Tracks(tr).BehaviorCode_Structure       = BehaviorCode;       
    Tracks(tr).BehaviorCode_HighLevel       = BehaviorCode.HighLevel.BehaviorVector;       
    Tracks(tr).BehaviorCode_LowLevel        = BehaviorCode.LowLevel.BehaviorVector;      
end
warning on;

TracksStats.BehaviorCode.LowLevel_CodeNames    = BehaviorCode.LowLevel.BehaviorCodeName;
TracksStats.BehaviorCode.LowLevel_CodeNumbers  = BehaviorCode.LowLevel.BehaviorCodeNumbers;
TracksStats.BehaviorCode.HighLevel_CodeNames   = BehaviorCode.HighLevel.BehaviorCodeName;
TracksStats.BehaviorCode.HighLevel_CodeNumbers = BehaviorCode.HighLevel.BehaviorCodeNumbers;

%% Generate data structure
%  File containing 'Data' structure with subfields of interesting locomotion/morphology/behavior features of worms. 'Data' STRUCTURE DOESN'T CONTAIN ANY MOVIE DISPLAY FEATURES such as perimeter coordinates. 
%           example for subfields:  Speed, AngularVelocity, Eccentricity, BehaviorCodeLowLevel, BehaviorCodeHighLevel, BehaviorProbLowLevel, BehaviorProbHighLevel, HeadDirection ...   
%  most matrices:                  data.(subfield_name) are (NxM) matrices; N=number of worms, M=number of frames. 
%  Behavior probability matrices:  data.(subfield_name) are (BxM) matrices; B=behavior code number, M=number of frames.  
%  Unavailable data-points will have a NaN value ONLY FOR SINGLE PRECISION MATRICES !!!! and '0' for logical, uint8 and uint16 variables 

[Data, File, Stats] = GenerateDataStructure (Tracks, File, TracksStats.BehaviorCode);

%% Update File structure and Save
load([File.TrackFile(1:end-4),'_DyePatterns.mat'],'FlowDelay')
load(File.BackgroundFile,'background','VignettingPattern')
File.GradientStartFrame                           = (TimeOfGradientStartValveSwitch*60+FlowDelay.DelayBetweenValveAndArenas)*File.FrameRate; 
File.GradientEndFrame                             = (TimeOfGradientEndValveSwitch*60+FlowDelay.DelayBetweenValveAndArenas)*File.FrameRate; 
File.GradientRiseTime_InitiationTo50Percent_InSec = FlowDelay.RiseTime_InitiationTo50Percent;
File.GradientRiseTime_InitiationTo90Percent_InSec = FlowDelay.RiseTime_InitiationTo90Percent;
File.FullGradientTime_InSec                       = File.GradientStartFrame/File.FrameRate+File.GradientRiseTime_InitiationTo90Percent_InSec;
File.FlowDelay                                    = FlowDelay;

disp([datestr(now),'Saving Data Matrices file of arena ',num2str(ArenaID),' of the movie: ',File.MovieName]);
save(File.FileNames.DataMatrices{ArenaID},'Data','File','background','VignettingPattern','TracksStats','FlowDelay','-v7.3');        
disp([datestr(now),'Saved']);

return

%% Internal functions
function ContinuousAngles = FindContinuousAngle(Angles, Mode, RangeEdges, plotfigures)      
% Inputs
% Angles             A vector of angles that may be trancated due to edge singularities (e.g. 180--> -180, etc)
% Mode               Optional string: 'degrees' if the range is [-180 180] (default) or radians if range is between [-pi pi].  
% RangeEdges         Optional vector of length 2. Manually specify the edges and ignore 'Mode'. 
% plotfigures        Optional logical. default = false;

% Output:
% ContinuousAngles   A vector of Angles without singularities, range may vary beyond total of 360 degrees (or 2*pi radians) 

%% find AnglesEdge
if ~exist('Mode','var')
    Mode = 'degrees'; % [-180 180] 
end
if ~exist('plotfigures','var')
    plotfigures = false; 
end
if strcmpi(Mode,'degrees')
    AnglesEdge = [-180 180];
elseif strcmpi(Mode,'radians')
    AnglesEdge = [-pi pi];
else
    if exist('RangeEdges','var')
        AnglesEdge = RangeEdges;
    else
        disp('Specify Mode or RangeEdges. Aborting ''FindContinuousAngle''... ');
        return
    end
end
    

%% Find when angles are changes by more than half the range: fix discontinuity by correcting all following angles. Repeat until end of vector
ContinuousAngles = Angles;
AngleLimit       = AnglesEdge(2);
for step = 1:1e4
    index       = find(abs(diff(ContinuousAngles))> AngleLimit ,1,'first'); % incontinuity between index and index+1
    if isempty(index)
        break
    end
    CorrectSign = sign(ContinuousAngles(index)-ContinuousAngles(index+1));
    ContinuousAngles((index+1):end) = ContinuousAngles((index+1):end) + 2*AngleLimit * CorrectSign;
end
if ~isempty(index)
    disp('Warning!! The function ''FindContinuousAngle'' was aborted before Angles were properly fixed')
end

if plotfigures
    figure('position',get(0,'screensize')); 
    subplot(2,2,1);
    plot(Angles,'b*-'); hold on;                  title ( 'Angles');
    line (get(gca,'xlim'),[AnglesEdge(1) AnglesEdge(1)],'linestyle',':','color','k');
    line (get(gca,'xlim'),[AnglesEdge(2) AnglesEdge(2)],'linestyle',':','color','k');
    subplot(2,2,3);
    plot(diff(Angles),'r*-');                     title ( 'diff(Angles)');
%     subplot(5,1,3); 
%     plot(abs(diff(Angles))>AnglesEdge(2),'r*-');  title ( 'abs(diff(Angles))>AnglesEdge(2)');
    subplot(2,2,2); 
    plot(ContinuousAngles,'b*-'); hold on;        title ( 'ContinuousAngles');
    subplot(2,2,4); 
    plot(diff(ContinuousAngles),'r*-'); hold on;  title ( 'diff(ContinuousAngles)');
end

return

function PlotLocomotion_And_Behavior(CurrentTrack, BehaviorCode)

Path         = CurrentTrack.Path;
Head         = squeeze(CurrentTrack.HeadTail(:,1,:));
Tail         = squeeze(CurrentTrack.HeadTail(:,2,:));
MidlineX     = CurrentTrack.Midline.X_coordinates_short;
MidlineY     = CurrentTrack.Midline.Y_coordinates_short;
PerimeterX   = CurrentTrack.WormPerimeter.Xcoordinate;
PerimeterY   = CurrentTrack.WormPerimeter.Ycoordinate;
FlagReliable = CurrentTrack.Midline.FlagTrueIfReliable;

limits = [ min([Head(:,1)' Tail(:,1)']),  max([Head(:,1)' Tail(:,1)']), min([Head(:,2)' Tail(:,2)']),  max([Head(:,2)' Tail(:,2)'])]; 
Xwidth = limits(2)-limits(1);
Ywidth = limits(4)-limits(3);
SCSZ = get(0,'ScreenSize');
ScreenHeightToWidth = SCSZ(4)/SCSZ(3);
PathHeightToWidth   = Ywidth/Xwidth;
if ScreenHeightToWidth > PathHeightToWidth
    SCSZ(4)=SCSZ(3)*PathHeightToWidth;
else
    SCSZ(3)=SCSZ(4)/PathHeightToWidth;
end   

if exist('BehaviorCode','var')
    Logical_Vectors.OutOfBounds     = BehaviorCode.LowLevel.BehaviorVector==0;
    Logical_Vectors.Forward         = BehaviorCode.LowLevel.BehaviorVector==1;
    Logical_Vectors.CurvingForward  = BehaviorCode.LowLevel.BehaviorVector==2;
    Logical_Vectors.Pause           = BehaviorCode.LowLevel.BehaviorVector==3;
    Logical_Vectors.Reverse         = BehaviorCode.LowLevel.BehaviorVector==4;
    Logical_Vectors.Omega           = BehaviorCode.LowLevel.BehaviorVector==5;
    Logical_Vectors.SharpTurns      = BehaviorCode.LowLevel.BehaviorVector==6;
 
    Logical_Vectors_HighLevel.OutOfBounds     = BehaviorCode.HighLevel.BehaviorVector==0;
    Logical_Vectors_HighLevel.Forward         = BehaviorCode.HighLevel.BehaviorVector==1;
    Logical_Vectors_HighLevel.CurvingForward  = BehaviorCode.HighLevel.BehaviorVector==2;
    Logical_Vectors_HighLevel.Pause           = BehaviorCode.HighLevel.BehaviorVector==3;
    Logical_Vectors_HighLevel.Reversal        = BehaviorCode.HighLevel.BehaviorVector==4;
    Logical_Vectors_HighLevel.OmegaWithoutPirouette = BehaviorCode.HighLevel.BehaviorVector==5;
    Logical_Vectors_HighLevel.SharpTurns      = BehaviorCode.HighLevel.BehaviorVector==6;
    Logical_Vectors_HighLevel.Pirouette       = (BehaviorCode.HighLevel.BehaviorVector==7)+(BehaviorCode.LowLevel.BehaviorVector==8);

    f_behavior = figure('position',[50 floor(SCSZ(4)/2) SCSZ(3)-100 floor(SCSZ(4)/2)-100],'name','Low and High Level Behavior'); 
    plot(Logical_Vectors.Forward,       'b.','markersize',6); hold on;  % Not extended
    plot(Logical_Vectors.CurvingForward,'bo','markersize',6,'markerfacecolor','b'); hold on;  % Not extended
    plot(Logical_Vectors.Reverse,     'r.','markersize',6); hold on;  % Not extended
    plot(Logical_Vectors.SharpTurns,  'ks','markerfacecolor','k','markersize',10); hold on;
    plot(Logical_Vectors.Omega ,      'gv','markerfacecolor','g'); 
    plot(Logical_Vectors.Pause ,      'y.'); hold on;
    plot(Logical_Vectors.OutOfBounds ,'ms','markerfacecolor','m'); hold on;
    plot(1.01*Logical_Vectors_HighLevel.Pirouette,              'rs','markerfacecolor','g');
    plot(1.01*Logical_Vectors_HighLevel.Reversal,               'r.');
    plot(1.01*Logical_Vectors_HighLevel.OmegaWithoutPirouette,  'gv','markerfacecolor','g');
    plot(1.01*Logical_Vectors_HighLevel.SharpTurns,             'ks','markerfacecolor','k');
    plot(1.01*Logical_Vectors_HighLevel.Forward,                'b.');
    plot(1.01*Logical_Vectors_HighLevel.Pause,                  'y.');
    plot(1.01*Logical_Vectors_HighLevel.OutOfBounds,            'ms','markerfacecolor','m'); hold on;
    ylimits = [0.955 1.06];
    ylim(ylimits)
    xlim([0 CurrentTrack.TrackLength]);
    behavior_plot = gca;
end

f_worm = figure('position',SCSZ);
for frame_ind = 1:CurrentTrack.TrackLength
    figure(f_worm);
    plot(Path(frame_ind,1),Path(frame_ind,2),'r*','markersize',25);  hold on;
    plot(PerimeterX{frame_ind},PerimeterY{frame_ind},'k:'); 
    if FlagReliable(frame_ind)
        plot(MidlineX{frame_ind},  MidlineY{frame_ind},'k');  hold on;
        plot(Head(frame_ind,1),Head(frame_ind,2),'gx','markersize',15);
        plot(Tail(frame_ind,1),Tail(frame_ind,2),'go','markersize',15);
    end
    axis(limits);
    title(frame_ind);
    hold off;
    
    if exist('BehaviorCode','var')
        figure(f_behavior);
        line_h = line([frame_ind frame_ind],ylimits,'color','k','parent',behavior_plot);                               
    end       
%     pause
    pause(0.1)
    delete(line_h);
end
    

return

function [Data, File, Stats] = GenerateDataStructure (Tracks, File, BehaviorCode)

NumOfWorms   = length(Tracks);
NumOfFrames  = File.NumberOfFrames;

NaN_Values   = true(NumOfWorms, NumOfFrames);
for tr_ind = 1:NumOfWorms    
    NaN_Values(tr_ind,Tracks(tr_ind).Frames) = false;    
end
Data.NaN    = NaN_Values;
Real_Values = ~NaN_Values;
Stats.NumberOfWorms = uint8(sum(Real_Values,1));

for tr_ind = 1:NumOfWorms    
    Tracks(tr_ind).Coordinates_X                 = Tracks(tr_ind).Path(:,1);    
    Tracks(tr_ind).Coordinates_Y                 = Tracks(tr_ind).Path(:,2);    
    Tracks(tr_ind).Coordinates_X_Smoothed        = Tracks(tr_ind).Path_smoothed_X_coordinate;    
    Tracks(tr_ind).Coordinates_Y_Smoothed        = Tracks(tr_ind).Path_smoothed_Y_coordinate;   
    Tracks(tr_ind).OutOfBounds_WithCollisions    = Tracks(tr_ind).OutOfBounds;  
    Tracks(tr_ind).FlagTrueIfMidlineIsReliable   = Tracks(tr_ind).Midline.FlagTrueIfReliable;      
    Tracks(tr_ind).NumOfBends_HighRes            = Tracks(tr_ind).Midline.NumOfBends_HighRes;      
    Tracks(tr_ind).NumOfBends_LowRes             = Tracks(tr_ind).Midline.NumOfBends_LowRes;      
end
Tracks = rmfield(Tracks,{'Path','Path_smoothed_X_coordinate','Path_smoothed_Y_coordinate','OutOfBounds'});
for tr_ind = 1:NumOfWorms    
    Tracks(tr_ind).HeadCoordinates_X        = Tracks(tr_ind).Head(:,1);    
    Tracks(tr_ind).HeadCoordinates_Y        = Tracks(tr_ind).Head(:,2);    
    Tracks(tr_ind).TailCoordinates_X        = Tracks(tr_ind).Tail(:,1);    
    Tracks(tr_ind).TailCoordinates_Y        = Tracks(tr_ind).Tail(:,2);   
    Tracks(tr_ind).NeuronCoordinates_X      = Tracks(tr_ind).Neuron.CoordinatesMatrix(:,1);    
    Tracks(tr_ind).NeuronCoordinates_Y      = Tracks(tr_ind).Neuron.CoordinatesMatrix(:,2);   
    Tracks(tr_ind).NeuronValue              = Tracks(tr_ind).Neuron.Value;   
%     Tracks(tr_ind).deltaFOverF              = Tracks(tr_ind).Neuron.deltaFOverF;                 % Not used, this will be redefined in later functions 
%     Tracks(tr_ind).deltaFOverF_Interpolated = Tracks(tr_ind).Neuron.deltaFOverF_Interpolated;   
%     Tracks(tr_ind).deltaFOverF_Filtered     = Tracks(tr_ind).Neuron.deltaFOverF_Filtered;   
end
Tracks = rmfield(Tracks,{'Head','Tail'});

%% Unit conversions !!!!
% Speed, Velocity_X, Velocity_Y   [pixels/frame]   --> [micrometer/second] 
% Acceleration                    [pixels/frame^2] --> [micrometer/second^2] 
% Jerk                            [pixels/frame^3] --> [micrometer/second^3]
% Curvature                       [angle/pixel]    --> [angle/micrometer]
% AngularVelocity                 [angle/frame]    --> [angle/second] 
%                                  keep [angle/frame] as a new field: 'AngularVelocity_AnglePerFrame'  
% HeadDirectionChange             [angle/frame]    --> [angle/second] 
%                                  keep [angle/frame] as a new field: 'HeadDirectionChange_AnglePerFrame'  
FramesPerSecond    = File.FrameRate ;
MicroMeterPerPixel = 1000/File.PixelSize;   % File.PixelSize == how many pixels in 1mm ....

for tr_ind = 1:NumOfWorms    
    Tracks(tr_ind).Speed        = Tracks(tr_ind).Speed        * MicroMeterPerPixel * FramesPerSecond   ;  % [pixels/frame]   --> [micrometer/second]   
    Tracks(tr_ind).Velocity_X   = Tracks(tr_ind).Velocity_X   * MicroMeterPerPixel * FramesPerSecond   ;  % [pixels/frame]   --> [micrometer/second]   
    Tracks(tr_ind).Velocity_Y   = Tracks(tr_ind).Velocity_Y   * MicroMeterPerPixel * FramesPerSecond   ;  % [pixels/frame]   --> [micrometer/second]   
    Tracks(tr_ind).Acceleration = Tracks(tr_ind).Acceleration * MicroMeterPerPixel * FramesPerSecond^2 ;  % [pixels/frame^2] --> [micrometer/second^2] 
    Tracks(tr_ind).Jerk         = Tracks(tr_ind).Jerk         * MicroMeterPerPixel * FramesPerSecond^3 ;  % [pixels/frame^3] --> [micrometer/second^3] 
    Tracks(tr_ind).Curvature    = Tracks(tr_ind).Curvature    / MicroMeterPerPixel ;                      % [angle/pixel]    --> [angle/micrometer]   

    Tracks(tr_ind).AngularVelocity_AnglePerFrame     = Tracks(tr_ind).AngularVelocity ;           
    Tracks(tr_ind).AngularVelocity                   = Tracks(tr_ind).AngularVelocity      * FramesPerSecond ;    % [angle/frame]   --> [angle/second]   
    Tracks(tr_ind).HeadDirectionChange_AnglePerFrame = Tracks(tr_ind).HeadDirectionChange ;           
    Tracks(tr_ind).HeadDirectionChange               = Tracks(tr_ind).HeadDirectionChange  * FramesPerSecond ;    % [angle/frame]   --> [angle/second]   
end

%% Each variable will be modified, if needed, to a single/uint16,uint8 or logical matrix. Most variables are already in the right format.        
SinglePrecisionFields  = {'Coordinates_X','Coordinates_Y','Coordinates_X_Smoothed','Coordinates_Y_Smoothed',...
                          'Eccentricity','MajorAxes','MinorAxes','Orientation',...
                          'Speed','AngularDisplacement','Velocity_X','Velocity_Y','AngularDisplacement_ContinuousAngle','Acceleration','Jerk','AngularVelocity','Curvature',...
                          'PerimeterOverSize','NormalizedSize','NormalizedSkeletonLength','NormalizedPerimeterLength','NormalizedPerimeterOverSize',...
                          'HeadDirection','HeadDirectionChange','AngularVelocity_AnglePerFrame','HeadDirectionChange_AnglePerFrame','HeadAngleRelativeToDirection',...
                          'HeadCoordinates_X','HeadCoordinates_Y','TailCoordinates_X','TailCoordinates_Y',...         
                          'NeuronCoordinates_X','NeuronCoordinates_Y','NeuronValue',...
                          'Stimulus_Head','Stimulus_Neuron','Stimulus_Centroid','FlowAxis_Head','FlowAxis_Neuron','FlowAxis_Centroid'};
Uint16Fields           = {'Size','SkeletonLength','PerimeterLength'};      
Uint8Fields            = {'NumOfBends_HighRes','NumOfBends_LowRes','BehaviorCode_HighLevel','BehaviorCode_LowLevel'};
LogicalFields          = {'OutOfBounds_WithCollisions','FlagTrueIfMidlineIsReliable','HeadIsSqueezed','PossibleCollision'};  

% Initialization
for f_ind = 1:length(SinglePrecisionFields)
    fieldname = SinglePrecisionFields{f_ind};
    Data.(fieldname) = zeros(NumOfWorms,NumOfFrames,'single')*NaN;
end
for f_ind = 1:length(Uint16Fields)
    fieldname = Uint16Fields{f_ind};
    Data.(fieldname) = zeros(NumOfWorms,NumOfFrames,'uint16');
end
for f_ind = 1:length(Uint8Fields)
    fieldname = Uint8Fields{f_ind};
    Data.(fieldname) = zeros(NumOfWorms,NumOfFrames,'uint8');
end
for f_ind = 1:length(LogicalFields)
    fieldname = LogicalFields{f_ind};
    Data.(fieldname) = false(NumOfWorms,NumOfFrames);
end
Data.FrameNumber = zeros(NumOfWorms,NumOfFrames,'single')*NaN;

% Assigning Tracks information
for tr_ind = 1:NumOfWorms
    CurrentTrackFrames = Real_Values(tr_ind,:);
    
    for f_ind = 1:length(SinglePrecisionFields)
        fieldname = SinglePrecisionFields{f_ind};
        Data.(fieldname)(tr_ind,CurrentTrackFrames) = single(Tracks(tr_ind).(fieldname));
    end
    for f_ind = 1:length(Uint16Fields)
        fieldname = Uint16Fields{f_ind};
        Data.(fieldname)(tr_ind,CurrentTrackFrames) = uint16(Tracks(tr_ind).(fieldname));
    end
    for f_ind = 1:length(Uint8Fields)
        fieldname = Uint8Fields{f_ind};
        Data.(fieldname)(tr_ind,CurrentTrackFrames) = uint8(Tracks(tr_ind).(fieldname));
    end
    for f_ind = 1:length(LogicalFields)
        fieldname = LogicalFields{f_ind};
        Data.(fieldname)(tr_ind,CurrentTrackFrames) = logical(Tracks(tr_ind).(fieldname));
    end 
    FramesIndices = find(CurrentTrackFrames);
    Data.FrameNumber(tr_ind,CurrentTrackFrames) = single(FramesIndices);
end

% Interpolation_Factor = 1e-4;
% BehaviorCodeNumbers  = BehaviorCode.LowLevel_CodeNumbers;  
% BehaviorCodeMAT      = Data.BehaviorCode_LowLevel;
% [BehaviorProbability, BehaviorProbability_Smoothed, DetectionAndSegmentationProbabilityVector] = CalculateBehaviorProbability (BehaviorCodeMAT, BehaviorCodeNumbers, Interpolation_Factor);
% Data.BehaviorProbability_LowLevel                        = BehaviorProbability;
% Data.BehaviorProbability_Smoothed_LowLevel               = BehaviorProbability_Smoothed;
% STATS.DetectionAndSegmentationProbabilityVector_LowLevel = DetectionAndSegmentationProbabilityVector;
% 
% Interpolation_Factor = 1e-4;
% BehaviorCodeNumbers  = BehaviorCode.HighLevel_CodeNumbers;  
% BehaviorCodeMAT      = Data.BehaviorCode_HighLevel;
% [BehaviorProbability, BehaviorProbability_Smoothed, DetectionAndSegmentationProbabilityVector] = CalculateBehaviorProbability (BehaviorCodeMAT, BehaviorCodeNumbers, Interpolation_Factor);
% Data.BehaviorProbability_HighLevel                        = BehaviorProbability;
% Data.BehaviorProbability_Smoothed_HighLevel               = BehaviorProbability_Smoothed;
% STATS.DetectionAndSegmentationProbabilityVector_HighLevel = DetectionAndSegmentationProbabilityVector;

File.FileNames.DataMatrices                    = cell(1,File.NumArenas);
for ar=1:File.NumArenas
    File.FileNames.DataMatrices{ar} = [File.TrackFile(1:end-4),'_DataMatrices_Arena',num2str(ar),'.mat'];
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

DetectionAndSegmentationProbabilityVector = NumberOfSegmentedTracks/NumOfWorms;
% Probability of each (not out of bound) behavior relative to the total number of tracks that were segmented for behaviors (detected + not out of bound) 
BehaviorProb          = behhist_SegmentedTracks ./ repmat(NumberOfSegmentedTracks,size(behhist_SegmentedTracks,1),1);   

BehaviorProb_Smoothed = zeros(size(BehaviorProb));
for row = 1:size(BehaviorProb,1)
    Vec                          = BehaviorProb(row,:);
    x                            = 1:length(Vec);
    BehaviorProb_Smoothed(row,:) = csaps(x,Vec,Interpolation_Factor, x);  
end

return

function Tracks = ExtendTracksVectorsToMovieLength (Tracks_in, MovieLength) 
% Tracks will be similar to Tracks_in, but all double precision variables will be converted to single precision variables.  
Tracks     = Tracks_in;
FIELDNAMES = fieldnames(Tracks);

RelevantFieldsToCorrect = false(1,length(FIELDNAMES));
tr=1;
for f_ind = 1:length(FIELDNAMES)
    fieldname = FIELDNAMES{f_ind};
    CurrentNumberOfFrames = length(Tracks(tr).Frames);
    RelevantFieldsToCorrect(f_ind) = ismember(CurrentNumberOfFrames,size(Tracks(tr).(fieldname)));
end
FIELDNAMES_FirstField    = FIELDNAMES(RelevantFieldsToCorrect);
FIELDNAMES_midline       = fieldnames(Tracks(tr).Midline);
FIELDNAMES_WormPerimeter = fieldnames(Tracks(tr).WormPerimeter);

for tr=1:length(Tracks)
    Frames                                         = Tracks(tr).Frames;
    CurrentNumberOfFrames                          = length(Frames);
    Tracks(tr).Frames_BeforeExtendingToMovieLength = Frames;
    Tracks(tr).Frames                              = 1:MovieLength;
    NewOutOfBounds                                 = true(1,MovieLength);
    NewOutOfBounds(Frames)                         = false;
    for f_ind = 1:length(FIELDNAMES_FirstField)
        fieldname  = FIELDNAMES_FirstField{f_ind};
        if strcmpi(fieldname,'Frames')
            continue
        end
        CurrentVec = Tracks(tr).(fieldname);
        if length(size(CurrentVec))<3            
            Y=whos('CurrentVec'); CurrentClass = Y.class;
            NewSize = Y.size; NewSize(NewSize==CurrentNumberOfFrames)=MovieLength;
            if islogical(CurrentVec)            
                NewVec = false(NewSize);
            elseif iscell(CurrentVec)
                NewVec = cell(NewSize);
            else
                NewVec = zeros(NewSize,CurrentClass)*NaN; % NaN is applied only to 'single', uint8 & uint16 still have zeros!!
            end
            if min(size(CurrentVec))==1
                NewVec(Frames) = CurrentVec;
            else
                NewVec(Frames,:) = CurrentVec;
            end
            
        elseif strcmpi(fieldname,'HeadTail')   % [Framesx2x2 single]
            NewVec = zeros(size(CurrentVec),'single')*NaN;
            NewVec(Frames,:,:) = CurrentVec;            
        else
            disp(fieldname,' -- unnkown format')
        end
        Tracks(tr).(fieldname) = NewVec;                        
    end  
    
    % Midline fields
    for f_ind = 1:length(FIELDNAMES_midline)
        fieldname  = FIELDNAMES_midline{f_ind};
        CurrentVec = Tracks(tr).Midline.(fieldname);
        Y=whos('CurrentVec'); CurrentClass = Y.class;
        NewSize = Y.size; NewSize(NewSize==CurrentNumberOfFrames)=MovieLength;
        if islogical(CurrentVec)            
            NewVec = false(NewSize);
        elseif iscell(CurrentVec)
            NewVec = cell(NewSize);
        else
            NewVec = zeros(NewSize,CurrentClass)*NaN; % NaN is applied only to 'single', uint8 & uint16 still have zeros!!
        end
        if min(size(CurrentVec))==1
            NewVec(Frames) = CurrentVec;
        else
            NewVec(Frames,:) = CurrentVec;
        end
        Tracks(tr).Midline.(fieldname) = NewVec;    
    end
    
     % WormPerimeter fields
    for f_ind = 1:length(FIELDNAMES_WormPerimeter)
        fieldname  = FIELDNAMES_WormPerimeter{f_ind};
        CurrentVec = Tracks(tr).WormPerimeter.(fieldname);
        Y=whos('CurrentVec'); CurrentClass = Y.class;
        NewSize = Y.size; NewSize(NewSize==CurrentNumberOfFrames)=MovieLength;
        if islogical(CurrentVec)            
            NewVec = false(NewSize);
        elseif iscell(CurrentVec)
            NewVec = cell(NewSize);
        else
            NewVec = zeros(NewSize,CurrentClass)*NaN; % NaN is applied only to 'single', uint8 & uint16 still have zeros!!
        end
        if min(size(CurrentVec))==1
            NewVec(Frames) = CurrentVec;
        else
            NewVec(Frames,:) = CurrentVec;
        end
        Tracks(tr).WormPerimeter.(fieldname) = NewVec;    
    end
    
    % correct out-of-bounds:  add all NaN frames
    Tracks(tr).OutOfBounds(NewOutOfBounds) = true;    
end

return

function Tracks = CorrectHeadTailSegmentationBasedOnNeuronCoordinates(Tracks, plotme)

if ~exist('plotme','var')
    plotme = false;
end
% Correct Head-Tail based on NeuronCoordinates 
for tr = 1:length(Tracks)
    HeadCoordinates             = Tracks(tr).Head;
    TailCoordinates             = Tracks(tr).Tail;
    NeuronCoordinates           = Tracks(tr).Neuron.CoordinatesMatrix;
    HeadToNeuronDistance        = sum((NeuronCoordinates-HeadCoordinates).^2,2);
    TailToNeuronDistance        = sum((NeuronCoordinates-TailCoordinates).^2,2);
    FlipIndices                 = find(HeadToNeuronDistance>TailToNeuronDistance);
    
    if plotme
        figure('name',['Track ',num2str(tr),', Head and neuron coordinates']); 
        plot(HeadCoordinates(:,1),'b:'); hold on;  
        plot(HeadCoordinates(:,2),'r:'); hold on;  
        plot(NeuronCoordinates(:,1),'b'); hold on;  
        plot(NeuronCoordinates(:,2),'r'); hold on;  
        legend('Head 1','Head 2','neuron 1','neuron 2')
        ylabel('Coordinates [pix]'); xlabel('Frames'); 

        figure('name',['Track ',num2str(tr),', Head-Neuron and Tail-Neuron distances']); 
        plot(HeadToNeuronDistance,'k'); hold on;  
        plot(TailToNeuronDistance,'g'); hold on;  
        plot(FlipIndices,TailToNeuronDistance(FlipIndices),'r*'); hold on;  
        ylabel('distance [pix]'); xlabel('Frames'); 
        legend('Head-Neuron','Tail-Neuron','Tail Segmentation is actually Head')
    end
    
    HeadCoordinates_old = HeadCoordinates;
    TailCoordinates_old = TailCoordinates;
    HeadCoordinates(FlipIndices,:) = TailCoordinates_old(FlipIndices,:);
    TailCoordinates(FlipIndices,:) = HeadCoordinates_old(FlipIndices,:);

    Tracks(tr).Head = HeadCoordinates;
    Tracks(tr).Tail = TailCoordinates;
end

return

function [Tracks, GradientAxisLimits, FlowAxisLimits] = CorrectOOBBasedOnNeuronCoordinatesAndMask(Tracks, Mask, plotme)
% Free parameters, based on magnification and arena size
PixelsFromLeftEdge  = 122; % 122 until the middle of the 2nd column. 61 until the begining of the first column;
PixelsFromRightEdge = 90;  % 90 until the begining of the 'end' of the second column from the right;
if ~exist('plotme','var')
    plotme = false;
end

% Calculate limits
Xmin =  find(mean(Mask,1)>0.9,1,'first');
Xmax =  find(mean(Mask,1)>0.9,1,'last');
FlowAxisLimits = [Xmin Xmax];
OOBLimit.LowX  = Xmin + PixelsFromLeftEdge;
OOBLimit.HighX = Xmax - PixelsFromRightEdge;
% Find Ymin, Ymax. This is the gradient axis limits. Note!! it is the 'row' dimension in Mask but in imshow(Mask) it is the Y axis.    
Ymin =  find(mean(Mask,2)>=0.55,1,'first');
Ymax =  find(mean(Mask,2)>=0.55,1,'last');
GradientAxisLimits = [Ymin Ymax];

if plotme
    figure; Vec=mean(Mask,2); plot(1:length(Vec),Vec); hold on; plot([Ymin Ymax],Vec([Ymin Ymax]),'r*'); title('gradient axis limits')
    figure; Vec=mean(Mask,1); plot(1:length(Vec),Vec); hold on; plot([Xmin Xmax],Vec([Xmin Xmax]),'r*'); title('Flow axis limits')
    figure; imshow(Mask,[]); hold on; 
    rectangle('Position',[Xmin,Ymin,Xmax-Xmin,Ymax-Ymin],'edgecolor','r','linewidth',2)
    rectangle('Position',[Xmin,Ymin,OOBLimit.LowX-Xmin,Ymax-Ymin],           'edgecolor','b','linewidth',2)
    rectangle('Position',[OOBLimit.HighX,Ymin,Xmax-OOBLimit.HighX,Ymax-Ymin],'edgecolor','b','linewidth',2)
end

% Correct OOB based on NeuronCoordinates 
for tr = 1:length(Tracks)    
    NearBarrierEdges            = (Tracks(tr).Path(:,2)<OOBLimit.LowX) | (Tracks(tr).Path(:,2)>OOBLimit.HighX);
    NoNeuronCoordinates         = isnan(Tracks(tr).Neuron.CoordinatesMatrix(:,1));
    NoWormCoordinates           = isnan(Tracks(tr).Path(:,1));
    OutOfBounds_new             = NoWormCoordinates | (NoNeuronCoordinates & NearBarrierEdges);
    OutOfBounds                 = OutOfBounds_new';    
    Tracks(tr).OutOfBounds      = OutOfBounds;    
   
    if plotme
        figure('name',['Track ',num2str(tr),', Worm coordinates, red=OOB']); 
        plot(Tracks(tr).Path(:,2), Tracks(tr).Path(:,1),'b.'); hold on; 
        plot(Tracks(tr).Path(OutOfBounds,2), Tracks(tr).Path(OutOfBounds,1),'r.')
        ylabel('Coordinates [pix]'); xlabel('Coordinates [pix]'); 
        
        figure('name',['Track ',num2str(tr),', Head coordinates, red=OOB']); 
        plot(Tracks(tr).Head(:,2), Tracks(tr).Head(:,1),'b.'); hold on; 
        plot(Tracks(tr).Head(OutOfBounds,2), Tracks(tr).Head(OutOfBounds,1),'r.')        
        ylabel('Coordinates [pix]'); xlabel('Coordinates [pix]'); 
    end    
end

return

function Tracks = CalculateRelativeCoordinates(Tracks, GradientAxisLimits, FlowAxisLimits)
% For Stimulus:  -1 corresponds to bottom if the field of view (High values of Dimension 1)   
% For Flow:      -1 corresponds to near the outlet barriers    (High values of Dimension 2)    

for dimension = 1:2
    switch dimension
        case 1
            FieldOut  = 'Stimulus'; 
            Limits    = GradientAxisLimits;
        case 2
            FieldOut  = 'FlowAxis'; 
            Limits    = FlowAxisLimits;
    end            
    for tr = 1:length(Tracks)
        HeadCoordinates     = Tracks(tr).Head(:,dimension);
        NeuronCoordinates   = Tracks(tr).Neuron.CoordinatesMatrix(:,dimension);
        CentroidCoordinates = Tracks(tr).Path(:,dimension);

        HeadCoordinates(HeadCoordinates<Limits(1))=Limits(1);
        HeadCoordinates(HeadCoordinates>Limits(2))=Limits(2);
        HeadCoordinates = -((HeadCoordinates-Limits(1))/(Limits(2)-Limits(1))*2 - 1); % between (-1) and 1, where -1 is high values 

        NeuronCoordinates(NeuronCoordinates<Limits(1))=Limits(1);
        NeuronCoordinates(NeuronCoordinates>Limits(2))=Limits(2);
        NeuronCoordinates = -((NeuronCoordinates-Limits(1))/(Limits(2)-Limits(1))*2 - 1); % between (-1) and 1, where -1 is high values 

        CentroidCoordinates(CentroidCoordinates<Limits(1))=Limits(1);
        CentroidCoordinates(CentroidCoordinates>Limits(2))=Limits(2);
        CentroidCoordinates = -((CentroidCoordinates-Limits(1))/(Limits(2)-Limits(1))*2 - 1); % between (-1) and 1, where -1 is high values 

        figure('name',[FieldOut,' , track ',num2str(tr)]); 
        plot(HeadCoordinates,'r-'); hold on; plot(NeuronCoordinates,'b-'); plot(CentroidCoordinates,'k-');

        Tracks(tr).([FieldOut,'_Head'])     = HeadCoordinates;
        Tracks(tr).([FieldOut,'_Neuron'])   = NeuronCoordinates;
        Tracks(tr).([FieldOut,'_Centroid']) = CentroidCoordinates;
    end
end

return




