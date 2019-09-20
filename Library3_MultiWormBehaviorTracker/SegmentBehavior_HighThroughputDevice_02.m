function [Tracks, File, TracksStats] = SegmentBehavior_HighThroughputDevice_02(File, ArenaID)
% Note: The Tracks(tr).Path variables has a switched X-Y coordinate relative to all other parameters
% This function uses three external functions:
%      (1) AnalyzeLocomotionProperties_SL_HighThroughputDevice_02 >> locomotion analysis
%      (2) AnalyzeMorphologyProperties_HighThroughputDevice_02    >> morphology analysis
%      (3) SegmentBehavior_SingleTrack_HighThroughputDevice_02    >> behavior analysis
%
%   September 2018, Sagi Levy
%
%% Initialize Settings parameters 
Settings           = File.VariablesInformation.SegmentationSettings;           % Call External Function with segmentation parameters. 
Settings.FrameRate = File.FrameRate; 
Settings.PixelSize = File.PixelSize;

%% Segment Tracks in each file
disp(['Segmenting behavior in arena ',num2str(ArenaID),' of the movie: ',File.MovieName]);
disp([datestr(now),' -- Loading tracks file']);
load(File.FileNames.SafeStitchedTracks{ArenaID},'Tracks','background','TracksStats');        
File.DefineOutOfBoundPolygons = true;

AllowVeryLongTracks = false; % Use true when only few worms are found within in each arena
if isfield(Settings,'AllowVeryLongTracks')
    AllowVeryLongTracks = Settings.AllowVeryLongTracks;
end

%% LOOP OVER ALL TRACKS AND CALCULATE LOCOMOTION AND MORPHOLOGY DYNAMICS
LEN_Tracks = length(Tracks);
if isfield(File,'OutOfBoundPolygons') || ( (isfield(File,'DefineOutOfBoundPolygons'))&&(File.DefineOutOfBoundPolygons))  
    disp('correcting out-of-bound field to include polygons where worms thrash');
    File     = DefineOutOfBoundsForBehaviorSegmentation(File, ArenaID);
    Polygon1 = File.OutOfBoundPolygons{ArenaID,1};
    Polygon2 = File.OutOfBoundPolygons{ArenaID,2};    
    AllPATHs = vertcat(Tracks.Path);

    for tr = 1:LEN_Tracks 
        CurrentPath = Tracks(tr).Path;
        in1         = inpolygon(CurrentPath(:,1),CurrentPath(:,2),Polygon1(:,1),Polygon1(:,2));
        in2         = inpolygon(CurrentPath(:,1),CurrentPath(:,2),Polygon2(:,1),Polygon2(:,2));
        in_outofbound_area = (in1 | in2)'; 
        Tracks(tr).OutOfBounds =   Tracks(tr).OutOfBounds | in_outofbound_area;
    end          
end
disp([datestr(now),' -- Analyzing locomotion data in ',num2str(LEN_Tracks),' tracks']);
h_waitbar  = waitbar(0,['Analyzing locomotion data in ',num2str(LEN_Tracks),' tracks']);

% Analyze locomotion and morphology properties
Sizes                               = [Tracks.Size];
TracksStats.MedianSize              = nanmedian(single(Sizes));
Midlines                            = [Tracks.Midline];
ReliableSkeletons                   = [Midlines.FlagTrueIfReliable];
SkeletonLengths                     = [Tracks.SkeletonLength];
PerimeterLengths                    = [Tracks.PerimeterLength];
TracksStats.MedianSkeleton          = nanmedian(single(SkeletonLengths(ReliableSkeletons)));
TracksStats.MedianPerimeter         = nanmedian(single(PerimeterLengths(ReliableSkeletons)));
TracksStats.MedianPerimeterOverSize = nanmedian(single(PerimeterLengths(ReliableSkeletons))./single(Sizes(ReliableSkeletons)));

for tr = 1:LEN_Tracks         
    waitbar(tr/LEN_Tracks , h_waitbar);  
    Tracks(tr).Path = Tracks(tr).Path(:,[2 1]);                % NOTE !!! FLIPPING TO FIT THE REAL X-Y COORDINATE!!!    
    Head            = squeeze(Tracks(tr).HeadTail(:,1,:));
    Tail            = squeeze(Tracks(tr).HeadTail(:,2,:));    
   
    % Calculate locomotion and morphology dynamics using EXTERNAL FUNCTIONS 
    TrackLocomotionData = AnalyzeLocomotionProperties_SL_HighThroughputDevice_02(Tracks(tr), Settings, Head, Tail, AllowVeryLongTracks);               % Tracks(tr).Path =  Worm Coordinates XY(t), time in units of frames   
    TrackMorphologyData = AnalyzeMorphologyProperties_HighThroughputDevice_02(Tracks(tr), Settings, TracksStats);           
   
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
    Tracks(tr).Head                                  = Head;
    Tracks(tr).Tail                                  = Tail;    
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
        
        HeadAngleRelativeToDirection   = HeadDirection - Tracks(tr).AngularDisplacement; 
        HeadAngleRelativeToDirection   = asin(sin(HeadAngleRelativeToDirection/180*pi))*180/pi;     % Angle of the head relative to movement direction between [-90 to 90] degrees == [left right]
        Tracks(tr).HeadAngleRelativeToDirection = single(HeadAngleRelativeToDirection);             % Angle [degree]
    end
end      

%% Segment behavior properties
disp([datestr(now),'Segmenting behavior in ',num2str(LEN_Tracks),' tracks']);
waitbar(0,h_waitbar,['Segmenting behavior in ',num2str(LEN_Tracks),' tracks']);
warning off
for tr = 1:LEN_Tracks  
    waitbar(tr/LEN_Tracks , h_waitbar);      
   [BehaviorCode, Segments_HighLevelBehavior, Segments_LowLevelBehavior] = SegmentBehavior_SingleTrack_HighThroughputDevice_02(Tracks(tr), tr, Settings);
       
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

%% Save Tracks and behavior matrices to 'SegmentedTracksFileName'
waitbar(0,h_waitbar,'Saving Tracks file');
disp([datestr(now),'Saving Tracks file of arena ',num2str(ArenaID),' of the movie: ',File.MovieName]);

if File.VariablesInformation.Background_params.UseOnlyFramesInFragement 
    save(File.FileNames.SegmentedTracks{ArenaID},'Tracks','File','background','TracksStats','AllBackgrounds','-v7.3');      
else
    save(File.FileNames.SegmentedTracks{ArenaID},'Tracks','File','background','TracksStats','-v7.3');        
end

%% %%%%%     Force Link Tracks together  %%%%%%      
if length(Tracks)>1
    disp([datestr(now),'Forcing Tracks linking']);
    [~ , NumOfDetectedWorms, TrackHistory, FlagTrueIfLinkingIsGood] = FindTrackLinks(Tracks, File);       % f=figure; plot(NumOfDetectedWorms,'b.'); hold on;

    TracksStats.NumOfWorms              = max(NumOfDetectedWorms);
    DetectionMode                       = File.VariablesInformation.DetectionMode;
    Tracks                              = LinkTracks_ByTrackHistory (Tracks,TrackHistory, DetectionMode);
    TracksStats.FlagTrueIfLinkingIsGood = FlagTrueIfLinkingIsGood;
    TracksStats.NumOfDetectedWorms      = NumOfDetectedWorms;
else
    TracksStats.NumOfWorms              = 1;
    TracksStats.FlagTrueIfLinkingIsGood = true;
    TracksStats.NumOfDetectedWorms      = ones(1,File.NumberOfFrames,'single');
end

%%  save Forced-link tracks        
disp([datestr(now),'Saving Linked Tracks file of arena ',num2str(ArenaID),' of the movie: ',File.MovieName]);
if File.VariablesInformation.Background_params.UseOnlyFramesInFragement 
    save(File.FileNames.SegmentedAndForcedStitchedTracks{ArenaID},'Tracks','File','background','TracksStats','AllBackgrounds','-v7.3');      
else
    save(File.FileNames.SegmentedAndForcedStitchedTracks{ArenaID},'Tracks','File','background','TracksStats','-v7.3');        
end

%% Generate and save data structure
%  File containing 'Data' structure with subfields of interesting locomotion/morphology/behavior features of worms. 'Data' STRUCTURE DOESN'T CONTAIN ANY MOVIE DISPLAY FEATURES such as perimeter coordinates. 
%           example for subfields:  Speed, AngularVelocity, Eccentricity, BehaviorCodeLowLevel, BehaviorCodeHighLevel, BehaviorProbLowLevel, BehaviorProbHighLevel, HeadDirection ...   
%  most matrices:      data.(subfield_name) are (NxM) matrices; N=number of worms, M=number of frames. 
%  Behavior matrices:  data.(subfield_name) are (BxM) matrices; B=behavior code number, M=number of frames.  
%  Unavailable data-points will have a NaN value ONLY FOR SINGLE PRECISION MATRICES !!!! and '0' for logical, uint8 and uint16 variables 

[Data, File, Stats] = GenerateDataStructure (Tracks, File, TracksStats.BehaviorCode);
disp([datestr(now),'Saving Data Matrices file of arena ',num2str(ArenaID),' of the movie: ',File.MovieName]);
save(File.FileNames.DataMatrices{ArenaID},'Data','File','background','TracksStats','-v7.3');        

%% Update status file
loadsuccess = 0;
while ~loadsuccess 
    try
        load(File.StatusFile,'File');
        loadsuccess = 1;
    catch
        disp('Error loading status file.  Retrying...');
        pause(5);
    end
end 
File.BehaviorSegmentation.CompletedArenas(ArenaID) = true;
save(File.StatusFile,'File','-append');
disp('*** Status file updated ***');

quit

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
    subplot(2,2,2); 
    plot(ContinuousAngles,'b*-'); hold on;        title ( 'ContinuousAngles');
    subplot(2,2,4); 
    plot(diff(ContinuousAngles),'r*-'); hold on;  title ( 'diff(ContinuousAngles)');
end

return

function File = DefineOutOfBoundsForBehaviorSegmentation(File, ArenaID)

if isfield(File,'OutOfBoundPolygons')
    disp('Out of bounds polygons were already defined');
    return
end

load (File.BackgroundFile,'background');
MovSize = File.FrameSize;
NumArenas = File.NumArenas;
if exist('ArenaID','var')
    ArenasToCorrect = ArenaID;
else
    ArenasToCorrect = 1:NumArenas;
end

figure; imshow(background,[]);

for ar=ArenasToCorrect
    h     = rectangle('Position',File.Arena(ar).TrackBox); set(h,'EdgeColor',[1,0,0]);    
    for polygon_number=1:2
        if polygon_number==1
            side = 'left';
        else
            side = 'right';
        end
        txt = ['Select ' side,' polygon area (arena ',num2str(ar),'). Enter to end.'];
        title(txt); 
        txt = ['Select ' side,' polygon area (arena ',num2str(ar),...
               ')',char(10),'This area is out of bound for behaviour segmentation. Enter to end.'];
        label = text(MovSize(2)/2,MovSize(1)*0.4,txt,'FontSize',18,'HorizontalAlignment','center');
        success = false; 
        while success == 0        
            [X,Y,button] = ginput; 
            X = [X(1:end); X(1)];
            Y = [Y(1:end); Y(1)];
            line_handle(ar) = line(X,Y,'color','b');
            txt = 'Click inside to confirm first polygon or outside to redefine it';
            set(label,'String',txt);
            [x,y]=ginput(1);
            in = inpolygon(x,y,X,Y);
            if in
                success = true;
                set(line_handle(ar),'color','r'); 
                File.OutOfBoundPolygons{ar,polygon_number} = [X,Y];
                delete(label);
            else
                delete(line_handle(ar));
                txt = ['Select ', side,' polygon area (arena ',num2str(ar),...
                       ')',char(10),'This area is out of bound for behaviour segmentation. Enter to end.'];
                set(label,'String',txt);
            end
        end
    end
    delete(h);   
    clear h line_handle       
end

return

function [TrackLinks, NumOfDetectedWorms, TrackHistory, FlagTrueIfLinkingIsGood] = FindTrackLinks(Tracks, File)

MaxFramesForLinkingBasedOnLocation_vec  = [6:40 45:5:90 100:20:700 1000 1500 2000 5000 inf];
    
NumOfTracks                             = length(Tracks);
NumOfFrames                             = File.NumberOfFrames;
TracksFramesMatrix                      = false(NumOfTracks, NumOfFrames);
clear TracksEdgeParameters
TracksEdgeParameters.Start.XYcoordinates  = single(zeros(NumOfTracks,2));
TracksEdgeParameters.End.XYcoordinates    = single(zeros(NumOfTracks,2));
TracksEdgeParameters.Start.Frame          = single(zeros(1,NumOfTracks));
TracksEdgeParameters.End.Frame            = single(zeros(1,NumOfTracks));

for tr_ind = 1:NumOfTracks
    TracksFramesMatrix(tr_ind,Tracks(tr_ind).Frames) = true;
    TracksEdgeParameters.Start.XYcoordinates(tr_ind,:) = Tracks(tr_ind).Path(1,:);
    TracksEdgeParameters.End.XYcoordinates(tr_ind,:)   = Tracks(tr_ind).Path(end,:);
    TracksEdgeParameters.Start.Frame(tr_ind)           = Tracks(tr_ind).Frames(1);
    TracksEdgeParameters.End.Frame(tr_ind)             = Tracks(tr_ind).Frames(end);    
end
NumOfDetectedWorms                 = sum(single(TracksFramesMatrix),1);
TracksLinked_ToLaterTracks         = false(1,NumOfTracks);
TracksLinked_ToPreviousTracks      = false(1,NumOfTracks);
TrackLinks                         = [];


for frame = NumOfFrames:(-1):1   %  1:NumOfFrames  
    % find tracks that ends at this frame
    truncated_tracks   = find(TracksEdgeParameters.End.Frame==frame);    
    if isempty(truncated_tracks)
        continue
    end
    for MaxFramesForLinkingBasedOnLocation = MaxFramesForLinkingBasedOnLocation_vec
        if MaxFramesForLinkingBasedOnLocation==26
            'ma';
        end
        % find tracks that ends at this frame    
        PossibleNextTracks = find((TracksEdgeParameters.Start.Frame>frame)&(TracksEdgeParameters.Start.Frame<(frame+MaxFramesForLinkingBasedOnLocation)));  % Search for possible tracks to link them
        PossibleNextTracks = PossibleNextTracks(~TracksLinked_ToPreviousTracks(PossibleNextTracks));                                                        % avoid using trakcs that were already linked    

        if isempty(PossibleNextTracks)
            continue
        end    
        [LinkageMatrix, LinkageStats] = CheckLinkageBasedOnLocation (truncated_tracks,PossibleNextTracks, TracksEdgeParameters, TracksFramesMatrix(truncated_tracks,:), TracksFramesMatrix(PossibleNextTracks,:)); % LinkageMatrix: column 1 = truncated tracks, column 2 = Next Tracks;
        TrackLinks                     = [TrackLinks; LinkageMatrix];

        TracksLinked_ToLaterTracks(LinkageStats.truncated_tracks_linked)        = true;
        TracksLinked_ToPreviousTracks(LinkageStats.Possible_Next_tracks_linked) = true;      
        
        truncated_tracks = setdiff(truncated_tracks, LinkageStats.truncated_tracks_linked);
        if isempty(truncated_tracks)
            break
        end
    end
end

[TrackHistory, FlagTrueIfLinkingIsGood] = FindTrackHistory(Tracks, TrackLinks, NumOfDetectedWorms);

return

function [LinkageMatrix, LinkageStats] = CheckLinkageBasedOnLocation (truncated_tracks,PossibleNextTracks, TracksEdgeParameters, truncated_TracksFramesMatrix, PossibleNext_TracksFramesMatrix) % LinkageMatrix: column 1 = truncated tracks, column 2 = Next Tracks;
if isempty(truncated_tracks)||isempty(PossibleNextTracks)
    disp('DEBUG LINKAGE MATRIX !!!');
    LinkageStats.truncated_tracks_linked     = [];
    LinkageStats.Possible_Next_tracks_linked = [];      
    LinkageMatrix = [];
    return
end
    
Coordinates.TruncatedTracks    = TracksEdgeParameters.End.XYcoordinates(truncated_tracks,:);
Coordinates.PossibleNextTracks = TracksEdgeParameters.Start.XYcoordinates(PossibleNextTracks,:);

DistanceMatrix         = single(zeros(length(truncated_tracks),length(PossibleNextTracks)));

for tr_ind1 = 1:length(truncated_tracks)
    for tr_ind2 = 1:length(PossibleNextTracks)
        
        OverlappingFrames = sum(truncated_TracksFramesMatrix(tr_ind1,:).* PossibleNext_TracksFramesMatrix(tr_ind2,:));
        if OverlappingFrames
            Distance = inf;  % Don't allow concatination of tracks with overlapping frames
        else
            DistanceX = Coordinates.TruncatedTracks(tr_ind1,1) - Coordinates.PossibleNextTracks(tr_ind2,1);
            DistanceY = Coordinates.TruncatedTracks(tr_ind1,2) - Coordinates.PossibleNextTracks(tr_ind2,2);
            Distance  = sqrt(DistanceX.^2 + DistanceY.^2);    % Distance in pixels b/w each the end of the first track and the beginning of the next one. 
        end
        DistanceMatrix(tr_ind1,tr_ind2) = Distance;        
    end
end                  

if length(truncated_tracks) <= length(PossibleNextTracks)
    LinkageMatrix(:,1) = truncated_tracks';
    IndexOfNextTracks  = zeros(1,length(truncated_tracks));
    for tr_ind = 1:length(truncated_tracks)
        [MinValue, CurrentIndexOfMinimalDistanceNextTracks] = min(DistanceMatrix(tr_ind,:)); 
        if isfinite(MinValue)            
            IndexOfNextTracks(tr_ind)                    = CurrentIndexOfMinimalDistanceNextTracks;
            LinkageMatrix(tr_ind,2)                      = PossibleNextTracks(CurrentIndexOfMinimalDistanceNextTracks);
            DistanceMatrix(:,CurrentIndexOfMinimalDistanceNextTracks) = inf;  
        else
            LinkageMatrix = LinkageMatrix(1:(tr_ind-1),:);
        end
    end  

else    % length(truncated_tracks) > length(PossibleNextTracks)
    LinkageMatrix(:,2) = PossibleNextTracks';
    for tr_ind2 = 1:length(PossibleNextTracks)
        [MinValue, CurrentIndexOfMinimalDistanceTruncatedTracks] = min(DistanceMatrix(:,tr_ind2),[],1); 
        if isfinite(MinValue)            
            LinkageMatrix(tr_ind2,1)                      = truncated_tracks(CurrentIndexOfMinimalDistanceTruncatedTracks);
            DistanceMatrix(CurrentIndexOfMinimalDistanceTruncatedTracks,:) = inf;  
        else
            LinkageMatrix = LinkageMatrix(1:(tr_ind-1),:);
        end
    end  
end

LinkageStats.truncated_tracks_linked     = LinkageMatrix(:,1);
LinkageStats.Possible_Next_tracks_linked = LinkageMatrix(:,2);     

return

function [TrackHistory, FlagTrueIfLinkingIsGood] = FindTrackHistory(Tracks, TrackLinks, NumOfDetectedWorms)

TruncatedTracks = TrackLinks(:,1);
NextTracks      = TrackLinks(:,2);
NumOfWorms      = max(NumOfDetectedWorms);
TrackHistory    = cell(1,NumOfWorms*2);
Overlap = false(1,length(TruncatedTracks));
CurrentNumOfWorms = 0;

for link_num =  1:length(TruncatedTracks)
    truncated_track = TruncatedTracks(link_num);
    next_track      = NextTracks(link_num);
    Overlap(link_num) = ~isempty(intersect(Tracks(truncated_track).Frames, Tracks(next_track).Frames));    
    
    CurrentWorm = 0;
    for worm_num = 1:length(TrackHistory)
        if find(ismember([truncated_track  next_track], TrackHistory{worm_num}))
            CurrentWorm = worm_num;
            break;
        end
    end
    if CurrentWorm==0
        CurrentNumOfWorms = CurrentNumOfWorms+1;
        CurrentWorm       = CurrentNumOfWorms;
    end
    TrackHistory{CurrentWorm} = [TrackHistory{CurrentWorm}, truncated_track, next_track];               
end

for CurrentWorm =  1:length(TrackHistory)
    if isempty(TrackHistory{CurrentWorm})
        break
    else
        TrackHistory{CurrentWorm} = sort(unique( TrackHistory{CurrentWorm} ));
    end
end
TrackHistory = TrackHistory(1:CurrentWorm-1);

% Find overlap
Overlap = zeros(length(TrackHistory),length(TrackHistory));
DeletedWorms = [];
for CurrentWorm1 =  1:length(TrackHistory)
    for CurrentWorm2 =  (CurrentWorm1+1):length(TrackHistory)        
        Overlap(CurrentWorm1,CurrentWorm2) = ~isempty(intersect(TrackHistory{CurrentWorm1},TrackHistory{CurrentWorm2}));
    end
    for CurrentWorm2 = find(Overlap(CurrentWorm1,:))
        TrackHistory{CurrentWorm1} = unique([TrackHistory{CurrentWorm1} TrackHistory{CurrentWorm2}]);
        TrackHistory{CurrentWorm2} = [];
        DeletedWorms  = [DeletedWorms CurrentWorm2];
    end
end
if isempty(DeletedWorms)
    FlagTrueIfLinkingIsGood = true;
else
    FlagTrueIfLinkingIsGood = false;
    disp('WARNING- Some tracks may have been deleted due to bad linking of tracks');
end
    
TrackHistory = TrackHistory(setdiff(1:length(TrackHistory),DeletedWorms));

return

function Tracks = LinkTracks_ByTrackHistory (Tracks,TrackHistory, DetectionMode)
if isempty(TrackHistory)
    return
end

if strcmpi(DetectionMode, 'AddAllProperties')
    AddProps = 3;
    disp('Stitching file for pattern recognition needs to be debugged, or user older version. Aborting...');
    return
elseif strcmpi(DetectionMode, 'AddAdvancedMorphologyProperties')
    AddProps = 2;
elseif strcmpi(DetectionMode, 'AddBasicMorphologyProperties')
    AddProps = true; 
else
    AddProps = false;
end

RowVectorBasicFields            = {'Frames';'Size';'Eccentricity';'MajorAxes';'MinorAxes';'Orientation';'OutOfBounds'};
RowVectorMidlineFields          = {'X_coordinates_short';'Y_coordinates_short';'Angle';'AngleFirstPoint';'AngleLastPoint';'NumOfBends_HighRes';'NumOfBends_LowRes';'FlagTrueIfReliable'};
ColVectorAndMatricesBasicFields = {'Path';'Box'};

RowVectorFields = {'OutOfBounds';'Speed';'AngularDisplacement';'Velocity_X';'Velocity_Y';'AngularDisplacement_ContinuousAngle';'Acceleration';'Jerk';'AngularVelocity';'Path_smoothed_X_coordinate';'Path_smoothed_Y_coordinate';...
                   'Curvature';'PerimeterOverSize';'NormalizedSize';'NormalizedSkeletonLength';'NormalizedPerimeterLength';'NormalizedPerimeterOverSize';'HeadDirection';'HeadDirectionChange';'HeadAngleRelativeToDirection';...
                   'BehaviorCode_HighLevel';'BehaviorCode_LowLevel'};

TracksToKeep = zeros(1,length(TrackHistory));
for Worm_num = 1:length(TrackHistory)
%     Worm_num
    CurrentTrackHistory = TrackHistory{Worm_num};
    tr_ind1             = CurrentTrackHistory(1);
    TracksToKeep(Worm_num) = tr_ind1;
    
    for link_num =  2:length(CurrentTrackHistory)     
        link_num;      
        tr_ind2 = CurrentTrackHistory(link_num);    

        for f_ind=1:length(RowVectorBasicFields)
            fieldname = RowVectorBasicFields{f_ind};
            Tracks(tr_ind1).(fieldname) = [Tracks(tr_ind1).(fieldname), Tracks(tr_ind2).(fieldname)];
            Tracks(tr_ind2).(fieldname) = [];
        end
        for f_ind=1:length(ColVectorAndMatricesBasicFields)
            fieldname = ColVectorAndMatricesBasicFields{f_ind};
            Tracks(tr_ind1).(fieldname) = [Tracks(tr_ind1).(fieldname); Tracks(tr_ind2).(fieldname)];
            Tracks(tr_ind2).(fieldname) = [];
        end   
        Tracks(tr_ind1).LastCoordinates   = Tracks(tr_ind2).LastCoordinates;
        Tracks(tr_ind1).LastSize          = Tracks(tr_ind2).LastSize;
        TrackFrameNum                     = length(Tracks(tr_ind1).Size);
        Tracks(tr_ind1).TrackLength       = single(TrackFrameNum);

         if AddProps

            Tracks(tr_ind1).SkeletonLength              = [Tracks(tr_ind1).SkeletonLength,                  Tracks(tr_ind2).SkeletonLength] ;
            Tracks(tr_ind1).PerimeterLength             = [Tracks(tr_ind1).PerimeterLength,                 Tracks(tr_ind2).PerimeterLength] ;    
            Tracks(tr_ind1).WormPerimeter.Xcoordinate   = [Tracks(tr_ind1).WormPerimeter.Xcoordinate,       Tracks(tr_ind2).WormPerimeter.Xcoordinate];          % vector of indices per worm per frame
            Tracks(tr_ind1).WormPerimeter.Ycoordinate   = [Tracks(tr_ind1).WormPerimeter.Ycoordinate,       Tracks(tr_ind2).WormPerimeter.Ycoordinate];          % vector of indices per worm per frame                    
            Tracks(tr_ind2).SkeletonLength              = [] ;
            Tracks(tr_ind2).PerimeterLength             = [] ;    
            Tracks(tr_ind2).WormPerimeter.Xcoordinate   = [];         
            Tracks(tr_ind2).WormPerimeter.Ycoordinate   = [];                      

            if AddProps>=2
                for f_ind=1:length(RowVectorMidlineFields)
                    fieldname = RowVectorMidlineFields{f_ind};
                    Tracks(tr_ind1).Midline.(fieldname) = [Tracks(tr_ind1).Midline.(fieldname), Tracks(tr_ind2).Midline.(fieldname)];
                    Tracks(tr_ind2).Midline.(fieldname) = [];
                end
                if AddProps == 3
                    Tracks(tr_ind1).PatternMatrix         = [Tracks(tr_ind1).PatternMatrix;                   Tracks(tr_ind2).PatternMatrix];       
                    Tracks(tr_ind2).PatternMatrix         = [];       
                end

            else   % If midline was not calculated, at least store the basic 'skeleton' information.
                Tracks(tr_ind1).WormSkeleton.Xcoordinate = [Tracks(tr_ind1).WormSkeleton.Xcoordinate,     Tracks(tr_ind2).WormSkeleton.Xcoordinate];        
                Tracks(tr_ind1).WormSkeleton.Ycoordinate = [Tracks(tr_ind1).WormSkeleton.Ycoordinate,     Tracks(tr_ind2).WormSkeleton.Ycoordinate];      
                Tracks(tr_ind2).WormSkeleton.Xcoordinate = [];        
                Tracks(tr_ind2).WormSkeleton.Ycoordinate = [];      
            end
        end

        %% Head tail and behavior segmentation
         Tracks(tr_ind1).Head = [Tracks(tr_ind1).Head ;                       Tracks(tr_ind2).Head];      
         Tracks(tr_ind1).Tail = [Tracks(tr_ind1).Tail ;                       Tracks(tr_ind2).Tail];      
         Tracks(tr_ind2).Head = [];      
         Tracks(tr_ind2).Tail = [];      
         for f_ind=1:length(RowVectorFields)
             fieldname = RowVectorFields{f_ind};
             Tracks(tr_ind1).(fieldname) = [Tracks(tr_ind1).(fieldname), Tracks(tr_ind2).(fieldname)];
             Tracks(tr_ind2).(fieldname) = [];
         end                    
    end
end

Tracks = Tracks(TracksToKeep);

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
    Tracks(tr_ind).OutOfBounds_WithoutCollisions = Tracks(tr_ind).OutOfBounds;  
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
                          'HeadCoordinates_X','HeadCoordinates_Y','TailCoordinates_X','TailCoordinates_Y'};         
Uint16Fields           = {'Size','SkeletonLength','PerimeterLength'};      
Uint8Fields            = {'NumOfBends_HighRes','NumOfBends_LowRes','BehaviorCode_HighLevel','BehaviorCode_LowLevel'};
LogicalFields          = {'OutOfBounds_WithoutCollisions','FlagTrueIfMidlineIsReliable'};  

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

%% Compute behavior probability
Interpolation_Factor = 1e-4;
BehaviorCodeNumbers  = BehaviorCode.LowLevel_CodeNumbers;  
BehaviorCodeMAT      = Data.BehaviorCode_LowLevel;
[BehaviorProbability, BehaviorProbability_Smoothed, DetectionAndSegmentationProbabilityVector] = CalculateBehaviorProbability (BehaviorCodeMAT, BehaviorCodeNumbers, Interpolation_Factor);
Data.BehaviorProbability_LowLevel                        = BehaviorProbability;
Data.BehaviorProbability_Smoothed_LowLevel               = BehaviorProbability_Smoothed;
STATS.DetectionAndSegmentationProbabilityVector_LowLevel = DetectionAndSegmentationProbabilityVector;

Interpolation_Factor = 1e-4;
BehaviorCodeNumbers  = BehaviorCode.HighLevel_CodeNumbers;  
BehaviorCodeMAT      = Data.BehaviorCode_HighLevel;
[BehaviorProbability, BehaviorProbability_Smoothed, DetectionAndSegmentationProbabilityVector] = CalculateBehaviorProbability (BehaviorCodeMAT, BehaviorCodeNumbers, Interpolation_Factor);
Data.BehaviorProbability_HighLevel                        = BehaviorProbability;
Data.BehaviorProbability_Smoothed_HighLevel               = BehaviorProbability_Smoothed;
STATS.DetectionAndSegmentationProbabilityVector_HighLevel = DetectionAndSegmentationProbabilityVector;

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






