function StitchFiles_MultipleSessions_v02(File, ar)

%%  Stitching Tracks from all Fragments (only "safe" stitching is allowed in this stage) 
DetectionMode = File.VariablesInformation.DetectionMode;
if strcmpi(DetectionMode, 'AddAllProperties')
    AddProps = 3;
    disp('The code for stitching tracks with animal pattern is not available in this function. Aborting');
    return
elseif strcmpi(DetectionMode, 'AddAdvancedMorphologyProperties')
    AddProps = 2;
elseif strcmpi(DetectionMode, 'AddBasicMorphologyProperties')
    AddProps = true; 
else
    AddProps = false;
end    
NumFragments  = File.Fragments;  

%% STEP 1: load and concaenate all Tracks data
% Combine necessary Tracks information from all fragments to one structure. Use only the small Tracks structures
disp([datestr(now),' -- LOADING AND STITCHING TRACKS FILE -- Arena ',num2str(ar)]);    
AllFragmentsTracksFileName = File.FileNames.ConcatinatedTracks{ar};

if (File.VariablesInformation.Background_params.UseOnlyFramesInFragement) && (ar==1)      % backgrounds are the same for different arenas
    AllBackgrounds.FragmentFrames  = File.FragmentFrames;
    AllBackgrounds.Matrix          = cell(1,NumFragments);
end

for Fragment = 1:NumFragments
    FileName = [File.FragmentSaveNames{Fragment}(1:end-4),'_Arena',num2str(ar),'.mat'];
    disp([datestr(now),' -- Loading Track File ',int2str(Fragment),'... ']);

    loadsuccess = 0;
    while ~loadsuccess 
        try
            if (File.VariablesInformation.Background_params.UseOnlyFramesInFragement) && (ar==1)  % backgrounds are the same for different arenas
                load(FileName,'Tracks','background');  
                AllBackgrounds.Matrix{Fragment} = background;
            else
                load(FileName,'Tracks');          
            end 
            loadsuccess = 1;
        catch
            disp(['Error loading Track File ',int2str(Fragment),'.  Retrying in 20 seconds...']);
            pause(20);
        end
    end                      

    if Fragment==1
        AllTracks  = Tracks;
    else
        AllTracks  = [AllTracks, Tracks];
    end
end

Tracks = AllTracks;
clear AllTracks;
disp('Saving variables to files');
if File.VariablesInformation.Background_params.UseOnlyFramesInFragement 
    save(AllFragmentsTracksFileName,'Tracks','File','background','AllBackgrounds','-v7.3');      
else
    load(File.FragmentSaveNames{Fragment},'background');  
    save(AllFragmentsTracksFileName,'Tracks','File','background','-v7.3');        
end

%% STEP 2: Linking Tracks based on location in trivial non-conflicting events 
disp([datestr(now),' -- SAFE LINKING TRACKS BASED ON LOCATION. ONLY WHEN LINKAGE IS OBVIOUS.']);
%%%%%%    Initialize linking parameters    %%%%%%
SafeStitchedTracksFileName          = File.FileNames.SafeStitchedTracks{ar};    
MaxSpeedForTrackLinking_mm_sec      = File.VariablesInformation.MaxSpeedForTrackLinking_mm_sec;
MaxTimeForLinkingBasedOnLocation1   = File.VariablesInformation.MaxTimeForLinkingBasedOnLocation1;
MaxTimeForLinkingBasedOnLocation2   = File.VariablesInformation.MaxTimeForLinkingBasedOnLocation2;      
MaxPixelsPerFrame                   = MaxSpeedForTrackLinking_mm_sec * File.PixelSize / File.FrameRate;    % File.PixelSize is actually How many pixels per mm. MaxDistance has units of pixels/frame 
MaxFramesForLinkingBasedOnLocation1 = MaxTimeForLinkingBasedOnLocation1 * File.FrameRate;                  % frames
MaxFramesForLinkingBasedOnLocation2 = MaxTimeForLinkingBasedOnLocation2 * File.FrameRate;                  % frames
TracksStats.BeforeStitching         = CalculateTracksStats (Tracks);

%%%%%%    safe- linking tracks    %%%%%%
GoodLinks1                       = FindGoodLinks  (Tracks, File, MaxFramesForLinkingBasedOnLocation1,  MaxPixelsPerFrame);   % first link 'close' tracks (frame-wise)
if ~isempty(GoodLinks1)
    Tracks                       = LinkGoodTracks (Tracks,GoodLinks1, AddProps);
end
[GoodLinks2, NumOfDetectedWorms] = FindGoodLinks  (Tracks, File, MaxFramesForLinkingBasedOnLocation2,  MaxPixelsPerFrame);   % then allow 'far' tracks
if ~isempty(GoodLinks2)
    Tracks                       = LinkGoodTracks (Tracks,GoodLinks2, AddProps);
end
TracksStats.AfterStitching                    = CalculateTracksStats (Tracks);
TracksStats.AfterStitching.GoodLinks1         = GoodLinks1;
TracksStats.AfterStitching.GoodLinks2         = GoodLinks2;
TracksStats.AfterStitching.NumOfDetectedWorms = NumOfDetectedWorms;
TracksStats.AfterStitching.MaxNumOfWorms      = max(NumOfDetectedWorms);    
File.NumberOfTracks(ar) = length(Tracks);
File.MaxNumOfWorms(ar)  = max(NumOfDetectedWorms);    

%%%%%   Correct midline based on deviation from the median midline lengths of all worms within the arena  
disp([datestr(now),' -- FLAGGING ADDITIONAL MIDLINE CALCULATION ERRORS']);
Tracks = IdentifyMidlineErrors (Tracks,File);

%% Save safe-stitched data file 
disp([datestr(now),' -- SAVING STITCHED TRACK FILE']);
if File.VariablesInformation.Background_params.UseOnlyFramesInFragement 
    save(SafeStitchedTracksFileName,'Tracks','File','background','TracksStats','AllBackgrounds','-v7.3');      
else
    load(File.FragmentSaveNames{Fragment},'background');  
    save(SafeStitchedTracksFileName,'Tracks','File','background','TracksStats','-v7.3');        
end       

load (File.StatusFile,'Stitched_arenas');
Stitched_arenas(ar)=1;
save(File.StatusFile,'Stitched_arenas','-append');

quit;

return

function STATS = CalculateTracksStats (Tracks)    
    STATS.TrackLengths.mean   = nanmean([Tracks.TrackLength]);
    STATS.TrackLengths.median = nanmedian([Tracks.TrackLength]);
    STATS.TrackLengths.std    = nanstd([Tracks.TrackLength],1);
    [N,X]                     = hist([Tracks.TrackLength],100);
    STATS.TrackLengths.N      = N/sum(N);
    STATS.TrackLengths.X      = X;
    STATS.NumberOfTracks      = length(Tracks);
return

function [LinkageMatrix, LinkageStats, ProximalTracksMatrix] = CheckLinkageBasedOnLocation (truncated_tracks,PossibleNextTracks, TracksEdgeParameters, MaxPixelsPerFrame) % LinkageMatrix: column 1 = truncated tracks, column 2 = Next Tracks;

Coordinates.TruncatedTracks    = TracksEdgeParameters.End.XYcoordinates(truncated_tracks,:);
Coordinates.PossibleNextTracks = TracksEdgeParameters.Start.XYcoordinates(PossibleNextTracks,:);
Frames.TruncatedTracks         = TracksEdgeParameters.End.Frame(truncated_tracks);
Frames.PossibleNextTracks      = TracksEdgeParameters.Start.Frame(PossibleNextTracks);

DistanceMatrix         = single(zeros(length(truncated_tracks),length(PossibleNextTracks)));
DistancePerFrameMatrix = single(zeros(length(truncated_tracks),length(PossibleNextTracks)));

for tr_ind1 = 1:length(truncated_tracks)
    CurrentTruncatedTrack = truncated_tracks(tr_ind1);
    for tr_ind2 = 1:length(PossibleNextTracks)
        CurrentNextTrack = PossibleNextTracks(tr_ind2);
        
        DistanceX = Coordinates.TruncatedTracks(tr_ind1,1) - Coordinates.PossibleNextTracks(tr_ind2,1);
        DistanceY = Coordinates.TruncatedTracks(tr_ind1,2) - Coordinates.PossibleNextTracks(tr_ind2,2);
        Distance  = sqrt(DistanceX.^2 + DistanceY.^2);    % Distance in pixels b/w each the end of the first track and the beginning of the next one. 
        
        FrameInterval    = Frames.PossibleNextTracks(tr_ind2) - Frames.TruncatedTracks (tr_ind1);
        DistancePerFrame =  Distance/FrameInterval;
        
        DistanceMatrix(tr_ind1,tr_ind2)         = Distance;
        DistancePerFrameMatrix(tr_ind1,tr_ind2) = DistancePerFrame;
    end
end                  
ProximalTracksMatrix      = DistancePerFrameMatrix<=MaxPixelsPerFrame;

NextTracksCollisions      = sum(ProximalTracksMatrix,1)>1;  % more than one truncated track is proximal to its beginning
TruncatedTracksCollisions = sum(ProximalTracksMatrix,2)>1;  % more then one possible next track is proximal to its ending
LinkageMatrix = [];

for tr_ind = 1:length(truncated_tracks)
    CurrentTruncatedTrack = truncated_tracks(tr_ind);
    if ~TruncatedTracksCollisions(tr_ind) % up to one possible next track is available
        CurrentNextTrack_Ind = find(ProximalTracksMatrix(tr_ind,:)==1);
        if (~isempty(CurrentNextTrack_Ind)) && (~NextTracksCollisions(CurrentNextTrack_Ind))  % if there is one next track that is not proximal to any additional truncated tracks
            CurrentNextTrack = PossibleNextTracks(CurrentNextTrack_Ind);
            LinkageMatrix = [LinkageMatrix; CurrentTruncatedTrack CurrentNextTrack];
        end               
    end
end    
if isempty(LinkageMatrix)
    LinkageStats.truncated_tracks_linked     = [];
    LinkageStats.Possible_Next_tracks_linked = [];      
else
    LinkageStats.truncated_tracks_linked     = LinkageMatrix(:,1);
    LinkageStats.Possible_Next_tracks_linked = LinkageMatrix(:,2);      
end    

return

function [GoodLinks, NumOfDetectedWorms] = FindGoodLinks(Tracks, File, MaxFramesForLinkingBasedOnLocation,  MaxPixelsPerFrame)

NumOfTracks                             = length(Tracks);
NumOfFrames                             = File.NumberOfFrames;
TracksFramesMatrix                      = false(NumOfTracks);
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
NumOfDetectedWorms  = sum(single(TracksFramesMatrix),1);
% Fix truncated tracks without conflicting events
TracksLinked_ToLaterTracks       = false(1,NumOfTracks);
TracksLinked_ToPreviousTracks    = false(1,NumOfTracks);
GoodLinks                        = [];

for frame = 1:NumOfFrames  
    % find tracks that ends at this frame
    truncated_tracks   = find(TracksEdgeParameters.End.Frame==frame);    
    if isempty(truncated_tracks)
        continue
    end
    % find tracks that ends at this frame    
    PossibleNextTracks = find((TracksEdgeParameters.Start.Frame>frame)&(TracksEdgeParameters.Start.Frame<(frame+MaxFramesForLinkingBasedOnLocation)));  % Search for possible tracks to link them
    PossibleNextTracks = PossibleNextTracks(~TracksLinked_ToPreviousTracks(PossibleNextTracks));                                                        % avoid using trakcs that were already linked
    if isempty(PossibleNextTracks)
        continue
    end    
    
    [LinkageMatrix, LinkageStats] = CheckLinkageBasedOnLocation (truncated_tracks,PossibleNextTracks, TracksEdgeParameters, MaxPixelsPerFrame); % LinkageMatrix: column 1 = truncated tracks, column 2 = Next Tracks;
    GoodLinks                     = [GoodLinks; LinkageMatrix];
    
    TracksLinked_ToLaterTracks(LinkageStats.truncated_tracks_linked)        = true;
    TracksLinked_ToPreviousTracks(LinkageStats.Possible_Next_tracks_linked) = true;             
end

return

function Tracks = LinkGoodTracks (Tracks,GoodLinks, AddProps)

NumOfLinks = size(GoodLinks,1);

for link_num =  NumOfLinks:(-1):1   % from the end to the beginning
    tr_ind1 = GoodLinks(link_num,1);
    tr_ind2 = GoodLinks(link_num,2);
    
    Tracks(tr_ind1).Path              = [Tracks(tr_ind1).Path;        Tracks(tr_ind2).Path];
    Tracks(tr_ind1).LastCoordinates   = Tracks(tr_ind2).LastCoordinates;
    Tracks(tr_ind1).Frames            = [Tracks(tr_ind1).Frames,      Tracks(tr_ind2).Frames];
    Tracks(tr_ind1).Size              = [Tracks(tr_ind1).Size,        Tracks(tr_ind2).Size];
    Tracks(tr_ind1).LastSize          = Tracks(tr_ind2).LastSize;
    Tracks(tr_ind1).Eccentricity      = [Tracks(tr_ind1).Eccentricity, Tracks(tr_ind2).Eccentricity];
    Tracks(tr_ind1).MajorAxes         = [Tracks(tr_ind1).MajorAxes   , Tracks(tr_ind2).MajorAxes ];
    Tracks(tr_ind1).MinorAxes         = [Tracks(tr_ind1).MinorAxes   , Tracks(tr_ind2).MinorAxes];
    Tracks(tr_ind1).Orientation       = [Tracks(tr_ind1).Orientation , Tracks(tr_ind2).Orientation ];
    Tracks(tr_ind1).Box               = [Tracks(tr_ind1).Box;          Tracks(tr_ind2).Box];
    TrackFrameNum                     = length(Tracks(tr_ind1).Size);
    Tracks(tr_ind1).TrackLength       = single(TrackFrameNum);
    
    Tracks(tr_ind1).OutOfBounds       = [Tracks(tr_ind1).OutOfBounds   , Tracks(tr_ind2).OutOfBounds];
   
    if AddProps
        
        Tracks(tr_ind1).SkeletonLength              = [Tracks(tr_ind1).SkeletonLength,                  Tracks(tr_ind2).SkeletonLength] ;
        Tracks(tr_ind1).PerimeterLength             = [Tracks(tr_ind1).PerimeterLength,                 Tracks(tr_ind2).PerimeterLength] ;    
        Tracks(tr_ind1).WormPerimeter.Xcoordinate   = [Tracks(tr_ind1).WormPerimeter.Xcoordinate,       Tracks(tr_ind2).WormPerimeter.Xcoordinate];          % vector of indices per worm per frame
        Tracks(tr_ind1).WormPerimeter.Ycoordinate   = [Tracks(tr_ind1).WormPerimeter.Ycoordinate,       Tracks(tr_ind2).WormPerimeter.Ycoordinate];          % vector of indices per worm per frame                    

        if AddProps>=2
            
            Tracks(tr_ind1).HeadTail                    = [Tracks(tr_ind1).HeadTail ;                       Tracks(tr_ind2).HeadTail];              
            Tracks(tr_ind1).Midline.X_coordinates_short = [Tracks(tr_ind1).Midline.X_coordinates_short,     Tracks(tr_ind2).Midline.X_coordinates_short];        
            Tracks(tr_ind1).Midline.Y_coordinates_short = [Tracks(tr_ind1).Midline.Y_coordinates_short,     Tracks(tr_ind2).Midline.Y_coordinates_short];        
            Tracks(tr_ind1).Midline.Angle               = [Tracks(tr_ind1).Midline.Angle,                   Tracks(tr_ind2).Midline.Angle];        
            Tracks(tr_ind1).Midline.AngleFirstPoint     = [Tracks(tr_ind1).Midline.AngleFirstPoint,         Tracks(tr_ind2).Midline.AngleFirstPoint];        
            Tracks(tr_ind1).Midline.AngleLastPoint      = [Tracks(tr_ind1).Midline.AngleLastPoint,          Tracks(tr_ind2).Midline.AngleLastPoint];        
            Tracks(tr_ind1).Midline.NumOfBends_HighRes  = [Tracks(tr_ind1).Midline.NumOfBends_HighRes,      Tracks(tr_ind2).Midline.NumOfBends_HighRes];        
            Tracks(tr_ind1).Midline.NumOfBends_LowRes   = [Tracks(tr_ind1).Midline.NumOfBends_LowRes,       Tracks(tr_ind2).Midline.NumOfBends_LowRes];        
            Tracks(tr_ind1).Midline.FlagTrueIfReliable  = [Tracks(tr_ind1).Midline.FlagTrueIfReliable,      Tracks(tr_ind2).Midline.FlagTrueIfReliable];                           

            if AddProps == 3
                Tracks(tr_ind1).PatternMatrix         = [Tracks(tr_ind1).PatternMatrix;                   Tracks(tr_ind2).PatternMatrix];       
            end

        else   % If midline was not calculated, at least store the basic 'skeleton' information.
            Tracks(tr_ind1).WormSkeleton.Xcoordinate = [Tracks(tr_ind1).WormSkeleton.Xcoordinate,     Tracks(tr_ind2).WormSkeleton.Xcoordinate];        
            Tracks(tr_ind1).WormSkeleton.Ycoordinate = [Tracks(tr_ind1).WormSkeleton.Ycoordinate,     Tracks(tr_ind2).WormSkeleton.Ycoordinate];      
        end
    end
end

TracksToRemove                 = false(1,length(Tracks));
TracksToRemove(GoodLinks(:,2)) = true;
Tracks                         = Tracks(~TracksToRemove);

return

function Tracks = IdentifyMidlineErrors (Tracks, File)
plotme = true;

% Maximum deviation allowed (in [%]) from the median midline lengths that was found in all worms within the arena.                                                                                               
% This is used to flag out all Midlines that were 'too short' due to error in morphology calculations.  
MaximumDeviationFromMedianLength = File.VariablesInformation.MidlineCalculationParams.MaximumDeviationFromMedianLength;  

SkeletonLengths   = double([Tracks.SkeletonLength]);
Midlines          = [Tracks.Midline];
ReliableSkeletons = [Midlines.FlagTrueIfReliable];
SkeletonLengths   = SkeletonLengths(ReliableSkeletons);
MinSkeletonLength = median(SkeletonLengths)*(1-MaximumDeviationFromMedianLength/100);

if plotme
    figure; 
    hist(SkeletonLengths,200); 
    line([MinSkeletonLength MinSkeletonLength],get(gca,'ylim'));
    xlabel('Midline length')
    ylabel('Count')
    title('Histogram before corrections. Midline lengths below the line cutoff are flagged as unreliable')
end

for tr_ind = 1:length(Tracks)
    BadIndices = find(Tracks(tr_ind).SkeletonLength < MinSkeletonLength);
    BadIndices = setdiff(BadIndices,find(~Tracks(tr_ind).Midline.FlagTrueIfReliable));
    Tracks(tr_ind).Midline.FlagTrueIfMidlineWasDetected   = Tracks(tr_ind).Midline.FlagTrueIfReliable;
    Tracks(tr_ind).Midline.FlagTrueIfReliable(BadIndices) = false;    
end

return















