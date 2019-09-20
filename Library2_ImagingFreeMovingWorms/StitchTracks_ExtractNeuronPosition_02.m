function File = StitchTracks_ExtractNeuronPosition_02(File)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Load and concatenate all Tracks data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Combine information from all fragments to one structure 
disp([datestr(now),' -- LOADING AND STITCHING TRACKS FILE']);    
NumFragments  = File.Fragments;  

if (File.VariablesInformation.Background_params.UseOnlyFramesInFragement)       % backgrounds are the same for different arenas
    AllBackgrounds.FragmentFrames  = File.FragmentFrames;
    AllBackgrounds.Matrix          = cell(1,NumFragments);
end
NumberOfTracksPerFragment = zeros(1,NumFragments);
for Fragment = 1:NumFragments
    FileName = [File.FragmentSaveNames{Fragment}(1:end-4),'.mat'];
    disp([datestr(now),' -- Loading Track File ',int2str(Fragment),'... ']);
    loadsuccess = 0;
    while ~loadsuccess 
        try
            if (File.VariablesInformation.Background_params.UseOnlyFramesInFragement)   % backgrounds are the same for different arenas
                load(FileName,'Tracks','background');  
                AllBackgrounds.Matrix{Fragment} = background;
            else
                load(FileName,'Tracks');          
            end 
            NumberOfTracksPerFragment(Fragment) = length(Tracks);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Link Tracks automatically based on location when stitching is safe, then allow manual linkage for non-safe track stitching 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Initialization
arena=1;
MaxWormPerimeterInMicroMeters = 2e3;
MaxMeanPerimeterInPixels      = MaxWormPerimeterInMicroMeters/ (1e3/File.PixelSize);  % in our case ~320
MaxNumOfWorms                 = File.VariablesInformation.MaxNumOfWorms;
NumberOfWormsInFOVIsConstant  = false;
if isfield(File.VariablesInformation,'NumberOfWormsInFOVIsConstant')
    if File.VariablesInformation.NumberOfWormsInFOVIsConstant
        NumberOfWormsInFOVIsConstant = true;
    end
end
%%%%%%    Initialize linking parameters    %%%%%%
MaxFramesForLinkingBasedOnLocationVec = File.VariablesInformation.MaxTimeForLinkingBasedOnLocationVec * File.FrameRate;
MaxPixelsPerFrameVec_Link1            = File.VariablesInformation.MaxSpeedForTrackLinking_mm_sec_Link1 * File.PixelSize / File.FrameRate;    % File.PixelSize is actually How many pixels per mm. MaxDistance has units of pixels/frame 
MaxPixelsPerFrameVec_Link2            = File.VariablesInformation.MaxSpeedForTrackLinking_mm_sec_Link2 * File.PixelSize / File.FrameRate;    % File.PixelSize is actually How many pixels per mm. MaxDistance has units of pixels/frame 
MaxPixelsForTrackLinkingVec           = File.VariablesInformation.MaxAbsoluteDistanceForTrackLinking_mm * File.PixelSize;    % File.PixelSize is actually How many pixels per mm. MaxDistance has units of pixels/frame 

%% Link Tracks automatically based on location when stitching is safe
disp([datestr(now),' -- SAFE LINKING TRACKS BASED ON LOCATION. ONLY WHEN LINKAGE IS SAFE.']);
TracksStats.BeforeStitching         = CalculateTracksStats (Tracks);

% Exclude tracks that have too large perimeter, indicating of wrong segmentaion
Non_Deleted_Tracks_BeforeStitching  = true(1,length(Tracks));
AveragePerimeter                    = zeros(1,length(Tracks));
for tr=1:length(Tracks)
    AveragePerimeter(tr) = mean(double(Tracks(tr).PerimeterLength));
end
Non_Deleted_Tracks_BeforeStitching(AveragePerimeter > MaxMeanPerimeterInPixels) = false;
% assuming long tracks (>1 second) are more probable to be correctly segmenting the worm perimeter  
MaxAveragePerimeterForShortTracks = prctile(AveragePerimeter([AllTracks.TrackLength]>File.FrameRate),99); 
Non_Deleted_Tracks_BeforeStitching( (AveragePerimeter > MaxAveragePerimeterForShortTracks)& ([Tracks.TrackLength]<File.FrameRate)) = false;
if MaxNumOfWorms==1 % use all tracks if only 1 worm exists in the arena
    Non_Deleted_Tracks_BeforeStitching  = true(1,length(Tracks));
end    
Tracks = Tracks(Non_Deleted_Tracks_BeforeStitching);

%%%%%%    safe- linking tracks    %%%%%%
GoodLinks_ALL_ind = 0;
if NumberOfWormsInFOVIsConstant 
    GoodLinks = FindGoodLinks_BasedOnConstantNumberOfWormsInFOV  (Tracks, File);   % first link 'close' tracks (frame-wise)
    if ~isempty(GoodLinks)
        Tracks = LinkGoodTracks (Tracks,GoodLinks); 
        GoodLinks_ALL_ind                = GoodLinks_ALL_ind + 1;
        GoodLinks_ALL{GoodLinks_ALL_ind} = GoodLinks;
        FramesAreUnique = CheckUniqueFrames(Tracks);
    end          
end
% add PossibleCollision field 
for tr = 1:length(Tracks)
    Tracks(tr).PossibleCollision = false(1,File.NumberOfFrames);
end

% continue safe-linking: from lower to higher 'risk' linking
for MaxFrames_ind = 1:length(MaxFramesForLinkingBasedOnLocationVec)    
    CurrentMaxFramesForLinkingBasedOnLocation = MaxFramesForLinkingBasedOnLocationVec(MaxFrames_ind);
    CurrentMaxPixelsForTrackLinking           = MaxPixelsForTrackLinkingVec(MaxFrames_ind);
    for speed_ind=1:length(MaxPixelsPerFrameVec_Link1)
        switch speed_ind
            case 1
                CurrentMaxPixelsPerFrame = MaxPixelsPerFrameVec_Link1(MaxFrames_ind);                
            case 2
                CurrentMaxPixelsPerFrame = MaxPixelsPerFrameVec_Link2(MaxFrames_ind);
        end

        [GoodLinks, NumOfDetectedWorms] = FindGoodLinks  (Tracks, File, CurrentMaxFramesForLinkingBasedOnLocation,  CurrentMaxPixelsPerFrame, CurrentMaxPixelsForTrackLinking);   % first link 'close' tracks (frame-wise)
        if ~isempty(GoodLinks)
            PossibleCollision                = true;
            Tracks                           = LinkGoodTracks (Tracks,GoodLinks, PossibleCollision); 
            GoodLinks_ALL_ind                = GoodLinks_ALL_ind + 1;
            GoodLinks_ALL{GoodLinks_ALL_ind} = GoodLinks;
            FramesAreUnique = CheckUniqueFrames(Tracks);
        end        
        if NumberOfWormsInFOVIsConstant
            GoodLinks = FindGoodLinks_BasedOnConstantNumberOfWormsInFOV  (Tracks, File);   % first link 'close' tracks (frame-wise)
            if ~isempty(GoodLinks)
                PossibleCollision                = false;
                Tracks                           = LinkGoodTracks (Tracks,GoodLinks, PossibleCollision); 
                GoodLinks_ALL_ind                = GoodLinks_ALL_ind + 1;
                GoodLinks_ALL{GoodLinks_ALL_ind} = GoodLinks;
                FramesAreUnique = CheckUniqueFrames(Tracks);
            end          
        end
    end
end
TracksStats.AfterSafeStitching                    = CalculateTracksStats (Tracks);
TracksStats.AfterSafeStitching.GoodLinks_ALL      = GoodLinks_ALL;
TracksStats.AfterSafeStitching.NumOfDetectedWorms = NumOfDetectedWorms;
TracksStats.AfterSafeStitching.MaxNumOfWorms      = max(NumOfDetectedWorms);    

File.NumberOfTracks = length(Tracks);
File.MaxNumOfWorms  = max(NumOfDetectedWorms);    

%% GUI allowing manual linkage for non-safe track stitching 
Non_Deleted_Tracks_BeforeManualInspection = [Tracks.TrackLength]>=File.FrameRate/2;  % delete tracks that are shorter than 0.5 second
% Do not analyze Tracks that are most of the time out of bounds 
NumberOOBFrames  = zeros(1,length(Tracks));
for tr=1:length(Tracks)
    NumberOOBFrames(tr)  = length(find(Tracks(tr).OutOfBounds));
end
FractionOutOfBounds = NumberOOBFrames./[Tracks.TrackLength];
NumberInBoundFrames = [Tracks.TrackLength] - NumberOOBFrames;
MostlyOutOfBoundsFrames = find((FractionOutOfBounds>0.5)&(NumberInBoundFrames<3*File.FrameRate));
Non_Deleted_Tracks_BeforeManualInspection(MostlyOutOfBoundsFrames) = false;

Tracks                               = Tracks(Non_Deleted_Tracks_BeforeManualInspection);
[EstimatedLinks, NumOfDetectedWorms] = FindLinks_NonExactEstimationByDistances(Tracks, File);

if MaxNumOfWorms~=1 % use all tracks if only 1 worm exists in the arena
    GoodLinks  =  ManualCorrectionOfTrackLinkingEstimation(File,Tracks, EstimatedLinks);  % manual corrections using GUI
    Tracks     =  LinkGoodTracks (Tracks,GoodLinks);     
    GoodLinks_ALL_ind                             = GoodLinks_ALL_ind + 1;
    GoodLinks_ALL{GoodLinks_ALL_ind}              = GoodLinks;
end

TracksStats.AfterStitching                    = CalculateTracksStats (Tracks);
TracksStats.AfterStitching.GoodLinks_ALL      = GoodLinks_ALL;
TracksStats.AfterStitching.NumOfDetectedWorms = NumOfDetectedWorms;
TracksStats.AfterStitching.MaxNumOfWorms      = max(NumOfDetectedWorms);    
TracksStats.AfterStitching.Non_Deleted_Tracks_BeforeManualInspection = Non_Deleted_Tracks_BeforeManualInspection;    

TracksCorrectionParameters = [];
% In case an error accured during manual inspection and needs to be corrected, you can still correct it using the function below in debug mode  
% disp('Entering debugging mode to allow visualization of linking. Run ''dbcont'' to continue');
% keyboard;
% [Tracks, TracksCorrectionParameters] = CorrectLinkingBetweenTwoTracks_inline (Tracks, tr1, tr2, frame); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Extract midline, Head-Tail, Pattern matrix, and neuron coordinates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% remove fields with too much data. Store 'WormArea' seperately for further usage    
WormArea = [Tracks.WormArea];
Tracks   = rmfield(Tracks,{'WormArea','Frame'});

%% Correct midline based on deviation from the median midline lengths. Use all midlines of the same worm in different frames  
disp([datestr(now),' -- FLAGGING ADDITIONAL MIDLINE CALUCLATION ERRORS']);
Tracks = IdentifyMidlineErrors_IndividualWorms (Tracks, File);

%% Update pattern matrix with Tracks Linking. 
% PatternMatrix{tr}(f,W,L) is the pixel value of worm 'tr' at frame 'f' in position (W,L) with respect to the worm Width/Length coordinate system   
VariableName  = 'PatternMatrix';
PatternMatrix = LinkPatternMatrix (File, VariableName, GoodLinks_ALL, Non_Deleted_Tracks_BeforeStitching, Non_Deleted_Tracks_BeforeManualInspection, TracksCorrectionParameters, NumberOfTracksPerFragment);

%% Estimate Head/Tail position based on pattern matrix: automatic estimation followed by Manual correction (of flagged possible errors) GUI
Tracks  = EstimatePositions_Head_Tail_Neuron_inline(PatternMatrix,Tracks);  % Automated estimation
Tracks  = ManualCorrectionOfHeadTailPositions_inline(File,Tracks);          % manual corrections using GUI

%% Correct PatternMatrix direction to fit the correct head/tail segmentation, and then save it
PatternMatrix                   = CorrectPatternMatrix_UsingHeadTailCoordinates    (Tracks, PatternMatrix);
StitchedPatternMatricesFileName = [File.TrackFile(1:end-4),'_StitchedPatternMatrices_Arena',num2str(arena),'.mat'];
save(StitchedPatternMatricesFileName,'PatternMatrix','-v7.3');

%% Find estimated neuron position in worm coordinates and real coordinates and add useful midline properties
plotme = false;
Tracks  = Estimate_Neuron_Position_inline(PatternMatrix,Tracks, plotme);
Tracks  = CorrectMidlineForHeadDirection_AndCalculateProperties(Tracks, File.VariablesInformation.MidlineCalculationParams);  
Tracks  = ExtractEstimatedNeuronPositionInRealCoordinates(Tracks,File, GoodLinks_ALL, Non_Deleted_Tracks_BeforeStitching, Non_Deleted_Tracks_BeforeManualInspection, TracksCorrectionParameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Save Tracks and WormArea to file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
disp([datestr(now),' -- SAVING STITCHED TRACK FILE']);
StitchedTracksFileName = [File.TrackFile(1:end-4),'_StitchedAndManuallyCorrected_Arena',num2str(arena),'.mat'];
if File.VariablesInformation.Background_params.UseOnlyFramesInFragement 
    save(StitchedTracksFileName,'Tracks','File','background','TracksStats','AllBackgrounds','WormArea','-v7.3');      
else
    load(File.FragmentSaveNames{1},'background');  
    save(StitchedTracksFileName,'Tracks','File','background','TracksStats','WormArea','-v7.3');        
end  
if File.VariablesInformation.Background_params.CalcMaskAndVignetting_fromDye
    load(File.BackgroundFile,'background','Mask','VignettingPattern');  
    save(StitchedTracksFileName,'background','Mask','VignettingPattern','-append');              
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Extract Neural activity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
File = ExtractNeuralActivityFromTracksAndOriginalMovie_inline(File, Tracks, WormArea);

return

%% Stitchig functions
function STATS = CalculateTracksStats (Tracks)    
    STATS.TrackLengths.mean   = nanmean([Tracks.TrackLength]);
    STATS.TrackLengths.median = nanmedian([Tracks.TrackLength]);
    STATS.TrackLengths.std    = nanstd([Tracks.TrackLength],1);
    [N,X]                     = hist([Tracks.TrackLength],100);
    STATS.TrackLengths.N      = N/sum(N);
    STATS.TrackLengths.X      = X;
    STATS.NumberOfTracks      = length(Tracks);
return

function [LinkageMatrix, LinkageStats, ProximalTracksMatrix] = CheckLinkageBasedOnLocation (truncated_tracks,PossibleNextTracks, TracksEdgeParameters, MaxPixelsPerFrame, MaxPixels) % LinkageMatrix: column 1 = truncated tracks, column 2 = Next Tracks;

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
        AbsoluteDistanceBetweenTracksIsSmallEnough =  DistanceMatrix(tr_ind,CurrentNextTrack_Ind);
        if (~isempty(CurrentNextTrack_Ind)) && (~NextTracksCollisions(CurrentNextTrack_Ind))&&(AbsoluteDistanceBetweenTracksIsSmallEnough)
            % if there is one next track that is not proximal to any additional truncated tracks, and is also close enough (absolute distance, not per frame)    
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

function [GoodLinks, NumOfDetectedWorms] = FindGoodLinks(Tracks, File, MaxFramesForLinkingBasedOnLocation,  MaxPixelsPerFrame, MaxPixels)

NumOfTracks                             = length(Tracks);
NumOfFrames                             = File.NumberOfFrames;
TracksFramesMatrix                      = false(NumOfTracks);
% TracksEdgeParameters.Start.Size         = single(zeros(1,NumOfTracks);
% TracksEdgeParameters.Start.Eccentricity = single(zeros(1,NumOfTracks);
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
% figure; imagesc(~TracksFramesMatrix); colormap gray
NumOfDetectedWorms  = sum(single(TracksFramesMatrix),1);
% MaxNumOfWorms       = max(NumOfDetectedWorms);
% MinNumOfWorms       = min(NumOfDetectedWorms);
% N                   = histc(NumOfDetectedWorms,((MinNumOfWorms-0.5):(MaxNumOfWorms+0.5)));
% figure; bar(MinNumOfWorms:MaxNumOfWorms,N(1:end-1)) 

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
    
    [LinkageMatrix, LinkageStats] = CheckLinkageBasedOnLocation (truncated_tracks,PossibleNextTracks, TracksEdgeParameters, MaxPixelsPerFrame, MaxPixels); % LinkageMatrix: column 1 = truncated tracks, column 2 = Next Tracks;
    GoodLinks                     = [GoodLinks; LinkageMatrix];
    
    TracksLinked_ToLaterTracks(LinkageStats.truncated_tracks_linked)        = true;
    TracksLinked_ToPreviousTracks(LinkageStats.Possible_Next_tracks_linked) = true;             
end

% figure; [N,X]=hist([Tracks.TrackLength],100); bar(X, N/sum(N)), title('before')

return

function [EstimatedLinks, NumOfDetectedWorms] = FindLinks_NonExactEstimationByDistances(Tracks, File)

NumOfTracks                             = length(Tracks);
NumOfFrames                             = File.NumberOfFrames;
TracksFramesMatrix                      = false(NumOfTracks);
% TracksEdgeParameters.Start.Size         = single(zeros(1,NumOfTracks);
% TracksEdgeParameters.Start.Eccentricity = single(zeros(1,NumOfTracks);
clear TracksEdgeParameters
TracksEdgeParameters.Start.XYcoordinates  = single(zeros(NumOfTracks,2));
TracksEdgeParameters.End.XYcoordinates    = single(zeros(NumOfTracks,2));
TracksEdgeParameters.Start.Frame          = single(zeros(1,NumOfTracks));
TracksEdgeParameters.End.Frame            = single(zeros(1,NumOfTracks));
TracksEdgeParameters.overlapping_tracks   = cell(1,NumOfTracks);

for tr_ind = 1:NumOfTracks
    TracksFramesMatrix(tr_ind,Tracks(tr_ind).Frames) = true;
    TracksEdgeParameters.Start.XYcoordinates(tr_ind,:) = Tracks(tr_ind).Path(1,:);
    TracksEdgeParameters.End.XYcoordinates(tr_ind,:)   = Tracks(tr_ind).Path(end,:);
    TracksEdgeParameters.Start.Frame(tr_ind)           = Tracks(tr_ind).Frames(1);
    TracksEdgeParameters.End.Frame(tr_ind)             = Tracks(tr_ind).Frames(end);   
    TracksEdgeParameters.overlapping_tracks{tr_ind}    = find(any(TracksFramesMatrix & repmat(TracksFramesMatrix(tr_ind,:),NumOfTracks,1),2));
end
% figure; imagesc(~TracksFramesMatrix); colormap gray
NumOfDetectedWorms  = sum(single(TracksFramesMatrix),1);
% MaxNumOfWorms       = max(NumOfDetectedWorms);
% MinNumOfWorms       = min(NumOfDetectedWorms);
% N                   = histc(NumOfDetectedWorms,((MinNumOfWorms-0.5):(MaxNumOfWorms+0.5)));
% figure; bar(MinNumOfWorms:MaxNumOfWorms,N(1:end-1)) 

% Fix truncated tracks without conflicting events
TracksLinked_ToLaterTracks       = false(1,NumOfTracks);
TracksLinked_ToPreviousTracks    = false(1,NumOfTracks);
EstimatedLinks                   = [];

for frame = NumOfFrames:-1:1  
    % find tracks that ends at this frame
    tracks_initiated   = find(TracksEdgeParameters.Start.Frame==frame);    
    if isempty(tracks_initiated)
        continue
    end
    
    for tr_ind = tracks_initiated
        % find tracks that ends before this frame    
        PossiblePreviousTracks = find(TracksEdgeParameters.End.Frame<frame);                        % Search for possible tracks to link them
        PossiblePreviousTracks = setdiff(PossiblePreviousTracks, find(TracksLinked_ToLaterTracks)); % avoid using trakcs that were already linked
        PossiblePreviousTracks = setdiff(PossiblePreviousTracks, TracksEdgeParameters.overlapping_tracks{tr_ind}); % avoid using trakcs with overlapping frames
        if isempty(PossiblePreviousTracks)
            continue
        end   
        % Most recent tracks [closest in time]
        EndTimes                 = TracksEdgeParameters.End.Frame(PossiblePreviousTracks);
        MostRecentPreviousTracks = PossiblePreviousTracks(EndTimes == max(EndTimes));
        NumOfRecentTracks        = length(MostRecentPreviousTracks);
        if NumOfRecentTracks==1
            PreviousTrack = MostRecentPreviousTracks;
        else   % if several options exists: choose the one with the closest distance
            ReferenceXYcoordinates = TracksEdgeParameters.Start.XYcoordinates(tr_ind,:);
            DistanceToReference    = zeros(1,NumOfRecentTracks);
            for ind = 1:NumOfRecentTracks
                DistanceToReference(ind) = sum((TracksEdgeParameters.End.XYcoordinates(MostRecentPreviousTracks(ind),:) - ReferenceXYcoordinates).^2);            
            end

            [~,ind_of_closest_track] = min(DistanceToReference);
            PreviousTrack = MostRecentPreviousTracks(ind_of_closest_track);
        end
   
        EstimatedLinks                             = [EstimatedLinks; [PreviousTrack  tr_ind] ];    
        TracksLinked_ToLaterTracks(PreviousTrack)  = true;
        TracksLinked_ToPreviousTracks(tr_ind)      = true;   
    end
end

% figure; [N,X]=hist([Tracks.TrackLength],100); bar(X, N/sum(N)), title('before')

return

function [GoodLinks, NumOfDetectedWorms] = FindGoodLinks_BasedOnConstantNumberOfWormsInFOV(Tracks, File)

NumOfTracks                             = length(Tracks);
NumOfFrames                             = File.NumberOfFrames;
TracksFramesMatrix                      = false(NumOfTracks);
% TracksEdgeParameters.Start.Size         = single(zeros(1,NumOfTracks);
% TracksEdgeParameters.Start.Eccentricity = single(zeros(1,NumOfTracks);
clear TracksEdgeParameters
TracksEdgeParameters.Start.XYcoordinates  = single(zeros(NumOfTracks,2));
TracksEdgeParameters.End.XYcoordinates    = single(zeros(NumOfTracks,2));
TracksEdgeParameters.Start.Frame          = single(zeros(1,NumOfTracks));
TracksEdgeParameters.End.Frame            = single(zeros(1,NumOfTracks));

for tr_ind = 1:NumOfTracks
    TracksFramesMatrix(tr_ind,Tracks(tr_ind).Frames) = true;
    TracksEdgeParameters.Start.XYcoordinates(tr_ind,:) = Tracks(tr_ind).Path(1,:);
    TracksEdgeParameters.End.XYcoordinates(tr_ind,:)   = Tracks(tr_ind).Path(end,:);
    TracksEdgeParameters.Start.Frame(tr_ind)           = min(Tracks(tr_ind).Frames);
    TracksEdgeParameters.End.Frame(tr_ind)             = max(Tracks(tr_ind).Frames(end));    
end
% figure; imagesc(~TracksFramesMatrix); colormap gray
NumOfDetectedWorms  = sum(single(TracksFramesMatrix),1);
NumOfWorms          = prctile(NumOfDetectedWorms,99); % don't take 100% in case there are a few detection errors
FramesWithOnlyOneWormMissing = NumOfDetectedWorms==(NumOfWorms-1);
FramesWithAllWormsDetected   = NumOfDetectedWorms== NumOfWorms;

% that accidently add "worm-like" objects   
% MaxNumOfWorms       = max(NumOfDetectedWorms);
% MinNumOfWorms       = min(NumOfDetectedWorms);
% N                   = histc(NumOfDetectedWorms,((MinNumOfWorms-0.5):(MaxNumOfWorms+0.5)));
% figure; bar(MinNumOfWorms:MaxNumOfWorms,N(1:end-1)) 

% Fix truncated tracks without conflicting events
TracksLinked_ToPreviousTracks    = false(1,NumOfTracks);
GoodLinks                        = [];

for frame = 1:NumOfFrames  
    % find tracks that ends at this frame
    truncated_tracks   = find(TracksEdgeParameters.End.Frame==frame);    
    if (length(truncated_tracks)==1) && (FramesWithAllWormsDetected(frame))   % Only one track is truncated and at that frame all worms are detected
        CurrentActiveTracks           = find(TracksFramesMatrix(:,frame));
        NextFrameWithAllWormsDetected = frame+find(FramesWithAllWormsDetected((frame+1):end),1,'first');
        FramesInBetween               = (frame+1):(NextFrameWithAllWormsDetected-1);
        InBetweenFramesOnlyThisTrackIsMissing = all(FramesWithOnlyOneWormMissing(FramesInBetween));
        if InBetweenFramesOnlyThisTrackIsMissing
            NextTrackIndex = setdiff(find(TracksFramesMatrix(:,NextFrameWithAllWormsDetected)),CurrentActiveTracks);
            if (length(NextTrackIndex)==1) && (NextTrackIndex<NumOfTracks) && (~TracksLinked_ToPreviousTracks(NextTrackIndex))
                GoodLinks = [GoodLinks; [truncated_tracks, NextTrackIndex]];
                TracksLinked_ToPreviousTracks(NextTrackIndex) = true; 
            end
        end        
        
%         CurrentActiveTracks = find(TracksFramesMatrix(:,frame));
%         NextTrackIndex = max(CurrentActiveTracks)+1;
%         if NextTrackIndex>NumOfTracks
%             break
%         end
%         if TracksLinked_ToPreviousTracks(NextTrackIndex)  % make sure that tracks are not linked more than once !
%             continue
%         end
%         NextTrackStartFrame                      = TracksEdgeParameters.Start.Frame(NextTrackIndex);
%         FramesBetweenTracks                      = (frame+1):(NextTrackStartFrame-1);
%         BetweenTheTwoTracks_OnlyOneWormIsMissing = all(FramesWithOnlyOneWormMissing(FramesBetweenTracks));
%         
%         if BetweenTheTwoTracks_OnlyOneWormIsMissing
%             GoodLinks = [GoodLinks; [truncated_tracks, NextTrackIndex]];
%             TracksLinked_ToPreviousTracks(NextTrackIndex) = true; 
%         end
    end    
end

return

function Tracks = LinkGoodTracks (Tracks,GoodLinks, PossibleCollision)

NumOfLinks = size(GoodLinks,1);

for link_num =  NumOfLinks:(-1):1   % from the end to the beginning
    tr_ind1 = GoodLinks(link_num,1);
    tr_ind2 = GoodLinks(link_num,2);
    currentLastFrame                  = Tracks(tr_ind1).Frames(end);
    
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
    
    Tracks(tr_ind1).SkeletonLength              = [Tracks(tr_ind1).SkeletonLength,                  Tracks(tr_ind2).SkeletonLength] ;
    Tracks(tr_ind1).PerimeterLength             = [Tracks(tr_ind1).PerimeterLength,                 Tracks(tr_ind2).PerimeterLength] ;    
    Tracks(tr_ind1).WormPerimeter.Xcoordinate   = [Tracks(tr_ind1).WormPerimeter.Xcoordinate,       Tracks(tr_ind2).WormPerimeter.Xcoordinate];          % vector of indices per worm per frame
    Tracks(tr_ind1).WormPerimeter.Ycoordinate   = [Tracks(tr_ind1).WormPerimeter.Ycoordinate,       Tracks(tr_ind2).WormPerimeter.Ycoordinate];          % vector of indices per worm per frame                    
    Tracks(tr_ind1).WormPerimeter.LinearIndices = [Tracks(tr_ind1).WormPerimeter.LinearIndices,     Tracks(tr_ind2).WormPerimeter.LinearIndices];          % vector of indices per worm per frame
    if isfield(Tracks(tr_ind1),'WormArea')        
        Tracks(tr_ind1).WormArea.Xcoordinate        = [Tracks(tr_ind1).WormArea.Xcoordinate,            Tracks(tr_ind2).WormArea.Xcoordinate];          % vector of indices per worm per frame
        Tracks(tr_ind1).WormArea.Ycoordinate        = [Tracks(tr_ind1).WormArea.Ycoordinate,            Tracks(tr_ind2).WormArea.Ycoordinate];          % vector of indices per worm per frame                    
        Tracks(tr_ind1).WormArea.LinearIndices      = [Tracks(tr_ind1).WormArea.LinearIndices,          Tracks(tr_ind2).WormArea.LinearIndices];          % vector of indices per worm per frame                    
    end
    Tracks(tr_ind1).HeadTail                    = [Tracks(tr_ind1).HeadTail ;                       Tracks(tr_ind2).HeadTail];              
    Tracks(tr_ind1).Midline.X_coordinates_short = [Tracks(tr_ind1).Midline.X_coordinates_short,     Tracks(tr_ind2).Midline.X_coordinates_short];        
    Tracks(tr_ind1).Midline.Y_coordinates_short = [Tracks(tr_ind1).Midline.Y_coordinates_short,     Tracks(tr_ind2).Midline.Y_coordinates_short];        
    Tracks(tr_ind1).Midline.LinearIndices       = [Tracks(tr_ind1).Midline.LinearIndices,           Tracks(tr_ind2).Midline.LinearIndices];        
    Tracks(tr_ind1).Midline.Pattern             = [Tracks(tr_ind1).Midline.Pattern,                 Tracks(tr_ind2).Midline.Pattern];        
    Tracks(tr_ind1).Midline.X_coordinates       = [Tracks(tr_ind1).Midline.X_coordinates,           Tracks(tr_ind2).Midline.X_coordinates];        
    Tracks(tr_ind1).Midline.Y_coordinates       = [Tracks(tr_ind1).Midline.Y_coordinates,           Tracks(tr_ind2).Midline.Y_coordinates];        
    Tracks(tr_ind1).Midline.Angle               = [Tracks(tr_ind1).Midline.Angle,                   Tracks(tr_ind2).Midline.Angle];        
    Tracks(tr_ind1).Midline.FlagTrueIfReliable  = [Tracks(tr_ind1).Midline.FlagTrueIfReliable,      Tracks(tr_ind2).Midline.FlagTrueIfReliable];                           

    if isfield(Tracks(1).Midline,'AngleFirstPoint')
        Tracks(tr_ind1).Midline.AngleFirstPoint     = [Tracks(tr_ind1).Midline.AngleFirstPoint,         Tracks(tr_ind2).Midline.AngleFirstPoint];        
        Tracks(tr_ind1).Midline.AngleLastPoint      = [Tracks(tr_ind1).Midline.AngleLastPoint,          Tracks(tr_ind2).Midline.AngleLastPoint];        
        Tracks(tr_ind1).Midline.NumOfBends_HighRes  = [Tracks(tr_ind1).Midline.NumOfBends_HighRes,      Tracks(tr_ind2).Midline.NumOfBends_HighRes];        
        Tracks(tr_ind1).Midline.NumOfBends_LowRes   = [Tracks(tr_ind1).Midline.NumOfBends_LowRes,       Tracks(tr_ind2).Midline.NumOfBends_LowRes];        
    end
    
    if isfield(Tracks(tr_ind1),'Head')        
        Tracks(tr_ind1).Head              = [Tracks(tr_ind1).Head;        Tracks(tr_ind2).Head];
        Tracks(tr_ind1).Tail              = [Tracks(tr_ind1).Tail;        Tracks(tr_ind2).Tail];
    end
    if isfield(Tracks(tr_ind1),'NeuronInWormCoordinates')        
        Tracks(tr_ind1).NeuronInWormCoordinates = [Tracks(tr_ind1).NeuronInWormCoordinates;        Tracks(tr_ind2).NeuronInWormCoordinates];
    end
    if isfield(Tracks(tr_ind1),'FramesWithReliableHeadSide')        
        Tracks(tr_ind1).FramesWithReliableHeadSide = [Tracks(tr_ind1).FramesWithReliableHeadSide,  Tracks(tr_ind2).FramesWithReliableHeadSide];
%         Tracks(tr_ind1).HeadDirectionWasSwitched   = [Tracks(tr_ind1).HeadDirectionWasSwitched,    Tracks(tr_ind2).HeadDirectionWasSwitched];
    end
    if isfield(Tracks(tr_ind1),'FramesWithEstimatedHeadSide')    
%         Tracks(tr_ind1).FramesWithEstimatedHeadSide= [Tracks(tr_ind1).FramesWithEstimatedHeadSide, Tracks(tr_ind2).FramesWithEstimatedHeadSide];
    end
    
    if exist('PossibleCollision','var')
        if PossibleCollision
            Tracks(tr_ind1).PossibleCollision(currentLastFrame:(Tracks(tr_ind2).Frames(1))) = true;
        end
    end

end

TracksToRemove                 = false(1,length(Tracks));
TracksToRemove(GoodLinks(:,2)) = true;
Tracks                         = Tracks(~TracksToRemove);

return

function Cell_out = StitchCell_BasedOnGoodLinks (Cell_in,GoodLinks_ALL)

Cell_out                  = Cell_in;
NumberOfLinkingIterations = length(GoodLinks_ALL);

for link_iteration_ind = 1:NumberOfLinkingIterations
    GoodLinks  = GoodLinks_ALL{link_iteration_ind};    
    NumOfLinks = size(GoodLinks,1);

    for link_num =  NumOfLinks:(-1):1   % from the end to the beginning
        tr_ind1 = GoodLinks(link_num,1);
        tr_ind2 = GoodLinks(link_num,2);

        Cell_out{tr_ind1} = [Cell_out{tr_ind1}, Cell_out{tr_ind2}];    
    end

    TracksToRemove                 = false(1,length(Cell_out));
    TracksToRemove(GoodLinks(:,2)) = true;
    Cell_out                       = Cell_out(~TracksToRemove);
end

return

function Cell_out = StitchMatrices_BasedOnGoodLinks (Cell_in,GoodLinks_ALL, Non_Deleted_Tracks_BeforeStitching, Non_Deleted_Tracks_BeforeManualInspection, TracksCorrectionParameters, MaxNumOfWorms)

Cell_out                  = Cell_in(Non_Deleted_Tracks_BeforeStitching);
NumberOfLinkingIterations = length(GoodLinks_ALL);

for link_iteration_ind = 1:NumberOfLinkingIterations
    if (link_iteration_ind == NumberOfLinkingIterations) && (MaxNumOfWorms>1) % last link is manual inspection. before that all short tracks are deleted 
        Cell_out = Cell_out(Non_Deleted_Tracks_BeforeManualInspection);
    end
    
    GoodLinks  = GoodLinks_ALL{link_iteration_ind}; 
    if isempty(GoodLinks)
        disp('Empty ''GoodLinks'' was detected'); 
        continue;
    end
    NumOfLinks = size(GoodLinks,1);
%     disp([datestr(now), num2str(link_iteration_ind)])
    for link_num =  NumOfLinks:(-1):1   % from the end to the beginning
        tr_ind1 = GoodLinks(link_num,1);
        tr_ind2 = GoodLinks(link_num,2);

        Cell_out{tr_ind1} = [Cell_out{tr_ind1}; Cell_out{tr_ind2}];   
        Cell_out{tr_ind2} = [];
    end

    TracksToRemove                 = false(1,length(Cell_out));
    TracksToRemove(GoodLinks(:,2)) = true;
    Cell_out                       = Cell_out(~TracksToRemove);
end

if ~isempty(TracksCorrectionParameters)  % ONE SWITCH - If determined by manual inspection
    tr_ind1      = TracksCorrectionParameters(1);
    tr_ind2      = TracksCorrectionParameters(2);
    index_in_tr1 = TracksCorrectionParameters(3);
    index_in_tr2 = TracksCorrectionParameters(4);
%     frame        = TracksCorrectionParameters(5);

    Cell_out_OLD = Cell_out;
    Cell_out{tr_ind1} = [Cell_out_OLD{tr_ind1}(1:index_in_tr1,:,:); Cell_out_OLD{tr_ind2}((index_in_tr2+1):end,:,:)]; 
    Cell_out{tr_ind2} = [Cell_out_OLD{tr_ind2}(1:index_in_tr2,:,:); Cell_out_OLD{tr_ind1}((index_in_tr1+1):end,:,:)];     
end

return

function [Tracks, TracksCorrectionParameters] = CorrectLinkingBetweenTwoTracks_inline (Tracks, tr1, tr2, frame)
% [Tracks, TracksCorrectionParameters] = CorrectLinkingBetweenTwoTracks_inline (Tracks, 1, 3, 40597)

Tracks_in = Tracks;

index_in_tr1 = find(Tracks(tr1).Frames<frame,1,'last');
index_in_tr2 = find(Tracks(tr2).Frames<frame,1,'last');
TracksCorrectionParameters = [tr1 tr2 index_in_tr1 index_in_tr2 frame];

FIELDNAMES = fieldnames(Tracks);
RelevantFieldsToCorrect = false(1,length(FIELDNAMES));
for f_ind = 1:length(FIELDNAMES)
    fieldname = FIELDNAMES{f_ind};
    CurrentNumberOfFrames = length(Tracks(tr1).Frames);
    RelevantFieldsToCorrect(f_ind) = ismember(CurrentNumberOfFrames,size(Tracks(tr1).(fieldname)));
end
FIELDNAMES_FirstField           = FIELDNAMES(RelevantFieldsToCorrect);
FIELDNAMES_midline              = fieldnames(Tracks(tr1).Midline);
FIELDNAMES_WormPerimeterAndArea = fieldnames(Tracks(tr1).WormPerimeter);

for f_ind = 1:length(FIELDNAMES_FirstField)
    fieldname  = FIELDNAMES_FirstField{f_ind};    
    CurrentVec = Tracks(tr1).(fieldname);
    if min(size(CurrentVec))==1
        Tracks(tr1).(fieldname) = [Tracks_in(tr1).(fieldname)(1:index_in_tr1), Tracks_in(tr2).(fieldname)((index_in_tr2+1):end)]; 
        Tracks(tr2).(fieldname) = [Tracks_in(tr2).(fieldname)(1:index_in_tr2), Tracks_in(tr1).(fieldname)((index_in_tr1+1):end)]; 
    elseif length(size(CurrentVec))<3   
        Tracks(tr1).(fieldname) = [Tracks_in(tr1).(fieldname)(1:index_in_tr1,:); Tracks_in(tr2).(fieldname)((index_in_tr2+1):end,:)]; 
        Tracks(tr2).(fieldname) = [Tracks_in(tr2).(fieldname)(1:index_in_tr2,:); Tracks_in(tr1).(fieldname)((index_in_tr1+1):end,:)];         
    elseif strcmpi(fieldname,'HeadTail')
        Tracks(tr1).(fieldname) = [Tracks_in(tr1).(fieldname)(1:index_in_tr1,:,:); Tracks_in(tr2).(fieldname)((index_in_tr2+1):end,:,:)]; 
        Tracks(tr2).(fieldname) = [Tracks_in(tr2).(fieldname)(1:index_in_tr2,:,:); Tracks_in(tr1).(fieldname)((index_in_tr1+1):end,:,:)];         
    else
        disp(fieldname,' -- unnkown format')
    end
end  

% Midline fields
for f_ind = 1:length(FIELDNAMES_midline)
    fieldname  = FIELDNAMES_midline{f_ind};
    Tracks(tr1).Midline.(fieldname) = [Tracks_in(tr1).Midline.(fieldname)(1:index_in_tr1), Tracks_in(tr2).Midline.(fieldname)((index_in_tr2+1):end)]; 
    Tracks(tr2).Midline.(fieldname) = [Tracks_in(tr2).Midline.(fieldname)(1:index_in_tr2), Tracks_in(tr1).Midline.(fieldname)((index_in_tr1+1):end)];  
end

% WormPerimeter 
for f_ind = 1:length(FIELDNAMES_WormPerimeterAndArea)
    fieldname  = FIELDNAMES_WormPerimeterAndArea{f_ind};
    Tracks(tr1).WormPerimeter.(fieldname) = [Tracks_in(tr1).WormPerimeter.(fieldname)(1:index_in_tr1), Tracks_in(tr2).WormPerimeter.(fieldname)((index_in_tr2+1):end)]; 
    Tracks(tr2).WormPerimeter.(fieldname) = [Tracks_in(tr2).WormPerimeter.(fieldname)(1:index_in_tr2), Tracks_in(tr1).WormPerimeter.(fieldname)((index_in_tr1+1):end)];  
end
% WormArea
for f_ind = 1:length(FIELDNAMES_WormPerimeterAndArea)
    fieldname  = FIELDNAMES_WormPerimeterAndArea{f_ind};
    Tracks(tr1).WormArea.(fieldname) = [Tracks_in(tr1).WormArea.(fieldname)(1:index_in_tr1), Tracks_in(tr2).WormArea.(fieldname)((index_in_tr2+1):end)]; 
    Tracks(tr2).WormArea.(fieldname) = [Tracks_in(tr2).WormArea.(fieldname)(1:index_in_tr2), Tracks_in(tr1).WormArea.(fieldname)((index_in_tr1+1):end)];  
end

return

function FramesAreUnique = CheckUniqueFrames(Tracks)

FramesAreUnique = true;
for tr=1:length(Tracks)
    Frames = Tracks(tr).Frames;
    if length(Frames)~=length(unique(Frames))
        FramesAreUnique = false;
        disp('Frames are not unique!!')
        return
    end
end

return

% Manual inspection
function GoodLinks = ManualCorrectionOfTrackLinkingEstimation (File, Tracks, EstimatedLinks)

% Filename:
FileNameAfterManualLinking       = [File.TrackFile(1:end-4),'_GoodLinksByManualInspection.mat'];

if exist(FileNameAfterManualLinking,'file')
    qstring = '''__GoodLinksByManualInspection.mat'' file was found. Load existing ''GoodLinks'' matrix or re-inspect Tracks linkage?';
    default = 'Load existing file';
    button = questdlg(qstring,'Manual inspection of tracks linkage','Load existing file','re-inspect Tracks linkage',default); 
    if strcmpi(button,'Load existing file')
        disp('loading existing manually inspected file');
        load(FileNameAfterManualLinking,'GoodLinks');
        return        
    end
end

EstimatedLinksToCheck      = EstimatedLinks;
uicontrol_xy               = [0.06 0.025]; 
uicontrol_xy_small         = [0.04 0.025]; 
done                       = false;
GoodLinks                  = [];

while ~done    
    NumberOfTrackLinksLeft       = size(EstimatedLinksToCheck,1);
    NumOfLinksToShow             = min([6, NumberOfTrackLinksLeft]);       
    LinksToShow                  = EstimatedLinksToCheck(1:NumOfLinksToShow,:);
    PreviousTracks               = LinksToShow(:,1);
    NextTracks                   = LinksToShow(:,2);   
    
    LastPreviousTrackFrames = zeros(1,NumOfLinksToShow);
    FirstNextTrackFrames    = zeros(1,NumOfLinksToShow);
    for link_num = 1:NumOfLinksToShow
        LastPreviousTrackFrames(link_num) = Tracks(PreviousTracks(link_num)).Frames(end); 
        FirstNextTrackFrames(link_num)    = Tracks(NextTracks(link_num)).Frames(1);         
    end
    Fig_handle=figure('position',get(0,'screensize'),'name','Track linkage correction. Correct track linkage if necessary then finish review');  
    for sp_ind = 1:NumOfLinksToShow
        subplot(2,3,sp_ind)
        ax=gca;
        CurrentFrame = LastPreviousTrackFrames(sp_ind);
        WormTrackingMovie([],[],File, CurrentFrame, Tracks, ax, PreviousTracks(sp_ind), NextTracks(sp_ind)); 
        P       = get(gca,'position');
        PstartX = P(1)+P(3)/2-uicontrol_xy(1);
        PstartY = P(2)-uicontrol_xy(2);
        ax=gca;
%         buttons(sp_ind)           = uicontrol(Fig_handle,'units','normalized','BackgroundColor','r','position',[PstartX-uicontrol_xy(1)/2  PstartY uicontrol_xy],'Style','popupmenu','String',{'Correct','Incorrect','not clear'});
        uicontrol(Fig_handle,'units','normalized','BackgroundColor','w','position',[PstartX+uicontrol_xy(1)*0.6  PstartY uicontrol_xy],'Style','text','String','Previous track:');
        uicontrol(Fig_handle,'units','normalized','BackgroundColor','w','position',[PstartX+uicontrol_xy(1)*0.6  PstartY-uicontrol_xy(2) uicontrol_xy],'Style','text','String','Next track:');
        buttons_previous(sp_ind)  = uicontrol(Fig_handle,'units','normalized','BackgroundColor','w','position',[PstartX+uicontrol_xy(1)*3/2  PstartY uicontrol_xy_small],'Style','edit','String',num2str(PreviousTracks(sp_ind)),'Value',PreviousTracks(sp_ind));
        buttons_next(sp_ind)      = uicontrol(Fig_handle,'units','normalized','BackgroundColor','w','position',[PstartX+uicontrol_xy(1)*3/2  PstartY-uicontrol_xy(2) uicontrol_xy_small],'Style','edit','String',num2str(NextTracks(sp_ind)),'Value',NextTracks(sp_ind));
        buttons2(sp_ind)          = uicontrol(Fig_handle,'units','normalized','BackgroundColor','r','position',[PstartX-uicontrol_xy(1)/2  PstartY-uicontrol_xy(2)/2 uicontrol_xy],'Style','pushbutton','String','movie',...
            'Callback', {@WormTrackingMovie,File, [], Tracks, ax, [], [], Fig_handle, sp_ind});
    end
    data.buttons_previous = buttons_previous;
    data.buttons_next     = buttons_next;
    guidata(Fig_handle,data);                                    % assign the data structure to the GUI
    button_finish = uicontrol(Fig_handle,'units','normalized','BackgroundColor','g','position',[0.0370    0.5040    0.0600    0.0250],'Style','popupmenu','String',{'Review in process','Finish review'});
    % Add here user input
    % if corrections are needed, use only corrected frames for further Head/Tail estimations  
    success = button_finish.Value==2;
    while ~success
        k = waitforbuttonpress;        
        success = button_finish.Value==2;
    end
    CorrectedPreviousTracks = zeros(1,NumOfLinksToShow);
    CorrectedNextTracks     = zeros(1,NumOfLinksToShow);
    for sp_ind = 1:NumOfLinksToShow    % correct linkage
        CorrectedPreviousTracks(sp_ind)   = str2double(buttons_previous(sp_ind).String);  
        CorrectedNextTracks(sp_ind)       = str2double(buttons_next(sp_ind).String);      
    end
    [CorrectedPreviousTracks' CorrectedNextTracks']
    disp(char(10))
    GoodLinks = [GoodLinks; [CorrectedPreviousTracks' CorrectedNextTracks']];
    
    if NumberOfTrackLinksLeft > NumOfLinksToShow
        EstimatedLinksToCheck = EstimatedLinksToCheck((NumOfLinksToShow+1):end,:);
    else
        done = true;
    end            
    close; 
    clear buttons  buttons2  buttons_previous  buttons_next    
end

[~,Indices]=sort(GoodLinks(:,2));
GoodLinks  = GoodLinks(Indices,:);

save(FileNameAfterManualLinking,'GoodLinks');
return

%% Midline
function Tracks = IdentifyMidlineErrors_IndividualWorms (Tracks, File)
% "clean" data for IDing and neuron tracking
plotme = true;

% Maximum deviation allowed (in [%]) from the median midline lengths that was found in all worms within the arena.                                                                                               
% This is used to flag out all Midlines that were 'too short' due to error in morphology calculations.  
MaximumDeviationFromMedianLength = File.VariablesInformation.MidlineCalculationParams.MaximumDeviationFromMedianLength;  
if plotme
    figure('name','Histograms of midline before corrections. Midline lengths below the line cutoff are flagged as unreliable'); 
end
for tr_ind = 1:length(Tracks)
    SkeletonLengths   = single([Tracks(tr_ind).SkeletonLength]);
    WormSizes         = single([Tracks(tr_ind).Size]);
    Midlines          = [Tracks(tr_ind).Midline];
    ReliableSkeletons = [Midlines.FlagTrueIfReliable];
    SkeletonLengths   = SkeletonLengths(ReliableSkeletons);
    WormSizes         = WormSizes(ReliableSkeletons);

    MaxSize           = prctile(WormSizes,99.8);
    MaxSkeletonLength = prctile(SkeletonLengths,99.8);
    MinSkeletonLength = prctile(SkeletonLengths,60)*(1-MaximumDeviationFromMedianLength/100);
    
    BadIndices = (Tracks(tr_ind).SkeletonLength < MinSkeletonLength)|(Tracks(tr_ind).SkeletonLength > MaxSkeletonLength)|...
                 (Tracks(tr_ind).Size > MaxSize);
    Tracks(tr_ind).Midline.FlagTrueIfMidlineWasDetected   = Tracks(tr_ind).Midline.FlagTrueIfReliable;
    Tracks(tr_ind).Midline.FlagTrueIfReliable(BadIndices) = false;    
    Tracks(tr_ind).TrackWormMidlineLength = median(single(Tracks(tr_ind).SkeletonLength(Tracks(tr_ind).Midline.FlagTrueIfReliable)));    
    Tracks(tr_ind).TrackWormSize          = median(single(Tracks(tr_ind).Size(Tracks(tr_ind).Midline.FlagTrueIfReliable)));
        
    if plotme
        subplot(length(Tracks),1,tr_ind)
        title(['Track ',num2str(tr_ind)]); 
        Xedges = 50:5:250;
        X      = (Xedges(1:end-1)+Xedges(2:end))/2;        
        N = histc(SkeletonLengths,Xedges); bar(X,N(1:end-1));
        line([MinSkeletonLength MinSkeletonLength],get(gca,'ylim'));
        line([MaxSkeletonLength MaxSkeletonLength],get(gca,'ylim'));
        xlabel('Midline length')
        ylabel('Count')
    end

end

return

function WormTrackingMovie(stam, stam2, File, SpecificFrames, Tracks, gca_handle, PreviousTrack, NextTrack, Fig_handle, sp_ind)

if exist('Fig_handle','var')
    data          = guidata(Fig_handle);
    PreviousTrack = str2double(data.buttons_previous(sp_ind).String);
    NextTrack     = str2double(data.buttons_next(sp_ind).String);
end

if isempty(SpecificFrames)
    LastPreviousTrackFrame = Tracks(PreviousTrack).Frames(end); 
    FirstNextTrackFrame    = Tracks(NextTrack).Frames(1); 
    NumOfFramesToShow      = 100; 

    FirstFrameToShow      = max([1,                   (LastPreviousTrackFrame-2*File.FrameRate)]); 
    LastFrameToShow       = min([File.NumberOfFrames, (FirstNextTrackFrame   +2*File.FrameRate)]); 
    FrameInterval         = max([1,(LastFrameToShow-FirstFrameToShow)/NumOfFramesToShow]);
    SpecificFrames        = round(FirstFrameToShow:FrameInterval:LastFrameToShow); 
    if isempty(SpecificFrames)
        disp('No frames to show');
        return
    end
end



axes(gca_handle); 
% origin: WormTracking_MovieDisplay_Imaging_tracker_v01  

%% Parameters:

MovieOptions.ZoomOnWorm                   = false; %%% Assuming there is only One Worm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
MovieOptions.ShowCTXposition              = false;                     
MovieOptions.ShowBehaviorFrom_Data_Matrices = false;
MovieOptions.ShowTracksOnBackground       = false;     % Default=false; Tracks are shown by default on the original movie frames that are read from the path given in 'File'

% Other paremeters
MovieOptions.deltaT_Collision             = 0.05;      % Display time in regular frames
MovieOptions.deltaT                       = 0.05;      % Display time in regular frames
% MovieOptions.arena                        = 1;     
if exist('ArenaID','var');
    MovieOptions.arena                      = ArenaID;     
end
if exist('SpecificFrames','var')
    if ~isempty(SpecificFrames)
        MovieOptions.ShowSpecificFrames = SpecificFrames;
    end
end

MovieOptions.TracksLine                   = false;      % true/false/Num = show/(don't show)/(Show last Num of frames) of track line preceeding to each frame.                                                         
MovieOptions.TracksText.TrackNum          = true;      % Track text diplay settings
MovieOptions.TracksText.ID1               = false;     
MovieOptions.TracksText.ID2               = false;     
MovieOptions.TracksText.Behavior         = false; 
% MovieOptions.TracksText.OnlyBehaviorColor_NoText = false;   % true or false

% MovieOptions.TracksText.ColorByID         = 0;     % No, same colors for all worms. different colors based on features 
MovieOptions.TracksText.ColorByID         = 2;     % Color set by track ID

% MovieOptions.ShowDetectedPixels           = {'Area','Skeleton','HeadTail'}; % A cell with the name of ields to display
% MovieOptions.ShowDetectedPixelsColors     = {'b','r','g'};                  % Color for display
MovieOptions.ShowDetectedPixels           = {'Perimeter'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'b'};                  % Color for display
MovieOptions.ShowDetectedPixels           = {'Perimeter','Midline'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'b','g'};                  % Color for display
MovieOptions.ShowDetectedPixels           = {'Perimeter','Midline','Head','Tail'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'b','g','r','m'};                      % Color for display
MovieOptions.ShowDetectedPixels           = {'Midline','Head','Tail'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'g','r','m'};                      % Color for display
MovieOptions.ShowDetectedPixels           = {'Perimeter','Midline','HeadTail'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'g','r','m'};                      % Color for display

BehaviorStructure = [];
            
% MovieOptions.ConstantScale                = true;        % If true, the function will look for the best scaling and scale each frame similarly. If false --> imshow(Mov,[]) 
MovieOptions.ConstantScale                = 2;        % If true, the function will look for the best scaling and scale each frame similarly. If false --> imshow(Mov,[]) 
% MovieOptions.ConstantScale                = false;     
MovieOptions.CreateMovie                  = false;       % If true create a movie
Data_Matrices = [];
background = [];

SimpleDisplay_FluorescenceMovieWithTracks1(Tracks, File, MovieOptions, Data_Matrices, BehaviorStructure, background, PreviousTrack, NextTrack);  % Optional Movie generation

return

function SimpleDisplay_FluorescenceMovieWithTracks1 (Tracks, File, MovieOptions, Data_Matrices, BehaviorStructure, background, PreviousTrack, NextTrack)
% Display a movie with Tracks
% Information about the Tracks and the movie is given in the Tracks and File input variables  
% MovieOptions-  A stucture with infomation about what and how to display, with the following fields:   
%    deltaT:             Pause between frames in seconds (Default = 5msec).    
%    deltaT_collision:   Pause between frames near end of tracks in seconds (Default = 1 sec).    
%    arena:              Zoom on one specific arena (Default = no Zoom).       
%                        This field is NOT used if specific tracks are shown (arena zoom will be set by the tracks) 
%    EnlargeArena:       Enlarge Arena zoom by the specified number of pixels (Default = no enlargement).       
%    ShowSpecificTracks: A vector with Tracks of interest to display. Frames will be set by Tracks (Default = Show all tracks within the frames of interest).   
%                        Other Tracks within the same region of interest will be shown as well. 
%    ShowEndOfTracks:    A structure that, if exists, specifiy how many frames to display before the end of the Track (Default= show all the track)
%                        This is convenient to see worms before and after collisions  
%         FramesBefore   Number of frames before the first track ends (before first collision).  
%         FramesAfter    Number of frames after the start of the last detected track (after last collision).            
%    ShowSpecificFrames: A vector of Frames of interest to display. All relevant tracks will be shown(Default = All frames in the movie). 
%    TracksText:         Structure defining what fields to display. Values are true or false (Default = show only track number). 
%         TrackNum       if exists and true, show the track number (#)
%         ID1            Worm ID == original track calculated from linking analysis by Dirk (DL).  
%         ID2            Worm ID == original track calculated from my linking analysis (SL).  
%         Behavior      Worm behavior as calculated in the segmentation functions.   
% MovieOptions.TracksLine

warning off; 
MoviePath      = File.MoviePath;
MovieFileNames = File.MovieFileNames;

%% Movie format
if     strcmpi(File.VideoFormat,'tiff')
    VideoFormat   = 1;                                  % read new file every loop iteration 
elseif strcmpi(File.VideoFormat,'multi-tiff')
    VideoFormat   = 2;
    FileFullName  = [MoviePath,'\',MovieFileNames{1}];  % one file
elseif strcmpi(File.VideoFormat,'avi')
    VideoFormat   = 3;
    FileFullName  = [MoviePath,'\',MovieFileNames{1}];  % one file
    MovieObj      = VideoReader(FileFullName);          % one object
elseif strcmpi(File.VideoFormat,'multi-avi')
    VideoFormat   = 4;
    PreviousMovieObjFile = 0;                           % read new file object if and when it is needed 
end


%% Find the frames of each track
LEN             = length([Tracks.Active]);
FirstFrames     = zeros(1,LEN);
LastFrames      = zeros(1,LEN);
ALLFRAMES       = false(File.FragmentFrames(end,end),LEN);
for tr=1:LEN
    CurrentFrames               = Tracks(tr).Frames;
    FirstFrames(tr)             = CurrentFrames(1);
    LastFrames(tr)              = CurrentFrames(end);
    ALLFRAMES(CurrentFrames,tr) = 1;
end

%% Set parameters for the movie:   show specific tracks or specific frames,  show a single arena or all arenas,  frame rate
if ~exist('MovieOptions','var')
    deltaT             = 0.005;  
    arena              = NaN;     
    Frames2Show        = 1:size(ALLFRAMES,1);  
    ShowTracksText     = false;  
    ShowDetectedPixels = false;
    ShowDetectedROIs   = false;
    ConstantScale      = false;
    CreateMovie        = false;
    ColorByID          = false;
    ShowTracksOnBackground = false;
    Add_CTXdata_And_SmoothPosition    = false;
    ShowBehaviorFrom_Data_Matrices    = false;
    OnlyBehaviorColor_NoText          = false;
else
    
    ShowTracksText = false;
    ColorByID      = false;
    OnlyBehaviorColor_NoText = false;
    if isfield(MovieOptions,'TracksText')
        field_names    = fieldnames(MovieOptions.TracksText);
        field_names    = field_names(~strcmpi(field_names,'ColorByID'));
        for f_ind = 1:length(field_names)
            if islogical(MovieOptions.TracksText.(field_names{f_ind}))
                ShowTracksText = ShowTracksText || MovieOptions.TracksText.(field_names{f_ind});  % at least one text field = true
            else
                ShowTracksText = true;
            end
        end
        if isfield(MovieOptions.TracksText,'ColorByID')
            ColorByID = MovieOptions.TracksText.ColorByID;
        end 
       
        if isfield(MovieOptions.TracksText,'OnlyBehaviorColor_NoText')   
            if MovieOptions.TracksText.OnlyBehaviorColor_NoText
                OnlyBehaviorColor_NoText = true;
            end
        end  
    end

    % delta T = [1/frame rate], in seconds 
    if isfield(MovieOptions,'deltaT')    
        deltaT = MovieOptions.deltaT;
    else
        deltaT = 0.005;    
    end
    if isfield(MovieOptions,'deltaT_Collision')    
        deltaT_collision = MovieOptions.deltaT_Collision;
    else
        deltaT_collision = 1;    
    end
    
    % Zoom on a specific arena
    if isfield(MovieOptions,'arena')     
        arena = MovieOptions.arena;
    else
        arena = NaN;    
    end
    
    % This field below exists if the user is interested only in specific tracks. Only relevant Frames and arenas will be shown.  
    if isfield(MovieOptions,'ShowSpecificTracks')  
        Tracks2Show    = MovieOptions.ShowSpecificTracks; 
%         if ~ MovieOptions.ZoomOnWorm
%             % which arena
%             TrackArena = [Tracks.Arena];
%             RelevantArenas = unique(TrackArena(Tracks2Show));
%             if length(RelevantArenas)==1
%                 arena = RelevantArenas;
%                 disp(['All relevant tracks are in arena number ',num2str(arena),char(10)]);
%             end
%         end
        
        % which frames              
        if isfield(MovieOptions,'ShowEndOfTracks')
            FramesBefore = MovieOptions.ShowEndOfTracks.FramesBefore;
            FramesAfter  = MovieOptions.ShowEndOfTracks.FramesAfter;            
            Framefirst   = max( min(LastFrames(Tracks2Show)) - FramesBefore , 1);                             % FramesBefore collisions or frame number one; 
            Framelast    = min( max(FirstFrames(Tracks2Show))+ FramesAfter  , File.FragmentFrames(end,end));  % FramesAfter collisions or last;             
            Frames2Show  = Framefirst:Framelast;
            if isempty(Frames2Show)  % No collisions 
                disp('no collisions ??. Showing All frames in the track ');
                Frames2Show = min(FirstFrames(Tracks2Show)) : max(LastFrames(Tracks2Show)); 
            end
        else
            try 
                Frames2Show = min(FirstFrames(Tracks2Show)) : max(LastFrames(Tracks2Show));      
            catch
                Frames2Show = FirstFrames(1):LastFrames(end);
            end
        end
        
    end
    
    % This field exists if user is interested only in Frames. All relevant Tracks will be shown.  
    if isfield(MovieOptions,'ShowSpecificFrames')  
        if isfield(MovieOptions,'ShowSpecificTracks')  
            disp(['Both specific tracks and specific frames are defined. Showing the movie by the frames definition',char(10)]);            
        end
        Frames2Show = MovieOptions.ShowSpecificFrames;  % 
    elseif ~isfield(MovieOptions,'ShowSpecificFrames') && ~isfield(MovieOptions,'ShowSpecificTracks')
        All_Frames = unique([Tracks.Frames]);
        Frames2Show = All_Frames;    
    end
    
    % The relevant arena Zoom area: 
    if ~isnan(arena)
        ArenaBox =  File.TrackBoxAxis(arena,:);
        if isfield(MovieOptions,'EnlargeArena')
            Enlarge = MovieOptions.EnlargeArena;
            ArenaBox = [ArenaBox(1)-Enlarge, ArenaBox(2)+Enlarge, ArenaBox(3)-Enlarge, ArenaBox(4)+Enlarge];
        end
    end
    
    if isfield(MovieOptions,'TracksLine')  
        TracksLine = MovieOptions.TracksLine;
    else
        TracksLine = false;
    end   
   
    if isfield(MovieOptions,'ShowDetectedPixels')  
        ShowDetectedPixels = true;
        WormFeatures       = MovieOptions.ShowDetectedPixels;        % A cell with the name of ields to display
        WormFeaturesCOLORS = MovieOptions.ShowDetectedPixelsColors;  % Color for display
        OriginalFrameSize  = File.FrameSize;

    else
        ShowDetectedPixels = false;
    end
    
    if isfield(MovieOptions,'ShowDetectedROI')  
        ShowDetectedROIs   = true;
        load(File.FluorescenceFiles{1}, 'Perimeter_of_ROI', 'CorrelationWithReference', 'Frames', 'NonReliableIndices','GCaMP_Summary');
        MinCorrelation     = GCaMP_Summary(1).MinCorrelationThreshold;
        ROI.Frames         = Frames;
        ROI.Reliable       = (CorrelationWithReference > MinCorrelation) & (~NonReliableIndices);
        [~ , index_in_Frames2Show, index_in_ROI_Matrices] = intersect(Frames2Show, ROI.Frames );        %  [C,ia,ib] = intersect(A,B);
        
        [~,ReliablePart,~] = intersect(index_in_ROI_Matrices, find(ROI.Reliable));
        
        index_in_ROI_Matrices = index_in_ROI_Matrices(ReliablePart);
        index_in_Frames2Show  = index_in_Frames2Show(ReliablePart);
        
        % This vectors are of length (Frames2Show) and relate to the indices in Frames2Show (frame_number) 
        ROI.Frames2Show.ShowOrNot                              = false(1,length(Frames2Show)); 
        ROI.Frames2Show.ShowOrNot       (index_in_Frames2Show) = true;
        ROI.Frames2Show.index_in_ROI_Mat                       = zeros(1,length(Frames2Show));
        ROI.Frames2Show.index_in_ROI_Mat(index_in_Frames2Show) = index_in_ROI_Matrices;
        
        clear NonReliableIndices  Frames  CorrelationWithReference
        
        NumberOfROIs       = length(Perimeter_of_ROI);        % A cell with the name of ields to display
        ROI_COLORS         = MovieOptions.ShowDetectedROIColors;  % Color for display

    else
        ShowDetectedROIs   = false;
    end

    if isfield(MovieOptions,'ConstantScale')
        ConstantScale = MovieOptions.ConstantScale;
    end
    if ConstantScale
        Max = 0;
        Min = inf;
        for Frame = Frames2Show(1:end/20:end)
            switch VideoFormat
            case 2                             % A single multiple-tiff file
                Mov           = imread(FileFullName, Frame);     
            case 1                             % A sequence of tiff files
                FileFullName  = [MoviePath,'\',MovieFileNames{Frame}]; 
                Mov           = imread(FileFullName);                          
            case 3                             % A single avi movie file
                Mov           = read(MovieObj, Frame);     
                Mov           = Mov(:,:,1);                %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!
            case 4                             % A sequence of avi movie files
                CurrentMovieObjFile = File.MultiAviFrameConversion.MovieFileNumber(Frame);
                if (PreviousMovieObjFile ~= CurrentMovieObjFile)        % create a new file object if necessary
                    FileFullName        = [MoviePath,'\',MovieFileNames{CurrentMovieObjFile}]; 
                    PreviousMovieObjFile = CurrentMovieObjFile;
                    MovieObj             = VideoReader(FileFullName);
                end
                CurrentFrameNumberInFile = File.MultiAviFrameConversion.FrameNumberInFile(Frame);
                Mov                      = read(MovieObj, CurrentFrameNumberInFile);     
                Mov                      = Mov(:,:,1);          %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!                                
            end 
                                  
            Max = max([Max, max(Mov(:))]);
            Min = min([Min, min(Mov(:))]);
        end
        
        BestScale = [Min Max];    
        if isnumeric(ConstantScale)   
            if (ConstantScale<=100)&&(ConstantScale>=0)   % Scale based on percentage given by the user
%                 disp(['re-setting scale (',num2str(ConstantScale),'%)']);
%                 HalfRange = (Max-Min)/2;
%                 Middle    = (Max+Min)/2;                           
%                 BestScale = [Middle-(HalfRange*ConstantScale/100) Middle+(HalfRange*ConstantScale/100)];        
%                 BestScale = prctile(double(Mov(:)),[ConstantScale 100]);
%                 BestScale = prctile(double(Mov(:)),[ConstantScale 100-ConstantScale]);
                BestScale = prctile(double(Mov(:)),[ConstantScale 99.99]);
            end 
        end
    end
    
    CreateMovie = false;
    if isfield(MovieOptions,'CreateMovie')   
        if MovieOptions.CreateMovie
            CreateMovie = true;
            writerObj = VideoWriter('K:\testvideowriter1.avi','Uncompressed AVI'); 
            writerObj.FrameRate = 6;
            open(writerObj);
        end
    end
    
    ZoomOnWorm = false;
    if isfield(MovieOptions,'ZoomOnWorm')   
        if MovieOptions.ZoomOnWorm
            ZoomOnWorm   = true;
        end
    end 
    
    ShowTracksOnBackground=false;
    if isfield(MovieOptions,'ShowTracksOnBackground')   
        if MovieOptions.ShowTracksOnBackground
            ShowTracksOnBackground   = true;
        end
    end  
    
    if exist('BehaviorStructure','var')
        if ~isempty(BehaviorStructure)
            BehaviorFieldName = BehaviorStructure.fieldname_inTracks;
        end
    end
    
    Add_CTXdata_And_SmoothPosition = false;
    if isfield(MovieOptions,'ShowCTXposition')   
        if MovieOptions.ShowCTXposition
            if ~isempty(Data_Matrices)
                Add_CTXdata_And_SmoothPosition = true;
            end
        end
    end
    if isfield(MovieOptions,'ShowValueOfThisField')   
        if ~isempty(MovieOptions.ShowValueOfThisField)
            if ~isempty(Data_Matrices)
                Add_CTXdata_And_SmoothPosition = true;
                FramesField = 'FrameNumber';
                X_field     = 'Coordinates_X_Smoothed';
                Y_field     = 'Coordinates_Y_Smoothed';
                CTX_field = MovieOptions.ShowValueOfThisField;                
                disp(['Showing values of ', CTX_field]); 
             end
        end
    end    
    
    ShowBehaviorFrom_Data_Matrices = false;
    if isfield(MovieOptions,'ShowBehaviorFrom_Data_Matrices')   
        if MovieOptions.ShowBehaviorFrom_Data_Matrices
            if ~isempty(Data_Matrices)
                ShowBehaviorFrom_Data_Matrices = true;
                if ~exist('CTX_field','var') % if it was not previously defined in 'MovieOptions.ShowValueOfThisField' then show 'CTX_position'
                    if isfield(Data_Matrices,'FrameNumber_Smoothed2')
                        FramesField = 'FrameNumber_Smoothed2';
                        X_field     = 'X_Smoothed2';
                        Y_field     = 'Y_Smoothed2';
                        CTX_field   = 'CTX_position_Smoothed2';
                    elseif isfield(Data_Matrices,'FrameNumber')
                        disp('Showing ''raw'' CTX data --> The shown CTX  value are not yet corrected to the butanone odor side at that experiment')
                        FramesField = 'FrameNumber';
                        X_field     = 'Coordinates_X_Smoothed';
                        Y_field     = 'Coordinates_Y_Smoothed';
                        CTX_field   = 'CTX_position';
                    end    
                end
            end
        end
    end  
       
    
end

% if ColorByID
%     COLORS_FOR_TRACKS   = {'w','g','b','m','c','y',[0.8 0.5 0.2],[0.8 0.5 1]}; % it will start from the green
%     if ColorByID==2
%         TracksColorsIndices = mod(1:length(Tracks),length(COLORS_FOR_TRACKS))+1;  
%     elseif ColorByID==1
%         TracksColorsIndices = mod(1:length(unique([Tracks.ID])),length(COLORS_FOR_TRACKS))+1;  
%     end    
% end

if ColorByID==2
    COLORS_FOR_TRACKS   = {'g','r','w','b','m','c','y',[0.8 0.5 0.2],[0.8 0.5 1]}; % it will start from the green
    TracksColorsIndices                             = ones(1,length(Tracks));  % all tracks are the same color
    TracksColorsIndices([PreviousTrack, NextTrack]) = 2;                       % except of the one that are checked for linkage
end

if ShowTracksOnBackground
    if isempty(background)
        load(File.BackgroundFile,'background');
    end
    Mov = background;
end


%% Show Movie with tracks and their linking
% f=figure('position',get(0,'screensize'));  
% tic
for frame_number = 1:length(Frames2Show)
    Frame = Frames2Show(frame_number);
    
    if ~ShowTracksOnBackground        
        % Get Frame
        switch VideoFormat
            case 2                             % A single multiple-tiff file
                Mov           = imread(FileFullName, Frame);     
            case 1                             % A sequence of tiff files
                FileFullName  = [MoviePath,'\',MovieFileNames{Frame}]; 
                Mov           = imread(FileFullName);                          
            case 3                             % A single avi movie file
                Mov           = read(MovieObj, Frame);     
                Mov           = Mov(:,:,1);                %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!
            case 4                             % A sequence of avi movie files
                CurrentMovieObjFile = File.MultiAviFrameConversion.MovieFileNumber(Frame);
                if (PreviousMovieObjFile ~= CurrentMovieObjFile)        % create a new file object if necessary
                    FileFullName        = [MoviePath,'\',MovieFileNames{CurrentMovieObjFile}]; 
                    PreviousMovieObjFile = CurrentMovieObjFile;
                    MovieObj             = VideoReader(FileFullName);
                end
                CurrentFrameNumberInFile = File.MultiAviFrameConversion.FrameNumberInFile(Frame);
                Mov                      = read(MovieObj, CurrentFrameNumberInFile);     
                Mov                      = Mov(:,:,1);          %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!                                
        end 
    end
    
    if ConstantScale
        imshow(Mov,BestScale,'initialmagnification',600);
    else
        imshow(Mov,[],'initialmagnification',600);
    end
%     imshow(Mov,[],'initialmagnification',600);
     title(['Tracks linking: ',num2str(PreviousTrack),' - ', num2str(NextTrack),', Frame ',num2str(Frame)]);
   
%     if ~ Add_CTXdata_And_SmoothPosition  
%         if ~isnan(arena)        
%             axis(ArenaBox); 
%             title(['Arena ',num2str(arena),', Frame ',num2str(Frame)]);
%         else
%             title(['Frame ',num2str(Frame)]);
%         end
%     else
%         if ~isnan(arena)        
%             axis(ArenaBox); 
%         end
%     end
    if ZoomOnWorm   %%% Assuming there is only One Worm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
        FrameIndexInTrack = find([Tracks.Frames]==Frame,1);
        if ~isempty(FrameIndexInTrack)   
%             load(File.TrackingVariablesFile,'MaxWorm_size')
%             HalfWormMaxSize  = File.VariablesInformation.MaxWorm_size/0.05 * File.PixelSize / 2;      % File.VariablesInformation.MaxWorm_size ~ 1mm*0.05 mm
%             WormCentroid = Tracks.Path(FrameIndexInTrack,:);
%             Yaxis(1)   = max(1,           WormCentroid(1)-HalfWormMaxSize*1.5); 
%             Yaxis(2)   = min(size(Mov,1),WormCentroid(1)+HalfWormMaxSize*1.5); 
%             Xaxis(1)   = max(1,           WormCentroid(2)-HalfWormMaxSize*1.5); 
%             Xaxis(2)   = min(size(Mov,2),WormCentroid(2)+HalfWormMaxSize*1.5); 
%             Xaxis = round(Xaxis);
%             Yaxis = round(Yaxis);
%             axis([Xaxis Yaxis]);
%             axis([Yaxis Xaxis]);
            Xaxis = [min(Tracks.WormPerimeter.Xcoordinate{FrameIndexInTrack}) max(Tracks.WormPerimeter.Xcoordinate{FrameIndexInTrack})];
            Yaxis = [min(Tracks.WormPerimeter.Ycoordinate{FrameIndexInTrack}) max(Tracks.WormPerimeter.Ycoordinate{FrameIndexInTrack})];
            axis([Yaxis Xaxis]);
            
        end
    end
    
    if exist('Tracks2Show','var')
        RelevantTracks = Tracks2Show;
    else
        if Add_CTXdata_And_SmoothPosition
            RelevantTracks = 1:length(Tracks);
        else
            RelevantTracks = find(ALLFRAMES(Frame,:));
        end
    end

    hold on;
    All_CTX_position = [];
    for t_ind = 1:length(RelevantTracks);        
        t_num       = RelevantTracks(t_ind);
        IndexInTracks = find(Tracks(t_num).Frames==Frame,1);
        if isempty(IndexInTracks)
            Coordinates = [];
        else
            Coordinates = double(Tracks(t_num).Path(IndexInTracks,:)); 
            % Flip X-Y for imshow
            Coordinates = Coordinates([2 1]);
        end
        FrameCount  = Frame - FirstFrames(t_num)+1;
        
        % Add Track line of previous frames
        if TracksLine
            RelevantPathFrames = find(Tracks(t_num).Frames<=Frame);                      % Show all previous track path            
            if ~islogical(TracksLine) && (TracksLine < length(RelevantPathFrames)) % Show only the specified number of previous frames)
                RelevantPathFrames = RelevantPathFrames((end-TracksLine+1):end);       
            end
            Path = Tracks(t_num).Path(RelevantPathFrames,:);           
            line(Path(:,2), Path(:,1),'linewidth',1,'color','k','linestyle','-');
        end
        % ColumnIndex in Data_Matrices
        if Add_CTXdata_And_SmoothPosition  || ShowBehaviorFrom_Data_Matrices
             ColumnIndex       = find(Data_Matrices.(FramesField)(t_num,:) == Frame,1);   
            if isempty(ColumnIndex) 
%                 ColumnIndex       = find(Data_Matrices.FrameNumber(t_num,:) == (Frame-1),1);   
                ColumnIndex       = find(Data_Matrices.(FramesField)(t_num,:) == (Frame-1),1);   
            end     
        end

        % Defining display text for each track
        if ShowTracksText && ~isempty(Coordinates)
            text_for_Movie = '';
            if MovieOptions.TracksText.TrackNum
                text_for_Movie = [text_for_Movie,'(#',num2str(t_num),')'];
            end
            if MovieOptions.TracksText.ID1
                text_for_Movie = [text_for_Movie,',L',num2str(Tracks(t_num).OriginalTrack)];
            end
            if MovieOptions.TracksText.ID2
                text_for_Movie = [text_for_Movie,',LL',num2str(Tracks(t_num).Worm_ID)];
            end
            if MovieOptions.TracksText.Behavior & ~OnlyBehaviorColor_NoText
                if ShowBehaviorFrom_Data_Matrices &&  ~isempty(ColumnIndex)
                    CurrentBehCode =  Data_Matrices.(BehaviorFieldName)(t_num,ColumnIndex);
                else
                    CurrentBehCode =  Tracks(t_num).(BehaviorFieldName)(IndexInTracks);
                end
%                 CurrentBehColor=  BehaviorStructure.BehaviorCodeColors{BehaviorStructure.BehaviorCodeNumbers==CurrentBehCode};
                CurrentBehName =  BehaviorStructure.BehaviorCodeName{BehaviorStructure.BehaviorCodeNumbers==CurrentBehCode};
                if isempty(text_for_Movie)
                    text_for_Movie = CurrentBehName;
                else                    
                    text_for_Movie = [text_for_Movie,char(10),CurrentBehName];
                end
            
            end 
        else
%             text_for_Movie = num2str(t_num);
            text_for_Movie = '';
        end
        
        if Add_CTXdata_And_SmoothPosition  &&  ~isempty(ColumnIndex)
            Current_Coordinates_X_Smoothed = double(Data_Matrices.(X_field)(t_num,ColumnIndex));
            Current_Coordinates_Y_Smoothed = double(Data_Matrices.(Y_field)(t_num,ColumnIndex));
            Current_CTX_position           = Data_Matrices.(CTX_field)(t_num,ColumnIndex);
            All_CTX_position = [All_CTX_position, Current_CTX_position];

            CTXtext  = num2str(Current_CTX_position,2);
%                     % From G --> R through brown
%                     R        = abs((Current_CTX_position+1)/2); 
%                     G        = abs((Current_CTX_position-1)/2); 
%                     B        = 0; 
            % From G --> R through black
            if Current_CTX_position>0
                R = Current_CTX_position;
                G = 0;
            else
                G = -Current_CTX_position;
                R = 0;
            end
            B        = 0; 
            CTXcolor = [R G B];            
        else
            CTXtext  = '';
            CTXcolor = 'b';                               
        end                         
    
        
                
        CurrentColor = 'b';
        if ~isempty(Coordinates)
            if ColorByID ==2      % Color set by track number
                CurrentColor = COLORS_FOR_TRACKS{TracksColorsIndices(t_num)};
            elseif ColorByID ==1  % Color set by worm ID
                Worm_ID = Tracks(t_num).ID;
                CurrentColor = COLORS_FOR_TRACKS{TracksColorsIndices(Worm_ID)};
            elseif ColorByID ==3  % Color set by Chemotaxis index
                CurrentColor = CTXcolor;
            elseif ColorByID ==4  % Color set by Behavior code
                if ShowBehaviorFrom_Data_Matrices &&  ~isempty(ColumnIndex)
                    CurrentBehCode =  Data_Matrices.(BehaviorFieldName)(t_num,ColumnIndex);
                else
                    CurrentBehCode =  Tracks(t_num).(BehaviorFieldName)(IndexInTracks);
                end
                CurrentBehColor=  BehaviorStructure.BehaviorCodeColors{BehaviorStructure.BehaviorCodeNumbers==CurrentBehCode};
                CurrentColor   = CurrentBehColor;
                CTXcolor       = CurrentColor;
            end
        else
            if ColorByID ==4
                CTXcolor = 'y';
            end                
        end  
        
        if Add_CTXdata_And_SmoothPosition && ~isempty(CTXtext)
            text(Current_Coordinates_Y_Smoothed-30,Current_Coordinates_X_Smoothed-30, CTXtext,'FontWeight','bold','color',CTXcolor);             % Track count
            plot(Current_Coordinates_Y_Smoothed   ,Current_Coordinates_X_Smoothed,'o','markersize',6,'color',CTXcolor,'markerfacecolor',CTXcolor); 
        end
        if ~isempty(Coordinates)
            text(Coordinates(2)+10,Coordinates(1),text_for_Movie,'FontWeight','bold','color',CurrentColor);             % Track count
            plot(Coordinates(2)   ,Coordinates(1),'.','markersize',12,'color',CurrentColor); 
        end
        
%         if Add_CTXdata_And_SmoothPosition  
%             if ~isnan(arena)                       
%                 title(['Arena ',num2str(arena),', Frame ',num2str(Frame),', CTX ',num2str(mean(All_CTX_position),2)]);
%             else
%                 title(['Frame ',num2str(Frame),', CTX ',num2str(mean(All_CTX_position),2)]);
%             end
%         end        
        %% Plot detected pixels
        if  ShowDetectedPixels  && ~isempty(Coordinates)
            for f_num = 1:length(WormFeatures)
                CurrentField = WormFeatures{f_num};        % name of field to display
                CurrentColor = WormFeaturesCOLORS{f_num};  % color
                if ColorByID ==2      % Color set by track number
                    CurrentColor = COLORS_FOR_TRACKS{TracksColorsIndices(t_num)};
                elseif ColorByID ==1  % Color set by worm ID
                    Worm_ID = Tracks(t_num).ID;
                    CurrentColor = COLORS_FOR_TRACKS{TracksColorsIndices(Worm_ID)};
                elseif ColorByID ==3  % Color set by Chemotaxis index
                    CurrentColor = CTXcolor;                    
                elseif ColorByID ==4  % Color set by Behavior code
                    CurrentColor = CurrentBehColor;                    
                end                         
                Plot_WormPixels_on_Frame1(Tracks(t_num), Frame, CurrentField, CurrentColor, ZoomOnWorm)              
                
            end
        end       
        %% Plot detected ROI
        if  ShowDetectedROIs                       
            if ROI.Frames2Show.ShowOrNot(frame_number) 
                Current_frame_ind_in_ROI = ROI.Frames2Show.index_in_ROI_Mat(frame_number);
                for ROI_ind = 1:NumberOfROIs                
                    plot(Perimeter_of_ROI{ROI_ind}(Current_frame_ind_in_ROI,:,2), ...
                         Perimeter_of_ROI{ROI_ind}(Current_frame_ind_in_ROI,:,1),...
                         '.','color', ROI_COLORS{ROI_ind}, 'markersize',0.2)  ;  hold on
                end
            end
        end       
    end
    
    if Add_CTXdata_And_SmoothPosition  
        if ~isnan(arena)                       
            title(['Arena ',num2str(arena),', Frame ',num2str(Frame),', CTX ',num2str(mean(All_CTX_position),2)]);
        else
            title(['Frame ',num2str(Frame),', CTX ',num2str(mean(All_CTX_position),2)]);
        end
    end
        
    % Allow long pauses before a track ends and when it begins (within one frame)  
    if ismember(Frame,[LastFrames(RelevantTracks) FirstFrames(RelevantTracks)])
% pause;
        pause(deltaT_collision);
    elseif ismember(Frame+1,LastFrames(RelevantTracks)) || ismember(Frame-1, FirstFrames(RelevantTracks))
        pause(deltaT_collision);
    else
        pause(deltaT);
    end
    
    if CreateMovie    
        set(f,'position',get(0,'screensize'));  
%         writeVideo(writerObj,getframe);
        writeVideo(writerObj,getframe(f));
    end
        
    hold off;
end

if CreateMovie    
    close(writerObj);
end

% toc

return

function Plot_WormPixels_on_Frame1(CurrentTrack, Frame_Num, PixelsField, COLOR, ZoomOnWorm )

FrameIndex = find( CurrentTrack.Frames == Frame_Num, 1, 'first'); 
if ZoomOnWorm 
    MARKERSIZE  = 1;        
else   
    MARKERSIZE = 0.5;
end

if ~exist('COLOR','var')
    COLOR   = 'r';
end
LINESTYLE   = 'none'; % for PixelIdxList and HeadTail_PixelIdxList
MARKER      = '.';    % for PixelIdxList

if strcmpi(PixelsField,'Midline') 
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        hold on;
        plot(CurrentTrack.Midline.Y_coordinates_short{FrameIndex},CurrentTrack.Midline.X_coordinates_short{FrameIndex},'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);   
    end
elseif strcmpi(PixelsField,'Perimeter')
    MARKERSIZE = 2;
%     LINESTYLE   = '-';    
%     MARKER      = 'none'; 
    hold on;
    plot(single(CurrentTrack.WormPerimeter.Ycoordinate{FrameIndex}),single(CurrentTrack.WormPerimeter.Xcoordinate{FrameIndex}),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);    
elseif strcmpi(PixelsField,'HeadTail')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = '*';  
        if ZoomOnWorm 
            MARKERSIZE  = 8;
        else            
            MARKERSIZE  = 4; 
        end
        hold on;
        HeadTail_Matrix = squeeze(CurrentTrack.HeadTail(FrameIndex,:,:));
        plot(HeadTail_Matrix(:,2),HeadTail_Matrix(:,1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);           
    end
elseif strcmpi(PixelsField,'Head')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = '*';  
        if ZoomOnWorm 
            MARKERSIZE  = 8;
        else
            MARKERSIZE  = 4;
        end
        hold on;
        Head_Vec = CurrentTrack.Head(FrameIndex,:);
        plot(Head_Vec(2),Head_Vec(1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);         
    end
elseif strcmpi(PixelsField,'Tail')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = 'o';  
        if ZoomOnWorm 
            MARKERSIZE  = 8;
        else   
            MARKERSIZE  = 4; 
        end
        hold on;
        Tail_Vec = CurrentTrack.Tail(FrameIndex,:);
        plot(Tail_Vec(2),Tail_Vec(1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);     
    end
elseif strcmpi(PixelsField,'NeuronCoordinates')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = 'o';  
        if ZoomOnWorm 
            MARKERSIZE  = 8;
        else   
            MARKERSIZE  = 4; 
        end
        hold on;
        NeuronCoordinates_Vec = CurrentTrack.NeuronCoordinates(FrameIndex,:);
        plot(NeuronCoordinates_Vec(2),NeuronCoordinates_Vec(1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);     
    end   
    
end

return

function Tracks = CorrectMidlineForHeadDirection_AndCalculateProperties(Tracks, MidlineCalculationParams)  
plotme = false;

%% first correct all 'Midline' fields for head-tail switching:
NumOfTracks = length(Tracks);
Allfields   = {'X_coordinates_short','Y_coordinates_short','LinearIndices','Pattern','X_coordinates','Y_coordinates','Angle'};
NumOfFields = length(Allfields);
for tr = 1:NumOfTracks
    disp([datestr(now), ' -- Correcting midline fields to fit head-tail switching- Track ',num2str(tr)]);
    HeadDirectionWasSwitched = Tracks(tr).HeadDirectionWasSwitched;  
    Midline                  = Tracks(tr).Midline;
    RelevantSwitchIndices    = find(HeadDirectionWasSwitched);
    for ind = RelevantSwitchIndices
        for field_ind = 1:NumOfFields
            Midline.(Allfields{field_ind}){ind} = Midline.(Allfields{field_ind}){ind}(end:(-1):1);
        end        
    end  
    
    disp([datestr(now), ' -- Adding additional midline properties fields - Track ',num2str(tr)]);    
    InterpolationFactorBodyAngle = MidlineCalculationParams.InterpolationFactorBodyAngle  ;  % for csaps spline function of midline angles. More smoothed!
    long_vec                     = 1 : (MidlineCalculationParams.WormMatrixLength-1)      ;
    NumOfFrames                  = length(Tracks(tr).Frames);
    Midline.Angle                = cell(1,NumOfFrames);
    Midline.AngleFirstPoint      = single(zeros(1,NumOfFrames)*NaN);
    Midline.AngleLastPoint       = single(zeros(1,NumOfFrames)*NaN);
    Midline.NumOfBends_HighRes   = uint8(zeros(1,NumOfFrames));
    Midline.NumOfBends_LowRes    = uint8(zeros(1,NumOfFrames));

    for frame_ind=1:NumOfFrames
        if Midline.FlagTrueIfReliable(frame_ind)
            X = Midline.X_coordinates{frame_ind};
            Y = Midline.Y_coordinates{frame_ind};    

            %% Calculating the angle of the normals to the worms skeleton (using increased sampling rate)  
            THETA       = cart2pol(diff(X),diff(Y));   % THETA is a counterclockwise angular displacement in radians from the positive x-axis
            % FIXING SINGULAR POINTS
            THETA_PEAKS = find(abs(diff(THETA))> pi); % indices BEFORE the 2pi 'jump'
            for p_ind = 1:length(THETA_PEAKS)
                StartIndex = THETA_PEAKS(p_ind);    
                if THETA(StartIndex) < THETA(StartIndex+1) % Expect a sudden increase of 2*pi 
                    THETA((StartIndex+1):end) =  THETA((StartIndex+1):end) - 2*pi;
                else                                       % Expect a sudden deccrease of 2*pi 
                    THETA((StartIndex+1):end) =  THETA((StartIndex+1):end) + 2*pi;
                end
            end
            THETA_PEAKS_OK = isempty(find(abs(diff(THETA))> pi , 1, 'first')); % indices BEFORE the 2pi 'jump'
            if ~ THETA_PEAKS_OK
                disp('WARNING: singularities were found in angles calculation');
                Midline.FlagTrueIfReliable(frame_ind)=false;
                Flag=false
            end

            % ADDITIONAL SMOOTHING OF THE ANGLE VECTOR        
            pp                          = csaps(long_vec', double(THETA)', InterpolationFactorBodyAngle); 
            THETA_smoothed              = single(fnval(pp,long_vec));                
            pp                          = csaps(long_vec', double(THETA)', InterpolationFactorBodyAngle/100);                                              
            THETA_ForBendsCount_HighRes = single(fnval(pp,long_vec));        
            pp                          = csaps(long_vec', double(THETA)', InterpolationFactorBodyAngle/10000);                                              
            THETA_ForBendsCount_LowRes  = single(fnval(pp,long_vec));

    %         figure; plot(long_vec,THETA,'b.'); hold on; plot(long_vec,THETA_smoothed,'r:'); plot(long_vec,THETA_ForBendsCount_HighRes,'r-'); plot(long_vec,THETA_ForBendsCount_LowRes,'g-'); 

            %% Angle Output
            % The angles in 'Midline.Angle' are defined as the COUNTERCLOCKWISE DISPLACEMENT FROM THE +X AXIS. i.e. [0 90 180 270] corresponds to [+X +Y -X -Y], respectively.  
            THETA_ang                     = THETA_smoothed/pi*180;
            THETA_ang_ForBendCount_HighRes= THETA_ForBendsCount_HighRes/pi*180;
            THETA_ang_ForBendCount_LowRes = THETA_ForBendsCount_LowRes/pi*180;
            Midline.Angle{frame_ind}             = single(THETA_ang);   
            Midline.AngleFirstPoint(frame_ind)   = THETA_ang(1);          % The head direction is: 360- Midline.AngleFirstPoint
            Midline.AngleLastPoint(frame_ind)    = THETA_ang(end);        % The tail direction is:      Midline.AngleLastPoint

            %% bending output: 
            %   one bend == the curvature of the body is always to the same direction    
            %   2 bends  == the curvature of the body switches direction 1 time   
            %   3 bends  == the curvature of the body switches direction 2 times   
            %   4 or more bends are NOT ALLOWED and considered as errors. Midline.FlagTrueIfReliable will be corrected appropriately.     

            % HIGH RESOLUTION
            THETA_RelativeToHead                = THETA_ang_ForBendCount_HighRes - THETA_ang_ForBendCount_HighRes(1);       % This is the angle relative to the head orientation;
            CurvatureDirection                  = sign(diff(THETA_RelativeToHead));                         % (1) for counterclockwise, (-1) for clockwise.
            CurvatureDirectionSwitchingPoints   = find(diff(CurvatureDirection)~=0);
            NumOfBends                          = 1+length(CurvatureDirectionSwitchingPoints);      
            Midline.NumOfBends_HighRes(frame_ind)      = NumOfBends;    

            % LOW RESOLUTION
            THETA_RelativeToHead                = THETA_ang_ForBendCount_LowRes - THETA_ang_ForBendCount_LowRes(1);       % This is the angle relative to the head orientation;
            CurvatureDirection                  = sign(diff(THETA_RelativeToHead));                         % (1) for counterclockwise, (-1) for clockwise.
            CurvatureDirectionSwitchingPoints   = find(diff(CurvatureDirection)~=0);
            NumOfBends                          = 1+length(CurvatureDirectionSwitchingPoints);      
            Midline.NumOfBends_LowRes(frame_ind) = NumOfBends;   

            if plotme
                figure; plot(long_vec,THETA,'b.'); hold on; plot(long_vec,THETA_smoothed,'r:'); plot(long_vec,THETA_ForBendsCount_HighRes,'r-'); plot(long_vec,THETA_ForBendsCount_LowRes,'g-'); 
                figure; 
                plot(X,Y,'r'); hold on; plot(X(1),Y(1),'r*')
                title(['(',num2str(Midline.NumOfBends_LowRes(frame_ind)),',',num2str(Midline.NumOfBends_HighRes(frame_ind)), ') bends at (low/high) resolution.  First Angle= ',num2str(mod(360+Midline.AngleFirstPoint(frame_ind),360)),'.  Last Angle= ',num2str(mod(360+Midline.AngleLastPoint(frame_ind),360))]);
                pause;  
                close; close;
            end
        end                                                      
    end

    for frame_ind=find(Midline.FlagTrueIfReliable==0)
        Midline.Angle{frame_ind}               = [];   
    end    
    Tracks(tr).Midline = Midline;
end

return

function PatternMatrix = CorrectPatternMatrix_UsingHeadTailCoordinates (Tracks, PatternMatrix)
NumOfTracks = length(Tracks);
for tr=1:NumOfTracks
    HeadDirectionWasSwitched                        = Tracks(tr).HeadDirectionWasSwitched;
    FramesWithEstimatedHeadSide                              = Tracks(tr).FramesWithEstimatedHeadSide;
    PatternMatrix{tr}(HeadDirectionWasSwitched,:,:) = PatternMatrix{tr}(HeadDirectionWasSwitched,:,end:(-1):1);
    PatternMatrix{tr}(~FramesWithEstimatedHeadSide,:,:)      = NaN;
end
return

function LinkedMatrices = LinkPatternMatrix (File, VariableName, GoodLinks_ALL, Non_Deleted_Tracks_BeforeStitching, Non_Deleted_Tracks_BeforeManualInspection, TracksCorrectionParameters, NumberOfTracksPerFragment)
%% Update pattern matrix with TracksLinking. Do not store information of unreliable midlines 

LinkedMatrices = []; 
NumFragments   = File.Fragments;
NumberOfMatrixRowsPerFragment = zeros(1,NumFragments);
for Fragment = 1:NumFragments
    FileName = [File.FragmentSaveNames{Fragment}(1:end-4),'.mat'];
    disp([datestr(now),' -- Loading Track File ',int2str(Fragment),'... ']);
    loadsuccess = 0;
    while ~loadsuccess 
        try
            load(FileName,VariableName); 
            eval(['NumberOfMatrixRowsPerFragment(Fragment) = length(',VariableName,');']);
%             loadsuccess = 1;
        catch
            disp([VariableName, 'was not found for Fragment ',int2str(Fragment),'.  Skipping...']);
            eval([VariableName,' = [];']);
%             pause(20);
        end
        loadsuccess = 1;
    end                     
    if Fragment==1
        eval(['LinkedMatrices = ',VariableName,';']);
    else
        eval(['LinkedMatrices  = [LinkedMatrices, ',VariableName,'];']);
    end
    eval(['clear ',VariableName,';']);
end
disp([datestr(now),' -- Stitching ', VariableName]);

MaxNumOfWorms  = File.VariablesInformation.MaxNumOfWorms;
LinkedMatrices = StitchMatrices_BasedOnGoodLinks(LinkedMatrices, GoodLinks_ALL, Non_Deleted_Tracks_BeforeStitching, Non_Deleted_Tracks_BeforeManualInspection, TracksCorrectionParameters, MaxNumOfWorms);


return

%% Extract head/tail position and estimated neuron position
function Tracks = EstimatePositions_Head_Tail_Neuron_inline(Tracks, PatternMatrix)
% NOTE:   I do NOT switch midline direction and do NOT switch patterns direction here
% FREE PARAMETERS
MaximalNeuronWidth      = 15;   
MinimalNeuronProminence = 0.6;  

plotme      = false;
NumOfTracks = length(Tracks);

for tr=1:NumOfTracks
    disp([datestr(now), ' -- Estimating head, tail and neuron positions for track ',num2str(tr)]);
    CurrentPatterns  = PatternMatrix{tr};
    NumOfFrames      = length(Tracks(tr).Frames);
    FramesIndices    = find(Tracks(tr).Midline.FlagTrueIfReliable);
    Tracks(tr).FramesWithReliableHeadSide = false(1,NumOfFrames);
    Tracks(tr).HeadDirectionWasSwitched   = false(1,NumOfFrames);
    Tracks(tr).NeuronInWormCoordinates    = zeros(NumOfFrames,2,'single')*NaN;
    Tracks(tr).Head                       = zeros(NumOfFrames,2,'single')*NaN;   
    Tracks(tr).Tail                       = zeros(NumOfFrames,2,'single')*NaN;   
    for ind = FramesIndices            
        MAT = squeeze(CurrentPatterns(ind,:,:));
        [success, CurrentHeadInStart, CurrentNeuronInWormCoordinates] = CheckPatternPeaks(MAT,MaximalNeuronWidth, MinimalNeuronProminence, plotme);
        if success
            Tracks(tr).FramesWithReliableHeadSide(ind) = true;
            Tracks(tr).NeuronInWormCoordinates(ind,:) = CurrentNeuronInWormCoordinates;
            if CurrentHeadInStart
                Tracks(tr).Head(ind,:) = squeeze(Tracks(tr).HeadTail(ind,1,:));
                Tracks(tr).Tail(ind,:) = squeeze(Tracks(tr).HeadTail(ind,2,:));
            else
                Tracks(tr).HeadDirectionWasSwitched(ind)= true;
                Tracks(tr).Head(ind,:) = squeeze(Tracks(tr).HeadTail(ind,2,:));
                Tracks(tr).Tail(ind,:) = squeeze(Tracks(tr).HeadTail(ind,1,:));
            end
        end
    end    
end

% ShowEstimatedPosition(PatternMatrix, ReliableFrameIndexPerTrack, GoodFrameExists, HeadInStart, NeuronCoordinates);

return

function [success, HeadInStart, NeuronCoordinates] = CheckPatternPeaks(MAT, MaximalNeuronWidth, MinimalNeuronProminence, plotme) 

success           = false;
HeadInStart       = false;
NeuronCoordinates = [NaN NaN];

LeftSideEnd  = size(MAT,2)/5;
RightSideEnd = size(MAT,2)*4/5;

Vec                     = mean(MAT,1);
[pks,locs,widths,proms] = findpeaks(double(Vec));
if length(locs)>1
    proms((locs>LeftSideEnd)& (locs<RightSideEnd))=0;
else
    return
end

[~,sortorder]          = sort(proms);
LargestPeaksIndices    = sortorder([end,(end-1)]);
LargestPeaksLocations  = locs(LargestPeaksIndices);   

HeadSideIsDefined = false;
if all(LargestPeaksLocations < LeftSideEnd)   % Both peaks are on left side >> high probability that head side is left
    if LargestPeaksLocations(2)<LargestPeaksLocations(1) % Peaks need to be switched since neuron comes before autofluorescence
        LargestPeaksIndices = LargestPeaksIndices([2 1]);
        LargestPeaksLocations  = locs(LargestPeaksIndices);
    end
    HeadSideIsDefined = true;
    HeadInStart       = true;
elseif all(LargestPeaksLocations > RightSideEnd)   % Both peaks are on right side >> high probability that head side is right
    if LargestPeaksLocations(1)<LargestPeaksLocations(2) % Peaks need to be switched since neuron comes before autofluorescence
        LargestPeaksIndices = LargestPeaksIndices([2 1]);
        LargestPeaksLocations  = locs(LargestPeaksIndices);
    end
    HeadSideIsDefined = true;
    HeadInStart       = false; 
end
        
LargestPeaksValues     = pks(LargestPeaksIndices);
LargestPeaksProminence = proms(LargestPeaksIndices);
LargestPeaksWidth      = widths(LargestPeaksIndices);

ProminenceRatio      = LargestPeaksProminence(1)/LargestPeaksProminence(2);
DynamicRange         = max(Vec)-min(Vec); 
NormalizedProminence = LargestPeaksProminence / DynamicRange;

if HeadSideIsDefined && (LargestPeaksWidth(1)<MaximalNeuronWidth) && (NormalizedProminence(1) > MinimalNeuronProminence)
    success = true;
end
if success
    [~, NeuronWidthCoordinate] = max(MAT(:,LargestPeaksLocations(1)));
    NeuronCoordinates = [LargestPeaksLocations(1) NeuronWidthCoordinate];
end    
    
% if plotme     
if plotme && success     
    figure('position',get(0,'ScreenSize'))
    subplot(2,1,1);    
        imshow(MAT,[]);        
        if success
            if HeadInStart
                TitleStr = 'Head In Start';
            else
                TitleStr = 'Head In End';
            end   
            hold on; plot(NeuronCoordinates(1),NeuronCoordinates(2),'ro');
        else
            TitleStr = 'NO success';
        end            
        title(TitleStr);
            
    subplot(2,1,2);         
        plot(Vec);   hold on; 
        plot(LargestPeaksLocations(1),LargestPeaksValues(1),'r*')
        plot(LargestPeaksLocations(2),LargestPeaksValues(2),'b*')
        xlim([0 length(Vec)])
        title(['Ratio= ',num2str(ProminenceRatio,2),...
            ', W1= ',num2str(LargestPeaksWidth(1),2),     ', W2= ',num2str(LargestPeaksWidth(2),2),...
            ', Prom= ',num2str(LargestPeaksProminence(1)),', Prom/Range=',num2str(NormalizedProminence(1),2)])  
end   

return

function ShowEstimatedPosition(PatternMatrix, ReliableFrameIndexPerTrack, GoodFrameExists, HeadInStart, NeuronCoordinates)
cols = 3;
rows = 8;
MaxPlotsPerFigure = cols*rows;
NumOfTracks                = length(PatternMatrix);

subplot_ind = 0;
figure('position',get(0,'ScreenSize'))
for tr=1:NumOfTracks
    subplot_ind = subplot_ind+1;
    if subplot_ind > MaxPlotsPerFigure
        subplot_ind = 1;
        figure('position',get(0,'ScreenSize'))
    end
    MAT     = squeeze(PatternMatrix{tr}(ReliableFrameIndexPerTrack(tr),:,:));
    if GoodFrameExists(tr)
        CurrentNeuronCoordinates = NeuronCoordinates(tr,:);
        if ~HeadInStart(tr)           
            MAT = MAT(:,end:-1:1);
            CurrentNeuronCoordinates(1) = size(MAT,2)-CurrentNeuronCoordinates(1);
        end
    else
        MAT = zeros(size(MAT));
    end
    subplot(rows,cols,subplot_ind);    
    imshow(MAT,[]);           
    hold on; plot(CurrentNeuronCoordinates(1),CurrentNeuronCoordinates(2),'ro');
    title(tr)
end 

return

%% If needed, functions for manual correction of head/tail position
function Tracks = ManualCorrectionOfHeadTailPositions_inline(File,Tracks)
% parameters
ThresholdForManualInspection                = 0.7;     % minimal fraction of "hard segmentation" area needed for checking the track manually 
ThresholdForManualInspection2               = 0.3; 
MaximalFrameDistanceForAutoHeadEstimation   = 1/3* File.FrameRate; % not more than 1/3 second  
MaximalDistanceForAutoHeadEstimation        = File.PixelSize / 5;  % not far from 200 micrometer  
Tracks_BeforeManualCorrection               = Tracks;

% Filename:
FileNameAfterInspection       = [File.TrackFile(1:end-4),'_TracksAfterManualInspection.mat'];

if exist(FileNameAfterInspection,'file')
    qstring = '''_TracksAfterManualInspection.mat'' file was found. Load existing file or re-inspect Tracks file?';
    default = 'Load existing file';
    button = questdlg(qstring,'Manual inspection of head/tail segmentation','Load existing file','re-inspect Tracks',default); 
    if strcmpi(button,'Load existing file')
        disp('loading existing manually inspected file');
        load(FileNameAfterInspection,'Tracks');
        return        
    end
end
% Manual Inspection based on head segmentation quality 
NumOfTracks  = length(Tracks);
MidlineReliabilityMatrix = zeros(NumOfTracks,File.NumberOfFrames,'single');
HeadReliabilityMatrix    = zeros(NumOfTracks,File.NumberOfFrames,'single');
for tr = 1:NumOfTracks
    Frames                     = Tracks(tr).Frames;
    FlagTrueIfReliableHeadSide = Tracks(tr).FramesWithReliableHeadSide;
    FlagTrueIfReliableMidline  = Tracks(tr).Midline.FlagTrueIfReliable;
    FramesWithReliableHeadSide = Frames(FlagTrueIfReliableHeadSide);
    FramesWithReliableMidline  = Frames(FlagTrueIfReliableMidline);
    HeadReliabilityMatrix(tr,FramesWithReliableHeadSide)   = 1;
    MidlineReliabilityMatrix(tr,FramesWithReliableMidline) = 1;       
end
HardSegmentationMatrix = (HeadReliabilityMatrix==0)&(MidlineReliabilityMatrix==1);
windowSize = File.FrameRate*2;
b          = (1/windowSize)*ones(1,windowSize);
a          = 1;
delay      = round(windowSize/2);
HardSegmentationMatrix_Filtered = filter(b,a,single(HardSegmentationMatrix),[],2);
HardSegmentationMatrix_Filtered = [HardSegmentationMatrix_Filtered(:,delay:end), ones(NumOfTracks,delay-1)*NaN];
HeadReliabilityMatrix_Filtered  = filter(b,a,HeadReliabilityMatrix,[],2);
HeadReliabilityMatrix_Filtered  = [HeadReliabilityMatrix_Filtered(:,delay:end), ones(NumOfTracks,delay-1)*NaN];

% FramesToInspect = (HardSegmentationMatrix_Filtered > ThresholdForManualInspection);
FramesToInspect = (HardSegmentationMatrix_Filtered > ThresholdForManualInspection)|(HeadReliabilityMatrix_Filtered<ThresholdForManualInspection2);

% reject head segmentation in frames for which the head was frequently not properly segmented within a close time interval   
for tr = 1:NumOfTracks
    Frames                                = Tracks(tr).Frames;
    Additional_Unreliable_Frames          = FramesToInspect(tr,Frames);
    Tracks(tr).FramesWithReliableHeadSide(Additional_Unreliable_Frames) = false;   
    Tracks(tr).HeadDirectionWasSwitched(Additional_Unreliable_Frames)   = false;   
    Tracks(tr).Head(Additional_Unreliable_Frames,:)                     = NaN;   
    Tracks(tr).Tail(Additional_Unreliable_Frames,:)                     = NaN;   
    Tracks(tr).NeuronInWormCoordinates(Additional_Unreliable_Frames,:)  = NaN;   
end
for tr = 1:NumOfTracks
    Frames                        = Tracks(tr).Frames;
    FlagTrueIfReliableHeadSide    = Tracks(tr).FramesWithReliableHeadSide;
    FramesWithNONReliableHeadSide = Frames(~FlagTrueIfReliableHeadSide);
    HeadReliabilityMatrix(tr,FramesWithNONReliableHeadSide)   = 0;
end    

%% Estimate head/tail position based on location and distance from head at reliable frames         
% Loop over each track and find frames without head/tail estimation
% try to estimate head/tail at frames with proper midline detected based on distance from coordinates at reliable frames
% Check consistency (between previous and following reliable frames). If consistent, except estimation, store flag-decision in vector 
% Find frames that are still without head/tail estimation. Manually ask for estimated head position

Tracks = EstimateHeadTailByPositionUsingReliableFrames (Tracks, MaximalFrameDistanceForAutoHeadEstimation, MaximalDistanceForAutoHeadEstimation);

disp('iteration 1 of 4');
MaximalAllowedFrameIntervalWithoutHeadSegmentation = 1000;
for tr = 1:NumOfTracks
    Tracks = ManualInspectHeadTail_DenseSampling (File, Tracks, tr, MaximalAllowedFrameIntervalWithoutHeadSegmentation, Tracks_BeforeManualCorrection);
end
Tracks = EstimateHeadTailByPositionUsingReliableFrames (Tracks, MaximalFrameDistanceForAutoHeadEstimation, MaximalDistanceForAutoHeadEstimation);

disp('iteration 2 of 4');
MaximalAllowedFrameIntervalWithoutHeadSegmentation = 300;
for tr = 1:NumOfTracks
    Tracks = ManualInspectHeadTail_DenseSampling (File, Tracks, tr, MaximalAllowedFrameIntervalWithoutHeadSegmentation, Tracks_BeforeManualCorrection);
end
Tracks = EstimateHeadTailByPositionUsingReliableFrames (Tracks, MaximalFrameDistanceForAutoHeadEstimation, MaximalDistanceForAutoHeadEstimation);

disp('iteration 3 of 4');
MaximalAllowedFrameIntervalWithoutHeadSegmentation = 90;
for tr = 1:NumOfTracks
    Tracks = ManualInspectHeadTail_DenseSampling (File, Tracks, tr, MaximalAllowedFrameIntervalWithoutHeadSegmentation, Tracks_BeforeManualCorrection);
end
Tracks = EstimateHeadTailByPositionUsingReliableFrames (Tracks, MaximalFrameDistanceForAutoHeadEstimation, MaximalDistanceForAutoHeadEstimation);

disp('iteration 4 of 4');
MaximalAllowedFrameIntervalWithoutHeadSegmentation = 30;
for tr = 1:NumOfTracks
    Tracks = ManualInspectHeadTail_DenseSampling (File, Tracks, tr, MaximalAllowedFrameIntervalWithoutHeadSegmentation, Tracks_BeforeManualCorrection);
end
Tracks = EstimateHeadTailByPositionUsingReliableFrames (Tracks, MaximalFrameDistanceForAutoHeadEstimation, MaximalDistanceForAutoHeadEstimation);

%% Save Manually corrected tracks
save(FileNameAfterInspection,'Tracks','File','-v7.3')

% save('L:\ImagingData\20160211_2\20160211_2_TracksAfterManualInspection2.mat','Tracks','-v7.3')

return

function Tracks = ManualInspectHeadTail_DenseSampling (File, Tracks, tr, MaximalAllowedFrameIntervalWithoutHeadSegmentation, Tracks_BeforeManualCorrection)

AllFrames                   = Tracks(tr).Frames;
FlagTrueIfEstimatedHeadSide = Tracks(tr).FramesWithEstimatedHeadSide;
FlagTrueIfReliableMidline   = Tracks(tr).Midline.FlagTrueIfReliable;
FramesWithEstimatedHeadSide = AllFrames(FlagTrueIfEstimatedHeadSide);

FlagTrueIfReliableMidlineWithoutHeadSegmentation = FlagTrueIfReliableMidline & (~FlagTrueIfEstimatedHeadSide);


IndicesWhereNoHeadSegmentationStarts = find(diff(FlagTrueIfReliableMidlineWithoutHeadSegmentation)==1)+1;
FramesWhereNoHeadSegmentationStarts  = AllFrames(IndicesWhereNoHeadSegmentationStarts);

IndicesToCheckManually = [];
for ind = 1:length(IndicesWhereNoHeadSegmentationStarts)
    CurrentIndexWithNoHeadSegmentation = IndicesWhereNoHeadSegmentationStarts(ind);
    CurrentFrameWithNoHeadSegmentation = FramesWhereNoHeadSegmentationStarts(ind);
    NextFrameWithHeadSegmentation      = FramesWithEstimatedHeadSide(find(FramesWithEstimatedHeadSide>CurrentFrameWithNoHeadSegmentation,1,'first'));
    if isempty(NextFrameWithHeadSegmentation)
        NextFrameWithHeadSegmentation = File.NumberOfFrames+1;
    end
    FramesUntilNextHeadSegmentedFrame  = NextFrameWithHeadSegmentation - CurrentFrameWithNoHeadSegmentation;
    if FramesUntilNextHeadSegmentedFrame > MaximalAllowedFrameIntervalWithoutHeadSegmentation  
        IndicesToCheckManually = [IndicesToCheckManually CurrentIndexWithNoHeadSegmentation];                      
    end
end
done = isempty(IndicesToCheckManually);
IndicesOfFramesWithInitialGuessAvailable  = find(Tracks_BeforeManualCorrection(tr).FramesWithReliableHeadSide);

uicontrol_xy               = [0.06 0.025]; 
if ~done
    IndicesToCheckManually = ...
        [IndicesToCheckManually(diff(AllFrames(IndicesToCheckManually)) > MaximalAllowedFrameIntervalWithoutHeadSegmentation/2), ....
         IndicesToCheckManually(end) ];
     IndicesToCheckManually = unique(IndicesToCheckManually);
end

while ~done        
    NumOfFramesToShow            = min([8, length(IndicesToCheckManually)]);       
    IndicesOfFramesToShow        = IndicesToCheckManually(1:NumOfFramesToShow);  
    FramesToShow                 = AllFrames(IndicesOfFramesToShow);
    
    Tracks(tr).Head(IndicesOfFramesToShow,:) = squeeze(Tracks(tr).HeadTail(IndicesOfFramesToShow,1,:)); % initial random guess
    Tracks(tr).Tail(IndicesOfFramesToShow,:) = squeeze(Tracks(tr).HeadTail(IndicesOfFramesToShow,2,:));
    
    CurrentIndicesWithGuessAvailable = intersect(IndicesOfFramesWithInitialGuessAvailable,IndicesOfFramesToShow);
    if ~isempty(CurrentIndicesWithGuessAvailable)
        Tracks(tr).Head(CurrentIndicesWithGuessAvailable,:) = Tracks_BeforeManualCorrection(tr).Head(CurrentIndicesWithGuessAvailable,:); 
        Tracks(tr).Tail(CurrentIndicesWithGuessAvailable,:) = Tracks_BeforeManualCorrection(tr).Tail(CurrentIndicesWithGuessAvailable,:);
        Tracks(tr).HeadDirectionWasSwitched(CurrentIndicesWithGuessAvailable) = Tracks_BeforeManualCorrection(tr).HeadDirectionWasSwitched(CurrentIndicesWithGuessAvailable);
    end
    
    Head = Tracks(tr).Head;
    Tail = Tracks(tr).Tail;
    
    Fig_handle=figure('position',get(0,'screensize'),'name',...
              ['Track ',num2str(tr),'. RIGHT click if OK (* designates Head). LEFT click on subplots with WRONG head-tail segmentation']);  

    for sp_ind = 1:NumOfFramesToShow
        subplot(2,4,sp_ind)
        CurrentFrame = FramesToShow(sp_ind);
        WormTrackingImage(File, CurrentFrame, Tracks(tr)); title(CurrentFrame);
        P       = get(gca,'position');
        PstartX = P(1)+P(3)/2-uicontrol_xy(1)/2;
        PstartY = P(2)-uicontrol_xy(2);
        if IndicesOfFramesToShow(sp_ind) < 2*File.FrameRate+1;
            MovieFrames = CurrentFrame:min([CurrentFrame+2*File.FrameRate max(AllFrames)]);            
        else
            MovieFrames  = (CurrentFrame-2*File.FrameRate):CurrentFrame;
        end
        ax=gca;
        buttons(sp_ind)  = uicontrol(Fig_handle,'units','normalized','BackgroundColor','r','position',[PstartX  PstartY uicontrol_xy],'Style','popupmenu','String',{'Correct','Incorrect','not clear'});
        buttons2(sp_ind) = uicontrol(Fig_handle,'units','normalized','BackgroundColor','r','position',[PstartX  PstartY-uicontrol_xy(2) uicontrol_xy],'Style','pushbutton','String','1 second movie','Callback', {@WormTrackingMovie2,File, MovieFrames, Tracks(tr), ax});
    end
    button_finish = uicontrol(Fig_handle,'units','normalized','BackgroundColor','g','position',[0.0370    0.5040    0.0600    0.0250],'Style','popupmenu','String',{'Review in process','Finish review'});
    % Add here user input
    % if corrections are needed, use only corrected frames for further Head/Tail estimations  
    success = button_finish.Value==2;
    while ~success
        k = waitforbuttonpress;        
        success = button_finish.Value==2;
    end
    Correct   = [buttons.Value]==1; % correct segmentation
    Incorrect = [buttons.Value]==2; % need to be switched
    NotClear  = [buttons.Value]==3; % not clear which side is the head
    
    CorrectIndices             = IndicesOfFramesToShow(Correct);
    IncorrectIndices           = IndicesOfFramesToShow(Incorrect);
    NotClearIndices            = IndicesOfFramesToShow(NotClear);        
    
    if ~isempty(IncorrectIndices)
        Tracks(tr).HeadDirectionWasSwitched(IncorrectIndices) = ~Tracks(tr).HeadDirectionWasSwitched(IncorrectIndices);
        Tracks(tr).Head(IncorrectIndices,:)                   = Tail(IncorrectIndices,:);
        Tracks(tr).Tail(IncorrectIndices,:)                   = Head(IncorrectIndices,:);
    end
    if ~isempty(NotClearIndices)
        Tracks(tr).Head(NotClearIndices,:)                    = zeros(length(NotClearIndices),2)*NaN;
        Tracks(tr).Tail(NotClearIndices,:)                    = zeros(length(NotClearIndices),2)*NaN;
        Tracks(tr).HeadDirectionWasSwitched(NotClearIndices)  = false;
    end
    
    Tracks(tr).FramesWithReliableHeadSide([CorrectIndices IncorrectIndices])  = true;
    Tracks(tr).FramesWithEstimatedHeadSide([CorrectIndices IncorrectIndices]) = true;
    IndicesToCheckManually = setdiff(IndicesToCheckManually,IndicesOfFramesToShow);
    
    done = isempty(IndicesToCheckManually);
    close; 
    clear buttons buttons2
    
end


return

function Tracks = EstimateHeadTailByPositionUsingReliableFrames (Tracks, MaximalFrameDistanceForAutoHeadEstimation, MaximalDistanceForAutoHeadEstimation);
NumOfTracks = length(Tracks);
for tr=1:NumOfTracks
    FramesWithReliableHeadSide  = Tracks(tr).FramesWithReliableHeadSide;
    HeadDirectionWasSwitched    = Tracks(tr).HeadDirectionWasSwitched;
    FlagTrueIfReliableMidline   = Tracks(tr).Midline.FlagTrueIfReliable;
    FramesNumberInMovie         = Tracks(tr).Frames;  
    if isfield(Tracks(tr),'FramesWithEstimatedHeadSide') && ~isempty(Tracks(tr).FramesWithEstimatedHeadSide)
        FramesWithEstimatedHeadSide  = Tracks(tr).FramesWithEstimatedHeadSide;
    else
        FramesWithEstimatedHeadSide = FramesWithReliableHeadSide;   % initialization
    end
    HeadTail_Coordinates        = Tracks(tr).HeadTail;        
    Head_Coordinates            = Tracks(tr).Head;        
    Tail_Coordinates            = Tracks(tr).Tail;    
    MinimalDistance             = zeros(1,length(FramesWithReliableHeadSide),'single')*NaN;    
    DistanceDifferences         = zeros(1,length(FramesWithReliableHeadSide),'single')*NaN;    
    
    FramesToAutoAnalyzeBasedOnLocation = find((~FramesWithEstimatedHeadSide)&(FlagTrueIfReliableMidline));
    
    % Track preferentially from the beginning to the end of each track
    for ind = FramesToAutoAnalyzeBasedOnLocation        
        distance_to_closest_previous_estimation = find(FramesWithEstimatedHeadSide(ind:(-1):1),1,'first')-1;
        ind_closest_previous_estimation         = ind - distance_to_closest_previous_estimation;      
        FrameDistance                           = FramesNumberInMovie(ind) - FramesNumberInMovie(ind_closest_previous_estimation);
        if isempty(FrameDistance)
            continue
        elseif FrameDistance==1
            ind_closest = ind_closest_previous_estimation;
        else
            ind_closest = ind_closest_previous_estimation;
            distance_to_closest_next_estimation     = find(FramesWithEstimatedHeadSide(ind:end),1,'first')-1;
            ind_closest_next_estimation             = ind + distance_to_closest_next_estimation;    
            FrameDistance_next                      = FramesNumberInMovie(ind_closest_next_estimation) - FramesNumberInMovie(ind);
            if FrameDistance>FrameDistance_next  % compare to next frame rather than to previous one
                ind_closest   = ind_closest_next_estimation;
                FrameDistance = FrameDistance_next;
            end
        end
        Current_Head_Coordinates   = squeeze(HeadTail_Coordinates(ind,1,:))';
        Current_Tail_Coordinates   = squeeze(HeadTail_Coordinates(ind,2,:))';
        Reference_Head_Coordinates = Head_Coordinates(ind_closest,:);

        Distance.Head_ReferenceHead = sum((Current_Head_Coordinates-Reference_Head_Coordinates).^2);
        Distance.Tail_ReferenceHead = sum((Current_Tail_Coordinates-Reference_Head_Coordinates).^2);
        CurrentDistance             = min ([Distance.Head_ReferenceHead, Distance.Tail_ReferenceHead]);
        MinimalDistance(ind)        = CurrentDistance;
        DistanceDifferences(ind)    = abs (Distance.Head_ReferenceHead -Distance.Tail_ReferenceHead);
        DistanceMakeSense           = (CurrentDistance<=MaximalDistanceForAutoHeadEstimation) || ...
                                      ((CurrentDistance<=3*MaximalDistanceForAutoHeadEstimation)&&(DistanceDifferences(ind)>3*MaximalDistanceForAutoHeadEstimation));

        if (FrameDistance <= MaximalFrameDistanceForAutoHeadEstimation) && DistanceMakeSense
            FramesWithEstimatedHeadSide(ind) = true;
            
            if Distance.Head_ReferenceHead <= Distance.Tail_ReferenceHead
                Head_Coordinates(ind,:) = Current_Head_Coordinates;
                Tail_Coordinates(ind,:) = Current_Tail_Coordinates;
            else
                HeadDirectionWasSwitched(ind) = true;
                Head_Coordinates(ind,:) = Current_Tail_Coordinates;
                Tail_Coordinates(ind,:) = Current_Head_Coordinates;                     
                
            end
        end        
    end 
    % Now track from the end to the beginning of each track
    FramesToAutoAnalyzeBasedOnLocation = find((~FramesWithEstimatedHeadSide)&(FlagTrueIfReliableMidline));
    for ind = FramesToAutoAnalyzeBasedOnLocation(end:-1:1)        
        distance_to_closest_previous_estimation = find(FramesWithEstimatedHeadSide(ind:(-1):1),1,'first')-1;
        ind_closest_previous_estimation         = ind - distance_to_closest_previous_estimation;      
        FrameDistance                           = FramesNumberInMovie(ind) - FramesNumberInMovie(ind_closest_previous_estimation);
        if isempty(FrameDistance)
            continue
        elseif FrameDistance==1
            ind_closest = ind_closest_previous_estimation;
        else
            ind_closest = ind_closest_previous_estimation;
            distance_to_closest_next_estimation     = find(FramesWithEstimatedHeadSide(ind:end),1,'first')-1;
            ind_closest_next_estimation             = ind + distance_to_closest_next_estimation;    
            FrameDistance_next                      = FramesNumberInMovie(ind_closest_next_estimation) - FramesNumberInMovie(ind);
            if FrameDistance>FrameDistance_next  % compare to next frame rather than to previous one
                ind_closest   = ind_closest_next_estimation;
                FrameDistance = FrameDistance_next;
            end
        end
        Current_Head_Coordinates   = squeeze(HeadTail_Coordinates(ind,1,:))';
        Current_Tail_Coordinates   = squeeze(HeadTail_Coordinates(ind,2,:))';
        Reference_Head_Coordinates = Head_Coordinates(ind_closest,:);

        Distance.Head_ReferenceHead = sum((Current_Head_Coordinates-Reference_Head_Coordinates).^2);
        Distance.Tail_ReferenceHead = sum((Current_Tail_Coordinates-Reference_Head_Coordinates).^2);
        CurrentDistance             = min ([Distance.Head_ReferenceHead, Distance.Tail_ReferenceHead]);
        MinimalDistance(ind)        = CurrentDistance;
        DistanceDifferences(ind)    = abs (Distance.Head_ReferenceHead -Distance.Tail_ReferenceHead);

        DistanceMakeSense           = (CurrentDistance<=MaximalDistanceForAutoHeadEstimation) || ...
                                      ((CurrentDistance<=3*MaximalDistanceForAutoHeadEstimation)&&(DistanceDifferences(ind)>3*MaximalDistanceForAutoHeadEstimation));
                                  
        if (FrameDistance <= MaximalFrameDistanceForAutoHeadEstimation) && DistanceMakeSense
            FramesWithEstimatedHeadSide(ind) = true;
            
            if Distance.Head_ReferenceHead <= Distance.Tail_ReferenceHead
                Head_Coordinates(ind,:) = Current_Head_Coordinates;
                Tail_Coordinates(ind,:) = Current_Tail_Coordinates;
            else
                HeadDirectionWasSwitched(ind) = true;
                Head_Coordinates(ind,:) = Current_Tail_Coordinates;
                Tail_Coordinates(ind,:) = Current_Head_Coordinates;                     
                
            end
        end        
    end 
    
    
    Tracks(tr).Head                        = Head_Coordinates;    
    Tracks(tr).Tail                        = Tail_Coordinates;    
    Tracks(tr).HeadDirectionWasSwitched    = HeadDirectionWasSwitched;    
    Tracks(tr).FramesWithEstimatedHeadSide = FramesWithEstimatedHeadSide;                
%     figure; plot(Head_Coordinates(:,1),Head_Coordinates(:,2),'b*'); hold on; 
%             plot(Tracks(tr).Head(:,1),Tracks(tr).Head(:,2),'bo'); 
%             plot(Tail_Coordinates(:,1),Tail_Coordinates(:,2),'g*');
%             plot(Tracks(tr).Tail(:,1),Tracks(tr).Tail(:,2),'go');
end

return 

function WormTrackingImage(File, SpecificFrames, Tracks, ArenaID)

%% Parameters:
MovieOptions.ZoomOnWorm                   = true; %%% Assuming there is only One Worm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
MovieOptions.ShowCTXposition              = false;                     
MovieOptions.ShowBehaviorFrom_Data_Matrices = false;
MovieOptions.ShowTracksOnBackground       = false;     % Default=false; Tracks are shown by default on the original movie frames that are read from the path given in 'File'

% Other paremeters
MovieOptions.deltaT_Collision             = 0.05;      % Display time in regular frames
MovieOptions.deltaT                       = 0.05;      % Display time in regular frames
% MovieOptions.arena                        = 1;     
if exist('ArenaID','var')
    MovieOptions.arena                      = ArenaID;     
end
if exist('SpecificFrames','var')
    if ~isempty(SpecificFrames)
        MovieOptions.ShowSpecificFrames = SpecificFrames;
    end
end

MovieOptions.TracksLine                   = false;      % true/false/Num = show/(don't show)/(Show last Num of frames) of track line preceeding to each frame.                                                         
MovieOptions.TracksText.TrackNum          = true;      % Track text diplay settings
MovieOptions.TracksText.ID1               = false;     
MovieOptions.TracksText.ID2               = false;     
MovieOptions.TracksText.Behavior          = false; 
MovieOptions.TracksText.ColorByID         = 0;     % No, same colors for all worms. different colors based on features 

% MovieOptions.ShowDetectedPixels           = {'Area','Skeleton','HeadTail'}; % A cell with the name of ields to display
% MovieOptions.ShowDetectedPixelsColors     = {'b','r','g'};                  % Color for display
MovieOptions.ShowDetectedPixels           = {'Perimeter'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'b'};                  % Color for display
MovieOptions.ShowDetectedPixels           = {'Perimeter','Midline'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'b','g'};                  % Color for display
MovieOptions.ShowDetectedPixels           = {'Perimeter','Midline','Head','Tail'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'b','g','r','m'};                      % Color for display
MovieOptions.ShowDetectedPixels           = {'Midline','Head','Tail'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'g','r','m'};                      % Color for display

BehaviorStructure = [];
            
% MovieOptions.ConstantScale                = true;   % If true, the function will look for the best scaling and scale each frame similarly. If false --> imshow(Mov,[]) 
MovieOptions.ConstantScale                = 2;        % If true, the function will look for the best scaling and scale each frame similarly. If false --> imshow(Mov,[]) 
% MovieOptions.ConstantScale                = false;     
MovieOptions.CreateMovie                  = false;    % If true create a movie
Data_Matrices = [];
background = [];

SimpleDisplay_FluorescenceMovieWithTracks1(Tracks, File, MovieOptions, Data_Matrices, BehaviorStructure, background);  % Optional Movie generation

return

function SimpleDisplay_FluorescenceMovieWithTracks2 (Tracks, File, MovieOptions, Data_Matrices, BehaviorStructure, background)
% Display a movie with Tracks
% Information about the Tracks and the movie is given in the Tracks and File input variables  
% MovieOptions-  A stucture with infomation about what and how to display, with the following fields:   
%    deltaT:             Pause between frames in seconds (Default = 5msec).    
%    deltaT_collision:   Pause between frames near end of tracks in seconds (Default = 1 sec).    
%    arena:              Zoom on one specific arena (Default = no Zoom).       
%                        This field is NOT used if specific tracks are shown (arena zoom will be set by the tracks) 
%    EnlargeArena:       Enlarge Arena zoom by the specified number of pixels (Default = no enlargement).       
%    ShowSpecificTracks: A vector with Tracks of interest to display. Frames will be set by Tracks (Default = Show all tracks within the frames of interest).   
%                        Other Tracks within the same region of interest will be shown as well. 
%    ShowEndOfTracks:    A structure that, if exists, specifiy how many frames to display before the end of the Track (Default= show all the track)
%                        This is convenient to see worms before and after collisions  
%         FramesBefore   Number of frames before the first track ends (before first collision).  
%         FramesAfter    Number of frames after the start of the last detected track (after last collision).            
%    ShowSpecificFrames: A vector of Frames of interest to display. All relevant tracks will be shown(Default = All frames in the movie). 
%    TracksText:         Structure defining what fields to display. Values are true or false (Default = show only track number). 
%         TrackNum       if exists and true, show the track number (#)
%         ID1            Worm ID == original track calculated from linking analysis by Dirk (DL).  
%         ID2            Worm ID == original track calculated from my linking analysis (SL).  
%         Behavior      Worm behavior as calculated in the segmentation functions.   
% MovieOptions.TracksLine

warning off; 
MoviePath      = File.MoviePath;
MovieFileNames = File.MovieFileNames;

%% Movie format
if     strcmpi(File.VideoFormat,'tiff')
    VideoFormat   = 1;                                  % read new file every loop iteration 
elseif strcmpi(File.VideoFormat,'multi-tiff')
    VideoFormat   = 2;
    FileFullName  = [MoviePath,'\',MovieFileNames{1}];  % one file
elseif strcmpi(File.VideoFormat,'avi')
    VideoFormat   = 3;
    FileFullName  = [MoviePath,'\',MovieFileNames{1}];  % one file
    MovieObj      = VideoReader(FileFullName);          % one object
elseif strcmpi(File.VideoFormat,'multi-avi')
    VideoFormat   = 4;
    PreviousMovieObjFile = 0;                           % read new file object if and when it is needed 
end


%% Find the frames of each track
LEN             = length([Tracks.Active]);
FirstFrames     = zeros(1,LEN);
LastFrames      = zeros(1,LEN);
ALLFRAMES       = false(File.FragmentFrames(end,end),LEN);
for tr=1:LEN
    CurrentFrames               = Tracks(tr).Frames;
    FirstFrames(tr)             = CurrentFrames(1);
    LastFrames(tr)              = CurrentFrames(end);
    ALLFRAMES(CurrentFrames,tr) = 1;
end

%% Set parameters for the movie:   show specific tracks or specific frames,  show a single arena or all arenas,  frame rate
if ~exist('MovieOptions','var')
    deltaT             = 0.005;  
    arena              = NaN;     
    Frames2Show        = 1:size(ALLFRAMES,1);  
    ShowTracksText     = false;  
    ShowDetectedPixels = false;
    ShowDetectedROIs   = false;
    ConstantScale      = false;
    CreateMovie        = false;
    ColorByID          = false;
    ShowTracksOnBackground = false;
    Add_CTXdata_And_SmoothPosition    = false;
    ShowBehaviorFrom_Data_Matrices    = false;
    OnlyBehaviorColor_NoText          = false;
else
    
    ShowTracksText = false;
    ColorByID      = false;
    OnlyBehaviorColor_NoText = false;
    if isfield(MovieOptions,'TracksText')
        field_names    = fieldnames(MovieOptions.TracksText);
        field_names    = field_names(~strcmpi(field_names,'ColorByID'));
        for f_ind = 1:length(field_names)
            if islogical(MovieOptions.TracksText.(field_names{f_ind}))
                ShowTracksText = ShowTracksText || MovieOptions.TracksText.(field_names{f_ind});  % at least one text field = true
            else
                ShowTracksText = true;
            end
        end
        if isfield(MovieOptions.TracksText,'ColorByID')
            ColorByID = MovieOptions.TracksText.ColorByID;
        end 
       
        if isfield(MovieOptions.TracksText,'OnlyBehaviorColor_NoText')   
            if MovieOptions.TracksText.OnlyBehaviorColor_NoText
                OnlyBehaviorColor_NoText = true;
            end
        end  
    end

    % delta T = [1/frame rate], in seconds 
    if isfield(MovieOptions,'deltaT')    
        deltaT = MovieOptions.deltaT;
    else
        deltaT = 0.005;    
    end
    if isfield(MovieOptions,'deltaT_Collision')    
        deltaT_collision = MovieOptions.deltaT_Collision;
    else
        deltaT_collision = 1;    
    end
    
    % Zoom on a specific arena
    if isfield(MovieOptions,'arena')     
        arena = MovieOptions.arena;
    else
        arena = NaN;    
    end
    
    % This field below exists if the user is interested only in specific tracks. Only relevant Frames and arenas will be shown.  
    if isfield(MovieOptions,'ShowSpecificTracks')  
        Tracks2Show    = MovieOptions.ShowSpecificTracks; 
%         if ~ MovieOptions.ZoomOnWorm
%             % which arena
%             TrackArena = [Tracks.Arena];
%             RelevantArenas = unique(TrackArena(Tracks2Show));
%             if length(RelevantArenas)==1
%                 arena = RelevantArenas;
%                 disp(['All relevant tracks are in arena number ',num2str(arena),char(10)]);
%             end
%         end
        
        % which frames              
        if isfield(MovieOptions,'ShowEndOfTracks')
            FramesBefore = MovieOptions.ShowEndOfTracks.FramesBefore;
            FramesAfter  = MovieOptions.ShowEndOfTracks.FramesAfter;            
            Framefirst   = max( min(LastFrames(Tracks2Show)) - FramesBefore , 1);                             % FramesBefore collisions or frame number one; 
            Framelast    = min( max(FirstFrames(Tracks2Show))+ FramesAfter  , File.FragmentFrames(end,end));  % FramesAfter collisions or last;             
            Frames2Show  = Framefirst:Framelast;
            if isempty(Frames2Show)  % No collisions 
                disp('no collisions ??. Showing All frames in the track ');
                Frames2Show = min(FirstFrames(Tracks2Show)) : max(LastFrames(Tracks2Show)); 
            end
        else
            try 
                Frames2Show = min(FirstFrames(Tracks2Show)) : max(LastFrames(Tracks2Show));      
            catch
                Frames2Show = FirstFrames(1):LastFrames(end);
            end
        end
        
    end
    
    % This field exists if user is interested only in Frames. All relevant Tracks will be shown.  
    if isfield(MovieOptions,'ShowSpecificFrames')  
        if isfield(MovieOptions,'ShowSpecificTracks')  
            disp(['Both specific tracks and specific frames are defined. Showing the movie by the frames definition',char(10)]);            
        end
        Frames2Show = MovieOptions.ShowSpecificFrames;  % 
    elseif ~isfield(MovieOptions,'ShowSpecificFrames') && ~isfield(MovieOptions,'ShowSpecificTracks')
        All_Frames = unique([Tracks.Frames]);
        Frames2Show = All_Frames;    
    end
    
    % The relevant arena Zoom area: 
    if ~isnan(arena)
        ArenaBox =  File.TrackBoxAxis(arena,:);
        if isfield(MovieOptions,'EnlargeArena')
            Enlarge = MovieOptions.EnlargeArena;
            ArenaBox = [ArenaBox(1)-Enlarge, ArenaBox(2)+Enlarge, ArenaBox(3)-Enlarge, ArenaBox(4)+Enlarge];
        end
    end
    
    if isfield(MovieOptions,'TracksLine')  
        TracksLine = MovieOptions.TracksLine;
    else
        TracksLine = false;
    end   
   
    if isfield(MovieOptions,'ShowDetectedPixels')  
        ShowDetectedPixels = true;
        WormFeatures       = MovieOptions.ShowDetectedPixels;        % A cell with the name of ields to display
        WormFeaturesCOLORS = MovieOptions.ShowDetectedPixelsColors;  % Color for display
        OriginalFrameSize  = File.FrameSize;

    else
        ShowDetectedPixels = false;
    end
    
    if isfield(MovieOptions,'ShowDetectedROI')  
        ShowDetectedROIs   = true;
        load(File.FluorescenceFiles{1}, 'Perimeter_of_ROI', 'CorrelationWithReference', 'Frames', 'NonReliableIndices','GCaMP_Summary');
        MinCorrelation     = GCaMP_Summary(1).MinCorrelationThreshold;
        ROI.Frames         = Frames;
        ROI.Reliable       = (CorrelationWithReference > MinCorrelation) & (~NonReliableIndices);
        [~ , index_in_Frames2Show, index_in_ROI_Matrices] = intersect(Frames2Show, ROI.Frames );        %  [C,ia,ib] = intersect(A,B);
        
        [~,ReliablePart,~] = intersect(index_in_ROI_Matrices, find(ROI.Reliable));
        
        index_in_ROI_Matrices = index_in_ROI_Matrices(ReliablePart);
        index_in_Frames2Show  = index_in_Frames2Show(ReliablePart);
        
        % This vectors are of length (Frames2Show) and relate to the indices in Frames2Show (frame_number) 
        ROI.Frames2Show.ShowOrNot                              = false(1,length(Frames2Show)); 
        ROI.Frames2Show.ShowOrNot       (index_in_Frames2Show) = true;
        ROI.Frames2Show.index_in_ROI_Mat                       = zeros(1,length(Frames2Show));
        ROI.Frames2Show.index_in_ROI_Mat(index_in_Frames2Show) = index_in_ROI_Matrices;
        
        clear NonReliableIndices  Frames  CorrelationWithReference
        
        NumberOfROIs       = length(Perimeter_of_ROI);        % A cell with the name of ields to display
        ROI_COLORS         = MovieOptions.ShowDetectedROIColors;  % Color for display

    else
        ShowDetectedROIs   = false;
    end

    if isfield(MovieOptions,'ConstantScale')
        ConstantScale = MovieOptions.ConstantScale;
    end
    if ConstantScale
        Max = 0;
        Min = inf;
        for Frame = Frames2Show(1:end/20:end)
            switch VideoFormat
            case 2                             % A single multiple-tiff file
                Mov           = imread(FileFullName, Frame);     
            case 1                             % A sequence of tiff files
                FileFullName  = [MoviePath,'\',MovieFileNames{Frame}]; 
                Mov           = imread(FileFullName);                          
            case 3                             % A single avi movie file
                Mov           = read(MovieObj, Frame);     
                Mov           = Mov(:,:,1);                %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!
            case 4                             % A sequence of avi movie files
                CurrentMovieObjFile = File.MultiAviFrameConversion.MovieFileNumber(Frame);
                if (PreviousMovieObjFile ~= CurrentMovieObjFile)        % create a new file object if necessary
                    FileFullName        = [MoviePath,'\',MovieFileNames{CurrentMovieObjFile}]; 
                    PreviousMovieObjFile = CurrentMovieObjFile;
                    MovieObj             = VideoReader(FileFullName);
                end
                CurrentFrameNumberInFile = File.MultiAviFrameConversion.FrameNumberInFile(Frame);
                Mov                      = read(MovieObj, CurrentFrameNumberInFile);     
                Mov                      = Mov(:,:,1);          %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!                                
            end 
                                  
            Max = max([Max, max(Mov(:))]);
            Min = min([Min, min(Mov(:))]);
        end
        
        BestScale = [Min Max];    
        if isnumeric(ConstantScale)   
            if (ConstantScale<=100)&&(ConstantScale>=0)   % Scale based on percentage given by the user
%                 disp(['re-setting scale (',num2str(ConstantScale),'%)']);
%                 HalfRange = (Max-Min)/2;
%                 Middle    = (Max+Min)/2;                           
%                 BestScale = [Middle-(HalfRange*ConstantScale/100) Middle+(HalfRange*ConstantScale/100)];        
%                 BestScale = prctile(double(Mov(:)),[ConstantScale 100]);
%                 BestScale = prctile(double(Mov(:)),[ConstantScale 100-ConstantScale]);
                BestScale = prctile(double(Mov(:)),[ConstantScale 99.999]);
            end 
        end
    end
    
    CreateMovie = false;
    if isfield(MovieOptions,'CreateMovie')   
        if MovieOptions.CreateMovie
            CreateMovie = true;
            writerObj = VideoWriter('K:\testvideowriter1.avi','Uncompressed AVI'); 
            writerObj.FrameRate = 6;
            open(writerObj);
        end
    end
    
    ZoomOnWorm = false;
    if isfield(MovieOptions,'ZoomOnWorm')   
        if MovieOptions.ZoomOnWorm
            ZoomOnWorm   = true;
        end
    end 
    
    ShowTracksOnBackground=false;
    if isfield(MovieOptions,'ShowTracksOnBackground')   
        if MovieOptions.ShowTracksOnBackground
            ShowTracksOnBackground   = true;
        end
    end  
    
    if exist('BehaviorStructure','var')
        if ~isempty(BehaviorStructure)
            BehaviorFieldName = BehaviorStructure.fieldname_inTracks;
        end
    end
    
    Add_CTXdata_And_SmoothPosition = false;
    if isfield(MovieOptions,'ShowCTXposition')   
        if MovieOptions.ShowCTXposition
            if ~isempty(Data_Matrices)
                Add_CTXdata_And_SmoothPosition = true;
            end
        end
    end
    if isfield(MovieOptions,'ShowValueOfThisField')   
        if ~isempty(MovieOptions.ShowValueOfThisField)
            if ~isempty(Data_Matrices)
                Add_CTXdata_And_SmoothPosition = true;
                FramesField = 'FrameNumber';
                X_field     = 'Coordinates_X_Smoothed';
                Y_field     = 'Coordinates_Y_Smoothed';
                CTX_field = MovieOptions.ShowValueOfThisField;                
                disp(['Showing values of ', CTX_field]); 
             end
        end
    end    
    
    ShowBehaviorFrom_Data_Matrices = false;
    if isfield(MovieOptions,'ShowBehaviorFrom_Data_Matrices')   
        if MovieOptions.ShowBehaviorFrom_Data_Matrices
            if ~isempty(Data_Matrices)
                ShowBehaviorFrom_Data_Matrices = true;
                if ~exist('CTX_field','var') % if it was not previously defined in 'MovieOptions.ShowValueOfThisField' then show 'CTX_position'
                    if isfield(Data_Matrices,'FrameNumber_Smoothed2')
                        FramesField = 'FrameNumber_Smoothed2';
                        X_field     = 'X_Smoothed2';
                        Y_field     = 'Y_Smoothed2';
                        CTX_field   = 'CTX_position_Smoothed2';
                    elseif isfield(Data_Matrices,'FrameNumber')
                        disp('Showing ''raw'' CTX data --> The shown CTX  value are not yet corrected to the butanone odor side at that experiment')
                        FramesField = 'FrameNumber';
                        X_field     = 'Coordinates_X_Smoothed';
                        Y_field     = 'Coordinates_Y_Smoothed';
                        CTX_field   = 'CTX_position';
                    end    
                end
            end
        end
    end  
       
    
end

if ColorByID
%     COLORS_FOR_TRACKS   = {'w','g','b','m','c','y','k'}; % it will start from the green
    COLORS_FOR_TRACKS   = {'w','g','b','m','c','y',[0.8 0.5 0.2],[0.8 0.5 1]}; % it will start from the green
    if ColorByID==2
        TracksColorsIndices = mod(1:length(Tracks),length(COLORS_FOR_TRACKS))+1;  
    elseif ColorByID==1
        TracksColorsIndices = mod(1:length(unique([Tracks.ID])),length(COLORS_FOR_TRACKS))+1;  
    end    
end

if ShowTracksOnBackground
    if isempty(background)
        load(File.BackgroundFile,'background');
    end
    Mov = background;
end


%% Show Movie with tracks and their linking
% f=figure('position',get(0,'screensize'));  
% tic
for frame_number = 1:length(Frames2Show)
    Frame = Frames2Show(frame_number);
    
    if ~ShowTracksOnBackground        
        % Get Frame
        switch VideoFormat
            case 2                             % A single multiple-tiff file
                Mov           = imread(FileFullName, Frame);     
            case 1                             % A sequence of tiff files
                FileFullName  = [MoviePath,'\',MovieFileNames{Frame}]; 
                Mov           = imread(FileFullName);                          
            case 3                             % A single avi movie file
                Mov           = read(MovieObj, Frame);     
                Mov           = Mov(:,:,1);                %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!
            case 4                             % A sequence of avi movie files
                CurrentMovieObjFile = File.MultiAviFrameConversion.MovieFileNumber(Frame);
                if (PreviousMovieObjFile ~= CurrentMovieObjFile)        % create a new file object if necessary
                    FileFullName        = [MoviePath,'\',MovieFileNames{CurrentMovieObjFile}]; 
                    PreviousMovieObjFile = CurrentMovieObjFile;
                    MovieObj             = VideoReader(FileFullName);
                end
                CurrentFrameNumberInFile = File.MultiAviFrameConversion.FrameNumberInFile(Frame);
                Mov                      = read(MovieObj, CurrentFrameNumberInFile);     
                Mov                      = Mov(:,:,1);          %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!                                
        end 
    end
    
    if ConstantScale
        imshow(Mov,BestScale,'initialmagnification',600);
    else
        imshow(Mov,[],'initialmagnification',600);
    end
%     imshow(Mov,[],'initialmagnification',600);
    
    if ~ Add_CTXdata_And_SmoothPosition  
        if ~isnan(arena)        
            axis(ArenaBox); 
            title(['Arena ',num2str(arena),', Frame ',num2str(Frame)]);
        else
            title(['Frame ',num2str(Frame)]);
        end
    else
        if ~isnan(arena)        
            axis(ArenaBox); 
        end
    end
    if ZoomOnWorm   %%% Assuming there is only One Worm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
        FrameIndexInTrack = find([Tracks.Frames]==Frame,1);
        if ~isempty(FrameIndexInTrack)   
%             load(File.TrackingVariablesFile,'MaxWorm_size')
%             HalfWormMaxSize  = File.VariablesInformation.MaxWorm_size/0.05 * File.PixelSize / 2;      % File.VariablesInformation.MaxWorm_size ~ 1mm*0.05 mm
%             WormCentroid = Tracks.Path(FrameIndexInTrack,:);
%             Yaxis(1)   = max(1,           WormCentroid(1)-HalfWormMaxSize*1.5); 
%             Yaxis(2)   = min(size(Mov,1),WormCentroid(1)+HalfWormMaxSize*1.5); 
%             Xaxis(1)   = max(1,           WormCentroid(2)-HalfWormMaxSize*1.5); 
%             Xaxis(2)   = min(size(Mov,2),WormCentroid(2)+HalfWormMaxSize*1.5); 
%             Xaxis = round(Xaxis);
%             Yaxis = round(Yaxis);
%             axis([Xaxis Yaxis]);
%             axis([Yaxis Xaxis]);
            Xaxis = [min(Tracks.WormPerimeter.Xcoordinate{FrameIndexInTrack}) max(Tracks.WormPerimeter.Xcoordinate{FrameIndexInTrack})];
            Yaxis = [min(Tracks.WormPerimeter.Ycoordinate{FrameIndexInTrack}) max(Tracks.WormPerimeter.Ycoordinate{FrameIndexInTrack})];
            axis([Yaxis Xaxis]);
            
        end
    end
    
    if exist('Tracks2Show','var')
        RelevantTracks = Tracks2Show;
    else
        if Add_CTXdata_And_SmoothPosition
            RelevantTracks = 1:length(Tracks);
        else
            RelevantTracks = find(ALLFRAMES(Frame,:));
        end
    end

    hold on;
    All_CTX_position = [];
    for t_ind = 1:length(RelevantTracks);        
        t_num       = RelevantTracks(t_ind);
        IndexInTracks = find(Tracks(t_num).Frames==Frame,1);
        if isempty(IndexInTracks)
            Coordinates = [];
        else
            Coordinates = double(Tracks(t_num).Path(IndexInTracks,:)); 
            % Flip X-Y for imshow
            Coordinates = Coordinates([2 1]);
        end
        FrameCount  = Frame - FirstFrames(t_num)+1;
        
        % Add Track line of previous frames
        if TracksLine
            RelevantPathFrames = find(Tracks(t_num).Frames<=Frame);                      % Show all previous track path            
            if ~islogical(TracksLine) && (TracksLine < length(RelevantPathFrames)) % Show only the specified number of previous frames)
                RelevantPathFrames = RelevantPathFrames((end-TracksLine+1):end);       
            end
            Path = Tracks(t_num).Path(RelevantPathFrames,:);           
            line(Path(:,2), Path(:,1),'linewidth',1,'color','k','linestyle','-');
        end
        % ColumnIndex in Data_Matrices
        if Add_CTXdata_And_SmoothPosition  || ShowBehaviorFrom_Data_Matrices
             ColumnIndex       = find(Data_Matrices.(FramesField)(t_num,:) == Frame,1);   
            if isempty(ColumnIndex) 
%                 ColumnIndex       = find(Data_Matrices.FrameNumber(t_num,:) == (Frame-1),1);   
                ColumnIndex       = find(Data_Matrices.(FramesField)(t_num,:) == (Frame-1),1);   
            end     
        end

        % Defining display text for each track
        if ShowTracksText && ~isempty(Coordinates)
            text_for_Movie = '';
            if MovieOptions.TracksText.TrackNum
                text_for_Movie = [text_for_Movie,'(#',num2str(t_num),')'];
            end
            if MovieOptions.TracksText.ID1
                text_for_Movie = [text_for_Movie,',L',num2str(Tracks(t_num).OriginalTrack)];
            end
            if MovieOptions.TracksText.ID2
                text_for_Movie = [text_for_Movie,',LL',num2str(Tracks(t_num).Worm_ID)];
            end
            if MovieOptions.TracksText.Behavior & ~OnlyBehaviorColor_NoText
                if ShowBehaviorFrom_Data_Matrices &&  ~isempty(ColumnIndex)
                    CurrentBehCode =  Data_Matrices.(BehaviorFieldName)(t_num,ColumnIndex);
                else
                    CurrentBehCode =  Tracks(t_num).(BehaviorFieldName)(IndexInTracks);
                end
%                 CurrentBehColor=  BehaviorStructure.BehaviorCodeColors{BehaviorStructure.BehaviorCodeNumbers==CurrentBehCode};
                CurrentBehName =  BehaviorStructure.BehaviorCodeName{BehaviorStructure.BehaviorCodeNumbers==CurrentBehCode};
                if isempty(text_for_Movie)
                    text_for_Movie = CurrentBehName;
                else                    
                    text_for_Movie = [text_for_Movie,char(10),CurrentBehName];
                end
            
            end 
        else
%             text_for_Movie = num2str(t_num);
            text_for_Movie = '';
        end
        
        if Add_CTXdata_And_SmoothPosition  &&  ~isempty(ColumnIndex)
            Current_Coordinates_X_Smoothed = double(Data_Matrices.(X_field)(t_num,ColumnIndex));
            Current_Coordinates_Y_Smoothed = double(Data_Matrices.(Y_field)(t_num,ColumnIndex));
            Current_CTX_position           = Data_Matrices.(CTX_field)(t_num,ColumnIndex);
            All_CTX_position = [All_CTX_position, Current_CTX_position];

            CTXtext  = num2str(Current_CTX_position,2);
%                     % From G --> R through brown
%                     R        = abs((Current_CTX_position+1)/2); 
%                     G        = abs((Current_CTX_position-1)/2); 
%                     B        = 0; 
            % From G --> R through black
            if Current_CTX_position>0
                R = Current_CTX_position;
                G = 0;
            else
                G = -Current_CTX_position;
                R = 0;
            end
            B        = 0; 
            CTXcolor = [R G B];            
        else
            CTXtext  = '';
            CTXcolor = 'b';                               
        end                         
    
        
                
        CurrentColor = 'b';
        if ~isempty(Coordinates)
            if ColorByID ==2      % Color set by track number
                CurrentColor = COLORS_FOR_TRACKS{TracksColorsIndices(t_num)};
            elseif ColorByID ==1  % Color set by worm ID
                Worm_ID = Tracks(t_num).ID;
                CurrentColor = COLORS_FOR_TRACKS{TracksColorsIndices(Worm_ID)};
            elseif ColorByID ==3  % Color set by Chemotaxis index
                CurrentColor = CTXcolor;
            elseif ColorByID ==4  % Color set by Behavior code
                if ShowBehaviorFrom_Data_Matrices &&  ~isempty(ColumnIndex)
                    CurrentBehCode =  Data_Matrices.(BehaviorFieldName)(t_num,ColumnIndex);
                else
                    CurrentBehCode =  Tracks(t_num).(BehaviorFieldName)(IndexInTracks);
                end
                CurrentBehColor=  BehaviorStructure.BehaviorCodeColors{BehaviorStructure.BehaviorCodeNumbers==CurrentBehCode};
                CurrentColor   = CurrentBehColor;
                CTXcolor       = CurrentColor;
            end
        else
            if ColorByID ==4
                CTXcolor = 'y';
            end                
        end  
        
        if Add_CTXdata_And_SmoothPosition && ~isempty(CTXtext)
            text(Current_Coordinates_Y_Smoothed-30,Current_Coordinates_X_Smoothed-30, CTXtext,'FontWeight','bold','color',CTXcolor);             % Track count
            plot(Current_Coordinates_Y_Smoothed   ,Current_Coordinates_X_Smoothed,'o','markersize',6,'color',CTXcolor,'markerfacecolor',CTXcolor); 
        end
        if ~isempty(Coordinates)
            text(Coordinates(2)+10,Coordinates(1),text_for_Movie,'FontWeight','bold','color',CurrentColor);             % Track count
            plot(Coordinates(2)   ,Coordinates(1),'.','markersize',12,'color',CurrentColor); 
        end
        
%         if Add_CTXdata_And_SmoothPosition  
%             if ~isnan(arena)                       
%                 title(['Arena ',num2str(arena),', Frame ',num2str(Frame),', CTX ',num2str(mean(All_CTX_position),2)]);
%             else
%                 title(['Frame ',num2str(Frame),', CTX ',num2str(mean(All_CTX_position),2)]);
%             end
%         end        
        %% Plot detected pixels
        if  ShowDetectedPixels  && ~isempty(Coordinates)
            for f_num = 1:length(WormFeatures)
                CurrentField = WormFeatures{f_num};        % name of field to display
                CurrentColor = WormFeaturesCOLORS{f_num};  % color
                if ColorByID ==2      % Color set by track number
                    CurrentColor = COLORS_FOR_TRACKS{TracksColorsIndices(t_num)};
                elseif ColorByID ==1  % Color set by worm ID
                    Worm_ID = Tracks(t_num).ID;
                    CurrentColor = COLORS_FOR_TRACKS{TracksColorsIndices(Worm_ID)};
                elseif ColorByID ==3  % Color set by Chemotaxis index
                    CurrentColor = CTXcolor;                    
                elseif ColorByID ==4  % Color set by Behavior code
                    CurrentColor = CurrentBehColor;                    
                end                         
                Plot_WormPixels_on_Frame2(Tracks(t_num), Frame, CurrentField, CurrentColor, ZoomOnWorm)              
                
            end
        end       
        %% Plot detected ROI
        if  ShowDetectedROIs                       
            if ROI.Frames2Show.ShowOrNot(frame_number) 
                Current_frame_ind_in_ROI = ROI.Frames2Show.index_in_ROI_Mat(frame_number);
                for ROI_ind = 1:NumberOfROIs                
                    plot(Perimeter_of_ROI{ROI_ind}(Current_frame_ind_in_ROI,:,2), ...
                         Perimeter_of_ROI{ROI_ind}(Current_frame_ind_in_ROI,:,1),...
                         '.','color', ROI_COLORS{ROI_ind}, 'markersize',0.2)  ;  hold on
                end
            end
        end       
    end
    
    if Add_CTXdata_And_SmoothPosition  
        if ~isnan(arena)                       
            title(['Arena ',num2str(arena),', Frame ',num2str(Frame),', CTX ',num2str(mean(All_CTX_position),2)]);
        else
            title(['Frame ',num2str(Frame),', CTX ',num2str(mean(All_CTX_position),2)]);
        end
    end
        
    % Allow long pauses before a track ends and when it begins (within one frame)  
    if ismember(Frame,[LastFrames(RelevantTracks) FirstFrames(RelevantTracks)])
% pause;
        pause(deltaT_collision);
    elseif ismember(Frame+1,LastFrames(RelevantTracks)) || ismember(Frame-1, FirstFrames(RelevantTracks))
        pause(deltaT_collision);
    else
        pause(deltaT);
    end
    
    if CreateMovie    
        set(f,'position',get(0,'screensize'));  
%         writeVideo(writerObj,getframe);
        writeVideo(writerObj,getframe(f));
    end
        
    hold off;
end

if CreateMovie    
    close(writerObj);
end

% toc

return

function Plot_WormPixels_on_Frame2(CurrentTrack, Frame_Num, PixelsField, COLOR, ZoomOnWorm )

FrameIndex = find( CurrentTrack.Frames == Frame_Num, 1, 'first'); 
if ZoomOnWorm 
    MARKERSIZE  = 1;        
else   
    MARKERSIZE = 0.5;
end

if ~exist('COLOR','var')
    COLOR   = 'r';
end
LINESTYLE   = 'none'; % for PixelIdxList and HeadTail_PixelIdxList
MARKER      = '.';    % for PixelIdxList

if strcmpi(PixelsField,'Midline') 
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        hold on;
        plot(CurrentTrack.Midline.Y_coordinates_short{FrameIndex},CurrentTrack.Midline.X_coordinates_short{FrameIndex},'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);   
    end
elseif strcmpi(PixelsField,'Perimeter')
    MARKERSIZE = 2;
%     LINESTYLE   = '-';    
%     MARKER      = 'none'; 
    hold on;
    plot(single(CurrentTrack.WormPerimeter.Ycoordinate{FrameIndex}),single(CurrentTrack.WormPerimeter.Xcoordinate{FrameIndex}),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);    
elseif strcmpi(PixelsField,'HeadTail')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = '*';  
        if ZoomOnWorm 
            MARKERSIZE  = 8;
        else            
            MARKERSIZE  = 4; 
        end
        hold on;
        HeadTail_Matrix = squeeze(CurrentTrack.HeadTail(FrameIndex,:,:));
        plot(HeadTail_Matrix(:,2),HeadTail_Matrix(:,1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);           
    end
elseif strcmpi(PixelsField,'Head')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = '*';  
        if ZoomOnWorm 
            MARKERSIZE  = 8;
        else
            MARKERSIZE  = 4;
        end
        hold on;
        Head_Vec = CurrentTrack.Head(FrameIndex,:);
        plot(Head_Vec(2),Head_Vec(1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);         
    end
elseif strcmpi(PixelsField,'Tail')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = 'o';  
        if ZoomOnWorm 
            MARKERSIZE  = 8;
        else   
            MARKERSIZE  = 4; 
        end
        hold on;
        Tail_Vec = CurrentTrack.Tail(FrameIndex,:);
        plot(Tail_Vec(2),Tail_Vec(1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);     
    end
elseif strcmpi(PixelsField,'NeuronCoordinates')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = 'o';  
        if ZoomOnWorm 
            MARKERSIZE  = 8;
        else   
            MARKERSIZE  = 4; 
        end
        hold on;
        NeuronCoordinates_Vec = CurrentTrack.NeuronCoordinates(FrameIndex,:);
        plot(NeuronCoordinates_Vec(2),NeuronCoordinates_Vec(1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);     
    end   
    
end

return

function WormTrackingMovie2(stam, stam2, File, SpecificFrames, Tracks, gca_handle)

axes(gca_handle); 
% origin: WormTracking_MovieDisplay_Imaging_tracker_v01  

%% Parameters:

MovieOptions.ZoomOnWorm                   = true; %%% Assuming there is only One Worm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
MovieOptions.ShowCTXposition              = false;                     
MovieOptions.ShowBehaviorFrom_Data_Matrices = false;
MovieOptions.ShowTracksOnBackground       = false;     % Default=false; Tracks are shown by default on the original movie frames that are read from the path given in 'File'

% Other paremeters
MovieOptions.deltaT_Collision             = 0.05;      % Display time in regular frames
MovieOptions.deltaT                       = 0.05;      % Display time in regular frames
% MovieOptions.arena                        = 1;     
if exist('ArenaID','var');
    MovieOptions.arena                      = ArenaID;     
end
if exist('SpecificFrames','var')
    if ~isempty(SpecificFrames)
        MovieOptions.ShowSpecificFrames = SpecificFrames;
    end
end

MovieOptions.TracksLine                   = false;      % true/false/Num = show/(don't show)/(Show last Num of frames) of track line preceeding to each frame.                                                         
MovieOptions.TracksText.TrackNum          = true;      % Track text diplay settings
MovieOptions.TracksText.ID1               = false;     
MovieOptions.TracksText.ID2               = false;     
MovieOptions.TracksText.Behavior         = false; 
% MovieOptions.TracksText.OnlyBehaviorColor_NoText = false;   % true or false

MovieOptions.TracksText.ColorByID         = 0;     % No, same colors for all worms. different colors based on features 

% MovieOptions.ShowDetectedPixels           = {'Area','Skeleton','HeadTail'}; % A cell with the name of ields to display
% MovieOptions.ShowDetectedPixelsColors     = {'b','r','g'};                  % Color for display
MovieOptions.ShowDetectedPixels           = {'Perimeter'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'b'};                  % Color for display
MovieOptions.ShowDetectedPixels           = {'Perimeter','Midline'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'b','g'};                  % Color for display
MovieOptions.ShowDetectedPixels           = {'Perimeter','Midline','Head','Tail'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'b','g','r','m'};                      % Color for display
MovieOptions.ShowDetectedPixels           = {'Midline','Head','Tail'}; % A cell with the name of ields to display
MovieOptions.ShowDetectedPixelsColors     = {'g','r','m'};                      % Color for display

BehaviorStructure = [];
            
% MovieOptions.ConstantScale                = true;        % If true, the function will look for the best scaling and scale each frame similarly. If false --> imshow(Mov,[]) 
MovieOptions.ConstantScale                = 2;        % If true, the function will look for the best scaling and scale each frame similarly. If false --> imshow(Mov,[]) 
% MovieOptions.ConstantScale                = false;     
MovieOptions.CreateMovie                  = false;       % If true create a movie
Data_Matrices = [];
background = [];

SimpleDisplay_FluorescenceMovieWithTracks2(Tracks, File, MovieOptions, Data_Matrices, BehaviorStructure, background);  % Optional Movie generation

return

%% Extract neuron position
function Tracks  = Estimate_Neuron_Position_inline(PatternMatrix,Tracks, plotme)
MinimalNeuronProminence = 0.6;   % FREE PARAMETER
NeuronLocationInWormCoordinates        = [25 75];
DefaultNeuronLocationInWormCoordinates = [40 18];

if plotme
    FigureString    = {'Track 1, frame index ','Track 2, frame index ','Track 3, frame index '};
    FigureHandle(1) = figure('position',[32    50   1615   268]); 
    FigureHandle(2) = figure('position',[32    390  1615   268]); 
    FigureHandle(3) = figure('position',[32    730  1615   268]); 
end

for tr = 1:length(Tracks)
    Tracks(tr).NeuronInWormCoordinates_Estimations             = zeros(size(Tracks(tr).Path),'single')*NaN;
    Tracks(tr).NeuronInWormCoordinates_TrueIfGoodEstimation    = false(1,size(Tracks(tr).Path,1));
    Tracks(tr).NeuronLocationWasDetected                       = false(1,size(Tracks(tr).Path,1));
end

for tr = 1:length(Tracks)
    disp(['Searching for estimated neuron locations for Track ',num2str(tr)])
    for ind = 1:length(Tracks(tr).Frames)
        if plotme
            pause(0.2)
            CurrentFigureHandle = FigureHandle(tr);
            set(FigureHandle(tr),'name',[FigureString{tr},num2str(ind)]);
        else           
            CurrentFigureHandle = [];
        end           
        MAT = squeeze(PatternMatrix{tr}(ind,:,:)); 
        if any(~isnan(MAT(:)))
            [NeuronLocationWasDetected, NeuronCoordinates] = CheckPatternPeaks2...
               (MAT, NeuronLocationInWormCoordinates, DefaultNeuronLocationInWormCoordinates , MinimalNeuronProminence, CurrentFigureHandle); 
            % asign values to Tracks:
            Tracks(tr).NeuronInWormCoordinates_Estimations(ind,:) = NeuronCoordinates;
            Tracks(tr).NeuronLocationWasDetected(ind)             = NeuronLocationWasDetected;           
       end   
   end    
end

return

function [NeuronLocationWasDetected, NeuronCoordinates] = CheckPatternPeaks2(MAT, NeuronLocationInWormCoordinates, DefaultNeuronLocationInWormCoordinates , MinimalNeuronProminence, FigureHandle) 

NeuronCoordinates       = [NaN NaN];
LeftSideEnd             = size(MAT,2)/5;
Vec                     = mean(MAT(:,1:LeftSideEnd),1);
[~,locs,~,proms] = findpeaks(double(Vec));

DynamicRange         = max(Vec)-min(Vec); 
NormalizedProminence = proms / DynamicRange;
if length(locs)>=1  % disregard small peaks
    locs  = locs(NormalizedProminence>MinimalNeuronProminence);
    proms = proms(NormalizedProminence>MinimalNeuronProminence);
end

if length(locs)==1  
    if (locs>NeuronLocationInWormCoordinates(1)) && (locs<NeuronLocationInWormCoordinates(2))
        NeuronCoordinates(1) =  locs;
    end
elseif length(locs)>1    
    [~,sortorder]          = sort(proms);
    LargestPeaksIndices    = sort(sortorder(end:(-1):(end-1))); % [first largest peak, second largest peak];
    LargestPeaksLocations  = locs(LargestPeaksIndices);  
    if     (LargestPeaksLocations(1)>NeuronLocationInWormCoordinates(1)) && (LargestPeaksLocations(1)<NeuronLocationInWormCoordinates(2))
        % choose first peak
        NeuronCoordinates(1) =  LargestPeaksLocations(1);
    elseif (LargestPeaksLocations(2)>NeuronLocationInWormCoordinates(1)) && (LargestPeaksLocations(2)<NeuronLocationInWormCoordinates(2))
        % choose first peak
        NeuronCoordinates(1) =  LargestPeaksLocations(2); 
    end
end
NeuronLocationWasDetected = ~isnan(NeuronCoordinates(1));
       
if NeuronLocationWasDetected
    [~, NeuronWidthCoordinate] = nanmax(MAT(:,NeuronCoordinates(1)));
    NeuronCoordinates(2)       = NeuronWidthCoordinate;
else
    NeuronCoordinates = DefaultNeuronLocationInWormCoordinates;
end
    
% if FigureHandle exists     
if exist('FigureHandle','var') && ~isempty(FigureHandle)
    figure(FigureHandle)
    subplot(1,2,1);   
        hold off;
        imshow(MAT,[]); hold on;       
        if NeuronLocationWasDetected
            plot(NeuronCoordinates(1),NeuronCoordinates(2),'ro');
        else
            plot(NeuronCoordinates(1),NeuronCoordinates(2),'bo');
        end            
            
    subplot(1,2,2);
        hold off;
        plot(Vec,'k-');   hold on; 
        if NeuronLocationWasDetected
            plot(NeuronCoordinates(1),Vec(NeuronCoordinates(1)),'r*')
        else
            plot(NeuronCoordinates(1),Vec(NeuronCoordinates(1)),'b*')
        end        
        xlim([0 length(Vec)])
end   

return

function Tracks  = ExtractEstimatedNeuronPositionInRealCoordinates(Tracks,File, GoodLinks_ALL, Non_Deleted_Tracks_BeforeStitching, Non_Deleted_Tracks_BeforeManualInspection, TracksCorrectionParameters)
%% Update pattern matrix with TracksLinking. Do not store information of unreliable midlines 
VariableName    = 'PositionMatrixX';
PositionMatrixX = LinkPatternMatrix (File, VariableName, GoodLinks_ALL, Non_Deleted_Tracks_BeforeStitching, Non_Deleted_Tracks_BeforeManualInspection, TracksCorrectionParameters);
PositionMatrixX = CorrectPatternMatrix_UsingHeadTailCoordinates    (Tracks, PositionMatrixX);
VariableName    = 'PositionMatrixY';
PositionMatrixY = LinkPatternMatrix (File, VariableName, GoodLinks_ALL, Non_Deleted_Tracks_BeforeStitching, Non_Deleted_Tracks_BeforeManualInspection, TracksCorrectionParameters);
PositionMatrixY = CorrectPatternMatrix_UsingHeadTailCoordinates    (Tracks, PositionMatrixY);

for tr = 1:length(Tracks)
    Tracks(tr).NeuronCoordinates_Estimations     = zeros(size(Tracks(tr).Path),'single')*NaN;
    NeuronInWormCoordinates_Estimations          = Tracks(tr).NeuronInWormCoordinates_Estimations;
    relevant_indices = find(~isnan(NeuronInWormCoordinates_Estimations(:,1)));
    for ind = relevant_indices'
        CurrentWormCoordinates = NeuronInWormCoordinates_Estimations(ind,:);
        Tracks(tr).NeuronCoordinates_Estimations(ind,1) = PositionMatrixX{tr}(ind, CurrentWormCoordinates(2),CurrentWormCoordinates(1)); 
        Tracks(tr).NeuronCoordinates_Estimations(ind,2) = PositionMatrixY{tr}(ind, CurrentWormCoordinates(2),CurrentWormCoordinates(1)); 
    end    
end

return

%% Extract neuron fluorescence
function File = ExtractNeuralActivityFromTracksAndOriginalMovie_inline(File, Tracks, WormArea)

%% Initialization and parameters
t0             = clock;
MoviePath      = File.MoviePath;
MovieFileNames = File.MovieFileNames;
PixelSize      = File.PixelSize;              % PixelSize == # pixels per mm (1D length) 
FrameRate      = File.FrameRate;              % PixelSize == # pixels per mm (1D length) 
NumOfFrames    = File.NumberOfFrames;
FrameSize      = File.FrameSize;
NumOfTracks    = length(Tracks);
LimitedSearchDistanceInMicrometers  = 25;
LimitedSearchDistanceInPixels       = floor(LimitedSearchDistanceInMicrometers/(1000/PixelSize));
SearchDistanceFromPreviousFrameInMicrometers = 40;
SearchDistanceFromPreviousFrameInPixels      = floor(SearchDistanceFromPreviousFrameInMicrometers/(1000/PixelSize));
SearchDistanceFromPreviousFrameInMicrometersWhenHeadIsSqueezed= 15;
SearchDistanceFromPreviousFrameInPixelsWhenHeadIsSqueezed     = floor(SearchDistanceFromPreviousFrameInMicrometersWhenHeadIsSqueezed/(1000/PixelSize));
DisplaySquareSize                   = 25; % pixels
                               
%% FramesWithinTracks matrix
mismatch           = false;
ConvertMovieFrameToIndexWithinTracks = zeros(NumOfTracks,NumOfFrames,'single')*NaN;
for tr=1:NumOfTracks
    TotalNumOfFramesInTrack                  = length(Tracks(tr).Frames);
    ConvertMovieFrameToIndexWithinTracks(tr,Tracks(tr).Frames) = 1:TotalNumOfFramesInTrack;
    if TotalNumOfFramesInTrack~=length(WormArea(tr).LinearIndices) 
        mismatch = true;
    end
end
if mismatch==true
    disp('ERROR !! mismatch between ''WormArea'' and ''Tracks'' structures. Aborting...')
    return
end                 
    
if     strcmpi(File.VideoFormat,'tiff')
    VideoFormat   = 1;                                  % read new file every loop iteration 
elseif strcmpi(File.VideoFormat,'multi-tiff')
    VideoFormat   = 2;
    FileFullName  = [MoviePath,'\',MovieFileNames{1}];  % one file
elseif strcmpi(File.VideoFormat,'avi')
    VideoFormat   = 3;
    FileFullName  = [MoviePath,'\',MovieFileNames{1}];  % one file
    MovieObj      = VideoReader(FileFullName);          % one object
elseif strcmpi(File.VideoFormat,'multi-avi')
    VideoFormat   = 4;
    PreviousMovieObjFile = 0;                           % read new file object if and when it is needed 
end
           
%% Load background, Mask and vignetting pattern
load(File.BackgroundFile, 'Mask', 'VignettingPattern','background');
dilatedMask = ~imdilate(~Mask,strel('disk',1,0));
VignettingPattern(~dilatedMask)=NaN;    
%  figure; imagesc(VignettingPattern); set(gca,'clim',[0 1]); axis equal; colorbar

%% Find worm objects in each frame
NeuronCoordinatesSearchFlag                  = zeros(NumOfTracks,NumOfFrames,'uint8');  % 1=search near estimated position, 2=search near previous position 
ReliableNeuronPosition                       = false(NumOfTracks,NumOfFrames);
FlagTrueIfNeuronPositionMatchPreviousFrame   = false(NumOfTracks,NumOfFrames);
HeadIsSqueezedMAT                            = false(NumOfTracks,NumOfFrames);
NeuronCoordinatesMatrix                      = zeros(NumOfTracks,NumOfFrames,2,'single')*NaN;
NeuronCoordinatesMAT_NoSearch                = zeros(NumOfTracks,NumOfFrames,2,'single')*NaN;
NeuronCoordinatesMAT_AllowLimitedSearch      = zeros(NumOfTracks,NumOfFrames,2,'single')*NaN;
NeuronCoordinatesMAT_SearchFromPreviousFrame = zeros(NumOfTracks,NumOfFrames,2,'single')*NaN;
NeuronDisplayMatrix                          = zeros(NumOfTracks,NumOfFrames,DisplaySquareSize,DisplaySquareSize,'single')*NaN;
PreviousNeuronCoordinatesFrame = zeros(1,NumOfTracks)*NaN;
MinimumPosition                = 3;
MaximumPosition_vec            = FrameSize-3;

h = waitbar(0,['analyzing neural activity in ',num2str(NumOfFrames),' frames']); 
for frame = 1:NumOfFrames     
    waitbar(frame/NumOfFrames,h);
    %% Get Frame and subtract background       
    switch VideoFormat
        case 2                             % A single multiple-tiff file
            Mov           = imread(FileFullName, frame);     
        case 1                             % A sequence of tiff files
            FileFullName  = [MoviePath,'\',MovieFileNames{frame}]; 
            Mov           = imread(FileFullName);                          
        case 3                             % A single avi movie file
            Mov           = read(MovieObj, frame);     
            Mov           = Mov(:,:,1);                %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!
        case 4                             % A sequence of avi movie files
            CurrentMovieObjFile = File.MultiAviFrameConversion.MovieFileNumber(frame);
            if (PreviousMovieObjFile ~= CurrentMovieObjFile)        % create a new file object if necessary
                FileFullName        = [MoviePath,'\',MovieFileNames{CurrentMovieObjFile}]; 
                PreviousMovieObjFile = CurrentMovieObjFile;
                MovieObj             = VideoReader(FileFullName);
            end
            CurrentFrameNumberInFile = File.MultiAviFrameConversion.FrameNumberInFile(frame);
            Mov                      = read(MovieObj, CurrentFrameNumberInFile);     
            Mov                      = Mov(:,:,1);          %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!                                
    end            
    Mov = imsubtract(Mov,background);    
              
    %% Keep Only Worm Objects and correct for vignetting patterns
    NormFrame = zeros(FrameSize,'single')*NaN;
    for tr=1:NumOfTracks
        ind = ConvertMovieFrameToIndexWithinTracks(tr,frame);
        if ~isnan(ind)
            NormFrame(WormArea(tr).LinearIndices{ind}) = Mov(WormArea(tr).LinearIndices{ind});
        end
    end    
    NormFrame = NormFrame./VignettingPattern;
    
    %% Find neuron coordinates    
    for tr=1:NumOfTracks
        ind                                   = ConvertMovieFrameToIndexWithinTracks(tr,frame); 
        if isnan(ind)                    
            continue
        end
        NeuronLocationWasDetected             = Tracks(tr).NeuronLocationWasDetected(ind);
        NeuronCoordinates_Estimations         = Tracks(tr).NeuronCoordinates_Estimations(ind,:);        
        NeuronCoordinates_Estimations_rounded = round(NeuronCoordinates_Estimations);
        
        if ~isnan(NeuronCoordinates_Estimations(1))
            NeuronCoordinatesMAT_NoSearch(tr,frame,:)           = NeuronCoordinates_Estimations_rounded;
            [Position, Distance_LimitedSearch]                  = FindLocalMaximum (NormFrame,NeuronCoordinates_Estimations_rounded, LimitedSearchDistanceInPixels);
            Position(1)                                         = min([Position(1)  MaximumPosition_vec(1)]);
            Position(2)                                         = min([Position(2)  MaximumPosition_vec(2)]);          
            NeuronCoordinatesMAT_AllowLimitedSearch(tr,frame,:) = Position;
        else
            Position = []; 
        end
                        
        PreviousFrame  = PreviousNeuronCoordinatesFrame(tr);
        FramesInterval = frame-PreviousFrame;
        if FramesInterval < FrameRate/2
            PreviousPosition                    = squeeze(NeuronCoordinatesMatrix(tr,PreviousFrame,:))';             
            [Position2, Distance2]              = FindLocalMaximum (NormFrame,PreviousPosition, SearchDistanceFromPreviousFrameInPixels);
            Position2(Position2<MinimumPosition)= MinimumPosition;  
            Position2(1)                        = min([Position2(1)  MaximumPosition_vec(1)]);
            Position2(2)                        = min([Position2(2)  MaximumPosition_vec(2)]);          
            NeuronCoordinatesMAT_SearchFromPreviousFrame(tr,frame,:) = Position2;
        else
            Position2 = [];                      
        end
        
        if ~isempty(Position) && ~isempty(Position2)               
            TrueIfBothPositionEstimationMatches                  = sqrt(sum((Position - Position2).^2)) <= 2;  % less than 2 pixels away from each other
            FlagTrueIfNeuronPositionMatchPreviousFrame(tr,frame) = TrueIfBothPositionEstimationMatches;
            
            % default values, use new position -- unless the conditions below suggest to use the previous position 
            NeuronCoordinatesSearchFlag(tr,frame)     = 1; 
            NeuronCoordinatesMatrix(tr,frame,:)       = NeuronCoordinatesMAT_AllowLimitedSearch(tr,frame,:);
            if TrueIfBothPositionEstimationMatches
                ReliableNeuronPosition(tr,frame)      = true;
            else  % Which one is a better approximation?  check if one is reliable while another is not...               
                if ~ReliableNeuronPosition(PreviousFrame)  % previous position is not reliable >> use new position                    
                    if NeuronLocationWasDetected && (Distance_LimitedSearch<=2)
                        ReliableNeuronPosition(tr,frame)  = true;
                    end
                elseif ReliableNeuronPosition(PreviousFrame) && (Distance2<=4) && ~(NeuronLocationWasDetected && (Distance_LimitedSearch<=2))
                    % previous position is reliable and the new position is not reliable 
                    NeuronCoordinatesMatrix(tr,frame,:)   = NeuronCoordinatesMAT_SearchFromPreviousFrame(tr,frame,:);
                    NeuronCoordinatesSearchFlag(tr,frame) = 2; 
                    ReliableNeuronPosition(tr,frame)      = true;                    
                end                                       
            end
            
        elseif ~isempty(Position) % only estimation based on current position is available
            NeuronCoordinatesMatrix(tr,frame,:)   = NeuronCoordinatesMAT_AllowLimitedSearch(tr,frame,:);
            NeuronCoordinatesSearchFlag(tr,frame) = 1; 
            if NeuronLocationWasDetected && (Distance_LimitedSearch<=2)
                ReliableNeuronPosition(tr,frame)  = true;
            end
            
        elseif ~isempty(Position2) % only estimation based on previous position is available
            NeuronCoordinatesMatrix(tr,frame,:)   = NeuronCoordinatesMAT_SearchFromPreviousFrame(tr,frame,:);
            NeuronCoordinatesSearchFlag(tr,frame) = 2; 
            if ReliableNeuronPosition(tr,PreviousFrame) &&  (Distance2<=4)
                ReliableNeuronPosition(tr,frame)  = true;
            end            
        else            
            continue    % No neuron was found >> don't update any matrix       
        end
                     
        PreviousNeuronCoordinatesFrame(tr) = frame;
        DisplayMatrix                      = AssignToDisplay(NormFrame, NeuronCoordinatesMatrix(tr,frame,:), DisplaySquareSize); 
        
        %% Quality check: if head is squeezed at a "omega-trial" state then the segmentation is unreliable in the approach below is necessary: 
        HeadIsSqueezed = sum(isnan(DisplayMatrix(:)))/ (DisplaySquareSize^2) < 0.5; % 0.3; 

        if HeadIsSqueezed
           HeadIsSqueezedMAT(tr,frame) = true;
           % Use other methods to find neuron:
           % A. Do not use current midline etimations      
           % B. allow very small distance changes from previous frame location      
           if FramesInterval < FrameRate/2
                [Position3, Distance3]              = FindLocalMaximum (NormFrame,PreviousPosition, ...
                                                                        SearchDistanceFromPreviousFrameInPixelsWhenHeadIsSqueezed); % only 2 pixels
                Position3(Position3<MinimumPosition)= MinimumPosition;  
                Position3(1)                        = min([Position3(1)  MaximumPosition_vec(1)]);
                Position3(2)                        = min([Position3(2)  MaximumPosition_vec(2)]);          
                NeuronCoordinatesMAT_SearchFromPreviousFrame(tr,frame,:) = Position3;
                NeuronCoordinatesMatrix(tr,frame,:)                      = Position3;
                NeuronCoordinatesSearchFlag(tr,frame)                    = 2;      
                DisplayMatrix                               = AssignToDisplay(NormFrame, NeuronCoordinatesMatrix(tr,frame,:), DisplaySquareSize); 
           end                   
        end
        NeuronDisplayMatrix(tr,frame,:,:) = DisplayMatrix; 
    end
end
 
% Assign Position to Tracks:
for tr=1:NumOfTracks
   Tracks(tr).NeuronCoordinates        = squeeze(NeuronCoordinatesMatrix(tr,Tracks(tr).Frames,:));     
   Tracks(tr).Neuron.CoordinatesMatrix = squeeze(NeuronCoordinatesMatrix(tr,:,:));     
   Tracks(tr).Neuron.DisplayMatrix     = squeeze(NeuronDisplayMatrix(tr,:,:,:));     
   Tracks(tr).HeadIsSqueezed           = HeadIsSqueezedMAT(tr,:);     
end

%% SAVE
arena = 1;
NeuronsDataBeforeManualCorrectionFileName = [File.TrackFile(1:end-4),'_NeuronsDataBeforeManualCorrection_Arena',num2str(arena),'.mat'];
tic
save(NeuronsDataBeforeManualCorrectionFileName,'File','Tracks','WormArea','background','VignettingPattern',...
    'NeuronCoordinatesSearchFlag','ReliableNeuronPosition','FlagTrueIfNeuronPositionMatchPreviousFrame',...
    'NeuronCoordinatesMatrix','NeuronCoordinatesMAT_NoSearch','NeuronCoordinatesMAT_AllowLimitedSearch','NeuronCoordinatesMAT_SearchFromPreviousFrame',...
    'NeuronDisplayMatrix', 'ConvertMovieFrameToIndexWithinTracks', 'HeadIsSqueezedMAT', '-v7.3');      
toc

%% elapsed time timing
t1=clock; 
disp(['Total cpu time = ',num2str(etime(t1,t0)/60),' minutes',char(10)])

return

function [Position, Distance] = FindLocalMaximum (Mov,NeuronCoordinates_Estimations, MaximalDistanceInPixels)
MaximalDistanceInPixels1 = MaximalDistanceInPixels;
MaximalDistanceInPixels2 = MaximalDistanceInPixels;

while ((NeuronCoordinates_Estimations(1)-MaximalDistanceInPixels1)<1)||((NeuronCoordinates_Estimations(1)+MaximalDistanceInPixels1)>size(Mov,1))
    MaximalDistanceInPixels1 = MaximalDistanceInPixels1-1;
end
while ((NeuronCoordinates_Estimations(2)-MaximalDistanceInPixels2)<1)||((NeuronCoordinates_Estimations(2)+MaximalDistanceInPixels2)>size(Mov,2))
    MaximalDistanceInPixels2 = MaximalDistanceInPixels2-1;
end
    
MAT = Mov((NeuronCoordinates_Estimations(1)-MaximalDistanceInPixels1):(NeuronCoordinates_Estimations(1)+MaximalDistanceInPixels1), ...
          (NeuronCoordinates_Estimations(2)-MaximalDistanceInPixels2):(NeuronCoordinates_Estimations(2)+MaximalDistanceInPixels2));
[M,I]     = max(MAT(:));
[I,J]     = ind2sub(size(MAT),I);
MATcenter = [MaximalDistanceInPixels1  MaximalDistanceInPixels2]+1;

Distance = sqrt(sum( (MATcenter - [I,J]).^2 ));
Position = round(NeuronCoordinates_Estimations) + [I,J] - MATcenter;

return

function DisplayMatrix = AssignToDisplay(Mov, NeuronCoordinates, DisplaySquareSize)
DisplayMatrix           = zeros(DisplaySquareSize,DisplaySquareSize,'single')*NaN;
MaximalDistanceInPixels = ones(2)*(DisplaySquareSize-1)/2;
MidPixelInDisplay       = (DisplaySquareSize+1)/2;

while ((NeuronCoordinates(1)-MaximalDistanceInPixels(1,1))<1)
    MaximalDistanceInPixels(1,1) = MaximalDistanceInPixels(1,1)-1;
end
while ((NeuronCoordinates(2)-MaximalDistanceInPixels(2,1))<1)
    MaximalDistanceInPixels(2,1) = MaximalDistanceInPixels(2,1)-1;
end
while ((NeuronCoordinates(1)+MaximalDistanceInPixels(1,2))>size(Mov,1))
    MaximalDistanceInPixels(1,2) = MaximalDistanceInPixels(1,2)-1;
end
while ((NeuronCoordinates(2)+MaximalDistanceInPixels(2,2))>size(Mov,2))
    MaximalDistanceInPixels(2,2) = MaximalDistanceInPixels(2,2)-1;
end

MAT = Mov((NeuronCoordinates(1)-MaximalDistanceInPixels(1,1)):(NeuronCoordinates(1)+MaximalDistanceInPixels(1,2)), ...
          (NeuronCoordinates(2)-MaximalDistanceInPixels(2,1)):(NeuronCoordinates(2)+MaximalDistanceInPixels(2,2)));


DisplayMatrix((MidPixelInDisplay-MaximalDistanceInPixels(1,1)):(MidPixelInDisplay+MaximalDistanceInPixels(1,2)), ...
              (MidPixelInDisplay-MaximalDistanceInPixels(2,1)):(MidPixelInDisplay+MaximalDistanceInPixels(2,2))) = MAT;
         
return





