function [TrackMorphologyData, CurrentTrack] = AnalyzeMorphologyProperties_HighThroughputDevice_02(CurrentTrack, Settings, TracksStats)

% Input 
% coordinates = XY path of the track. column 1 and 2 corresponds to X and Y, respectively.  
% Settings    = parameters, smoothing information, and plot information. Defined by the user using the Settings function.   

% Output 
% TrackLocomotionData. A structure with the track locomotion information. See full decription of fields below.

%% Settings: parameters, smoothing information, and plot information
PlotInfo                              = Settings.PlotInfo;  % Information about which plots to show
CurrentTrack.PerimeterOverSize        = single(CurrentTrack.PerimeterLength) ./ single(CurrentTrack.Size);
TrackMorphologyData.Size              = CurrentTrack.Size;
TrackMorphologyData.SkeletonLength    = CurrentTrack.SkeletonLength;
TrackMorphologyData.PerimeterLength   = CurrentTrack.PerimeterLength;
TrackMorphologyData.Eccentricity      = CurrentTrack.Eccentricity;
TrackMorphologyData.PerimeterOverSize = CurrentTrack.PerimeterOverSize;

%% Normalized_sizes vectors: 
%   Short Tracks: Normalize the size relative to the median of all tracks. 
%   Long tracks:  Normalize with repect to the median in the track.
%   NOTE: ALLOW SIZE NORMALIZATION PER ARENA !!

MinTrackLengthForSizeNormalization = Settings.Thresholds.MinTrackLengthForSizeNormalization;         % seconds. A good guess is tracks of 20 seconds.
MinTrackFramesForSizeNormalization = MinTrackLengthForSizeNormalization * Settings.FrameRate;            % Frames 

if CurrentTrack.TrackLength < MinTrackFramesForSizeNormalization 
    TrackMorphologyData.NormalizedSize              = single(CurrentTrack.Size)              / TracksStats.MedianSize;        
    TrackMorphologyData.NormalizedSkeletonLength    = single(CurrentTrack.SkeletonLength)    / TracksStats.MedianSkeleton;        
    TrackMorphologyData.NormalizedPerimeterLength   = single(CurrentTrack.PerimeterLength)   / TracksStats.MedianPerimeter;        
    TrackMorphologyData.NormalizedPerimeterOverSize = single(CurrentTrack.PerimeterOverSize) / TracksStats.MedianPerimeterOverSize;        
else
    TrackMorphologyData.NormalizedSize              = single(CurrentTrack.Size)              / nanmedian(single(CurrentTrack.Size));        
    TrackMorphologyData.NormalizedSkeletonLength    = single(CurrentTrack.SkeletonLength)    / nanmedian(single(CurrentTrack.SkeletonLength));       
    TrackMorphologyData.NormalizedPerimeterLength   = single(CurrentTrack.PerimeterLength)   / nanmedian(single(CurrentTrack.PerimeterLength));    
    TrackMorphologyData.NormalizedPerimeterOverSize = single(CurrentTrack.PerimeterOverSize) / nanmedian(single(CurrentTrack.PerimeterOverSize));    
    if isfield(CurrentTrack,'HeadIsSqueezed')
        TrackMorphologyData.HeadIsSqueezed = CurrentTrack.HeadIsSqueezed(CurrentTrack.Frames); 
    end        
end            

%% Plots
ScreenSize = get(0,'screensize');
if PlotInfo.Morphology_All
    figure('name','Morphology All','position',ScreenSize); 
    plot(TrackMorphologyData.NormalizedPerimeterLength,'bo-'); hold on;  
    plot(TrackMorphologyData.NormalizedSkeletonLength,'rx-');
    plot(TrackMorphologyData.NormalizedSize,'g.-');
    plot(TrackMorphologyData.Eccentricity,'k.-');
    plot(TrackMorphologyData.NormalizedPerimeterOverSize,'m.-');
    if isfield(TrackMorphologyData,'HeadIsSqueezed')
        plot(TrackMorphologyData.HeadIsSqueezed,'c.');
        legend('normalized perimeter','normalized skeleton','normalized size','Eccentricity','normalized (perimeter/size)','Head is squeezed');
    else        
        legend('normalized perimeter','normalized skeleton','normalized size','Eccentricity','normalized (perimeter/size)');
    end       
    ylim([0.5 1.2])
end
if PlotInfo.Eccentricity_Perimeter
    figure('name','Morphology All','position',ScreenSize); 
    plot(TrackMorphologyData.NormalizedPerimeterLength,'bo-'); hold on;  
    plot(TrackMorphologyData.Eccentricity,'k.-');
    if isfield(TrackMorphologyData,'HeadIsSqueezed')
        plot(TrackMorphologyData.HeadIsSqueezed,'c.');
    end
    legend('normalized perimeter','Eccentricity');
    ylim([0.5 1.2])

    figure('name','Morphology All','position',ScreenSize); 
    plot(TrackMorphologyData.NormalizedPerimeterLength,'bo-'); hold on;  
    plot(TrackMorphologyData.Eccentricity,'k.-');
    legend('normalized perimeter','Eccentricity');
    ylim([0.2 1.2])
end

return








