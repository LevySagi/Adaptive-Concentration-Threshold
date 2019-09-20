function [TrackMorphologyData, CurrentTrack] = AnalyzeMorphologyProperties_ImagingSetup_v3(CurrentTrack, Settings)
% Input 
% CurrentTrack = Track structure of one worm.    
% Settings     = parameters, smoothing information, and plot information. Defined by the user using the Settings function.   

% Output 
% TrackLocomotionData. A structure with the track locomotion information. See full decription of fields below.

%% Settings: parameters, smoothing information, and plot information
PlotInfo                = Settings.PlotInfo;                                                          % Information about which plots to show, and whether to show the movie. 

CurrentTrack.PerimeterOverSize           = single(CurrentTrack.PerimeterLength) ./ single(CurrentTrack.Size);
CurrentTrack.NormalizedSize              = single(CurrentTrack.Size)              / nanmedian(single(CurrentTrack.Size));        
CurrentTrack.NormalizedSkeletonLength    = single(CurrentTrack.SkeletonLength)    / nanmedian(single(CurrentTrack.SkeletonLength));        
CurrentTrack.NormalizedPerimeterLength   = single(CurrentTrack.PerimeterLength)   / nanmedian(single(CurrentTrack.PerimeterLength));        
CurrentTrack.NormalizedPerimeterOverSize = single(CurrentTrack.PerimeterOverSize) / nanmedian(single(CurrentTrack.PerimeterOverSize));   

Fields = {'Size','SkeletonLength','PerimeterLength','Eccentricity','PerimeterOverSize',...    
          'NormalizedSize','NormalizedSkeletonLength','NormalizedPerimeterLength','NormalizedPerimeterOverSize',...
          'HeadIsSqueezed'};
for f_ind = 1:length(Fields)
    fieldname = Fields{f_ind};
    TrackMorphologyData.(fieldname) = single(CurrentTrack.(fieldname));
end

%% Plots
ScreenSize       = get(0,'screensize');

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
        legend('normalized perimeter','Eccentricity','Head is squeezed');        
    else
        legend('normalized perimeter','Eccentricity');        
    end
    ylim([0.5 1.2])

    figure('name','Morphology All','position',ScreenSize); 
    plot(TrackMorphologyData.NormalizedPerimeterLength,'bo-'); hold on;  
    plot(TrackMorphologyData.Eccentricity,'k.-');
    legend('normalized perimeter','Eccentricity');
    ylim([0.2 1.2])
end

return








