function [BehaviorCode, Segments_HighLevelBehavior, Segments_LowLevelBehavior, Segments] = SegmentBehavior_SingleTrack_ImagingSetup_v3(CurrentTrack, tr_num, Settings)

% September 2019, Sagi Levy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Parameters initialization
FrameRate          = Settings.FrameRate;                      % FrameRate = Frames per second.
PixelSize          = Settings.PixelSize;                      % PixelSize = # pixels per mm.  
PlotSummaryFigures = false;   
PlotConfiguration  = Settings.PlotInfo.PlotConfiguration;     % Internal choice: When ploting the results, which plot configurations to show?      
                                                              %   PlotConfiguration = 0; % No Plots
                                                              %   PlotConfiguration = 1; % 6 plots
                                                              %   PlotConfiguration = 2; % 4 plots
Thresholds                = Settings.Thresholds;

Frames                    = CurrentTrack.Frames ;
FramesOutOfBounds         = CurrentTrack.OutOfBounds ;
VelocityAmplitude         = CurrentTrack.Speed; 
VelocityAngle             = CurrentTrack.AngularDisplacement; 
ChangeInVelocityAngle     = CurrentTrack.AngularVelocity; 
Eccentricity              = CurrentTrack.Eccentricity; 
NormalizedSize            = CurrentTrack.NormalizedSize;
NormalizedPerimeterLength = CurrentTrack.NormalizedPerimeterLength;
NormalizedPerimeterOverSize = CurrentTrack.NormalizedPerimeterOverSize;
VelocityContinuousAngle   = CurrentTrack.AngularDisplacement_ContinuousAngle; 
ChangeInVelocityAmplitude = CurrentTrack.Acceleration; 
Path                      = CurrentTrack.Path; 
Head                      = CurrentTrack.Head; 
Tail                      = CurrentTrack.Tail; 

NumOfFrames               = length(Frames);

%%% Convert position and time units of all threholds given by the user (from the Settings.Threholds variable) to units of [pixels] and [frames]      

Thresholds.Pause.MinFramesForShortPause                        = round(Thresholds.Pause.MinTimeForShortPause                              * FrameRate);        
Thresholds.Pause.MinFramesForLongPause                         = round(Thresholds.Pause.MinTimeForLongPause                               * FrameRate);                                                                              
Thresholds.Pause.MinFramesForVeryLongPause                     = round(Thresholds.Pause.MinTimeForVeryLongPause                           * FrameRate);                                                                              
Thresholds.StraightRun.MinPixelPerFrame_ForRunDefinition       = Thresholds.StraightRun.MinVelocityAmplitude_ForRunDefinition       / FrameRate * PixelSize;    
Thresholds.StraightRun.MinPixelPerFrame_ForInermediateVelocity = Thresholds.StraightRun.MinVelocityAmplitude_ForInermediateVelocity / FrameRate * PixelSize;    
Thresholds.StraightRun.MinFramesForShortRun                    = round(Thresholds.StraightRun.MinTimeForShortRun                          * FrameRate);      % Used for Out of bounds extension  
Thresholds.StraightRun.AllowedExtensionFrameInterval           = round(Thresholds.StraightRun.AllowedExtensionTime                        * FrameRate);      
Thresholds.StraightRun.AllowedPausesFrameInterval              = round(Thresholds.StraightRun.AllowedPauseTimeInExtendedStraightRun       * FrameRate);      
Thresholds.StraightRun.AllowedPausesFrameInterval_VeryLong     = round(Thresholds.StraightRun.AllowedPauseTimeInLongExtendedStraightRun   * FrameRate);      
Thresholds.Omega.MinFramesForOmegaState_EccentricityCriterion  = round(Thresholds.Omega.MinTimeForOmegaState_EccentricityCriterion        * FrameRate);        
Thresholds.Omega.MinFramesForOmegaState_PerimeterCriterion     = round(Thresholds.Omega.MinTimeForOmegaState_PerimeterCriterion           * FrameRate);        
Thresholds.Omega.MinFramesForOmegaState                        = min([Thresholds.Omega.MinFramesForOmegaState_EccentricityCriterion  Thresholds.Omega.MinFramesForOmegaState_PerimeterCriterion]);         
Thresholds.Pirouette.MaxFramesBetweenReverseAndOmega           = round(Thresholds.Pirouette.MaxTimeIntervalBetweenReverseAndOmega         * FrameRate); 
Thresholds.Pirouette.MaxFramesBetweenLastReverse_And_Forward = round(Thresholds.Pirouette.MaxTimeIntervalBetweenLastReverse_And_Forward * FrameRate); 
Thresholds.Pirouette.MaxPauseFramesWithinPirouetteReversal     = round(Thresholds.Pirouette.MaxPauseTimeIntervalWithinPirouetteReversal   * FrameRate); 
Thresholds.Curve.MinFrames_ForCurveDefinition                  = round(Thresholds.Curve.MinTime_ForCurveDefinition                        * FrameRate); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Calculate Reverse versus Forward flags based on the worm's autofluorescence or gray pattern    
% forward_NumOfFrames  =  10;  % Good for 30Hz
forward_NumOfFrames    =  round(FrameRate/3);  

if (NumOfFrames >= forward_NumOfFrames+1) && (NumOfFrames >= Thresholds.MinimumTimeForBehaviorDefinition * FrameRate)  % The track is not too short for defining direction and for defining behavioral states.
    DirectionDuringRunsIsDefined = true;
    
    Distances.HeadToCentroid          =  sqrt((Head(:,1)-Path(:,1)).^2 + (Head(:,2)-Path(:,2)).^2);
    Distances.HeadToNextFrameCentroid =  sqrt((Head(1:end-forward_NumOfFrames,1)-Path(1+forward_NumOfFrames:end,1)).^2 + (Head(1:end-forward_NumOfFrames,2)-Path(1+forward_NumOfFrames:end,2)).^2);
    Distances.TailToCentroid          =  sqrt((Tail(:,1)-Path(:,1)).^2 + (Tail(:,2)-Path(:,2)).^2);
    Distances.TailToNextFrameCentroid =  sqrt((Tail(1:end-forward_NumOfFrames,1)-Path(1+forward_NumOfFrames:end,1)).^2 + (Tail(1:end-forward_NumOfFrames,2)-Path(1+forward_NumOfFrames:end,2)).^2);

    Distances.HeadToNextFrameCentroid(end:length(Frames))=NaN;
    Distances.TailToNextFrameCentroid(end:length(Frames))=NaN;

    ForwardByHead        = (Distances.HeadToNextFrameCentroid - Distances.HeadToCentroid)< 0;
    ForwardByTail        = (Distances.TailToNextFrameCentroid - Distances.TailToCentroid)> 0;  % (+forward_NumOfFrames*4)
    ForwardByHeadAndTail = ForwardByHead & ForwardByTail;
    ReverseByHead        = (Distances.HeadToNextFrameCentroid - Distances.HeadToCentroid)> 0;
    ReverseByTail        = (Distances.TailToNextFrameCentroid - Distances.TailToCentroid)< 0;  % (+forward_NumOfFrames*4)
    ReverseByHeadAndTail = ReverseByHead & ReverseByTail;

    if PlotSummaryFigures
        figure('name','Forward definition by head vs. Forward definition by tail'); 
        plot(ForwardByHead,'b.'); hold on;  plot(ForwardByTail,'bo'); 
        plot(0.995*(VelocityAmplitude>Thresholds.StraightRun.MinPixelPerFrame_ForRunDefinition ),'gx');
        ylim([0.98 1.01]); % xlim([0 202]);
        legend('ForwardByHead','ForwardByTail','High Velocity');
        figure('name','Forward and reverse definition by head and tail'); 
        plot(ForwardByHeadAndTail,'b.'); hold on;  plot(ReverseByHeadAndTail,'ro'); 
        plot(0.995*(VelocityAmplitude>Thresholds.StraightRun.MinPixelPerFrame_ForRunDefinition ),'gx');
        ylim([0.98 1.01]); % xlim([0 202]); 
        legend('ForwardByHeadAndTail','ReverseByHeadAndTail','High Velocity');
    end
    
else    
    DirectionDuringRunsIsDefined = false;
    FramesOutOfBounds            = true(size(Frames));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Find candidate sharp turns
if DirectionDuringRunsIsDefined
%     [locs_all, pks_all] = FindPeaks_Allow3FramesIntegration(ChangeInVelocityAngle, Thresholds.Turns.MinVelocityAngleChange, PlotSummaryFigures); 
    [locs_all, pks_all] = FindPeaks_Allow9FramesIntegration(ChangeInVelocityAngle, Thresholds.Turns.MinVelocityAngleChange, PlotSummaryFigures);   % Good for 30Hz
    PossibleTurns       = locs_all;

    FramesOutOfBoundsAndTheirNeighboringFrames  = find(  FramesOutOfBounds + [false FramesOutOfBounds(1:end-1)] + [FramesOutOfBounds(2:end) false]  );
    [PossibleTurns_NoNaNs,indices]              = setdiff(PossibleTurns,FramesOutOfBoundsAndTheirNeighboringFrames) ;          
    pks_all_NoNaNs                              = pks_all(indices);
    pks_all_NoNaNs_FullVec                      = single(zeros(1,length(Frames))*NaN);
    pks_all_NoNaNs_FullVec(PossibleTurns_NoNaNs)= pks_all_NoNaNs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Plot raw locomotion data
if PlotConfiguration == 1
    f=figure('name',['Analysis of track ',num2str(tr_num)]); 
    subplot(6,1,1); plot(Frames, ChangeInVelocityAngle); hold on; 
                    plot(Frames, SMOOTHED_ChangeInVelocityAngle_30deg,'r.-'); ylabel('Acc Angle'); 
                    Add_Xlines(gca, [-180 -120 -60 0 60 120 180] , 1); 
    subplot(6,1,2); plot(Frames, JerkAmplitude,'b.-');                      ylabel('Jerk amplitude');
    subplot(6,1,3); plot(Frames, VelocityAngle,'r.-');                      ylabel('Vel Angle');
                    Add_Xlines(gca, [-180 -120 -60 0 60 120 180] , 1) 
    subplot(6,1,4); plot(Frames, ChangeInVelocityAmplitude,'b.-');          ylabel('Acc amplitude');
    subplot(6,1,5); plot(Frames, VelocityAmplitude,'b.-');                  ylabel('Vel amplitude');
    subplot(6,1,6); plot(Frames, NormalizedSize,'b.-');                     ylabel('Normalized Size');
    subplot(6,1,6); [AX,H1,H2] = plotyy(Frames, Eccentricity,Frames, NormalizedSize);             
    set(get(AX(1),'Ylabel'),'String','Eccentricity','color','k');    set(AX(1),'ylim',[floor(nanmin(Eccentricity)*10) ceil(nanmax(Eccentricity)*10)]/10,'Ycolor','k');     set(H1,'marker','.','color','k');
    set(get(AX(2),'Ylabel'),'String','Normalized Size','color','g'); set(AX(2),'ylim',[floor(nanmin(NormalizedSize)*10-0.5) ceil(nanmax(NormalizedSize)*10)]/10,'Ycolor','g','ytick',[0.8 1 1.2]); set(H2,'marker','.','color','g');
    children_h = get(gcf,'children'); for h=children_h; set(h,'fontsize',6); end

elseif PlotConfiguration == 2
    f=figure('name',['Analysis of track ',num2str(tr_num)]); 
    subplot(4,1,1); plot(Frames, ChangeInVelocityAngle,'k.-'); hold on;  ylabel(['Acceleration',char(10),'Angle']); 
                    Add_Xlines(gca, [-180 -120 -60 0 60 120 180] , 1); 
    subplot(4,1,2); plot(Frames, VelocityContinuousAngle,'k.-');            ylabel(['Velocity',char(10),'Continuous Angle']);
                    set(gca,'ytick',[floor(min(VelocityContinuousAngle)/60)*60 : 60: ceil(max(VelocityContinuousAngle)/60)*60],'ygrid','on')
    subplot(4,1,3); plot(Frames, VelocityAmplitude,'k.-');                  ylabel(['Velocity',char(10),'Amplitude']);
    subplot(4,1,4); 
        [AX,H1,H2] = plotyy(Frames, Eccentricity,Frames, NormalizedSize);             
        set(get(AX(1),'Ylabel'),'String','Eccentricity','color','k');    
        set(AX(1),'ylim',[floor(nanmin(Eccentricity)*10)       ceil(nanmax(Eccentricity)*10)]/10,  'Ycolor','k');                     set(H1,'marker','.','color','k');
        set(get(AX(2),'Ylabel'),'String',['Normalized',char(10),'Size'],'color','g'); 
        set(AX(2),'ylim',[floor(nanmin(NormalizedSize)*10-0.5) ceil(nanmax(NormalizedSize)*10)]/10,'Ycolor','g','ytick',[0.8 1 1.2]); set(H2,'marker','.','color','g');
end    

IncludeSegmentsOfOneFrame = true;

%%  Behavior segmentation Strategy: 
if DirectionDuringRunsIsDefined
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %% 'Low Level Behavior Code'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % A. Define High and intermediate velocity frames.
    %    Define Out-of-bounds frames
    %
    % B. Define straight-runs and extend each segments If necessary.
    %    Define the direction and curvature of each run and sharp turns.
    %
    % C. Define Pause frames, short pause states and long pause states. Some of which will be considered as omegas.   
    %
    % D. Define omegas, some of the pauses may be omegas
    %    NOTE!! The shape criteria do not require any velocity requirement. In each segment, some of the frames may be considered a forward motion later on.     

    % D. Concatinate straight runs that are: proximal, without pause/omega/(sharp turns), same direction.   
    %    Correct the omegas and pauses states acordingly. 

    % E. Define curving frames within the straight runs.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% A. Define High and intermediate velocity frames. Define Out-of-bounds frames  
    % High velocity frames
    HighVelocity                                = VelocityAmplitude >= Thresholds.StraightRun.MinPixelPerFrame_ForRunDefinition;
    IntermediateVelocity                        = (VelocityAmplitude >= Thresholds.StraightRun.MinPixelPerFrame_ForInermediateVelocity)&(~HighVelocity);

    % Out of Bounds Segments and their extension
    Segments.OutOfBounds                                     = FindSegments      (find(FramesOutOfBounds),  IncludeSegmentsOfOneFrame);
    Logical_Vectors.OutOfBounds                              = false(1,length(Frames));
    if ~isempty(Segments.OutOfBounds)    
        Segments.OutOfBounds                                     = ExtendOutOfBound_OnlyBasedOnFrameInterval (Segments.OutOfBounds, NumOfFrames, Thresholds.StraightRun.MinFramesForShortRun);
        Indices_Vectors.OutOfBounds                              = [Segments.OutOfBounds.FramesIndices];
    else
        Indices_Vectors.OutOfBounds                              = [];
    end
    Logical_Vectors.OutOfBounds(Indices_Vectors.OutOfBounds) = true;
    FramesOutOfBounds                                        = Logical_Vectors.OutOfBounds;

    %% B. Define straight-runs and extend each segments if necessary.
    %     Find high velocity segments.
    %     Define the direction of each run. 
    %     Extend each straight run and evaluate real sharp turns as follows:
    %     Loop 1: For each high velocity segment:
    %       Check if there are adjacent SHARP TURNS   before/after   the segment. 
    %         If such a turn exist, add to this segment all frames that:  
    %            1. are after/before that turn to the straight run segment  ** AND **    
    %            2. have at least an intermediate velocity amplitude. 
    %         Is there an adjacent straight run of OPPOSITE DIRECTION before/after that turn?    
    %             If yes --> CONSIDER THIS AS A REAL SHARP TURN. others will be considered as jitters.    
    %     Loop 2: For each high velocity segment:
    %       Combine segments that are proximal and in the same direction. These small pauses will NOT be considered as a pause state.  

    % Straight runs (forward and reverse) and sharp turns.
    ForwardAndHighVelocity                      = ForwardByHeadAndTail;
    ForwardAndHighVelocity(~HighVelocity)       = false;
    ReverseAndHighVelocity                      = ReverseByHeadAndTail;
    ReverseAndHighVelocity(~HighVelocity)       = false;

    HighVelocity_And_InBound                    = HighVelocity;
    HighVelocity_And_InBound(FramesOutOfBounds) = false;
    HighVelocity_And_InBound                    = find(HighVelocity_And_InBound);  % logical representation --> frames indices

    % Straight runs: without extension
    Segments.HighVelocity                 = FindSegments        (HighVelocity_And_InBound,  IncludeSegmentsOfOneFrame); 
    [Segments.HighVelocity, DirectionFlag]= FindDirection       (Segments.HighVelocity ,    ForwardAndHighVelocity,   ReverseAndHighVelocity); 
    Segments.HighVelocity_not_extended    = Segments.HighVelocity;
    
    Logical_Vectors_not_extended.Forward                            = false(1,length(Frames));
    Logical_Vectors_not_extended.Reverse                            = false(1,length(Frames));
    Logical_Vectors.Forward                                         = false(1,length(Frames));
    Logical_Vectors.Reverse                                         = false(1,length(Frames));
    Logical_Vectors_AllowJitters.Forward                            = false(1,length(Frames));
    Logical_Vectors_AllowJitters.Reverse                            = false(1,length(Frames));
    Logical_Vectors_Long_extension.Forward                          = false(1,length(Frames));
    Logical_Vectors_Long_extension.Reverse                          = false(1,length(Frames));
    Logical_Vectors_Long_extension_AllowJitters.Forward             = false(1,length(Frames));
    Logical_Vectors_Long_extension_AllowJitters.Reverse             = false(1,length(Frames));
    Logical_Vectors.Forward                                         = false(1,length(Frames));
    Logical_Vectors.Reverse                                         = false(1,length(Frames));
    Logical_Vectors.Curving                                         = false(1,length(Frames));
    Logical_Vectors.CurvingForward                                  = false(1,length(Frames));
    Logical_Vectors.CurvingReverse                                  = false(1,length(Frames));
    Logical_Vectors.Forward_NoCurving                               = false(1,length(Frames));
    
    if ~isempty(Segments.HighVelocity) 
        DirectionOfSegment                                                          = [Segments.HighVelocity_not_extended.Direction];
        ForwardSegments                                                             = (DirectionOfSegment==0)|(DirectionOfSegment==1);
        Indices_Vectors_not_extended.Forward                                        = [Segments.HighVelocity_not_extended(ForwardSegments).FramesIndices];
        Indices_Vectors_not_extended.Reverse                                        = [Segments.HighVelocity_not_extended(~ForwardSegments).FramesIndices];
        Logical_Vectors_not_extended.Forward(Indices_Vectors_not_extended.Forward)  = true;
        Logical_Vectors_not_extended.Reverse(Indices_Vectors_not_extended.Reverse)  = true;


        % Straight runs: with extension. Don't Allow Jitters within a straight run   
        AllowJittersWithinStraightRuns = false;   
        [Segments.StraightRuns, ~, SharpTurns]= ExtendStraightRuns  (Segments.HighVelocity ,  DirectionFlag, FramesOutOfBounds, PossibleTurns_NoNaNs, ...
                                                                                              IntermediateVelocity, Thresholds.StraightRun,  AllowJittersWithinStraightRuns); 
        Segments.StraightRuns_Extended_NoJitters                = Segments.StraightRuns;
        Indices_Vectors.SharpTurns                              = find(SharpTurns);
        Logical_Vectors.SharpTurns                              = SharpTurns;
        SharpTurnsValues                                        = pks_all_NoNaNs_FullVec(Logical_Vectors.SharpTurns);

        DirectionOfSegment                               = [Segments.StraightRuns.Direction];
        ForwardSegments                                  = (DirectionOfSegment==0)|(DirectionOfSegment==1);
        Indices_Vectors.Forward                          = [Segments.StraightRuns(ForwardSegments).FramesIndices];
        Indices_Vectors.Reverse                          = [Segments.StraightRuns(~ForwardSegments).FramesIndices];
        Logical_Vectors.Forward(Indices_Vectors.Forward) = true;
        Logical_Vectors.Reverse(Indices_Vectors.Reverse) = true;

        % Straight runs: with extension (short). ALLOW Jitters within a straight run   
        % recalculate when allowing concatination of straight segments with a very
        % long pause time in between
        AllowJittersWithinStraightRuns = true;   
        Segments.StraightRuns = ExtendStraightRuns  (Segments.HighVelocity ,  DirectionFlag, FramesOutOfBounds, PossibleTurns_NoNaNs, ...
                                                                                                              IntermediateVelocity, Thresholds.StraightRun, AllowJittersWithinStraightRuns); 
        Segments.StraightRuns_Extended_AllowJitters                   = Segments.StraightRuns;
        DirectionOfSegment                                            = [Segments.StraightRuns.Direction];
        ForwardSegments                                               = (DirectionOfSegment==0)|(DirectionOfSegment==1);
        Indices_Vectors_AllowJitters.Forward                          = [Segments.StraightRuns(ForwardSegments).FramesIndices];
        Indices_Vectors_AllowJitters.Reverse                          = [Segments.StraightRuns(~ForwardSegments).FramesIndices];
        Logical_Vectors_AllowJitters.Forward(Indices_Vectors_AllowJitters.Forward) = true;
        Logical_Vectors_AllowJitters.Reverse(Indices_Vectors_AllowJitters.Reverse) = true;

        % Straight runs: with LONG extension. Don't Allow Jitters within a straight run   
        % recalculate when allowing concatination of straight segments with a very
        % long pause time in between
        AllowJittersWithinStraightRuns = false;   
        Thresholds.StraightRun.AllowedPausesFrameInterval = Thresholds.StraightRun.AllowedPausesFrameInterval_VeryLong ;
        [Segments.StraightRuns, DirectionFlag_StraightRuns_SuperExtended]= ExtendStraightRuns  (Segments.HighVelocity ,  DirectionFlag, FramesOutOfBounds, PossibleTurns_NoNaNs, ...
                                                                                                          IntermediateVelocity, Thresholds.StraightRun, AllowJittersWithinStraightRuns); 

        Segments.StraightRuns_Long_extension_NoJitters                  = Segments.StraightRuns;
        DirectionOfSegment                                              = [Segments.StraightRuns.Direction];
        ForwardSegments                                                 = (DirectionOfSegment==0)|(DirectionOfSegment==1);
        Indices_Vectors_Long_extension.Forward                          = [Segments.StraightRuns(ForwardSegments).FramesIndices];
        Indices_Vectors_Long_extension.Reverse                          = [Segments.StraightRuns(~ForwardSegments).FramesIndices];
        Logical_Vectors_Long_extension.Forward(Indices_Vectors_Long_extension.Forward) = true;
        Logical_Vectors_Long_extension.Reverse(Indices_Vectors_Long_extension.Reverse) = true;

        % Straight runs: with LONG extension. ALLOW JITTERS within a straight run   
        % recalculate when allowing concatination of straight segments with a very long pause time between them    
        AllowJittersWithinStraightRuns = true;   
        Thresholds.StraightRun.AllowedPausesFrameInterval = Thresholds.StraightRun.AllowedPausesFrameInterval_VeryLong ;
        if Thresholds.StraightRun.AllowedExtensionFrameInterval>=1
            Segments.StraightRuns = ExtendStraightRuns  (Segments.HighVelocity ,  DirectionFlag, FramesOutOfBounds, PossibleTurns_NoNaNs, ...
                                                                                                                  IntermediateVelocity, Thresholds.StraightRun, AllowJittersWithinStraightRuns); 
        end
        Segments.StraightRuns = Add_CurvingAndSpeed_Stats (Segments.StraightRuns, Thresholds, VelocityAngle, ChangeInVelocityAngle, VelocityAmplitude, ChangeInVelocityAmplitude, ...
                                                                                find(SharpTurns), SharpTurnsValues);

        Segments.StraightRuns_Long_extension_AllowJitters               = Segments.StraightRuns;
        DirectionOfSegment                                              = [Segments.StraightRuns.Direction];
        ForwardSegments                                                 = (DirectionOfSegment==0)|(DirectionOfSegment==1);
        Segments.Forward_Long_extension_AllowJitters                    = Segments.StraightRuns(ForwardSegments);
        Segments.Reverse_Long_extension_AllowJitters                    = Segments.StraightRuns(~ForwardSegments);
        Indices_Vectors_Long_extension_AllowJitters.Forward             = [Segments.StraightRuns(ForwardSegments).FramesIndices];
        Indices_Vectors_Long_extension_AllowJitters.Reverse             = [Segments.StraightRuns(~ForwardSegments).FramesIndices];
        Logical_Vectors_Long_extension_AllowJitters.Forward(Indices_Vectors_Long_extension_AllowJitters.Forward) = true;
        Logical_Vectors_Long_extension_AllowJitters.Reverse(Indices_Vectors_Long_extension_AllowJitters.Reverse) = true;

        %%%% CHOOSE DEFINITION OF STRAIGHT RUNS %%%% It is best to use 'StraightRuns_Extended_AllowJitters' for the definition. 
        % Only this one will be used for the low level forward and reverse runs
        % Segments.StraightRuns = Segments.StraightRuns_Long_extension_AllowJitters;        % This definition uses a too long extension parameters     
        Segments.StraightRuns = Segments.StraightRuns_Extended_AllowJitters;
        Segments.StraightRuns = Add_CurvingAndSpeed_Stats (Segments.StraightRuns, Thresholds, VelocityAngle, ChangeInVelocityAngle, VelocityAmplitude, ChangeInVelocityAmplitude, ...
                                                                                find(SharpTurns), SharpTurnsValues);
        StraightRunFrames     = [Segments.StraightRuns.FramesIndices];

        % Define forward and reverse frames using 'Segments.StraightRuns'
        DirectionOfSegment    = [Segments.StraightRuns.Direction];
        ForwardSegments       = (DirectionOfSegment==0)|(DirectionOfSegment==1);
        Segments.Forward      = Segments.StraightRuns(ForwardSegments);
        Segments.Reverse      = Segments.StraightRuns(~ForwardSegments);
        Indices_Vectors.Forward                          = [Segments.Forward.FramesIndices];
        Indices_Vectors.Reverse                          = [Segments.Reverse.FramesIndices];
        Logical_Vectors.Forward(Indices_Vectors.Forward) = true;
        Logical_Vectors.Reverse(Indices_Vectors.Reverse) = true;

        % Define forward curving frames using 'Segments.StraightRuns'
        CurvingFrames                                                  = [Segments.StraightRuns.CurvingFrames];
        Indices_Vectors.Curving                                        = CurvingFrames;
        Indices_Vectors.CurvingForward                                 = intersect(Indices_Vectors.Forward, CurvingFrames);  % This is a subset of the Forward indices !!!!!!!!!!!
        Indices_Vectors.CurvingReverse                                 = intersect(Indices_Vectors.Reverse, CurvingFrames);  % This is a subset of the Reverse indices !!!!!!!!!!!
        Logical_Vectors.Curving(Indices_Vectors.Curving)               = true;
        Logical_Vectors.CurvingForward(Indices_Vectors.CurvingForward) = true;
        Logical_Vectors.CurvingReverse(Indices_Vectors.CurvingReverse) = true;

        Logical_Vectors.Forward_NoCurving                                 = Logical_Vectors.Forward;
        Logical_Vectors.Forward_NoCurving(Indices_Vectors.CurvingForward) = false;

    else        
        SharpTurns                                              = false(1,length(Frames));        
        SharpTurnsValues                                        = [];
        
        Segments.StraightRuns_Extended_NoJitters                = []; 
        Segments.StraightRuns_Extended_AllowJitters             = []; 
        Segments.StraightRuns_Long_extension_NoJitters          = []; 
        Segments.StraightRuns_Long_extension_AllowJitters       = [];         
        Segments.StraightRuns                                   = []; 
        Segments.Forward_Long_extension_AllowJitters            = []; 
        Segments.Reverse_Long_extension_AllowJitters            = []; 
        Segments.Forward                                        = [];
        Segments.Reverse                                        = [];        
        StraightRunFrames                                       = [];

    end             

    %%%%%%%%        Define sharp turns  %%%%%%%%%%%%%%%
    SharpTurns(StraightRunFrames) = false;   % The rare cases of sharp turns during a forward movement is dealt within the 'Add_CurvingAndSpeed_Stats' function 

    Segments.SharpTurns              = FindSegments        (find(SharpTurns),            IncludeSegmentsOfOneFrame); 
    Logical_Vectors.SharpTurns       = false(1,length(Frames));
    if ~isempty(Segments.SharpTurns)
        Indices_Vectors.SharpTurns   = [Segments.SharpTurns.FramesIndices];
    else
        Indices_Vectors.SharpTurns   = [];
    end    
    Logical_Vectors.SharpTurns(Indices_Vectors.SharpTurns)  = true;

    for turn_ind = 1:length(SharpTurnsValues)
        Segments.SharpTurns(turn_ind).ChangeInAngle = SharpTurnsValues(turn_ind);
    end        
    
    %% C. Pause frames, short pause states and long pause states. Some of which will be considered as omegas  
    LowVelocity                          = true(1,length(Frames));
    LowVelocity(StraightRunFrames)       = false;
    LowVelocity(SharpTurns)              = false;
    LowVelocity(FramesOutOfBounds)       = false;
    Segments.LowVelocity                 = FindSegments        (find(LowVelocity),  IncludeSegmentsOfOneFrame); 
    Segments.LowVelocity_LongSegments    = DeleteShortSegments (Segments.LowVelocity,   Thresholds.Pause.MinFramesForLongPause); 

    %% D. Initial definition of omega state based on the worm shape and the proximity of a long pause state. 
    %  The 'OmegaShape' criteria do not require any velocity requirement. 
    %  NOTE!!!! Some Of the 'OmegaState' frames will be corrected to be part of straight runs (the shape starts to change when the worm is still at high speed).     
    %
    % THE OMEGA SHAPE CRITERIA 
    % Omega state can occur if one of 3 criteria accur:
    %   1. Omega by eccentricity:           low eccentricity for a long time
    %   2. Omega by perimeter:              low perimeter but high perimeter/size for a long time (longer than eccentricity criteria). This happens when worms squeeze between posts. 
    %   3. Omega by pause and eccentricity: low eccentricity for a short time but coupled with a long pause. 

    % Find frames with specific omega features
    OmegaSegmentedByPerimeter     = ( (NormalizedPerimeterLength  < Thresholds.Omega.MaxNormPerimeterForOmega).*(NormalizedPerimeterOverSize  > Thresholds.Omega.MinNormPerimeterOverSizeForOmega) );
    OmegaSegmentedByPerimeter(FramesOutOfBounds)    = false;
    OmegaSegmentedByPerimeter                       = find(OmegaSegmentedByPerimeter);
    OmegaSegmentedByEccentricity                    = (    Eccentricity    < Thresholds.Omega.MaxEccentricityForOmega);
    OmegaSegmentedByEccentricity(FramesOutOfBounds) = false;   
    OmegaSegmentedByEccentricity                    = find(OmegaSegmentedByEccentricity);

    % Omega Shape segments
    Segments.OmegaByEccentricity    = FindSegments        (OmegaSegmentedByEccentricity,   IncludeSegmentsOfOneFrame);
    Segments.OmegaByPerimeter       = FindSegments        (OmegaSegmentedByPerimeter,      IncludeSegmentsOfOneFrame);
    Segments.OmegaByEccentricity    = DeleteShortSegments (Segments.OmegaByEccentricity,   Thresholds.Omega.MinFramesForOmegaState_EccentricityCriterion );
    Segments.OmegaByPerimeter       = DeleteShortSegments (Segments.OmegaByPerimeter,      Thresholds.Omega.MinFramesForOmegaState_PerimeterCriterion    );

    % Find overlapping segments with both eccentricity and perimeter omega features     
    if (~isempty(Segments.OmegaByPerimeter))&&(~isempty(Segments.OmegaByEccentricity))
        % IMPORTANT NOTE!! The 'OR' operation is finding all fragements that pass either one of the two field thresholds      
        % BUT ONLY IN SEGMENTS WHERE SOME FRAME OVERLAP EXISTS !!!!       
        Segments    = CombineSegments     (Segments , 'OmegaByEccentricity', 'OmegaByPerimeter', {'OR','AND'});            % optional additional input: {'name1','name2'}
    end

    % ALL OMEGA FRAMES THAT ARE PART OF A LONG ENOUGH OMEGA SHAPE SEGMENTS
    if ~isempty(Segments.OmegaByEccentricity)
        Omega_Frames                    =  [Segments.OmegaByEccentricity.FramesIndices];
        if ~isempty(Segments.OmegaByPerimeter)
            Omega_Frames = union([Segments.OmegaByPerimeter.FramesIndices],     [Segments.OmegaByEccentricity.FramesIndices]);         % ALL OMEGA FRAMES THAT HAVE LONG ENOUGH SEGMENTS
        end
    else
        Omega_Frames = [];
    end
    Segments.OmegaShape   = FindSegments        (Omega_Frames,                   IncludeSegmentsOfOneFrame);
    Segments.OmegaShape   = DeleteShortSegments (Segments.OmegaShape ,           Thresholds.Omega.MinFramesForOmegaState );
    Segments              = CombineSegments_OR_if_minimal_overlap_AND_otherwise (Segments , 'OmegaShape', 'LowVelocity_LongSegments', Thresholds.Omega.MinOverlapWithPauseForOmegaExtension);   
    if ~isempty(Segments.OmegaShape_ExtendedBy_LowVelocity_LongSegments)
        OmegaState            = intersect(find(LowVelocity), [Segments.OmegaShape_ExtendedBy_LowVelocity_LongSegments.FramesIndices]);  % This 'cleans' the omega shape within the straight run
    else
        OmegaState = [];
    end
    Segments.OmegaState   = FindSegments  (OmegaState,  IncludeSegmentsOfOneFrame);
    Segments.Omega        = Segments.OmegaState;

    % Pauses are the low velocity frames that are not considered an omega state    
    Pause                     = LowVelocity;
    Pause(OmegaState)         = false;
    Segments.Pause            = FindSegments  (find(Pause),  IncludeSegmentsOfOneFrame);
    Segments.VeryLongPause    = DeleteShortSegments (Segments.Pause,  Thresholds.Pause.MinFramesForVeryLongPause);
    Logical_VeryLongPauses    = false(1, length(Frames));
    if ~isempty(Segments.VeryLongPause)
        Logical_VeryLongPauses([Segments.VeryLongPause.FramesIndices]) = true;
    end

    Logical_Vectors.Omega                           = false(1,length(Frames));
    Logical_Vectors.Pause                           = false(1,length(Frames));

    % Assign values
    if ~isempty(Segments.OmegaState)
        Indices_Vectors.Omega                           = [Segments.OmegaState.FramesIndices];
        Logical_Vectors.Omega(Indices_Vectors.Omega)    = true;
    end
    if ~isempty(Segments.Pause)           
        Indices_Vectors.Pause                           = [Segments.Pause.FramesIndices];           
        Logical_Vectors.Pause(Indices_Vectors.Pause)    = true;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%   End of low level behavior segmentation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Assign values to segment low level behavior structure   
    Segments_LowLevelBehavior.Forward     = Segments.Forward;    % INCLUDING CURVING FRAMES !!!!
    Segments_LowLevelBehavior.Reverse     = Segments.Reverse;    % INCLUDING CURVING FRAMES !!!!
    Segments_LowLevelBehavior.Pause       = Segments.Pause;      % INCLUDING Jitters !!!!
    Segments_LowLevelBehavior.Omega       = Segments.Omega;
    Segments_LowLevelBehavior.OutOfBounds = Segments.OutOfBounds;
    Segments_LowLevelBehavior.SharpTurns  = Segments.SharpTurns;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %% High level bahavior code. This includes Pirouettes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    [Segments_HighLevelBehavior, Logical_Vectors_HighLevel] = CalculateHighLevelBehavior  (Segments_LowLevelBehavior, Logical_Vectors, Segments.Reverse_Long_extension_AllowJitters, ...
                                                                                                Thresholds, Logical_VeryLongPauses);
    % Add curving code:                                                                                        
    Logical_Vectors_HighLevel.Forward_NoCurving   = Logical_Vectors_HighLevel.Forward & (~Logical_Vectors.CurvingForward);                          
    Logical_Vectors_HighLevel.CurvingForward      = Logical_Vectors_HighLevel.Forward &   Logical_Vectors.CurvingForward; 
    
else      % If the Track length is too short for behavior segmentation (probably due to collisions?) then behaviour will be defined as OutOfBounds    
    
    Logical_Vectors.OutOfBounds        = true(1,length(Frames));
    Logical_Vectors.Omega              = false(1,length(Frames));
    Logical_Vectors.Pause              = false(1,length(Frames));
    Logical_Vectors.Curving            = false(1,length(Frames));
    Logical_Vectors.CurvingForward     = false(1,length(Frames));
    Logical_Vectors.CurvingReverse     = false(1,length(Frames));
    Logical_Vectors.Forward_NoCurving  = false(1,length(Frames));
    Logical_Vectors.Forward            = false(1,length(Frames));
    Logical_Vectors.Reverse            = false(1,length(Frames));
    Logical_Vectors.SharpTurns         = false(1,length(Frames));
    
    Logical_Vectors_HighLevel.OutOfBounds                   = true(1,length(Frames));
    Logical_Vectors_HighLevel.Pirouette                     = false(1,length(Frames));
    Logical_Vectors_HighLevel.Pirouette_ExactInitiationTime = false(1,length(Frames));
    Logical_Vectors_HighLevel.Pirouette_NotAfterOutOfBound  = false(1,length(Frames));
    Logical_Vectors_HighLevel.Pirouette_UntilPause          = false(1,length(Frames));
    Logical_Vectors_HighLevel.Pirouette_InitialReversal     = false(1,length(Frames));
    Logical_Vectors_HighLevel.Pirouette_AfterPause          = false(1,length(Frames));
    Logical_Vectors_HighLevel.Pirouette_AfterReversal       = false(1,length(Frames));
    Logical_Vectors_HighLevel.Forward                       = false(1,length(Frames));
    Logical_Vectors_HighLevel.Reversal                      = false(1,length(Frames));
    Logical_Vectors_HighLevel.OmegaWithoutPirouette         = false(1,length(Frames));
    Logical_Vectors_HighLevel.SharpTurns                    = false(1,length(Frames));
    Logical_Vectors_HighLevel.Pause                         = false(1,length(Frames));
    Logical_Vectors_HighLevel.Forward_NoCurving             = false(1,length(Frames));
    Logical_Vectors_HighLevel.CurvingForward                = false(1,length(Frames));
    
    Segments_HighLevelBehavior = [];
    Segments_LowLevelBehavior  = [];
    Segments                   = [];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Generate behavior code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
BehaviorCode = GenerateBehaviorCode (Logical_Vectors, Logical_Vectors_HighLevel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if PlotSummaryFigures && DirectionDuringRunsIsDefined
    SCSZ = get(0,'ScreenSize');

    figure('position',[50 floor(SCSZ(4)/2) SCSZ(3)-100 floor(SCSZ(4)/2)-100],'name','compare Straight Run strategy'); 
    plot(Logical_Vectors_not_extended.Forward,         'b.','markersize',6); hold on;  % Not extended
    plot(Logical_Vectors_not_extended.Reverse,         'r.','markersize',6); hold on;  % Not extended
    plot(Logical_Vectors_AllowJitters.Forward,         'bo','markersize',6); hold on;  % Not extended
    plot(Logical_Vectors_AllowJitters.Reverse,         'ro','markersize',6); hold on;  % Not extended
    plot(Logical_Vectors.Forward,                      'bo','markersize',10); hold on;  % Not extended
    plot(Logical_Vectors.CurvingForward,               'bo','markersize',10); hold on;  % Not extended
    plot(Logical_Vectors.Reverse,                      'ro','markersize',10); hold on;  % Not extended
    plot(Logical_Vectors.SharpTurns,'ks','markerfacecolor','k','markersize',10); hold on;
    plot(Logical_Vectors.Omega ,'gv','markerfacecolor','g'); hold on;
    plot(Logical_Vectors.Pause ,'y.'); hold on;
    plot(Logical_Vectors.OutOfBounds ,'ms','markerfacecolor','m'); hold on;
    legend('F (w/oJ)','R (w/oJ)','F (withJ)','R (withJ)','Forward','Curve','Reverse','SharpTurn','Omega','Pause','OOB');
    ylim([0.955 1.06])

    SCSZ = get(0,'ScreenSize');
    figure('position',[50 floor(SCSZ(4)/2) SCSZ(3)-100 floor(SCSZ(4)/2)-100],'name','Low and High Level Behavior, Frame index'); 
    plot(Logical_Vectors.Forward,        'b.','markersize',6); hold on;  % Not extended
    plot(Logical_Vectors.CurvingForward, 'b.','markersize',6); hold on;  % Not extended
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
    ylim([0.955 1.06]);
    xlabel('Frame index')
    legend('Forward','Curve','Reverse','SharpTurn','Omega','Pause','OOB', 'Pirouette','Reversal','Omega (noPir)','SharpTurns','Forward','Pause','OOB');

    figure('position',[50 floor(SCSZ(4)/2) SCSZ(3)-100 floor(SCSZ(4)/2)-100],'name','Low and High Level Behavior, Frames in movie'); 
    plot(Frames, Logical_Vectors.Forward,        'b.','markersize',6); hold on;  % Not extended
    plot(Frames, Logical_Vectors.CurvingForward, 'b.','markersize',6); hold on;  % Not extended
    plot(Frames, Logical_Vectors.Reverse,     'r.','markersize',6); hold on;  % Not extended
    plot(Frames, Logical_Vectors.SharpTurns,  'ks','markerfacecolor','k','markersize',10); hold on;
    plot(Frames, Logical_Vectors.Omega ,      'gv','markerfacecolor','g'); 
    plot(Frames, Logical_Vectors.Pause ,      'y.'); hold on;
    plot(Frames, Logical_Vectors.OutOfBounds ,'ms','markerfacecolor','m'); hold on;
    plot(Frames, 1.01*Logical_Vectors_HighLevel.Pirouette,              'rs','markerfacecolor','g');
    plot(Frames, 1.01*Logical_Vectors_HighLevel.Reversal,               'r.');
    plot(Frames, 1.01*Logical_Vectors_HighLevel.OmegaWithoutPirouette,  'gv','markerfacecolor','g');
    plot(Frames, 1.01*Logical_Vectors_HighLevel.SharpTurns,             'ks','markerfacecolor','k');
    plot(Frames, 1.01*Logical_Vectors_HighLevel.Forward,                'b.');
    plot(Frames, 1.01*Logical_Vectors_HighLevel.Pause,                  'y.');
    plot(Frames, 1.01*Logical_Vectors_HighLevel.OutOfBounds,            'ms','markerfacecolor','m'); hold on;
    ylim([0.955 1.06]);
    xlabel('Frame in movie')
    legend('Forward','Curve','Reverse','SharpTurn','Omega','Pause','OOB', 'Pirouette','Reversal','Omega (noPir)','SharpTurns','Forward','Pause','OOB');

    if ~isempty(Segments_HighLevelBehavior.Pirouette)
        figure('position',[50 floor(SCSZ(4)/2) SCSZ(3)-100 floor(SCSZ(4)/2)-100],'name','Pirouette angles'); 
        plot([Segments_HighLevelBehavior.Pirouette.FirstFrameIndices],[Segments_HighLevelBehavior.Pirouette.AngleBeforeReversal], 'k*'); hold on;
        plot([Segments_HighLevelBehavior.Pirouette.LastFrameIndices],[Segments_HighLevelBehavior.Pirouette.AngleAfterPirouette],'b*');
        ylim([-180 180]);
        legend('Angle Before Reversal','Angle After Reversal');
    end
end

return

%% Inline functions
function Segments = FindSegments (Indices, IncludeSegmentsOfOneFrame, ForceSegmentsOfOneFrame)
if isempty(Indices)
    Segments=[];
    return
end
Indices           = single(Indices);
LongSegmentsExist = ~isempty(find(diff(Indices)==1, 1));

if exist('ForceSegmentsOfOneFrame','var')
    if ForceSegmentsOfOneFrame
        LongSegmentsExist=false;
    end
end        
        
if LongSegmentsExist
    IsSegment  = diff(Indices)==1;
    if IsSegment(end)== true
        IsSegment = [IsSegment false]; % Pad the ending with zero if necessary
    end
    FirstIndex       = find(IsSegment,1,'first');      
    Segment_index    = 0;
    for i=1:5000 
        Segment_index                   = Segment_index+1;
        SegmentLength                   = find(IsSegment(FirstIndex:end)==0,1,'first');    

        Segments(Segment_index).FirstFrameIndices  = Indices(FirstIndex);
        Segments(Segment_index).FramesIndices      = Indices(FirstIndex:(FirstIndex+SegmentLength-1));
        Segments(Segment_index).MiddleFrameIndices = floor(Indices(FirstIndex)+ (SegmentLength-1)/2);
        Segments(Segment_index).Length             = SegmentLength;

        IsSegment(FirstIndex:(FirstIndex+SegmentLength-1))  = false; 
        FirstIndex = find(IsSegment,1,'first');      
        if isempty(FirstIndex)
            break
        end
    end
    if i==5000
        disp('aborting segmentation. Too many segments were found')
        return
    end
else
    Segments.FramesIndices = [];
    if length(Indices)==1
        IsSegment = 0;
    end
end

%% Add segments of one frame
if IncludeSegmentsOfOneFrame
    % Segments with one frame
    FramesIncludedInLongSegments = ismember(Indices, [Segments.FramesIndices]);
    segment_index_short            = 0;
    for frame_index = Indices(~FramesIncludedInLongSegments)
        segment_index_short = segment_index_short + 1;
        ShortSegments(segment_index_short).FirstFrameIndices  = frame_index;
        ShortSegments(segment_index_short).FramesIndices      = frame_index;
        ShortSegments(segment_index_short).MiddleFrameIndices = frame_index;
        ShortSegments(segment_index_short).Length             = 1;
    end
    if ~isempty(frame_index)                             % If at least one SHORT segment was found    
        if ~isempty(find(FramesIncludedInLongSegments,1))  % If at least one LONG  segment was found               
            % Merge short segments and long segments based in their timing    
            TemporaryMergedSegments = Segments;
            TemporaryMergedSegments((length(Segments)+1):(length(ShortSegments) + length(Segments))) = ShortSegments;

            [~,Order] = sort([TemporaryMergedSegments.FirstFrameIndices]);

            for segment_index_new = 1 : (length(TemporaryMergedSegments))
                NewSegments(segment_index_new) = TemporaryMergedSegments(Order(segment_index_new));               
            end
            Segments = NewSegments;  
        else
            Segments = ShortSegments;
        end
    end
end

return

function Segments = DeleteShortSegments (SegmentsIn , MinSegmentLength)
if isempty(SegmentsIn)
    Segments=[];
    return
end
IncludeSegment = [SegmentsIn.Length] >= MinSegmentLength;  % indices of segments that are long enough

if isempty(IncludeSegment)
    Segments = [];
else
    segment_index_new = 0;
    for segment_index_old = find(IncludeSegment)   
        segment_index_new           = segment_index_new + 1;
        Segments(segment_index_new) = SegmentsIn(segment_index_old);               
    end    
end
if ~segment_index_new
    Segments = [];
end
return

function [SegmentsHighVelocity, DirectionFlag] = FindDirection (SegmentsHighVelocity , Forward, Reverse)
% Example: Segments.HighVelocity  = FindDirection  (Segments.HighVelocity ,  ForwardAndHighVelocity,   ReverseAndHighVelocity, Settings); 
if isempty(SegmentsHighVelocity)
    DirectionFlag = NaN;
    return
end
SegmentsHighVelocity(1).Direction = NaN;
DirectionFlag                     = zeros(1,length(Forward));

for seg_ind = 1:length(SegmentsHighVelocity)
    FramesIndices = SegmentsHighVelocity(seg_ind).FramesIndices;
    Length        = SegmentsHighVelocity(seg_ind).Length;
    ForwardCount  = sum(Forward(FramesIndices));   % Frames that are definitely Forward
    ReverseCount  = sum(Reverse(FramesIndices));   % Frames that are definitely Reverse
    
    CorrectDirection = ForwardCount;
    WrongDirection   = ReverseCount;
    if ForwardCount == ReverseCount          % Unknown direction. most probably forward
        Direction        = 0;      
    elseif ForwardCount>ReverseCount         % Forward
        Direction        = 1; 
    else                                     % Reverse
        Direction        = -1; 
        CorrectDirection = ReverseCount;
        WrongDirection   = ForwardCount;
    end
    DirectionFlag(FramesIndices) = ones(size(FramesIndices))*Direction;
    
    SegmentsHighVelocity(seg_ind).Direction                       = single(Direction);
    SegmentsHighVelocity(seg_ind).ForwardCount                    = single(ForwardCount);
    SegmentsHighVelocity(seg_ind).ReverseCount                    = single(ReverseCount);
    SegmentsHighVelocity(seg_ind).CorrectVsWrongDirectionFlags    = single(CorrectDirection / WrongDirection);    
    SegmentsHighVelocity(seg_ind).FractionOfDirectionFlagedFrames = single(CorrectDirection / Length);        
end

return

function Segments_StraightRuns = Add_CurvingAndSpeed_Stats (Segments_StraightRunsIn , Thresholds, ...
                                        VelocityAngle, ChangeInVelocityAngle, VelocityAmplitude, ChangeInVelocityAmplitude, SharpTurnsFrames, SharpTurnsValues)

MinPixelPerFrame_ForRunDefinition    = Thresholds.StraightRun.MinPixelPerFrame_ForRunDefinition;
MinVelocityAngleChange               = Thresholds.Curve.MinVelocityAngleChange;
MinFrames_ForCurveDefinition         = Thresholds.Curve.MinFrames_ForCurveDefinition; 
                   
Segments_StraightRuns      = Segments_StraightRunsIn;   % This is either Forward or Reverse segments
HighVelocityFrames         = VelocityAmplitude > MinPixelPerFrame_ForRunDefinition; 

CurvingFrames                                                      = false(1,length(VelocityAngle));
CurvingFrames(abs(ChangeInVelocityAngle)>MinVelocityAngleChange)   = true;  
CurvingFrames(~HighVelocityFrames)                                 = false;
CurvingDirection                                                   = int8(zeros(1,length(VelocityAngle)));
CurvingDirection(CurvingFrames)                                    = int8(sign(ChangeInVelocityAngle(CurvingFrames)));

FastRuns_NoCurving                = HighVelocityFrames;
FastRuns_NoCurving(CurvingFrames) = false; 

NumberOfStraightRuns  = length(Segments_StraightRuns);

for seg_ind = 1:NumberOfStraightRuns
    FramesIndices                                = Segments_StraightRuns(seg_ind).FramesIndices;
    FirstFrame                                   = Segments_StraightRuns(seg_ind).FirstFrameIndices;
    LastFrame                                    = Segments_StraightRuns(seg_ind).FramesIndices(end);
    Length                                       = Segments_StraightRuns(seg_ind).Length;      
    
    % Curving
    CurvingFramesInSegment                       = CurvingFrames(FramesIndices);  
    CurrentCurvingFrames                         = FramesIndices(CurvingFramesInSegment);  
    CurvingValuesInSegment                       = ChangeInVelocityAngle(FramesIndices);
    CurvingDirectionInSegment                    = CurvingDirection(FramesIndices);
    NumberOfCurvingFrames                        = length(find(CurvingFramesInSegment)); 
    FractionOfCurvingFrames                      = NumberOfCurvingFrames/Length; 
    CurvingSpeedAverage                          = nanmean(abs(CurvingValuesInSegment(CurvingFramesInSegment)));    % Average ONLY DURING THE CURVING FRAMES !!!!!!
    TotalCurvature                               = nansum(abs(CurvingValuesInSegment(CurvingFramesInSegment)));    
    
    % Angle
    FirstNoCurvingFrame                          = FramesIndices(find(FastRuns_NoCurving(FramesIndices),1,'first'));
    Angle_FirstNoCurvingFrame                    = VelocityAngle(FirstNoCurvingFrame);
    Angle_RunDirection                           = Angle_FirstNoCurvingFrame;    
    Angle_InitialFrame                           = VelocityAngle(FirstFrame);
    Angle_FinalFrame                             = VelocityAngle(LastFrame);
    Angle_TotalChangeDuringRun                   = Angle_FinalFrame - Angle_InitialFrame;    

    % Angle_RunDirection correction (if needed) 
    % if all the frames in the segment are curving ones (direction is constantly changed), the run direction will be the taken as the direction after a short time, as defined by the user.  
    FramesForCheckingInitialCurvature            = (FirstFrame):(min([FirstFrame+MinFrames_ForCurveDefinition-1, LastFrame]));
    if isempty(Angle_RunDirection) 
        Angle_RunDirection = VelocityAngle(FramesForCheckingInitialCurvature(end));
    end
    if nansum(CurvingFrames(FramesForCheckingInitialCurvature))>0.5  % at least half of the initial frames are curving
        SegmentStartsWithCurving = true;
    else
        SegmentStartsWithCurving = false;
    end 
    
    % Speed 
    CurrentHighVelocityFrames      = FramesIndices(HighVelocityFrames(FramesIndices));
    Speed_Average_including_Pauses = nanmean(VelocityAmplitude(FramesIndices));
    Speed_Average                  = nanmean(VelocityAmplitude(CurrentHighVelocityFrames));
    Speed_Max                      = nanmax(VelocityAmplitude(CurrentHighVelocityFrames));
    Speed_Min                      = nanmin(VelocityAmplitude(CurrentHighVelocityFrames));
    Speed_AverageChange            = nanmean(ChangeInVelocityAmplitude(CurrentHighVelocityFrames));
    
    % Sharp turn?? very rare. This may accur at the device edges
    [FramesWithSharpTurns, ia]     = intersect(SharpTurnsFrames, CurrentHighVelocityFrames);
    CurrentSharpTurnsValues        = SharpTurnsValues(ia);
    
    % Assign to structure
    Segments_StraightRuns(seg_ind).Angle_InitialFrame              = Angle_InitialFrame;
    Segments_StraightRuns(seg_ind).Angle_FinalFrame                = Angle_FinalFrame;
    Segments_StraightRuns(seg_ind).Angle_RunDirection              = Angle_RunDirection;
    Segments_StraightRuns(seg_ind).Angle_TotalChangeDuringRun      = Angle_TotalChangeDuringRun;
    
    Segments_StraightRuns(seg_ind).CurvingFramesInSegment          = CurvingFramesInSegment;
    Segments_StraightRuns(seg_ind).CurvingFrames                   = CurrentCurvingFrames;
    Segments_StraightRuns(seg_ind).CurvingDirectionInSegment       = CurvingDirectionInSegment; 
    Segments_StraightRuns(seg_ind).FractionOfCurvingFrames         = FractionOfCurvingFrames;     
    Segments_StraightRuns(seg_ind).CurvingSpeedAverage             = CurvingSpeedAverage;   % Average ONLY DURING THE CURVING FRAMES !!!!!!
    Segments_StraightRuns(seg_ind).TotalCurvature                  = TotalCurvature;
    Segments_StraightRuns(seg_ind).TotalCurvatureOverSegmentLength = TotalCurvature/Length;
    Segments_StraightRuns(seg_ind).SegmentStartsWithCurving        = SegmentStartsWithCurving;
    
    Segments_StraightRuns(seg_ind).Speed_Average_including_Pauses  = Speed_Average_including_Pauses;
    Segments_StraightRuns(seg_ind).Speed_Average                   = Speed_Average;
    Segments_StraightRuns(seg_ind).Speed_Max                       = Speed_Max;
    Segments_StraightRuns(seg_ind).Speed_Min                       = Speed_Min;
    Segments_StraightRuns(seg_ind).Speed_AverageChange             = Speed_AverageChange;    
    
    Segments_StraightRuns(seg_ind).SharpTurnsFrames                = FramesWithSharpTurns;
    Segments_StraightRuns(seg_ind).SharpTurnsValues                = CurrentSharpTurnsValues;          
end

return

function [SegmentsHighVelocity, DirectionFlag, SharpTurns] = ExtendStraightRuns (SegmentsHighVelocity, DirectionFlag, FramesOutOfBounds, PossibleTurns_NoNaNs, ...
                                                                                  IntermediateVelocity, ThresholdsStraightRun, AllowJittersWithinStraightRuns)
if ~exist('AllowJittersWithinStraightRuns','var')
    AllowJittersWithinStraightRuns = false;
end

% Correct the DirectionFlag vector to consist of only 1,-1 or NaNs. All 0 values (uncertain runs) will be considered as forward runs (1) 
DirectionFlag_NoZerosWithNaNs                                   = zeros(size(DirectionFlag))*NaN;
StraightRunFrames                                               = [SegmentsHighVelocity.FramesIndices];
DirectionFlag_NoZerosWithNaNs(StraightRunFrames)                = DirectionFlag(StraightRunFrames);
DirectionFlag_NoZerosWithNaNs(DirectionFlag_NoZerosWithNaNs==0) = 1;   % All uncertain straight runs are considered forward runs
DirectionFlag                                                   = DirectionFlag_NoZerosWithNaNs;

% Example:
% [Segments.HighVelocity, DirectionFlag]= ExtendStraightRuns  (Segments.HighVelocity ,  DirectionFlag, FramesOutOfBounds, PossibleTurns_NoNaNs, Thresholds.StraightRun.AllowedExtensionFrameInterval); 
AllowedExtensionFrameInterval    =  ThresholdsStraightRun.AllowedExtensionFrameInterval;
AllowedPausesFrameInterval       =  ThresholdsStraightRun.AllowedPausesFrameInterval;

CurrentTrackLength                          = length(FramesOutOfBounds);
SegmentsHighVelocity(1).DirectionWasChanged = false;                    
SegmentsHighVelocity(1).DirectionWillChange = false;
SharpTurns                                  = false(1,CurrentTrackLength);

for seg_ind = 1:length(SegmentsHighVelocity)                               % Loop over all straight runs 
    % initialization
    SegmentsHighVelocity(seg_ind).DirectionWasChanged = false;
    SegmentsHighVelocity(seg_ind).DirectionWillChange = false;

    % Find relevant flanking regions
    FrameIndices              = SegmentsHighVelocity(seg_ind).FramesIndices;
    FramesAlreadyDefinedAsStraightRun                    = false(1,length(DirectionFlag));
    FramesAlreadyDefinedAsStraightRun(StraightRunFrames) = true;
    FramesAlreadyDefinedAsStraightRun(FrameIndices)      = false;  % Only the current segment frames are allowed
    
    FirstFrame                = FrameIndices(1);
    LastFrame                 = FrameIndices(end);                    
    CurrentDirection          = SegmentsHighVelocity(seg_ind).Direction;           
    if CurrentDirection==0
        DirectionUncertaintyFlag = true;
        CurrentDirection         = 1;   % CurrentDirection is now either 1 or -1    
        SegmentsHighVelocity(seg_ind).Direction = 1;
    else
        DirectionUncertaintyFlag = false;
    end
    SegmentsHighVelocity(seg_ind).DirectionUncertaintyFlag = DirectionUncertaintyFlag;
    
    RelevantIndices_BeforeRun = (FirstFrame-AllowedExtensionFrameInterval):(FirstFrame-1);
    RelevantIndices_AfterRun  = (LastFrame+1):(LastFrame+AllowedExtensionFrameInterval);        
    if seg_ind == 1
        RelevantIndices_BeforeRun = RelevantIndices_BeforeRun(RelevantIndices_BeforeRun>0);
    end
    if seg_ind == length(SegmentsHighVelocity)  
        RelevantIndices_AfterRun = RelevantIndices_AfterRun(RelevantIndices_AfterRun<=CurrentTrackLength);
    end       
    % Avoid overlapping Straight Run segments
    LastOverlappingIndex = find(ismember(RelevantIndices_BeforeRun,StraightRunFrames),1,'last');
    if ~isempty(LastOverlappingIndex)
        RelevantIndices_BeforeRun = RelevantIndices_BeforeRun((LastOverlappingIndex+1):end);
    end
    FirstOverlappingIndex = find(ismember(RelevantIndices_AfterRun,StraightRunFrames),1,'first');
    if ~isempty(FirstOverlappingIndex)
        RelevantIndices_AfterRun = RelevantIndices_AfterRun(1:(FirstOverlappingIndex-1));
    end
 
    % Skip if flanking regions consist of out of bounds frames or if no turns were found within the range 
    if sum(FramesOutOfBounds(RelevantIndices_BeforeRun))>0  % includes out of bounds frame
        RelevantIndices_BeforeRun = [];
    elseif ~sum(ismember(RelevantIndices_BeforeRun,PossibleTurns_NoNaNs))
        RelevantIndices_BeforeRun = [];
    end
    if sum(FramesOutOfBounds(RelevantIndices_AfterRun))>0  % includes out of bounds frame
        RelevantIndices_AfterRun = [];
    elseif ~sum(ismember(RelevantIndices_AfterRun,PossibleTurns_NoNaNs))
        RelevantIndices_AfterRun = [];
    end

    if ~isempty(RelevantIndices_BeforeRun)  % If a turn was found BEFORE the segment
        last_turning_frame = RelevantIndices_BeforeRun(find(ismember(RelevantIndices_BeforeRun,PossibleTurns_NoNaNs),1,'last'));
        
        RelevantFrames     = (last_turning_frame+1):(FirstFrame-1);
        if ~isempty(RelevantFrames)            
            % Is the previous segment in the opposite direction?   Use this to decide whether the SharpTurn is 'real' and not just a jitter              
            RelevantFramesForOppositeDirectionCheck = (last_turning_frame-AllowedExtensionFrameInterval):(last_turning_frame-1);
            RelevantFramesForOppositeDirectionCheck = RelevantFramesForOppositeDirectionCheck(RelevantFramesForOppositeDirectionCheck>0);                
            OppositeDirection = false;
            if ~isempty(RelevantFramesForOppositeDirectionCheck)  
                DirectionBeforeLastTurn = sign(nansum(DirectionFlag(RelevantFramesForOppositeDirectionCheck)));
                OppositeDirection       = ( CurrentDirection == -(DirectionBeforeLastTurn) );
            end                          
            if OppositeDirection
                SegmentsHighVelocity(seg_ind).DirectionWasChanged = true;
                SharpTurns(last_turning_frame) = true;                                                                                                   
            end
            % extend the straight run segment if the velocity has at least an intermediate value  
            RelevantFrames                = RelevantFrames(IntermediateVelocity(RelevantFrames));    % add only frames with at least intermediate velocity            
            % Avoid overlapping with other Straight Run segments
            MemberVec = ismember(RelevantFrames,StraightRunFrames);
            if sum(MemberVec)
                RelevantFrames = RelevantFrames(~MemberVec);
            end            
            DirectionFlag(RelevantFrames) = CurrentDirection*ones(size(RelevantFrames));
        end
        RelevantIndices_BeforeRun = RelevantFrames;
    end
    
    if ~isempty(RelevantIndices_AfterRun)  % If a turn was found AFTER the segment
        first_turning_frame = RelevantIndices_AfterRun(find(ismember(RelevantIndices_AfterRun,PossibleTurns_NoNaNs),1,'first'));
        
        RelevantFrames      = (LastFrame+1):(first_turning_frame-1);
        if ~isempty(RelevantFrames)            
            % Is the next segment in the opposite direction?   Use this to decide whether the SharpTurn is 'real' and not just a jitter              
            RelevantFramesForOppositeDirectionCheck = (first_turning_frame+1):(first_turning_frame+AllowedExtensionFrameInterval);
            RelevantFramesForOppositeDirectionCheck = RelevantFramesForOppositeDirectionCheck(RelevantFramesForOppositeDirectionCheck<CurrentTrackLength);                
            OppositeDirection = false;
            if ~isempty(RelevantFramesForOppositeDirectionCheck)  
                DirectionAfterFirstTurn = sign(nansum(DirectionFlag(RelevantFramesForOppositeDirectionCheck)));
                OppositeDirection       = ( CurrentDirection == -(DirectionAfterFirstTurn) );
            end
            if OppositeDirection
                SegmentsHighVelocity(seg_ind).DirectionWillChange = true;
                SharpTurns(first_turning_frame) = true;                                                                                                  
            end
            % extend the straight run segment if the velocity has at least an intermediate value  
            RelevantFrames                = RelevantFrames(IntermediateVelocity(RelevantFrames)); % add only frames with at least intermediate velocity
            % Avoid overlapping with other Straight Run segments
            MemberVec = ismember(RelevantFrames,StraightRunFrames);
            if sum(MemberVec)
                RelevantFrames = RelevantFrames(~MemberVec);
            end            
            DirectionFlag(RelevantFrames) = CurrentDirection*ones(size(RelevantFrames));
        end
        RelevantIndices_AfterRun = RelevantFrames;
    end
    % Assign extended values  
    FrameIndices                                          = unique([FrameIndices  RelevantIndices_BeforeRun   RelevantIndices_AfterRun]);    
    SegmentsHighVelocity(seg_ind).FirstHighVelocityFrame  = FirstFrame;
    SegmentsHighVelocity(seg_ind).LastHighVelocityFrame   = LastFrame;
    SegmentsHighVelocity(seg_ind).FramesIndices           = FrameIndices;
    SegmentsHighVelocity(seg_ind).FirstFrameIndices       = FrameIndices(1);
    SegmentsHighVelocity(seg_ind).MiddleFrameIndices      = round((FrameIndices(1)+FrameIndices(end))/2);
    SegmentsHighVelocity(seg_ind).Length                  = FrameIndices(end)-FrameIndices(1)+1;
           
end

%% Concatinate proximal straight runs if they are in the same direction and without ANY POSSIBLE sharp turns 
if AllowedPausesFrameInterval 
    SegmentsHighVelocity(1).FirstPauseFrame    =  NaN;
    SegmentsHighVelocity(1).LastPauseFrame     =  NaN;
    SegmentsHighVelocity(1).FirstJitterFrame   =  NaN;
    SegmentsHighVelocity(1).LastJitterFrame    =  NaN;
    if length(SegmentsHighVelocity)>1 
        PossibleTurns_NoNaNs_logical                       = false(size(FramesOutOfBounds));
        PossibleTurns_NoNaNs_logical(PossibleTurns_NoNaNs) = true; 

        SegmentsHighVelocity_old                   = SegmentsHighVelocity;
        Direction                = [SegmentsHighVelocity.Direction];
        UncertaintyInDirection   = [SegmentsHighVelocity.DirectionUncertaintyFlag];
        Length                   = [SegmentsHighVelocity.Length];
        FirstFrame               = [SegmentsHighVelocity.FirstFrameIndices];
        LastFrame                = FirstFrame+Length-1;
        IntervalBetweenSegments  = FirstFrame(2:end)-LastFrame(1:end-1)-1;
        CloseSegments            = IntervalBetweenSegments <= AllowedPausesFrameInterval;

        KnownDirection                      = [SegmentsHighVelocity.Direction]~=0;
        KnownDirection_In2FollowingSegments = KnownDirection(1:end-1)& KnownDirection(2:end);
        NoChangeInDirection                 = diff(Direction)==0;

        SegmentsForConcatination = CloseSegments & KnownDirection_In2FollowingSegments & NoChangeInDirection;


        clear SegmentsHighVelocity;
        SegmentsHighVelocity(1) = SegmentsHighVelocity_old(1);
        segment_index_new = 1;

        if AllowJittersWithinStraightRuns           % Don't allow concatination when the inerval consist of frames out of bound 
            Frames_Not_Allowed_For_Extension = FramesOutOfBounds;
        else                                        % Don't allow concatination when the inerval consist of frames out of bound AND when interval includes sharp turns
            Frames_Not_Allowed_For_Extension                       = FramesOutOfBounds;
            Frames_Not_Allowed_For_Extension(PossibleTurns_NoNaNs) = true;
        end        

        for segment_index_old = 2:length(SegmentsHighVelocity_old)
            Criterion_For_Concatination = SegmentsForConcatination(segment_index_old-1);  % Is this in proximity and in similar direction to the previous segment?

            if Criterion_For_Concatination
                % Check that there are no sharp turns or out of bound frames in the interval between the segments    
                FramesInInterval = (LastFrame(segment_index_old-1)+1):(FirstFrame(segment_index_old)-1);
                if sum(Frames_Not_Allowed_For_Extension(FramesInInterval))
                    Criterion_For_Concatination = false;
                end
            end

            if ~ Criterion_For_Concatination            % If it should not be concatinated with the previous segment
                segment_index_new                                          = segment_index_new + 1;
                SegmentsHighVelocity(segment_index_new)                    = SegmentsHighVelocity_old(segment_index_old);
                SegmentsHighVelocity(segment_index_new).FirstPauseFrame    =  NaN;
                SegmentsHighVelocity(segment_index_new).LastPauseFrame     =  NaN;
                SegmentsHighVelocity(segment_index_new).FirstJitterFrame   =  NaN;
                SegmentsHighVelocity(segment_index_new).LastJitterFrame    =  NaN;
            else                                                 % Concatinate with the previous segment
                FirstPauseFrame      = FramesInInterval(1);
                FirstPauseFrame      = nanmin([SegmentsHighVelocity(segment_index_new).FirstPauseFrame   FirstPauseFrame]);     % if additional concatinations accured, take the first frame in series
                LastPauseFrame       = FramesInInterval(end);            
                FirstJitterFrame     = FramesInInterval(find(PossibleTurns_NoNaNs_logical(FramesInInterval),1,'first'));          
                FirstJitterFrame     = nanmin([SegmentsHighVelocity(segment_index_new).FirstJitterFrame   FirstJitterFrame]);   % if additional concatinations accured, take the first frame in series 
                LastJitterFrame      = FramesInInterval(find(PossibleTurns_NoNaNs_logical(FramesInInterval),1,'last'));          
                CurrentFramesIndices = (SegmentsHighVelocity(segment_index_new).FirstFrameIndices):(SegmentsHighVelocity_old(segment_index_old).FramesIndices(end));
                ForwardCount         = SegmentsHighVelocity(segment_index_new).ForwardCount +  SegmentsHighVelocity_old(segment_index_old).ForwardCount;
                ReverseCount         = SegmentsHighVelocity(segment_index_new).ReverseCount +  SegmentsHighVelocity_old(segment_index_old).ReverseCount;
                if     Direction(segment_index_old)==1,  
                    CorrectVsWrongDirectionFlags    = ForwardCount/ReverseCount;
                    FractionOfDirectionFlagedFrames = ForwardCount/length(CurrentFramesIndices); 
                elseif Direction(segment_index_old)==-1, 
                    CorrectVsWrongDirectionFlags    = ReverseCount/ForwardCount;
                    FractionOfDirectionFlagedFrames = ReverseCount/length(CurrentFramesIndices); 
               end

                % re-assign to fields. The fields 'Direction', 'DirectionWasChanged' and 'FirstHighVelocityFrame' are unchanged   
                SegmentsHighVelocity(segment_index_new).FirstFrameIndices               =  CurrentFramesIndices(1);
                SegmentsHighVelocity(segment_index_new).FramesIndices                   =  CurrentFramesIndices;
                SegmentsHighVelocity(segment_index_new).MiddleFrameIndices              =  round((CurrentFramesIndices(1)+CurrentFramesIndices(end))/2);
                SegmentsHighVelocity(segment_index_new).Length                          =  length(CurrentFramesIndices);
                SegmentsHighVelocity(segment_index_new).ForwardCount                    =  ForwardCount;
                SegmentsHighVelocity(segment_index_new).ReverseCount                    =  ReverseCount;
                SegmentsHighVelocity(segment_index_new).CorrectVsWrongDirectionFlags    =  CorrectVsWrongDirectionFlags;
                SegmentsHighVelocity(segment_index_new).FractionOfDirectionFlagedFrames =  FractionOfDirectionFlagedFrames;
                SegmentsHighVelocity(segment_index_new).DirectionWillChange             =  SegmentsHighVelocity_old(segment_index_old).DirectionWillChange;            
                SegmentsHighVelocity(segment_index_new).LastHighVelocityFrame           =  SegmentsHighVelocity_old(segment_index_old).LastHighVelocityFrame;
                DirectionFlag(CurrentFramesIndices)                                     =  Direction(segment_index_old);
                SegmentsHighVelocity(segment_index_new).DirectionUncertaintyFlag        =  ...
                       (SegmentsHighVelocity(segment_index_new).DirectionUncertaintyFlag) || (UncertaintyInDirection(segment_index_old));
                SegmentsHighVelocity(segment_index_new).FirstPauseFrame                 =  FirstPauseFrame;
                SegmentsHighVelocity(segment_index_new).LastPauseFrame                  =  LastPauseFrame;
                SegmentsHighVelocity(segment_index_new).FirstJitterFrame                =  FirstJitterFrame;
                SegmentsHighVelocity(segment_index_new).LastJitterFrame                 =  LastJitterFrame;
            end    
        end    
    end
end

return

function SegmentsOutOfBound  = ExtendOutOfBound_OnlyBasedOnFrameInterval (SegmentsOutOfBound, NumOfFrames, MinFramesForShortRun)
%% Concatinate proximal Out-of-bounds segments if there is less than 'MinFramesForShortRun' of HighVelocity frames within the interval between segments 
if isempty(SegmentsOutOfBound)
    return    
end
    
AllOutOfBoundFrames_logical                                     = false(1,NumOfFrames);
AllOutOfBoundFrames_logical([SegmentsOutOfBound.FramesIndices]) = true;
Length                   = [SegmentsOutOfBound.Length];
FirstFrame               = [SegmentsOutOfBound.FirstFrameIndices];
LastFrame                = FirstFrame+Length-1;

for segment_index = 2:length(SegmentsOutOfBound)
    FramesInInterval = (LastFrame(segment_index-1)+1):(FirstFrame(segment_index)-1);
    if ~isempty(FramesInInterval)
        Criterion_For_Concatination = length(FramesInInterval) < MinFramesForShortRun ;          
        if Criterion_For_Concatination
            AllOutOfBoundFrames_logical(FramesInInterval) = true;
        end
    end
end 

SegmentsOutOfBound  = FindSegments  ( find(AllOutOfBoundFrames_logical),  true);

return

function [Segments, MatchingMatrix] = CombineSegments (SegmentsIn , Field1, Field2, Operator, NewFieldName )
%% Find INTERSECTION between two fields and for each intersection combine the fields by either an OR or and AND operator (on the frames).    
%% This function will generate new fields: 'Field1_OR_Field2' and 'Field1_AND_Field2', respectively;    
% Inputs:
%   SegmentsIn        =  Segment structure
%   Field1 and Field2 = strings with the name of fields
%   Operator          = cell array which contains the operation of interest. 3 Possible values: {'OR'} / {'AND'} / {'OR','AND'}  
%   NewFieldName      = optional. cell array with strings for the output fields names. length(NewFieldName) == length(Operator)  is a must!   
%
% Outputs:
%   Segment           = Segment structure with the additional fields of interest  
%   MatchingMatrix    = a logical matrix of size (#SegmentsField1, #SegmentsField2, #Operators). Shows true when combined frames were detected.
%
%%
if isempty(SegmentsIn) 
    Segments       = [];
    MatchingMatrix = [];
    return;
end
 
Segments = SegmentsIn;

Output_Fields         = cell(1,length(Operator));
Output_SegmentIndices = zeros(1,length(Operator));
for f_ind = 1:length(Operator)
    Operation_name = Operator{f_ind};
    Output_Fields{f_ind} = [Field1,'_AND_',Field2,'_FramesOperation_',Operation_name];
end

if exist('NewFieldName','var')
    if length(NewFieldName)== length(Operator)
        Output_Fields = NewFieldName;
    else
        disp('ignoring input field names');
    end
end

SegmentsField1 = SegmentsIn.(Field1);
SegmentsField2 = SegmentsIn.(Field2);

%% In the case of at least one empty segment: there are no mutual states to analyze --> return empty segments 
if isempty(SegmentsField1) || isempty(SegmentsField2)
    for f_ind = 1:length(Output_Fields)           
        Segments.(Output_Fields{f_ind}) = [];                     
    end
    MatchingMatrix=[];
    return
end

%% The two segments are not empty
MatchingMatrix = false(length(SegmentsField1),length(SegmentsField2), length(Output_Fields));

for seg_index1 = 1:length(SegmentsField1)
    Frames1 = SegmentsField1(seg_index1).FramesIndices; 
    for seg_index2 = 1:length(SegmentsField2)
        Frames2            = SegmentsField2(seg_index2).FramesIndices; 
        FramesInBothFields = intersect(Frames1,Frames2);   
        
        if ~isempty(FramesInBothFields)   % Update a new field in Segments if mutual frames were found
        
            for f_ind = 1:length(Output_Fields)            

                    if strcmpi(Operator{f_ind},'AND')       % Looking for mutual frames                   --> intersect
                        Frames_NewField = FramesInBothFields;                
                    elseif strcmpi(Operator{f_ind},'OR')    % Looking for either one of the fields frames --> Union
                        Frames_NewField = union(Frames1,Frames2);     
                    end

                    Output_SegmentIndices(f_ind) = Output_SegmentIndices(f_ind) + 1;
                    new_segment_ind              = Output_SegmentIndices(f_ind);
                    
                    Segments.(Output_Fields{f_ind})(new_segment_ind).FirstFrameIndices  =  Frames_NewField(1);
                    Segments.(Output_Fields{f_ind})(new_segment_ind).FramesIndices      =  Frames_NewField;
                    Segments.(Output_Fields{f_ind})(new_segment_ind).MiddleFrameIndices =  floor(Frames_NewField(1)+(length(Frames_NewField)-1)/2);
                    Segments.(Output_Fields{f_ind})(new_segment_ind).Length             =  length(Frames_NewField);
                    
                    MatchingMatrix(seg_index1, seg_index2, f_ind) = true;
            end 
        end
    end
end

% Remove possible repeated segments
if isfield(Segments,Output_Fields{1})
    FirstFrameIndices              = [Segments.(Output_Fields{1}).FirstFrameIndices];
    RepeatedSegmentsIndices        = [0 diff(FirstFrameIndices)==0];
    
    if ~isempty(find(RepeatedSegmentsIndices,1))
        NoRepeatsSegments = SegmentsIn;
        for f_ind = 1:length(Output_Fields)   
            NoRepeatsSegments.(Output_Fields{f_ind}) = Segments.(Output_Fields{f_ind})(~RepeatedSegmentsIndices);
        end
        Segments = NoRepeatsSegments;
    end
    
else   % If no combination of fields was found then generate fields that are empty   
    for f_ind = 1:length(Output_Fields)   
        Segments.(Output_Fields{f_ind}) = [];        
    end
end

return

function Segments = CombineSegments_OR_if_minimal_overlap_AND_otherwise (SegmentsIn , Field1, Field2, MinOverlap_ForOROperator)

%% Find INTERSECTION of Field1 relative to Field2: combine the fields by either an OR or and AND operator (on the frames) DEPENDING ON THEIR OVERLAP.
% For example:   
%   Segments = CombineSegments_OR_if_minimal_overlap_AND_otherwise     (Segments , 'OmegaShape', 'LowVelocity_LongSegments', Thresholds.Omega.MinOverlapWithPauseForOmegaExtension);   
%   This function will make a new field called:    'OmegaShape_ExtendedBy_LowVelocity_LongSegments'  
%   The new field will extend the OmegaShape segments if there is a significant overlap with a 'LowVelocity_LongSegments' state.
%   Significance is set by: 'MinOverlapWithPauseForOmegaExtension'
%     for MinOverlapWithPauseForOmegaExtension == 50%,  an OmegaShape segment with 10 frames will be extended if it overlaps at least with 50% of the 'LowVelocity_LongSegments' frames. 
%     in this example, the criterion cannot be reached if the length of the 'LowVelocity_LongSegments' is more than twice of the 'OmegaShape' segment.        


%% This function will generate new fields: 'Field1_ExtendedBy_Field2'
% Inputs:
%   SegmentsIn        =  Segment structure
%   Field1 and Field2 = strings with the name of fields
%   MinOverlap_ForOROperator =  if the (overlap>MinOverlap_ForOROperator)  --> use OR operation on intersected states frames     
%                               if the (overlap<MinOverlap_ForOROperator)  --> do not extend the state           
% Outputs:
%   Segment           = Segment structure with the additional field of interest  
%
%%
if isempty(SegmentsIn) 
    Segments       = [];
    return;
end
 
Segments       = SegmentsIn;
SegmentsField1 = SegmentsIn.(Field1);
SegmentsField2 = SegmentsIn.(Field2);
Output_Field   = [Field1,'_ExtendedBy_',Field2];

%% In the case of at least one empty segment: there are no mutual states to analyze 
if isempty(SegmentsField1)|| isempty(SegmentsField2) 
    Segments.(Output_Field) = [];      
    return
end

%% The two segments are not empty
new_segment_ind = 0;
for seg_index1 = 1:length(SegmentsField1)
    Frames1 = SegmentsField1(seg_index1).FramesIndices; 
    for seg_index2 = 1:length(SegmentsField2)
        Frames2            = SegmentsField2(seg_index2).FramesIndices; 
        FramesInBothFields = intersect(Frames1,Frames2);         
        
        if ~isempty(FramesInBothFields)   % Update a new field in Segments if mutual frames were found            
            new_segment_ind = new_segment_ind+1; 
            
            FramesInEitherOneOfTheFields = union(Frames1,Frames2);
            CurrentOverlap  = length(FramesInBothFields)/length(FramesInEitherOneOfTheFields) * 100;   % Percentage of intersected frames versus union frames 

            if CurrentOverlap <  MinOverlap_ForOROperator       % No good overlap: use an AND operator (intersect)
                Frames_NewField = FramesInBothFields;                
            else                                                % Good overlap: use an OR operator (Union)
                Frames_NewField = FramesInEitherOneOfTheFields;     
            end

            Segments.(Output_Field)(new_segment_ind).FirstFrameIndices  =  Frames_NewField(1);
            Segments.(Output_Field)(new_segment_ind).FramesIndices      =  Frames_NewField;
            Segments.(Output_Field)(new_segment_ind).MiddleFrameIndices =  floor(Frames_NewField(1)+(length(Frames_NewField)-1)/2);
            Segments.(Output_Field)(new_segment_ind).Length             =  length(Frames_NewField);                             
        end
    end
end

% Remove possible repeated segments
if isfield(Segments,Output_Field)
    FirstFrameIndices              = [Segments.(Output_Field).FirstFrameIndices];
    RepeatedSegmentsIndices        = [0 diff(FirstFrameIndices)==0];
    
    if ~isempty(find(RepeatedSegmentsIndices,1))
        NoRepeatsSegments = SegmentsIn;
        NoRepeatsSegments.(Output_Field) = Segments.(Output_Field)(~RepeatedSegmentsIndices); 
        Segments         = NoRepeatsSegments;
    end
    
else   % If no combination of fields was found then generate fields that are empty   
    Segments.(Output_Field) = []; 
end

return

function [Segments, IndicesOfGoodSegments] = FilterSegment (SegmentsIn , FramesNotAllowed)
%%%% Return a Segment structure that is filtered from segments that INCLUDES Frames that are not allowed. Keep all fields within that structure.     
% Inputs:
%   SegmentsIn        = Segment structure
%   FramesNotAllowed  = logical vector. true for frames that are not allowed. 
%
% Outputs:
%   Segment               = Segment structure with the additional fields of interest  
%   IndicesOfGoodSegments = a logical vector of length (SegmentsIn). true for indices that pass the filter. 

if isempty(SegmentsIn)
    Segments       = [];
    IndicesOfGoodSegments = [];
    return;
end

IndicesOfGoodSegments = false(1,length(SegmentsIn));
for seg_index = 1:length(SegmentsIn)
    IndicesOfGoodSegments(seg_index) = ~ logical(sum(FramesNotAllowed([SegmentsIn(seg_index).FramesIndices])));  % Segment is good only if there is no 'FramesNotAllowed' within its range
end

Segments = SegmentsIn(IndicesOfGoodSegments);

return

function [Segments_HighLevelBehavior, Logical_Vectors_HighLevel] = CalculateHighLevelBehavior(Segments_LowLevelBehavior, Logical_Vectors, Segments_LongReverse, Thresholds, Logical_LongPauses, CurrentTrack) 
                                                        
MaxFramesBetweenReverseAndOmega             = Thresholds.Pirouette.MaxFramesBetweenReverseAndOmega;
MaxFramesBetweenLastReverse_And_Forward     = Thresholds.Pirouette.MaxFramesBetweenLastReverse_And_Forward; 
MaxPauseFramesWithinPirouetteReversal       = Thresholds.Pirouette.MaxPauseFramesWithinPirouetteReversal;
                                                        
Omega          = Logical_Vectors.Omega;
Forward        = Logical_Vectors.Forward;
Reverse        = Logical_Vectors.Reverse;
Pause          = Logical_Vectors.Pause;
OutOfBounds    = Logical_Vectors.OutOfBounds;

% Forward_RunDirection = vector. NaN if not a forward run. values equal the previously calculated forward segment RunDirection. 
%  THIS IS NOT THE INSTANTANOUS ANGLE VECTOR !!!
Forward_RunDirection             = single(zeros(size(Logical_Vectors.Forward))*NaN);
Forward_FractionOfCurvingFrames  = single(zeros(size(Logical_Vectors.Forward))*NaN);
Forward_SegmentStartsWithCurving = false(size(Logical_Vectors.Forward));

for FWD_ind = 1:length(Segments_LowLevelBehavior.Forward)
    Forward_RunDirection(Segments_LowLevelBehavior.Forward(FWD_ind).FramesIndices)             =  Segments_LowLevelBehavior.Forward(FWD_ind).Angle_RunDirection;
    Forward_FractionOfCurvingFrames(Segments_LowLevelBehavior.Forward(FWD_ind).FramesIndices)  =  Segments_LowLevelBehavior.Forward(FWD_ind).FractionOfCurvingFrames;
    Forward_SegmentStartsWithCurving(Segments_LowLevelBehavior.Forward(FWD_ind).FramesIndices) =  Segments_LowLevelBehavior.Forward(FWD_ind).SegmentStartsWithCurving;
end

UncertaintyFrames = OutOfBounds + Logical_LongPauses + Omega;  % If the worm has one of these states before Pirouette initiation --> starting time is not well defined.
ReverseOrOmega    = Reverse + Omega;
NumberOfFrames    = length(Omega);

% ReverseSegments                 = Segments_LowLevelBehavior.Reverse;
if ~isempty(Segments_LongReverse)
    ReverseSegments                 = Segments_LongReverse;

    NumberOfReversalSegments        = length(ReverseSegments);
    ReversalsStartTime              = [ReverseSegments.FirstFrameIndices];
    ReversalsLength                 = [ReverseSegments.Length];
    ReversalsEndTime                = ReversalsStartTime+ReversalsLength-1;
    ReversalsFirstPauseFrame        = [ReverseSegments.FirstPauseFrame];
    ReversalsFirstJitterFrame       = [ReverseSegments.FirstJitterFrame];
    ReversalsFirstHighVelocityFrame = [ReverseSegments.FirstHighVelocityFrame];
    ReversalsStartTime              = [ReverseSegments.FirstFrameIndices];
    ReversalsAngle_InitialFrame     = [ReverseSegments.Angle_InitialFrame];
    AngleBeforeReversal             = mod(ReversalsAngle_InitialFrame+360,360)-180  ;

    ReversalsEndTime_Until_Pause_Or_ReversalEnd                                    = ReversalsEndTime;
    ReversalsEndTime_Until_Pause_Or_ReversalEnd(~isnan(ReversalsFirstPauseFrame))  = ReversalsFirstPauseFrame(~isnan(ReversalsFirstPauseFrame));
    ReversalsEndTime_Until_Jitter_Or_ReversalEnd                                   = ReversalsEndTime;
    ReversalsEndTime_Until_Jitter_Or_ReversalEnd(~isnan(ReversalsFirstJitterFrame))= ReversalsFirstJitterFrame(~isnan(ReversalsFirstJitterFrame));

    PirouetteSegments               = false(1,NumberOfReversalSegments);
    Pirouette_seg_ind               = 0;
    FramesAlreadyUsedAsPirouettes   = false(1,NumberOfFrames);

    for seg_ind = 1:NumberOfReversalSegments
        FrameIntervalAfterReversal = (ReversalsEndTime(seg_ind)+1):(ReversalsEndTime(seg_ind)+MaxFramesBetweenReverseAndOmega);
        FrameIntervalAfterReversal = FrameIntervalAfterReversal(FrameIntervalAfterReversal<NumberOfFrames);   % to ovoid movie edge effect
        FirstOmegaFrameInRange     = FrameIntervalAfterReversal(find(Omega(FrameIntervalAfterReversal),1,'first'));

        if FirstOmegaFrameInRange  % Omega Exists
            % Check if there is a Forward motion RIGHT AFTER the omega segment    
            OmegaLength     = find(Omega(FirstOmegaFrameInRange:end)==0,1,'first')-1;
            LastOmegaFrame  = FirstOmegaFrameInRange+ OmegaLength -1; 
            FrameAfterOmega = LastOmegaFrame+1; 
            if FrameAfterOmega < NumberOfFrames
                ReverseAfterOmega                = Reverse(FrameAfterOmega); 
                NoForwardBeforeOmega             = ~logical(sum(Forward(FrameIntervalAfterReversal)));

                if ReverseAfterOmega && NoForwardBeforeOmega   % Concatination of additional reversals is required
                    % If there is a REVERSE immediately after the omega then look for omegas of the type: (Reverse)-->(omega)-->(reverse with short pauses and omegas)--> forward              
                    % i.e. allow the Forward run to take place a while after the LAST run ending.
                    
                    SkipForwardAfterOmegaConditionCheck = false;
                    ContinueReversalConcatination       = true;
                    CurrentFrame                        = FrameAfterOmega;
                    while ContinueReversalConcatination  
                        ContinueReversalConcatination    = false;
                        ConcatinatedReverseSegmentLength = find(ReverseOrOmega(CurrentFrame:end)==0,1,'first')-1;
                        LastFrameOfConcatinatedSegment   = CurrentFrame + ConcatinatedReverseSegmentLength -1; 
                        LastPirouetteFrame               = LastFrameOfConcatinatedSegment;                % Correct the Pirouette definition to include 
                        if Pause(LastPirouetteFrame+1) 
                            % if the next frame after this concatinated reversal is a Pause, check whether there is an additional reversal within the range. Is so continue concatination. 
                            FirstPauseFrame    = LastPirouetteFrame+1;
                            PauseSegmentLength = find(Pause((LastPirouetteFrame+1):end)==0,1,'first')-1;
                            FrameAfterPause    = FirstPauseFrame + PauseSegmentLength; 
                            if ~isempty(FrameAfterPause)
                                if (PauseSegmentLength <= MaxPauseFramesWithinPirouetteReversal) && (ReverseOrOmega(FrameAfterPause))
                                    % Add this pause segment to the Pirouette reversal and continue concatination of the next reversal segment   
                                    ContinueReversalConcatination  = true;
                                    CurrentFrame                   = FrameAfterPause;
                                else
                                    LastPirouetteFrame = FrameAfterPause - 1; % Include the last pause segment in the pirouette
                                end
                            else
                                ForwardAfterOmega                   = false;      % The pause doesn't end until the end of the track, i.e. there is no Forward motion after the omega
                                SkipForwardAfterOmegaConditionCheck = true;
                            end
                        end
                    end
                    % Find whether the forward run is within range (yes=true) 
                    if ~ SkipForwardAfterOmegaConditionCheck
                        LastFrameToCheck   = min([(LastPirouetteFrame+MaxFramesBetweenLastReverse_And_Forward) , length(Forward)]);
                        ForwardAfterOmega  = logical(sum(Forward((LastPirouetteFrame+1):LastFrameToCheck)));  
                    end

                else                     % If there is no additional reversals, the next frame must be a forward run, i.e. expect a classical piruoette sequence: (Reverse --> omega --> Forward)    
                    ForwardAfterOmega        = Forward(FrameAfterOmega);    
                    LastPirouetteFrame       = LastOmegaFrame;
                end                
                FramesIndices                    = (ReversalsStartTime(seg_ind)):LastPirouetteFrame;
                OutOfBound_WithinPirouette       = logical(sum(Logical_Vectors.OutOfBounds(FramesIndices)));   % make sure that there are no out of bound frames during the pirouette sequence
                FramesArePartOfPreviousPirouette = logical(sum(FramesAlreadyUsedAsPirouettes(FramesIndices))); % make sure that this is not part of a previously defined long pirouette  

                % Check the pirouette condition. The existance of a Reverse and an omega is already known here.    
                PirouetteSegments(seg_ind)       = ForwardAfterOmega && (~OutOfBound_WithinPirouette)&&(~FramesArePartOfPreviousPirouette);             
            end        
            if PirouetteSegments(seg_ind)  % If it is a piruoette sequence: (Reverse --> omega --> Forward) 
                Pirouette_seg_ind                            = Pirouette_seg_ind+1;   
                FramesAlreadyUsedAsPirouettes(FramesIndices) = true;
                PirouetteWithForwardBeforeOmega              = false;  % initialization. Criterion is checked below         

                % Properties of Forward motion AFTER the pirouette ends
                RelevantFrames                                   = Logical_Vectors.Forward;
                RelevantFrames(1:LastPirouetteFrame)             = false;                     % Only forward frames after this pirouette
                FirstForwardIndex                                = find(RelevantFrames,1,'first');

                AngleAfterPirouette                               = Forward_RunDirection(FirstForwardIndex);
                AfterPirouette_SegmentStartsWithCurving           = Forward_SegmentStartsWithCurving(FirstForwardIndex);
                AfterPirouette_FractionOfCurvingFrames            = Forward_FractionOfCurvingFrames(FirstForwardIndex);           

                TotalChangeInAngle = mod((AngleAfterPirouette - AngleBeforeReversal(seg_ind) + 360), 360) - 180;

                % Is it a Forward type piruoette? i.e. there is a short forward segment between the reverse and the omega  
                ForwardBeforeOmega_StartTime = FrameIntervalAfterReversal(find(Forward(FrameIntervalAfterReversal),1,'first'));           
                if ~isempty(ForwardBeforeOmega_StartTime) % a Forward type piruoette                
                    PirouetteWithForwardBeforeOmega = (ForwardBeforeOmega_StartTime < FirstOmegaFrameInRange);  % Check that this forard is indeed before the omega. If yes, criterion is met
                end                      
                if PirouetteWithForwardBeforeOmega % a Forward type piruoette                
                    ForwardBeforeOmega_Length  = find(Omega(ForwardBeforeOmega_StartTime:end)==0,1,'first')-1;
                    ForwardBeforeOmega_EndTime = ForwardBeforeOmega_StartTime+ ForwardBeforeOmega_Length -1; 
                else
                    ForwardBeforeOmega_Length  = [];
                    ForwardBeforeOmega_EndTime = [];
                end

                % evaluate how good is the estimated pirouette staring time
                FirstFrameIndex                    = ReversalsStartTime(seg_ind);
                InitiationTimeIsWellDefined        = false;
                InitiationTimeIsNotAfterOutOfBound = false;
                if FirstFrameIndex~=1                          % Inititation time is well defined only if the pirouettes doesn't start at the FIRST FRAME of the track
                    if ~(UncertaintyFrames(FirstFrameIndex-1)) % Inititation time is well defined only if the pirouettes doesn't start at an UNCERTAINTY FRAME like out of bounds. 
                        InitiationTimeIsWellDefined = true;
                    end
                    if ~OutOfBounds(FirstFrameIndex-1)         % If the worm was NOT (out of bounds) before the reversal initiated
                        InitiationTimeIsNotAfterOutOfBound = true;
                    end
                end

                % get here all relevant features and assign them into the  Piruoette Segment    
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).FirstFrameIndices               = FirstFrameIndex;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).LastFrameIndices                = LastPirouetteFrame;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).FramesIndices                   = FramesIndices;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).Length                          = (LastPirouetteFrame-ReversalsStartTime(seg_ind))+1;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).MiddleFrameIndices              = round((LastPirouetteFrame+ReversalsStartTime(seg_ind))/2);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).FirstReversalHighVelocityFrame  = ReversalsFirstHighVelocityFrame(seg_ind);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).ReversalsEndTime                = ReversalsEndTime(seg_ind);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).FirstReversalHighVelocityFrame  = ReversalsFirstHighVelocityFrame(seg_ind);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).ReversalsFirstPauseFrame        = ReversalsFirstPauseFrame(seg_ind);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).ReversalsFirstJitterFrame       = ReversalsFirstJitterFrame(seg_ind);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).ReversalsLength                 = ReversalsLength(seg_ind);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).ReversalsEndTime_Until_Pause_Or_ReversalEnd  = ReversalsEndTime_Until_Pause_Or_ReversalEnd(seg_ind);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).ReversalsEndTime_Until_Jitter_Or_ReversalEnd = ReversalsEndTime_Until_Jitter_Or_ReversalEnd(seg_ind);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).FramesUntilPause                = FirstFrameIndex:ReversalsEndTime_Until_Pause_Or_ReversalEnd(seg_ind);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).InitialReversalFrames           = FirstFrameIndex:ReversalsEndTime(seg_ind);

                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).PirouetteWithForwardBeforeOmega = PirouetteWithForwardBeforeOmega;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).ForwardBeforeOmega_Length       = ForwardBeforeOmega_Length;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).ForwardBeforeOmega_EndTime      = ForwardBeforeOmega_EndTime;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).InitiationTimeIsWellDefined     = InitiationTimeIsWellDefined;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).InitiationTimeIsNotAfterOutOfBound = InitiationTimeIsNotAfterOutOfBound;

                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).AngleBeforeReversal             = AngleBeforeReversal(seg_ind);
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).AngleAfterPirouette             = AngleAfterPirouette;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).TotalChangeInAngle              = TotalChangeInAngle;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).AfterPirouette_SegmentStartsWithCurving = AfterPirouette_SegmentStartsWithCurving;
                Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ind).AfterPirouette_FractionOfCurvingFrames  = AfterPirouette_FractionOfCurvingFrames;                      

            end
        end          
    end        
end

%% Filter Pirouettes if they were counted that were counted twice due to two multiple reversals before the omega 
Logical_Vectors_HighLevel.Pirouette                          = false(1,NumberOfFrames);
Logical_Vectors_HighLevel.Pirouette_ExactInitiationTime      = false(1,NumberOfFrames);
Logical_Vectors_HighLevel.Pirouette_NotAfterOutOfBound       = false(1,NumberOfFrames);
Logical_Vectors_HighLevel.Pirouette_UntilPause               = false(1,NumberOfFrames);
Logical_Vectors_HighLevel.Pirouette_AfterPause               = false(1,NumberOfFrames);
Logical_Vectors_HighLevel.Pirouette_InitialReversal          = false(1,NumberOfFrames);
Logical_Vectors_HighLevel.Pirouette_AfterReversal            = false(1,NumberOfFrames);

if exist('Segments_HighLevelBehavior','var')
    Logical_Vectors_HighLevel.Pirouette([Segments_HighLevelBehavior.Pirouette.FramesIndices]) = true;
    PirouetteFrames                                                                           = [Segments_HighLevelBehavior.Pirouette.FramesIndices];
%     Pirouettes_Are_Well_Defined_If_Value_Equals_Zero                                          = length(PirouetteFrames)-length(unique(PirouetteFrames))
    if length(PirouetteFrames) ~= length(unique(PirouetteFrames))
        disp('Pirouettes are not well defined')
    end
    
    Pirouette_seg_ExactTiming                                                                 = [Segments_HighLevelBehavior.Pirouette.InitiationTimeIsWellDefined];
    Logical_Vectors_HighLevel.Pirouette_ExactInitiationTime([Segments_HighLevelBehavior.Pirouette(Pirouette_seg_ExactTiming).FramesIndices])  = true;
    Pirouette_seg_NotAfterOutOfBound                                                         = [Segments_HighLevelBehavior.Pirouette.InitiationTimeIsNotAfterOutOfBound];
    Logical_Vectors_HighLevel.Pirouette_NotAfterOutOfBound([Segments_HighLevelBehavior.Pirouette(Pirouette_seg_NotAfterOutOfBound).FramesIndices])  = true;    
    Logical_Vectors_HighLevel.Pirouette_UntilPause([Segments_HighLevelBehavior.Pirouette.FramesUntilPause])       = true;
    Logical_Vectors_HighLevel.Pirouette_AfterPause    = (Logical_Vectors_HighLevel.Pirouette)&(~Logical_Vectors_HighLevel.Pirouette_UntilPause); 
    Logical_Vectors_HighLevel.Pirouette_InitialReversal([Segments_HighLevelBehavior.Pirouette.InitialReversalFrames])  = true;
    Logical_Vectors_HighLevel.Pirouette_AfterReversal    = (Logical_Vectors_HighLevel.Pirouette)&(~Logical_Vectors_HighLevel.Pirouette_InitialReversal); 
else
    Segments_HighLevelBehavior.Pirouette = [];
end

%% Assign Values to high level behavior structure
% Out of bounds and forward remain the same
Logical_Vectors_HighLevel.OutOfBounds                                                     = Logical_Vectors.OutOfBounds;
Segments_HighLevelBehavior.OutOfBounds                                                    = Segments_LowLevelBehavior.OutOfBounds;

% Forward. Filter forward segments that are within a pirouette sequence (Forward-type pirouette = (reversal)>(short FWD)>(omega)>(FWD) )   
Logical_Vectors_HighLevel.Forward                                                         = false(1,NumberOfFrames);
if ~isempty(Segments_LowLevelBehavior.Forward)
    [Segments_HighLevelBehavior.Forward , IndicesOfGoodSegments]                              = FilterSegment (Segments_LowLevelBehavior.Forward  , Logical_Vectors_HighLevel.Pirouette); 
    Logical_Vectors_HighLevel.Forward([Segments_HighLevelBehavior.Forward.FramesIndices])     = true;
else
    Segments_HighLevelBehavior.Forward = []; 
end

% Reversals. Filter Reversals segments that are within a pirouette sequence    
Logical_Vectors_HighLevel.Reversal                                                        = false(1,NumberOfFrames);
if ~isempty(Segments_LowLevelBehavior.Reverse)
    [Segments_HighLevelBehavior.Reversal , IndicesOfGoodSegments]                             = FilterSegment (Segments_LowLevelBehavior.Reverse  , Logical_Vectors_HighLevel.Pirouette); 
    Logical_Vectors_HighLevel.Reversal([Segments_HighLevelBehavior.Reversal.FramesIndices])   = true;
else
    Segments_HighLevelBehavior.Reverse = []; 
end

% Omegas without Pirouettes
Logical_Vectors_HighLevel.OmegaWithoutPirouette                                           = false(1,NumberOfFrames);
if ~isempty(Segments_LowLevelBehavior.Omega)
    [Segments_HighLevelBehavior.OmegaWithoutPirouette , IndicesOfGoodSegments]                = FilterSegment (Segments_LowLevelBehavior.Omega  ,    Logical_Vectors_HighLevel.Pirouette); 
    Logical_Vectors_HighLevel.OmegaWithoutPirouette([Segments_HighLevelBehavior.OmegaWithoutPirouette.FramesIndices]) = true;
else
    Segments_HighLevelBehavior.OmegaWithoutPirouette = [];
end

% Sharp Turns
Logical_Vectors_HighLevel.SharpTurns                                                      = false(1,NumberOfFrames);
if ~isempty(Segments_LowLevelBehavior.SharpTurns)
    [Segments_HighLevelBehavior.SharpTurns , IndicesOfGoodSegments]                           = FilterSegment (Segments_LowLevelBehavior.SharpTurns , Logical_Vectors_HighLevel.Pirouette); 
    Logical_Vectors_HighLevel.SharpTurns([Segments_HighLevelBehavior.SharpTurns.FramesIndices])= true;
else
    Segments_HighLevelBehavior.SharpTurns = [];
end

% Pauses --> All the rest. Some of the pauses are in the low level will be considered part of a Pirouette 
Logical_Vectors_HighLevel.Pause                                                          = false(1,NumberOfFrames);
if ~isempty(Segments_LowLevelBehavior.Pause)
    [Segments_HighLevelBehavior.Pause , IndicesOfGoodSegments]                               = FilterSegment (Segments_LowLevelBehavior.Pause , Logical_Vectors_HighLevel.Pirouette); 
    Logical_Vectors_HighLevel.Pause([Segments_HighLevelBehavior.Pause.FramesIndices])        = true;
else
    Segments_HighLevelBehavior.Pause = [];
end

return

function [locs_all, pks_all, locs, pks] = FindPeaks_Allow3FramesIntegration(vec,ThresholdValue, plotme)

% locs, pks          = The locations and peaks values without allowing integration with neighbouring frames
% locs_all, pks_all  = The locations and peaks allowing integration with neighbouring frames.     
%                      NOTE !!! the values in pks_all are the integrated values only for peaks that were found by integration!!  

% Example: [locs_all, pks_all] = FindPeaks_Allow3FramesIntegration(CurrentTrack.Locomotion_ChangeInVelocity_Angle, 60, true)
% find sudden angle changes

% peaks above threshold:
[pks1,locs1] = findpeaks(double(vec),'MINPEAKHEIGHT',ThresholdValue);
[pks2,locs2] = findpeaks(double(-vec),'MINPEAKHEIGHT',ThresholdValue);
locs         = [locs1  locs2];
pks          = [pks1  -pks2];
[locs,IX]    = sort(locs);
pks          = pks(IX);

% true if it's a peak location or its neighbor
locs_logical               = false(size(vec));
locs_logical(locs)         = true;
locs_logical_WithNeighbors = locs_logical;
locs_logical_WithNeighbors(1:end-1) = locs_logical_WithNeighbors(1:end-1)+ locs_logical(2:end);             
locs_logical_WithNeighbors(2:end)   = locs_logical_WithNeighbors(2:end)  + locs_logical(1:end-1);    

% find also peaks that are more shallow, but represent an overall large angle change (integral...)  
%  Algorithm: sum 2 angle changes: For each frame add the value of its two neighbours. Then find peaks...    
vec_sum          = vec;             
vec_sum(1:end-1) = vec_sum(1:end-1)+ vec(2:end);             
vec_sum(2:end)   = vec_sum(2:end)  + vec(1:end-1);     

[pks1,locs1] = findpeaks(double(vec_sum),'MINPEAKHEIGHT',ThresholdValue);
[pks2,locs2] = findpeaks(double(-vec_sum),'MINPEAKHEIGHT',ThresholdValue);
locs_sum     = [locs1  locs2];
pks_sum      = [pks1  -pks2];
[locs_sum,IX]= sort(locs_sum);
pks_sum      = pks_sum(IX);
         
locs_sum_logical            = false(size(vec));
locs_sum_logical(locs_sum)  = true;
pks_sum_values              = zeros(size(vec));
pks_sum_values(locs_sum)    = pks_sum;

% find the peaks that were found by integration and add them to the pks pool  
locs_sum_NotPreviouslyDiscovered = find(locs_sum_logical .* ~ locs_logical_WithNeighbors);
pks_sum_NotPreviouslyDiscovered  = pks_sum_values(locs_sum_NotPreviouslyDiscovered);

locs_all     = [locs  locs_sum_NotPreviouslyDiscovered];
pks_all      = [pks   pks_sum_NotPreviouslyDiscovered];
[locs_all,IX]= sort(locs_all);
pks_all      = pks_all(IX);

if plotme
    figure('name','Peaks of changes in velocity- with vs. without integration'); 
    plot(vec,'k.-'); hold on; plot(locs,pks,'r*'); 
    plot(vec_sum,'b.:'); hold on; plot(locs_sum,pks_sum,'go'); 

    figure('name','Peaks of changes in velocity- All'); 
    plot(vec,'k.-'); hold on; plot(locs_all,pks_all,'r*'); 
end

return

function [locs_all, pks_all, locs, pks] = FindPeaks_Allow9FramesIntegration(vec,ThresholdValue, plotme)
% HOW MANY FRAMES FOR INTEGRATION????  Allow 300 micro-seconds
% locs, pks          = The locations and peaks values without allowing integration with neighbouring frames
% locs_all, pks_all  = The locations and peaks allowing integration with neighbouring frames.     
%                      NOTE !!! the values in pks_all are the integrated values only for peaks that were found by integration!!  

% Example: [locs_all, pks_all] = FindPeaks_Allow3FramesIntegration(CurrentTrack.Locomotion_ChangeInVelocity_Angle, 60, true)
% find sudden angle changes

% peaks above threshold:
[pks1,locs1] = findpeaks(double(vec),'MINPEAKHEIGHT',ThresholdValue);
[pks2,locs2] = findpeaks(double(-vec),'MINPEAKHEIGHT',ThresholdValue);
locs         = [locs1  locs2];
pks          = [pks1  -pks2];
[locs,IX]    = sort(locs);
pks          = pks(IX);

% neighbors filter proeperties
NumberOfNeighborFrames = 9;
b                      = ones(1,NumberOfNeighborFrames); % this is summing up, not taking the mean
a                      = 1;

% true if it's a peak location or if by adding neighbor frames it becomes a peak 
locs_logical               = false(size(vec));
locs_logical(locs)         = true;
locs_logical_WithNeighbors = filter(b,a,locs_logical);
delay                      = floor(NumberOfNeighborFrames/2);
locs_logical_WithNeighbors = [locs_logical_WithNeighbors((delay+1):end), locs_logical((end-delay+1):end)];

% find also peaks that are more shallow, but represent an overall large angle change (integral...)  
%  Algorithm: sum 2 angle changes: For each frame add the value of its two neighbours. Then find peaks...    
vec_sum = filter(b,a,vec);
vec_sum = [vec_sum((delay+1):end), vec_sum((end-delay+1):end)];

[pks1,locs1] = findpeaks(double(vec_sum),'MINPEAKHEIGHT',ThresholdValue);
[pks2,locs2] = findpeaks(double(-vec_sum),'MINPEAKHEIGHT',ThresholdValue);
locs_sum     = [locs1  locs2];
pks_sum      = [pks1  -pks2];
[locs_sum,IX]= sort(locs_sum);
pks_sum      = pks_sum(IX);
         
locs_sum_logical            = false(size(vec));
locs_sum_logical(locs_sum)  = true;
pks_sum_values              = zeros(size(vec));
pks_sum_values(locs_sum)    = pks_sum;

% find the peaks that were found by integration and add them to the pks pool  
locs_sum_NotPreviouslyDiscovered = find(locs_sum_logical .* ~ locs_logical_WithNeighbors);
pks_sum_NotPreviouslyDiscovered  = pks_sum_values(locs_sum_NotPreviouslyDiscovered);

locs_all     = [locs  locs_sum_NotPreviouslyDiscovered];
pks_all      = [pks   pks_sum_NotPreviouslyDiscovered];
[locs_all,IX]= sort(locs_all);
pks_all      = pks_all(IX);

if plotme
    figure('name','Peaks of changes in velocity- with vs. without integration'); 
    plot(vec,'k.-'); hold on; plot(locs,pks,'r*'); 
    plot(vec_sum,'b.:'); hold on; plot(locs_sum,pks_sum,'go'); 

    figure('name','Peaks of changes in velocity- All'); 
    plot(vec,'k.-'); hold on; plot(locs_all,pks_all,'r*'); 
end

return

function BehaviorCode = GenerateBehaviorCode (Logical_Vectors, Logical_Vectors_HighLevel)

NumberOfFrames = length(Logical_Vectors.OutOfBounds);

%% Low level behavior
LowLevel_CodeNumber = [      0             1                  2               3       4         5        6      ];     
LowLevel_CodeName   = {'Out of bounds','Forward',          'Curve',         'Pause','Reverse','Omega','SharpTurn'};
LowLevel_FieldNames = {'OutOfBounds'  ,'Forward_NoCurving','CurvingForward','Pause','Reverse','Omega','SharpTurns'};

TotalFramesCount = 0;
BehaviorVector   = single(zeros(1,NumberOfFrames));
for f_ind = 1:length(LowLevel_FieldNames)
    name                       = LowLevel_FieldNames{f_ind};
    LogicalVec                 = Logical_Vectors.(name);
    BehaviorVector(LogicalVec) = LowLevel_CodeNumber(f_ind);    
    TotalFramesCount           = TotalFramesCount + length(find(LogicalVec));
end
if TotalFramesCount~=NumberOfFrames    % validation. This message should never be displayed
    disp('WARNING: OVERLAPPING SEGMENTS!!!  LOW LEVEL BEHAVIOR IS NOT UNIQUELY SEGMENTED!')
end
BehaviorCode.LowLevel.BehaviorVector      = BehaviorVector;
BehaviorCode.LowLevel.BehaviorCodeName    = LowLevel_CodeName;
BehaviorCode.LowLevel.BehaviorCodeNumbers = LowLevel_CodeNumber;

%% High level behavior
HighLevel_CodeNumber = [      0             1                  2               3       4              5                     6                  7                           8            ];     
HighLevel_CodeName   = {'Out of bounds','Forward',          'Curve',         'Pause','Reversal','Omega-Pause',          'SharpTurn', 'Pirouette. Initial reversal', 'Pirouette. After reversal'};
HighLevel_FieldNames = {'OutOfBounds'  ,'Forward_NoCurving','CurvingForward','Pause','Reversal','OmegaWithoutPirouette','SharpTurns','Pirouette_InitialReversal',   'Pirouette_AfterReversal'};

TotalFramesCount = 0;
BehaviorVector   = single(zeros(1,NumberOfFrames));
for f_ind = 1:length(HighLevel_FieldNames)
    name                       = HighLevel_FieldNames{f_ind};
    LogicalVec                 = Logical_Vectors_HighLevel.(name);
    BehaviorVector(LogicalVec) = HighLevel_CodeNumber(f_ind);    
    TotalFramesCount           = TotalFramesCount + length(find(LogicalVec));
end
if TotalFramesCount~=NumberOfFrames  % validation. This message should never be displayed
    disp('WARNING: OVERLAPPING SEGMENTS!!! HIGH LEVEL BEHAVIOR IS NOT UNIQUELY SEGMENTED!')
end
BehaviorCode.HighLevel.BehaviorVector      = BehaviorVector;
BehaviorCode.HighLevel.BehaviorCodeName    = HighLevel_CodeName;               
BehaviorCode.HighLevel.BehaviorCodeNumbers = HighLevel_CodeNumber;  
BehaviorCode.HighLevel.Pirouette_ExactInitiationTime = Logical_Vectors_HighLevel.Pirouette_ExactInitiationTime; % true for pirouettes frames that the initial pirouette timing is well defined.
BehaviorCode.HighLevel.Pirouette_NotAfterOutOfBound  = Logical_Vectors_HighLevel.Pirouette_NotAfterOutOfBound; % true for pirouettes that are not after out-of-bound frames. Subset of 'Pirouette_ExactInitiationTime'

return

function Add_Xlines(h, Values, ChangeTicks)  % Add horizontal lines only to the current axis.
  
xlimits = get(h,'Xlim');
for V = Values
    line(xlimits,[V V],'color','k','linestyle',':','parent',h); hold on;   
end
if exist('ChangeTicks','var')
    if ChangeTicks
        Ticks = Values;
        set(h,'ytick',Ticks);
    end
end

return
























