function [TrackLocomotionData, CurrentTrack] = AnalyzeLocomotionProperties_SL_HighThroughputDevice_02(CurrentTrack, Settings, Head, Tail, AllowVeryLongTracks)
% This function analyzes animal locomotion
% It will put NaN values whereever there are too many sequential frames missing.   
% The threshold for the maximal amount of missing frames allowed is set by the used in the Settings variable using the Settings function.  

% Input 
% coordinates = CurrentTrack.Path is the XY path of the worm track. column 1 and 2 corresponds to X and Y, respectively.  
% Settings    = parameters, smoothing information, and plot information. Defined by the user using the Settings function.   
% Head, Tail  = optional variables.  
% AllowVeryLongTracks = true/false. This defines the maximal number of calculations allowed when finding continous angles as 1e6 or 1e4, respectively.  

% Output 
% TrackLocomotionData. A structure with the track locomotion information. See full decription of fields below.
%
%
%   September 2018, Sagi Levy
%

%% Initialization
coordinates = single(CurrentTrack.Path);
TrackLength = size(coordinates,1);

% Settings: parameters, smoothing information, and plot information
PlotInfo                = Settings.PlotInfo;                                                          % Information about which plots to show, and whether to show the movie. 
InterpolationFactor     = Settings.SmoothingInfo.InterpolationFactor;                                 % InterpolationFactor    = 0.9; 
DeltaFrameSmoothed      = Settings.SmoothingInfo.DeltaFrameSmoothed;                                  % DeltaFrameSmoothed     = 1;
MaxAllowedMissingFrames = Settings.SmoothingInfo.MaxAllowedTimeForMissingFrames * Settings.FrameRate; % The threshold for the maximal amount of missing frames. 

% Find Long NaN segments where the worm was not properly detected
NaN_Frames       = isnan(coordinates(:,1));
StartNaNSegments = (find(diff(NaN_Frames)==1)+1)'; if NaN_Frames(1)==1,   StartNaNSegments = [1 StartNaNSegments]; end
EndNaNSegments   = (find(diff(NaN_Frames)==-1))';  if NaN_Frames(end)==1, EndNaNSegments   = [EndNaNSegments length(NaN_Frames)]; end
if length(StartNaNSegments) ~= length(EndNaNSegments)
    disp('error in NaN segments calculation in function: AnalyzeLocomotionProperties_SL_vXX');
end

LongNaNSegments_indices = find((EndNaNSegments - StartNaNSegments)>MaxAllowedMissingFrames);
if ~isempty(LongNaNSegments_indices)
    CorrectSmoothedVectors = true;
    FramesWithLongNaNSegments = false(1,length(coordinates));
    for ind = LongNaNSegments_indices
        FramesWithLongNaNSegments(StartNaNSegments(ind):EndNaNSegments(ind)) = true;    
    end
else
    CorrectSmoothedVectors = false;
end

%% Math
% Smooth
RawVectorIndices          = 1: TrackLength; 
SmoothedVectorsIndices    = single(1: DeltaFrameSmoothed: TrackLength);
Smoothed_coordinates      = single(zeros(length(SmoothedVectorsIndices),2));
Smoothed_coordinates(:,1) = single(csaps(RawVectorIndices,double(coordinates(:,1)'), InterpolationFactor, SmoothedVectorsIndices));
Smoothed_coordinates(:,2) = single(csaps(RawVectorIndices,double(coordinates(:,2)'), InterpolationFactor, SmoothedVectorsIndices));
if CorrectSmoothedVectors 
    Smoothed_coordinates(FramesWithLongNaNSegments,:)=NaN;
end

% Differentiate
Displacement                = diff(coordinates,1);                                  % Units = pixels.
VelocityVector              = Displacement;                                         % Units = pixels/frame.
Smoothed_Displacement       = diff(Smoothed_coordinates,1);                         % Units = pixels.
Smoothed_VelocityVector     = Smoothed_Displacement / DeltaFrameSmoothed;           % Units = pixels/frame.
if  CurrentTrack.TrackLength > 2
    AccelerationVector          = diff(VelocityVector);                                 % Units = pixels/frame^2.  
    Smoothed_AccelerationVector = diff(Smoothed_VelocityVector,1) / DeltaFrameSmoothed; % Units = pixels/frame^2.
end

% Velocity in polar coordinates
[Velocity.Angle,          Velocity.Amplitude]         = cart2pol(VelocityVector(:,1)             , VelocityVector(:,2));           % [THETA,RHO] = cart2pol(X,Y);
[Smoothed_Velocity.Angle, Smoothed_Velocity.Amplitude]= cart2pol(Smoothed_VelocityVector(:,1)    , Smoothed_VelocityVector(:,2));  % [THETA,RHO] = cart2pol(X,Y);
Velocity.Angle          = Velocity.Angle          /(2*pi)*360;
Smoothed_Velocity.Angle = Smoothed_Velocity.Angle /(2*pi)*360;
Velocity.ContinuousAngle              = FindContinuousAngle(Velocity.Angle,'degrees',[],false, AllowVeryLongTracks);               % FindContinuousAngle(Velocity.Angle,'degrees',[],true); % plotfigure = true                                                   
Smoothed_Velocity.ContinuousAngle     = FindContinuousAngle(Smoothed_Velocity.Angle,'degrees',[],false, AllowVeryLongTracks);      % FindContinuousAngle(Smoothed_Velocity.Angle,'degrees',[],true); % plotfigure = true                                                   

if  CurrentTrack.TrackLength > 2
    % Acceleration in polar coordinates
    [Acceleration.Angle,          Acceleration.Amplitude]         = cart2pol(AccelerationVector(:,1)             , AccelerationVector(:,2));           % [THETA,RHO] = cart2pol(X,Y);
    [Smoothed_Acceleration.Angle, Smoothed_Acceleration.Amplitude]= cart2pol(Smoothed_AccelerationVector(:,1)    , Smoothed_AccelerationVector(:,2));  % [THETA,RHO] = cart2pol(X,Y);
    Acceleration.Angle          = Acceleration.Angle          /(2*pi)*360;
    Smoothed_Acceleration.Angle = Smoothed_Acceleration.Angle /(2*pi)*360;
end

% Differentiation of the velocity polar coordinates  
if CurrentTrack.TrackLength > 2
    ChangesIn_Velocity.Amplitude          = diff(Velocity.Amplitude);                                    % d(|V(t)|)/dt, units= pixels/Frames^2 
    ChangesIn_Velocity.Angle              = diff(Velocity.ContinuousAngle) ;                             % d( U(t) )/dt, units= angle /Frames^2 
    ChangesIn_Smoothed_Velocity.Amplitude = diff(Smoothed_Velocity.Amplitude)/ DeltaFrameSmoothed;       % d(|V(t)|)/dt, units= pixels/Frames^2 
    ChangesIn_Smoothed_Velocity.Angle     = diff(Smoothed_Velocity.ContinuousAngle)/ DeltaFrameSmoothed; % d( U(t) )/dt, units= angle /Frames^2 
    ChangesIn_Smoothed_Velocity.ContinuousAngle  = FindContinuousAngle(ChangesIn_Smoothed_Velocity.Angle,'degrees',[],false, AllowVeryLongTracks);      % FindContinuousAngle(Smoothed_Velocity.Angle,'degrees',[],true); % plotfigure = true           
end

if CurrentTrack.TrackLength > 3
    Jerk.Amplitude                                  = diff(ChangesIn_Smoothed_Velocity.Amplitude)/ DeltaFrameSmoothed;       % d(|V(t)|)/dt^2, units= pixels/Frames^3 
    Jerk.Angle                                      = diff(ChangesIn_Smoothed_Velocity.ContinuousAngle)/ DeltaFrameSmoothed; % d( U(t) )/dt, units= angle /Frames^2 
end

%% KEEP SIMILAR DIMENSIONS BY ADDING NANS
% Add NaN at the last velocity value, and at the last and first acceleration values
VelocityVector                          = [VelocityVector;              [NaN NaN]];
Smoothed_VelocityVector                 = [Smoothed_VelocityVector;     [NaN NaN]];
Velocity.Angle                          = [Velocity.Angle;                    NaN];
Velocity.Amplitude                      = [Velocity.Amplitude;                NaN];
Velocity.ContinuousAngle                = [Velocity.ContinuousAngle;          NaN];
Smoothed_Velocity.Angle                 = [Smoothed_Velocity.Angle;           NaN];
Smoothed_Velocity.Amplitude             = [Smoothed_Velocity.Amplitude;       NaN];
Smoothed_Velocity.ContinuousAngle       = [Smoothed_Velocity.ContinuousAngle; NaN];

if CurrentTrack.TrackLength > 2
    AccelerationVector                      = [[NaN NaN]; AccelerationVector; [NaN NaN]];
    Acceleration.Angle                      = [NaN; Acceleration.Angle;            NaN];
    Acceleration.Amplitude                  = [NaN; Acceleration.Amplitude;        NaN];
    Smoothed_AccelerationVector             = [[NaN NaN]; Smoothed_AccelerationVector; [NaN NaN]];
    Smoothed_Acceleration.Angle             = [NaN; Smoothed_Acceleration.Angle;            NaN];
    Smoothed_Acceleration.Amplitude         = [NaN; Smoothed_Acceleration.Amplitude;        NaN];
    ChangesIn_Smoothed_Velocity.Angle       = [NaN; ChangesIn_Smoothed_Velocity.Angle;      NaN];
    ChangesIn_Smoothed_Velocity.Amplitude   = [NaN; ChangesIn_Smoothed_Velocity.Amplitude;  NaN];
else
    AccelerationVector                      = NaN*single(zeros(CurrentTrack.TrackLength,2));  
    Acceleration.Angle                      = NaN*single(zeros(CurrentTrack.TrackLength,1)); 
    Acceleration.Amplitude                  = NaN*single(zeros(CurrentTrack.TrackLength,1)); 
    Smoothed_AccelerationVector             = NaN*single(zeros(CurrentTrack.TrackLength,2)); 
    Smoothed_Acceleration.Angle             = NaN*single(zeros(CurrentTrack.TrackLength,1)); 
    Smoothed_Acceleration.Amplitude         = NaN*single(zeros(CurrentTrack.TrackLength,1)); 
    ChangesIn_Smoothed_Velocity.Angle       = NaN*single(zeros(CurrentTrack.TrackLength,1)); 
    ChangesIn_Smoothed_Velocity.Amplitude   = NaN*single(zeros(CurrentTrack.TrackLength,1));     
end
if CurrentTrack.TrackLength > 3
    Jerk.Angle                              = [NaN; NaN; Jerk.Angle;      NaN];
    Jerk.Amplitude                          = [NaN; NaN; Jerk.Amplitude;  NaN];
else
    Jerk.Angle                              = NaN*single(zeros(CurrentTrack.TrackLength,1)); 
    Jerk.Amplitude                          = NaN*single(zeros(CurrentTrack.TrackLength,1));        
end

%% Assign Vectors to TrackLocomotionData Output variable
TrackLocomotionData.FrameIndices            = SmoothedVectorsIndices;                  % 
TrackLocomotionData.Path.X                  = Smoothed_coordinates(:,1);               % Cartesian coordinate X
TrackLocomotionData.Path.Y                  = Smoothed_coordinates(:,2);               % Cartesian coordinate Y
TrackLocomotionData.Velocity                = Smoothed_Velocity;                       % Contains Polar coordinates fields: Amplitude, Angle, Continuous Angle. 
TrackLocomotionData.Velocity.X              = Smoothed_VelocityVector(:,1);            % Cartesian coordinate X
TrackLocomotionData.Velocity.Y              = Smoothed_VelocityVector(:,2);            % Cartesian coordinate Y 
TrackLocomotionData.Acceleration            = Smoothed_Acceleration;                   % Contains Polar coordinates fields: Amplitude, Angle. 
TrackLocomotionData.Acceleration.X          = Smoothed_AccelerationVector(:,1);        % Cartesian coordinate X
TrackLocomotionData.Acceleration.Y          = Smoothed_AccelerationVector(:,2);        % Cartesian coordinate Y 
TrackLocomotionData.ChangeInVelocity        = ChangesIn_Smoothed_Velocity;             % Second derivatives. Contains Polar coordinates fields: Amplitude, Angle. Direct differentiation of polar coordinates. 
TrackLocomotionData.Jerk                    = Jerk;                                    % Third derivatives.  Contains Polar coordinates fields: Amplitude, Angle. Direct differentiation of polar coordinates. 

TrackLocomotionData.SmoothingInfo.InterpolationFactor = InterpolationFactor;
TrackLocomotionData.SmoothingInfo.DeltaFrameSmoothed  = DeltaFrameSmoothed;

%% Plots
PlotTrackLocomotion(TrackLocomotionData, coordinates, PlotInfo, Head, Tail) ;
if PlotInfo.MOVIE_using_ChangesInVelocity                      % Movie of Path
    MovieOfTrackLocomotion(TrackLocomotionData, coordinates) ;
end

return

%% Calculation functions
function ContinuousAngles = FindContinuousAngle(Angles, Mode_in, RangeEdges, plotfigures_in, AllowVeryLongTracks_in)      
% Inputs
% Angles             A vector of angles that may be trancated due to edge singularities (e.g. 180--> -180, etc)
% Mode               Optional string: 'degrees' if the range is [-180 180] (default) or radians if range is between [-pi pi].  
% RangeEdges         Optional vector of length 2. Manually specify the edges and ignore 'Mode'. 
% plotfigures        Optional logical. default = false;

% Output:
% ContinuousAngles   A vector of Angles without singularities, range may vary beyond total of 360 degrees (or 2*pi radians) 

%% find AnglesEdge
Mode = 'degrees'; % [-180 180] 
if exist('Mode_in','var')
    Mode = Mode_in;    
end
plotfigures = false; 
if exist('plotfigures_in','var')
    if ~isempty(plotfigures_in)
        plotfigures = plotfigures_in; 
    end
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
AllowVeryLongTracks = false; 
if exist('AllowVeryLongTracks_in','var')
    if ~isempty(AllowVeryLongTracks_in)
        AllowVeryLongTracks = AllowVeryLongTracks_in; 
    end
end
if AllowVeryLongTracks
    MaxNumberOfCalculationSteps = 1e6;
else
    MaxNumberOfCalculationSteps = 1e4;
end   

%% Find when angles are changes by more than half the range: fix discontinuity by correcting all following angles. Repeat until end of vector
ContinuousAngles = Angles;
AngleLimit       = AnglesEdge(2);
for step = 1:MaxNumberOfCalculationSteps
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

%% Plot and movie functions
function Add_Ygrid_For_Angle_Plots(XLIMITS)

hold on; 
Angles = [-180 -90 0 90 180];
for Angle = Angles
    line([XLIMITS(1) XLIMITS(2)],[Angle Angle],'linestyle',':','color',[0.5 0.5 0.5]);
end
set(gca,'ytick',Angles);
return

function MyHistogram(Vec,XLABEL_STRING)
MinimumNumberOfBins = 10;
NumberOfBins        = max(MinimumNumberOfBins, round(length(Vec)/10));

hist(Vec, NumberOfBins); 
ylabel('Frame Count','fontsize',16);
title(['mean = ',num2str(nanmean(Vec),2),' , STD = ',num2str(nanstd(Vec),2)],'fontsize',16);
set(gca,'fontsize',20);
if exist('XLABEL_STRING','var')
    xlabel(XLABEL_STRING,'fontsize',16);
end    
return

function PlotTrackLocomotion(TrackLocomotionData, coordinates, PlotInfo, Head, Tail) 

SmoothedVectorsIndices          = TrackLocomotionData.FrameIndices;         % 
Smoothed_coordinates(:,1)       = TrackLocomotionData.Path.X;               % Cartesian coordinate X
Smoothed_coordinates(:,2)       = TrackLocomotionData.Path.Y;               % Cartesian coordinate Y
Smoothed_Velocity               = TrackLocomotionData.Velocity;             % Contains Polar coordinates fields: Amplitude, Angle, Continuous Angle. 
Smoothed_VelocityVector(:,1)    = TrackLocomotionData.Velocity.X;           % Cartesian coordinate X
Smoothed_VelocityVector(:,2)    = TrackLocomotionData.Velocity.Y;           % Cartesian coordinate Y 
Smoothed_Acceleration           = TrackLocomotionData.Acceleration;         % Contains Polar coordinates fields: Amplitude, Angle. 
Smoothed_AccelerationVector(:,1)= TrackLocomotionData.Acceleration.X;       % Cartesian coordinate X
Smoothed_AccelerationVector(:,2)= TrackLocomotionData.Acceleration.Y;       % Cartesian coordinate Y 
ChangesIn_Smoothed_Velocity     = TrackLocomotionData.ChangeInVelocity;     % Contains Polar coordinates fields: Amplitude, Angle. Direct differentiation of polar coordinates. 

RawVectorIndices = 1:length(coordinates(:,1));
ScreenSize       = get(0,'screensize');

%% Plots
if PlotInfo.Path   % Path
    figure('name','Path','position',ScreenSize); 
    subplot(1,3,1); plot(RawVectorIndices,coordinates(:,1),'k*');hold on; plot(SmoothedVectorsIndices,Smoothed_coordinates(:,1),'r.-'); title ('X')
    subplot(1,3,2); plot(RawVectorIndices,coordinates(:,2),'k*');hold on; plot(SmoothedVectorsIndices,Smoothed_coordinates(:,2),'r.-'); title ('Y') 
    subplot(1,3,3); plot(coordinates(:,1)         ,coordinates(:,2)         ,'k*',coordinates(1,1),         coordinates(1,2),         'gd'); hold on;
                    plot(Smoothed_coordinates(:,1),Smoothed_coordinates(:,2),'r.-',Smoothed_coordinates(1,1),Smoothed_coordinates(1,2),'go'); 
    xlabel ('X'); ylabel ('Y')
end
if PlotInfo.VelocityCartesian
    figure('name','Velocity- Cartesyen coordinates','position',ScreenSize); 
    subplot(1,3,1); plot(SmoothedVectorsIndices,Smoothed_VelocityVector(:,1),'r.-');  title ('Velocity X')
    subplot(1,3,2); plot(SmoothedVectorsIndices,Smoothed_VelocityVector(:,2),'r.-');  title ('Velocity Y')
    subplot(1,3,3); plot(Smoothed_VelocityVector(:,1) ,Smoothed_VelocityVector(:,2) ,'r.-', Smoothed_VelocityVector(1,1) ,Smoothed_VelocityVector(1,2), 'go'); hold on;
    xlabel ('Velocity X'); ylabel ('Velocity Y')
end
if PlotInfo.VelocityPolar   
    figure('name','Velocity- polar coordinates','position',ScreenSize); 
    subplot(1,3,1); plot(SmoothedVectorsIndices,Smoothed_Velocity.Amplitude,'r.-');  title ('Velocity amplitude')
    subplot(1,3,2); plot(SmoothedVectorsIndices,Smoothed_Velocity.Angle,'r.-');      title ('Velocity angle')
    subplot(1,3,3); plot(Smoothed_Velocity.Amplitude ,Smoothed_Velocity.Angle ,'r.-', Smoothed_Velocity.Amplitude(1) ,Smoothed_Velocity.Angle(1)   ,'go'); hold on;
    xlabel ('Velocity amplitude'); ylabel ('Velocity angle')
end
if PlotInfo.AccelerationCartesian   
    figure('name','Acceleration- Cartesyen coordinates','position',ScreenSize); 
    subplot(1,3,1); plot(SmoothedVectorsIndices,Smoothed_AccelerationVector(:,1),'r.-');  title ('Acceleration X')
    subplot(1,3,2); plot(SmoothedVectorsIndices,Smoothed_AccelerationVector(:,2),'r.-');  title ('Acceleration Y')
    subplot(1,3,3); plot(Smoothed_AccelerationVector(:,1) ,Smoothed_AccelerationVector(:,2) ,'r.-', Smoothed_AccelerationVector(1,1) ,Smoothed_AccelerationVector(1,2), 'go'); hold on;
    xlabel ('Acceleration X'); ylabel ('Acceleration Y')
end
if PlotInfo.AccelerationPolar  
    figure('name','Acceleration','position',ScreenSize); 
    subplot(1,3,1); plot(SmoothedVectorsIndices,Smoothed_Acceleration.Amplitude,'r.-');  title ('Acceleration amplitude')
    subplot(1,3,2); plot(SmoothedVectorsIndices,Smoothed_Acceleration.Angle,'r.-');      title ('Acceleration angle')
    subplot(1,3,3); plot(Smoothed_Acceleration.Amplitude ,Smoothed_Acceleration.Angle ,'r.-', Smoothed_Acceleration.Amplitude(1) ,Smoothed_Acceleration.Angle(1)   ,'go'); hold on;
    xlabel ('Acceleration amplitude'); ylabel ('Acceleration angle')
end
if PlotInfo.ChangesedInVelocityPolar 
    figure('name','Changes in Velocity (Not Acceleration)','position',ScreenSize); 
    subplot(1,3,1); plot(SmoothedVectorsIndices,ChangesIn_Smoothed_Velocity.Amplitude,'r.-');  title ('Changes In Velocity: amplitude')
    subplot(1,3,2); plot(SmoothedVectorsIndices,ChangesIn_Smoothed_Velocity.Angle,'r.-');      title ('Changes In Velocity: angle')
    subplot(1,3,3); plot(ChangesIn_Smoothed_Velocity.Amplitude ,ChangesIn_Smoothed_Velocity.Angle ,'r.-', ChangesIn_Smoothed_Velocity.Amplitude(1) ,ChangesIn_Smoothed_Velocity.Angle(1)   ,'go'); hold on;
    xlabel ('Changes In Velocity: amplitude'); ylabel ('Changes In Velocity: angle')
end
if PlotInfo.Path_Velcoty_Acceleration_AllCartesian % All variables in Cartesian coordinates
    figure('name','Path, Velocity and Acceleration- Cartesyen coordinates','position',ScreenSize); 
    subplot(3,2,1); plot(RawVectorIndices,coordinates(:,1),'k*');hold on; plot(SmoothedVectorsIndices,Smoothed_coordinates(:,1),'r.-'); ylabel ('X'); grid;
    subplot(3,2,3); plot(SmoothedVectorsIndices,Smoothed_VelocityVector(:,1),'r.-');  ylabel ('Velocity X'); grid;
    subplot(3,2,5); plot(SmoothedVectorsIndices,Smoothed_AccelerationVector(:,1),'r.-');  ylabel ('Acceleration X'); grid;
    subplot(3,2,2); plot(RawVectorIndices,coordinates(:,2),'k*');hold on; plot(SmoothedVectorsIndices,Smoothed_coordinates(:,2),'r.-'); ylabel ('Y'); grid;
    subplot(3,2,4); plot(SmoothedVectorsIndices,Smoothed_VelocityVector(:,2),'r.-');  ylabel ('Velocity Y'); grid;
    subplot(3,2,6); plot(SmoothedVectorsIndices,Smoothed_AccelerationVector(:,2),'r.-');  ylabel ('Acceleration Y'); grid;
    % set(get(gcf,'children'),'xtick',0:10:TrackLength)
end
if PlotInfo.Path_Velcoty_Acceleration_Polar  % Path in Cartesian coordinates, velocity and acceleration in Polar
    figure('name','Path (cartesyen coordinates) and Velocity and Acceleration (polar coordinates)','position',ScreenSize); 
    subplot(3,2,1); plot(RawVectorIndices,coordinates(:,1),'k*');hold on; plot(SmoothedVectorsIndices,Smoothed_coordinates(:,1),'r.-'); ylabel ('X'); grid;
    subplot(3,2,3); plot(SmoothedVectorsIndices,Smoothed_Velocity.Amplitude,'r.-');  ylabel ('Velocity Amplitude'); grid;
                    LIMITS.VelocityAmplitudeX = get(gca,'xlim');  LIMITS.VelocityAmplitudeY = get(gca,'ylim'); 
    subplot(3,2,5); plot(SmoothedVectorsIndices,Smoothed_Acceleration.Amplitude,'r.-');   ylabel ('Acceleration Amplitude'); grid;
                    LIMITS.AccelerationAmplitudeX = get(gca,'xlim');  LIMITS.AccelerationAmplitudeY = get(gca,'ylim'); 
    subplot(3,2,2); plot(RawVectorIndices,coordinates(:,2),'k*');hold on; plot(SmoothedVectorsIndices,Smoothed_coordinates(:,2),'r.-'); ylabel ('Y'); grid;
    subplot(3,2,4); plot(SmoothedVectorsIndices,Smoothed_Velocity.Angle,'r.-');      ylabel ('Velocity angle'); grid;
                    LIMITS.VelocityAngleX = get(gca,'xlim');  LIMITS.VelocityAngleY = get(gca,'ylim'); 
                    Add_Ygrid_For_Angle_Plots(LIMITS.VelocityAngleX); hold off;
    subplot(3,2,6); plot(SmoothedVectorsIndices,Smoothed_Acceleration.Angle,'r.-');  ylabel ('Acceleration angle'); grid;
                    LIMITS.AccelerationAngleX = get(gca,'xlim');  LIMITS.AccelerationAngleY = get(gca,'ylim'); 
                    Add_Ygrid_For_Angle_Plots(LIMITS.AccelerationAngleX); hold off;
    % set(get(gcf,'children'),'xtick',0:10:TrackLength)
end
if PlotInfo.Path_Velcoty_ChangesInVelocity_Polar % Path in Cartesian coordinates, velocity and changes in velocity in Polar coordinates   
    figure('name','Path (cartesyen coordinates) and Velocity and Changes in velocity (NOT ACCELERATION, polar coordinates)','position',ScreenSize); 
    subplot(3,2,1); plot(RawVectorIndices,coordinates(:,1),'k*');hold on; plot(SmoothedVectorsIndices,Smoothed_coordinates(:,1),'r.-'); ylabel ('X'); grid;
    subplot(3,2,3); plot(SmoothedVectorsIndices,Smoothed_Velocity.Amplitude,'r.-');  ylabel ('Velocity Amplitude'); grid;
                    LIMITS.VelocityAmplitudeX = get(gca,'xlim');  LIMITS.VelocityAmplitudeY = get(gca,'ylim'); 
    subplot(3,2,5); plot(SmoothedVectorsIndices,ChangesIn_Smoothed_Velocity.Amplitude,'r.-');   ylabel ('Changes In Velocity Amplitude'); grid;
                    LIMITS.ChangesInVelocityAmplitudeX = get(gca,'xlim');  LIMITS.ChangesInVelocityAmplitudeY = get(gca,'ylim'); 
    subplot(3,2,2); plot(RawVectorIndices,coordinates(:,2),'k*');hold on; plot(SmoothedVectorsIndices,Smoothed_coordinates(:,2),'r.-'); ylabel ('Y'); grid;
    subplot(3,2,4); plot(SmoothedVectorsIndices,Smoothed_Velocity.ContinuousAngle,'r.-');  ylabel ('Velocity CONTINUOUS angle'); grid;
                    LIMITS.VelocityAngleX = get(gca,'xlim');  LIMITS.VelocityContinuousAngleY = get(gca,'ylim'); 
    subplot(3,2,6); plot(SmoothedVectorsIndices,ChangesIn_Smoothed_Velocity.Angle,'r.-');  ylabel ('Changes In Velocity Angle'); grid;
                    LIMITS.ChangesInVelocityAngleX = get(gca,'xlim');  LIMITS.ChangesInVelocityAngleY = get(gca,'ylim'); 
                    Add_Ygrid_For_Angle_Plots(LIMITS.ChangesInVelocityAngleX); hold off;
    % set(get(gcf,'children'),'xtick',0:10:TrackLength)
end
if PlotInfo.histograms    % Histograms  
    figure('name','Histograms','position',ScreenSize); 
    subplot(2,2,1); MyHistogram(Smoothed_Velocity.Amplitude,          'Velocity Amplitude [pixels/frame]');
    subplot(2,2,2); MyHistogram(Smoothed_Velocity.Angle,              'Velocity Angle [deg/frame]');
                    set(gca,'xtick',[-180 -90 0 90 180]);
    subplot(2,2,3); MyHistogram(ChangesIn_Smoothed_Velocity.Amplitude,'Changes in Velocity Amplitude [pixels/(frame^2)]');  
    subplot(2,2,4); MyHistogram(ChangesIn_Smoothed_Velocity.Angle,    'Changes in Velocity Angle [deg/(frame^2)]');        % turning direction !!
                    set(gca,'xtick',[-180 -90 0 90 180])
end
if PlotInfo.Velocity_vs_ChangesInVelocity   % Plots NOT TIME KINETICS   
    figure('name','Plot combinations of Velocity and Changes in Velocity, Amplitude and Angles','position',ScreenSize); 
    subplot(3,2,1); plot(Smoothed_Velocity.Amplitude,           Smoothed_Velocity.Angle,                'r.');   xlabel ('Velocity Amplitude');            ylabel ('Velocity Angle');                  grid; 
                    set(gca,'ytick',[-180 -90 0 90 180])
    subplot(3,2,2); plot(ChangesIn_Smoothed_Velocity.Amplitude, ChangesIn_Smoothed_Velocity.Angle,      'r.');   xlabel ('Changes in Velocity Amplitude'); ylabel ('Changes in Velocity Angle');       grid;
                    set(gca,'ytick',[-180 -90 0 90 180])
    subplot(3,2,3); plot(Smoothed_Velocity.Amplitude,           ChangesIn_Smoothed_Velocity.Amplitude,  'r.');   xlabel ('Velocity Amplitude');            ylabel ('Changes in Velocity Amplitude');   grid;
    subplot(3,2,4); plot(Smoothed_Velocity.Amplitude,           ChangesIn_Smoothed_Velocity.Angle,      'r.');   xlabel ('Velocity Amplitude');            ylabel ('Changes in Velocity Angle');       grid;
                    set(gca,'ytick',[-180 -90 0 90 180])
    subplot(3,2,5); plot(Smoothed_Velocity.Angle,               ChangesIn_Smoothed_Velocity.Amplitude,  'r.');   xlabel ('Velocity Angle');                ylabel ('Changes in Velocity Amplitude');   grid;
                    set(gca,'xtick',[-180 -90 0 90 180])
    subplot(3,2,6); plot(Smoothed_Velocity.Angle,               ChangesIn_Smoothed_Velocity.Angle,      'r.');   xlabel ('Velocity Angle');                ylabel ('Changes in Velocity Angle');       grid;
                    set(gca,'xtick',[-180 -90 0 90 180],'ytick',[-180 -90 0 90 180])
end

if PlotInfo.HeadTailPath
    if exist('Tail','var')
        figure('name','Head, Tail and Path','position',ScreenSize); 
        subplot(2,1,1); 
            plot(RawVectorIndices,coordinates(:,1),'k*');hold on;   
            plot(RawVectorIndices,Head(:,1),'b*');                 
            plot(RawVectorIndices,Tail(:,1),'r*');                  
            ylabel ('X');
            legend('Centroid','Head','Tail');
        subplot(2,1,2); 
            plot(RawVectorIndices,coordinates(:,2),'k*');hold on; 
            plot(RawVectorIndices,Head(:,2),'b*');                
            plot(RawVectorIndices,Tail(:,2),'r*');                
            ylabel ('Y') 
            xlabel ('Time [frames]'); 
            legend('Centroid','Head','Tail');
    end    
end

return

%% Movies functions
function MovieOfTrackLocomotion(TrackLocomotionData, coordinates) 

SmoothedVectorsIndices          = TrackLocomotionData.FrameIndices;         % 
Smoothed_coordinates(:,1)       = TrackLocomotionData.Path.X;               % Cartesian coordinate X
Smoothed_coordinates(:,2)       = TrackLocomotionData.Path.Y;               % Cartesian coordinate Y
Smoothed_Velocity               = TrackLocomotionData.Velocity;             % Contains Polar coordinates fields: Amplitude, Angle, Continuous Angle. 
ChangesIn_Smoothed_Velocity     = TrackLocomotionData.ChangeInVelocity;     % Contains Polar coordinates fields: Amplitude, Angle. Direct differentiation of polar coordinates. 

DeltaFrameSmoothed              = TrackLocomotionData.SmoothingInfo.DeltaFrameSmoothed;

%%  Movie of Path: Using the derivative of the velocity angle (not Acceleration angle per se)    
ScaleAxisPositionToPixelSize    = true;
FramesPerSecond                 = 3;
% MOVIE_MODE                    = 'OnlyGenerateFrames';
% MOVIE_MODE                    = 'RunMovie';
MOVIE_MODE                      = 'CreateMovieFile';                                                  
AddAngleGrid                    = false;
AddToMovieTitle                 = 'Change in velocity = derivative of Velocity angle and of velocti amplitude (Not Acceleration)';

MinThresholdForTurning          = 90;  % degrees/frames^2
MaxThresholdForPausing          = 0.77; % pixel/frame   % NOTE!! Pauses in Dirk's code is defined as [ V < (0.022 mm/s) ~ (0.77 pixels/frame)].
High_ChangesInVelocityAngle     = find(abs(ChangesIn_Smoothed_Velocity.Angle) > MinThresholdForTurning);
Low_VelocityAmplitude           = find(    Smoothed_Velocity.Amplitude        < MaxThresholdForPausing);
BOTH_HighAngleChange_and_LowVelocityAmplitude = intersect(High_ChangesInVelocityAngle, Low_VelocityAmplitude);
High_ChangesInVelocityAngle     = setdiff(High_ChangesInVelocityAngle, BOTH_HighAngleChange_and_LowVelocityAmplitude);
Pause                           = setdiff(Low_VelocityAmplitude,       BOTH_HighAngleChange_and_LowVelocityAmplitude);  % Pause = Low_VelocityAmplitude that is NOT a turn

AddPauses.Condition1.Indices    = High_ChangesInVelocityAngle;  
AddPauses.Condition1.color      = 'm';   
AddPauses.Condition1.PauseTime  = 0;   % Frames to pause in the output movie  
AddPauses.Condition2.Indices    = Pause;  
AddPauses.Condition2.color      = 'y';   
AddPauses.Condition2.PauseTime  = 0;   % Frames to pause in the output movie  
AddPauses.Condition3.Indices    = BOTH_HighAngleChange_and_LowVelocityAmplitude;  
AddPauses.Condition3.color      = 'r';   
AddPauses.Condition3.PauseTime  = 0;   % Frames to pause in the output movie  

CreateTrackMovie(coordinates, Smoothed_coordinates, DeltaFrameSmoothed, SmoothedVectorsIndices, ...
                  Smoothed_Velocity.Amplitude, Smoothed_Velocity.ContinuousAngle, ChangesIn_Smoothed_Velocity.Amplitude, ChangesIn_Smoothed_Velocity.Angle, ...
                  AddAngleGrid, ScaleAxisPositionToPixelSize, FramesPerSecond, MOVIE_MODE, AddToMovieTitle, AddPauses);

return

function [MOV, FileName] = CreateTrackMovie(coordinates, Smoothed_coordinates, DeltaFrameSmoothed, ...
                                            SmoothedVectorsIndices, VelocityAmplitude, VelocityAngle, ChangeInVelocityAmplitude, ChangeInVelocityAngle, ...
                                            AddAngleGrid, ScaleAxisPositionToPixelSize, FramesPerSecond, MOVIE_MODE, AddToMovieTitle, AddPauses)
%% NOTE !!
% The input Xaxis will be VectorsIndices or SmoothedVectorsIndices depending if the Y data is smoothed or not (see below)  
% The input VelocityAmplitude may refer to Velocity.Amplitude (Raw data) or Smoothed_Velocity.Amplitude , etc....  
% The input VelocityAngle     may refer to Velocity.Angle     (Raw data) or Smoothed_Velocity.Angle      or Smoothed_Velocity.ContinuousAngle , etc....  
% The same for the acceleration inputs.

%% Initialization
FileName     = [];
NumOfFrames  = length(coordinates(:,1));

if ~exist('ScaleAxisPositionToPixelSize','var')
    ScaleAxisPositionToPixelSize = true;
end
if ~exist('FramesPerSecond','var')
    FramesPerSecond = 5;
end
if ~exist('AddAngleGrid','var')
    AddAngleGrid = false;    
end
if ~exist('MOVIE_MODE','var')
    MOVIE_MODE = 'OnlyGenerateFrames';
%     MOVIE_MODE = 'RunMovie';
%     MOVIE_MODE = 'CreateMovieFile';
end
if ~exist('AddToMovieTitle','var')
    MovieName = 'Track Analysis Movie';
else
    MovieName = ['Track Analysis Movie. ',AddToMovieTitle];
end

% Add pauses in interesting time points
AddPausesToMovie = false; 
if exist('AddPauses','var')
    AddPausesToMovie    = true;
    PauseConditionNames = fieldnames(AddPauses);
    PauseFramesInMovie  = zeros(1,NumOfFrames); % 0 if no pause, 1 for condition #1, 2 for condition #2, ...
    AdditionalFrames    = 0;
    for cond_num = 1:length(PauseConditionNames)
        Indices                     = AddPauses.(PauseConditionNames{cond_num}).Indices;
        PauseFramesInMovie(Indices) = cond_num;
        AdditionalFrames            = AdditionalFrames + length(Indices)* (AddPauses.(PauseConditionNames{cond_num}).PauseTime - 1);
    end
end

% Movie structure
Movie_Figure_handle   = figure('name',MovieName,'position',get(0,'screensize')); 
if ~AddPausesToMovie
    MOV(NumOfFrames)  = struct('cdata',[],'colormap',[]);        % Here I will load frames to the movie
else
    MovieLength       = NumOfFrames + AdditionalFrames;
    MOV(MovieLength)  = struct('cdata',[],'colormap',[]);        % Here I will load frames to the movie
end    

% Position of 'Path' plot if it needs to be rescales (if ScaleAxisPositionToPixelSize==true )
XYposition            = [0.0740    0.09    0.3907    0.8150];    % position of the 'path' plot, before rescaling it.
subplot(1,2,1); plot(coordinates(:,1),coordinates(:,2)); hold off; 
LIMITS.X              = get(gca,'xlim'); 
LIMITS.Y              = get(gca,'ylim'); 
LIMITS.YtoX_AxisRatio = diff(LIMITS.Y) / diff(LIMITS.X);    

if LIMITS.YtoX_AxisRatio > XYposition(4)/XYposition(3)  % [Changes in Y] are larger than [changes in X] and their ratio is larger than the maximal axis ratio. 
    % Shrink X axis to keep pixel size
    XYposition(3)= XYposition(4)/ LIMITS.YtoX_AxisRatio;
    XYposition(1)= XYposition(1) + (0.5 - XYposition(3))/2;    % Centralize the graph     
else  % [Changes in Y] are larger than [changes in X] but their ratio is smaller than the maximal axis ratio. 
    % Shrink Y axis to keep pixel size
    XYposition(4)= XYposition(3)* LIMITS.YtoX_AxisRatio; 
    XYposition(2)= 0.5 - XYposition(4)/2;    % Centralize the graph 
end

% Generate Movie File (if requested)
if strcmpi(MOVIE_MODE,'CreateMovieFile')
    writerObj           = VideoWriter([MovieName, '.avi']);
    writerObj.FrameRate = FramesPerSecond; 
    open(writerObj);
    FileName = writerObj.Filename;
    
end

% Get vector limits
LIMITS.subplot.X = [0 round(length(SmoothedVectorsIndices)*1.06)]; % add ~6% ; 
figure; 
plot(SmoothedVectorsIndices,VelocityAmplitude);          LIMITS.subplot.VelocityAmplitude         = get(gca,'ylim'); 
plot(SmoothedVectorsIndices,VelocityAngle);              LIMITS.subplot.VelocityAngle             = get(gca,'ylim'); 
plot(SmoothedVectorsIndices,ChangeInVelocityAmplitude);  LIMITS.subplot.ChangeInVelocityAmplitude = get(gca,'ylim'); 
plot(SmoothedVectorsIndices,ChangeInVelocityAngle);      LIMITS.subplot.ChangeInVelocityAngle     = get(gca,'ylim'); 
close;

%% Plot Prames
Movie_frame = 0;
for time = 1:NumOfFrames
    Movie_frame = Movie_frame+1;
    CurrentSmoothIndex = 1 + (time-1)/DeltaFrameSmoothed;    
    Indices.Raw    = 1:time;
    Indices.Smooth = 1:CurrentSmoothIndex;
    
    subplot(1,2,1);  % PATH
        plot(coordinates(Indices.Raw,1),                coordinates(Indices.Raw,2),                '*','color',[0.7 0.7 0.7]); hold on; 
        plot(Smoothed_coordinates(Indices.Smooth,1),    Smoothed_coordinates(Indices.Smooth,2)    ,'k.-','markersize',2); hold on;
        PathEndPointHandle = plot(Smoothed_coordinates(CurrentSmoothIndex,1),Smoothed_coordinates(CurrentSmoothIndex,2),'go');  hold off;
        title(['Track path. Frame =', num2str(time)], 'fontsize',20); 
        xlabel('X position', 'fontsize',20); 
        ylabel('Y position', 'fontsize',20);     
        set(gca,'xlim',LIMITS.X,'ylim',LIMITS.Y);
        if ScaleAxisPositionToPixelSize
            set(gca,'position',XYposition)
        end
    
    subplot(4,2,2);  
        plot(SmoothedVectorsIndices(Indices.Smooth),    VelocityAmplitude(Indices.Smooth)    ,'-','color',[0.33 0.33 1]); hold on;
        plot(SmoothedVectorsIndices(CurrentSmoothIndex),VelocityAmplitude(CurrentSmoothIndex),'bo');                      hold off;
        ylabel('Velocity Amplitude', 'fontsize',16);     
        set(gca,'xlim',LIMITS.subplot.X,'ylim',LIMITS.subplot.VelocityAmplitude);    
    subplot(4,2,4);  
        plot(SmoothedVectorsIndices(Indices.Smooth),    ChangeInVelocityAmplitude(Indices.Smooth)    ,'-','color',[0.33 0.33 1]); hold on;
        plot(SmoothedVectorsIndices(CurrentSmoothIndex),ChangeInVelocityAmplitude(CurrentSmoothIndex),'bo');                      hold off;
        ylabel(['Change in',char(10),'Velocity Amplitude'], 'fontsize',16);  
        set(gca,'xlim',LIMITS.subplot.X,'ylim',LIMITS.subplot.ChangeInVelocityAmplitude);
    subplot(4,2,6);      
        plot(SmoothedVectorsIndices(Indices.Smooth),    VelocityAngle(Indices.Smooth)    ,'-','color',[1 0.33 0.33]); hold on;
        plot(SmoothedVectorsIndices(CurrentSmoothIndex),VelocityAngle(CurrentSmoothIndex),'ro');                      hold off;
        if AddAngleGrid
            Add_Ygrid_For_Angle_Plots(LIMITS.subplot.X); hold off;
        end
        ylabel('Velocity Angle',     'fontsize',16); 
        set(gca,'xlim',LIMITS.subplot.X,'ylim',LIMITS.subplot.VelocityAngle);        
    subplot(4,2,8);      
        plot(SmoothedVectorsIndices(Indices.Smooth),    ChangeInVelocityAngle(Indices.Smooth)    ,'-','color',[1 0.33 0.33]); hold on;
        plot(SmoothedVectorsIndices(CurrentSmoothIndex),ChangeInVelocityAngle(CurrentSmoothIndex),'ro');                      hold off;    
        Add_Ygrid_For_Angle_Plots(LIMITS.subplot.X); hold off;          % 'Acceleration' Angle will always be shown with [-180 -90 0 90 180] grid    
        ylabel(['Change in',char(10),'Velocity Angle'],     'fontsize',16); 
        set(gca,'xlim',LIMITS.subplot.X,'ylim',LIMITS.subplot.ChangeInVelocityAngle);
    
    % PAUSE IF NEEDED
    if AddPausesToMovie             
        if PauseFramesInMovie(time)               % Add pauses only in the relevant frames
            cond_num  = PauseFramesInMovie(time);
            COLOR     = AddPauses.(PauseConditionNames{cond_num}).color;
            PauseTime = AddPauses.(PauseConditionNames{cond_num}).PauseTime - 1;  % One frame will anyway be added
            % Fixing the plot 
            set(PathEndPointHandle,'MarkerFaceColor',COLOR,'MarkerEdgeColor',COLOR);%,'markersize',12);
            
            for FramesToAdd = 1:PauseTime   % PAUSE: Insert several identical frames to the movie.
                MOV(Movie_frame)= getframe(Movie_Figure_handle);
                if strcmpi(MOVIE_MODE,'CreateMovieFile')
                    writeVideo(writerObj,MOV(Movie_frame))
                end
                Movie_frame = Movie_frame+1;
            end            
        end     
    end
    % Add the frame to the movie
    MOV(Movie_frame)= getframe(Movie_Figure_handle);
    if strcmpi(MOVIE_MODE,'CreateMovieFile')
        writeVideo(writerObj,MOV(Movie_frame))
    end

end

% Run Movie or close Movie file
if strcmpi(MOVIE_MODE,'RunMovie')
    Movie_Figure_handle = figure('name',MovieName,'position',get(0,'screensize')); 
    movie(Movie_Figure_handle,MOV,1,FramesPerSecond);
elseif strcmpi(MOVIE_MODE,'CreateMovieFile')
    close(writerObj);
end

return 









