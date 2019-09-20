function [Tracks, deltaFOverF_Filtered, NeuronCTX_index] = ExtractNeuralActivityFromTracks_AfterGUI_02(Tracks, File, plotme)

%% Notes
% 1) This function is called by GUI function: AnimalAndNeuronTrackerGUI_02.m, but can also run independently. 
%    To run this function after GUI:  load('MoviePath\MovieName_AfterNeuronPositionManualCorrection.mat','Tracks','File') 
% 2) The following fields are NOT USED IN FOLLOWING FUNCTIONS and will be redefined more carefully:     
%       Tracks(tr).Neuron.deltaFOverF      
%       Tracks(tr).Neuron.deltaFOverF_Interpolated 
%       Tracks(tr).Neuron.deltaFOverF_Filtered 
%    These fields are defined here just for allowing fast data plot

%% Free parameter - Please adjust it based on the neuron of interest, imaging quality and magnification
%  NumOfPixelsInNeuron is set to 4 for AWCon imaging line with 4x objective and 0.63 demagnification. 
%  Note that results are relatively robest to this parameter, e.g. very similar results are obtained when set to 6 and 8  
NumOfPixelsInNeuron = 4;

%% Initialization
% prctile_for_F0_definition = 20;   
prctile_for_F0_definition = 3;  % 
if ~exist('plotme','var')
    plotme = false;
end
DisplayMatrixSize             = size(Tracks(1).Neuron.DisplayMatrix,2);
MidPixel                      = round((DisplayMatrixSize+1)/2);
NumOfFrames                   = length(Tracks(1).Neuron.CoordinatesMatrix(:,1));
NumOfTracks                   = length(Tracks);
NeuronYCoordinate             = zeros(NumOfTracks,NumOfFrames,'single')*NaN;
NeuronValues_PixelSortedByAmp = zeros(NumOfTracks,NumOfFrames,DisplayMatrixSize,'single')*NaN;
NeuronValue                   = zeros(NumOfTracks,NumOfFrames,'single')*NaN;

%% Compute NeuronValue and NeuronYCoordinate
for tr=1:NumOfTracks
    NeuronYCoordinate(tr,:)        = Tracks(tr).Neuron.CoordinatesMatrix(:,1); % gradient axis
    FiveOnFiveDisplayMAT           = squeeze(Tracks(tr).Neuron.DisplayMatrix(:,(MidPixel-2):(MidPixel+2),(MidPixel-2):(MidPixel+2))); % (frame,1:5,1:5)
    FiveOnFiveDisplayMATForSorting = zeros(NumOfFrames,DisplayMatrixSize,'single')*NaN;
    for frame=1:NumOfFrames
        SmallMAT = squeeze(FiveOnFiveDisplayMAT(frame,:,:));
        FiveOnFiveDisplayMATForSorting(frame,:) = SmallMAT(:);        
    end    
    NaN_indices                                 = all(isnan(FiveOnFiveDisplayMATForSorting),2);
    FiveOnFiveDisplayMATForSorting(isnan(FiveOnFiveDisplayMATForSorting))=0;
    
    NeuronValues_PixelSortedByAmp(tr,:,:) = sort(FiveOnFiveDisplayMATForSorting,2,'descend');
    NeuronValue(tr,:)                     = mean(squeeze(NeuronValues_PixelSortedByAmp(tr,:,1:NumOfPixelsInNeuron)),2);
    NeuronValue(tr,NaN_indices)           = NaN;
    NeuronYCoordinate(tr,NaN_indices)     = NaN;
end
% CTX_index
minY            = min(NeuronYCoordinate(:));
maxY            = max(NeuronYCoordinate(:));
NeuronCTX_index = ((NeuronYCoordinate-minY)/(maxY-minY)*2)-1;
NeuronCTX_index = -NeuronCTX_index;  % 1==up==butanone, -1==down==pentanedione
% deltaFOverF
F0          = prctile(NeuronValue,prctile_for_F0_definition,2); 
F0_MAT      = repmat(F0,[1 NumOfFrames]);
deltaFOverF = (NeuronValue - F0_MAT)./F0_MAT;

%% Find all frames in deltaFOverF that are not within a LONG NaN segment
% Long will be considered above 0.25 sec
% For FrameRate of 30Hz --> the maximal allowed NaN segment is 7 frames 
MaxFramesWithNaN = floor(0.25*File.FrameRate);
NaNsMatrix       = isnan(deltaFOverF);
NaNSegmentStartsOrEnds = diff( [[0;0],NaNsMatrix] ,1,2);
AllowedFrames    = true(size(deltaFOverF));
InterpolatedDeltaFOverF = deltaFOverF; % initialization

for tr=1:NumOfTracks
    NaNSegmentsFirstFrames = find(NaNSegmentStartsOrEnds(tr,:)==1);
    NaNSegmentsLastFrames  = find(NaNSegmentStartsOrEnds(tr,:)==-1)-1;
    if NaNsMatrix(tr,end)==1
        NaNSegmentsLastFrames = [NaNSegmentsLastFrames, size(deltaFOverF,2)]; 
    end
    NaNSegmentLength = NaNSegmentsLastFrames-NaNSegmentsFirstFrames+1;
    NonAllowedSegments = find(NaNSegmentLength>MaxFramesWithNaN) ;
    NonAllowedFrames = [];
    for seg_ind = NonAllowedSegments
        NonAllowedFrames = [NonAllowedFrames, NaNSegmentsFirstFrames(seg_ind):NaNSegmentsLastFrames(seg_ind)];
    end
    AllowedFrames(tr,NonAllowedFrames) = false;
    X = find(~NaNsMatrix(tr,:));
    InterpolatedDeltaFOverF(tr,:) = interp1(X, deltaFOverF(tr,X), 1:size(deltaFOverF,2));
end
InterpolatedDeltaFOverF(~AllowedFrames)= NaN;

%% low pass filter
LPS = designfilt('lowpassfir', ...
  'PassbandFrequency',3,'StopbandFrequency',5, ...
  'PassbandRipple',1,'StopbandAttenuation',80, ...
  'DesignMethod','equiripple','SampleRate',30);
D = round(mean(grpdelay(LPS)));  
deltaFOverF_Filtered = zeros(size(InterpolatedDeltaFOverF),'single');
for tr=1:NumOfTracks
    y = filter(LPS,[InterpolatedDeltaFOverF(tr,:)'; zeros(D,1)]);y = y(D+1:end);
    deltaFOverF_Filtered(tr,:) = y;
end

%% Assign To Tracks
for tr=1:NumOfTracks
    Tracks(tr).Neuron.Value                = NeuronValue(tr,:);
    Tracks(tr).Neuron.F0                   = F0(tr);
    Tracks(tr).Neuron.deltaFOverF          = deltaFOverF(tr,:);
    Tracks(tr).Neuron.deltaFOverF_Interpolated = InterpolatedDeltaFOverF(tr,:);
    Tracks(tr).Neuron.deltaFOverF_Filtered = deltaFOverF_Filtered(tr,:);
    Tracks(tr).Neuron.CTX_index            = NeuronCTX_index(tr,:);       % [-1 1] over the gradient axis
end

%% Plots
if plotme
    TrackColors = {'g','b','m'};    
    % CTX_index_ToPlot = NeuronCTX_index;        % -1 to 1
    CTX_index_ToPlot   = NeuronCTX_index+1;      % 0 to 2  >>  Easier for visualization
    % TimeToPlot       = (1:NumOfFrames);  xlabelstring = 'Frames';
    TimeToPlot         = (1:NumOfFrames)/30;   xlabelstring = 'Time [sec]';
        
    ind=1;
    f(ind) = figure('name','Raw values'); ind=ind+1;
    for tr=1:NumOfTracks
        subplot(NumOfTracks,1,tr); 
        plot(TimeToPlot,NeuronValue(tr,:),'.','color',TrackColors{tr}); hold on;
        xlim([0 TimeToPlot(end)]);
        ylabel('Raw fluorescence');
    end
    xlabel(xlabelstring); 
    
    f(ind) = figure('name','deltaF/F'); ind=ind+1;
    for tr=1:NumOfTracks
        subplot(NumOfTracks,1,tr); 
        plot(TimeToPlot,deltaFOverF(tr,:),'.','color',TrackColors{tr}); hold on;
        plot(TimeToPlot,CTX_index_ToPlot(tr,:),'-','color','k'); hold on;   
        xlim([0 TimeToPlot(end)])
        ylabel('\DeltaF/F');
    end
    xlabel(xlabelstring); 

    f(ind) = figure('name','filtered vs. unfiltered deltaF/F'); ind=ind+1;
    for tr=1:NumOfTracks
        subplot(NumOfTracks,1,tr); 
        plot(TimeToPlot,deltaFOverF(tr,:),'.','color',TrackColors{tr}); hold on;
        plot(TimeToPlot,deltaFOverF_Filtered(tr,:),'-','color',TrackColors{tr}); hold on;
        plot(TimeToPlot,CTX_index_ToPlot(tr,:),'-','color','k'); hold on;   
        xlim([0 TimeToPlot(end)])
        ylim([0 5])
        ylabel('\DeltaF/F');
    end
    xlabel(xlabelstring); 

    f(ind) = figure('name','filtered deltaF/F'); ind=ind+1;
    for tr=1:NumOfTracks
        subplot(NumOfTracks,1,tr); 
        plot(TimeToPlot,deltaFOverF_Filtered(tr,:),'.','color',TrackColors{tr}); hold on;
        plot(TimeToPlot,CTX_index_ToPlot(tr,:),'-','color','k'); hold on;   
        xlim([0 TimeToPlot(end)])
        ylim([0 5])
        ylabel('\DeltaF/F');
    end
    xlabel(xlabelstring); 

    set(f,'position',get(0,'ScreenSize'));
    
    fvtool(LPS)   
end
    
return


