function [Dendrite, Axon, Soma] = ImmobilizedWorm_ProcessTracker(Soma, TrackAxon, TrackDendrite, Frame_For_ROI_Definition, AxonWidthWithBackground, DendriteWidthWithBackground)

%% This function extract basic fluorescence pattern of processes, plots results, and run a short movie to evaluate the segmentation performance 
%
%% Input descriptions:  
%  Soma                         Output Structure generated by the function 'ImmobilizedWorm_SomaTracker.m'  
%  TrackAxon, TrackDendrite     true/false for calculating Axon and dendrite fluorescence           
%  Frame_For_ROI_Definition     Movie frame for visualization during user interface process definition  
%  AxonWidthWithBackground      Width larger than the axon, which includes pixels for axon background calculation      
%  DendriteWidthWithBackground  Width larger than the dendrite, which includes pixels for dendrite background calculation    
%
%% Example inputs
%  TrackAxon = true; TrackDendrite = true; Frame_For_ROI_Definition = 1;

%% Free parameters. These choices fit imaging with 100x objective and demagnification of 0.63x
if ~exist('AxonWidthWithBackground','var') || isempty(AxonWidthWithBackground)
    AxonWidthWithBackground     = 13;%7;
end
if ~exist('DendriteWidthWithBackground','var') || isempty(DendriteWidthWithBackground)
    DendriteWidthWithBackground     = 13;%7;
end
ProcessWidth                = 2;
NumOfPixelsForZoom_Dendrite = 280;
NumOfPixelsForZoom_Axon     = 100;

%% Get movie information and soma tracking information
PathName          = Soma.MovieInfo.PathName;
FrameNames        = Soma.MovieInfo.FrameNames;
NeuronCoordinates = Soma.NeuronCoordinates;
MovieLength       = size(NeuronCoordinates,1);
FractionOfDendritePixels = ProcessWidth/DendriteWidthWithBackground;
FractionOfAxonPixels     = ProcessWidth/AxonWidthWithBackground;

%% Initialization with default parmeters, if needed 
if ~exist('TrackAxon','var') 
    TrackAxon     = true; 
end
if ~exist('TrackDendrite','var') 
    TrackDendrite = true; 
end
if ~exist('Frame_For_ROI_Definition','var') 
    Frame_For_ROI_Definition = 1; 
end

%% Define process lines and compute polygons
SomaReferenceCoordinates                        = NeuronCoordinates(Frame_For_ROI_Definition,:);
SomaCoordinatesChangeRelativeToReferenceFrame   = NeuronCoordinates - repmat(SomaReferenceCoordinates,MovieLength,1);

frame     = Frame_For_ROI_Definition;
FileName  = [PathName,'\',FrameNames{frame}];
f         = imread(FileName);    

if TrackDendrite
    TitleStr                       = 'dendrite';
    Dendrite                       = DefineProcessPolygon(TitleStr,f,DendriteWidthWithBackground, SomaCoordinatesChangeRelativeToReferenceFrame, NumOfPixelsForZoom_Dendrite);  
    DendriteMask                   = Dendrite.ProcessMask;
    Dendrite.NumOfPixels           = round(length(find(DendriteMask(:)))* FractionOfDendritePixels);    
    Dendrite.MovieInfo.PathName    = PathName;
    Dendrite.MovieInfo.FrameNames  = FrameNames;
else
    Dendrite = [];
end
if TrackAxon
    TitleStr                   = 'axon';
    Axon                       = DefineProcessPolygon(TitleStr,f,AxonWidthWithBackground, SomaCoordinatesChangeRelativeToReferenceFrame, NumOfPixelsForZoom_Axon);  
    AxonMask                   = Axon.ProcessMask;
    Axon.NumOfPixels           = round(length(find(AxonMask(:)))* FractionOfAxonPixels);    
    Axon.MovieInfo.PathName    = PathName;
    Axon.MovieInfo.FrameNames  = FrameNames;
else
    Axon = [];
end

%% Main loop: Run over all frames and extract neural activity
disp(['Extracting neural activity from the following path: ',PathName]);
waitbar_handle = waitbar(0, 'Extracting processes activity'); 
for frame = 1:MovieLength 
    waitbar(frame/MovieLength, waitbar_handle)
    FileName  = [PathName,'\',FrameNames{frame}];
    f         = imread(FileName);    
    
    if TrackDendrite    
        Xmin = Dendrite.InititalXminVec(frame);
        Xmax = Xmin + Dendrite.MatrixXlength-1;
        Ymin = Dendrite.InititalYminVec(frame);
        Ymax = Ymin + Dendrite.MatrixYlength-1;
        CroppedFrame                         = f(Ymin:Ymax, Xmin:Xmax);       
        Dendrite.RawMatrix(:,:,frame)        = CroppedFrame;
        CroppedFrame(~DendriteMask)          = NaN;
        Dendrite.RawMatrix_Masked(:,:,frame) = CroppedFrame;  
    end
    if TrackAxon    
        Xmin = Axon.InititalXminVec(frame);
        Xmax = Xmin + Axon.MatrixXlength-1;
        Ymin = Axon.InititalYminVec(frame);
        Ymax = Ymin + Axon.MatrixYlength-1;
        CroppedFrame                         = f(Ymin:Ymax, Xmin:Xmax);       
        Axon.RawMatrix(:,:,frame)            = CroppedFrame;
        CroppedFrame(~AxonMask)              = NaN;
        Axon.RawMatrix_Masked(:,:,frame)     = CroppedFrame;  
    end    
end
close(waitbar_handle) 

%% calculate process activity by taking the mean over all process pixels. Process pixels are found by their relative intensity 
if TrackDendrite 
    Dendrite = CalculateActivity(Dendrite);
end
if TrackAxon 
    Axon     = CalculateActivity(Axon);
end

%% Save Mat File
MatFileName = [PathName,'\','ProcessActivity.mat'];
save(MatFileName,'Dendrite','Axon','-v7.3');
disp(['Finished saving process activity. File name: ''',MatFileName,'''']) 

%% figure of neuron and processes activity
figure('name','neuron and processes activity','position', [102   125   560   878]); 
subplot(3,1,1); 
    plot(Soma.NeuronMeanRawValue,'k-'); hold on; 
    plot(Soma.NeuronBackground,'k:'); hold on; 
    LegendCell = {'soma activity','soma background'};
    if TrackDendrite 
        plot(Dendrite.ProcessMeanRawValue,'b-'); hold on; 
        plot(Dendrite.ProcessBackground,'b:'); hold on;         
        LegendCell(end+1:end+2) = {'dendrite activity','dendrite background'};
    end
    if TrackAxon 
        plot(Axon.ProcessMeanRawValue,'r-'); hold on; 
        plot(Axon.ProcessBackground,'r:'); hold on;         
        LegendCell(end+1:end+2) = {'axon activity','axon background'};
    end    
    ylabel('raw fluorescence'); 
    legend(LegendCell);
    xlim([0 MovieLength]);
    
subplot(3,1,2); 
    plot(Soma.NeuronValue_BackgroundSubtructed,'k-'); hold on; 
    line([0 MovieLength],[Soma.Neuron_F0, Soma.Neuron_F0],'color','k','linestyle',':');        
    LegendCell = {'soma activity','soma F0'};
    if TrackDendrite 
        plot(Dendrite.ProcessValue_BackgroundSubtructed,'b-'); hold on; 
        line([0 MovieLength],[Dendrite.Process_F0, Dendrite.Process_F0],'color','b','linestyle',':');        
        LegendCell(end+1:end+2) = {'dendrite activity','dendrite F0'};
    end
    if TrackAxon 
        plot(Axon.ProcessValue_BackgroundSubtructed,'r-'); hold on; 
        line([0 MovieLength],[Axon.Process_F0, Axon.Process_F0],'color','r','linestyle',':');        
        LegendCell(end+1:end+2) = {'axon activity','axon F0'};
    end            
    ylabel('raw fluorescence background corrected'); 
    legend(LegendCell);
    xlim([0 MovieLength]);    
subplot(3,1,3); 
    plot(Soma.Neuron_DeltaFOverF,'k-'); hold on; 
    LegendCell = {'soma activity'};
    if TrackDendrite 
        plot(Dendrite.Process_DeltaFOverF,'b-'); hold on; 
        LegendCell(end+1) = {'dendrite activity'};
    end
    if TrackAxon 
        plot(Axon.Process_DeltaFOverF,'r-'); hold on; 
        LegendCell(end+1) = {'axon activity'};
    end            
    ylabel('\DeltaF/F');
    legend(LegendCell);
    xlim([0 MovieLength]);
    
%% Run fast movie check
MinMaxSomaMovie     = prctile(single(Soma.NeuronRawMatrix(:)),[0.1 99.9]);
if TrackDendrite, MinMaxDendriteMovie = prctile(single(Dendrite.RawMatrix(:)),[0.1 99.9]); end
if TrackAxon,     MinMaxAxonMovie     = prctile(single(Axon.RawMatrix(:)),[0.01 99.99]);   end
pause(0.2);
figure; 
pause(0.5);
for frame = 1:10:MovieLength
    subplot(1,3,1);
    CroppedFrame = Soma.NeuronRawMatrix(:,:,frame);
    imshow(CroppedFrame,MinMaxSomaMovie,'initialmagnification',400); title(frame)
    if isfield(Soma,'Polygon')   
        hold on; plot(Soma.Polygon([1:end 1],1),Soma.Polygon([1:end 1],2),'r'); hold off;
    end
    
    if TrackDendrite 
        subplot(1,3,2);
        CroppedFrame = Dendrite.RawMatrix(:,:,frame);
        imshow(CroppedFrame,MinMaxDendriteMovie,'initialmagnification',400); title('dendrite')
        PolygonX = Dendrite.Polygon.PolygonXcoordinates-min(Dendrite.Polygon.PolygonXcoordinates)+1;    
        PolygonY = Dendrite.Polygon.PolygonYcoordinates-min(Dendrite.Polygon.PolygonYcoordinates)+1;    
        hold on; plot(PolygonX([1:end 1]),PolygonY([1:end 1]),'r'); hold off;
    end
    if TrackAxon 
        subplot(1,3,3);
        CroppedFrame = Axon.RawMatrix(:,:,frame);
        imshow(CroppedFrame,MinMaxAxonMovie,'initialmagnification',400); title('axon')
        PolygonX = Axon.Polygon.PolygonXcoordinates-min(Axon.Polygon.PolygonXcoordinates)+1;    
        PolygonY = Axon.Polygon.PolygonYcoordinates-min(Axon.Polygon.PolygonYcoordinates)+1;    
        hold on; plot(PolygonX([1:end 1]),PolygonY([1:end 1]),'r'); hold off;
    end
    
    pause(0.01);
end
keyboard;
% close all hidden;
return

function Process = DefineProcessPolygon(TitleStr,f,ProcessWidthWithBackground, SomaCoordinatesChangeRelativeToReferenceFrame, NumOfPixelsForZoom)
SomaCoordinatesChangeRelativeToReferenceFrame_X = SomaCoordinatesChangeRelativeToReferenceFrame(:,2);
SomaCoordinatesChangeRelativeToReferenceFrame_Y = SomaCoordinatesChangeRelativeToReferenceFrame(:,1);
MovieLength = length(SomaCoordinatesChangeRelativeToReferenceFrame_X);

figure; 
imshow(f,[]); 
title(['Click on middle of ',TitleStr,' for zoom']);   
[I, J] = ginput(1);
Xmin = max([1,         round(J-NumOfPixelsForZoom/2)]);
Xmax = min([size(f,1), round(J+NumOfPixelsForZoom/2)]);
Ymin = max([1,         round(I-NumOfPixelsForZoom/2)]);
Ymax = min([size(f,2), round(I+NumOfPixelsForZoom/2)]);
fzoomed = f(Xmin:Xmax, Ymin:Ymax);

imshow(fzoomed,[],'initialmagnification',400); 
title(['Click on ',TitleStr,' line. Then press Enter']);   
[I_line,J_line]= ginput;
I_line = I_line + Ymin-1;
J_line = J_line + Xmin-1;

[PolygonXcoordinates, PolygonYcoordinates] = DefinePolygonFromLine(I_line,J_line, ProcessWidthWithBackground);

imshow(f,[]); 
hold on;
plot(PolygonXcoordinates([1:end 1]),PolygonYcoordinates([1:end 1]),'b');
plot(I_line,J_line,'r-');

Polygon.Xline                = I_line;
Polygon.Yline                = J_line;
Polygon.PolygonXcoordinates  = PolygonXcoordinates;
Polygon.PolygonYcoordinates  = PolygonYcoordinates;
PolygonLength                = length(PolygonXcoordinates);

PolygonMatrix.PolygonXcoordinates = repmat(PolygonXcoordinates,MovieLength,1) + repmat(SomaCoordinatesChangeRelativeToReferenceFrame_X,1,PolygonLength);
PolygonMatrix.PolygonYcoordinates = repmat(PolygonYcoordinates,MovieLength,1) + repmat(SomaCoordinatesChangeRelativeToReferenceFrame_Y,1,PolygonLength);

Process.Polygon          = Polygon;
Process.PolygonMatrix    = PolygonMatrix;
Process.InititalXminVec  = floor(min(PolygonXcoordinates)) + SomaCoordinatesChangeRelativeToReferenceFrame_X;
Process.InititalYminVec  = floor(min(PolygonYcoordinates)) + SomaCoordinatesChangeRelativeToReferenceFrame_Y;   
Process.MatrixXlength    = ceil(max(PolygonXcoordinates)-min(PolygonXcoordinates));
Process.MatrixYlength    = ceil(max(PolygonYcoordinates)-min(PolygonYcoordinates));
Process.RawMatrix        = zeros(Process.MatrixYlength,Process.MatrixXlength,MovieLength,'uint16');    
Process.RawMatrix_Masked = zeros(Process.MatrixYlength,Process.MatrixXlength,MovieLength,'uint16');       
[X,Y]                    = meshgrid(1:Process.MatrixXlength,1:Process.MatrixYlength);
ProcessMask              = inpolygon(X,Y,PolygonXcoordinates-min(PolygonXcoordinates)+1,PolygonYcoordinates-min(PolygonYcoordinates)+1);
Process.ProcessMask      = ProcessMask;  

return

function  [PolygonXcoordinates, PolygonYcoordinates] = DefinePolygonFromLine(I_line,J_line, Width)
I_line = I_line'; J_line = J_line';
[THETA,~] = cart2pol(diff(I_line),diff(J_line));
THETA(end+1) = THETA(end);
NormTheta_up   = THETA + pi/2; % normal to line
NormTheta_down = THETA - pi/2; % normal to line

I_line_up   = I_line + Width*cos(NormTheta_up);
I_line_down = I_line + Width*cos(NormTheta_down);
J_line_up   = J_line + Width*sin(NormTheta_up);
J_line_down = J_line + Width*sin(NormTheta_down);
PolygonXcoordinates        = [I_line_up, I_line_down(end:-1:1)];
PolygonYcoordinates        = [J_line_up, J_line_down(end:-1:1)];

% THETA/pi*180

return

function  Process = CalculateActivity(Process)
RawMatrix   = Process.RawMatrix_Masked;
NewLength   = size(RawMatrix,1) * size(RawMatrix,2);
MovieLength = size(RawMatrix,3);
RawMatrix   = reshape(RawMatrix,NewLength, MovieLength);

% Calculate raw signal
MatrixForRawActivityCalculation = sort(RawMatrix,1,'descend');
MatrixForRawActivityCalculation = MatrixForRawActivityCalculation(1:Process.NumOfPixels,:);
ProcessMeanRawValue             = single(mean(MatrixForRawActivityCalculation,1));

% Calculate background
MatrixForBackgroundCalculation = sort(RawMatrix,1);
MatrixForBackgroundCalculation = MatrixForBackgroundCalculation(~(sum(MatrixForBackgroundCalculation,2)==0),:);
MatrixForBackgroundCalculation = MatrixForBackgroundCalculation(1:5,:);
ProcessBackground              = single(mean(MatrixForBackgroundCalculation,1));

% Calculate process activity
ProcessValue_BackgroundSubtructed = ProcessMeanRawValue - ProcessBackground;
Process_F0                        = min(ProcessValue_BackgroundSubtructed);
Process_DeltaFOverF               = (ProcessValue_BackgroundSubtructed - Process_F0)/Process_F0; 

% assign to structure 
Process.ProcessMeanRawValue                                 = ProcessMeanRawValue;
Process.ProcessBackground                                   = ProcessBackground;
Process.ProcessValue_BackgroundSubtructed                   = ProcessValue_BackgroundSubtructed;
Process.Process_F0                                          = Process_F0;
Process.Process_DeltaFOverF                                 = Process_DeltaFOverF;
return


















