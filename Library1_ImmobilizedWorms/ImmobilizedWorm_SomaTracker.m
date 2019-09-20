function OutputStructure = ImmobilizedWorm_SomaTracker(PathName, MovieLength, TrackUserDefinedNeuron)

%% This function extract basic soma fluorescence pattern, plots results, and run a short movie to evaluate the segmentation performance 
%
%% Input descriptions: 
%  PathName                 Path of movie folder 
%  MovieLength              Total number of frames relevant for image processing  
%  TrackUserDefinedNeuron   false: Auto detect soma. This assumes the soma fluoresnce of the relevant neuron always contains the brightest pixels within the frame 
%                           true:  Use user interface to define soma boundaries in a single frame and auto-analyze other frames       
%                                  Assumes the soma fluoresnce of the neuron of interest always contains the brightest pixels in the frame 
%
%% Free parameters 
%  Optimized for 100x objectives with demagnification of 0.63x
NeuronAreaInPixels = 200;  % Radious ~ 8 pixels
LengthOfSquare     = 51;
MaxAllowedDistanceFrameNeuronCenter = 25;
AllowJittersInTracking   = true;    % false may be used only if soma location is perfectly stable within the full movie length 

MaxDistanceBetweenFrames = 8; % This parameter is important only if TrackUserDefinedNeuron==true
if ~exist('TrackUserDefinedNeuron','var')
    TrackUserDefinedNeuron = false;
end

%% Get movie information
% Get path if input doesn't exist
if ~exist('PathName','var') || isempty(PathName)
    UserPrompt     = 'Select a parent directory that contains all relevant movies files';
    PathName = uigetdir('', UserPrompt);
end
% Get movie files, one frame per file
FrameNames     = dir(PathName);
FrameNames     = {FrameNames(3:end).name};
% Here movie files have initials: "StreamEXP_" and other files are excluded. Please modify initials to match your movie file names  
C              = strfind(FrameNames,'StreamEXP_');          
TrueMovieFile  = false(1,length(C));    
for file_ind = 1:length(C)
    TrueMovieFile(file_ind) = ~ isempty(C{file_ind}); 
end 
FrameNames = FrameNames(TrueMovieFile);
% Correct movie length if shorter than suggested in input variable  
if exist('MovieLength','var') && ~isempty(MovieLength)
    MovieLength = min([MovieLength length(FrameNames)]);
else
    MovieLength = length(FrameNames);
end

%% Use user interface to define soma boundaries if needed
if TrackUserDefinedNeuron 
    frame = 1;
    FileName  = [PathName,'\',FrameNames{frame}];
    f         = imread(FileName);
    figure; imshow(f,[]); title('Click on Neuron Centroid. Then Click Enter');   
    [I,J]= ginput(1);
    CurrentCoordinatesGuess = [J,I];
    NeuronRawMatrix_Masked  = zeros(LengthOfSquare,LengthOfSquare,MovieLength,'uint16');
end

%% Run over all frames and extract neural activity
% Initialization
NeuronCoordinates  = zeros(MovieLength,2,'single');
NeuronRawMatrix    = zeros(LengthOfSquare,LengthOfSquare,MovieLength,'uint16');

disp(['Extracting neural activity from the following path: ',PathName]);
waitbar_handle = waitbar(0, 'Extracting neural activity'); 
IndicesVector = ceil(-(LengthOfSquare/2)):floor((LengthOfSquare/2));

% Main loop
for frame = 1:MovieLength
    waitbar(frame/MovieLength, waitbar_handle)
    FileName  = [PathName,'\',FrameNames{frame}];
    f         = imread(FileName);
    if frame==1, FrameSize = size(f); end
    BlurryFrame = imfilter(f-min(f(:)),fspecial('average',15));   % Please modify filter window to match your movie resolution and quality  
       
    if TrackUserDefinedNeuron
        if frame ==1
            I = CurrentCoordinatesGuess(1);
            J = CurrentCoordinatesGuess(2);
            Xmin = max([1            round(I-(LengthOfSquare-1)/2)]);
            Xmax = min([FrameSize(1) round(I+(LengthOfSquare-1)/2)]);
            Ymin = max([1            round(J-(LengthOfSquare-1)/2)]);
            Ymax = min([FrameSize(2) round(J+(LengthOfSquare-1)/2)]);
            BlurryFrame(1:Xmin,:)   = 0;
            BlurryFrame(:,1:Ymin)   = 0;
            BlurryFrame(Xmax:end,:) = 0;
            BlurryFrame(:,Ymax:end) = 0;   
        else
            BlurryFrame(~RelevantPartOfFullImage) = 0;
        end
    end
    
    [a, ind] = max(BlurryFrame(:));
    [I,J] = ind2sub(FrameSize,ind);
    NeuronCoordinates(frame,:) = [I, J]';
    if (~AllowJittersInTracking)&&(frame>1)
        DistanceToPreviousCoordinate = sqrt(sum(NeuronCoordinates(frame,:)-NeuronCoordinates(frame-1,:)).^2);
        if DistanceToPreviousCoordinate>MaxDistanceBetweenFrames
            [frame, DistanceToPreviousCoordinate]
            % If segmentation is wrong segmentation (distance is too far) then use previous neuron location
            NeuronCoordinates(frame,:) = NeuronCoordinates(frame-1,:); 
            I                          = NeuronCoordinates(frame-1,1);
            J                          = NeuronCoordinates(frame-1,2);
        end        
    end
    if TrackUserDefinedNeuron
        CurrentCoordinatesGuess(1)=I;
        CurrentCoordinatesGuess(2)=J;
    end
        
    % figure; imshow(f,[]); hold on; plot(J,I,'r*'); figure; imshow(BlurryFrame,[]); hold on; plot(J,I,'r*')
    Xmin = max([1            round(I-(LengthOfSquare-1)/2)]);
    Xmax = min([FrameSize(1) round(I+(LengthOfSquare-1)/2)]);
    Ymin = max([1            round(J-(LengthOfSquare-1)/2)]);
    Ymax = min([FrameSize(2) round(J+(LengthOfSquare-1)/2)]);
    
    CroppedFrame = f(Xmin:Xmax, Ymin:Ymax);    
    NeuronRawMatrix(:,:,frame) = CroppedFrame;
    
    if TrackUserDefinedNeuron
        if frame==1   %% Use user interface to define soma boundaries if needed
            figure; imshow(CroppedFrame,[],'initialmagnification',400); 
            title('choose polygon that encompass the neuron, then press enter')
            [Xv,Yv] = ginput;
            if isempty(Yv)
                Xv = [ 1 LengthOfSquare LengthOfSquare 1]';
                Yv = [ 1 1 LengthOfSquare LengthOfSquare]';                
            end
            [X,Y] = meshgrid(1:LengthOfSquare,1:LengthOfSquare);
            NeuronMask = inpolygon(X,Y,Xv,Yv);            
            CroppedFrame(~NeuronMask) = 0;
            figure; imshow(CroppedFrame,[],'initialmagnification',400);            
        end
        CroppedFrame(~NeuronMask) = 0;
        NeuronRawMatrix_Masked(:,:,frame) = CroppedFrame;                
    end    
    RelevantPartOfFullImage = false(size(f));
    RelevantPartOfFullImage(NeuronCoordinates(frame,1)+IndicesVector, NeuronCoordinates(frame,2)+IndicesVector) = NeuronMask;
%     figure; imshow(f,[]);  figure; imshow(CroppedFrame,[]); 
end
close(waitbar_handle) 

%% Calculate neuron activity by taking the mean over all neuron pixels
ind1 = round((LengthOfSquare-1)/2) - MaxAllowedDistanceFrameNeuronCenter+1;
ind2 = round((LengthOfSquare+1)/2) + MaxAllowedDistanceFrameNeuronCenter;
if TrackUserDefinedNeuron
    MatrixForRawActivityCalculation = NeuronRawMatrix_Masked(ind1:ind2,ind1:ind2,:);
else
    MatrixForRawActivityCalculation = NeuronRawMatrix(ind1:ind2,ind1:ind2,:);
end
NewLength = length(ind1:ind2)^2;
MatrixForRawActivityCalculation = reshape(MatrixForRawActivityCalculation,NewLength, MovieLength);
MatrixForRawActivityCalculation = sort(MatrixForRawActivityCalculation,1,'descend');
MatrixForRawActivityCalculation = MatrixForRawActivityCalculation(1:NeuronAreaInPixels,:);
NeuronMeanRawValue              = single(mean(MatrixForRawActivityCalculation,1));
if ind1>1
    MatrixForBackgroundCalculation  = NeuronRawMatrix([1:ind1-1, ind2+2:end],[1:ind1-1, ind2+2:end],:);
    NewLength = length([1:ind1-1, ind2+2:LengthOfSquare])^2;
else
    MatrixForBackgroundCalculation  = NeuronRawMatrix;
    NewLength = LengthOfSquare^2;
end
    
MatrixForBackgroundCalculation = reshape(MatrixForBackgroundCalculation,NewLength, MovieLength);
NeuronBackground                 = single(min(MatrixForBackgroundCalculation,[],1)); 
NeuronValue_BackgroundSubtructed = NeuronMeanRawValue - NeuronBackground;

Neuron_F0          = min(NeuronValue_BackgroundSubtructed);
Neuron_DeltaFOverF = (NeuronValue_BackgroundSubtructed - Neuron_F0)/Neuron_F0; 

OutputStructure.FreeParameters.NeuronAreaInPixels                  = NeuronAreaInPixels;
OutputStructure.FreeParameters.LengthOfSquare                      = LengthOfSquare;
OutputStructure.FreeParameters.MaxAllowedDistanceFrameNeuronCenter = MaxAllowedDistanceFrameNeuronCenter;
OutputStructure.MovieInfo.PathName                                 = PathName;
OutputStructure.MovieInfo.FrameNames                               = FrameNames;
OutputStructure.NeuronCoordinates                                  = NeuronCoordinates;
OutputStructure.NeuronRawMatrix                                    = NeuronRawMatrix;

if TrackUserDefinedNeuron
    OutputStructure.TrackUserDefinedNeuron                         = true;
    OutputStructure.NeuronMaskDefinedByUser                        = NeuronMask;
    OutputStructure.Polygon                                        = [Xv Yv];
    OutputStructure.NeuronRawMatrix_Masked                         = NeuronRawMatrix_Masked;
end

OutputStructure.NeuronMeanRawValue                                 = NeuronMeanRawValue;
OutputStructure.NeuronBackground                                   = NeuronBackground;
OutputStructure.NeuronValue_BackgroundSubtructed                   = NeuronValue_BackgroundSubtructed;
OutputStructure.Neuron_F0                                          = Neuron_F0;
OutputStructure.Neuron_DeltaFOverF                                 = Neuron_DeltaFOverF;

% Save Mat File
MatFileName = [PathName,'\','NeuralActivity.mat'];
save(MatFileName,'OutputStructure','-v7.3');
disp(['Finished saving neural activity. File name: ''',MatFileName,'''']) 

%% figures
figure('name','neuron coordinates'); 
subplot(1,2,1); 
    plot(OutputStructure.NeuronCoordinates(:,1),OutputStructure.NeuronCoordinates(:,2),'k-'); hold on; 
    plot(OutputStructure.NeuronCoordinates(1,1),OutputStructure.NeuronCoordinates(1,2),'r*'); hold on; 
    plot(OutputStructure.NeuronCoordinates(end,1),OutputStructure.NeuronCoordinates(end,2),'bo'); hold on; 
    xlabel('X');ylabel('Y'); axis equal
subplot(2,2,2); 
    plot(OutputStructure.NeuronCoordinates(:,1),'k-'); hold on; 
    ylabel('X');
subplot(2,2,4); 
    plot(OutputStructure.NeuronCoordinates(:,2),'k-'); hold on; 
    ylabel('Y');
figure('name','neuron activity'); 
subplot(3,1,1); 
    plot(OutputStructure.NeuronMeanRawValue,'r-'); hold on; 
    plot(OutputStructure.NeuronBackground,'k-'); hold on; 
    ylabel('raw fluorescence'); legend('neuron activity','background');
    xlim([0 MovieLength]);
subplot(3,1,2); 
    plot(OutputStructure.NeuronValue_BackgroundSubtructed,'r-'); hold on; 
    line([0 MovieLength],[OutputStructure.Neuron_F0, OutputStructure.Neuron_F0],'color','r','linestyle',':');
    ylabel('raw fluorescence background corrected'); legend('neuron activity','F0');
    xlim([0 MovieLength]);
subplot(3,1,3); 
    plot(OutputStructure.Neuron_DeltaFOverF,'k-'); hold on; 
    ylabel('\DeltaF/F');
    xlim([0 MovieLength]);
    
%% Run fast movie check
MinMaxMovie = prctile(single(OutputStructure.NeuronRawMatrix(:)),[0.1 99.9]);
figure; 
for frame = 1:10:MovieLength
    CroppedFrame = OutputStructure.NeuronRawMatrix(:,:,frame);
    imshow(CroppedFrame,MinMaxMovie,'initialmagnification',400); title(frame)
    if TrackUserDefinedNeuron   
        hold on; plot(OutputStructure.Polygon([1:end 1],1),OutputStructure.Polygon([1:end 1],2),'r'); hold off;
    end
    pause(0.01);
end
% close all;
return


