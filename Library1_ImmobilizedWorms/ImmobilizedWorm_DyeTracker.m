function OutputStructure = ImmobilizedWorm_DyeTracker(PathName, DyeOnFrame, DyeOffFrame)
% PathName    >> Path of movie files
% DyeOnFrame  >> Frame when dye flow starts 
% DyeOffFrame >> Frame when dye flow ends   

%% Free Parameters Choice for 100x with demagnification of 0.63
LengthOfSquare     = 61;

%% Get movie information
% Get path if input doesn't exist
if ~exist('PathName','var') || isempty(PathName)
    UserPrompt     = 'Select a parent directory that contains all relevant movies files';
    PathName = uigetdir('', UserPrompt);
end
FrameNames     = dir(PathName);
FrameNames     = {FrameNames(3:end).name};
C              = strfind(FrameNames,'StreamEXP_');
TrueMovieFile  = false(1,length(C));    
for file_ind = 1:length(C)
    TrueMovieFile(file_ind) = ~ isempty(C{file_ind}); 
end 
FrameNames = FrameNames(TrueMovieFile);
MovieLength = length(FrameNames);


%% User interface for dye location
frame     = DyeOffFrame;
FileName  = [PathName,'\',FrameNames{frame}];
f         = imread(FileName);
figure; imshow(f,[]); title('Click on the middle of dye flow area'); 
FrameSize = size(f);
[J,I]= ginput(1);
Xmin = max([1            round(I-(LengthOfSquare-1)/2)]);
Xmax = min([FrameSize(1) round(I+(LengthOfSquare-1)/2)]);
Ymin = max([1            round(J-(LengthOfSquare-1)/2)]);
Ymax = min([FrameSize(2) round(J+(LengthOfSquare-1)/2)]);
hold on; 
rectangle('Position',[Ymin,Xmin,Ymax-Ymin,Xmax-Xmin],'edgecolor','r');

%% Run over all frames and extract DyePattern
disp(['Extracting dye pattern from the following path: ',PathName]);
waitbar_handle = waitbar(0, 'Extracting dye patten'); 
DyePattern    = zeros(1,MovieLength,'single');  % Initialization
for frame = 1:MovieLength
    waitbar(frame/MovieLength, waitbar_handle)
    FileName  = [PathName,'\',FrameNames{frame}];
    f         = imread(FileName);
%     figure; imshow(f,[]);        
    CroppedFrame      = f(Xmin:Xmax, Ymin:Ymax);
    DyePattern(frame) = mean(CroppedFrame(:));        
end
close(waitbar_handle) 
    
FirstFrameForBackground = max([1, DyeOnFrame-200]);
MeanBackground = mean(DyePattern(FirstFrameForBackground:DyeOnFrame));
StdBackground  = std(DyePattern(FirstFrameForBackground:DyeOnFrame),1);
DyePattern_BackgroundSubtracted = DyePattern - MeanBackground;
DyePattern_Normalized           = DyePattern_BackgroundSubtracted ./ max(DyePattern_BackgroundSubtracted);
% StdBackground_Normalized      = StdBackground / max(DyePattern_BackgroundSubtracted);
StdBackground_Normalized        = 2e-3;  % Free Parameter - It is robust, just needs initial adjustment to match the imaging quality

%% Find Delay time between switch and beginning of dye change and extract dye pattern ON and OFF vectors
DyeIncreasing      = DyePattern_Normalized(DyeOnFrame:DyeOffFrame);
DelayTimeInFrames  = find(DyeIncreasing<StdBackground_Normalized,1,'last');
DelayTimeInSeconds = DelayTimeInFrames/10;

DyePatternParams.MeanBackground            = MeanBackground;
DyePatternParams.StdBackground             = StdBackground;
DyePatternParams.StdBackground_Normalized  = StdBackground_Normalized;
DyePatternParams.DelayTimeInFrames         = DelayTimeInFrames;
DyePatternParams.DelayTimeInSeconds        = DelayTimeInSeconds;

DyeOnset                     = DyePattern_Normalized(DyeOnFrame:(DyeOffFrame+DelayTimeInFrames));
DyeOnset_TimingWithoutDelay  = -DelayTimeInFrames:(DyeOffFrame-DyeOnFrame);
DyeOffset                    = DyePattern_Normalized(DyeOffFrame:end);
DyeOffset_TimingWithoutDelay = -DelayTimeInFrames:(length(DyePattern_Normalized)-DyeOffFrame-DelayTimeInFrames);

% ComputeRiseTimes 
RiseTimePercentages = [10 25 50 75 90 95];
VecOn  = DyeOnset(DelayTimeInFrames:end);
VecOff = DyeOffset(DelayTimeInFrames:end);

DyePatternParams.RiseTimes_10_25_50_75_90_95  = zeros(1,length(RiseTimePercentages),'single')*NaN;
DyePatternParams.DecayTimes_10_25_50_75_90_95 = zeros(1,length(RiseTimePercentages),'single')*NaN;
for percentage_ind = 1:length(RiseTimePercentages)
    Percentage = RiseTimePercentages(percentage_ind);
    OnFrame  = find(VecOn>Percentage/100,1,'first');
    if ~isempty(OnFrame)
        DyePatternParams.RiseTimes_10_25_50_75_90_95(percentage_ind) = OnFrame;
    end
    OffFrame = find(VecOff<(1-Percentage/100),1,'first');
    if ~isempty(OffFrame)
        DyePatternParams.DecayTimes_10_25_50_75_90_95(percentage_ind)= OffFrame;
    end
end    

%% Assign to structure
OutputStructure.FreeParameters.DyeSquareCentroid    = [J,I];
OutputStructure.FreeParameters.LengthOfSquare       = LengthOfSquare;
OutputStructure.MovieInfo.PathName                  = PathName;
OutputStructure.MovieInfo.FrameNames                = FrameNames;
OutputStructure.DyePattern                          = DyePattern;
OutputStructure.DyePattern_BackgroundSubtracted     = DyePattern_BackgroundSubtracted;
OutputStructure.DyePattern_Normalized               = DyePattern_Normalized;
OutputStructure.DyeOnset                            = DyeOnset;
OutputStructure.DyeOnset_TimingWithoutDelay         = DyeOnset_TimingWithoutDelay;
OutputStructure.DyeOffset                           = DyeOffset;
OutputStructure.DyeOffset_TimingWithoutDelay        = DyeOffset_TimingWithoutDelay;
OutputStructure.DyePatternParams                    = DyePatternParams;

%% figure
figure('name','Dye pattern'); 
subplot(1,2,1); 
    plot(OutputStructure.DyePattern,'k-'); hold on; 
    ylabel('raw dye pattern'); 
    xlim([0 MovieLength]);
subplot(1,2,2); 
    plot(OutputStructure.DyeOnset_TimingWithoutDelay,OutputStructure.DyeOnset,'r-','linewidth',2); hold on; 
    plot(OutputStructure.DyeOffset_TimingWithoutDelay,OutputStructure.DyeOffset,'b-','linewidth',2); hold on; 
    line([0 0],[0 1],'color','k');
    ylabel('Normalized Dye Patterns'); legend('Onset','Offset'); 
    title(['Delay Time = ',num2str(OutputStructure.DyePatternParams.DelayTimeInFrames),' frames']);
    ylim([-0.01 1]); xlim([-OutputStructure.DyePatternParams.DelayTimeInFrames, max([OutputStructure.DyeOnset_TimingWithoutDelay, OutputStructure.DyeOffset_TimingWithoutDelay])]);

%% Save Mat File
MatFileName = [PathName,'\','DyePattern.mat'];
tic
save(MatFileName,'OutputStructure','-v7.3');
toc
disp(['Finished saving dye pattern. File name: ''',MatFileName,'''']) 

return