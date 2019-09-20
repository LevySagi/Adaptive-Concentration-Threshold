function [File, Data] = CalculateActivityAndStimulusFeatures_04(File, Data, ExpInfo, plotme)

% September 2019, Sagi Levy

% % Example for ExpInfo
% ExpInfo.ConcentrationAtArenaTop    = 1e-6;   % dilution
% ExpInfo.ConcentrationAtArenaBottom = 1e-8;   % dilution

%% Initialization and free parameters
MaxAllowedMissingFrames                       = 0.3*File.FrameRate;   % max frames for considering it a high confidence for activity interpolation 
MinAverageForDownGradientMovement_CTXperFrame = 0.01/File.FrameRate;  % max(CTX/frame) for considering it 'downgradient' movement, corresponds to 21.5 micrometer/sec (4300 micrometer in 2CTX)
% Activity and Peak parameters. Note that the following parameters may depend on imaging quality and typical signals of the neuron of interest 
PRCTILE_For_F0  = 20;     % 20% Only for display of full tracks. Not used in further analysis
ThresholdValue  = 0.5;
ProminenceValue = 0.1;    % minimum for the non-averaged signal
MaxTimeBetweenPeakAndAverageSlopeInitiation   = 15; % seconds
MinOfMaxSlopeVal = 0.01;
MinActivityGradientForSharpActivation         = 1/1.5;  % minimal activation change (of filtered activity) rate to be considered "Sharp". Units: (DeltaF/F)/second
MinForSharpActivation_DeltaFOverFperFrame     = MinActivityGradientForSharpActivation/ File.FrameRate;
MaxFramesBetweenPeakAndAverageSlopeInitiation = MaxTimeBetweenPeakAndAverageSlopeInitiation* File.FrameRate;
% Filter parameters
LPS = designfilt('lowpassfir', ...
  'PassbandFrequency',1,'StopbandFrequency',5, ...
  'PassbandRipple',1,'StopbandAttenuation',80, ...
  'DesignMethod','equiripple','SampleRate',30);
D_LPS = round(mean(grpdelay(LPS)));
% fvtool(LPS)  
% average windows parameters
windowSize      = File.FrameRate/2;
D_AverageWindow = round(windowSize/2);
b               = (1/windowSize)*ones(1,windowSize);
a               = 1;
windowSize_Activity = File.FrameRate*2;
D_AverageActivity   = round(windowSize_Activity/2);
b_Activity          = (1/windowSize_Activity)*ones(1,windowSize_Activity);

if ~exist('plotme','var')
    plotme = false;
end
NumOfTracks = size(Data.Coordinates_X,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Extraction of filtered stimulus and neuronal activity experienced by the animal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  Stimulus is measured in Chemotaxis Index units based on location of the head in the gradient(-1<CTX<1)
%  Neuronal activity is measured by NeuronValue (fluoresence- background fluoresence) or DeltaF/F     
%% Interpolate Stimulus and Neuron signals to estimate values at NaN frames   
%  High confidence = interpolate only fragments with up to several missing frames  
%  Low confidence  = interpolate all NaN frames  
fieldsIn = {'Stimulus_Head','Stimulus_Neuron','Stimulus_Centroid','FlowAxis_Head','FlowAxis_Neuron','FlowAxis_Centroid','NeuronValue'};
for f_ind = 1:length(fieldsIn)
    fieldnameIn = fieldsIn{f_ind};
    NewName_Low        = [fieldnameIn,'_InterpLowCondifence'];
    NewName_High       = [fieldnameIn,'_InterpHighCondifence'];
    Data.(NewName_Low) = zeros(size(Data.(fieldnameIn)),'single');
    Data.(NewName_High)= zeros(size(Data.(fieldnameIn)),'single');
    for tr =1:NumOfTracks
        Stimulus = Data.(fieldnameIn)(tr,:);
        % Low confidence interpolation: use data from non-NaN frames to calculate the expected data in all frames  
        NaN_Frames                   = isnan(Stimulus)';
        X                            = find(~NaN_Frames);
        InterpStimulus_LowConfidence = interp1(X, Stimulus(X), 1:length(Stimulus));
        if NaN_Frames(1)
            % Starts with NaN >> assign value of the first non-NaN frame
           FirstNonNaNFrame = find(NaN_Frames==0,1,'first');
           InterpStimulus_LowConfidence(1:(FirstNonNaNFrame-1))= Stimulus(FirstNonNaNFrame);            
        end
        if NaN_Frames(end)
            % Ends with NaN >> assign value of the last non-NaN frame
           LastNonNaNFrame = find(NaN_Frames==0,1,'last');
           InterpStimulus_LowConfidence((LastNonNaNFrame+1):end)= Stimulus(LastNonNaNFrame);            
        end
        % High confidence interpolation: reject interpolation estimation at times when data is too scarce.    
        FramesWithLongNaNSegments     = FindSegment(NaN_Frames, MaxAllowedMissingFrames); % Find Long NaN segments where the worm was not properly detected     
        InterpStimulus_HighConfidence = InterpStimulus_LowConfidence;
        InterpStimulus_HighConfidence(FramesWithLongNaNSegments) = NaN; 
        % assign to structure        
        Data.(NewName_Low)(tr,:)  = InterpStimulus_LowConfidence;
        Data.(NewName_High)(tr,:) = InterpStimulus_HighConfidence;        
    end
end

%% Stimulus signal vector
%  Stimulus = Head stimulus, except frames where Head is low confidence while neuron is high confidence --> use neuron stimulus at these specific frames. 
Data.Stimulus_LowConfidence = Data.Stimulus_Head_InterpLowCondifence;
Data.Stimulus_HighConfidence= Data.Stimulus_Head_InterpHighCondifence;
Data.FlowAxis_LowConfidence = Data.FlowAxis_Head_InterpLowCondifence;
Data.FlowAxis_HighConfidence= Data.FlowAxis_Head_InterpHighCondifence;
for tr =1:NumOfTracks
    NeuronHigh_HeadLow = find(isnan(Data.Stimulus_Head_InterpHighCondifence(tr,:)) & ...
                             ~isnan(Data.Stimulus_Neuron_InterpHighCondifence(tr,:))); 
    Data.Stimulus_HighConfidence(tr,NeuronHigh_HeadLow) = Data.Stimulus_Neuron_InterpHighCondifence(tr,NeuronHigh_HeadLow);
    Data.Stimulus_LowConfidence(tr,NeuronHigh_HeadLow)  = Data.Stimulus_Neuron_InterpHighCondifence(tr,NeuronHigh_HeadLow);
    
    Data.FlowAxis_HighConfidence(tr,NeuronHigh_HeadLow) = Data.FlowAxis_Neuron_InterpHighCondifence(tr,NeuronHigh_HeadLow);
    Data.FlowAxis_LowConfidence(tr,NeuronHigh_HeadLow)  = Data.FlowAxis_Neuron_InterpHighCondifence(tr,NeuronHigh_HeadLow);
    if plotme
        figure('name',['Stimulus, activity, and non-OOB behavior, Track ', num2str(tr)]); 
        plot(Data.Stimulus_LowConfidence(tr,:),'k.-'); hold on;
        plot(Data.Stimulus_HighConfidence(tr,:),'b.-'); 
        plot(Data.NeuronValue_InterpHighCondifence(tr,:)/max(Data.NeuronValue_InterpHighCondifence(tr,:))*2-1,'k.-'); hold on;
        NotOOB = find(Data.BehaviorCode_LowLevel(tr,:)>0);
        plot(NotOOB,zeros(1,length(NotOOB)),'g.'); hold on;
        legend('Stimulus (low conf.)','Stimulus (high conf.)','Activity (high conf.)','non-OOB behavior');
    end
end

%%  DeltaF/F
MAT     = Data.NeuronValue_InterpHighCondifence;
F0      = prctile(MAT,PRCTILE_For_F0,2);
Data.F0 = F0;
F0_Mat  = repmat(F0,1,size(MAT,2));
Data.DeltaFOverF = (MAT-F0_Mat)./F0_Mat;
if plotme
    figure('name','DeltaF/F dynamics'); 
    for tr=1:NumOfTracks
        subplot(NumOfTracks,1,tr); plot(Data.DeltaFOverF(tr,:)); title(tr)
    end
end
if plotme
    figure('name','DeltaF/F histograms'); 
    for tr=1:NumOfTracks
        subplot(NumOfTracks,1,tr); 
        Edges = -1:0.02:7; X = (Edges(1:end-1)+Edges(2:end))/2;
        N =histc(Data.DeltaFOverF(tr,:),Edges); N=N(1:end-1); N = N/sum(N);
        bar(X,N); title(['Track ',num2str(tr),', skewness = ', num2str(skewness(N)),', %(/DeltaF.F>1) = ', num2str(100*sum(N(X>=1)))]); ylim([0 0.05])
    end
end

%% Filter and analyze activity pattern 
%%% Initialization
InitializationMat         = zeros(size(Data.Stimulus_LowConfidence),'single');
Data.NeuronValueFiltered  = InitializationMat;
Data.DeltaFOverFFiltered  = InitializationMat;
Data.ActivationSegments.InitiationFrames = cell(1,2);
Data.ActivationSegments.PeakFrames       = cell(1,2);
Data.ActivationSegments.OOB              = cell(1,2);
Data.ActivationSegments.SharpActivation  = cell(1,2);

for tr = 1:NumOfTracks
    NeuronValue = Data.NeuronValue_InterpLowCondifence(tr,:);
    NeuronNaNs  = isnan(Data.NeuronValue_InterpHighCondifence(tr,:));
    
    % Filter activity
    y = filter(LPS,[NeuronValue'; zeros(D_LPS,1)]);y = y(D_LPS+1:end);
    NeuronValueFiltered     = y / nanmedian(y)*nanmedian(NeuronValue);
    F0                      = prctile(NeuronValueFiltered,PRCTILE_For_F0);
    Data.F0_Filtered(tr)    = F0;
    DeltaFOverFFiltered     = (NeuronValueFiltered-F0)/F0;
    DiffActivityFiltered    = [NaN diff(DeltaFOverFFiltered')]; % diff filter
    
    % average the filter
    y = filter(b_Activity,a,DeltaFOverFFiltered);     
    y = y(D_AverageActivity+1:end); y(end+1:length(DeltaFOverFFiltered))=NaN;
    AverageActivity = y;
    DiffAverageActivity = [NaN diff(AverageActivity)'];   % diff average filter   

    DeltaFOverFFiltered(NeuronNaNs) = NaN;
    AverageActivity(NeuronNaNs)     = NaN;
    DiffActivityFiltered(NeuronNaNs) = NaN;
    DiffAverageActivity(NeuronNaNs) = NaN;    

    %% findpeaks
    [StartFrames, PeakFrames] = FindActivityPeaks(DeltaFOverFFiltered, AverageActivity, DiffActivityFiltered, DiffAverageActivity, ...
         ThresholdValue, ProminenceValue, MaxFramesBetweenPeakAndAverageSlopeInitiation, MinOfMaxSlopeVal, plotme);
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Assign to structure
    Data.NeuronValueFiltered(tr,:)               = NeuronValueFiltered;
    Data.DeltaFOverFFiltered(tr,:)               = DeltaFOverFFiltered;
    Data.ActivationSegments.InitiationFrames{tr} = StartFrames;
    Data.ActivationSegments.PeakFrames{tr}       = PeakFrames;
end

if plotme
    figure('name','Filtered(\DeltaF/F)'); 
    for tr=1:NumOfTracks
        subplot(NumOfTracks,1,tr)
        plot(Data.DeltaFOverF(tr,:),'k.');         hold on; 
        plot(Data.DeltaFOverFFiltered(tr,:),'ro');
        legend('\DeltaF/F','Filtered(\DeltaF/F)'); title(['Track ',num2str(tr)]);
    end
end

%% Correction after manual inspection of activity events results
load('D:\ImagingInGradientDevices\ManualInspection_Activity.mat','num','txt','raw');
MovieColumns  = find(strcmpi(txt(1,:),File.MovieName));
TracksColumns = num(1,MovieColumns);
FirstRowInNum = find(~isnan(num(:,1)),1,'first');
for tr = 1:NumOfTracks
   DeltaFOverFFiltered  = Data.DeltaFOverFFiltered(tr,:);
   CurrentThreeColumns  = MovieColumns(TracksColumns==tr);
   InitiationFrames     = num(FirstRowInNum:end,CurrentThreeColumns(1));
   PeakFrames           = num(FirstRowInNum:end,CurrentThreeColumns(2));
   OOB                  = num(FirstRowInNum:end,CurrentThreeColumns(3));
   NoNaNs               = ~isnan(InitiationFrames);
   InitiationFrames     = InitiationFrames(NoNaNs);
   PeakFrames           = PeakFrames(NoNaNs);
   OOB                  = logical(OOB(NoNaNs));
   InitiationToPeak     = PeakFrames - InitiationFrames;
   DELTA_DeltaFOverF    = DeltaFOverFFiltered(PeakFrames) - DeltaFOverFFiltered(InitiationFrames);
   SharpActivation      = (DELTA_DeltaFOverF./InitiationToPeak') > MinForSharpActivation_DeltaFOverFperFrame;
   
   Data.ActivationSegments.InitiationFrames{tr} = InitiationFrames;
   Data.ActivationSegments.PeakFrames{tr}       = PeakFrames;
   Data.ActivationSegments.OOB{tr}              = OOB;
   Data.ActivationSegments.SharpActivation{tr}  = SharpActivation;
  
   if plotme
        figure('name',['Track ',num2str(tr),' Activation initiation and peaks']); 
        ActivationInitiationFrames = Data.ActivationSegments.InitiationFrames{tr};
        ActivationPeakFrames       = Data.ActivationSegments.PeakFrames{tr};
        ActivationOOB              = Data.ActivationSegments.OOB{tr};
        Stimulus                   = Data.Stimulus_LowConfidence(tr,:);
        OOB                        = Data.BehaviorCode_LowLevel(tr,:)==0 ;
        Stimulus_OOB               = Stimulus; Stimulus_OOB(~OOB)=NaN;
        subplot(2,1,1)
        plot(DeltaFOverFFiltered,'k.-'); hold on; 
        plot(ActivationInitiationFrames(ActivationOOB), DeltaFOverFFiltered(ActivationInitiationFrames(ActivationOOB)),'ro'); hold on; 
        plot(ActivationPeakFrames(ActivationOOB), DeltaFOverFFiltered(ActivationPeakFrames(ActivationOOB)),'go'); hold on; 
        plot(ActivationInitiationFrames(~ActivationOOB), DeltaFOverFFiltered(ActivationInitiationFrames(~ActivationOOB)),'ro','markerfacecolor','r'); hold on; 
        plot(ActivationPeakFrames(~ActivationOOB), DeltaFOverFFiltered(ActivationPeakFrames(~ActivationOOB)),'go','markerfacecolor','g'); hold on; 
        ylabel('\DeltaF/F'); xlabel('Frame');
        subplot(2,1,2)
        plot(Stimulus,'k.-'); hold on; 
        plot(Stimulus_OOB,'.-','color',[0.5 0.5 0.5]); hold on; 
        % plot(Data.Stimulus_HighConfidence(tr,:),'k.-'); hold on; 
        plot(ActivationInitiationFrames(ActivationOOB), Stimulus(ActivationInitiationFrames(ActivationOOB)),'ro'); hold on; 
        plot(ActivationPeakFrames(ActivationOOB), Stimulus(ActivationPeakFrames(ActivationOOB)),'go'); hold on; 
        plot(ActivationInitiationFrames(~ActivationOOB), Stimulus(ActivationInitiationFrames(~ActivationOOB)),'ro','markerfacecolor','r'); hold on; 
        plot(ActivationPeakFrames(~ActivationOOB), Stimulus(ActivationPeakFrames(~ActivationOOB)),'go','markerfacecolor','g'); hold on; 
        ylabel('Stimulus'); xlabel('Frame');
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Filter locomotion data (e.g. Stimulus), find downgradient navigation segments, find activity and locomotion features during each segment  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%% Initialization
InitializationMat         = zeros(size(Data.Stimulus_LowConfidence),'single');
Data.FlowAxisFiltered     = InitializationMat;
Data.StimulusFiltered     = InitializationMat;
Data.DiffStimulusFiltered = InitializationMat;
Data.AverageStimulus      = InitializationMat;
Data.DiffAverageStimulus  = InitializationMat;
Data.Downgradient         = false(size(Data.Stimulus_LowConfidence));
Data.DowngradientSegments.StartFrames = cell(1,NumOfTracks);
Data.DowngradientSegments.EndFrames   = cell(1,NumOfTracks);

for tr=1:NumOfTracks  % Loop over Tracks
    NaNs         = Data.BehaviorCode_LowLevel(tr,:)==0;
    Stimulus     = Data.Stimulus_LowConfidence(tr,:);
    FlowAxis     = Data.FlowAxis_LowConfidence(tr,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%% Step 1: Define two smoothed coordinates vectors 'StimulusFiltered' and 'AverageStimulus'    
    %           Find two gradient vectors by differentiating them 
    % Filter FlowAxis
    y = filter(LPS,[FlowAxis'; zeros(D_LPS,1)]);y = y(D_LPS+1:end);
    FlowAxisFiltered     = y / prctile(abs(y),99.9);
    % Filter stimulus
    y = filter(LPS,[Stimulus'; zeros(D_LPS,1)]);y = y(D_LPS+1:end);
    StimulusFiltered     = y / prctile(abs(y),99.9);
    DiffStimulusFiltered = [NaN diff(StimulusFiltered')]; % diff filter
    % average the filter
    y = filter(b,a,StimulusFiltered);     
    y = y(D_AverageWindow+1:end); y(end+1:length(Stimulus))=NaN;
    AverageStimulus = y;
    DiffAverageStimulus = [NaN diff(AverageStimulus')];   % diff average filter   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%% Step 2: Define 'Downgradient' by finding high enough AVERAGE downgradient speed:   
    %               DiffAverageStimulus < - MinAverageForDownGradientMovement_CTXperFrame  
    %           Let's take an average of 50 microM/sec
    %           arena: 14*300+100=4300 micrometer = 2 CTX units --> 21.5 microM = 0.01 CTX units 
    %           MinAverageForDownGradientMovement_CTXperFrame = 0.01/File.FrameRate;   

    Downgradient = (DiffAverageStimulus<-MinAverageForDownGradientMovement_CTXperFrame)& ~NaNs ;
    [~, DowngradientStartFrames, DowngradientEndFrames]= FindSegment(Downgradient, 1);           
    if plotme
        ActivationInitiationFrames = Data.ActivationSegments.InitiationFrames{tr};
        ActivationPeakFrames       = Data.ActivationSegments.PeakFrames{tr};
        ActivationOOB              = Data.ActivationSegments.OOB{tr};
        Stimulus_OOB = Stimulus;          Stimulus_OOB(~NaNs)=NaN;
        StimulusDowngradient = Stimulus;  StimulusDowngradient(~Downgradient)= NaN;
        figure('name','filtering and averaging stimulus vector'); 
            plot(Stimulus,'k.'); hold on; plot(StimulusFiltered,'b.'); hold on; plot(AverageStimulus,'r.');
            legend('Stimulus','filtered','smooth filter');
        figure('name','Gradient vectors'); 
            plot(DiffStimulusFiltered,'b.-'); hold on; plot(DiffAverageStimulus,'k.'); hold on; 
            plot(find(Downgradient),DiffAverageStimulus(Downgradient),'r.'); 
            legend('diff(filtered)','diff(smooth filtered)','diff(smooth filtered), Estimated DownGradient');
        figure('name','downgradient estimation and activity estimations'); 
            plot(Stimulus,'k.-'); hold on; 
            plot(StimulusDowngradient,'b.-'); 
            plot(DowngradientStartFrames,Stimulus(DowngradientStartFrames),'bx'); 
            plot(DowngradientEndFrames,Stimulus(DowngradientEndFrames),'b*');         
            plot(Stimulus_OOB,'.-','color',[0.5 0.5 0.5]); 
            plot(ActivationInitiationFrames(ActivationOOB), Stimulus(ActivationInitiationFrames(ActivationOOB)),'ro'); hold on; 
            plot(ActivationPeakFrames(ActivationOOB), Stimulus(ActivationPeakFrames(ActivationOOB)),'go'); hold on; 
            plot(ActivationInitiationFrames(~ActivationOOB), Stimulus(ActivationInitiationFrames(~ActivationOOB)),'ro','markerfacecolor','r'); hold on; 
            plot(ActivationPeakFrames(~ActivationOOB), Stimulus(ActivationPeakFrames(~ActivationOOB)),'go','markerfacecolor','g'); hold on; 
            legend('Stimulus','Estimated down-gradient','start','end','OOB',...
                   'ActivationInitiation-OOB','ActivationPeak-OOB','ActivationInitiation','ActivationPeak')
            %hold on; plot(AverageStimulus,'r.');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Assign to structure
    Data.FlowAxisFiltered(tr,:)     = FlowAxisFiltered;
    Data.StimulusFiltered(tr,:)     = StimulusFiltered;
    Data.DiffStimulusFiltered(tr,:) = DiffStimulusFiltered;
    Data.AverageStimulus(tr,:)      = AverageStimulus;
    Data.DiffAverageStimulus(tr,:)  = DiffAverageStimulus;
    Data.Downgradient(tr,:)         = Downgradient;
    Data.DowngradientSegments.StartFrames{tr} = DowngradientStartFrames;
    Data.DowngradientSegments.EndFrames{tr}   = DowngradientEndFrames;
end

%% Activity pattern during downgradient navigation
Data = FindActivityDuringDowngradient(Data, File.FrameRate, ThresholdValue);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Concentration And Adaptive-Threshold
Data = ExtractConcentrationAndAdaptiveThreshold(Data, File.FrameRate, ExpInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% SAVE
FileName     = [File.MoviePath,'\',File.MovieName,'_DataMatrices_WithGradient.mat'];
File.ExpInfo = ExpInfo;
save(FileName,'File','Data','-v7.3');
disp(['File saved: ', FileName]);

return

function  [VecTrueIfFrameWithinSegment, StartFrames, EndFrames]= FindSegment(LogicalVec, MinimumSegmentLength)
% Example:
% [FramesWithLongNaNSegments, StartNaNSegments, EndNaNSegments]= FindSegment(NaN_Frames, MaxAllowedMissingFrames);

StartFrames = (find(diff(LogicalVec)==1)+1)'; if LogicalVec(1)==1,   StartFrames = [1 StartFrames]; end
EndFrames   = (find(diff(LogicalVec)==-1))';  if LogicalVec(end)==1, EndFrames   = [EndFrames length(LogicalVec)]; end
if length(StartFrames) ~= length(EndFrames)
    disp('error in segments calculation. keyboard...');
    keyboard;
end
LongSegments_indices = find((EndFrames - StartFrames)>MinimumSegmentLength);
if ~isempty(LongSegments_indices)
    VecTrueIfFrameWithinSegment = false(1,length(LogicalVec));
    for ind = LongSegments_indices
        VecTrueIfFrameWithinSegment(StartFrames(ind):EndFrames(ind)) = true;    
    end
end        


return

function  [StartFrames, PeakFrames] = FindActivityPeaks(DeltaFOverFFiltered, AverageActivity, DiffActivityFiltered, DiffAverageActivity, ...
         ThresholdValue, ProminenceValue, MaxFramesBetweenPeakAndAverageSlopeInitiation, MinOfMaxSlopeVal, plotme)
% Find candidate peaks
[pks1,locs1] = findpeaks(double(AverageActivity),'MINPEAKHEIGHT',ThresholdValue,'MinPeakProminence',0);
[pks2,locs2] = findpeaks(double(DeltaFOverFFiltered),'MINPEAKHEIGHT',ThresholdValue,'MinPeakProminence',ProminenceValue);
    
%     figure; plot(DeltaFOverFFiltered,'k.'); hold on; plot(AverageActivity,'g*'); plot(locs1,pks1,'ro'); plot(locs2,pks2,'ko');
%     figure; plot(DiffAverageActivity,'g.-'); hold on;  plot(locs1,DiffAverageActivity(locs1),'ro');
%     figure; plot(DiffActivityFiltered,'k.-'); hold on; plot(DiffAverageActivity,'g*-'); plot(locs1,DiffAverageActivity(locs1),'ro');
    
StartFrames = zeros(1,length(locs1));
PeakFrames  = zeros(1,length(locs1));
RealPeak    = false(1,length(locs1));
for pk_ind = 1:length(locs1)
    loc = locs1(pk_ind);
    pk  = pks1(pk_ind);

    % Find lower bound for activity initiation: the time when differential of AVERAGE activity begins to rise to peak
    %   First, find the time before peak, when 70% peak is reached. Second, go back until the first negative value 
    Ind70percent        = find(AverageActivity(1:loc)<0.7*pk,1,'last');
    if isempty(Ind70percent)
        continue
    end
    CurrentLowerBound = find(DiffAverageActivity(1:Ind70percent)<0,1,'last');
    if isempty(CurrentLowerBound)
        continue
    end
    MaxSlopeVal       = max(DiffAverageActivity(CurrentLowerBound:loc));

    % Find whether the peak is 'real'
    %   (1) Maximal average slope is large enough.    
    %   (2) Activity initiation occur within a reasonable time from peak time ('easy' constraint: 10 second limit)  
    %   (3) No previous peaks between current activity initiation and peak time   
    if (MaxSlopeVal > MinOfMaxSlopeVal) && ...
       ((loc-CurrentLowerBound)<MaxFramesBetweenPeakAndAverageSlopeInitiation) && ...
       ~any(ismember(locs1(1:(pk_ind-1)),CurrentLowerBound:loc))

        RealPeak(pk_ind) = true;

        % Find exact peak characteristics using the non-averaged vectors: DeltaFOverFFiltered, DiffActivityFiltered    
        %   (1) Find dynamic range from DeltaFOverF at the time when average vector is at lower bound and at peak time 
        %   (2) Activity initiation: go back from 25% dynamic range to when the differential is zero   
        %   (3) Peak time: First prominent peak from 75% dynamic range 
        MinVal       = DeltaFOverFFiltered(CurrentLowerBound);
        MaxVal       = DeltaFOverFFiltered(loc);
        DynamicRange = MaxVal-MinVal;

        Value25Percent = MinVal + DynamicRange*0.25;
        Ind25percent   = find(DeltaFOverFFiltered(CurrentLowerBound:loc)<Value25Percent,1,'last');
        IndStart       = find(DiffActivityFiltered(CurrentLowerBound:(CurrentLowerBound+Ind25percent))<0,1,'last');
        if isempty(IndStart) % NEED TO BE DEBUGGED...
            RealPeak(pk_ind) = false;
            continue
        end            
        StartFrame     = CurrentLowerBound+IndStart;
        StartFrames(pk_ind) = StartFrame;

        Value75Percent  = MinVal + DynamicRange*0.75;
        Ind75percent    = find(DeltaFOverFFiltered(CurrentLowerBound:loc)>Value75Percent,1,'first');
        Frame75percent  = CurrentLowerBound+Ind75percent;
        OptionalPeaksLocations = locs2; 
        OptionalPeaksLocations(OptionalPeaksLocations<Frame75percent) = 0;
        PeakFrame = OptionalPeaksLocations(find(OptionalPeaksLocations>0,1,'first'));            
        PeakFrames(pk_ind) = PeakFrame;        
    end     
end
% Keep only real peaks
StartFrames = StartFrames(RealPeak);
PeakFrames  = PeakFrames(RealPeak);

if plotme        
    figure('name','peak and deflection point detection'); 
    plot(DeltaFOverFFiltered,'k.-'); hold on;
    plot(StartFrames,DeltaFOverFFiltered(StartFrames),'ro','markerfacecolor','r'); 
    plot(PeakFrames,DeltaFOverFFiltered(PeakFrames),'go','markerfacecolor','g'); 
    plot(get(gca,'xlim'),[ThresholdValue ThresholdValue],'k:');
    plot(get(gca,'xlim'),[0 0],'k-');
end
    
return

function Data = FindActivityDuringDowngradient(Data, FrameRate, ThresholdValue)
NumOfTracks = size(Data.Stimulus_LowConfidence,1);
%%% Initialization
InitializationMat = zeros(size(Data.Stimulus_LowConfidence),'single')*NaN;
Data.DeltaFOverF_RelativeToDowngradientStart = InitializationMat;  % Activity Relative To DownGradient Start Frame 

InitializationMat  = false(size(Data.Stimulus_LowConfidence));
Data.DuringActivation                = InitializationMat;
Data.DuringActivation_OnlySharp      = InitializationMat;
Data.DuringActivation_NotOOB         = InitializationMat;
Data.Downgradient_ActivityStartsLow        = InitializationMat;
Data.Downgradient_ActivityStartsHigh       = InitializationMat;  
Data.Downgradient_NoActivation             = InitializationMat;
Data.Downgradient_InitialActivationFrame   = InitializationMat;
Data.Downgradient_DuringActivation         = InitializationMat;
Data.Downgradient_PeakActivationFrame      = InitializationMat;
Data.Downgradient_AfterPeakActivation      = InitializationMat;

Data.ActivationSegments.NotDuringDowngradient = cell(1,NumOfTracks); % average 0.5 second  after downgradient ends

Data.DowngradientSegments.InitialActivity        = cell(1,NumOfTracks); % average 0.5 second before downgradient begins
Data.DowngradientSegments.FinalActivity          = cell(1,NumOfTracks); % average 0.5 second  after downgradient ends
Data.DowngradientSegments.FinalRelativeActivity  = cell(1,NumOfTracks); % average 0.5 second  after downgradient ends
Data.DowngradientSegments.Activity_2sec_AfterSegment         = cell(1,NumOfTracks); % average 0.5 second  after downgradient ends
Data.DowngradientSegments.RelativeActivity_2sec_AfterSegment = cell(1,NumOfTracks); % average 0.5 second  after downgradient ends
Data.DowngradientSegments.Activation             = cell(1,NumOfTracks); % average 0.5 second  after downgradient ends
Data.DowngradientSegments.SharpActivation        = cell(1,NumOfTracks); % average 0.5 second  after downgradient ends

for tr=1:NumOfTracks  % Loop over Tracks
    
    Downgradient               = Data.Downgradient(tr,:);
    
    %% Calculate "during activation" frames
    ActivationInitiationFrames = Data.ActivationSegments.InitiationFrames{tr};
    ActivationPeakFrames       = Data.ActivationSegments.PeakFrames{tr};
    ActivationOOB              = Data.ActivationSegments.OOB{tr};
    ActivationSharp            = Data.ActivationSegments.SharpActivation{tr};
    Data.ActivationSegments.NotDuringDowngradient{tr} = false(size(ActivationInitiationFrames));
    
    for f_ind=1:length(ActivationInitiationFrames)
        Start = ActivationInitiationFrames(f_ind);
        End   = ActivationPeakFrames(f_ind);
        Data.DuringActivation(tr,Start:End)= true;
        if ~ActivationOOB(f_ind)
            Data.DuringActivation_NotOOB(tr,Start:End)= true;
        end
        if ActivationSharp(f_ind)
            Data.DuringActivation_OnlySharp(tr,Start:End)= true;
        end        
        if ~Downgradient(Start)
            Data.ActivationSegments.NotDuringDowngradient{tr}(f_ind) = true;
        end
    end
    %% events during downgradients
    Data.Downgradient_InitialActivationFrame(tr,ActivationInitiationFrames(~ActivationOOB)) = true;
    Data.Downgradient_PeakActivationFrame(tr,ActivationPeakFrames(~ActivationOOB))          = true;    
    Data.Downgradient_InitialActivationFrame(tr,:) = Downgradient & Data.Downgradient_InitialActivationFrame(tr,:);
    Data.Downgradient_PeakActivationFrame(tr,:)    = Downgradient & Data.Downgradient_PeakActivationFrame(tr,:);
    
    %% Characterisctics on the Downgeradient Segment level    
    % Initialization: input
    NeuronValueFiltered        = Data.NeuronValueFiltered(tr,:); 
    DeltaFOverFFiltered        = Data.DeltaFOverFFiltered(tr,:);
    DowngradientStartFrames    = Data.DowngradientSegments.StartFrames{tr};
    DowngradientEndFrames      = Data.DowngradientSegments.EndFrames{tr}  ;
    
    % Initialization: output
    Data.DowngradientSegments.InitialActivity{tr}         = zeros(1,length(DowngradientStartFrames))*NaN; % average 0.5 second before downgradient begins
    Data.DowngradientSegments.FinalActivity{tr}           = zeros(1,length(DowngradientStartFrames))*NaN; % average 0.5 second  after downgradient ends
    Data.DowngradientSegments.FinalRelativeActivity{tr}   = zeros(1,length(DowngradientStartFrames))*NaN; % average 0.5 second  after downgradient ends
    Data.DowngradientSegments.Activity_2sec_AfterSegment{tr}         = zeros(1,length(DowngradientStartFrames))*NaN; % average b/w 2-2.5 second  after downgradient ends
    Data.DowngradientSegments.RelativeActivity_2sec_AfterSegment{tr} = zeros(1,length(DowngradientStartFrames))*NaN; % average b/w 2-2.5 second  after downgradient ends
    Data.DowngradientSegments.Activation{tr}              = false(1,length(DowngradientStartFrames)); 
    Data.DowngradientSegments.SharpActivation{tr}         = false(1,length(DowngradientStartFrames)); 
    
        
    for f_ind = 1:length(DowngradientStartFrames)
        CurrentStartFrame = DowngradientStartFrames(f_ind);
        CurrentEndFrame   = DowngradientEndFrames(f_ind);
        
        InitialRawNeuronValue       = NeuronValueFiltered(CurrentStartFrame);
        CurrentRelativeDeltaFOverF  = (NeuronValueFiltered - InitialRawNeuronValue) / InitialRawNeuronValue;
        Data.DeltaFOverF_RelativeToDowngradientStart(tr,CurrentStartFrame:CurrentEndFrame) = ...
                ( NeuronValueFiltered(CurrentStartFrame:CurrentEndFrame) - InitialRawNeuronValue) / InitialRawNeuronValue;   
               
        % Mean 0.5 second before segment starts
        FirstFrame            = max([ 1  round(CurrentStartFrame-FrameRate*0.5)]);
        InitialActivity       = nanmean(DeltaFOverFFiltered(FirstFrame:CurrentStartFrame));
        % Mean 0.5 second after segment ends        
        LastFrame             = min([ length(DeltaFOverFFiltered)  round(CurrentStartFrame+FrameRate*0.5)]);
        FinalActivity         = nanmean(DeltaFOverFFiltered(CurrentEndFrame:LastFrame));
        FinalRelativeActivity = nanmean(CurrentRelativeDeltaFOverF(CurrentEndFrame:LastFrame));
        % Mean 0.5 second, starting from 2 seconds after segment ends        
        FirstFrame        = min([ length(DeltaFOverFFiltered)  round(CurrentStartFrame+FrameRate*2)]);
        LastFrame         = min([ length(DeltaFOverFFiltered)  round(CurrentStartFrame+FrameRate*2.5)]);        
        Activity_2sec_AfterSegment         = nanmean(DeltaFOverFFiltered(FirstFrame:LastFrame));
        RelativeActivity_2sec_AfterSegment = nanmean(CurrentRelativeDeltaFOverF(FirstFrame:LastFrame));
        
        Data.DowngradientSegments.InitialActivity{tr}(f_ind)            = InitialActivity;
        Data.DowngradientSegments.FinalActivity{tr}(f_ind)              = FinalActivity;
        Data.DowngradientSegments.FinalRelativeActivity{tr}(f_ind)      = FinalRelativeActivity;
        Data.DowngradientSegments.Activity_2sec_AfterSegment{tr}(f_ind)         = Activity_2sec_AfterSegment; 
        Data.DowngradientSegments.RelativeActivity_2sec_AfterSegment{tr}(f_ind) = RelativeActivity_2sec_AfterSegment; 
        
        if any(Data.DuringActivation_NotOOB(tr,CurrentStartFrame:CurrentEndFrame))
             Data.DowngradientSegments.Activation{tr}(f_ind) = true;
        end
        if any(Data.DuringActivation_OnlySharp(tr,CurrentStartFrame:CurrentEndFrame))
             Data.DowngradientSegments.SharpActivation{tr}(f_ind) = true;
        end
        
        %% characterisctics on the frame level
        AllFrames                                        = CurrentStartFrame:CurrentEndFrame;
        ActivationInSegment                              = Data.DuringActivation_NotOOB(tr,AllFrames);        
        Data.Downgradient_DuringActivation(tr,AllFrames) = ActivationInSegment;
        
        InitialActivationIndex = find(ActivationInSegment==1,1,'first');
        if isempty(InitialActivationIndex);
            Data.Downgradient_NoActivation(tr,AllFrames) = true;
        else
            Data.Downgradient_NoActivation(tr,AllFrames(1:(InitialActivationIndex-1)))     = true;
        end
            
        LastActivationIndex = find(ActivationInSegment==1,1,'last');
        if ~isempty(LastActivationIndex);
            Data.Downgradient_AfterPeakActivation(tr,AllFrames((LastActivationIndex+1):end)) = true;
        end  
        
        if InitialActivity > ThresholdValue   % segment starts with a high deltaFoverF value
            Data.Downgradient_ActivityStartsHigh(tr,AllFrames)= true ;
        else
            Data.Downgradient_ActivityStartsLow(tr,AllFrames) = true ;
        end
    end
    
end

return

function Data = ExtractConcentrationAndAdaptiveThreshold(Data, FrameRate, ExpInfo)
ConcentrationAtArenaTop     = ExpInfo.ConcentrationAtArenaTop;      % dilution
ConcentrationAtArenaBottom  = ExpInfo.ConcentrationAtArenaBottom;   % dilution

MicroM_ScalingFactor  = 11.16;  % dilution 1 = 11.16 Molars
% Soma parameters
Tao       = 17;                           % seconds
Beta      = 1/Tao;                        % =0.059 1/sec
K         = 5.5e-6/MicroM_ScalingFactor;  % = 4.93e-7 dilution  (=5.5 microM)

NumOfTracks = size(Data.Stimulus_LowConfidence,1);
NumOfFrames = size(Data.Stimulus_LowConfidence,2);
TimeVector  = (1:NumOfFrames)/FrameRate;

for tr = 1:NumOfTracks
    StimulusFiltered     = Data.StimulusFiltered(tr,:);
    AverageStimulus      = Data.AverageStimulus(tr,:);
    DiffStimulusFiltered = Data.DiffStimulusFiltered(tr,:);
    DiffAverageStimulus  = Data.DiffAverageStimulus(tr,:);
    
    C                     = (StimulusFiltered+1)/2 * (ConcentrationAtArenaTop - ConcentrationAtArenaBottom) + ConcentrationAtArenaBottom; % dC/dt = (dS/dt)*(dC/dS)
    C_Average             = (AverageStimulus+1)/2  * (ConcentrationAtArenaTop - ConcentrationAtArenaBottom) + ConcentrationAtArenaBottom; % dC/dt = (dS/dt)*(dC/dS)
    dCdFrame              = (DiffStimulusFiltered+1)/2  * (ConcentrationAtArenaTop - ConcentrationAtArenaBottom);                         % dC/dt = (dS/dt)*(dC/dS)
    dCdFrame_Average      = (DiffAverageStimulus+1)/2  * (ConcentrationAtArenaTop - ConcentrationAtArenaBottom);                          % dC/dt = (dS/dt)*(dC/dS) 
    dCdt                  = dCdFrame * FrameRate;
    dCdt_Average          = dCdFrame_Average * FrameRate;
    FoldChange            = dCdt_Average ./ C;
    
    AdaptiveThreshold = CalculateAdaptiveThreshold_inline(TimeVector, C, K, Beta);
    
    Data.C(tr,:)                 = C;
    Data.C_Average(tr,:)         = C_Average;
    Data.dCdt(tr,:)              = dCdt;
    Data.dCdt_Average(tr,:)      = dCdt_Average;
    Data.FoldChange(tr,:)        = FoldChange;      
    Data.AdaptiveThreshold(tr,:) = AdaptiveThreshold;    
end

return

function AdaptiveThreshold = CalculateAdaptiveThreshold_inline (Time, Concentration, K, Beta)

TimeDiffPerStep = diff(Time([1 2])); % seconds
C0              = Concentration(1);
VectorLength    = length(Time);

AdaptiveThreshold = K*(1-exp(-C0/K))*ones(1,VectorLength,'single');
for ind = 2:VectorLength
    CurrentConcentration   = Concentration(ind);
    AdaptiveThreshold(ind) = (AdaptiveThreshold(ind-1)*exp(-Beta*TimeDiffPerStep) + K*(1-exp(-CurrentConcentration/K))*(1- exp(-Beta*TimeDiffPerStep)));
end 

return

