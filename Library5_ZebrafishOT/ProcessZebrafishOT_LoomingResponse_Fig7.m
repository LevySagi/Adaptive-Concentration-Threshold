function ProcessZebrafishOT_LoomingResponse_Fig7
%  
%  Analyze zebrafish OT calcium responses to looming visual stimuli 
% 
%  Written by Sagi Levy, September 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic data processing  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data published in Figure 3E of Dunn et al. at Neuron, 
% Stimulus 1 is 'slow', stimulus 2 is 'fast', and stimulus 3 is 'medium' expansion speed.
load('D:\Fish Looming Data\Data_ZebrafishImaging_Raw.mat');
Stimulus  = stim_timecourse;            % looming stimulus time course for fish 1-8
Stimulus2 = stim_timecourse_fish_9_10;  % looming stimulus time course for fish 9 and 10
Time      = ((1:85)-16)/3;              % Data acquired @3Hz 
ZeroIndex = find(Time==0,1,'first');    % Stimulus initiation
EndIndex  = find(Time==15,1,'first');   % last stimulus time point @15sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate stimulus derivative and fold-change
FrameRate   = 3;
Derivative  = cell(1,3);
Derivative2 = cell(1,3);
FoldChange  = cell(1,3);
FoldChange2 = cell(1,3);
for speed_ind = 1:3
    Stimulus{speed_ind}(EndIndex:end) = NaN;    % Analysis of OT responses to visual stimulus (t<15sec), not to stimulus removal. 
    Stimulus2{speed_ind}(EndIndex:end)= NaN;    % Analysis of OT responses to visual stimulus (t<15sec), not to stimulus removal.  
    Derivative{speed_ind}  = [NaN diff(Stimulus{speed_ind}) * FrameRate]; 
    Derivative2{speed_ind} = [NaN diff(Stimulus2{speed_ind})* FrameRate]; 
    FoldChange{speed_ind}  = Derivative{speed_ind}  ./ Stimulus{speed_ind};  
    FoldChange2{speed_ind} = Derivative2{speed_ind} ./ Stimulus2{speed_ind};       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Matrices with all OT neuronal responses 
%  RawMAT/RawMAT2 corresponding to fish 1-8, and fish 9-10, respectively
RawMAT              = [];
FishIndexVec        = [];   % Vector with fish ID
CurrentNeuronNumber = 1;
for FishNumber = 1:8    
    CurrentMAT = all_resps{FishNumber};   % 3D matrix (time,neuron,stimulus speed)
    RawMAT     = [RawMAT, CurrentMAT];    
    FishIndexVec(CurrentNeuronNumber:(CurrentNeuronNumber+size(CurrentMAT,2)-1)) = FishNumber;
    CurrentNeuronNumber = CurrentNeuronNumber + size(CurrentMAT,2);
end
RawMAT = shiftdim(RawMAT,1);

RawMAT2             = [];
FishIndexVec2       = [];    % Vector with fish ID
CurrentNeuronNumber = 1;
for FishNumber = 9:10    
    CurrentMAT = all_resps{FishNumber};
    RawMAT2    = [RawMAT2, CurrentMAT];   
    FishIndexVec2(CurrentNeuronNumber:(CurrentNeuronNumber+size(CurrentMAT,2)-1)) = FishNumber;
    CurrentNeuronNumber = CurrentNeuronNumber + size(CurrentMAT,2);
end
RawMAT2 = shiftdim(RawMAT2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Delicate smoothing of calcium activity with LPS Filter (filter at @1.5Hz for sampling rate of 3Hz)
LPS = designfilt('lowpassfir', ...
  'PassbandFrequency',0.1,'StopbandFrequency',1.5, ...
  'PassbandRipple',1,'StopbandAttenuation',80, ...
  'DesignMethod','equiripple','SampleRate',3);
VectorLength = size(RawMAT,3);
% filter
SmoothedMAT  = shiftdim(filter(LPS,shiftdim(RawMAT,2)),1); 
SmoothedMAT2 = shiftdim(filter(LPS,shiftdim(RawMAT2,2)),1); 

% correction of filter delay
Delay        = mean(grpdelay(LPS)); % indices
RoundedDelay = round(Delay); 
TimeWithDelay = Time - Delay/3;
for speed_ind = 1:3
    CurrentMAT      = squeeze(SmoothedMAT(:,speed_ind,:));
    DelayedFixedMAT = interp1(TimeWithDelay, CurrentMAT', Time)';   % Interpolating values to the correct time indices                     
    SmoothedMAT(:,speed_ind,:) = DelayedFixedMAT;

    CurrentMAT      = squeeze(SmoothedMAT2(:,speed_ind,:));
    DelayedFixedMAT = interp1(TimeWithDelay, CurrentMAT', Time)';                        
    SmoothedMAT2(:,speed_ind,:) = DelayedFixedMAT;                
end        
SmoothedMAT(:,:,(VectorLength-RoundedDelay+1):VectorLength) = NaN;
SmoothedMAT2(:,:,(VectorLength-RoundedDelay+1):VectorLength)= NaN;

%%%%%% Examples of data filtering  %%%%%
% figure; 
% subplot(2,1,1);
% hold on; plot(Time,squeeze(RawMAT(100,2,:)),'k.'); hold on; plot(Time,squeeze(SmoothedMAT(100,2,:)),'k-'); 
% hold on; plot(Time,squeeze(RawMAT(100,3,:)),'.','color',[0.5 0 0]); hold on; plot(Time,squeeze(SmoothedMAT(100,3,:)),'-','color',[0.5 0 0]); 
% hold on; plot(Time,squeeze(RawMAT(100,1,:)),'r.'); hold on; plot(Time,squeeze(SmoothedMAT(100,1,:)),'r-'); 
% ylim([-0.3 3.2]); xlim([-1 15])
% line(get(gca,'xlim'),[0.1 0.1])
% legend({' Raw Activity, fast stimulus',' Smoothed Activity, fast stimulus',' Raw Activity, slow stimulus',' Smoothed Activity, slow stimulus'})
% subplot(2,1,2);
% hold on; plot(Time,squeeze(RawMAT(1100,2,:)),'k.'); hold on; plot(Time,squeeze(SmoothedMAT(1100,2,:)),'k-'); 
% hold on; plot(Time,squeeze(RawMAT(1100,3,:)),'.','color',[0.5 0 0]); hold on; plot(Time,squeeze(SmoothedMAT(1100,3,:)),'-','color',[0.5 0 0]); 
% hold on; plot(Time,squeeze(RawMAT(1100,1,:)),'r.'); hold on; plot(Time,squeeze(SmoothedMAT(1100,1,:)),'r-'); 
% ylim([-0.2 1.5]); xlim([-1 15])
% line(get(gca,'xlim'),[0.1 0.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Extract outliers: Neurons that are active prior to stimulus initiation or are non-responsive to any given stimuli
%%% (A) Neurons that are active prior to stimulus initiation %%%
%       strategy: Extract noise levels from raw matrices (standard deviation before stimulus initiates (index<=15))  
%                 Find traces that start significantly higher than the noise level
%   Fish 1-8
BeforeStimulusMAT        = RawMAT(:,:,1:15);  
SecondBeforeStimulusMAT  = RawMAT(:,:,13:15); % 1 second
BaseLines                = mean(SecondBeforeStimulusMAT,3);
BaseLinesFluctuations    = std(BeforeStimulusMAT,1,3);  
STDestimation            = mean(BaseLinesFluctuations(:));
% figure; hist(BaseLinesFluctuations(:),100)
TracesThatStartHigh      = BaseLines > 2 * STDestimation;        % corresponds to ~97 percentile of the Baseline fluctuations
PercentageStartHigh      = 100 * length(find(TracesThatStartHigh(:)))/length(TracesThatStartHigh(:));
NeuronWithTracesThatStartHigh = sum(TracesThatStartHigh,2)>0;    % High activity prior to stimulus initiation at any speed 

%   Fish 9-10
BeforeStimulusMAT2        = RawMAT2(:,:,1:15);  
SecondBeforeStimulusMAT2  = RawMAT2(:,:,13:15); % 1 second
BaseLines2                = mean(SecondBeforeStimulusMAT2,3);
BaseLinesFluctuations2    = std(BeforeStimulusMAT2,1,3);  
STDestimation2            = mean(BaseLinesFluctuations2(:));
% figure; hist(BaseLinesFluctuations(:),100)
TracesThatStartHigh2      = BaseLines2 > 2 * STDestimation2;     % corresponds to ~97 percentile of the Baseline fluctuations
PercentageStartHigh2      = 100 * length(find(TracesThatStartHigh2(:)))/length(TracesThatStartHigh2(:));
NeuronWithTracesThatStartHigh2 = sum(TracesThatStartHigh2,2)>0;  % High activity prior to stimulus initiation at any speed 

%%% (B) Non-reponsive neurons %%%
%       strategy: Extract max(SNR) of each neuron from smoothed matrices  
%                 Find non-reponsive neurons (SNR < 5)
%   Fish 1-8
BaseLines_Smoothed    = mean(SmoothedMAT(:,:,13:15),3);    % 1 second prior to stimulus initiation
MaxOfTrace            = max(SmoothedMAT(:,:,16:end),[],3); % Max after stimulus initiation
MaxOfNeuron           = max(MaxOfTrace,[],2);
BaseLinesFluctuations_Smoothed = std(SmoothedMAT(:,:,1:15),1,3);  
STDestimation_Smoothed         = mean(BaseLinesFluctuations_Smoothed(:));
DynamicRange          = repmat(MaxOfNeuron,1,3) - BaseLines_Smoothed;
MaxSNRperTrace        = DynamicRange ./ STDestimation_Smoothed;
MaxSNRperNeuron       = max(MaxSNRperTrace,[],2);
NeuronWithLowSNR      = MaxSNRperNeuron < 5;     % Threshold of SNR<5
PercentageLowSNR      = 100 * length(find(NeuronWithLowSNR(:)))/length(NeuronWithLowSNR(:));
% figure; hist(MaxSNRperNeuron(:),200); xlabel('SNR'); ylabel('# of neurons')

%   Fish 9-10
BaseLines_Smoothed2    = mean(SmoothedMAT2(:,:,13:15),3);    % 1 second
MaxOfTrace2            = max(SmoothedMAT2(:,:,16:end),[],3); % Max after stimulus initiation
MaxOfNeuron2           = max(MaxOfTrace2,[],2);
BaseLinesFluctuations_Smoothed2 = std(SmoothedMAT2(:,:,1:15),1,3);  
STDestimation_Smoothed2         = mean(BaseLinesFluctuations_Smoothed2(:));
DynamicRange2          = repmat(MaxOfNeuron2,1,3) - BaseLines_Smoothed2;
MaxSNRperTrace2        = DynamicRange2 ./ STDestimation_Smoothed2;
MaxSNRperNeuron2       = max(MaxSNRperTrace2,[],2);
NeuronWithLowSNR2      = MaxSNRperNeuron2 < 5;     % Threshold of SNR<5
PercentageLowSNR2      = 100 * length(find(NeuronWithLowSNR2(:)))/length(NeuronWithLowSNR2(:));

%%% (C) Neurons activated prior to stimulus initiation (A) or non-responsive (B) %%%  
NeuronToExclude       = NeuronWithTracesThatStartHigh | NeuronWithLowSNR;
NeglectMatrix         = repmat(NeuronToExclude,1,3,85);
PercentageToNeglect   = 100 * length(find(NeuronToExclude(:))) /length(NeuronToExclude(:));

NeuronToExclude2       = NeuronWithTracesThatStartHigh2 | NeuronWithLowSNR2;
NeglectMatrix2         = repmat(NeuronToExclude2,1,3,85);
PercentageToNeglect2   = 100 * length(find(NeuronToExclude2(:))) /length(NeuronToExclude2(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Normalization (relative to dynamic range) and outliers exclusion 
NormSmoothedMAT  = SmoothedMAT./repmat(MaxOfNeuron,1,3,85);   % Normalized by dynamic range [0 1]
NormSmoothedMAT2 = SmoothedMAT2./repmat(MaxOfNeuron2,1,3,85); % Normalized by dynamic range [0 1]
ProcessedMAT     = NormSmoothedMAT(~NeuronToExclude,:,:);     % ProcessedMat = smoothed + normalized by dynamic range + exclude outliers
ProcessedMAT2    = NormSmoothedMAT2(~NeuronToExclude2,:,:);   % ProcessedMat = smoothed + normalized by dynamic range + exclude outliers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Stimulus, Derivative and FoldChange 3D matrices (neuron,speed,time)
%  fish 1-8
TotalNumberOfNeurons = size(ProcessedMAT,1);  
StimulusBigMat       = ones(TotalNumberOfNeurons,3,size(ProcessedMAT,3));
DerivativeBigMat     = ones(TotalNumberOfNeurons,3,size(ProcessedMAT,3));
FoldChangeBigMat     = ones(TotalNumberOfNeurons,3,size(ProcessedMAT,3));
for speed_ind = 1:3
    StimulusBigMat(:,speed_ind,:)   = repmat(Stimulus{speed_ind},   TotalNumberOfNeurons, 1);    
    DerivativeBigMat(:,speed_ind,:) = repmat(Derivative{speed_ind}, TotalNumberOfNeurons, 1);    
    FoldChangeBigMat(:,speed_ind,:) = repmat(FoldChange{speed_ind}, TotalNumberOfNeurons, 1);    
end

%  fish 9-10
TotalNumberOfNeurons2 = size(ProcessedMAT2,1);  
StimulusBigMat2       = ones(TotalNumberOfNeurons2,3,size(ProcessedMAT2,3));
DerivativeBigMat2     = ones(TotalNumberOfNeurons2,3,size(ProcessedMAT2,3));
FoldChangeBigMat2     = ones(TotalNumberOfNeurons2,3,size(ProcessedMAT2,3));
for speed_ind = 1:3
    StimulusBigMat2(:,speed_ind,:)   = repmat(Stimulus2{speed_ind},   TotalNumberOfNeurons2, 1);    
    DerivativeBigMat2(:,speed_ind,:) = repmat(Derivative2{speed_ind}, TotalNumberOfNeurons2, 1);    
    FoldChangeBigMat2(:,speed_ind,:) = repmat(FoldChange2{speed_ind}, TotalNumberOfNeurons2, 1);    
end

% All fish
ProcessedMAT_All     = [ProcessedMAT;     ProcessedMAT2];
StimulusBigMat_All   = [StimulusBigMat;   StimulusBigMat2];
DerivativeBigMat_All = [DerivativeBigMat; DerivativeBigMat2];
FoldChangeBigMat_All = [FoldChangeBigMat; FoldChangeBigMat2];

%%%%%  Plot basline stats   %%%%%
% BeforeStimulusInitiation           = ZeroIndex-1;
% BaseLinesFluctuations_All_Smoothed = std(ProcessedMAT_All(:,:,1:BeforeStimulusInitiation),1,3);  
% STDestimation_All_Smoothed         = mean(BaseLinesFluctuations_All_Smoothed(:));
% STDestimation_All_Smoothed_STD     = std(BaseLinesFluctuations_All_Smoothed(:));
% SNR4_All = 4* STDestimation_All_Smoothed;
% BaseLines_All_Smoothed = mean(ProcessedMAT_All(:,:,BeforeStimulusInitiation-2:BeforeStimulusInitiation),3);  

% figure; hist(BaseLines_All_Smoothed(:),300); xlim([0 0.4]); title('baseline amplitude')
% figure; hist(BaseLinesFluctuations_All_Smoothed(:),300); xlim([0 0.4]); title('baseline fluctuations')
% figure; hist((BaseLines_All_Smoothed(:)+BaseLinesFluctuations_All_Smoothed(:)),300); xlim([0 0.4]); title('baseline amplitude + fluctuations')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Increased temporal resolution (30Hz) of stimulus, derivative and fold-change
% Regular temporal resolution (3Hz) of stimulus, derivative and fold-change, from stimulus initiation to stimulus termination 
LengthOfStimulusVec     = (EndIndex-ZeroIndex+1);
StimulusMat1            = ones(3,LengthOfStimulusVec);
StimulusMat2            = ones(3,LengthOfStimulusVec);
DerivativeMat1          = ones(3,LengthOfStimulusVec);
DerivativeMat2          = ones(3,LengthOfStimulusVec);
% Increased temporal resolution (30Hz) of stimulus, derivative and fold-change, from stimulus initiation to stimulus termination 
LengthOfLongStimulusVec = (EndIndex-ZeroIndex)*10+1;
DeltaT                  = 1/30;         % x10 increase temporal resolution 
TimeLong                = 0:DeltaT:15;
StimulusMatLong1        = ones(3,LengthOfLongStimulusVec);
StimulusMatLong2        = ones(3,LengthOfLongStimulusVec);
DerivativeMatLong1      = ones(3,LengthOfLongStimulusVec);
DerivativeMatLong2      = ones(3,LengthOfLongStimulusVec);

for speed_ind = 1:3        
    StimulusMat1(speed_ind,:)     = Stimulus{speed_ind}([ZeroIndex:EndIndex-1, EndIndex-1] );
    StimulusMat2(speed_ind,:)     = Stimulus2{speed_ind}([ZeroIndex:EndIndex-1, EndIndex-1]);    
    StimulusMatLong1(speed_ind,:) = interp1( StimulusMat1(speed_ind,:), 1:0.1:LengthOfStimulusVec,'pchip'); % x10 increase temporal resolution 
    StimulusMatLong2(speed_ind,:) = interp1( StimulusMat2(speed_ind,:), 1:0.1:LengthOfStimulusVec,'pchip'); % x10 increase temporal resolution 
    
    DerivativeMat1(speed_ind,:)     = Derivative{speed_ind}([ZeroIndex:EndIndex-1, EndIndex-1]);
    DerivativeMat2(speed_ind,:)     = Derivative2{speed_ind}([ZeroIndex:EndIndex-1, EndIndex-1]);    
    DerivativeMatLong1(speed_ind,:) = interp1( DerivativeMat1(speed_ind,:), 1:0.1:LengthOfStimulusVec,'pchip'); % x10 increase temporal resolution 
    DerivativeMatLong2(speed_ind,:) = interp1( DerivativeMat2(speed_ind,:), 1:0.1:LengthOfStimulusVec,'pchip'); % x10 increase temporal resolution 
end
FoldChangeLong1 = DerivativeMatLong1 ./ StimulusMatLong1;
FoldChangeLong2 = DerivativeMatLong2 ./ StimulusMatLong2;

% Generate 3D matrices of increased temporal resolution stimulus, derivative and fold-change. Matrix dimensions = (neuron, speed, time)  
TotalNumberOfNeurons = size(ProcessedMAT,1);  
StimulusBigMat1      = ones(TotalNumberOfNeurons,3,LengthOfLongStimulusVec,'single');
DerivativeBigMat1    = ones(TotalNumberOfNeurons,3,LengthOfLongStimulusVec,'single');
FoldChangeBigMat1    = ones(TotalNumberOfNeurons,3,LengthOfLongStimulusVec,'single');
for speed_ind = 1:3
    StimulusBigMat1(:,speed_ind,:)   = repmat(StimulusMatLong1(speed_ind,:),   TotalNumberOfNeurons, 1);    
    DerivativeBigMat1(:,speed_ind,:) = repmat(DerivativeMatLong1(speed_ind,:), TotalNumberOfNeurons, 1);    
    FoldChangeBigMat1(:,speed_ind,:) = repmat(FoldChangeLong1(speed_ind,:),    TotalNumberOfNeurons, 1);    
end

TotalNumberOfNeurons2 = size(ProcessedMAT2,1);  
StimulusBigMat2       = ones(TotalNumberOfNeurons2,3,LengthOfLongStimulusVec,'single');
DerivativeBigMat2     = ones(TotalNumberOfNeurons2,3,LengthOfLongStimulusVec,'single');
FoldChangeBigMat2     = ones(TotalNumberOfNeurons2,3,LengthOfLongStimulusVec,'single');
for speed_ind = 1:3
    StimulusBigMat2(:,speed_ind,:)   = repmat(StimulusMatLong2(speed_ind,:),   TotalNumberOfNeurons2, 1);    
    DerivativeBigMat2(:,speed_ind,:) = repmat(DerivativeMatLong2(speed_ind,:), TotalNumberOfNeurons2, 1);    
    FoldChangeBigMat2(:,speed_ind,:) = repmat(FoldChangeLong2(speed_ind,:),    TotalNumberOfNeurons2, 1);    
end

StimulusBigMat_All   = [StimulusBigMat1;   StimulusBigMat2];
DerivativeBigMat_All = [DerivativeBigMat1; DerivativeBigMat2];
FoldChangeBigMat_All = [FoldChangeBigMat1; FoldChangeBigMat2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Increased temporal resolution (30Hz) of neuronal activity
MAT                   = ProcessedMAT_All(:,:,ZeroIndex:EndIndex);
ProcessedMAT_Long_All = ones(size(MAT,1),3,LengthOfLongStimulusVec,'single');

for neuron_ind = 1:size(MAT,1)
    for speed_ind = 1:3        
        ProcessedMAT_Long_All(neuron_ind,speed_ind,:) = interp1( squeeze(MAT(neuron_ind,speed_ind,:)), 1:0.1:LengthOfStimulusVec,'spline');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Compute responsivenese, latency and stimulus features at the time of neuronal activation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
[Responsiveness, ActivationTimes, ActivationIndices] = CalculateActivationTime(TimeLong, ProcessedMAT_Long_All); 
StimulusAtActivationTime   = CalculateSignalPropertiesAtActivationTime(ActivationIndices, StimulusBigMat_All);
DerivativeAtActivationTime = CalculateSignalPropertiesAtActivationTime(ActivationIndices, DerivativeBigMat_All);
FoldChangeAtActivationTime = CalculateSignalPropertiesAtActivationTime(ActivationIndices, FoldChangeBigMat_All);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% ACT model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Explore ACT model parameters
ACTmodelStructure = FishACT_screen_inline(TimeLong, StimulusMatLong1, StimulusMatLong2) ;  % Screen over a small range of parameters, for parameter tuning   

%% Fitted adaptive threshold model parameters
%  B = (Alpha-(2*A)/X0);  C = (A/X0^2);  
X0 = 70;  A = 31;  Alpha = 1;  tao = 2;   
% Find corresponding indices in screen structure
A_vec     = ACTmodelStructure.ScreenParams.A_vec;
Alpha_vec = ACTmodelStructure.ScreenParams.Alpha_vec;
X0_vec    = ACTmodelStructure.ScreenParams.X0_vec;
Tao_vec   = ACTmodelStructure.ScreenParams.Tao_vec;

X0_ind    = find(X0_vec>(X0-1e-5),1,'first');
A_ind     = find(A_vec>(A-1e-5),1,'first');
Tao_ind   = find(Tao_vec>(tao-1e-5),1,'first');
Alpha_ind = find(Alpha_vec>(Alpha-1e-5),1,'first');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% save processed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
save('D:\Fish Looming Data\Data_ZebrafishImaging_Processed.mat','-v7.3');

return

%% Basic processing 
function [Responsiveness, ActivationTimes, ActivationIndices] = CalculateActivationTime(Time, MAT)

ThresholdForFindingPeak                   = 0.2;  % if empty, use minimum of 0.1
ThresholdForFindingPeak_Try2              = 0.1;  % if 'ThresholdForFindingPeak' give no result, use minimum of 0.1
ThresholdForFindingActivationTimeFromPeak = 0.05; 

ZeroIndex   = find(Time==0,1,'first');
EndIndex    = find(Time==15,1,'first');

ThresholdMAT_PeakDetection                        = MAT >= ThresholdForFindingPeak;
ThresholdMAT_PeakDetection(:,:,1:ZeroIndex)       = false; % Neglect events before stimulus begins
ThresholdMAT_PeakDetection(:,:,(EndIndex+1):end)  = false; % Neglect events after stimulus ends
Diff_ThresholdMAT_PeakDetection                   = diff(ThresholdMAT_PeakDetection,1,3);

ThresholdMAT_PeakDetection_Try2                        = MAT >= ThresholdForFindingPeak_Try2;
ThresholdMAT_PeakDetection_Try2(:,:,1:ZeroIndex)       = false; % Neglect events before stimulus begins
ThresholdMAT_PeakDetection_Try2(:,:,(EndIndex+1):end)  = false; % Neglect events after stimulus ends
Diff_ThresholdMAT_PeakDetection_Try2                   = diff(ThresholdMAT_PeakDetection_Try2,1,3);

ThresholdMAT_ActivationTime                       = MAT >= ThresholdForFindingActivationTimeFromPeak;
ThresholdMAT_ActivationTime(:,:,1:ZeroIndex)      = false; % Neglect events before stimulus begins
Diff_ThresholdMAT_ActivationTime                  = diff(ThresholdMAT_ActivationTime,1,3);

NumberOfNeurons        = size(MAT,1);
ActivationTimes        = zeros(NumberOfNeurons,3)*NaN;
ActivationIndices      = zeros(NumberOfNeurons,3)*NaN;
for n_ind = 1:NumberOfNeurons
    for speed_ind = 1:3
        Index_MiddleOfRisingToPeak = find(squeeze(Diff_ThresholdMAT_PeakDetection(n_ind,speed_ind,:))==1,1,'first')+1; 
        if isempty(Index_MiddleOfRisingToPeak)
            Index_MiddleOfRisingToPeak = find(squeeze(Diff_ThresholdMAT_PeakDetection_Try2(n_ind,speed_ind,:))==1,1,'first')+1;            
        end
        if ~isempty(Index_MiddleOfRisingToPeak) && Index_MiddleOfRisingToPeak~= size(MAT,3)
            Index_ActivationTime   = find(squeeze(Diff_ThresholdMAT_ActivationTime(n_ind,speed_ind,1:Index_MiddleOfRisingToPeak))==1,1,'last');        
            ActivationTimes(n_ind,speed_ind)           = Time(Index_ActivationTime);                       
            ActivationIndices(n_ind,speed_ind)         = Index_ActivationTime;                       
        end
    end    
end

Responsiveness = ~isnan(ActivationTimes); 

return

function SignalPropertyAtActivationTime = CalculateSignalPropertiesAtActivationTime(ActivationIndices, SignalPropertyMAT)

NumberOfNeurons                = size(ActivationIndices,1);
SignalPropertyAtActivationTime = ActivationIndices * NaN; % initialization
for n_ind = 1:NumberOfNeurons
    for speed_ind = 1:3
        index                 = ActivationIndices(n_ind,speed_ind);
        if ~isnan(index)
            SignalPropertyAtActivationTime(n_ind, speed_ind) = SignalPropertyMAT(n_ind, speed_ind, index);
        end
    end    
end

return

%% ACT model analysis
function StructureOut = FishACT_screen_inline(TimeLong, StimulusMatLong1, StimulusMatLong2) 

% Polynomial steady state function is inspired by previous studies of photoreceptor responses and behavior responses to visual stimulus,    
%   e.g. Burkhardt, 1994; Donner et al., 1990; Frishman and Sieving, 1995; Rieke and Rudd, 2009  
%   Y = A + (Alpha-(2*A)/X0)*X + (A/X0^2)*X^2 ;

%% Screen parameters
% A_vec     = [20 23 25:36 40 50];
% Alpha_vec = [1:0.01:1.05 1.07 1.1 1.12 1.15 1.2 1.3 1.5 2 3 5 10];
% X0_vec    = [10 20 40 70 100 150 200 300];
% Tao_vec   = [0.1 0.4 0.7 1 1.2 1.5 2 3 5 7 10 15 20 50 200];
% 
% A_vec     = [20 23 25:36 40 50];
% Alpha_vec = [1:0.01:1.1 1.12 1.15 1.2 1.3 1.4 1.5 2 5 10];
% X0_vec    = [10 20 40 70 100 150 200 300 2000 1e4];
% Tao_vec   = [0.1 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.5 3 10 15 20 200];

A_vec     = [20 23 25:36 40 50];
Alpha_vec = [1:0.02:1.14 1.17 1.2 1.3 1.4 1.5 2 5 10];
X0_vec    = [10 20 40 50 60 70 80 100 200 300 1e4];
Tao_vec   = [0.1 0.5 1 1.1 1.2 1.5 1.7 2 2.3 2.6 3 10 15 20 200];

ScreenParams.TimeLong         = TimeLong;
ScreenParams.StimulusMatLong1 = StimulusMatLong1;
ScreenParams.StimulusMatLong2 = StimulusMatLong2;
ScreenParams.A_vec            = A_vec;
ScreenParams.Alpha_vec        = Alpha_vec;
ScreenParams.X0_vec           = X0_vec;
ScreenParams.Tao_vec          = Tao_vec;
ScreenParams.SteadyStatefunction = 'Y = A + (Alpha-(2*A)/X0)*X + (A/X0^2)*X^2'; 

LengthOfLongStimulusVec = length(TimeLong);
DeltaT                  = diff(TimeLong([1 2])); % seconds


%% Run screen
ACT_MAT1                      = ones(length(A_vec),length(Alpha_vec),length(X0_vec),length(Tao_vec),3,LengthOfLongStimulusVec,'single')*NaN; % protocol 1
ACT_MAT2                      = ones(length(A_vec),length(Alpha_vec),length(X0_vec),length(Tao_vec),3,LengthOfLongStimulusVec,'single')*NaN; % protocol 2
PredictedActivationIndex_MAT1 = ones(length(A_vec),length(Alpha_vec),length(X0_vec),length(Tao_vec),3,'single')*NaN;     % protocol 1
PredictedActivationIndex_MAT2 = ones(length(A_vec),length(Alpha_vec),length(X0_vec),length(Tao_vec),3,'single')*NaN;     % protocol 2
PredictedActivationTime_MAT1  = ones(length(A_vec),length(Alpha_vec),length(X0_vec),length(Tao_vec),3,'single')*NaN;     % protocol 1
PredictedActivationTime_MAT2  = ones(length(A_vec),length(Alpha_vec),length(X0_vec),length(Tao_vec),3,'single')*NaN;     % protocol 2 
StimulusAtActivationTime_MAT1 = ones(length(A_vec),length(Alpha_vec),length(X0_vec),length(Tao_vec),3,'single')*NaN;     % protocol 1
StimulusAtActivationTime_MAT2 = ones(length(A_vec),length(Alpha_vec),length(X0_vec),length(Tao_vec),3,'single')*NaN;     % protocol 2 

plotme = false;
for A_ind = 1:length(A_vec)    
    A  = A_vec(A_ind);
    T0 = A;
    for Alpha_ind = 1:length(Alpha_vec)
        Alpha = Alpha_vec(Alpha_ind);            
        for X0_ind = 1:length(X0_vec)                
            X0 = X0_vec(X0_ind);       
            for Tao_ind = 1:length(Tao_vec)            
                Tao = Tao_vec(Tao_ind);            
                for speed_ind = 1:3
                    % Compute ACT for stimulus in protocol 1
                    CurrentACT = ComputeFishACT_inline (StimulusMatLong1(speed_ind,:), DeltaT, A, Alpha, X0, Tao, T0, plotme);   
                    ACT_MAT1(A_ind, Alpha_ind, X0_ind, Tao_ind, speed_ind, :) = CurrentACT; 
                    Index      = find( CurrentACT <= StimulusMatLong1(speed_ind,:) ,1,'first');
                    if ~isempty(Index)
                        PredictedActivationIndex_MAT1(A_ind, Alpha_ind, X0_ind, Tao_ind, speed_ind) = Index;
                        PredictedActivationTime_MAT1( A_ind, Alpha_ind, X0_ind, Tao_ind, speed_ind) = TimeLong(Index);
                        StimulusAtActivationTime_MAT1(A_ind, Alpha_ind, X0_ind, Tao_ind, speed_ind) = StimulusMatLong1(speed_ind,Index);                    
                    end

                    % Compute ACT for stimulus in protocol 2
                    CurrentACT = ComputeFishACT_inline (StimulusMatLong2(speed_ind,:), DeltaT, A, Alpha, X0, Tao, T0, plotme);                
                    ACT_MAT2(A_ind, Alpha_ind, X0_ind, Tao_ind, speed_ind, :) = CurrentACT;             
                    Index      = find( CurrentACT <= StimulusMatLong2(speed_ind,:) ,1,'first');
                    if ~isempty(Index)
                        PredictedActivationIndex_MAT2(A_ind, Alpha_ind, X0_ind, Tao_ind, speed_ind) = Index;
                        PredictedActivationTime_MAT2( A_ind, Alpha_ind, X0_ind, Tao_ind, speed_ind) = TimeLong(Index);
                        StimulusAtActivationTime_MAT2(A_ind, Alpha_ind, X0_ind, Tao_ind, speed_ind) = StimulusMatLong2(speed_ind,Index);                    
                    end

                end        
            end
        end
    end    
end

%% Assign to output structure
StructureOut.ScreenParams                  = ScreenParams;
StructureOut.PredictedActivationTime_MAT1  = PredictedActivationTime_MAT1;
StructureOut.PredictedActivationTime_MAT2  = PredictedActivationTime_MAT2;
StructureOut.PredictedActivationIndex_MAT1 = PredictedActivationIndex_MAT1;
StructureOut.PredictedActivationIndex_MAT2 = PredictedActivationIndex_MAT2;
StructureOut.StimulusAtActivationTime_MAT1 = StimulusAtActivationTime_MAT1;
StructureOut.StimulusAtActivationTime_MAT2 = StimulusAtActivationTime_MAT2;
StructureOut.ACT_MAT1                      = ACT_MAT1;
StructureOut.ACT_MAT2                      = ACT_MAT2;

return

function Y = ComputeFishACT_inline (X, TimeDiffPerStep, A, Alpha, X0, Tao, Y0, plotme)
% X = Stimulus, Y = Adaptive-threshold
% SteadyStateFunction:   Y = A + (Alpha-(2*A)/X0)*X + (A/X0^2)*X^2 ;
if ~exist('Y0','var')|| isempty(Y0)
    Y0 = A + (Alpha-(2*A)/X0)*X(1) + (A/X0^2)*X(1)^2 ;  % If initial conditions are not given, assume that X(-inf<t<0)=X(1)
end

% TimeDiffPerStep = diff(Time([1 2])); % seconds
Beta         = 1/Tao;
VectorLength = length(X);

Y = Y0 * ones(1,VectorLength,'single');
for ind = 2:VectorLength
    CurrentX = X(ind);
    Y(ind)   = Y(ind-1)*exp(-Beta*TimeDiffPerStep) + ( A + (Alpha-(2*A)/X0)*CurrentX + (A/X0^2)*CurrentX^2 )* (1- exp(-Beta*TimeDiffPerStep));
end 

if exist('plotme','var')&& plotme
    figure; plot(X,Y);        
end

return




