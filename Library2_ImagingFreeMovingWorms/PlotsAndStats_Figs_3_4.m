function PlotsAndStats_Figs_3_4(Data)

% September 2019, Sagi Levy

%% Load Data, if needed
PathName       = 'F:\';
SaveFileName   = [PathName,'GradientImagingMetaAnalysis.mat'];
if ~exist('Data','var')
    load(SaveFileName,'Data')  
end
FrameRate = 30;

%% Number of repeats/frames
Y = sum(Data.GradientFramesMatrix,2);
WormIndices = [ 5     6     7     8    19    20]; 
NumOfWorms  = length(WormIndices); Minutes = sum(Y(WormIndices))/FrameRate/60;
disp(['1e-5 repeats: ',num2str(NumOfWorms), ' worms, ',num2str(Minutes),' minutes, ',num2str(sum(Y(WormIndices))),' frames']);
WormIndices = 11:18; 
NumOfWorms  = length(WormIndices); Minutes = sum(Y(WormIndices))/FrameRate/60;
disp(['1e-6 repeats: ',num2str(NumOfWorms), ' worms, ',num2str(Minutes),' minutes, ',num2str(sum(Y(WormIndices))),' frames']);
WormIndices = 1:4; 
NumOfWorms  = length(WormIndices); Minutes = sum(Y(WormIndices))/FrameRate/60;
disp(['1e-7 repeats: ',num2str(NumOfWorms), ' worms, ',num2str(Minutes),' minutes, ',num2str(sum(Y(WormIndices))),' frames']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Figure 3E and Figures S4C-E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
[Stats, HistogramStructure] = Plot_ModelPredictorsComparison(Data);    % Figure 3E, S4C-E, including SEM by bootstrap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Plot activity and behavior matrices and CDFs when aligned by events-of-choice
%  Events-of-choice includes:
%       Measured activity initiation events and measured no-initiation events 
%           a) Regardless of behavior 
%           b) When behavior = Forward (Figures 4A and 4B)
%           c) When behavior = Reverse (Figure S4F)
%       Predicted activity initiation events and predicted no-initiation events  
%           a) Regardless of behavior (Figure 3D)
%           b) When behavior = Forward (Figures 4A and 4B)
%           c) When behavior = Reverse (Figure S4F)
%       Control: all downgradient frames activity initiation events and predicted no-initiation events  
%           a) Regardless of behavior (Figure 3D)
%           b) When behavior = Forward (Figures 4A and 4B)
%           c) When behavior = Reverse (Figure S4F)
%       Measured reversals initiation
%           Plot mean activity as a function of ACT predictor (Figures 4C and S4H)
%       Measured downgradient navigation termination    
%           Plot mean activity as a function of ACT predictor (Figures 4C and S4H)
%
%% Initialization
FramesWindow = (-5*FrameRate):(5*FrameRate); % [-5,5] seconds
ZeroIndex    = find(FramesWindow==0);
PlotInfo.NeuronClimit   = [-0.8 2];  % for actual activity full range [0 5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Matrices aligned by actual or predicted activity (Figures 3D, 4A, 4B and S4F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
PlotInfo.XLIMITS        = [0 5];

%% PREDICTED ACTIVITY:  Align data by events of crossing (or not crossing) the adaptive threshold
%%% Generate matrices and save processed data in structures  %%%%
StructureOut = CalculateMatricesAlignedByThresholdCrossingEvents(Data, FramesWindow);
% Assign F0 to be the activity at t=0 (instead of average before the event), before plotting activity figures 
StructureOut.CrossingThreshold.EventPastNeuronValue                      = StructureOut.CrossingThreshold.NeuronValue(:,ZeroIndex);
StructureOut.CrossingAboveThreshold.EventPastNeuronValue                 = StructureOut.CrossingAboveThreshold.NeuronValue(:,ZeroIndex);
StructureOut.DowngradientAboveThreshold.EventPastNeuronValue             = StructureOut.DowngradientAboveThreshold.NeuronValue(:,ZeroIndex);
StructureOut.DowngradientAboveThresholdAtLeast5sec.EventPastNeuronValue  = StructureOut.DowngradientAboveThresholdAtLeast5sec.NeuronValue(:,ZeroIndex);
StructureOut.AllFrames_RegardlessOfThreshold.EventPastNeuronValue        = StructureOut.AllFrames_RegardlessOfThreshold.NeuronValue(:,ZeroIndex);
StructureOut.AllFramesAboveThreshold.EventPastNeuronValue                = StructureOut.AllFramesAboveThreshold.NeuronValue(:,ZeroIndex);
StructureOut.AllFramesAboveThresholdAtLeast5sec.EventPastNeuronValue     = StructureOut.AllFramesAboveThresholdAtLeast5sec.NeuronValue(:,ZeroIndex);
StructureOut.DowngradientAboveThresholdAtLeast5sec_Binning10.EventPastNeuronValue = StructureOut.DowngradientAboveThresholdAtLeast5sec_Binning10.NeuronValue(:,(5*FrameRate/10+1));
StructureOut.DowngradientAboveThresholdAtLeast5sec_Binning30.EventPastNeuronValue = StructureOut.DowngradientAboveThresholdAtLeast5sec_Binning30.NeuronValue(:,(5*FrameRate/30+1));

%%%% Figure 3D %%%%
PlotInfo.PlotActivity          = true;
PlotInfo.PlotBehavior          = false;
MinFramesInBoundFromEventStart = 0;           % Analyze all events
PlotInfo.FigureTitle = 'Predict activation';    PlotMatrices_AlignByActivityOrPrediction(StructureOut.CrossingThreshold, MinFramesInBoundFromEventStart, PlotInfo);                    % Predicted activity
PlotInfo.FigureTitle = 'Predict no-activation'; PlotMatrices_AlignByActivityOrPrediction(StructureOut.DowngradientAboveThresholdAtLeast5sec, MinFramesInBoundFromEventStart, PlotInfo);% Predicted no-activity

%%%% Figures 4A, 4B and S4F %%%%
PlotInfo.PlotActivity          = false;
PlotInfo.PlotBehavior          = true;
MinFramesInBoundFromEventStart = 5*FrameRate; % Do not analyze behavior of animals that are out of bounds

% Predicted activity
PlotInfo.FigureTitle = 'Predict activation';    
[CDF_Fwd, CDF_Rev, SortedTimes_Fwd, SortedTimes_Rev] = PlotMatrices_AlignByActivityOrPrediction(StructureOut.CrossingThreshold, MinFramesInBoundFromEventStart, PlotInfo);
CDFs.Predicted.Fwd = CDF_Fwd;
CDFs.Predicted.Rev = CDF_Rev;
Repeats.Predicted.Fwd = length(SortedTimes_Fwd);
Repeats.Predicted.Rev = length(SortedTimes_Rev);
Times.Predicted.Fwd   = SortedTimes_Fwd;
Times.Predicted.Rev   = SortedTimes_Rev;

% Predicted no-activity
PlotInfo.FigureTitle = 'Predict no-activation'; 
[CDF_Fwd, CDF_Rev, SortedTimes_Fwd, SortedTimes_Rev] = PlotMatrices_AlignByActivityOrPrediction(StructureOut.DowngradientAboveThresholdAtLeast5sec, MinFramesInBoundFromEventStart, PlotInfo);
CDFs.NotPredicted.Fwd    = CDF_Fwd;
CDFs.NotPredicted.Rev    = CDF_Rev;
Repeats.NotPredicted.Fwd = length(SortedTimes_Fwd);
Repeats.NotPredicted.Rev = length(SortedTimes_Rev);
Times.NotPredicted.Fwd   = SortedTimes_Fwd;
Times.NotPredicted.Rev   = SortedTimes_Rev;

% Control- all frames (regardless of the adaptive threshold or locomotion)
PlotInfo.FigureTitle = 'Control (all frames)'; 
[CDF_Fwd, CDF_Rev, SortedTimes_Fwd, SortedTimes_Rev] = PlotMatrices_AlignByActivityOrPrediction(StructureOut.AllFrames_RegardlessOfThreshold, MinFramesInBoundFromEventStart, PlotInfo);
CDFs.Baseline.Fwd    = CDF_Fwd;
CDFs.Baseline.Rev    = CDF_Rev;
Repeats.Baseline.Fwd = length(SortedTimes_Fwd);
Repeats.Baseline.Rev = length(SortedTimes_Rev);
Times.Baseline.Fwd   = SortedTimes_Fwd;
Times.Baseline.Rev   = SortedTimes_Rev;

% Control #2: see that analysis of no-activity prediction are not sensitive to the frame rate, use 10x binning 
PlotInfo.PlotActivity          = true;
PlotInfo.PlotBehavior          = false;
MinFramesInBoundFromEventStart = 0;           % Analyze all events
BinningFactor = 10; PlotInfo.FigureTitle = ['Predict no-activation, binning x',num2str(BinningFactor)]; 
PlotMatrices_AlignByActivityOrPrediction(StructureOut.DowngradientAboveThresholdAtLeast5sec_Binning10, MinFramesInBoundFromEventStart/BinningFactor, PlotInfo);
BinningFactor = 30; PlotInfo.FigureTitle = ['Predict no-activation, binning x',num2str(BinningFactor)]; 
PlotMatrices_AlignByActivityOrPrediction(StructureOut.DowngradientAboveThresholdAtLeast5sec_Binning30, MinFramesInBoundFromEventStart/BinningFactor, PlotInfo);

PlotInfo.PlotActivity          = false;
PlotInfo.PlotBehavior          = true;
MinFramesInBoundFromEventStart = 5*FrameRate; % Do not analyze behavior of animals that are out of bounds
BinningFactor = 10; PlotInfo.FigureTitle = ['Predict no-activation, binning x',num2str(BinningFactor)]; 
PlotMatrices_AlignByActivityOrPrediction(StructureOut.DowngradientAboveThresholdAtLeast5sec_Binning10, MinFramesInBoundFromEventStart/BinningFactor, PlotInfo);
BinningFactor = 30; PlotInfo.FigureTitle = ['Predict no-activation, binning x',num2str(BinningFactor)]; 
PlotMatrices_AlignByActivityOrPrediction(StructureOut.DowngradientAboveThresholdAtLeast5sec_Binning30, MinFramesInBoundFromEventStart/BinningFactor, PlotInfo);

%% MEASURED ACTIVITY:  Align data by measured activity (or no-activity) events 
%%% Generate matrices and save processed data in structures  %%%%
StructureOut = CalculateMatricesAlignedByActivity(Data, FramesWindow);

%%%% Figures 4A, 4B and S4F %%%%
PlotInfo.PlotActivity          = false;
PlotInfo.PlotBehavior          = true;
MinFramesInBoundFromEventStart = 5*FrameRate; % Do not analyze behavior of animals that are out of bounds

% Measured activity
PlotInfo.FigureTitle = 'Measured activity';    
[CDF_Fwd, CDF_Rev, SortedTimes_Fwd, SortedTimes_Rev] = PlotMatrices_AlignByActivityOrPrediction(StructureOut.ActivationInitiation, MinFramesInBoundFromEventStart, PlotInfo);
CDFs.Activated.Fwd    = CDF_Fwd;
CDFs.Activated.Rev    = CDF_Rev;
Repeats.Activated.Fwd = length(SortedTimes_Fwd);
Repeats.Activated.Rev = length(SortedTimes_Rev);
Times.Activated.Fwd   = SortedTimes_Fwd;
Times.Activated.Rev   = SortedTimes_Rev;

% Measured no-activity
PlotInfo.FigureTitle = 'Measured no-activity';    
[CDF_Fwd, CDF_Rev, SortedTimes_Fwd, SortedTimes_Rev] = PlotMatrices_AlignByActivityOrPrediction(StructureOut.DowngradientNoActivation, MinFramesInBoundFromEventStart, PlotInfo);
CDFs.NotActivated.Fwd    = CDF_Fwd;
CDFs.NotActivated.Rev    = CDF_Rev;
Repeats.NotActivated.Fwd = length(SortedTimes_Fwd);
Repeats.NotActivated.Rev = length(SortedTimes_Rev);
Times.NotActivated.Fwd   = SortedTimes_Fwd;
Times.NotActivated.Rev   = SortedTimes_Rev;

%%% Control #2: see that analysis of no-activity are not sensitive to the frame rate
PlotInfo.PlotActivity = false;
PlotInfo.PlotBehavior = true;
% Use 10x binning 
BinningFactor         = 10;   
MinFramesInBoundFromEventStart = 5*FrameRate/BinningFactor; % Do not analyze behavior of animals that are out of bounds
PlotInfo.FigureTitle = ['Measured no-activation, binning x',num2str(BinningFactor)]; 
PlotMatrices_AlignByActivityOrPrediction(StructureOut.DowngradientNoActivation_Binning10, MinFramesInBoundFromEventStart, PlotInfo);

% Use 30x binning
BinningFactor         = 30;   
MinFramesInBoundFromEventStart = 5*FrameRate/BinningFactor; % Do not analyze behavior of animals that are out of bounds
PlotInfo.FigureTitle = ['Measured no-activation, binning x',num2str(BinningFactor)]; 
PlotMatrices_AlignByActivityOrPrediction(StructureOut.DowngradientNoActivation_Binning30, MinFramesInBoundFromEventStart, PlotInfo);

%% CDF figures 
CDF_Time = (1:length(CDFs.Activated.Fwd))/FrameRate;
figure('name','CDF for initiating aversive behavior','position',[680 558 201 420]);  % Figure 4B
subplot(2,1,1); 
    plot(CDF_Time,CDFs.Activated.Fwd,'k-'); hold on; 
    plot(CDF_Time,CDFs.NotActivated.Fwd,'k--'); 
    plot(CDF_Time,CDFs.Baseline.Fwd,'b-'); 
    axis([0 5 0 1]); set(gca,'box','off','xtick',0:5,'ytick',0:0.5:1); ylabel('CDF') 
subplot(2,1,2); 
    plot(CDF_Time,CDFs.Predicted.Fwd,'k-'); hold on; 
    plot(CDF_Time,CDFs.NotPredicted.Fwd,'k--'); 
    plot(CDF_Time,CDFs.Baseline.Fwd,'b-'); 
    axis([0 5 0 1]); set(gca,'box','off','xtick',0:5,'ytick',0:0.5:1); ylabel('CDF') 
    xlabel('Time [sec]')
figure('name','CDF for initiating forward behavior','position',[680 558 201 420]); % Figure S4D
subplot(2,1,1); 
    plot(CDF_Time,CDFs.Activated.Rev,'k-'); hold on; 
    plot(CDF_Time,CDFs.NotActivated.Rev,'k--'); 
    plot(CDF_Time,CDFs.Baseline.Rev,'b-'); 
    axis([0 5 0 1]); set(gca,'box','off','xtick',0:5,'ytick',0:0.5:1); ylabel('CDF')   
subplot(2,1,2); 
    plot(CDF_Time,CDFs.Predicted.Rev,'k-'); hold on; 
    plot(CDF_Time,CDFs.NotPredicted.Rev,'k--'); 
    plot(CDF_Time,CDFs.Baseline.Rev,'b-'); 
    axis([0 5 0 1]); set(gca,'box','off','xtick',0:5,'ytick',0:0.5:1); ylabel('CDF')  
    xlabel('Time [sec]')

%% Two-sample Kolmogorov-Smirnov test
[Pvalue.Predicted_Fwd,    KSstatistic.Predicted_Fwd]   = kstest2_inline(CDFs.Predicted.Fwd,    CDFs.Baseline.Fwd, Repeats.Predicted.Fwd,    Repeats.Baseline.Fwd, 0);
[Pvalue.NotPredicted_Fwd, KSstatistic.NotPredicted_Fwd]= kstest2_inline(CDFs.NotPredicted.Fwd, CDFs.Baseline.Fwd, Repeats.NotPredicted.Fwd, Repeats.Baseline.Fwd, 0);
[Pvalue.Predicted_Rev,    KSstatistic.Predicted_Rev]   = kstest2_inline(CDFs.Predicted.Rev,    CDFs.Baseline.Rev, Repeats.Predicted.Rev,    Repeats.Baseline.Rev, 0);
[Pvalue.NotPredicted_Rev, KSstatistic.NotPredicted_Rev]= kstest2_inline(CDFs.NotPredicted.Rev, CDFs.Baseline.Rev, Repeats.NotPredicted.Rev, Repeats.Baseline.Rev, 0);
[Pvalue.Activated_Fwd,    KSstatistic.Activated_Fwd]   = kstest2_inline(CDFs.Activated.Fwd,    CDFs.Baseline.Fwd, Repeats.Activated.Fwd,    Repeats.Baseline.Fwd, 0);
[Pvalue.NotActivated_Fwd, KSstatistic.NotActivated_Fwd]= kstest2_inline(CDFs.NotActivated.Fwd, CDFs.Baseline.Fwd, Repeats.NotActivated.Fwd, Repeats.Baseline.Fwd, 0);
[Pvalue.Activated_Rev,    KSstatistic.Activated_Rev]   = kstest2_inline(CDFs.Activated.Rev,    CDFs.Baseline.Rev, Repeats.Activated.Rev,    Repeats.Baseline.Rev, 0);
[Pvalue.NotActivated_Rev, KSstatistic.NotActivated_Rev]= kstest2_inline(CDFs.NotActivated.Rev, CDFs.Baseline.Rev, Repeats.NotActivated.Rev, Repeats.Baseline.Rev, 0);

EffectSize.MaxOfControls = max([KSstatistic.NotPredicted_Fwd, KSstatistic.NotPredicted_Rev, KSstatistic.NotActivated_Fwd, KSstatistic.NotActivated_Rev]);
EffectSize.MinOfForwards = min([KSstatistic.Predicted_Fwd, KSstatistic.Activated_Fwd]);
EffectSize.Predicted_Rev = KSstatistic.Predicted_Rev;
EffectSize.Activated_Rev = KSstatistic.Activated_Rev;

% save('CDFs_Figure4.mat','CDFs','Repeats','Pvalue','KSstatistic','EffectSize','Times','-v7.3') 
% load('CDFs_Figure4.mat','CDFs','Repeats','Pvalue','KSstatistic','EffectSize','Times') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Matrices aligned by behavioral events (Figures 4C and S4H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
PlotInfo.XLIMITS               = [-5 5];
MinFramesInBoundFromEventStart = 0;   

%% Align Matrices by behavior state
% %%% Generate matrices and save processed data in structures  %%%%
StructureOut = CalculateMatricesAlignedByBehaviorState(Data, FramesWindow);
% Assign F0 to be the activity at t=0 (instead of average before the event)
StructureOut.RevToFwd.EventPastNeuronValue           = StructureOut.RevToFwd.NeuronValue(:,ZeroIndex);
StructureOut.FwdToRev.EventPastNeuronValue           = StructureOut.FwdToRev.NeuronValue(:,ZeroIndex);
StructureOut.ReversalInitiation.EventPastNeuronValue = StructureOut.ReversalInitiation.NeuronValue(:,ZeroIndex);

PlotInfo.FigureTitle = 'Aligned by reversal initiation'; 
PlotMatrices_AlignedByBehavior(StructureOut.FwdToRev, MinFramesInBoundFromEventStart, PlotInfo);
% PlotMatrices_AlignedByBehavior(StructureOut.RevToFwd, MinFramesInBoundFromEventStart, PlotInfo);
% PlotMatrices_AlignedByBehavior(StructureOut.ReversalInitiation, MinFramesInBoundFromEventStart, PlotInfo);

%% Align Matrices by downgradient navigation events
StructureOut = CalculateMatricesAlignedByNavigation(Data, FramesWindow);
% Assign F0 to be the activity at t=0 (instead of average before the event)
StructureOut.DowngradientEnd_FastChange.EventPastNeuronValue   = StructureOut.DowngradientEnd_FastChange.NeuronValue(:,ZeroIndex);
StructureOut.DowngradientEnd_NotOOB.EventPastNeuronValue       = StructureOut.DowngradientEnd_NotOOB.NeuronValue(:,ZeroIndex);
StructureOut.DowngradientFramesNotNearEnd.EventPastNeuronValue = StructureOut.DowngradientFramesNotNearEnd.NeuronValue(:,ZeroIndex);

PlotInfo.FigureTitle = 'Aligned by downgradient termination'; 
PlotMatrices_AlignedByBehavior(StructureOut.DowngradientEnd_NotOOB, MinFramesInBoundFromEventStart, PlotInfo);
% PlotMatrices_AlignedByBehavior(StructureOut.DowngradientEnd_FastChange, MinFramesInBoundFromEventStart, PlotInfo);
% PlotMatrices_AlignedByBehavior(StructureOut.DowngradientFramesNotNearEnd, MinFramesInBoundFromEventStart, PlotInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Traces Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
PlotStimulusInsteadOfConcentration = false;
% tr = 16; Plot_Environment_Activity_Behavior(Data,tr, [ 1150 1165], 1.2, PlotStimulusInsteadOfConcentration,[1 5.5]*1e-7, true); % Forward example
% set(gcf,'position',[825   289   117   652]);
tr = 18; Plot_Environment_Activity_Behavior(Data,tr, [ 1364 1379], 1.8, PlotStimulusInsteadOfConcentration,[1 5.5]*1e-7, true);   % Figure 3C, Forward example
set(gcf,'position',[825   289   117   652]);
% tr = 19; Plot_Environment_Activity_Behavior(Data,tr, [172 187], 4, PlotStimulusInsteadOfConcentration,[0 50]*1e-7, true);       % Forward example
% set(gcf,'position',[825   289   117   652]);
tr = 12; Plot_Environment_Activity_Behavior(Data,tr, [18 33], 6, PlotStimulusInsteadOfConcentration,[1 5]*1e-7, true);            % Figure S4B, Reverse example
set(gcf,'position',[825   289   117   652]);

return

%% Data sorting by events
function StructureOut = CalculateMatricesAlignedByThresholdCrossingEvents(Data, FramesWindow)
MinimumFramesForThresholdPassing    = 0.2*30; % 0.2 second
MinimumFramesBeforeThresholdPassing = 5*30  ; % 5 seconds

%% Initialization
MatrixSize                   = size(Data.NaN);
NumberOfFrames               = MatrixSize(2);
NumberOfWorms                = MatrixSize(1);
C_Minus_T                    = Data.C - Data.AdaptiveThreshold; 
Sign_C_Minus_T               = sign(C_Minus_T);
BelowThreshold               = Sign_C_Minus_T == -1;
AboveThreshold               = Sign_C_Minus_T == 1;
OOB                          = Data.BehaviorCode_LowLevel==0;                        
CrossingThresholdFrames      = false(size(BelowThreshold));
CrossingAboveThresholdFrames = false(size(BelowThreshold));
for tr=1:NumberOfWorms
    if any(BelowThreshold(tr,:))
        [~, LongThresholdPassing_StartFrames, ~] = FindSegment(BelowThreshold(tr,:), MinimumFramesForThresholdPassing);   
        for frame = LongThresholdPassing_StartFrames'
            PreviousFrames = (frame-MinimumFramesBeforeThresholdPassing):(frame-1);
            PreviousFrames = PreviousFrames(PreviousFrames>=1);
            PreviousFrames = PreviousFrames(PreviousFrames<=NumberOfFrames);
            if all(AboveThreshold(tr,PreviousFrames)) && ~OOB(tr,frame)
                CrossingThresholdFrames(tr,frame) = true;
            end
        end
    end
    if any(AboveThreshold(tr,:))
        [~, LongThresholdPassing_StartFrames, ~] = FindSegment(AboveThreshold(tr,:), MinimumFramesForThresholdPassing);    
       for frame = LongThresholdPassing_StartFrames'
            PreviousFrames = (frame-MinimumFramesBeforeThresholdPassing):(frame-1);
            PreviousFrames = PreviousFrames(PreviousFrames>=1);
            PreviousFrames = PreviousFrames(PreviousFrames<=NumberOfFrames);
            if all(BelowThreshold(tr,PreviousFrames)) && ~OOB(tr,frame)
                CrossingAboveThresholdFrames(tr,frame) = true;
            end
       end
    end
end

%% Remove frames before and after gradient is established, and remove control non-gradient experiments 
GradientFramesMatrix = Data.GradientFramesMatrix;
CrossingThresholdFrames_DuringGradient                             = CrossingThresholdFrames;
CrossingThresholdFrames_DuringGradient(~GradientFramesMatrix)      = false;
CrossingAboveThresholdFrames_DuringGradient                        = CrossingAboveThresholdFrames;
CrossingAboveThresholdFrames_DuringGradient(~GradientFramesMatrix) = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% %%%%%%%  Aligning matrices by crossing BELOW the threshold
StructureOut_CrossingThreshold      = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, CrossingThresholdFrames_DuringGradient, FramesWindow);

%% %%%%%%%  Aligning matrices by crossing ABOVE the threshold
StructureOut_CrossingAboveThreshold = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, CrossingAboveThresholdFrames_DuringGradient, FramesWindow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% %%%%%%%  CONTROL: All downgradient frames for which C > T 
DownGradient                         = Data.Downgradient;   
DownGradient(~AboveThreshold)        = false;
DownGradient(~GradientFramesMatrix)  = false;
DownGradient(:,[1:301, 62699:63000]) = false; 
DownGradientAboveThreshold           = DownGradient;
StructureOut_DowngradientAboveThreshold = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, DownGradient, FramesWindow);

% eliminate frames that are up to 5 seconds before predicted activation 
DownGradient = DownGradientAboveThreshold;
Window = (-5*30):0; 
for row=1:size(DownGradient,1)
    CurrentCrossingThresholdFrames = find(CrossingThresholdFrames(row,:));
    for frame=CurrentCrossingThresholdFrames
        CurrentFrameWindow = Window + frame;
        CurrentFrameWindow = CurrentFrameWindow(CurrentFrameWindow>=1);
        DownGradient(row,CurrentFrameWindow) = false;
    end
end
StructureOut_DowngradientAboveThresholdForAtLeast5sec = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, DownGradient, FramesWindow);

%% %%%%%%%  CONTROL: All frames for which C > T regardless of their direction
AllFrames      = true(size(Data.C));   
OOB            = Data.BehaviorCode_LowLevel==0 ; 
AllFrames(OOB) = false;

AllFrames(~GradientFramesMatrix)  = false;
AllFrames(:,[1:301, 62699:63000]) = false; 
AllFrames_RegardlessOfThreshold   = AllFrames;
StructureOut_AllFrames_RegardlessOfThreshold = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, AllFrames_RegardlessOfThreshold, FramesWindow);

AllFrames(~AboveThreshold)        = false;
AllFramesOriginal                 = AllFrames;
StructureOut_AllFramesAboveThreshold = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, AllFrames, FramesWindow);

% eliminate frames that are up to 5 seconds before activation 
AllFrames = AllFramesOriginal;
Window = (-5*30):0; % eliminate frames that are up to 5 seconds before activation 
for row=1:size(AllFrames,1)
    CurrentCrossingThresholdFrames = find(CrossingThresholdFrames(row,:));
    for frame=CurrentCrossingThresholdFrames
        CurrentFrameWindow = Window + frame;
        CurrentFrameWindow = CurrentFrameWindow(CurrentFrameWindow>=1);
        AllFrames(row,CurrentFrameWindow) = false;
    end
end
StructureOut_AllFramesAboveThresholdForAtLeast5sec = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, AllFrames, FramesWindow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% CONTROL: validate that our analysis and conclusions are not sensitive to the movie frame rate
%  In Figures 3F and 4A we align matrices by "no activation" or "no predicted activation"  
%  Intuition: 
%    There are many frames like that, and it scales with the movie frame rate since an initiation event is a single frame
%    However, the matrices allow overlapping frames and is also scaling with the frame rate.   
%       For example, if frame 'n' and frame 'n+1' are both "no activation" frames then the matrix will include
%       a row with data from [n,n+FW] and an additional overlapping row with data from [n+1,n+1+FW], where FW is the matrix FramesWindow   
%    When using the matrices to calculate probabilities or averages, both denominator and numerator scales    
%    with frame rate and therefore this parameter cancels out.   
%%%%% All downgradient frames for which C > T  
%%%%% binning x10 %%%%%
Binning_NumOfFrames  = 10;
FramesWindowBinning = (FramesWindow(1)/Binning_NumOfFrames):(FramesWindow(end)/Binning_NumOfFrames);
[Data_Binned, ~, CrossingThresholdFrames_BinnedOrGate]    = BinData(Data,CrossingThresholdFrames,Binning_NumOfFrames);  % Outputs: [Data_Binned, AndGateBinned, OrGateBinned]
 
% Above threshold, in gradient
DownGradient_Binned    = Data_Binned.Downgradient;   
AboveThreshold_Binned  = Data_Binned.C > Data_Binned.AdaptiveThreshold;
DownGradient_Binned(~AboveThreshold_Binned)            = false;
DownGradient_Binned(~Data_Binned.GradientFramesMatrix) = false;
DownGradient_Binned(:,[1:round(301/Binning_NumOfFrames), (end-round(301/Binning_NumOfFrames)):end]) = false; 
DownGradientAboveThreshold_Binned = DownGradient_Binned;

% Eliminate frames that are up to 5 seconds before predicted activation 
DownGradient = DownGradientAboveThreshold_Binned;
Window = (-5*30/Binning_NumOfFrames):0; 
for row=1:size(DownGradient,1)
    CurrentCrossingThresholdFrames = find(CrossingThresholdFrames_BinnedOrGate(row,:));
    for frame=CurrentCrossingThresholdFrames
        CurrentFrameWindow = Window + frame;
        CurrentFrameWindow = CurrentFrameWindow(CurrentFrameWindow>=1);
        DownGradient(row,CurrentFrameWindow) = false;
    end
end
DownGradient_PredictedNoActivation_Binned = DownGradient;
% Align
StructureOut_DowngradientAboveThresholdForAtLeast5sec_Binning10      = Align_Concentration_Activity_Behavior_By_EventMatrix(Data_Binned, DownGradient_PredictedNoActivation_Binned, FramesWindowBinning);
StructureOut_DowngradientAboveThresholdForAtLeast5sec_Binning10.Time = StructureOut_DowngradientAboveThresholdForAtLeast5sec_Binning10.Time*Binning_NumOfFrames; % Correction from the function above

%%%%% binning x30 %%%%%
Binning_NumOfFrames  = 30;
FramesWindowBinning = (FramesWindow(1)/Binning_NumOfFrames):(FramesWindow(end)/Binning_NumOfFrames);
[Data_Binned, ~, CrossingThresholdFrames_BinnedOrGate]    = BinData(Data,CrossingThresholdFrames,Binning_NumOfFrames);  % Outputs: [Data_Binned, AndGateBinned, OrGateBinned]
 
% Above threshold, in gradient
DownGradient_Binned    = Data_Binned.Downgradient;   
AboveThreshold_Binned  = Data_Binned.C > Data_Binned.AdaptiveThreshold;
DownGradient_Binned(~AboveThreshold_Binned)            = false;
DownGradient_Binned(~Data_Binned.GradientFramesMatrix) = false;
DownGradient_Binned(:,[1:round(301/Binning_NumOfFrames), (end-round(301/Binning_NumOfFrames)):end]) = false; 
DownGradientAboveThreshold_Binned = DownGradient_Binned;

% Eliminate frames that are up to 5 seconds before predicted activation 
DownGradient = DownGradientAboveThreshold_Binned;
Window = (-5*30/Binning_NumOfFrames):0; 
for row=1:size(DownGradient,1)
    CurrentCrossingThresholdFrames = find(CrossingThresholdFrames_BinnedOrGate(row,:));
    for frame=CurrentCrossingThresholdFrames
        CurrentFrameWindow = Window + frame;
        CurrentFrameWindow = CurrentFrameWindow(CurrentFrameWindow>=1);
        DownGradient(row,CurrentFrameWindow) = false;
    end
end
DownGradient_PredictedNoActivation_Binned = DownGradient;
% Align
StructureOut_DowngradientAboveThresholdForAtLeast5sec_Binning30      = Align_Concentration_Activity_Behavior_By_EventMatrix(Data_Binned, DownGradient_PredictedNoActivation_Binned, FramesWindowBinning);
StructureOut_DowngradientAboveThresholdForAtLeast5sec_Binning30.Time = StructureOut_DowngradientAboveThresholdForAtLeast5sec_Binning30.Time*Binning_NumOfFrames; % Correction from the function above

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Assign to structure
StructureOut.CrossingThreshold                      = StructureOut_CrossingThreshold;
StructureOut.CrossingAboveThreshold                 = StructureOut_CrossingAboveThreshold;
StructureOut.DowngradientAboveThreshold             = StructureOut_DowngradientAboveThreshold;
StructureOut.DowngradientAboveThresholdAtLeast5sec  = StructureOut_DowngradientAboveThresholdForAtLeast5sec;
StructureOut.AllFrames_RegardlessOfThreshold        = StructureOut_AllFrames_RegardlessOfThreshold;
StructureOut.AllFramesAboveThreshold                = StructureOut_AllFramesAboveThreshold;
StructureOut.AllFramesAboveThresholdAtLeast5sec     = StructureOut_AllFramesAboveThresholdForAtLeast5sec;

StructureOut.DowngradientAboveThresholdAtLeast5sec_Binning10 = StructureOut_DowngradientAboveThresholdForAtLeast5sec_Binning10;
StructureOut.DowngradientAboveThresholdAtLeast5sec_Binning30 = StructureOut_DowngradientAboveThresholdForAtLeast5sec_Binning30;

return

function StructureOut = CalculateMatricesAlignedByActivity(Data, FramesWindow)
%% Initialization
ActivationInitiationFrames  = Data.ActivationStartFrames;  
ActivationPeakFrames        = Data.ActivationPeakFrames;  
ActivationOOB               = Data.ActivationFrameOOB;  

ActivationInitiation_NotOOB                = ActivationInitiationFrames; 
ActivationInitiation_NotOOB(ActivationOOB) = false;
ActivationPeak_NotOOB                      = ActivationPeakFrames; 
ActivationPeak_NotOOB(ActivationOOB)       = false;

%% Remove frames before and after gradient is established, and remove control non-gradient experiments 
GradientFramesMatrix = Data.GradientFramesMatrix;

ActivationInitiation_NotOOB_DuringGradient                        = ActivationInitiation_NotOOB;
ActivationInitiation_NotOOB_DuringGradient(~GradientFramesMatrix) = false;
ActivationPeak_NotOOB_DuringGradient                              = ActivationPeak_NotOOB;
ActivationPeak_NotOOB_DuringGradient(~GradientFramesMatrix)       = false;

%% %%%%%%%  CONTROL: Aligning matrices by NO-Initiation of activity
DownGradient_NoActivation_StartLow = Data.Downgradient_NoActivation & Data.Downgradient_ActivityStartsLow;   
DownGradient_NoActivation_StartLow(~GradientFramesMatrix)       = false;
DownGradient_NoActivation_StartLow(:,[1:301, 62699:63000])      = false; 
Window = (-5*30):0; % eliminate frames that are up to 5 seconds before activation 
for row=1:size(DownGradient_NoActivation_StartLow,1)
    CurrentActivationFrames = find(ActivationInitiationFrames(row,:));
    for frame=CurrentActivationFrames
        CurrentFrameWindow = Window + frame;
        DownGradient_NoActivation_StartLow(row,CurrentFrameWindow) = false;
    end
end
DownGradient_NoActivation_StartLow(Data.deltaFOverF_Interpolated(:)>1) = false;

StructureOut_DowngradientNoActivation = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, DownGradient_NoActivation_StartLow, FramesWindow);

%% %%%%%%%  Aligning matrices by peak of activity
StructureOut_ActivationPeak           = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, ActivationPeak_NotOOB_DuringGradient, FramesWindow);

%% %%%%%%%  Aligning matrices by initiation of activity
StructureOut_ActivationInitiation     = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, ActivationInitiation_NotOOB_DuringGradient, FramesWindow);


%% CONTROL: validate that our analysis and conclusions are not sensitive to the movie frame rate
%  In Figures 3F and 4A we align matrices by "no activation" or "no predicted activation"  
%  Intuition: 
%    There are many frames like that, and it scales with the movie frame rate since an initiation event is a single frame
%    However, the matrices allow overlapping frames and is also scaling with the frame rate.   
%       For example, if frame 'n' and frame 'n+1' are both "no activation" frames then the matrix will include
%       a row with data from [n,n+FW] and an additional overlapping row with data from [n+1,n+1+FW], where FW is the matrix FramesWindow   
%    When using the matrices to calculate probabilities or averages, both denominator and numerator scales    
%    with frame rate and therefore this parameter cancels out.  
%%%%% Binning x10 %%%%%%
Binning_NumOfFrames  = 10;
FramesWindowBinning = (FramesWindow(1)/Binning_NumOfFrames):(FramesWindow(end)/Binning_NumOfFrames);
[Data_Binned, DownGradient_NoActivation_StartLow_Binned] = BinData(Data, DownGradient_NoActivation_StartLow, Binning_NumOfFrames);  % Outputs: [Data_Binned, AndGateBinned]
StructureOut_DowngradientNoActivation_Binning10      = Align_Concentration_Activity_Behavior_By_EventMatrix(Data_Binned, DownGradient_NoActivation_StartLow_Binned, FramesWindowBinning);
StructureOut_DowngradientNoActivation_Binning10.Time = StructureOut_DowngradientNoActivation_Binning10.Time*Binning_NumOfFrames; % Correction from the function above

%%%%% Binning x30 %%%%%%
Binning_NumOfFrames  = 30;
FramesWindowBinning = (FramesWindow(1)/Binning_NumOfFrames):(FramesWindow(end)/Binning_NumOfFrames);
[Data_Binned, DownGradient_NoActivation_StartLow_Binned] = BinData(Data, DownGradient_NoActivation_StartLow, Binning_NumOfFrames);  % Outputs: [Data_Binned, AndGateBinned]
StructureOut_DowngradientNoActivation_Binning30      = Align_Concentration_Activity_Behavior_By_EventMatrix(Data_Binned, DownGradient_NoActivation_StartLow_Binned, FramesWindowBinning);
StructureOut_DowngradientNoActivation_Binning30.Time = StructureOut_DowngradientNoActivation_Binning30.Time*Binning_NumOfFrames; % Correction from the function above

%% Assign to structure
StructureOut.ActivationInitiation             = StructureOut_ActivationInitiation;
StructureOut.ActivationPeak                   = StructureOut_ActivationPeak;
StructureOut.DowngradientNoActivation         = StructureOut_DowngradientNoActivation;
StructureOut.DowngradientNoActivation_Binning10 = StructureOut_DowngradientNoActivation_Binning10;
StructureOut.DowngradientNoActivation_Binning30 = StructureOut_DowngradientNoActivation_Binning30;

return

function StructureOut = CalculateMatricesAlignedByBehaviorState(Data, FramesWindow)
%% Initialization
SwitchingRevToFwd        = Data.SwitchingRevToFwd;  
SwitchingFwdToRev        = Data.SwitchingFwdToRev;  
ReversalInitiationNotOOB = Data.ReversalInitiationNotOOB;  

%% Remove frames before and after gradient is established, and remove control non-gradient experiments 
GradientFramesMatrix = Data.GradientFramesMatrix;
SwitchingRevToFwd_DuringGradient                                = SwitchingRevToFwd;
SwitchingRevToFwd_DuringGradient(~GradientFramesMatrix)         = false;
SwitchingRevToFwd_DuringGradient(:,[1:301, 62699:63000])        = false;
SwitchingFwdToRev_DuringGradient                                = SwitchingFwdToRev;
SwitchingFwdToRev_DuringGradient(~GradientFramesMatrix)         = false;
SwitchingFwdToRev_DuringGradient(:,[1:301, 62699:63000])        = false;
ReversalInitiationNotOOB_DuringGradient                         = ReversalInitiationNotOOB;
ReversalInitiationNotOOB_DuringGradient(~GradientFramesMatrix)  = false;
ReversalInitiationNotOOB_DuringGradient(:,[1:301, 62699:63000]) = false;

%% Aligning matrices 
StructureOut_FwdToRev           = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, SwitchingFwdToRev_DuringGradient, FramesWindow);
StructureOut_RevToFwd           = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, SwitchingRevToFwd_DuringGradient, FramesWindow);
StructureOut_ReversalInitiation = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, ReversalInitiationNotOOB_DuringGradient, FramesWindow);

%% Assign to structure
StructureOut.RevToFwd           = StructureOut_RevToFwd;
StructureOut.FwdToRev           = StructureOut_FwdToRev;
StructureOut.ReversalInitiation = StructureOut_ReversalInitiation;

return

function StructureOut = CalculateMatricesAlignedByNavigation(Data, FramesWindow)
%% Initialization
MatrixSize                  = size(Data.NaN);
NumberOfFrames              = MatrixSize(2);
DowngradientNavigationEndFrames_WithoutOOB          = Data.DowngradientNavigationEndFrames_WithoutOOB;  
DowngradientNavigationEndFrames_FastDirectionChange = Data.DowngradientNavigationEndFrames_FastDirectionChange; 

%% Remove frames before and after gradient is established, and remove control non-gradient experiments 
GradientFramesMatrix = Data.GradientFramesMatrix;

DowngradientEnd_NotOOB_DuringGradient                            = DowngradientNavigationEndFrames_WithoutOOB;
DowngradientEnd_NotOOB_DuringGradient(~GradientFramesMatrix)     = false;
DowngradientEnd_NotOOB_DuringGradient(:,[1:301, 62699:63000])    = false;
DowngradientEnd_FastChange_DuringGradient                        = DowngradientNavigationEndFrames_FastDirectionChange;
DowngradientEnd_FastChange_DuringGradient(~GradientFramesMatrix) = false;
DowngradientEnd_FastChange_DuringGradient(:,[1:301, 62699:63000])= false;

%% %%%%%%%  CONTROL: Aligning matrices by downgradient behavior, NOT near downgradient ending frame 
DowngradientEnd_All_DuringGradient                                 = Data.DowngradientNavigationEndFrames_All;
DowngradientFramesNotNearEnd_DuringGradient                        = Data.DowngradientNavigation; 
DowngradientFramesNotNearEnd_DuringGradient(~GradientFramesMatrix) = false;
DowngradientFramesNotNearEnd_DuringGradient(:,[1:301, 62699:63000])= false; 

Window = (-5*30):(5*30); % eliminate frames that are up to 5 seconds before or after navigation change 
for row=1:size(DowngradientFramesNotNearEnd_DuringGradient,1)
    CurrentDowngradientEndFrames = find(DowngradientEnd_All_DuringGradient(row,:));
    for frame=CurrentDowngradientEndFrames
        CurrentFrameWindow = Window + frame;
        CurrentFrameWindow = CurrentFrameWindow((CurrentFrameWindow>0)&(CurrentFrameWindow<NumberOfFrames));
        DowngradientFramesNotNearEnd_DuringGradient(row,CurrentFrameWindow) = false;
    end
end

StructureOut_DowngradientFramesNotNearEnd = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, DowngradientFramesNotNearEnd_DuringGradient, FramesWindow);

%% %%%%%%%  Aligning matrices by ending of downgradient behavior: only events resulting in fast direction change
StructureOut_DowngradientEnd_FastChange   = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, DowngradientEnd_FastChange_DuringGradient, FramesWindow);

%% %%%%%%%  Aligning matrices by ending of downgradient behavior: all events that are not OOB
StructureOut_DowngradientEnd_NotOOB = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, DowngradientEnd_NotOOB_DuringGradient, FramesWindow);

%% Asign to structure
StructureOut.DowngradientEnd_NotOOB       = StructureOut_DowngradientEnd_NotOOB;
StructureOut.DowngradientEnd_FastChange   = StructureOut_DowngradientEnd_FastChange;
StructureOut.DowngradientFramesNotNearEnd = StructureOut_DowngradientFramesNotNearEnd;

return

function  [VecTrueIfFrameWithinSegment, StartFrames, EndFrames]= FindSegment(LogicalVec, MinimumSegmentLength)
% Example:
% [FramesWithLongNaNSegments, StartNaNSegments, EndNaNSegments]= FindSegment(NaN_Frames, MaxAllowedMissingFrames);

StartFrames = (find(diff(LogicalVec)==1)+1)'; if LogicalVec(1)==1,   StartFrames = [1; StartFrames]; end
EndFrames   = (find(diff(LogicalVec)==-1))';  if LogicalVec(end)==1, EndFrames   = [EndFrames; length(LogicalVec)]; end
if length(StartFrames) ~= length(EndFrames)
    disp('error in segments calculation. keyboard...');
    keyboard;
end
LongSegments_indices = find((EndFrames - StartFrames+1)>MinimumSegmentLength);
if ~isempty(LongSegments_indices)
    VecTrueIfFrameWithinSegment = false(1,length(LogicalVec));
    for ind = LongSegments_indices'
        VecTrueIfFrameWithinSegment(StartFrames(ind):EndFrames(ind)) = true;    
    end
end        
% Correct to keep only long segments
LogicalVec = VecTrueIfFrameWithinSegment;
StartFrames = (find(diff(LogicalVec)==1)+1)'; if LogicalVec(1)==1,   StartFrames = [1; StartFrames]; end
EndFrames   = (find(diff(LogicalVec)==-1))';  if LogicalVec(end)==1, EndFrames   = [EndFrames; length(LogicalVec)]; end
if length(StartFrames) ~= length(EndFrames)
    disp('error in segments calculation. keyboard...');
    keyboard;
end

return

%% Data Binning for control
function [Data_Binned, LogicalVector_AndGateBinned, LogicalVector_OrGateBinned] = BinData(Data, LogicalVector, Binning_NumOfFrames)

% This function bin SOME of the matrices in the structure Data 
% While binning is trivial for most numeric variables, it need to be done carefully for logical vectors.    
% Here we will bin logical vectors be two types possible binning strategies: 
%   OR  gate: bin 'n' is true if AT LEAST ONE frame in the interval [(n-1)*BinSize, n*BinSize] is true 
%   AND gate: bin 'n' is true if     ALL     frames in the interval [(n-1)*BinSize, n*BinSize] are true 
%   
%% Initialization
TotalNumberOfFrames   = size(Data.Coordinates_X_Smoothed,2);
% BinFrameVec           = 1:Binning_NumOfFrames:TotalNumberOfFrames;
BinFrameVec           = (Binning_NumOfFrames:Binning_NumOfFrames:TotalNumberOfFrames)-round(Binning_NumOfFrames/2);
LogicalFields_OrGate  = {'NaN','OutOfBounds_WithCollisions'}; 
LogicalFields_ANDGate = {'FlagTrueIfMidlineIsReliable','Downgradient','Downgradient_NoActivation','DowngradientNavigation','UpgradientNavigation','Dwell','GradientFramesMatrix'}; 
Fields_Behavior       = {'BehaviorCode_HighLevel','BehaviorCode_LowLevel'}; % Demand same behavior within the bin, else give an out-of-bounds value  
  
% Filter
a_MA = 1; 
b_MA = ones(1,Binning_NumOfFrames)/Binning_NumOfFrames;
D_MA = round(mean(grpdelay(b_MA,a_MA)));  % Filter delay, for float number vectors

%% Bin Data
FieldNames = fieldnames(Data);
for f_ind = 1:length(FieldNames)
    CurrentField = FieldNames{f_ind};
    X = Data.(CurrentField);
    
    if isfloat(X)
        X_filter = filter(b_MA,a_MA,X')'; 
        X_filter = X_filter(:,D_MA+1:end);
        X_filter(:,(TotalNumberOfFrames-D_MA):TotalNumberOfFrames) = NaN;
        X_Bin    = X_filter(:,BinFrameVec);        
        
    elseif ~isempty(find(strcmpi(CurrentField,Fields_Behavior),1)) % Behavior field
        X_Bin    = X(:,BinFrameVec);
%         %% Correct behavior to be out-of-bounds (OOB) if beahvior is not the same code throughout the bin 
%         %  The problem: it generates a bias- all non-OOB behaviors tend to terminate with an OOB state 
%         NoChangeInBehavior             = [NaN*ones(size(X,1),1), diff(double(X),[],2)]==0; 
%         NoChangeInBehavior_filter      = filter(b_MA,a_MA,double(NoChangeInBehavior)')'; 
%         NoChangeInBehavior_ANDgate     = NoChangeInBehavior_filter > (1-eps);
%         NoChangeInBehavior_Bin         = NoChangeInBehavior_ANDgate(:,BinFrameVec);             
%         X_Bin(~NoChangeInBehavior_Bin) = 0; 
        
%     elseif islogical(X)
    elseif ~isempty(find(strcmpi(CurrentField,LogicalFields_OrGate),1)) % logical OR gate field
        X_filter = filter(b_MA,a_MA,double(X)')'; 
        X_ORgate = X_filter > 0;
        X_Bin    = X_ORgate(:,BinFrameVec);
        
%         figure; plot(1:size(X,2),X(1,:),'k*'); hold on; plot(1:size(X,2),X_filter(1,:),'ro'); ylim([-0.1 1.1])  
%         figure; plot(1:size(X,2),X(1,:),'k*'); hold on; plot(1:size(X,2),X_ORgate(1,:),'ro'); ylim([-0.1 1.1])  
%         figure; plot(1:size(X,2),X(1,:),'k*'); hold on; plot(BinFrameVec,X_Bin(1,:),'ro'); ylim([-0.1 1.1])    
    elseif ~isempty(find(strcmpi(CurrentField,LogicalFields_ANDGate),1)) % logical AND gate field
        X_filter  = filter(b_MA,a_MA,double(X)')'; 
        X_ANDgate = X_filter > (1-eps);
        X_Bin     = X_ANDgate(:,BinFrameVec);  
        
%         figure; plot(1:size(X,2),X(1,:),'k*'); hold on; plot(1:size(X,2),X_filter(1,:),'ro'); ylim([-0.1 1.1])  
%         figure; plot(1:size(X,2),X(1,:),'k*'); hold on; plot(1:size(X,2),X_ANDgate(1,:),'ro'); ylim([-0.1 1.1])  
%         figure; plot(1:size(X,2),X(1,:),'k*'); hold on; plot(BinFrameVec,X_Bin(1,:),'ro'); ylim([-0.1 1.1])    
    end        
    Data_Binned.(CurrentField) = X_Bin;    
end    

%% Bin Logical Vectors
X = LogicalVector;
if exist('LogicalVector','var')
    X_filter = filter(b_MA,a_MA,double(X)')'; 
    X_ORgate = X_filter > 0;
    X_Bin    = X_ORgate(:,BinFrameVec);
    LogicalVector_OrGateBinned = X_Bin;

    X_filter  = filter(b_MA,a_MA,double(X)')'; 
    X_ANDgate = X_filter > (1-eps);
    X_Bin     = X_ANDgate(:,BinFrameVec);  
    LogicalVector_AndGateBinned = X_Bin;

%         figure; plot(1:size(X,2),X(1,:),'k*'); hold on; plot(1:size(X,2),X_ORgate(1,:),'ro'); ylim([-0.1 1.1])  
%         figure; plot(1:size(X,2),X(1,:),'k*'); hold on; plot(1:size(X,2),X_ANDgate(1,:),'ro'); ylim([-0.1 1.1])  
%         figure; plot(1:size(X,2),X(1,:),'k*'); hold on; plot(BinFrameVec,LogicalVector_OrGateBinned(1,:),'bo');     
%                                                         plot(BinFrameVec,LogicalVector_AndGateBinned(1,:),'rd'); ylim([-0.1 1.1])   
end    

return

%% Plots and related indices sorting and stats
function [Stats, HistogramStructure] = Plot_ModelPredictorsComparison(Data)
%% predictors histogram resolution
Edges_Predictor_HighRes = [-1 -0.4 -0.2 0 0.2 0.4 1]; % For all predictors, New definition: (C-T)/(C+T) 
Edges_Predictor         = [-1 -0.2 0 0.2 1];          % For all predictors, New definition: (C-T)/(C+T) 
% Edges_Predictor         = [-1 -0.3 0 0.3 1];          % For all predictors, New definition: (C-T)/(C+T) 

%% Initialization
MatrixSize                  = size(Data.NaN);
NumberOfFrames              = MatrixSize(2);
NumberOfWorms               = MatrixSize(1);
ActivationInitiationFrames  = Data.ActivationStartFrames;  
ActivationOOB               = Data.ActivationFrameOOB;  
ActivationInitiation_NotOOB                 = ActivationInitiationFrames; 
ActivationInitiation_NotOOB(ActivationOOB)  = false;
Downgradient_NoActivation_ActivityStartsLow = Data.Downgradient_NoActivation & Data.Downgradient_ActivityStartsLow;

%% Remove frames before and after gradient is established, and remove control non-gradient experiments 
GradientFramesMatrix = false(MatrixSize);
for tr=1:NumberOfWorms
    FirstFrameOfOdor     = Data.ExpInfo.GradientStartFrame(tr);   % Gradient STARTS to be established at time 0  
    FirstFrameOfGradient = FirstFrameOfOdor + Data.ExpInfo.GradientRiseTime_InitiationTo90Percent_InSec(tr)*30; % Gradient is established   
    LastFrameOfGradient  = Data.ExpInfo.GradientEndFrame(tr);     % Gradient STARTS to vanish at this frame 
    LastFrameOfGradient  = min([LastFrameOfGradient NumberOfFrames]);
    GradientFrames       = false(1,NumberOfFrames);
    GradientFrames(FirstFrameOfGradient:LastFrameOfGradient) = true;
    GradientFramesMatrix(tr,:) = GradientFrames;
end
ExperimentsWithOdorGradient = Data.ExpInfo.ConcentrationAtArenaTop ~= Data.ExpInfo.ConcentrationAtArenaBottom;
GradientFramesMatrix(~ExperimentsWithOdorGradient, :) = false;
ActivationInitiation_NotOOB(~GradientFramesMatrix)    = false; 

%% Model Predictors 
C                        = Data.C; 
C(~GradientFramesMatrix) = NaN;
Threshold                = Data.AdaptiveThreshold;
dCdt                              = Data.dCdt; 
dCdt(~GradientFramesMatrix)       = NaN;
FoldChange                        = Data.FoldChange; 
FoldChange(~GradientFramesMatrix) = NaN;
% Measured thresholds from previous (immobilzed worms) imaging experiments (Soma) 
Measured_C_Threshold               = 4.915e-7;  % units = dilution.  (=5.4852 microM)
Measured_dCdt_Threshold            = -2.5e-8;   % units = dilution/sec (=-0.2788 microM/s)
Measured_FoldChange_Threshold      = -0.0604;   % units = 1/sec 
% Predictors
Predictor.AdaptiveThreshold        = (C - Threshold) ./ ( C + Threshold );   
Predictor.NonAdaptiveThreshold     = (C - Measured_C_Threshold)./ ( C + Measured_C_Threshold );   
Predictor.Derivative               = -(dCdt - Measured_dCdt_Threshold)  ./ ( dCdt +  Measured_dCdt_Threshold );   
Predictor.Derivative(dCdt>0)       = 1; 
Predictor.FoldChange               = -(FoldChange - Measured_FoldChange_Threshold) ./ ( FoldChange +  Measured_FoldChange_Threshold );   
Predictor.FoldChange(FoldChange>0) = 1; 

%% Figure S4E: 2D histogram (colormaps). Activation probability for 2 predictors
Edges2 = Edges_Predictor;        
plotme = 1;

StructureOut = ComputeAndPlot2DActivationProbability(Edges2, Edges2, Predictor.NonAdaptiveThreshold, Predictor.AdaptiveThreshold, ...
    Downgradient_NoActivation_ActivityStartsLow, ActivationInitiation_NotOOB, plotme);
GreenColormap = CreateWhiteToGreenColormap_inline(true, [0.8 0.8 0.8]);
set(gcf,'colormap',GreenColormap,'name','Figure S4E. ACT versus absolute predictors')
ylabel('\lambda_A_C_T'); xlabel('\lambda_C')

StructureOut = ComputeAndPlot2DActivationProbability(Edges2, Edges2, Predictor.NonAdaptiveThreshold, Predictor.Derivative, ...
    Downgradient_NoActivation_ActivityStartsLow, ActivationInitiation_NotOOB, plotme);
GreenColormap = CreateWhiteToGreenColormap_inline(false);
set(gcf,'colormap',GreenColormap,'name','Figure S4E. Derivative versus absolute predictors')
ylabel('\lambda_d_C_/_d_t'); xlabel('\lambda_C')

StructureOut = ComputeAndPlot2DActivationProbability(Edges2, Edges2, Predictor.NonAdaptiveThreshold, Predictor.FoldChange, ...
    Downgradient_NoActivation_ActivityStartsLow, ActivationInitiation_NotOOB, plotme);
GreenColormap = CreateWhiteToGreenColormap_inline(false);
set(gcf,'colormap',GreenColormap,'name','Figure S4E. Fold-change versus absolute predictors')
ylabel('\lambda_F_o_l_d_C_h_a_n_g_e'); xlabel('\lambda_C')

%% Figure 3E (upper), 1D histograms: DeltaF/F at different predictors ranges
DeltaFOverF                        = Data.DeltaFOverFFiltered;
DeltaFOverF(~GradientFramesMatrix) = NaN;
EdgesDeltaFoverF                   = 0:0.05:5;

[AverageValuesMat, BinVec1, BinVec2, Histograms, BinVec_values] = ComputeAverageWithin2DHistogram_inline ...
                        (Predictor.AdaptiveThreshold, Predictor.NonAdaptiveThreshold, DeltaFOverF, Edges_Predictor_HighRes , [-1.1 1.1], EdgesDeltaFoverF);
CrelativetoT_vec = AverageValuesMat;

[AverageValuesMat, BinVec1, BinVec2, Histograms, BinVec_values] = ComputeAverageWithin2DHistogram_inline ...
                        (Predictor.AdaptiveThreshold, Predictor.NonAdaptiveThreshold, DeltaFOverF, [-1.1 1.1] , Edges_Predictor_HighRes, EdgesDeltaFoverF);
C_vec = AverageValuesMat;

[AverageValuesMat, BinVec1, BinVec2, Histograms, BinVec_values] = ComputeAverageWithin2DHistogram_inline ...
                        (Predictor.AdaptiveThreshold, Predictor.Derivative, DeltaFOverF, [-1.1 1.1] , Edges_Predictor_HighRes, EdgesDeltaFoverF);
dCdt_vec = AverageValuesMat;

[AverageValuesMat, BinVec1, BinVec2, Histograms, BinVec_values] = ComputeAverageWithin2DHistogram_inline ...
                        (Predictor.AdaptiveThreshold, Predictor.FoldChange, DeltaFOverF, [-1.1 1.1] , Edges_Predictor_HighRes, EdgesDeltaFoverF);
FC_vec = AverageValuesMat;

YLIMITS = [0 2];
figure('name','Figure 1E (upper). Average deltaF/F as a function of model predictors','position',[680   192   214   786]); 
subplot(4,1,1)
    bar(CrelativetoT_vec,'edgecolor','k','facecolor','k');
    set(gca,'xlim',[0.25 length(Edges_Predictor_HighRes)-0.25],'Xtick',(0:length(Edges_Predictor_HighRes)-1)+0.5,'xticklabel',Edges_Predictor_HighRes,'XticklabelRotation',90);        
    ylim(YLIMITS);
    ylabel('\DeltaF/F'); xlabel('\lambda_A_C_T')
subplot(4,1,2)
    bar(C_vec,'edgecolor','k','facecolor','k');
    set(gca,'xlim',[0.25 length(Edges_Predictor_HighRes)-0.25],'Xtick',(0:length(Edges_Predictor_HighRes)-1)+0.5,'xticklabel',Edges_Predictor_HighRes,'XticklabelRotation',90);        
    ylim(YLIMITS);
    ylabel('\DeltaF/F'); xlabel('\lambda_C')
subplot(4,1,3)
    bar(dCdt_vec,'edgecolor','k','facecolor','k');
    set(gca,'xlim',[0.25 length(Edges_Predictor_HighRes)-0.25],'Xtick',(0:length(Edges_Predictor_HighRes)-1)+0.5,'xticklabel',Edges_Predictor_HighRes,'XticklabelRotation',90);        
    ylim(YLIMITS);
    ylabel('\DeltaF/F'); xlabel('\lambda_d_C_/_d_T')
subplot(4,1,4)
    bar(FC_vec,'edgecolor','k','facecolor','k');
    set(gca,'xlim',[0.25 length(Edges_Predictor_HighRes)-0.25],'Xtick',(0:length(Edges_Predictor_HighRes)-1)+0.5,'xticklabel',Edges_Predictor_HighRes,'XticklabelRotation',90);        
    ylim(YLIMITS);
    ylabel('\DeltaF/F'); xlabel('\lambda_F_o_l_d_C_h_a_n_g_e')

%% Figure 3E (lower), 1D histograms: Activation probability at different predictors ranges
Edges2 = Edges_Predictor_HighRes;
PredictorFields = {'AdaptiveThreshold','NonAdaptiveThreshold','Derivative','FoldChange'};
h1 = figure('name','Activation Probability Per Second','position',[680   379   245   599]); 
h5 = figure('name','Activation Probability Per 5 Seconds','position',[680   379   245   599]); 

for predictor_ind = 1:length(PredictorFields)
    CurrentFieldName        = PredictorFields{predictor_ind};
    CurrentMat              = Predictor.(CurrentFieldName);
    CurrentMat_NoActivation = CurrentMat(Downgradient_NoActivation_ActivityStartsLow);
    CurrentMat_Activation   = CurrentMat(ActivationInitiation_NotOOB);
    N = histc(CurrentMat_NoActivation(:),Edges2); N_NoActivation = N(1:end-1);
    N = histc(CurrentMat_Activation(:),Edges2);   N_Activation   = N(1:end-1);
    ActivationProbabilityPerFrame         = N_Activation./(N_Activation + N_NoActivation);
    NonActivationProbabilityPerFrame      = 1-ActivationProbabilityPerFrame;
    t= 30;   NonActivationPerTimeInterval = (NonActivationProbabilityPerFrame).^t; ActivationProbability_PerSecond   = 1 - NonActivationPerTimeInterval;
    t= 30*5; NonActivationPerTimeInterval = (NonActivationProbabilityPerFrame).^t; ActivationProbability_Per5Seconds = 1 - NonActivationPerTimeInterval;
    
    figure(h1);
    subplot(length(PredictorFields),1,predictor_ind); bar(ActivationProbability_PerSecond,'edgecolor','k','facecolor','k'); set(gca,'Xtick',(1:length(Edges2))-0.5,'xticklabel',Edges2,'xlim',[0.25 length(Edges2)-0.25],'XticklabelRotation',90); ylim([0 0.55])
    xlabel(['\lambda(',CurrentFieldName,')']);
    ylabel('Prob (1s)');
    figure(h5);
    subplot(length(PredictorFields),1,predictor_ind); bar(ActivationProbability_Per5Seconds,'edgecolor','k','facecolor','k'); set(gca,'Xtick',(1:length(Edges2))-0.5,'xticklabel',Edges2,'xlim',[0.25 length(Edges2)-0.25],'XticklabelRotation',90); ylim([0 1])    
    xlabel(['\lambda(',CurrentFieldName,')']);
    ylabel('Prob (5s)');
end

%% 2D plots (colormaps)
%  probability(high(DeltaF/F)) as a function of ACT and C predictors
CutOff = 2; % DeltaFOverF
GreenColormap = CreateWhiteToGreenColormap_inline(true, [0.8 0.8 0.8]); 
[AverageValuesMat, BinVec1, BinVec2, Histograms, BinVec_values] = ComputeFractionAboveValueWithin2DHistogram_inline ...
                        (Predictor.AdaptiveThreshold, Predictor.NonAdaptiveThreshold, DeltaFOverF, Edges_Predictor , Edges_Predictor, EdgesDeltaFoverF, CutOff);
figure('name','Fraction above cutoff deltaF/F'); 
imagesc_withoutNaNs(1:size(AverageValuesMat,2), 1:size(AverageValuesMat,1),AverageValuesMat,...
                    [min(AverageValuesMat(:)) max(AverageValuesMat(:))],GreenColormap, [1 1 1])
set(gca,'Xtick',(1:length(Edges_Predictor))-0.5,'xticklabel',Edges_Predictor,'Ytick',(1:length(Edges_Predictor))-0.5,'Yticklabel',Edges_Predictor,'XticklabelRotation',90)            
axis xy; axis image; colorbar;
set(gca,'clim',[0.0018 0.3476])
set(gcf,'colormap',GreenColormap)
ylabel('\lambda_A_C_T'); xlabel('\lambda_C')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Predictor values at the time of activity initiation (Figures S4C and S4D) 
ValuesAtActivation.AdaptiveThreshold = Data.AdaptiveThreshold(ActivationInitiation_NotOOB(:));
ValuesAtActivation.C                 = Data.C(ActivationInitiation_NotOOB(:));
ValuesAtActivation.dCdt              = Data.dCdt(ActivationInitiation_NotOOB(:));
ValuesAtActivation.FoldChange        = Data.FoldChange(ActivationInitiation_NotOOB(:));

ValuesAtActivation.Predictor.AdaptiveThreshold    = (ValuesAtActivation.C - ValuesAtActivation.AdaptiveThreshold) ./ (ValuesAtActivation.C + ValuesAtActivation.AdaptiveThreshold);   
ValuesAtActivation.Predictor.NonAdaptiveThreshold = (ValuesAtActivation.C - Measured_C_Threshold) ./ (ValuesAtActivation.C + Measured_C_Threshold);   
ValuesAtActivation.Predictor.Derivative           = -(ValuesAtActivation.dCdt - Measured_dCdt_Threshold)  ./ (ValuesAtActivation.dCdt +  Measured_dCdt_Threshold);   
ValuesAtActivation.Predictor.Derivative(ValuesAtActivation.dCdt>0) = 1; 
ValuesAtActivation.Predictor.FoldChange           = -(ValuesAtActivation.FoldChange - Measured_FoldChange_Threshold) ./ (ValuesAtActivation.FoldChange +  Measured_FoldChange_Threshold);   
ValuesAtActivation.Predictor.FoldChange(ValuesAtActivation.FoldChange>0) = 1; 

binranges  = -1:0.13333:1.01;
bincenters = (binranges(1:end-1)+binranges(2:end))/2;
% Figure S4C
figure('name','Figure S4C: Model predictors values at neural activation time','position',[680    98   386   880]); 
ind = 0;
Stats.binranges      = binranges;
Stats.bincenters     = bincenters;
Stats.N              = zeros(4,length(bincenters))*NaN; 
Stats.Mean           = zeros(1,4)*NaN; 
Stats.Std            = zeros(1,4)*NaN; 
Stats.SEM            = zeros(1,4)*NaN; 
Stats.p              = zeros(1,4)*NaN; 
Stats.Skewness       = zeros(1,4)*NaN; 
Stats.ExcessKurtosis = zeros(1,4)*NaN; 

Vec = ValuesAtActivation.Predictor.AdaptiveThreshold;     
    ind=ind+1;  Stats = BarsAndStats(binranges, Vec, Stats, ind);
Vec = ValuesAtActivation.Predictor.NonAdaptiveThreshold;    
    ind=ind+1;  Stats = BarsAndStats(binranges, Vec, Stats, ind);
Vec = ValuesAtActivation.Predictor.Derivative;   
    ind=ind+1;  Stats = BarsAndStats(binranges, Vec, Stats, ind);
Vec = ValuesAtActivation.Predictor.FoldChange;    
    ind=ind+1;  Stats = BarsAndStats(binranges, Vec, Stats, ind);

Stats.Control.N = ones(1,size(Stats.N,2))/size(Stats.N,2);
% Percentage around Lambda=0.5
Vec = ValuesAtActivation.Predictor.AdaptiveThreshold;     
PercentageActivatedAtThreshold = length(find(Vec<=0.333 & Vec>=-0.333))./length(Vec)*100;
PercentageActivatedAtThreshold

%% Figure S4D
% Here the Jehnsen-Shannon divergence can be used to see the difference from a flat randomly spread statistics       
% This will show that some other model predictors are random  
% figure; bar(Stats.bincenters,Stats.Control.N)
DjsVec    = zeros(1,size(Stats.N,1))*NaN;
for ind = 1:4
    DjsVec(ind) = JehnsenShannonDivergence_inline (Stats.N(ind,:),Stats.Control.N);
end

% Control Djs (flat distribution)
SampleSize      = length(Vec);
NumOfIterations = 1e3;
RandMat         = 2*rand(NumOfIterations,SampleSize)-1;
Nvec            = ones(NumOfIterations,length(bincenters))*NaN; 
DjsVecControl   = ones(NumOfIterations,1)*NaN; 
for ind=1:1e3
    N          = histcounts(RandMat(ind,:), binranges);    
    Nvec(ind,:)= N/sum(N); 
    DjsVecControl(ind) = JehnsenShannonDivergence_inline (Nvec(ind,:),Stats.Control.N);
end

% Control Djs (Gaussian distribution)
Stats.GaussianControl.N = normpdf(bincenters,0,Stats.Std(1));  
Stats.GaussianControl.N = Stats.GaussianControl.N / sum(Stats.GaussianControl.N);
DjsVecGaussianControl   = JehnsenShannonDivergence_inline (Stats.GaussianControl.N,Stats.Control.N);

% Plot Figure S4D
figure('name','bars=Models, low line= lower significance limit (2 std from random uniform sample), high line=gaussian (same std as ACT model)','position', [680   558   335   420]); 
bar(DjsVec,'facecolor',[0.5 0.5 0.5]);         
hold on; plot(get(gca,'xlim'),(mean(DjsVecControl)+2*std(DjsVecControl))*ones(1,2),'k-');
hold on; plot(get(gca,'xlim'),mean(DjsVecGaussianControl)*ones(1,2),'k-');
ylabel('Djs'); title('Distance of predictors distribution from uniform random distribution')
xlim([-0.25 length(Stats.Mean)+0.25]+0.5)  
set(gca,'xtick',1:4,'xticklabel',{'ACT','C','dC/dt','(dC/dt)/C'})


figure('name','Mean and SEM','position', [680   558   335   420]); 
bar(1:length(Stats.Mean),Stats.Mean,'facecolor',[0.5 0.5 0.5],'BaseValue',-1); hold on;  
errorbar(1:length(Stats.Mean),Stats.Mean,Stats.SEM,'k.');  
xlim([-0.25 length(Stats.Mean)+0.25]+0.5)        
line(get(gca,'xlim'),[0 0],'color',[0.2 0.2 0.2])
ylabel('mean(\lambda)');
set(gca,'xtick',1:4,'xticklabel',{'ACT','C','dC/dt','(dC/dt)/C'})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Histograms with SEM calculated by bootstrap for Figures 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Initialization and resizing Matrices - omit non gradient experiments
PredictorMatrix                             = Predictor.AdaptiveThreshold(ExperimentsWithOdorGradient,:);
DeltaFOverF                                 = DeltaFOverF(ExperimentsWithOdorGradient,:);
Downgradient_NoActivation_ActivityStartsLow = Downgradient_NoActivation_ActivityStartsLow(ExperimentsWithOdorGradient,:);
ActivationInitiation_NotOOB                 = ActivationInitiation_NotOOB(ExperimentsWithOdorGradient,:);

NumOfIterations  = 1e3;
NumOfExperiments = size(PredictorMatrix,1);

binranges  = -1:0.133333:1.01;
bincenters = (binranges(1:end-1)+binranges(2:end))/2;
HistogramStructure.PredictorValuesAtActivation.Edges   = binranges;
HistogramStructure.PredictorValuesAtActivation.Centers = bincenters;

%% Figure S4C
% No bootstrap
PredictorValuesAtActivation = PredictorMatrix(ActivationInitiation_NotOOB(:));   
N = histcounts(PredictorValuesAtActivation, binranges);    
N = N/sum(N); 
HistogramStructure.PredictorValuesAtActivation.UsingAllData = N;

% Bootstrap
tic
rng('default'); rng(1);
RandMatrix = randi(NumOfExperiments, NumOfExperiments, NumOfIterations); % random integers between 1 and NumOfExperiments in matrix(NumOfExperiments x NumOfIterations) 
PDF_bootstrap_mat = zeros(length(bincenters),NumOfIterations)*NaN;
for n=1:NumOfIterations
    CurrentExperimentNumbers           = RandMatrix(:,n);
    CurrentPredictorMatrix             = PredictorMatrix(CurrentExperimentNumbers,:);
    CurrentActivationInitiation_NotOOB = ActivationInitiation_NotOOB(CurrentExperimentNumbers,:);
    CurrentPredictorValuesAtActivation = CurrentPredictorMatrix(CurrentActivationInitiation_NotOOB(:));   
    N                      = histcounts(CurrentPredictorValuesAtActivation, binranges);    
    N                      = N/sum(N);     
    PDF_bootstrap_mat(:,n) = N;
end
toc
HistogramStructure.PredictorValuesAtActivation.bootstrap_mat = PDF_bootstrap_mat;
HistogramStructure.PredictorValuesAtActivation.MeanBootstrap = mean(PDF_bootstrap_mat,2);
HistogramStructure.PredictorValuesAtActivation.STDBootstrap  = std(PDF_bootstrap_mat,[],2);
HistogramStructure.PredictorValuesAtActivation.SEMBootstrap  = std(PDF_bootstrap_mat,[],2) / sqrt(NumOfExperiments);

figure('name','Figure S4C: Model predictors values at neural activation time');
bar(HistogramStructure.PredictorValuesAtActivation.UsingAllData,'edgecolor','k','facecolor','k'); hold on;
errorbar(HistogramStructure.PredictorValuesAtActivation.UsingAllData, ...
         HistogramStructure.PredictorValuesAtActivation.SEMBootstrap, 'linestyle','none', 'color','k'); 
set(gca,'xlim',[0.25 length(binranges)-0.25],'Xtick',(0:length(binranges)-1)+0.5,'xticklabel',binranges,'XticklabelRotation',90);        

%% Figure 3E, 1D histograms: average deltaF/F as a function of ACT predictor
HistogramStructure.ActivationProb5sec.Edges = Edges_Predictor_HighRes;

% No bootstrap
[AverageValuesMat, BinVec1, BinVec2, Histograms, BinVec_values] = ComputeAverageWithin2DHistogram_inline ...
                        (PredictorMatrix, PredictorMatrix, DeltaFOverF, Edges_Predictor_HighRes , [-1.1 1.1]);
HistogramStructure.MeanDeltaFOverF.UsingAllData = AverageValuesMat;

% Bootstrap
tic
rng('default'); rng(1);
RandMatrix = randi(NumOfExperiments, NumOfExperiments, NumOfIterations); % random integers between 1 and NumOfExperiments in matrix(NumOfExperiments x NumOfIterations) 
DeltaFOverF_bootstrap_mat = zeros(length(AverageValuesMat),NumOfIterations)*NaN;
for n=1:NumOfIterations
    CurrentExperimentNumbers = RandMatrix(:,n);
    CurrentPredictorMatrix   = PredictorMatrix(CurrentExperimentNumbers,:);
    CurrentDeltaFOverF       = DeltaFOverF(CurrentExperimentNumbers,:);
    CurrentAverageValuesMat  = ComputeAverageWithin2DHistogram_inline ...
                 (CurrentPredictorMatrix, CurrentPredictorMatrix, CurrentDeltaFOverF, Edges_Predictor_HighRes , [-1.1 1.1]);
    DeltaFOverF_bootstrap_mat(:,n)   = CurrentAverageValuesMat;
end
toc
HistogramStructure.MeanDeltaFOverF.bootstrap_mat = DeltaFOverF_bootstrap_mat;
HistogramStructure.MeanDeltaFOverF.MeanBootstrap = mean(DeltaFOverF_bootstrap_mat,2);
HistogramStructure.MeanDeltaFOverF.STDBootstrap  = std(DeltaFOverF_bootstrap_mat,[],2);
HistogramStructure.MeanDeltaFOverF.SEMBootstrap  = std(DeltaFOverF_bootstrap_mat,[],2) / sqrt(NumOfExperiments);

YLIMITS = [0 2];
figure('name','Average deltaF/F as a function of model predictors');%,'position',[680   192   214   786]); 
bar(HistogramStructure.MeanDeltaFOverF.UsingAllData,'edgecolor','k','facecolor','k'); hold on;
errorbar(HistogramStructure.MeanDeltaFOverF.UsingAllData, ...
         HistogramStructure.MeanDeltaFOverF.SEMBootstrap, 'linestyle','none', 'color','k'); 
set(gca,'xlim',[0.25 length(Edges_Predictor_HighRes)-0.25],'Xtick',(0:length(Edges_Predictor_HighRes)-1)+0.5,'xticklabel',Edges_Predictor_HighRes,'XticklabelRotation',90);        
ylim(YLIMITS);
xlabel('\lambda'); ylabel('\DeltaF/F')


%% Figure 3E, 1D histogram: Activation probability as a function of ACT predictor 
TimeInterval = 5; % sec. Activation probability in 5 seconds
h5 = figure('name','Activation Probability Per 5 Seconds');%,'position',[680   379   245   599]); 
HistogramStructure.ActivationProb5sec.TimeInterval = TimeInterval;
HistogramStructure.ActivationProb5sec.Edges        = Edges_Predictor_HighRes;

% No bootstrap
CurrentMat_NoActivation = PredictorMatrix(Downgradient_NoActivation_ActivityStartsLow);
CurrentMat_Activation   = PredictorMatrix(ActivationInitiation_NotOOB);
ActivityProbability     = CalculateActivationProbability(CurrentMat_NoActivation, CurrentMat_Activation, Edges_Predictor_HighRes, TimeInterval);
HistogramStructure.ActivationProb5sec.UsingAllData = ActivityProbability;

% Bootstrap
tic
rng('default'); rng(1);
NumOfIterations  = 1e3;
NumOfExperiments = size(PredictorMatrix,1);
RandMatrix = randi(NumOfExperiments, NumOfExperiments, NumOfIterations); % random integers between 1 and NumOfExperiments in matrix(NumOfExperiments x NumOfIterations) 
PDF_bootstrap_mat = zeros(length(ActivityProbability),NumOfIterations)*NaN;
for n=1:NumOfIterations
    CurrentExperimentNumbers = RandMatrix(:,n);
    CurrentPredictorMatrix                     = PredictorMatrix(CurrentExperimentNumbers,:);
    CurrentDowngradientNoActivationStartsLow   = Downgradient_NoActivation_ActivityStartsLow(CurrentExperimentNumbers,:);
    CurrentActivityInitiation                  = ActivationInitiation_NotOOB(CurrentExperimentNumbers,:);
    CurrentMat_NoActivation  = CurrentPredictorMatrix(CurrentDowngradientNoActivationStartsLow);
    CurrentMat_Activation    = CurrentPredictorMatrix(CurrentActivityInitiation);
    CurrentPDF               = CalculateActivationProbability(CurrentMat_NoActivation, CurrentMat_Activation, Edges_Predictor_HighRes, TimeInterval);
    PDF_bootstrap_mat(:,n)   = CurrentPDF;
end
toc
HistogramStructure.ActivationProb5sec.bootstrap_mat = PDF_bootstrap_mat;
HistogramStructure.ActivationProb5sec.MeanBootstrap = mean(PDF_bootstrap_mat,2);
HistogramStructure.ActivationProb5sec.STDBootstrap  = std(PDF_bootstrap_mat,[],2);
HistogramStructure.ActivationProb5sec.SEMBootstrap  = std(PDF_bootstrap_mat,[],2) / sqrt(NumOfExperiments);

figure(h5);
bar(HistogramStructure.ActivationProb5sec.UsingAllData ,'edgecolor','k','facecolor','k'); hold on;
errorbar(HistogramStructure.ActivationProb5sec.UsingAllData, ...
         HistogramStructure.ActivationProb5sec.SEMBootstrap, 'linestyle','none', 'color','k'); 
set(gca,'Xtick',(1:length(Edges_Predictor_HighRes))-0.5,'xticklabel',Edges_Predictor_HighRes,'xlim',[0.25 length(Edges_Predictor_HighRes)-0.25],'XticklabelRotation',90); ylim([0 1])
xlabel('\lambda'); ylabel('Prob (5s)')

return

function PDF = CalculateActivationProbability(CurrentMat_NoActivation, CurrentMat_Activation, Edges2, TimeInterval)

N = histc(CurrentMat_NoActivation(:),Edges2); N_NoActivation = N(1:end-1);
N = histc(CurrentMat_Activation(:),Edges2);   N_Activation   = N(1:end-1);
ActivationProbabilityPerFrame         = N_Activation./(N_Activation + N_NoActivation);
NonActivationProbabilityPerFrame      = 1-ActivationProbabilityPerFrame;
t= 30*TimeInterval; 
NonActivationPerTimeInterval = (NonActivationProbabilityPerFrame).^t; 
PDF = 1 - NonActivationPerTimeInterval;  % ActivationProbability_Per5Seconds

return

function [pValue, KSstatistic]= kstest2_inline(sampleCDF1, sampleCDF2, n1, n2, tail)

%
% Compute the test statistic of interest.
%

switch tail
   case  0      %  2-sided test: T = max|F1(x) - F2(x)|.
      deltaCDF  =  abs(sampleCDF1 - sampleCDF2);

   case -1      %  1-sided test: T = max[F2(x) - F1(x)].
      deltaCDF  =  sampleCDF2 - sampleCDF1;

   case  1      %  1-sided test: T = max[F1(x) - F2(x)].
      deltaCDF  =  sampleCDF1 - sampleCDF2;
end

KSstatistic   =  max(deltaCDF);

%
% Compute the asymptotic P-value approximation and accept or
% reject the null hypothesis on the basis of the P-value.
%

% n1     =  length(x1);
% n2     =  length(x2);
n      =  n1 * n2 /(n1 + n2);
lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic , 0);

if tail ~= 0        % 1-sided test.

   pValue  =  exp(-2 * lambda * lambda);

else                % 2-sided test (default).
%
%  Use the asymptotic Q-function to approximate the 2-sided P-value.
%
   j       =  (1:101)';
   pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
   pValue  =  min(max(pValue, 0), 1);

end

% H  =  (alpha >= pValue);
return

function [CDF_ForwardTermination, CDF_ReverseTermination, SortedTimes_ForwardTermination, SortedTimes_ReverseTermination] = PlotMatrices_AlignByActivityOrPrediction(StructureOut, MinFramesInBoundFromEventStart, PlotInfo)

%% Initialization
Concentration       = StructureOut.Concentration;
FoldChange          = StructureOut.FoldChange;
Behavior            = StructureOut.Behavior;
NeuronValue         = StructureOut.NeuronValue;
F0_mat              = repmat(StructureOut.EventPastNeuronValue,1,size(NeuronValue,2));
RelativeDeltaFOverF = (NeuronValue-F0_mat)./F0_mat;
ACT_predictor       = (StructureOut.EventConcentration - StructureOut.EventThreshold)./(StructureOut.EventConcentration + StructureOut.EventThreshold); % 
Time                = StructureOut.Time;
ZeroIndex           = find(Time==0);
XLIMITS             = PlotInfo.XLIMITS;
NeuronClimit        = PlotInfo.NeuronClimit;
FigureTitleStr      = PlotInfo.FigureTitle;

load('ColormapLowLevelBehavior.mat','BehaviorColormapLowLevel');

FramesInBoundFromEventStart = zeros(size(Behavior,1),1)*NaN;
logicalMat = Behavior==0;
for row=1:length(FramesInBoundFromEventStart)
    Frame = find(logicalMat(row,ZeroIndex:end),1,'first')-1;
    if~isempty(Frame)
        FramesInBoundFromEventStart(row)= Frame;
    end
end
EventsInBound = isnan(FramesInBoundFromEventStart)|(FramesInBoundFromEventStart>MinFramesInBoundFromEventStart);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%% Plot activity sort by ACT_predictor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
if PlotInfo.PlotActivity 

        [~,SortByDistanceToThreshold] = sort(ACT_predictor);
    SortByDistanceToThreshold     = SortByDistanceToThreshold((StructureOut.EventConcentration(SortByDistanceToThreshold)>0)&EventsInBound);
    [~, ~, ~, RelativeDeltaFOverF_ACTsorted] = ...
        SortMatricesAndPlot(Concentration, FoldChange, Behavior, RelativeDeltaFOverF, SortByDistanceToThreshold, false, '', Time, BehaviorColormapLowLevel, XLIMITS, NeuronClimit);

    % colormap
    FigureTitle = [FigureTitleStr, ', DeltaF/F (sorted By ACT_predictor)'];
    Mat         = RelativeDeltaFOverF_ACTsorted;
    figure('name',FigureTitle,'position',[680   615   325   355]);
        imagesc(Time,1:size(Mat,1),Mat);       
        colorbar; set(gca,'clim',NeuronClimit); xlim(XLIMITS);
        set(gca,'xtick',XLIMITS(1):5:XLIMITS(end))

    % Average
    FigureTitle = [FigureTitleStr, ', DeltaF/F'];
    figure('name',FigureTitle,'position',[680   340   325   200]); 
    MEAN  = nanmean(Mat,1);
    SEM   = nanstd(Mat,1,1)./sqrt(size(Mat,1));
    upper = MEAN + SEM;
    lower = MEAN - SEM;
    ylim([0 1]);
    jbfill(Time,upper,lower,[0.5 0.5 0.5],[0.2 0.2 0.2]); hold on;
    plot(Time,MEAN,'k-','linewidth',2);
    xlim(XLIMITS); 
    set(gca,'xtick',XLIMITS(1):5:XLIMITS(end))

    % CDF of traces that reach deltaF/F=ThresholdValue 
    ThresholdValue         = 0.6;
    [~, SortedTimes, ~]    = FindReachThresholdIndices(Mat, ZeroIndex, ThresholdValue);

    NumOfTraces            = length(SortedTimes);
    CumulativeDistribution = zeros(1,(size(Mat,2)-ZeroIndex));
    LastNumberOfTracks = 0; Last_ind = 1;
    for frame=SortedTimes
        if frame==0
            CumulativeDistribution(1)=1;
            continue
        end
        if ~isnan(frame)
            CumulativeDistribution(Last_ind:frame-1)= LastNumberOfTracks;
            CumulativeDistribution(frame)           = LastNumberOfTracks+1;
            LastNumberOfTracks                      = LastNumberOfTracks+1;
            Last_ind                                = frame;
        else
            % no more tracks terminated
             CumulativeDistribution(Last_ind:end)   = LastNumberOfTracks;
             break
        end
    end
    CumulativeDistribution = CumulativeDistribution/NumOfTraces;

    FigureTitle = [FigureTitleStr, ', CDF(DeltaF/F)>0.6'];
    figure('name',FigureTitle,'position',[680   100   325   200]); 
    plot(Time(ZeroIndex+1:end), CumulativeDistribution, 'k-','linewidth',2); hold on;
    ylabel('Cumulative distribution'); xlabel('Time [sec]'); ylim([0 0.6])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%% Split Matrices to behavior start code, and then sort by time of behavior termination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
if PlotInfo.PlotBehavior

    NonZeroOdorConcentration     = find((StructureOut.EventConcentration>0)&EventsInBound);
    [Concentration_NonZero, FoldChange_NonZero, Behavior_NonZero, RelativeDeltaFOverF_NonZero] = ...
        SortMatricesAndPlot(Concentration, FoldChange, Behavior, RelativeDeltaFOverF, NonZeroOdorConcentration, false, '', Time, BehaviorColormapLowLevel, XLIMITS, NeuronClimit);
    EventBehavior_NonZero = StructureOut.EventBehavior(NonZeroOdorConcentration);

    %% Matrices starting in Forward state %%  
    % Event starting in forward state
    StartingForward                        = find(EventBehavior_NonZero==1);
    [~, ~, ~, TerminationBehavior] = FindBehaviorTerminationIndices(Behavior_NonZero(StartingForward,:), ZeroIndex); % eliminate behaviors that stop due to locomtion out of bounds. 
    StartingForward  = StartingForward(TerminationBehavior~=0);

    % Sort by forward termination
    [Concentration_NonZero_StartFwd, FoldChange_NonZero_StartFwd, Behavior_NonZero_StartFwd, RelativeDeltaFOverF_NonZero_StartFwd] = ...
        SortMatricesAndPlot(Concentration_NonZero, FoldChange_NonZero, Behavior_NonZero, RelativeDeltaFOverF_NonZero, StartingForward, false, '', Time, BehaviorColormapLowLevel, XLIMITS, NeuronClimit);
    [SortByForwardTermination, SortedTimes] = FindBehaviorTerminationIndices(Behavior_NonZero_StartFwd, ZeroIndex);

    % Plot behavior
    [~, ~, Behavior_sorted] = SortMatricesAndPlot(Concentration_NonZero_StartFwd, FoldChange_NonZero_StartFwd, Behavior_NonZero_StartFwd, RelativeDeltaFOverF_NonZero_StartFwd, SortByForwardTermination, false, '', Time, BehaviorColormapLowLevel, XLIMITS, NeuronClimit);
    FigureTitle = [FigureTitleStr, ', Behavior matrix. Only event starting in FORWARD state and sorted by FORWARD termination'];
    figure('name',FigureTitle); PlotBehaviorMatrix(Behavior_sorted, BehaviorColormapLowLevel, Time, 1:size(Behavior_sorted,1)); xlim(XLIMITS);

    % If matrices are huge, sub-sample 1% of the data only for behavior plot display
    TotalNumberOfEvent = length(SortByForwardTermination);
    if TotalNumberOfEvent> 1e3
        NumberOfEventsToDisplay = 100;     
        Indices           = round(linspace(1,TotalNumberOfEvent,NumberOfEventsToDisplay)); 
        [~, ~, Behavior_sorted] = SortMatricesAndPlot(Concentration_NonZero_StartFwd, FoldChange_NonZero_StartFwd, Behavior_NonZero_StartFwd, RelativeDeltaFOverF_NonZero_StartFwd, SortByForwardTermination(Indices), false, '', Time, BehaviorColormapLowLevel, XLIMITS, NeuronClimit);
        FigureTitle = [FigureTitleStr, ', Subsampled behavior matrix. Only event starting in FORWARD state and sorted by FORWARD termination'];
        figure('name',FigureTitle); PlotBehaviorMatrix(Behavior_sorted, BehaviorColormapLowLevel, Time, 1:size(Behavior_sorted,1)); xlim(XLIMITS);
    end

    % Cumulative distribution
    NumOfTraces            = length(SortedTimes);
    CumulativeDistribution = zeros(1,(size(Behavior_NonZero,2)-ZeroIndex));
    LastNumberOfTracks = 0; Last_ind = 1;
    for frame=SortedTimes
        if frame==0
            CumulativeDistribution(1)=1;
            continue
        end
        if ~isnan(frame)
            CumulativeDistribution(Last_ind:frame-1)= LastNumberOfTracks;
            CumulativeDistribution(frame)           = LastNumberOfTracks+1;
            LastNumberOfTracks                      = LastNumberOfTracks+1;
            Last_ind                                = frame;
        else
            % no more tracks terminated
             CumulativeDistribution(Last_ind:end)   = LastNumberOfTracks;
             break
        end
    end
    CDF_ForwardTermination         = CumulativeDistribution/NumOfTraces;
    SortedTimes_ForwardTermination = SortedTimes;

    %% Matrices starting in Reverse state %%  
    % Event starting in reverse state
    StartingReverse       = find((EventBehavior_NonZero==4));
    [~, ~, ~, TerminationBehavior] = FindBehaviorTerminationIndices(Behavior_NonZero(StartingReverse,:), ZeroIndex); % eliminate behaviors that stop due to locomotion out of bounds.
    StartingReverse  = StartingReverse(TerminationBehavior~=0);

    % Sort by reverse termination
    [Concentration_NonZero_StartRev, FoldChange_NonZero_StartRev, Behavior_NonZero_StartRev, RelativeDeltaFOverF_NonZero_StartRev] = ...
        SortMatricesAndPlot(Concentration_NonZero, FoldChange_NonZero, Behavior_NonZero, RelativeDeltaFOverF_NonZero, StartingReverse, false, '', Time, BehaviorColormapLowLevel, XLIMITS, NeuronClimit);
    [SortByForwardInitiation, SortedTimes] = FindBehaviorInitiationIndices(Behavior_NonZero_StartRev, ZeroIndex, [1 2]);

    % Plot behavior
    [~, ~, Behavior_sorted] = SortMatricesAndPlot(Concentration_NonZero_StartRev, FoldChange_NonZero_StartRev, Behavior_NonZero_StartRev, RelativeDeltaFOverF_NonZero_StartRev, SortByForwardInitiation, false, '', Time, BehaviorColormapLowLevel, XLIMITS, NeuronClimit);
    FigureTitle = [FigureTitleStr, ', Behavior matrix. Only event starting in REVERSE state and sorted by REVERSE termination'];
    figure('name',FigureTitle); PlotBehaviorMatrix(Behavior_sorted, BehaviorColormapLowLevel, Time, 1:size(Behavior_sorted,1)); xlim(XLIMITS);

    % If matrices are huge, sub-sample 1% of the data only for behavior plot display
    TotalNumberOfEvent = length(SortByForwardInitiation);
    if TotalNumberOfEvent> 1e3
        NumberOfEventsToDisplay = 100;          
        Indices           = round(linspace(1,TotalNumberOfEvent,NumberOfEventsToDisplay)); 
        [~, ~, Behavior_sorted] = SortMatricesAndPlot(Concentration_NonZero_StartRev, FoldChange_NonZero_StartRev, Behavior_NonZero_StartRev, RelativeDeltaFOverF_NonZero_StartRev, SortByForwardInitiation(Indices), false, '', Time, BehaviorColormapLowLevel, XLIMITS, NeuronClimit);
        FigureTitle = [FigureTitleStr, ', subsampled behavior matrix. Only event starting in REVERSE state and sorted by REVERSE termination'];
        figure('name',FigureTitle); PlotBehaviorMatrix(Behavior_sorted, BehaviorColormapLowLevel, Time, 1:size(Behavior_sorted,1)); xlim(XLIMITS);
    end

    % Cumulative distribution
    NumOfTraces            = length(SortedTimes);
    CumulativeDistribution = zeros(1,(size(Behavior_NonZero,2)-ZeroIndex));
    LastNumberOfTracks = 0; Last_ind = 1;
    for frame=SortedTimes
        if frame==0
            CumulativeDistribution(1)=1;
            continue
        end
        if ~isnan(frame)
            CumulativeDistribution(Last_ind:frame-1)= LastNumberOfTracks;
            CumulativeDistribution(frame)           = LastNumberOfTracks+1;
            LastNumberOfTracks                      = LastNumberOfTracks+1;
            Last_ind                                = frame;
        else
            % no more tracks terminated
             CumulativeDistribution(Last_ind:end)   = LastNumberOfTracks;
             break
        end
    end
    CDF_ReverseTermination         = CumulativeDistribution/NumOfTraces;
    SortedTimes_ReverseTermination = SortedTimes;

    %% Plot CDFs for forward and reversal initiation %%       
    FigureTitle = [FigureTitleStr, ', CDFs for forward and reverse termination'];        
    figure('name',FigureTitle); 
    plot(Time(ZeroIndex+1:end), CDF_ForwardTermination, 'k-','linewidth',2); hold on;
    plot(Time(ZeroIndex+1:end), CDF_ReverseTermination, 'r-','linewidth',2);
    ylabel('CDF'); xlabel('Time [sec]'); legend('forward termination','reverse termination')
else
    CDF_ForwardTermination =[];
    CDF_ReverseTermination =[];
    SortedTimes_ForwardTermination =[];
    SortedTimes_ReverseTermination =[];      
end

return

function PlotMatrices_AlignedByBehavior(StructureOut, MinFramesInBoundFromEventStart, PlotInfo)

%% Initialization
Behavior            = StructureOut.Behavior;
NeuronValue         = StructureOut.NeuronValue;
F0_mat              = repmat(StructureOut.EventPastNeuronValue,1,size(NeuronValue,2));
RelativeDeltaFOverF = (NeuronValue-F0_mat)./F0_mat;
ACT_predictor       = (StructureOut.EventConcentration - StructureOut.EventThreshold)./(StructureOut.EventConcentration + StructureOut.EventThreshold); % 
Time                = StructureOut.Time;
ZeroIndex           = find(Time==0,1);
XLIMITS             = PlotInfo.XLIMITS;
FigureTitleStr      = PlotInfo.FigureTitle;

load('ColormapLowLevelBehavior.mat','BehaviorColormapLowLevel');

FramesInBoundFromEventStart = zeros(size(NeuronValue,1),1)*NaN;
logicalMat = Behavior==0;
for row=1:length(FramesInBoundFromEventStart)
    Frame = find(logicalMat(row,ZeroIndex:end),1,'first')-1;
    if~isempty(Frame)
        FramesInBoundFromEventStart(row)= Frame;
    end
end
EventsInBound = isnan(FramesInBoundFromEventStart)|(FramesInBoundFromEventStart>MinFramesInBoundFromEventStart);

% Keep only in bound events 
RelativeDeltaFOverF = RelativeDeltaFOverF(EventsInBound,:);
ACT_predictor       = ACT_predictor(EventsInBound);

%% Figure 4C
for FigType=1:3
    switch FigType
        case 1                
            FigureTitle    = [FigureTitleStr, ', \lambda < -0.2'];
            CurrentIndices = find(ACT_predictor<-0.2);
        case 2
            FigureTitle    = [FigureTitleStr, '-0.2 <= \lambda  <= 0.2'];
            CurrentIndices = find((ACT_predictor>=-0.2) & (ACT_predictor<=0.2));
        case 3
            FigureTitle    = [FigureTitleStr, '\lambda > 0.2'];
            CurrentIndices = find(ACT_predictor>0.2);
    end               
    Mat = RelativeDeltaFOverF(CurrentIndices,:);
    
%     figure('name',[FigureTitle,'. Relative DeltaF/F']);
%         imagesc(Time,1:length(CurrentIndices),Mat);       
%         colorbar; set(gca,'clim',NeuronClimit); xlim(XLIMITS);
%         set(gca,'xtick',XLIMITS(1):5:XLIMITS(end))

    figure('name',FigureTitle,'position',[680   300   239   200]); 
        MEAN  = nanmean(Mat,1);
        SEM   = nanstd(Mat,1,1)./sqrt(size(Mat,1));
        upper = MEAN + SEM;
        lower = MEAN - SEM;
        ylim([-0.31 0.5]);
        jbfill(Time,upper,lower,[0.5 0.5 0.5],[0.2 0.2 0.2]); hold on;
        plot(Time,MEAN,'k-','linewidth',2);
        xlim(XLIMITS); 
        line(XLIMITS,[0 0],'color','k','linestyle',':')
        set(gca,'xtick',XLIMITS(1):5:XLIMITS(end))
end    

%% Figure S4F:  Fraction increasing: average(deltaF/F(t<0)) < average(deltaF/F(t>0))    
TimeIndices_0to2sec          = ZeroIndex:(find(Time>=2,1,'first'));
TimeIndices_minus2to0sec     = (find(Time<=-2,1,'last')):ZeroIndex;
AverageChangeInDeltaFoverF_AfterVsBefore = nanmean(RelativeDeltaFOverF(:,TimeIndices_0to2sec),2) - ...
                                           nanmean(RelativeDeltaFOverF(:,TimeIndices_minus2to0sec),2);

AverageChanges{1} = AverageChangeInDeltaFoverF_AfterVsBefore(ACT_predictor<=-0.2);
AverageChanges{2} = AverageChangeInDeltaFoverF_AfterVsBefore((ACT_predictor>-0.2)&(ACT_predictor<0.2));
AverageChanges{3} = AverageChangeInDeltaFoverF_AfterVsBefore(ACT_predictor>=0.2);
MinimumChangeForAnalysis = 0.1; % for baseline near 0
% bootstrap on repeats 
Bootstrap_NumOfIterations = 1000;
RandomDrawsMatrices       = cell(1,length(AverageChanges));
for vec_ind = 1:length(AverageChanges)
    Vec         = AverageChanges{vec_ind};
    LENGTH      = length(Vec);    
    Indices     = 1:LENGTH;
    RandomDraws = datasample(Indices,Bootstrap_NumOfIterations*LENGTH);      % vector of indices
    RandomDraws = reshape(RandomDraws,[LENGTH, Bootstrap_NumOfIterations]);  % Matrix of indices, each column is one bootstrap
    RandomDrawsMatrices{vec_ind} = RandomDraws;
end

Fraction_bootstrap = zeros(length(AverageChanges),Bootstrap_NumOfIterations,1)*NaN;   
LENs               = zeros(length(AverageChanges),1)*NaN;
for vec_ind = 1:length(AverageChanges)  
    Vec = AverageChanges{vec_ind};
    LEN = length(Vec);
    LENs(vec_ind) = LEN;
    for n = 1:Bootstrap_NumOfIterations  
        CurrentIndices = RandomDrawsMatrices{vec_ind}(:,n);
        CurrentVec     = Vec(CurrentIndices); 
        Fraction_bootstrap(vec_ind,n) = length(find(CurrentVec>MinimumChangeForAnalysis)) / LEN;
    end
end

MEANs = mean(Fraction_bootstrap,2);
STDs  = std(Fraction_bootstrap,1,2);
% SEMs  = STDs ./ sqrt(LENs);

% Figure S4F
figure('name','Bars: How ACT predictor affects fraction of increasing deltaF/F'); 
bar(1:3,MEANs,'edgecolor','k','facecolor',[0.7 0.7 0.7]); hold on; % ylim([0.5 1]);
errorbar(1:3,MEANs,STDs,'k','linestyle','none'); set(gca,'xtick',1:3,'xticklabel',{'[-1, -0.2]','(-0.2, 0.2)','[0.2, 1]'});
ylabel('Activation probability')

return

function [N, Vec1, Vec2] = Compute2DHistogram_inline (Mat1, Mat2, Edges1, Edges2, normalization)

% Mat1 is two dimensional: rows=worms, columns = single frames. 
% Mat2 is two dimensional: rows=worms, columns = single frames. 

Vec1       = (Edges1(1:end-1)+Edges1(2:end))/2;
Vec2       = (Edges2(1:end-1)+Edges2(2:end))/2;
NumOfBins1 = length(Vec1);
NumOfBins2 = length(Vec2);

N = hist3([Mat1(:) Mat2(:)],'Edges',{Edges1, Edges2});
N = N(1:end-1,1:end-1);

if     normalization == 1    % P(X|Y)= For a given Y range, the sum over all X range will be 1. i.e. each ROW in the images shown with 'Plot2Dstats' will sum up to 1.  
        TotalNumber      = sum(N,1);
        TotalNumber_Mat  = repmat(TotalNumber,[NumOfBins1 1]);    
elseif normalization == 2   % P(Y|X)= For a given X range, the sum over all Y range will be 1. i.e. each Column in the images shown with 'Plot2Dstats' will sum up to 1.  
        TotalNumber      = sum(N,2);   
        TotalNumber_Mat  = repmat(TotalNumber,[1 NumOfBins2]);    
elseif normalization == 3   % P(Y&X)= The sum over ALL numbers will be 1. This is the true probability of being at state (x,y). i.e. the sum over all values in the images shown with 'Plot2Dstats' will be 1.   
        TotalNumber      = sum(N(:));
        TotalNumber_Mat  = TotalNumber*ones(size(N));
else
    TotalNumber_Mat = ones(size(N)); % No Normalization
end
N = N./TotalNumber_Mat;


return

function [AverageValuesMat, BinVec1, BinVec2, Histograms, BinVec_values] = ComputeAverageWithin2DHistogram_inline ...
              (Mat1, Mat2, ValueMat, Edges1, Edges2, EdgesForValuesHistograms)

% Mat1 is two dimensional: rows=worms, columns = single frames. 
% Mat2 is two dimensional: rows=worms, columns = single frames. 

BinVec1    = (Edges1(1:end-1)+Edges1(2:end))/2;
BinVec2    = (Edges2(1:end-1)+Edges2(2:end))/2;
NumOfBins1 = length(BinVec1);
NumOfBins2 = length(BinVec2);

if exist('EdgesForValuesHistograms','var')
    BinVec_values    = (EdgesForValuesHistograms(1:end-1)+EdgesForValuesHistograms(2:end))/2;
    NumOfBins_values = length(BinVec_values);
    Histograms       = zeros(NumOfBins1,NumOfBins2,NumOfBins_values,'single')*NaN;
else
    BinVec_values    = [];
    Histograms       = [];    
end

Vec1     = Mat1(:);
Vec2     = Mat2(:);
ValueVec = ValueMat(:);

Indices_Mat1 = discretize(Vec1,Edges1);
Indices_Mat2 = discretize(Vec2,Edges2);

AverageValuesMat = zeros(NumOfBins1,NumOfBins2,'single')*NaN;
for ind1 = 1:NumOfBins1
    Current_Edge1_ind = Indices_Mat1==ind1;
    for ind2 = 1:NumOfBins2
         Current_Edge2_ind      = Indices_Mat2==ind2;  
         ValuesIn_Bin1_and_Bin2 = ValueVec(Current_Edge1_ind & Current_Edge2_ind);
         MeanValue              = nanmean(ValuesIn_Bin1_and_Bin2);
         AverageValuesMat(ind1, ind2) = MeanValue;

         if exist('EdgesForValuesHistograms','var')
             N = histc(ValuesIn_Bin1_and_Bin2,EdgesForValuesHistograms);    
             N = N(1:end-1);
             N = N/sum(N);
             Histograms(ind1,ind2,:) = N;
         end
    end
end                       
 
return

function [AverageValuesMat, BinVec1, BinVec2, Histograms, BinVec_values] = ComputeFractionAboveValueWithin2DHistogram_inline ...
              (Mat1, Mat2, ValueMat, Edges1, Edges2, EdgesForValuesHistograms, CutOffValue)

% Mat1 is two dimensional: rows=worms, columns = single frames. 
% Mat2 is two dimensional: rows=worms, columns = single frames. 

BinVec1    = (Edges1(1:end-1)+Edges1(2:end))/2;
BinVec2    = (Edges2(1:end-1)+Edges2(2:end))/2;
NumOfBins1 = length(BinVec1);
NumOfBins2 = length(BinVec2);

if exist('EdgesForValuesHistograms','var')
    BinVec_values    = (EdgesForValuesHistograms(1:end-1)+EdgesForValuesHistograms(2:end))/2;
    NumOfBins_values = length(BinVec_values);
    Histograms       = zeros(NumOfBins1,NumOfBins2,NumOfBins_values,'single')*NaN;
else
    BinVec_values    = [];
    Histograms       = [];    
end

Vec1     = Mat1(:);
Vec2     = Mat2(:);
ValueVec = ValueMat(:);

Indices_Mat1 = discretize(Vec1,Edges1);
Indices_Mat2 = discretize(Vec2,Edges2);

AverageValuesMat = zeros(NumOfBins1,NumOfBins2,'single')*NaN;
for ind1 = 1:NumOfBins1
    Current_Edge1_ind = Indices_Mat1==ind1;
    for ind2 = 1:NumOfBins2
         Current_Edge2_ind      = Indices_Mat2==ind2;  
         ValuesIn_Bin1_and_Bin2 = ValueVec(Current_Edge1_ind & Current_Edge2_ind);
         Fraction              = length(find(ValuesIn_Bin1_and_Bin2>CutOffValue)) / length(find(~isnan(ValuesIn_Bin1_and_Bin2)));
         AverageValuesMat(ind1, ind2) = Fraction;

         if exist('EdgesForValuesHistograms','var')
             N = histc(ValuesIn_Bin1_and_Bin2,EdgesForValuesHistograms);    
             N = N(1:end-1);
             N = N/sum(N);
             Histograms(ind1,ind2,:) = N;
         end
    end
end                       
 
return

function StructureOut = ComputeAndPlot2DActivationProbability(Edges1, Edges2, C, C_minus_Threshold_over_C, ...
    Downgradient_NoActivation_ActivityStartsLow, ActivationInitiation_NotOOB, plotme)
BlueToGreenColormap      = zeros(64,3); 
for ind=1:64
    BlueToGreenColormap(ind,2) = (ind-1)/(64-1);     
    BlueToGreenColormap(ind,3) = (64-ind)/(64-1);     
end

XTickVec = Edges1;
YTickVec = Edges2;
normalization  = 0; % Total Number of events

Mat1           = C(Downgradient_NoActivation_ActivityStartsLow);
Mat2           = C_minus_Threshold_over_C(Downgradient_NoActivation_ActivityStartsLow);
N_NotActivated = Compute2DHistogram_inline (Mat1, Mat2, Edges1, Edges2, normalization);
% figure; imagesc(N_NotActivated'); colorbar; set(gca,'xticklabel', XTICKLABELS, 'yticklabel', YTICKLABELS);

Mat1           = C(ActivationInitiation_NotOOB);
Mat2           = C_minus_Threshold_over_C(ActivationInitiation_NotOOB);
N_Activated    = Compute2DHistogram_inline (Mat1, Mat2, Edges1, Edges2, normalization);
% figure; imagesc(N_Activated'); colorbar; set(gca,'xticklabel', XTICKLABELS, 'yticklabel', YTICKLABELS);

ActivationProbability = N_Activated ./ (N_Activated + N_NotActivated);

% The Probability of NOT being activated within t = 1 second=30 frames, assuming a Poisson process   
% NotActivatedPerSecond           = exp(-ActivationProbability*t);
NotActivatedPerFrame            = 1-ActivationProbability;          
t                               = 30; % frames
NotActivatedPerSecond           = (NotActivatedPerFrame).^t;
ActivationProbability_PerSecond = 1 - NotActivatedPerSecond;
t                               = 30*5; % frames
NotActivatedPerSecond           = (NotActivatedPerFrame).^t;
ActivationProbability_Per5Second = 1 - NotActivatedPerSecond;

if plotme    
%     figure('name','Number of frames WITHOUT or WITH activation','position',get(0,'ScreenSize')) ; 
%     subplot(1,2,1); 
%         log10_N_NotActivated = N_NotActivated;
%         log10_N_NotActivated(log10_N_NotActivated==0)=1;
%         log10_N_NotActivated = log10(log10_N_NotActivated);
%         imagesc(log10_N_NotActivated'); colorbar; 
%         set(gca,'Xtick',(1:length(XTickVec))-0.5,'xticklabel',XTickVec,'Ytick',(1:length(YTickVec))-0.5,'Yticklabel',YTickVec)
%         axis xy; axis image; 
% %         title('Not activated')
%         title('Not activated, log10(#)');
%         ylabel('no normalization');
%     subplot(1,2,2); 
%         imagesc(N_Activated'); colorbar; 
%         set(gca,'Xtick',(1:length(XTickVec))-0.5,'xticklabel',XTickVec,'Ytick',(1:length(YTickVec))-0.5,'Yticklabel',YTickVec)
%         axis xy; axis image; title('Activated')
    %    
    figure('name','Activation probability','position',get(0,'ScreenSize')) ; 
    subplot(1,3,1); 
        imagesc_withoutNaNs(1:size(ActivationProbability,1), 1:size(ActivationProbability,2), ActivationProbability', ...
                            [0 0.12], BlueToGreenColormap, [1 1 1]);  colorbar; 
        set(gca,'Xtick',(1:length(XTickVec))-0.5,'xticklabel',XTickVec,'Ytick',(1:length(YTickVec))-0.5,'Yticklabel',YTickVec)
        axis xy; axis image; title('Probability of activation per frame')
    subplot(1,3,2); 
        imagesc_withoutNaNs(1:size(ActivationProbability_PerSecond,1), 1:size(ActivationProbability_PerSecond,2), ActivationProbability_PerSecond', ...
                            [0 1], BlueToGreenColormap, [1 1 1]);  colorbar; 
        set(gca,'Xtick',(1:length(XTickVec))-0.5,'xticklabel',XTickVec,'Ytick',(1:length(YTickVec))-0.5,'Yticklabel',YTickVec)
        axis xy; axis image; title('Probability of activation per second')
    subplot(1,3,3); 
        imagesc_withoutNaNs(1:size(ActivationProbability_Per5Second,1), 1:size(ActivationProbability_Per5Second,2), ActivationProbability_Per5Second', ...
                            [0 1], BlueToGreenColormap, [1 1 1]);  colorbar; 
        set(gca,'Xtick',(1:length(XTickVec))-0.5,'xticklabel',XTickVec,'Ytick',(1:length(YTickVec))-0.5,'Yticklabel',YTickVec)
        axis xy; axis image; title('Probability of activation per 5 seconds')
end

% assign to structure
StructureOut.N_NotActivated                  = N_NotActivated;
StructureOut.N_Activated                     = N_Activated;
StructureOut.ActivationProbabilityPerFrame   = ActivationProbability;
StructureOut.ActivationProbability_PerSecond = ActivationProbability_PerSecond;

return

function StructureOut = Align_Concentration_Activity_Behavior_By_EventMatrix(Data, EventLogicalNatrix, FramesWindow)
NumberOfFrames = size(Data.NaN,2);
NumberOfWorms  = size(Data.NaN,1);
NumberOfEvents = length(find(EventLogicalNatrix));

EventTracks          = zeros(NumberOfEvents,1,'single')*NaN;
EventInitialFrames   = zeros(NumberOfEvents,1,'single')*NaN;
EventConcentration   = zeros(NumberOfEvents,1,'single')*NaN;
EventThreshold       = zeros(NumberOfEvents,1,'single')*NaN;
EventPastDeltaFOverF = zeros(NumberOfEvents,1,'single')*NaN;
EventPastNeuronValue = zeros(NumberOfEvents,1,'single')*NaN;
EventBehavior        = zeros(NumberOfEvents,1,'single')*NaN;
Concentration        = zeros(NumberOfEvents,length(FramesWindow),'single')*NaN;
FoldChange           = zeros(NumberOfEvents,length(FramesWindow),'single')*NaN;
DeltaFOverF          = zeros(NumberOfEvents,length(FramesWindow),'single')*NaN;
NeuronValue          = zeros(NumberOfEvents,length(FramesWindow),'single')*NaN;
Behavior             = zeros(NumberOfEvents,length(FramesWindow),'single')*NaN;

event_ind = 0;
for tr= 1:NumberOfWorms
    CurrentEventsFrames = find(EventLogicalNatrix(tr,:));
    for frame = CurrentEventsFrames
       event_ind = event_ind +1;
       EventTracks(event_ind)          = tr;
       EventInitialFrames(event_ind)   = frame;
       EventConcentration(event_ind)   = Data.C(tr,frame);
       EventThreshold(event_ind)       = Data.AdaptiveThreshold(tr,frame);
       EventPastDeltaFOverF(event_ind) = nanmean(Data.DeltaFOverFFiltered(tr,(frame-15):frame));
       EventPastNeuronValue(event_ind) = nanmean(Data.NeuronValueFiltered(tr,(frame-15):frame));
       EventBehavior(event_ind)        = Data.BehaviorCode_LowLevel(tr,frame-1);
       
       CurrentFrames = frame + FramesWindow;
       if (CurrentFrames(end)>NumberOfFrames)||(CurrentFrames(1)<1)
           disp('debuging needed...')
           keyboard;
       end       
       Concentration(event_ind,:) = Data.C(tr,CurrentFrames);
       FoldChange(event_ind,:)    = Data.FoldChange(tr,CurrentFrames);
       DeltaFOverF(event_ind,:)   = Data.DeltaFOverFFiltered(tr,CurrentFrames);
       NeuronValue(event_ind,:)   = Data.NeuronValueFiltered(tr,CurrentFrames);
       Behavior(event_ind,:)      = Data.BehaviorCode_LowLevel(tr,CurrentFrames);       
    end    
end

StructureOut.EventTracks          = EventTracks;
StructureOut.EventInitialFrames   = EventInitialFrames;
StructureOut.EventConcentration   = EventConcentration;
StructureOut.EventThreshold       = EventThreshold;
StructureOut.EventPastDeltaFOverF = EventPastDeltaFOverF;
StructureOut.EventPastNeuronValue = EventPastNeuronValue;
StructureOut.EventBehavior        = EventBehavior;

StructureOut.Concentration        = Concentration;
StructureOut.FoldChange           = FoldChange;
StructureOut.DeltaFOverF          = DeltaFOverF;
StructureOut.NeuronValue          = NeuronValue;
StructureOut.Behavior             = Behavior;

StructureOut.Time                 = FramesWindow/30;
return

function PlotBehaviorMatrix(BehaviorMatrix, BehaviorColormapLowLevel, X, Y)
if ~exist('BehaviorColormapLowLevel','var')||isempty(BehaviorColormapLowLevel)
    load('ColormapLowLevelBehavior.mat','BehaviorColormapLowLevel');
end
if exist('X','var')
    imagesc(X,Y,BehaviorMatrix);     
else    
    imagesc(BehaviorMatrix);
end
colorbar; colormap(BehaviorColormapLowLevel)
set(gca,'clim',[0 6]);
return

function BehaviorProbability = PlotAverageBehaviorFromMatrix(BehaviorMatrix, Time)
%% Behavioe code
% OutOfBounds                 = (BehaviorCode_LowLevel_Vec==0);
% ForwardOrCurve              = (BehaviorCode_LowLevel_Vec==1)|(BehaviorCode_LowLevel_Vec==2);
% Pause                       = (BehaviorCode_LowLevel_Vec==3);
% Reverse                     = (BehaviorCode_LowLevel_Vec==4);
% Omega                       = (BehaviorCode_LowLevel_Vec==5);
% SharpTurns                  = (BehaviorCode_LowLevel_Vec==6);

ColorPerBehaviorCode = {'','b','y','r','g'}; % Codes: 1-5
ColorPerBehaviorCode = {'',[0.5 0.5 0.5],'k','r','m'}; % Codes: 1-5
ColorPerBehaviorCode = {'','','k','r','m'}; % Codes: 1-5
BehaviorCodeToInclude_ByOrder = [4 5 3 2 1]; % exclude out of bound, sharp turns are with pauses, and don't show the forward or curve
%% Calculate BehaviorProbability
TotalNumberOfInBoundTracks = sum(BehaviorMatrix>0,1);
BehaviorProbability        = zeros(6,size(BehaviorMatrix,2));
for b_ind = 1:6
    CurrentCodeVec = sum(BehaviorMatrix==b_ind,1);
    BehaviorProbability(b_ind,:) = CurrentCodeVec ./ TotalNumberOfInBoundTracks;
end
%% Add sharp turns to pauses for probabillity plots
BehaviorProbabilityForDisplay = BehaviorProbability(1:5,:);
BehaviorProbabilityForDisplay(3,:)= BehaviorProbability(3,:)+ BehaviorProbability(6,:); 

%% Plot 
BaseLine = zeros(1,size(BehaviorMatrix,2));
for b_ind = 1:length(BehaviorCodeToInclude_ByOrder)
    CurrentCode          = BehaviorCodeToInclude_ByOrder(b_ind);
    CurrentProbabilities = BehaviorProbabilityForDisplay(CurrentCode,:);
    CurrentColor         = ColorPerBehaviorCode{CurrentCode};
    if isempty(CurrentColor)||isempty(find(CurrentProbabilities>eps,1))
        continue
    end
    lower      = BaseLine;
    upper      = BaseLine + CurrentProbabilities;
    jbfill(Time,upper,lower,CurrentColor,'k',1,1); hold on;
    BaseLine = upper;
end


return

function [Concentration_sorted, FoldChange_sorted, Behavior_sorted, RelativeDeltaFOverF_sorted] = ...
    SortMatricesAndPlot(Concentration, FoldChange, Behavior, RelativeDeltaFOverF, SortingIndices, ...
    plotme, FigureTitle, Time, BehaviorColormapLowLevel, XLIMITS, NeuronClimit)

Concentration_sorted        = Concentration(SortingIndices,:);
FoldChange_sorted           = FoldChange(SortingIndices,:);
Behavior_sorted             = Behavior(SortingIndices,:);
RelativeDeltaFOverF_sorted  = RelativeDeltaFOverF(SortingIndices,:);
EventIndices                = 1:length(SortingIndices);
if plotme
    figure('name',[FigureTitle,'. log10(C)']);         
        imagesc(Time,EventIndices,log10(Concentration_sorted));      colorbar; xlim(XLIMITS);
    figure('name',[FigureTitle,'. Fold Change']);      
        imagesc(Time,1:size(FoldChange_sorted,1),FoldChange_sorted); colorbar; set(gca,'clim',[-1 1]); xlim(XLIMITS);
    figure('name',[FigureTitle,'. Relative DeltaF/F']);
        imagesc(Time,EventIndices,RelativeDeltaFOverF_sorted);       colorbar; set(gca,'clim',NeuronClimit); xlim(XLIMITS);
    figure('name',[FigureTitle,'. Behavior']); 
        PlotBehaviorMatrix(Behavior_sorted, BehaviorColormapLowLevel, Time, EventIndices); xlim(XLIMITS);
    figure('name',[FigureTitle,'. Behavior probability']); 
        PlotAverageBehaviorFromMatrix(Behavior_sorted, Time); xlim(XLIMITS); ylim([0 1]);           
end

return

function [SortByReachThresholdValue, SortedTimes, OriginalReachThresholdFrames] = FindReachThresholdIndices(DeltaFOverFMatrix, ZeroIndex, ThresholdValue)
NumOfTracks          = size(DeltaFOverFMatrix,1);
ReachThresholdFrames = ones(1,NumOfTracks)*NaN;
for row = 1:NumOfTracks
    EndIndex   = find(DeltaFOverFMatrix(row,ZeroIndex:end)>=ThresholdValue,1,'first')-1;
    if isempty(EndIndex)
        continue
    end
    ReachThresholdFrames(row)  = EndIndex;          
end
% figure; hist(ReachThresholdFrames,10);
OriginalReachThresholdFrames             = ReachThresholdFrames;
[SortedTimes, SortByReachThresholdValue] = sort(ReachThresholdFrames);

return

function [SortByBehaviorTermination, SortedTimes, OriginalTerminationFrames, OriginalTerminationBehavior] = FindBehaviorTerminationIndices(BehaviorMatrix, ZeroIndex)
InitialBehaviorCodes = BehaviorMatrix(:,ZeroIndex-1);
if ~isempty(find(diff(InitialBehaviorCodes)~=0,1))
    disp('Not all rows start with the same behavior code')
end
NumOfTracks         = size(BehaviorMatrix,1);
TerminationFrames   = ones(1,NumOfTracks)*NaN;
TerminationBehavior = zeros(1,NumOfTracks)*NaN;
for row = 1:NumOfTracks
    LogicalVec = BehaviorMatrix(row,:) == InitialBehaviorCodes(row);
    EndIndex   = find(LogicalVec(ZeroIndex:end)==0,1,'first')-1;
    if isempty(EndIndex)
        continue
    end
    TerminationFrames(row)  = EndIndex;          
    TerminationBehavior(row)= BehaviorMatrix(row, ZeroIndex+EndIndex);          
end
% figure; hist(TerminationFrames,10);
OriginalTerminationFrames                = TerminationFrames;
OriginalTerminationBehavior              = TerminationBehavior;
[SortedTimes, SortByBehaviorTermination] = sort(TerminationFrames);

return

function [SortByBehaviorInitiation, SortedTimes]  = FindBehaviorInitiationIndices(BehaviorMatrix, ZeroIndex, BehaviorCodes)

LogicalMatrix = BehaviorMatrix==BehaviorCodes(1);
if length(BehaviorCodes)>1
    for code_ind = 2:length(BehaviorCodes)
        LogicalMatrix = LogicalMatrix | (BehaviorMatrix==BehaviorCodes(code_ind));
    end    
end
InitialBehaviorCodes = BehaviorMatrix(:,ZeroIndex-1);
if ~isempty(find(diff(InitialBehaviorCodes)~=0,1))
    disp('Not all rows start with the same behavior code')
end
NumOfTracks = size(BehaviorMatrix,1);
StartFrames = ones(1,NumOfTracks)*NaN;
for row = 1:NumOfTracks
    LogicalVec = LogicalMatrix(row,:);
    StartIndex = find(LogicalVec(ZeroIndex:end)==1,1,'first')-1;
    if isempty(StartIndex)
        continue
    end
    StartFrames(row)=StartIndex;          
end
% figure; hist(TerminationFrames,10);
[SortedTimes, SortByBehaviorInitiation] = sort(StartFrames);

return

function Djs = JehnsenShannonDivergence_inline (P,Q)
% P and Q are distribution vectors of the same size
% Djs is a scalar

% avoid numerical problems
P(P==0)=eps;
Q(Q==0)=eps;

M   = (P+Q)/2;
Djs = sum(1/2* (P.*log2(P./M)) + 1/2* (Q.*log2(Q./M)));

return

function Stats = BarsAndStats(binranges, Vec, Stats, ind)

bincenters = (binranges(1:end-1)+binranges(2:end))/2;
N = histcounts(Vec, binranges);    
N = N/sum(N); 
[~,p]          = ttest(Vec);  
Skewness       = skewness(Vec,0); 
ExcessKurtosis = kurtosis(Vec,0)-3; 
Stats.Mean(ind)= mean(Vec); 
Stats.Std(ind) = std(Vec); 
Stats.SEM(ind) = std(Vec)/sqrt(length(Vec)); 
Stats.N(ind,:) = N; 
Stats.p(ind)   = p; 
Stats.Skewness(ind)       = Skewness; 
Stats.ExcessKurtosis(ind) = ExcessKurtosis; 
subplot(4,1,ind); bar(bincenters,N); title(['mean=',num2str(mean(Vec),3),' , std=',num2str(std(Vec),3),' , p=',num2str(p,2),' , ExKu=',num2str(ExcessKurtosis,2)]);     

return

function Plot_Environment_Activity_Behavior(Data,tr, Xlimits, NeuronYLimit, PlotStimulusInsteadOfConcentration, ConcentrationLimitsIN, StartXaxisFromZero)

if ~exist('StartXaxisFromZero','var') || isempty(StartXaxisFromZero)|| isempty(Xlimits)
    StartXaxisFromZero = false;  % False: time=0 at the full gradient initiation. True: time=0 at Xlimits(1)
end

ScalingFactorDilutionToMicroM = 11.16e6;
PlotConcentrationInMicroM = true;
PlotHeadAngle    = false;
NumberOfFrames   = size(Data.NaN,2);

FirstFrameOfOdor                = Data.ExpInfo.GradientStartFrame(tr); % Gradient STARTS to be established at time 0  
FirstFrameOfEstablishedGradient = FirstFrameOfOdor + Data.ExpInfo.GradientRiseTime_InitiationTo90Percent_InSec(tr)*30; % Gradient is established   
LastFrameOfEstablishedGradient  = Data.ExpInfo.GradientEndFrame(tr);   % Gradient STARTS to vanish at this frame 
LastFrameOfEstablishedGradient  = min([LastFrameOfEstablishedGradient NumberOfFrames]);

GradientFrames                                                                 = false(1,NumberOfFrames);
GradientFrames(FirstFrameOfEstablishedGradient:LastFrameOfEstablishedGradient) = true;

if ~exist('PlotStimulusInsteadOfConcentration','var') || isempty(PlotStimulusInsteadOfConcentration)
    PlotStimulusInsteadOfConcentration = false;
end

FramesVec = (1:size(Data.NaN,2))-FirstFrameOfEstablishedGradient;
% Time = FramesVec;         xlabelstr = 'Frames';
Time = FramesVec/30;      xlabelstr = 'Time [sec]';
if ~exist('Xlimits','var') || isempty(Xlimits)
%     Xlimits = [Time(1) Time(end)];
    Xlimits = [0 Time(end)]; 
end
if StartXaxisFromZero
    ReferenceTime = Xlimits(1);
    Time    = Time   - ReferenceTime;
    Xlimits = Xlimits- ReferenceTime;
    FrameRate = 30; % Hz
    disp(['t=0 corresponds to frame ',num2str(FirstFrameOfEstablishedGradient + FrameRate*ReferenceTime)])
else    
    disp(['t=0 corresponds to frame ',num2str(FirstFrameOfEstablishedGradient)]);
end
disp(['Movie name: ',Data.ExpInfo.MovieNames{Data.ExpInfo.MovieNumber(tr)}, ' , worm number ',num2str(Data.ExpInfo.TrackNumber(tr))]);

Stimulus       = Data.StimulusFiltered(tr,:);
Concentration  = Data.C(tr,:);
Threshold      = Data.AdaptiveThreshold(tr,:);
NeuralActivity = Data.DeltaFOverFFiltered(tr,:);  % use unfiltered data 

ActivationInitiationFrames = find(Data.ActivationStartFrames(tr,:));  
ActivationPeakFrames       = find(Data.ActivationPeakFrames(tr,:));  
ActivationOOB              = Data.ActivationFrameOOB(tr,ActivationInitiationFrames);  
SharpActivation            = Data.ActivationSharp(tr,ActivationInitiationFrames);  

NotOOB_Sharp    = ~ActivationOOB &  SharpActivation;
NotOOB_NotSharp = ~ActivationOOB & ~SharpActivation;

% Low level Behavior:
BehaviorCode_LowLevel_Vec   = Data.BehaviorCode_LowLevel(tr,:);
OutOfBounds                 = (BehaviorCode_LowLevel_Vec==0);
ForwardOrCurve              = (BehaviorCode_LowLevel_Vec==1)|(BehaviorCode_LowLevel_Vec==2);
Pause                       = (BehaviorCode_LowLevel_Vec==3);
Reverse                     = (BehaviorCode_LowLevel_Vec==4);
Omega                       = (BehaviorCode_LowLevel_Vec==5);
SharpTurns                  = (BehaviorCode_LowLevel_Vec==6);

Stimulus_InGradient                       = Stimulus;        
Stimulus_InGradient(~GradientFrames)      = NaN;
Stimulus_InGradient(OutOfBounds)          = NaN;
Concentration_InGradient                  = Concentration;   
Concentration_InGradient(~GradientFrames) = NaN;
Concentration_InGradient(OutOfBounds)     = NaN;

maxNeuralActivity   = max(NeuralActivity((Time>=Xlimits(1))&(Time<=Xlimits(2))));
if exist('NeuronYLimit','var')&&~isempty(NeuronYLimit)
    maxNeuralActivity = NeuronYLimit;
end
ConcentrationLimits = [Data.ExpInfo.ConcentrationAtArenaBottom(tr)   Data.ExpInfo.ConcentrationAtArenaTop(tr)];
if exist('ConcentrationLimitsIN','var')&&~isempty(ConcentrationLimitsIN)
    ConcentrationLimits = ConcentrationLimitsIN;
end

%% Relative Head Direction
if PlotHeadAngle
    NumOfSubplots = 5; 
    HeadAngleRelativeToDirection         = Data.HeadAngleRelativeToDirection(tr,:);
    HeadAngleRelativeToDirection_Forward = HeadAngleRelativeToDirection;
    HeadAngleRelativeToDirection_Reverse = HeadAngleRelativeToDirection;
    HeadAngleRelativeToDirection_Pause   = HeadAngleRelativeToDirection;
    HeadAngleRelativeToDirection_Forward(~ForwardOrCurve)= NaN;
    HeadAngleRelativeToDirection_Reverse(~Reverse)       = NaN;
    HeadAngleRelativeToDirection_Pause(~Pause)           = NaN;
else
    NumOfSubplots = 4;    
end

%% Plot
figure('name',['Tracks number ' , num2str(tr)],'position',get(0,'ScreenSize')); 
NumOfSubplots = NumOfSubplots-1;

EndAversiveBehaviorIndex   = find(diff(ForwardOrCurve)==1);
StartAversiveBehaviorIndex = find(diff(ForwardOrCurve)==-1)+1;
EndAversiveBehaviorTime    = Time(EndAversiveBehaviorIndex);
StartAversiveBehaviorTime  = Time(StartAversiveBehaviorIndex);
EndAversiveBehaviorTime    = EndAversiveBehaviorTime(   (EndAversiveBehaviorTime  >Xlimits(1)) & (EndAversiveBehaviorTime  <Xlimits(2))   );
StartAversiveBehaviorTime  = StartAversiveBehaviorTime( (StartAversiveBehaviorTime>Xlimits(1)) & (StartAversiveBehaviorTime<Xlimits(2)) );

subplot(NumOfSubplots,1,1); 
if PlotStimulusInsteadOfConcentration
    plot(Time, Stimulus,'-','color',[0.5 0.5 0.5]); hold on;               
    plot(Time, Stimulus_InGradient,'k-'); hold on;      
        
    plot(Time(ActivationInitiationFrames(ActivationOOB)),   Stimulus(ActivationInitiationFrames(ActivationOOB)),'rx'); hold on; 
    plot(Time(ActivationPeakFrames(ActivationOOB)),         Stimulus(ActivationPeakFrames(ActivationOOB)),'gx'); hold on; 
    plot(Time(ActivationInitiationFrames(NotOOB_NotSharp)), Stimulus(ActivationInitiationFrames(NotOOB_NotSharp)),'ro','markerfacecolor','r'); hold on; 
    plot(Time(ActivationPeakFrames(NotOOB_NotSharp)),       Stimulus(ActivationPeakFrames(NotOOB_NotSharp)),'go','markerfacecolor','g'); hold on; 
    plot(Time(ActivationInitiationFrames(NotOOB_Sharp)),    Stimulus(ActivationInitiationFrames(NotOOB_Sharp)),'r^','markerfacecolor','r'); hold on; 
    plot(Time(ActivationPeakFrames(NotOOB_Sharp)),          Stimulus(ActivationPeakFrames(NotOOB_Sharp)),'g^','markerfacecolor','g'); hold on; 
    
    xlim(Xlimits); ylim([-1 1]); ylabel('Gradient axis');
%     legend('Stimulus, OOB','Stimulus');
    
else
    if PlotConcentrationInMicroM
        Concentration            = ScalingFactorDilutionToMicroM * Concentration;
        Concentration_InGradient = ScalingFactorDilutionToMicroM * Concentration_InGradient;
        Threshold                = ScalingFactorDilutionToMicroM * Threshold;
        ConcentrationLimits      = ScalingFactorDilutionToMicroM * ConcentrationLimits;
    end
    plot(Time, Concentration,'-','color',[0.5 0.5 0.5]);  hold on;
    plot(Time, Concentration_InGradient,'k-');  hold on;
    plot(Time, Threshold,'-','color',[0.5 0 0 ]);  hold on;    
    
    plot(Time(ActivationInitiationFrames(ActivationOOB)),   Concentration(ActivationInitiationFrames(ActivationOOB)),'rx'); hold on; 
    plot(Time(ActivationPeakFrames(ActivationOOB)),         Concentration(ActivationPeakFrames(ActivationOOB)),'gx'); hold on; 
    plot(Time(ActivationInitiationFrames(NotOOB_NotSharp)), Concentration(ActivationInitiationFrames(NotOOB_NotSharp)),'ro','markerfacecolor','r'); hold on; 
    plot(Time(ActivationPeakFrames(NotOOB_NotSharp)),       Concentration(ActivationPeakFrames(NotOOB_NotSharp)),'go','markerfacecolor','g'); hold on; 
    plot(Time(ActivationInitiationFrames(NotOOB_Sharp)),    Concentration(ActivationInitiationFrames(NotOOB_Sharp)),'r^','markerfacecolor','r'); hold on; 
    plot(Time(ActivationPeakFrames(NotOOB_Sharp)),          Concentration(ActivationPeakFrames(NotOOB_Sharp)),'g^','markerfacecolor','g'); hold on; 

    xlim(Xlimits); ylim(ConcentrationLimits); ylabel('Concentration'); 
%     legend('Concentration, OOB','Concentration','Adaptive Threshold');
    
    line(ones(1,2)*Time(FirstFrameOfEstablishedGradient),get(gca,'ylim'),'color','b'); 
    line(ones(1,2)*Time(LastFrameOfEstablishedGradient), get(gca,'ylim'),'color','b'); 
    
    for line_ind=1:length(EndAversiveBehaviorTime)
        CurrentTime = EndAversiveBehaviorTime(line_ind);
        line(ones(1,2)*CurrentTime,get(gca,'ylim'),'color',[0.5 0.5 0.5]);         
    end    
    for line_ind=1:length(StartAversiveBehaviorTime)
        CurrentTime = StartAversiveBehaviorTime(line_ind);
        line(ones(1,2)*CurrentTime,get(gca,'ylim'),'color','r');         
    end
end
    
subplot(NumOfSubplots,1,2); 
    plot(Time, NeuralActivity,'k-');     hold on;           
    plot(Time(ActivationInitiationFrames(ActivationOOB)),   NeuralActivity(ActivationInitiationFrames(ActivationOOB)),'rx'); hold on; 
    plot(Time(ActivationPeakFrames(ActivationOOB)),         NeuralActivity(ActivationPeakFrames(ActivationOOB)),'gx'); hold on; 
    plot(Time(ActivationInitiationFrames(NotOOB_NotSharp)), NeuralActivity(ActivationInitiationFrames(NotOOB_NotSharp)),'ro','markerfacecolor','r'); hold on; 
    plot(Time(ActivationPeakFrames(NotOOB_NotSharp)),       NeuralActivity(ActivationPeakFrames(NotOOB_NotSharp)),'go','markerfacecolor','g'); hold on; 
    plot(Time(ActivationInitiationFrames(NotOOB_Sharp)),    NeuralActivity(ActivationInitiationFrames(NotOOB_Sharp)),'r^','markerfacecolor','r'); hold on; 
    plot(Time(ActivationPeakFrames(NotOOB_Sharp)),          NeuralActivity(ActivationPeakFrames(NotOOB_Sharp)),'g^','markerfacecolor','g'); hold on; 
    xlim(Xlimits); ylim([-0.2 maxNeuralActivity]); ylabel('\DeltaF/F');

    line(ones(1,2)*Time(FirstFrameOfEstablishedGradient),get(gca,'ylim'),'color','b'); 
    line(ones(1,2)*Time(LastFrameOfEstablishedGradient), get(gca,'ylim'),'color','b'); 
    for line_ind=1:length(EndAversiveBehaviorTime)
        CurrentTime = EndAversiveBehaviorTime(line_ind);
        line(ones(1,2)*CurrentTime,get(gca,'ylim'),'color',[0.5 0.5 0.5]);         
    end    
    for line_ind=1:length(StartAversiveBehaviorTime)
        CurrentTime = StartAversiveBehaviorTime(line_ind);
        line(ones(1,2)*CurrentTime,get(gca,'ylim'),'color','r');         
    end
    
subplot(NumOfSubplots,1,3);
    %% low level behavior
    C_Pause          = ones(length(Concentration),1); C_Pause(~Pause)=NaN;
    C_SharpTurns     = ones(length(Concentration),1); C_SharpTurns(~SharpTurns)=NaN;
    C_ForwardOrCurve = ones(length(Concentration),1); C_ForwardOrCurve(~ForwardOrCurve)=NaN;
    C_Reverse        = ones(length(Concentration),1); C_Reverse(~Reverse)=NaN;
    C_Omega          = ones(length(Concentration),1); C_Omega(~Omega)=NaN;
    C_OutOfBounds    = ones(length(Concentration),1); C_OutOfBounds(~OutOfBounds)=NaN;

    hold on; 
    plot(Time, C_ForwardOrCurve, '-','color',[0.58 0.58 0.58]); hold on;  
    plot(Time, C_OutOfBounds ,   'y-','linewidth',2); 
    plot(Time, C_SharpTurns,     'k*','markerfacecolor','k','markersize',5); 
    plot(Time, C_Pause ,         'k-','linewidth',2);                              
    plot(Time, C_Reverse,        'r-','linewidth',2);                   
    plot(Time, C_Omega,          'm-','linewidth',2); 
    xlim(Xlimits); ylabel('Behavior');
    ylim([0.9 1.1]);
    xlabel(xlabelstr)
%     legend('Forward','OOB','SharpTurn','Pause','Reverse','Omega');
    line(ones(1,2)*Time(FirstFrameOfEstablishedGradient),get(gca,'ylim'),'color','b'); 
    line(ones(1,2)*Time(LastFrameOfEstablishedGradient), get(gca,'ylim'),'color','b'); 
    
if PlotHeadAngle
    
subplot(NumOfSubplots,1,4); 
%     plot(Time, StimulusCentroid,'-','color',[0.5 0.5 0.5]); hold on; 
    plot(Time, HeadAngleRelativeToDirection_Pause,'y.'); hold on; 
    plot(Time, HeadAngleRelativeToDirection_Forward,'b.'); hold on; 
    plot(Time, -HeadAngleRelativeToDirection_Reverse,'r.'); hold on; 
    xlim(Xlimits); ylim([-90 90]); ylabel('Relative head angle');
    line(ones(1,2)*Time(FirstFrameOfEstablishedGradient),get(gca,'ylim'),'color','b'); 
    line(ones(1,2)*Time(LastFrameOfEstablishedGradient), get(gca,'ylim'),'color','b'); 
end   


return

function GreenColormap = CreateWhiteToGreenColormap_inline(SetLowestValueColor, LowestValueColor_RGBvec)
if ~exist('SetLowestValueColor','var')
    SetLowestValueColor = false;
end
%% black to green
% GreenColormap = zeros(64,3);          % initialization
% if SetLowestValueColor
%     GreenColormap(2:end,2) = linspace(0,1,63);
%     GreenColormap(1,:)     = LowestValueColor_RGBvec;
% else
%     GreenColormap(:,2)     = linspace(0,1,64);
% end

% %% white to green
% GreenColormap = ones(64,3);          % initialization
% if SetLowestValueColor
%     GreenColormap(2:end,1) = linspace(1,0,63);
%     GreenColormap(2:end,3) = linspace(1,0,63);
%     GreenColormap(1,:)     = LowestValueColor_RGBvec;
% else
%     GreenColormap(:,1)     = linspace(1,0,64);
%     GreenColormap(:,3)     = linspace(1,0,64);
% end

%% white to dark green
GreenColormap = ones(64,3);          % initialization
if SetLowestValueColor
    GreenColormap(2:end,1) = linspace(1,0,63);
    GreenColormap(2:end,3) = linspace(1,0,63);
    GreenColormap(2:end,2) = linspace(1,0.5,63);
    GreenColormap(1,:)     = LowestValueColor_RGBvec;
else
    GreenColormap(:,1)     = linspace(1,0,64);
    GreenColormap(:,3)     = linspace(1,0,64);
    GreenColormap(:,2)     = linspace(1,0.5,64);
end


return



