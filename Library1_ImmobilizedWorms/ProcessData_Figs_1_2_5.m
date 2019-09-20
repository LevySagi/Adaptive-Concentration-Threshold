function ProcessData_Figs_1_2_5
%
%  Written by Sagi Levy, September 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%  Input:  'Data_AWConImaging_Raw.mat'
%          Raw data for imaging of different compartments of AWC(ON) in animals immobilized in microfluidic devices     
%          Related to Figures 1, 2, 5, S1, S2, S3 and S5
%  Output: 'Data_AWConImaging_Processed.mat'
%               Processing of neuronal activity, 
%               error estimation
%               evaluation of all models' parameters, including ACT model and linear kernel for LN models  
%               Extract models' predictions  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%  Load raw data
load ('H:\Data_AWConImaging_Raw.mat','Data'); 

%%  Basic processing: Compute DeltaF/F and concentration dynamics features, smoothing, and basic statistics                 
Data = Process_Data_AWConImaging(Data);

%%  Find error estimation by bootstrap and evaluate ACT model parameters   OR  load intermediate processed data .mat file
%%% Threshold constant calculation
NumOfIterations     = 1750;   % Bootstrap is slow. Run parallel sessions to increase speed. See details in the function below  
[Data, ModelParams] = Calculate_Thresholds_And_K_BootstrapForErrorbars(Data, NumOfIterations); 
%%% Adaptation time calculation
CalculateACTModelAdaptationTimeAndPredictions_WithBootstrap(Data, ModelParams);  % Compute adaptation times with bootstrapping. Load output at: 'ProcessedData_ACT_ModelParams.mat'
%%% load intermediate processed data, until this point in the code.
% load('H:\ProcessedData_ACT_ModelParams.mat','Data','ModelParams','ACT_CX17256','ACT_CX17255','ACT_Egl4','ACT_Egl4_ad450','ACT_VaryInitialC'); 
 
%% LN model 
Data = ComputeImpulseFunctions(Data);
Data = PredictResponse_LNmodel(Data);     
Data = OptimizeRectifierForLNmodel(Data);

%% Compute all models' predicted responses to slow decrease of odor concentration, including responsiveness and latency. Compile also actual data statistics
CurrentStrain_Slow = 'str2_CX17256_FromD6_Slow';
CurrentStrain_Fast = 'str2_CX17256_FromD6_Fast';
[ModelPredictions_VaryFinalC,   MeasuredStats_VaryFinalC]   = ComputeAlternativeModelParams (Data, CurrentStrain_Slow, CurrentStrain_Fast);
CurrentStrain_Slow = 'str2_CX17256_ToBuffer_Slow';
BestFit            = ModelPredictions_VaryFinalC.Information.BestFit;  % Predictions for 'str2_CX17256_ToBuffer_Slow' using data from 'str2_CX17256_FromD6_Slow'
[ModelPredictions_VaryInitialC, MeasuredStats_VaryInitialC] = ComputeAlternativeModelParams (Data, CurrentStrain_Slow, CurrentStrain_Fast, BestFit);

%% save output structure
save('H:\Data_AWConImaging_Processed.mat',...
             'Data','ModelParams',...
             'ACT_CX17256','ACT_CX17255','ACT_Egl4','ACT_Egl4_ad450','ACT_VaryInitialC',...
             'ModelPredictions_VaryFinalC',  'MeasuredStats_VaryFinalC',...
             'ModelPredictions_VaryInitialC','MeasuredStats_VaryInitialC','-v7.3');  % save processed data

return


function readme
%% Most names of variables and fields are self explanatory. Read this part for more detailed descriptions

%% RawData 
%      This structure contains raw GCaMP5a fluorescence data measured at the soma, dendrite and axon of the AWC(ON) sensory neuron   
%      Data was recorded in immobalized animals and was used in figures 1, 2, and 5, and in their respective supplementary figures   
%  First layer field names  
%                     1st field contains information about the genotype and type of experiment.  
%                     e.g. Fast/Slow transitions, naive/desensitized, change final/initial concentrations   
%     'str2_CX17256'  The strain CX17256 was used for all experiments except of egl-4 experiments (Figure 5) 
%     'str2_CX17255'  The strain CX17255 was used for egl-4 control experiments (Figure 5), since this is the background strain of the egl-4 mutants imaging lines 
%     '_Egl4'         egl-4 (lf) 
%     'Egl4_ad450'    egl-4 (gf) 
%     '_FromDX_Fast'  Fast transitions from Dilution of 10^(-X) to various lower concentrations  
%     '_FromDX_Slow'  Slow transitions from Dilution of 10^(-X) to various lower concentrations  
%     'Desensitized'  Desensitized animals 
%     'ToBuffer_Slow' Slow transitions from various initial concentrations to buffer 
%  Second layer field Names, and their subfields
%     'PulsesInfo'                   Structure containing information about the experimental protocol
%             'PulseFromButanone'    contains a cell array with names of all experiment from this type (including in other genotypes/training conditions)       
%             'FinalConcentrations'  Final   Concentration, in microM. If this is a vector, each value corresponds to the respective index in 'PulseFromButanone'   
%             'InitialConcentration' Initial Concentration, in microM. If this is a vector, each value corresponds to the respective index in 'PulseFromButanone'   
%     Butanone_dX_To_Butanone_dY     Experiment Name. Transition from dilution of 10^(-X) to 10^(-Y). See 'PulsesInfo' for concentration information       
%             'Soma_RawValue_BackgroundSubtructed'      Fluorescence measured in the soma    
%             'Dendrite_RawValue_BackgroundSubtructed'  Fluorescence measured in the dendrite  
%             'Axon_RawValue_BackgroundSubtructed'      Fluorescence measured in the axon   
%             'TimeVec'                                 Time from stimulus initiation [sec]    
%             'FrameVec'                                Frames from stimulus initiation. Experiment was monitored at 10Hz
%        Additional subfields only for slow transitions:  
%             'Dye'                                     Dye test fluorescence values
%             'Dye_TimeVec'                             Dye test time in seconds 
%  
%% ProcessedData 
%  This structure contains raw and processed data.
%  Processing involved running the attached scripts on the 'RawData' structure (above)   
%  First layer fields in the structure are similar to the 'RawData' structure
%  Second layer fields with Experiment names (see 'RawData' structure) has the following additional processed data fields: 
%       'DeltaFOverF'               Measured response. F0 is defined as the fluorescence within 1 second before stimulus initiation    
%       'DeltaFOverF_MAX'           max (DeltaFOverF)    
%       'DeltaFOverF_Smoothed'      Smoothed (lowpass filtered) response
%       'DeltaFOverF_MA_Smoothed'   Smoothed (Moving Average) response
%       'Responsiveness'            [%] responders
%       'Latency'                   Response latency [sec], calculated from the smoothed (lowpass filtered) response 
%       'Latency_MA'                Response latency [sec], calculated from the smoothed (Moving Average) response 
%        Notes:
%           Activity measures were calculated seperately for the soma, dendrites and axon 
%           Other fields are calculation intermediates.
%
%  Additional second layer fields and their subfields
%       'Stats'             This structure gather information for all (second layer) experiments, for each compartment    
%                           'Information'       includes:   'NumberOfRepeats','PulseTypeNames','FinalConcentrations'. 
%                                               Shows only information regarding the specific experiment type (first field layer)   
%                           'DeltaFOverF_MAX'  
%                           'Responsiveness'  
%                           'Latency_InSec'  
%                           'Latency_MA_InSec'  
%       'DyeInformation'    This structure gather information about dye and concentration dynamics properties in slow transition experiments     
%                           'InitialConcentrationPerPulse'  matches the fields in Stats.Information  
%                           'FinalConcentrationPerPulse'    matches the fields in Stats.Information  
%                           'TimeVectorInSec'               for dye measurements 
%                           3D matrices (experiments/repeats/time):  Dye_Raw, Dye_Smoothed, C, C_Smoothed, dCdt, dCdt_Smoothed, d2Cdt2, DeltaCoverC 
%                           'Stats'
%                               2D matrices (experiments/repeats):  dCdt_Peak, dCdt_TimeToPeak, d2Cdt2_Peak, d2Cdt2_TimeToPeak, d2Cdt2_TimeToSignChange, DeltaCoverC_Peak, DeltaCoverC_TimeToPeak 
%                               Structures (experiments/repeats):  
%                                  Computed from lowpass filtered curves: C_AtNeuronActivation,    dCdt_AtNeuronActivation,    d2Cdt2_AtNeuronActivation,    DeltaCoverC_AtNeuronActivation    
%                                  Computed from Moving Average curves:   C_AtNeuronActivation_MA, dCdt_AtNeuronActivation_MA, d2Cdt2_AtNeuronActivation_MA, DeltaCoverC_AtNeuronActivation_MA    
%                                  Each structure contains subfields for each compartment with a 2D matrix (experiments/repeats)  
%     
%   
return

%% Basic Processing of 'RawData' obtained for imaging of immobilized animals
function Data = Process_Data_AWConImaging(Data)

%% Process response curves, concentration features and their corresponding statistics
Data = ComputeDeltaFOverF(Data);                      
Data = SmoothActivity_CalculateResponsesStats(Data);  % Process responses
Data = CalculateStatsOverRepeats (Data);              % Process responses statistics
Data = SmoothDyeAndCalculateStats(Data);              % Process dye curves, dynamic concentration features and their statistics

return

function Data = ComputeDeltaFOverF(Data) 
% Compute DeltaFOverF, Calculate F0 as the mean over one second before pulse initiation  

FieldNames  = fieldnames(Data);
for f_ind = 1:length(FieldNames)
    CurrentField   = FieldNames{f_ind};
    PulseTypeNames = Data.(CurrentField).PulsesInfo.PulseFromButanone;
%     if isfield(ProcessedData.(CurrentField).(PulseTypeNames{1}),'Dye')
%         DyeDataExists = true;
%     else
%         DyeDataExists = false;
%     end        
    
    for pulsetype_ind = 1:length(PulseTypeNames)
        CurrentPulseTypeName = PulseTypeNames{pulsetype_ind}; 
        if ~isfield(Data.(CurrentField),CurrentPulseTypeName)
            continue
        end
        ZeroIndex            = find(Data.(CurrentField).(CurrentPulseTypeName).TimeVec==0);
        F0_indices           = (ZeroIndex-9):ZeroIndex;
        
        % Soma
        SomaRawValue_BackgroundSubtructed = Data.(CurrentField).(CurrentPulseTypeName).Soma_RawValue_BackgroundSubtructed;                
        F0_beginning         = nanmean(SomaRawValue_BackgroundSubtructed(:,F0_indices),2);
        F0_beginning_MAT     = repmat(F0_beginning,1,size(SomaRawValue_BackgroundSubtructed,2));
        SomaDeltaFOverF      = (SomaRawValue_BackgroundSubtructed-F0_beginning_MAT)./F0_beginning_MAT;
        
        Data.(CurrentField).(CurrentPulseTypeName).Soma_DeltaFOverF       = SomaDeltaFOverF;
        
        % Axon and dendrite
        AxonRawValue_BackgroundSubtructed      = Data.(CurrentField).(CurrentPulseTypeName).Axon_RawValue_BackgroundSubtructed;
        DendriteRawValue_BackgroundSubtructed  = Data.(CurrentField).(CurrentPulseTypeName).Dendrite_RawValue_BackgroundSubtructed;
        ProcessActivityWasMonitored            = ~all(isnan(AxonRawValue_BackgroundSubtructed(:)));
        if ~ProcessActivityWasMonitored
            disp([CurrentField,'  ,  ', CurrentPulseTypeName]);
        end
        if ProcessActivityWasMonitored
            F0_beginning         = nanmean(AxonRawValue_BackgroundSubtructed(:,F0_indices),2);
            F0_beginning_MAT     = repmat(F0_beginning,1,size(AxonRawValue_BackgroundSubtructed,2));
            AxonDeltaFOverF      = (AxonRawValue_BackgroundSubtructed-F0_beginning_MAT)./F0_beginning_MAT;
                        
            Data.(CurrentField).(CurrentPulseTypeName).Axon_DeltaFOverF       = AxonDeltaFOverF;
                     
            F0_beginning             = nanmean(DendriteRawValue_BackgroundSubtructed(:,F0_indices),2);
            F0_beginning_MAT         = repmat(F0_beginning,1,size(DendriteRawValue_BackgroundSubtructed,2));
            DendriteDeltaFOverF      = (DendriteRawValue_BackgroundSubtructed-F0_beginning_MAT)./F0_beginning_MAT;
            
            Data.(CurrentField).(CurrentPulseTypeName).Dendrite_DeltaFOverF   = DendriteDeltaFOverF;
            
        else
            NaNMat = ones(size(SomaRawValue_BackgroundSubtructed),'single');
            Data.(CurrentField).(CurrentPulseTypeName).Axon_DeltaFOverF       = NaNMat;
            Data.(CurrentField).(CurrentPulseTypeName).Dendrite_DeltaFOverF   = NaNMat;
        end  
    end              
end

return

function Data = SmoothActivity_CalculateResponsesStats(Data)

%% Initialization
FieldNames = fieldnames(Data);

% Free parameters
DeflectionPoint_LowDiffThreshold            = 5e-3;  % until you reach this lower threshold
DynamicRangeThresholdForResponsiveness      = 0.2;   % Minimal deltaF/F dynamic range for safely considering it as a real response rather than fluctuation  
MinimalChangeRequiredForActivityInitiation  = 0.05;  % 5% deltaF/F

% Moving-average (MA) filter for dynamic range calculation. It has a very short delay, but doesn't fix ripples as good as LPF.
a_MA = 1; b_MA = ones(1,5)/5;
D_MA = round(mean(grpdelay(b_MA,a_MA)));  % Filter delay

% Low-pass filter (LPF) for activity. It fixes ripples to some extent but has a long delay (~3 sec). Not appropriate for quantification of step response latency.
LPF = designfilt('lowpassfir', ...
  'PassbandFrequency',0.1,'StopbandFrequency',1, ...
  'PassbandRipple',1,'StopbandAttenuation',80, ...
  'DesignMethod','equiripple','SampleRate',10);
D = round(mean(grpdelay(LPF)));           % Filter delay

NeuronFieldNames_Raw  = {'Soma_DeltaFOverF','Dendrite_DeltaFOverF','Axon_DeltaFOverF'};
NeuronFieldNames      = {'Soma','Dendrite','Axon'};

%% Calculate latency and responsiveness, and smoothed responses
for f_ind = 1:length(FieldNames)
    CurrentField    = FieldNames{f_ind};               
    PulseTypeNames  = Data.(CurrentField).PulsesInfo.PulseFromButanone;    
    
    for p_ind = 1:length(PulseTypeNames)
        CurrentPulseName = PulseTypeNames{p_ind};
        
        if isfield(Data.(CurrentField),CurrentPulseName)
            CurrentStructure = Data.(CurrentField).(CurrentPulseName);                           
            CurrentTimeVec   = CurrentStructure.TimeVec;  
            NumOfRepeats     = size(CurrentStructure.Soma_DeltaFOverF,1); 
            VectorLength     = size(CurrentStructure.Soma_DeltaFOverF,2);  
            TimeZeroIndex    = find(CurrentTimeVec==0); 

            for n_ind = 1:length(NeuronFieldNames_Raw)
                CurrentNeuronFieldName_Raw   = NeuronFieldNames_Raw{n_ind};
                MAT                          = CurrentStructure.(CurrentNeuronFieldName_Raw);
                
                % Initializtion
                SmoothedMAT_ForSlowResponses   = MAT * NaN;  
                SmoothedMAT_MA                 = MAT * NaN;  
                DynamicRange                   = zeros(NumOfRepeats,1,'single')*NaN;
                Responsiveness                 = zeros(NumOfRepeats,1,'single')*NaN;
                Latency_ByValue                = zeros(NumOfRepeats,1,'single')*NaN; % Using the average window vector! Good also for fast responses.
                Latency_ByDiffValue            = zeros(NumOfRepeats,1,'single')*NaN; % Using the LPF vector
                Latency                        = zeros(NumOfRepeats,1,'single')*NaN; % First timing from the two above
                AbsoluteValueIsDefiningLatency = zeros(NumOfRepeats,1,'single')*NaN;
                Latency_ByDiffValue_MA              = zeros(NumOfRepeats,1,'single')*NaN;
                Latency_MA                          = zeros(NumOfRepeats,1,'single')*NaN;
                AbsoluteValueIsDefiningLatency_MA   = zeros(NumOfRepeats,1,'single')*NaN;
                
                % Process curves and calculate responsiveness and latency 
                for rep_ind = 1:NumOfRepeats
                    CurrentVec = MAT(rep_ind,:);
                    if all(isnan(CurrentVec))
                        continue
                    end
                    %%%  Smooth Non-normalized deltaF/F
                    CurrentSmoothedVec_ForDR = filter(b_MA,a_MA,CurrentVec'); 
                    CurrentSmoothedVec_ForDR = CurrentSmoothedVec_ForDR(D_MA+1:end)';
                    CurrentSmoothedVec_ForDR((VectorLength-D_MA):VectorLength) = NaN;

                    %%%  Smooth Non-normalized deltaF/F
                    CurrentSmoothedVec = filter(LPF,CurrentVec'); 
                    CurrentSmoothedVec = CurrentSmoothedVec(D+1:end)';
                    CurrentSmoothedVec((VectorLength-D):VectorLength) = NaN;
                    FrameBeforeNaNs = find(~isnan(CurrentSmoothedVec),1,'last');
                    IndicesForAmplitudeNormalization = (FrameBeforeNaNs-10):(FrameBeforeNaNs);
                    NormlizingFactor         = mean(CurrentSmoothedVec(IndicesForAmplitudeNormalization)) / mean(CurrentVec(IndicesForAmplitudeNormalization));
                    CurrentSmoothedVec       = CurrentSmoothedVec / NormlizingFactor;

                    %%%  Differentiate smoothed activity
                    CurrentDiffVec     = [NaN, diff(CurrentSmoothedVec)]/(CurrentTimeVec(2)-CurrentTimeVec(1));

                    % Assign to Matrix      
                    SmoothedMAT_ForSlowResponses(rep_ind,:) = CurrentSmoothedVec;
                    SmoothedMAT_MA(rep_ind,:)               = CurrentSmoothedVec_ForDR;
                    
                    DynamicRange(rep_ind)   = max([0 CurrentSmoothedVec_ForDR(TimeZeroIndex:end)]);             
                    Responsiveness(rep_ind) = DynamicRange(rep_ind)>DynamicRangeThresholdForResponsiveness;   
                    
                    if Responsiveness(rep_ind) % Neuron Is Responsive                       

                        %% Find Latency by deflection point and by threshold value
                        % Deflection point by differential
                        IndexAfterDeflectionPoint = find(CurrentSmoothedVec_ForDR>DynamicRangeThresholdForResponsiveness,1,'first');
                        IndexOfDeflectionPoint    = find(CurrentDiffVec(1:IndexAfterDeflectionPoint)<=DeflectionPoint_LowDiffThreshold,1,'last');
                        if TimeZeroIndex > IndexOfDeflectionPoint  
                           IndexOfDeflectionPoint = TimeZeroIndex;
                        end 
                        if ~isempty(IndexOfDeflectionPoint)
                            Latency_ByDiffValue(rep_ind) = CurrentTimeVec(IndexOfDeflectionPoint);
                        end
                        
                        % Deflection point by max value with at least minimal differential            
                        IndexAfterDeflectionPoint = find(CurrentSmoothedVec_ForDR>DynamicRangeThresholdForResponsiveness,1,'first');
                        IndexOfHighActivity       = find(CurrentSmoothedVec_ForDR(1:IndexAfterDeflectionPoint)<=MinimalChangeRequiredForActivityInitiation,1,'last')+1;
                        if IndexOfHighActivity< TimeZeroIndex
                            % assume error, assign maximum length
                             IndexOfHighActivity = length(CurrentSmoothedVec_ForDR);
                        end
                        Latency_ByValue(rep_ind) = CurrentTimeVec(IndexOfHighActivity);                             

                        %%%% Deflection point, by activity or by diff(activity)?
                        if isempty(IndexOfDeflectionPoint) || (IndexOfHighActivity < IndexOfDeflectionPoint) 
                            IndexForLatency = IndexOfHighActivity;
                            AbsoluteValueIsDefiningLatency(rep_ind) = 1;
                        else
                            IndexForLatency = IndexOfDeflectionPoint;
                            AbsoluteValueIsDefiningLatency(rep_ind) = 0;
                        end
                        Latency(rep_ind) = CurrentTimeVec(IndexForLatency);    
                                                
                        %% Latency By deflection point (differential) using Average Window vectors. Good for Fast responses
                        IndexAfterDeflectionPoint = find(CurrentSmoothedVec_ForDR>DynamicRangeThresholdForResponsiveness,1,'first');
                        IndexOfDeflectionPoint    = find(diff(CurrentSmoothedVec_ForDR(1:IndexAfterDeflectionPoint))<=DeflectionPoint_LowDiffThreshold,1,'last');
                        if TimeZeroIndex > IndexOfDeflectionPoint  
                           IndexOfDeflectionPoint = TimeZeroIndex;
                        end 
                        if ~isempty(IndexOfDeflectionPoint)
                            Latency_ByDiffValue_MA(rep_ind) = CurrentTimeVec(IndexOfDeflectionPoint);
                        end
                        %%%% Deflection point, by activity or by diff(activity)?
                        if isempty(IndexOfDeflectionPoint) || (IndexOfHighActivity < IndexOfDeflectionPoint) 
                            IndexForLatency = IndexOfHighActivity;
                            AbsoluteValueIsDefiningLatency_MA(rep_ind) = 1;
                        else
                            IndexForLatency = IndexOfDeflectionPoint;
                            AbsoluteValueIsDefiningLatency_MA(rep_ind) = 0;
                        end
                        Latency_MA(rep_ind) = CurrentTimeVec(IndexForLatency);                                                 
                        
                    end
                end                
                Responsiveness = Responsiveness * 100; % Change from [0,1] scale to [0,100] scale
                
                % Assign to CurrentStructure
                CurrentStructure.([NeuronFieldNames{n_ind},'_DeltaFOverF_Smoothed'])         = SmoothedMAT_ForSlowResponses;
                CurrentStructure.([NeuronFieldNames{n_ind},'_DeltaFOverF_MA_Smoothed'])      = SmoothedMAT_MA;
                CurrentStructure.([NeuronFieldNames{n_ind},'_DeltaFOverF_MAX'])              = DynamicRange;
                CurrentStructure.([NeuronFieldNames{n_ind},'_Responsiveness'])               = Responsiveness;
                CurrentStructure.([NeuronFieldNames{n_ind},'_LatencyByValue_InSec'])         = Latency_ByValue;
                CurrentStructure.([NeuronFieldNames{n_ind},'_LatencyByDiffValue_InSec'])     = Latency_ByDiffValue;
                CurrentStructure.([NeuronFieldNames{n_ind},'_Latency_InSec'])                = Latency;
                CurrentStructure.([NeuronFieldNames{n_ind},'_LatencyDefinedByValue'])        = AbsoluteValueIsDefiningLatency;                  
                CurrentStructure.([NeuronFieldNames{n_ind},'_LatencyByDiffValue_MA_InSec'])  = Latency_ByDiffValue_MA;
                CurrentStructure.([NeuronFieldNames{n_ind},'_Latency_MA_InSec'])             = Latency_MA;
                CurrentStructure.([NeuronFieldNames{n_ind},'_LatencyDefinedByValue_MA'])     = AbsoluteValueIsDefiningLatency_MA;                  
                
            end                                                 
            % Assign CurrentStructure back to Data 
            Data.(CurrentField).(CurrentPulseName) = CurrentStructure;                        
        end    
    end       
end

return

function Data = CalculateStatsOverRepeats (Data)

StrainNames                     = fieldnames(Data);
CompartmentFieldNames           = {'Soma','Dendrite','Axon'};
MaxCurrentNumberOfRepeats       = 40;  % free parameter. Not relevant for calculation- just make sure it's higher than the maximal number of repeats

% Calculate Averages Responses over all repeats of the same strain and pulse   
for s_ind = 1:length(StrainNames)
    CurrentStrain       = StrainNames{s_ind};   
    PulseTypeNames      = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;   % All possible pulses in the current experiment type
    FinalConcentrations = Data.(CurrentStrain).PulsesInfo.FinalConcentration;  % Final concentration in microM, corresponding to indices in 'PulseTypeNames' 
    
    NumberOfRepeats      = squeeze(FindNumberOfRepeats(Data, s_ind));
    NumberOfRepeats      = NumberOfRepeats(1:length(FinalConcentrations),:);
    
    RelevantPulseIndices = ~isnan(NumberOfRepeats);
    
    clear Stats
    for c_ind = 1:length(CompartmentFieldNames)
        CurrentCompartment = CompartmentFieldNames{c_ind};
        NotNaNs            = RelevantPulseIndices(:,c_ind);
        CurrentNumberOfRepeats     = NumberOfRepeats(NotNaNs,c_ind);
        CurrentPulses              = PulseTypeNames(NotNaNs);
        CurrentFinalConcentrations = FinalConcentrations(NotNaNs);
        
        CurrentField                      = [CurrentCompartment,'_Information'];
        Stats.(CurrentField).NumberOfRepeats     = CurrentNumberOfRepeats;
        Stats.(CurrentField).PulseTypeNames      = CurrentPulses;
        Stats.(CurrentField).FinalConcentrations = CurrentFinalConcentrations;
        Stats.([CurrentCompartment,'_Responsiveness']).Matrix      = zeros(length(CurrentPulses),MaxCurrentNumberOfRepeats,'single')*NaN; 
        Stats.([CurrentCompartment,'_Latency_InSec']).Matrix       = zeros(length(CurrentPulses),MaxCurrentNumberOfRepeats,'single')*NaN; 
        Stats.([CurrentCompartment,'_Latency_MA_InSec']).Matrix    = zeros(length(CurrentPulses),MaxCurrentNumberOfRepeats,'single')*NaN; 
        Stats.([CurrentCompartment,'_DeltaFOverF_MAX']).Matrix     = zeros(length(CurrentPulses),MaxCurrentNumberOfRepeats,'single')*NaN; 
        
        
        for p_ind = 1:length(CurrentPulses)
            CurrentPulseName = CurrentPulses{p_ind};
            CurrentStructure = Data.(CurrentStrain).(CurrentPulseName);
            
            % Responsivenss
            CurrentField                                    = [CurrentCompartment,'_Responsiveness'];
            Vec                                             = CurrentStructure.(CurrentField);
            Stats.(CurrentField).Matrix(p_ind,1:length(Vec))= Vec;
            Stats.(CurrentField).Median(p_ind)              = nanmedian(Vec);
            Stats.(CurrentField).Mean(p_ind)                = nanmean(Vec);
            Stats.(CurrentField).STD(p_ind)                 = nanstd(Vec,[],1);
            if length(find(~isnan(Vec)))~=CurrentNumberOfRepeats(p_ind)
                disp('WTF')
                keyboard
            end
            
            % Latency
            CurrentField                                    = [CurrentCompartment,'_Latency_InSec'];            
            Vec                                             = CurrentStructure.(CurrentField);
            Stats.(CurrentField).Matrix(p_ind,1:length(Vec))= Vec;
            Stats.(CurrentField).Median(p_ind)              = nanmedian(Vec);
            Stats.(CurrentField).Mean(p_ind)                = nanmean(Vec);
            Stats.(CurrentField).STD(p_ind)                 = nanstd(Vec,[],1);
            
            CurrentField                                    = [CurrentCompartment,'_Latency_MA_InSec'];            
            Vec                                             = CurrentStructure.(CurrentField);
            Stats.(CurrentField).Matrix(p_ind,1:length(Vec))= Vec;
            Stats.(CurrentField).Median(p_ind)              = nanmedian(Vec);
            Stats.(CurrentField).Mean(p_ind)                = nanmean(Vec);
            Stats.(CurrentField).STD(p_ind)                 = nanstd(Vec,[],1);
            
            % Max deltaF/F
            CurrentField                                    = [CurrentCompartment,'_DeltaFOverF_MAX'];            
            Vec                                             = CurrentStructure.(CurrentField);
            Stats.(CurrentField).Matrix(p_ind,1:length(Vec))= Vec;
            Stats.(CurrentField).Median(p_ind)              = nanmedian(Vec);
            Stats.(CurrentField).Mean(p_ind)                = nanmean(Vec);
            Stats.(CurrentField).STD(p_ind)                 = nanstd(Vec,[],1);                          
        end
        
        CurrentField = [CurrentCompartment,'_Responsiveness'];     Stats.(CurrentField).SEM = Stats.(CurrentField).STD ./ sqrt(CurrentNumberOfRepeats)';
        CurrentField = [CurrentCompartment,'_Latency_InSec'];      Stats.(CurrentField).SEM = Stats.(CurrentField).STD ./ sqrt(CurrentNumberOfRepeats)';
        CurrentField = [CurrentCompartment,'_Latency_MA_InSec'];   Stats.(CurrentField).SEM = Stats.(CurrentField).STD ./ sqrt(CurrentNumberOfRepeats)';                    
        CurrentField = [CurrentCompartment,'_DeltaFOverF_MAX'];    Stats.(CurrentField).SEM = Stats.(CurrentField).STD ./ sqrt(CurrentNumberOfRepeats)';
                
        Data.(CurrentStrain).Stats = Stats;                
    end        
end

return

function [NumberOfRepeats, CellToPrint] = FindNumberOfRepeats(Data, RelevantStrainIndices) 

%% NumberOfRepeats
StrainNames           = fieldnames(Data);
MaxNumberOfPulseTypes = 20; 
NeuronFieldNames      = {'soma','dendrite','axon'};
NumberOfRepeats       = zeros(length(RelevantStrainIndices), MaxNumberOfPulseTypes, length(NeuronFieldNames), 'single') * NaN; % Initialization
CellToPrint           = {'Strain Name','Pulse Name','n (Soma)','n (Axon)','n (Dendrite)'};
CellRow               = 1;

for s_count = 1:length(RelevantStrainIndices)
    s_ind             = RelevantStrainIndices(s_count);
    CurrentStrain     = StrainNames{s_ind};
    CurrentStructure  = Data.(CurrentStrain);
%     CurrentNumberOfPulses =  
    CurrentPulseNames = CurrentStructure.PulsesInfo.PulseFromButanone; % some of which may be empty ! 
    NumOfPulses       = length(CurrentPulseNames);                     % some of which may be empty ! 
    
    for p_ind = 1:NumOfPulses
        PulseName = CurrentPulseNames{p_ind};
        if isfield(CurrentStructure,PulseName)
%             NumberOfRepeats(s_count,p_ind,1) = size(CurrentStructure.(PulseName).NeuronValue_BackgroundSubtructed,1);
%             NumberOfRepeats(s_count,p_ind,2) = length(find(all(~isnan(CurrentStructure.(PulseName).DendriteValue_BackgroundSubtructed),2)));
%             NumberOfRepeats(s_count,p_ind,3) = length(find(all(~isnan(CurrentStructure.(PulseName).AxonValue_BackgroundSubtructed),2)));            
            NumberOfRepeats(s_count,p_ind,1) = size(CurrentStructure.(PulseName).Soma_DeltaFOverF,1);
            NumberOfRepeats(s_count,p_ind,2) = length(find(~all(isnan(CurrentStructure.(PulseName).Dendrite_DeltaFOverF),2)));
            NumberOfRepeats(s_count,p_ind,3) = length(find(~all(isnan(CurrentStructure.(PulseName).Axon_DeltaFOverF),2)));            
            
            CellRow = CellRow+1;
            CellToPrint(CellRow,1:5) = {CurrentStrain, PulseName, NumberOfRepeats(s_count,p_ind,1), NumberOfRepeats(s_count,p_ind,2), NumberOfRepeats(s_count,p_ind,3)};
        end
    end    
end

% MaxNumberOfRepeats = max(NumberOfRepeats(:));
% xlswrite('D:\Paper\Repeats_Test.xlsx',CellToPrint)


return

function Data = SmoothDyeAndCalculateStats(Data)

CompartmentFieldNames = {'Soma','Dendrite','Axon'};
StrainNames           = fieldnames(Data);
warning('off')
DyeVectorLength    = length(Data.str2_CX17256_ToBuffer_Slow.Butanone_d5_To_Buffer.Dye_TimeVec); % LenDyeVec = 1100; 
MaxNumberOfRepeats = 40;  % free parameter. It just needs to be higher than the maximal number of repeats

%% low pass filters - including free parameters
% dye
LPS_DyeSmoothing1 = designfilt('lowpassfir', ...
  'PassbandFrequency',0.01,'StopbandFrequency',3, ...
  'PassbandRipple',1,'StopbandAttenuation',30, ...
  'DesignMethod','equiripple','SampleRate',10);
D_Dye1 = round(mean(grpdelay(LPS_DyeSmoothing1)));

% fitting functions:
ft_cubic = fittype('1-a*x-b*x^2', 'dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'}); % initial dye dynamics smoothing
ft_poly3 = fittype('a*x+b*x^2+c*x^3', 'dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c'}); % initital dye derivative smoothing
ft_poly9 = fittype('poly9'); % late dye derivative smoothing

%% Calculate dye dynamics properties
for s_ind = 1:length(StrainNames)
    CurrentStrain              = StrainNames{s_ind};  
    PulsesNames                = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;
    InitialConcentration       = Data.(CurrentStrain).PulsesInfo.InitialConcentration;      
    NumberOfPulses             = length(PulsesNames);
    if ~isfield(Data.(CurrentStrain).(PulsesNames{1}),'Dye_TimeVec')
        continue
    end
    
    if length(InitialConcentration)>1
        % Pulses To buffer from different initial concentrations to similar final concentrations 
        InitialConcentrationPerPulse = InitialConcentration;
    else
        % Pulses from the same initial concentration to different final concentrations 
        InitialConcentrationPerPulse = InitialConcentration * ones(1,NumberOfPulses);
    end
    
    clear DyeInformation
    DyeInformation.InitialConcentrationPerPulse = InitialConcentrationPerPulse;
    DyeInformation.FinalConcentrationPerPulse   = FinalConcentrationPerPulse;
    DyeInformation.Dye_Raw         = zeros(NumberOfPulses,MaxNumberOfRepeats,DyeVectorLength,'single')*NaN;
    DyeInformation.Dye_Smoothed    = zeros(NumberOfPulses,MaxNumberOfRepeats,DyeVectorLength,'single')*NaN;
    DyeInformation.C               = zeros(NumberOfPulses,MaxNumberOfRepeats,DyeVectorLength,'single')*NaN;
    DyeInformation.C_Smoothed      = zeros(NumberOfPulses,MaxNumberOfRepeats,DyeVectorLength,'single')*NaN;  % Concentration in units of dilution, corrected by the known initial and final concentrations 
    DyeInformation.dCdt            = zeros(NumberOfPulses,MaxNumberOfRepeats,DyeVectorLength,'single')*NaN;  % dilution/sec
    DyeInformation.dCdt_Smoothed   = zeros(NumberOfPulses,MaxNumberOfRepeats,DyeVectorLength,'single')*NaN;  % dilution/sec
    DyeInformation.d2Cdt2          = zeros(NumberOfPulses,MaxNumberOfRepeats,DyeVectorLength,'single')*NaN;  % dilution/sec^2
    DyeInformation.DeltaCoverC     = zeros(NumberOfPulses,MaxNumberOfRepeats,DyeVectorLength,'single')*NaN;  % 1/sec
    
    DyeInformation.Stats.dCdt_Peak               = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;  % MINIMUM !!
    DyeInformation.Stats.dCdt_TimeToPeak         = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;  % MINIMUM !!
    DyeInformation.Stats.d2Cdt2_Peak             = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;  % MINIMUM !!
    DyeInformation.Stats.d2Cdt2_TimeToPeak       = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;  % MINIMUM !!
    DyeInformation.Stats.d2Cdt2_TimeToSignChange = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;  % MINIMUM !!
    DyeInformation.Stats.DeltaCoverC_Peak        = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;  % MINIMUM !!
    DyeInformation.Stats.DeltaCoverC_TimeToPeak  = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;  % MINIMUM !!
    DyeInformation.Stats.C_AtNeuronActivation.Soma           = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;
    DyeInformation.Stats.C_AtNeuronActivation.Dendrite       = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;  
    DyeInformation.Stats.C_AtNeuronActivation.Axon           = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;  
    DyeInformation.Stats.dCdt_AtNeuronActivation             = DyeInformation.Stats.C_AtNeuronActivation;
    DyeInformation.Stats.d2Cdt2_AtNeuronActivation           = DyeInformation.Stats.C_AtNeuronActivation;
    DyeInformation.Stats.DeltaCoverC_AtNeuronActivation      = DyeInformation.Stats.C_AtNeuronActivation;    
    DyeInformation.Stats.C_AtNeuronActivation_MA             = DyeInformation.Stats.C_AtNeuronActivation ;
    DyeInformation.Stats.dCdt_AtNeuronActivation_MA          = DyeInformation.Stats.C_AtNeuronActivation;
    DyeInformation.Stats.d2Cdt2_AtNeuronActivation_MA        = DyeInformation.Stats.C_AtNeuronActivation;
    DyeInformation.Stats.DeltaCoverC_AtNeuronActivation_MA   = DyeInformation.Stats.C_AtNeuronActivation;

    for p_ind = 1:NumberOfPulses
        CurrentPulseName = PulsesNames{p_ind}; 
        InitialC      = InitialConcentrationPerPulse(p_ind);
        FinalC        = FinalConcentrationPerPulse(p_ind);
        DyeOnset      = double(Data.(CurrentStrain).(CurrentPulseName).Dye);        
        Dye_Raw       = 1 - DyeOnset;                       % Odor removal pulses BEFORE BASELINE CORRECTION   
        NumOfRepeats  = size(DyeOnset,1);               
        DyeTimeVec    = Data.(CurrentStrain).(CurrentPulseName).Dye_TimeVec;   % in seconds  
        if p_ind==1                
            DyeInformation.TimeVectorInSec = DyeTimeVec; % assign this only once
        end        
        TimeZeroIndex = find(DyeTimeVec==0);         
            
        for rep_ind = 1:NumOfRepeats
            CurrentDye_Raw                  = Dye_Raw(rep_ind,:);    
            LastNoNaN                       = find(~isnan(CurrentDye_Raw),1,'last');
            Factor                          = nanmean(CurrentDye_Raw((LastNoNaN-10*2):LastNoNaN));  % Normlize the dye zero using 2 last seconds
            CurrentDye_Raw                  = (CurrentDye_Raw - Factor)/(1-Factor);
            CurrentDye_Raw(1:TimeZeroIndex) = 1;            
            DyeInformation.Dye_Raw(p_ind,rep_ind,:) =  CurrentDye_Raw;   
            
            Concentration_Raw                 = FinalC + (InitialC - FinalC)* CurrentDye_Raw;             
            DyeInformation.C(p_ind,rep_ind,:) = Concentration_Raw;        

            %% Smooth dye pattern and its differential      
            %% A. Smooth dye pattern
            % Strategy for dye smoothing:
            %   at t<D_Dye1      >> cubic approximation
            %   at D_Dye1<t      >> LPS_DyeSmoothing1
            %   at t>peak slope  >> poly9 approximation

            % at t<D_Dye1  >> cubic approximation
            SmoothedDye0               = CurrentDye_Raw; 
            indices_LPS0               = TimeZeroIndex:(TimeZeroIndex+D_Dye1);
            x                          = (0:length(indices_LPS0)-1)';
            y                          = double(CurrentDye_Raw(indices_LPS0)');            
            [cf,~,~]                   = fit(x,y,ft_cubic);
            y_fit                      = 1-cf.a*x-cf.b*x.^2;            
            SmoothedDye0(indices_LPS0) = y_fit;

            % at D_Dye1<t<D_DiffDye  >> LPS_DyeSmoothing1
            SmoothedDye1   = filter(LPS_DyeSmoothing1,SmoothedDye0); 
            SmoothedDye1   = SmoothedDye1(D_Dye1+1:end)';
            SmoothedDye1   = SmoothedDye1/max(SmoothedDye1(TimeZeroIndex-20:TimeZeroIndex));        
            indices        = 1:indices_LPS0(end);                       
            SmoothedDye1(indices) = SmoothedDye0(indices);
            LastNoNaNIndex = find(~isnan(SmoothedDye0),1,'last');
            indices        = (LastNoNaNIndex-D_Dye1):LastNoNaNIndex;    
            SmoothedDye1(indices) = NaN;
            SmoothedDye = SmoothedDye1';
            LastNoNaNIndex = find(~isnan(SmoothedDye),1,'last');
            indices = (LastNoNaNIndex+1):DyeVectorLength;    
            SmoothedDye(indices) = NaN; 

            %%%%%   late dye pattern ripples correction  %%%% 
            FirstIndex = TimeZeroIndex + 10;
%             Indices    = FirstIndex+(0:find(isnan(SmoothedDye(FirstIndex:end)),1,'first')-2);
%             x          = (0:length(Indices)-1)';
            Indices    = find(~isnan(SmoothedDye));
            Indices    = Indices(Indices>FirstIndex);
            y          = double(SmoothedDye(Indices)');            
            x          = Indices'-min(Indices);
            [cf_,~,~]  = fit(x,y,ft_poly9);
            y_fit      = feval(cf_,x);  
            SmoothedDye_LateTimeCorrection                     = SmoothedDye;
            SmoothedDye_LateTimeCorrection(Indices)            = y_fit;
            SmoothedDye_LateTimeCorrection(Indices(end)+1:end) = NaN;

            % find stitching index
            DiffVec = (SmoothedDye_LateTimeCorrection- SmoothedDye);
            minIndex = FirstIndex+1;
            maxIndex = FirstIndex+20;
            [~,relative_ind] = min(abs(DiffVec(minIndex:maxIndex)));
            StitchIndex = minIndex+relative_ind-1;
            SmoothedDye_LateTimeCorrection(1:StitchIndex) = SmoothedDye(1:StitchIndex);

            SmoothedDye = SmoothedDye_LateTimeCorrection;  % Use only the corrected vector for further analysis
            LastNoNaN   = find(~isnan(SmoothedDye),1,'last');
            NaNIndices  = find(isnan(SmoothedDye(1:LastNoNaN)));
            
            if ~isempty(NaNIndices)
                RelevantIndices = find(~isnan(SmoothedDye(1:LastNoNaN)));                                  
                InterpolatedValues              = interp1(RelevantIndices, SmoothedDye(RelevantIndices),NaNIndices,'cubic');                    
                SmoothedDye(NaNIndices) = InterpolatedValues;                               
            end            
            SmoothedDye(SmoothedDye<0)= NaN;
            C_Smoothed                                   = FinalC + (InitialC - FinalC)* SmoothedDye;             
            DyeInformation.Dye_Smoothed(p_ind,rep_ind,:) = SmoothedDye;                                  
            DyeInformation.C_Smoothed(p_ind,rep_ind,:)   = C_Smoothed;        

            %% B. Smooth diff dye pattern
            dCdt                                 = [NaN, diff(C_Smoothed)]/(DyeTimeVec(2)-DyeTimeVec(1));     
            dCdt(TimeZeroIndex+1)                = NaN;
            RelevantIndices         = TimeZeroIndex:TimeZeroIndex+2;
            NotNaNs                 = setdiff(find(~isnan(dCdt)), RelevantIndices);
            InterpolatedValues      = interp1(NotNaNs, dCdt(NotNaNs),RelevantIndices,'cubic'); 
            dCdt(RelevantIndices)   = InterpolatedValues;
            
            DyeInformation.dCdt(p_ind,rep_ind,:) = dCdt;              
            [CurrentSlopePeakValue, CurrentSlopePeakIndex]           = min(dCdt); 
            
            if CurrentSlopePeakValue == 0 % constant concentration
                PeakDetected = false;
                dCdt_Smoothed = dCdt;
            else            
                PeakDetected = true;
                % at t<= slope peak  --> use poly3 approximation
                dCdt_Smoothed_0      = dCdt; 
%                 MaxIndexForSmoothing = CurrentSlopePeakIndex-1*10; % up to 1 seconds before peak 
                MaxIndexForSmoothing = CurrentSlopePeakIndex+2*10; % up to 2 seconds AFTER peak 
                indices_LPS0         = TimeZeroIndex:MaxIndexForSmoothing;
                x                    = (0:length(indices_LPS0)-1)';
                y                    = double(dCdt(indices_LPS0)');  
                RelevantIndices      = find(~isnan(dCdt(indices_LPS0)));
                x                    = x(RelevantIndices);
                y                    = y(RelevantIndices);     
                [cf,~,~]             = fit(x,y,ft_poly3);
                y_fit                = feval(cf,x);  
                dCdt_Smoothed_0(indices_LPS0(RelevantIndices)) = y_fit;
                
                % Stitch to the rest of the vector using cubic interpolation                  
                [~, BeforeNaNIndex] = min(dCdt_Smoothed_0);  
                AfterNaNIndex       = MaxIndexForSmoothing+1;   
                AfterNaNIndex       = max(AfterNaNIndex, BeforeNaNIndex+2);
                dCdt_Smoothed       = dCdt_Smoothed_0;
                RelevantIndices     = BeforeNaNIndex+1:AfterNaNIndex-1;
                dCdt_Smoothed(RelevantIndices)  = NaN; 
                NotNaNs                         = find(~isnan(dCdt_Smoothed));
                InterpolatedValues              = interp1(NotNaNs, dCdt_Smoothed(NotNaNs),RelevantIndices,'cubic'); 
                dCdt_Smoothed(RelevantIndices)  = InterpolatedValues;                                  
                
                DerivativeStitchingIndex1 = RelevantIndices(1);
                DerivativeStitchingIndex2 = RelevantIndices(end)+1;             
            end
%             
            [CurrentSlopePeakValue, CurrentSlopePeakIndex] = min(dCdt);  
            DyeInformation.dCdt_Smoothed(p_ind,rep_ind,:)  = dCdt_Smoothed;   

            %% Calculation of deltaCoverC == dCdt / C
            % Use ripple corrected dye and diffdye for calculation of deltaCoverC 
            DeltaCoverC  = dCdt_Smoothed./C_Smoothed;                  % [dCdt/C]    
            LimitForReliableDeltaCoverC = find(SmoothedDye<0.05,1);   % less than 5% dye signal is too noisy for deltaC/C measurements
            DeltaCoverC(LimitForReliableDeltaCoverC:end) = NaN;               
            
            DyeInformation.DeltaCoverC(p_ind,rep_ind,:)                = DeltaCoverC;    % [dCdt/C]                      
            [CurrentDeltaCoverCPeakValue, CurrentDeltaCoverCPeakIndex] = min(DeltaCoverC);

            %% Dye second derivative           
            d2Cdt2  = [NaN, diff(dCdt_Smoothed)]/(DyeTimeVec(2)-DyeTimeVec(1));   
            if PeakDetected
                NotNaN  = find(~isnan(d2Cdt2));
                d2Cdt2(DerivativeStitchingIndex1-2:DerivativeStitchingIndex1+2) = NaN;
                d2Cdt2(DerivativeStitchingIndex2-2:DerivativeStitchingIndex2+2) = NaN;
                RelevantIndices        = [DerivativeStitchingIndex1-2:DerivativeStitchingIndex1+2  DerivativeStitchingIndex2-2:DerivativeStitchingIndex2+2]; 
                InterpolatedValues     = interp1(NotNaN, d2Cdt2(NotNaN),RelevantIndices,'cubic'); 
                d2Cdt2(RelevantIndices)= InterpolatedValues;                                       
            end
            DyeInformation.d2Cdt2(p_ind,rep_ind,:)               = d2Cdt2;    % [d^2C/dt^2]                      
            [Current_d2Cdt2_PeakValue, Current_d2Cdt2_PeakIndex] = min(d2Cdt2);
            Current_d2Cdt2_SignSwitchIndex                       = TimeZeroIndex + find(d2Cdt2(TimeZeroIndex+1:end)>0,1,'first');
                       
            % STATS   
            if PeakDetected
                DyeInformation.Stats.dCdt_Peak(p_ind,rep_ind)               = CurrentSlopePeakValue;
                DyeInformation.Stats.dCdt_TimeToPeak(p_ind,rep_ind)         = DyeTimeVec(CurrentSlopePeakIndex);
                DyeInformation.Stats.d2Cdt2_Peak(p_ind,rep_ind)             = Current_d2Cdt2_PeakValue;
                DyeInformation.Stats.d2Cdt2_TimeToPeak(p_ind,rep_ind)       = DyeTimeVec(Current_d2Cdt2_PeakIndex);
                DyeInformation.Stats.d2Cdt2_TimeToSignChange(p_ind,rep_ind) = DyeTimeVec(Current_d2Cdt2_SignSwitchIndex);
                DyeInformation.Stats.DeltaCoverC_Peak(p_ind,rep_ind)        = CurrentDeltaCoverCPeakValue;
                DyeInformation.Stats.DeltaCoverC_TimeToPeak(p_ind,rep_ind)  = DyeTimeVec(CurrentDeltaCoverCPeakIndex);            
            end
            
            for c_ind = 1:length(CompartmentFieldNames)
                CurrentCompartment = CompartmentFieldNames{c_ind};
                LatencyInSec       = Data.(CurrentStrain).Stats.([CurrentCompartment,'_Latency_InSec']).Matrix(p_ind,rep_ind);
                if ~isnan(LatencyInSec)
                    CurrentIndex = find(DyeTimeVec >= LatencyInSec-eps,1,'first');                
                    DyeInformation.Stats.C_AtNeuronActivation.(CurrentCompartment)(p_ind,rep_ind)           = C_Smoothed(CurrentIndex);
                    DyeInformation.Stats.dCdt_AtNeuronActivation.(CurrentCompartment)(p_ind,rep_ind)        = dCdt_Smoothed(CurrentIndex);
                    DyeInformation.Stats.d2Cdt2_AtNeuronActivation.(CurrentCompartment)(p_ind,rep_ind)      = d2Cdt2(CurrentIndex);
                    DyeInformation.Stats.DeltaCoverC_AtNeuronActivation.(CurrentCompartment)(p_ind,rep_ind) = DeltaCoverC(CurrentIndex);      
                end
                LatencyInSec       = Data.(CurrentStrain).Stats.([CurrentCompartment,'_Latency_MA_InSec']).Matrix(p_ind,rep_ind);
                if ~isnan(LatencyInSec)
                    CurrentIndex = find(DyeTimeVec >= LatencyInSec-eps,1,'first');                
                    DyeInformation.Stats.C_AtNeuronActivation_MA.(CurrentCompartment)(p_ind,rep_ind)           = C_Smoothed(CurrentIndex);
                    DyeInformation.Stats.dCdt_AtNeuronActivation_MA.(CurrentCompartment)(p_ind,rep_ind)        = dCdt_Smoothed(CurrentIndex);
                    DyeInformation.Stats.d2Cdt2_AtNeuronActivation_MA.(CurrentCompartment)(p_ind,rep_ind)      = d2Cdt2(CurrentIndex);
                    DyeInformation.Stats.DeltaCoverC_AtNeuronActivation_MA.(CurrentCompartment)(p_ind,rep_ind) = DeltaCoverC(CurrentIndex);      
                end
            end
                       
            %% Plots
%             figure;
%             subplot(2,2,1); 
%             plot(DyeTimeVec,squeeze(DyeInformation.C_SuperSmoothedForDiff(p_ind,rep_ind,:)),'r-'); hold on;
%             plot(DyeTimeVec,squeeze(DyeInformation.C_Smoothed(p_ind,rep_ind,:)),'b:'); 
%             plot(DyeTimeVec,squeeze(DyeInformation.C(p_ind,rep_ind,:)),'k.'); 
%             title('C')
%             xlim([-10 90])
%             
%             subplot(2,2,2); 
%             plot(DyeTimeVec,squeeze(DyeInformation.dCdt_Smoothed(p_ind,rep_ind,:)),'b:'); hold on;
%             plot(DyeTimeVec,squeeze(DyeInformation.dCdt(p_ind,rep_ind,:)),'k.'); 
%             title('dCdt')
%             xlim([-10 90])
%             
%             subplot(2,2,3); 
%             plot(DyeTimeVec,squeeze(DyeInformation.d2Cdt2(p_ind,rep_ind,:)),'r-'); hold on;
%             title('d^2C/dt^2')
%             xlim([-10 90])
%             
%             subplot(2,2,4); 
%             plot(DyeTimeVec,squeeze(DyeInformation.dCdtOverC(p_ind,rep_ind,:)),'r-'); hold on;
%             xlim([-10 90])
%             title('\DeltaC/C')            
                         
        end
    end    
    Data.(CurrentStrain).DyeInformation = DyeInformation;           

end

warning('on')

return

%% ACT model parameters calculation and bootstrap for evaluation of responsiveness and latency errors and predictions
% Calculate Thresholds and threshold constants (K)
function [Data, ModelParams] = Calculate_Thresholds_And_K_BootstrapForErrorbars (Data, NumOfIterations)

StrainNames           = fieldnames(Data);
RelevantStrainIndices = setdiff(1:length(StrainNames), find(strcmpi(StrainNames,'str2_CX17256_ToBuffer_Slow')));
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
ConcentrationScalingFactor = 11.16*1e6;  % Concentration units in microM

% % disp([datestr(now),'  --  loading existing results for T, n, and K with (1750 iterations) and without boostrapping'])
% % load('H:\BootstrapThresholdForStepResponse_Combined.mat','Data')

%% No Bootstrapping
disp([datestr(now),'  --  Calculating T, n, and K using all data (no boostrapping)'])
[Data, ModelParams] = CalculateThresholdFor50PercentResponse_CalculateK (Data, RelevantStrainIndices); % Using all data, without bootsrtapping

%% With Bootstrapping
%% IMPORTANT NOTE !!
%%%%% The bootstrap loop is slow: ~9min/iteration (iteration= one random bootstrap over all experiment types and compartments).      
%     It is recommended to run it on several windows for parallel computation and then combine the data before further analysis. 
%%%   For combining parallel sessions use the function: 'Combine_BootstrapThresholdForStepResponse'  
if NumOfIterations>5
    disp([datestr(now),'  --  NOTE!! The bootstrap loop is slow. Parallel computation is highly recommended. See details within script.' ])  
end

disp([datestr(now),'  --  Calculating T, n, and K with boostrapping. Number of iterations = ',num2str(NumOfIterations)])
% Based on the function 'CalculateThresholdFor50PercentResponse_CalculateK' 

% Fit and Fit Options  - Concentration units in microM
ft_T_and_n         = fittype('T^n/(T^n+C^n) *100', 'dependent',{'y'},'independent',{'C'},'coefficients',{'T','n'});
opts               = fitoptions( 'Method', 'NonlinearLeastSquares','robust','Bisquare' );
opts.Display       = 'Off';
opts.MaxFunEvals   = 10000;
opts.MaxIter       = 10000;
% opts.TolFun        = 1e-10;
% opts.TolX          = 1e-10;
opts.TolFun        = 1e-14;
opts.TolX          = 1e-14;
opts.DiffMinChange = 1e-10;

rng('shuffle');
% The bootstrap loop- SLOW
for s_ind = 1:length(RelevantStrainIndices)
    CurrentStrainIndex   = RelevantStrainIndices(s_ind);
    CurrentStrainName    = StrainNames{CurrentStrainIndex};    
    InitialConcentration = Data.(CurrentStrainName).PulsesInfo.InitialConcentration; 
    Data.(CurrentStrainName).Stats.InitialConcentration = InitialConcentration;        % Units of microM
    InitialConcentration_Dilution = InitialConcentration /ConcentrationScalingFactor;  % Units of Dilution 
    
    for c_ind = 1:length(CompartmentFieldNames)        
        CurrentComparment     = CompartmentFieldNames{c_ind};
        CurrentField          = [CurrentComparment,'_Responsiveness']; 
        ResponsivenessMatrix  = Data.(CurrentStrainName).Stats.(CurrentField).Matrix;
        RepeatsPerCondition   = Data.(CurrentStrainName).Stats.([CurrentComparment,'_Information']).NumberOfRepeats;
        
        Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_InitialThresholdGuess']).Vector = zeros(1,NumOfIterations)*NaN;
        Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_Threshold']).Vector             = zeros(1,NumOfIterations)*NaN;
        Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_n']).Vector                     = zeros(1,NumOfIterations)*NaN;
        
        disp([newline, newline, newline, datestr(now),'  --  ',CurrentStrainName,' , ',CurrentComparment, newline]);
                       
        c1 = clock;
        for it = 1:NumOfIterations
            if it==1
                disp([datestr(now),'  --  Iteration #',num2str(it),' / ',num2str(NumOfIterations)])
            else
                c2 = clock;
                EstimatedTimeLeft = (NumOfIterations-it+1)*(etime(c2,c1)/(it-1)/60/60);
                disp([datestr(now),'  --  Iteration #',num2str(it),' / ',num2str(NumOfIterations),...
                                   ' . EstimatedTimeLeft: ',num2str(EstimatedTimeLeft),' hours']);
            end
            % iterate over experiments and find current responsivness vector for bootstrap threshold calculation  
            % Make sure to keep a constant number of repeats per experimental condition   
            FinalConcentrations          = Data.(CurrentStrainName).Stats.([CurrentComparment,'_Information']).FinalConcentrations;                                          
            FinalConcentrations_Dilution = FinalConcentrations / ConcentrationScalingFactor;  % Units of Dilution            
            
            CurrentResponsivenessMatrix = ResponsivenessMatrix * NaN; % Initialization
            for fc_ind = 1:length(FinalConcentrations)
                ResponsivenessVector = ResponsivenessMatrix(fc_ind,:);
                ResponsivenessVector = ResponsivenessVector(~isnan(ResponsivenessVector));
                Bootstrap_indices = randi(RepeatsPerCondition(fc_ind),1,RepeatsPerCondition(fc_ind));
                CurrentResponsivenessMatrix(fc_ind,1:RepeatsPerCondition(fc_ind)) = ResponsivenessVector(Bootstrap_indices);   
            end
%             CurrentResponsiveness = Data.(CurrentStrainName).Stats.(CurrentField).Mean;  % Without bootstrap 
            CurrentResponsiveness = nanmean(CurrentResponsivenessMatrix,2)';               % With bootstrap 
            if any(isnan(CurrentResponsiveness))
                keyboard
            end
    
            [CurrentThreshold, Current_n, Tguess] = CalculateThreshold_MM(double(FinalConcentrations_Dilution), double(CurrentResponsiveness), double(InitialConcentration_Dilution), ft_T_and_n, opts);  
            
            Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_InitialThresholdGuess']).Vector(it) = Tguess           * ConcentrationScalingFactor; % Units of microM         ;
            Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_Threshold']).Vector(it)             = CurrentThreshold * ConcentrationScalingFactor; % Units of microM         ;
            Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_n']).Vector(it)                     = Current_n;  
        end
    end
end
toc

% Stats over bootstrap results
for s_ind = 1:length(RelevantStrainIndices)
    CurrentStrainIndex = RelevantStrainIndices(s_ind);
    CurrentStrainName  = StrainNames{CurrentStrainIndex};    
    for c_ind = 1:length(CompartmentFieldNames)        
        CurrentComparment     = CompartmentFieldNames{c_ind};
        RepeatsPerCondition   = Data.(CurrentStrainName).Stats.([CurrentComparment,'_Information']).NumberOfRepeats;               
        MeanNumOfRepeats      = mean(RepeatsPerCondition);
        
        CurrentStructure         = Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_Threshold']);
        Vector                   = CurrentStructure.Vector;        
        CurrentStructure.Mean    = mean(Vector);
        CurrentStructure.Median  = median(Vector);
        CurrentStructure.STD     = std(Vector,1);
        CurrentStructure.SEM     = std(Vector,1)./ sqrt(MeanNumOfRepeats);
        CurrentStructure.MeanNumOfRepeats = MeanNumOfRepeats;
        Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_Threshold']) = CurrentStructure;

        CurrentStructure         = Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_n']);
        Vector                   = CurrentStructure.Vector;        
        CurrentStructure.Mean    = mean(Vector);
        CurrentStructure.Median  = median(Vector);
        CurrentStructure.STD     = std(Vector,1);
        CurrentStructure.SEM     = std(Vector,1)./ sqrt(MeanNumOfRepeats);
        CurrentStructure.MeanNumOfRepeats = MeanNumOfRepeats;                
        Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_n']) = CurrentStructure;
    end
end          

%% IMPORTANT NOTE !!
%%%%% The bootstrap loop is slow: ~9min/iteration (iteration= one random bootstrap over all experiment types and compartments).      
%     It is recommended to run it on several windows for parallel computation and then combine the data before further analysis. 
%%%   For combining parallel sessions use the function: 'Combine_BootstrapThresholdForStepResponse'  
% save('H:\BootstrapThresholdForStepResponse_1.mat','Data','ModelParams','-v7.3')
% save('H:\BootstrapThresholdForStepResponse_2.mat','Data','ModelParams','-v7.3') ...
% After all functions finished, run: 'Combine_BootstrapThresholdForStepResponse2'   
% Then:   load('H:\BootstrapThresholdForStepResponse_Combined.mat','Data','ModelParams')
% save BCK2
save('H:\BootstrapThresholdForStepResponse20190813_05.mat','Data','-v7.3')

%% Generate ModelParams
ListOfStrains            = {'str2_CX17256_FromD7_Fast','str2_CX17256_FromD6_Fast','str2_CX17256_FromD5_Fast','str2_CX17256_FromD4_Fast'};
CurrentTrainingCondition = 'Naive';
ModelParams.(CurrentTrainingCondition) = GenerateThresholdStructure_WithBootstrap(Data, ListOfStrains, CompartmentFieldNames);
ListOfStrains            = {'str2_CX17256_Desensitized_FromD6_Fast','str2_CX17256_Desensitized_FromD5_Fast','str2_CX17256_Desensitized_FromD4_Fast','str2_CX17256_Desensitized_FromD3_Fast'};
CurrentTrainingCondition = 'Desensitized';
ModelParams.(CurrentTrainingCondition) = GenerateThresholdStructure_WithBootstrap(Data, ListOfStrains, CompartmentFieldNames);

ListOfStrains            = {'str2_CX17255_AllFast'};
CurrentTrainingCondition = 'CX17255';
ModelParams.(CurrentTrainingCondition) = GenerateThresholdStructure_WithBootstrap(Data, ListOfStrains, CompartmentFieldNames);
ListOfStrains            = {'str2_Egl4_AllFast'};
CurrentTrainingCondition = 'Egl4';
ModelParams.(CurrentTrainingCondition) = GenerateThresholdStructure_WithBootstrap(Data, ListOfStrains, CompartmentFieldNames);
ListOfStrains            = {'str2_Egl4_ad450_AllFast'};
CurrentTrainingCondition = 'Egl4_ad450';
ModelParams.(CurrentTrainingCondition) = GenerateThresholdStructure_WithBootstrap(Data, ListOfStrains, CompartmentFieldNames);

%% Calculate K using exponential function
%%%%%%% Initialization of exponential fitting function: T = K*(1-exp(-x/K)) %%%%%
ft    = fittype('K*(1-exp(-x/K))','dependent',{'y'},'independent',{'x'},'coefficients',{'K'});
opts  = fitoptions( 'Method', 'NonlinearLeastSquares','robust','Bisquare' );
% opts.Display     = 'Off';
opts.MaxFunEvals = 10000;
opts.MaxIter     = 10000;
opts.TolFun      = 1e-15;
opts.TolX        = 1e-15;
opts.TolFun      = 1e-17;
opts.TolX        = 1e-17;
opts.DiffMinChange = 1e-9;
opts.DiffMaxChange = 1e-3;
% opts.StartPoint  = [1e-6 1 ]; 
opts.Upper       = 1e-3 ; 
opts.Lower       = 1e-9 ; 

warning('off','curvefit:fit:noStartPoint')
TrainingConditions = {'Naive','Desensitized'};

% No Bootstrap
for t_ind = 1:length(TrainingConditions)
    CurrentTrainingCondition = TrainingConditions{t_ind};
    for c_ind = 1:length(CompartmentFieldNames)
        CurrentComparment = CompartmentFieldNames{c_ind};    
        cf = fit(ModelParams.(CurrentTrainingCondition).(CurrentComparment).InitialConcetrations' / ConcentrationScalingFactor, ...
                 ModelParams.(CurrentTrainingCondition).(CurrentComparment).Threshold' / ConcentrationScalingFactor, ...
                 ft, opts);  
        ModelParams.(CurrentTrainingCondition).(CurrentComparment).K = cf.K * ConcentrationScalingFactor;         
    end
end    
% With Bootstrap
tic
NumOfIterations = size(ModelParams.Naive.Soma.Bootstrap.Threshold.Mat,1);
for t_ind = 1:length(TrainingConditions)
    CurrentTrainingCondition = TrainingConditions{t_ind}
    for c_ind = 1:length(CompartmentFieldNames)
        CurrentComparment = CompartmentFieldNames{c_ind};    
        InitialConcentrations = ModelParams.(CurrentTrainingCondition).(CurrentComparment).InitialConcetrations;
        for it = 1:NumOfIterations
            Thresholds = ModelParams.(CurrentTrainingCondition).(CurrentComparment).Bootstrap.Threshold.Mat(it,:);
            cf = fit(InitialConcentrations'/ConcentrationScalingFactor, Thresholds'/ConcentrationScalingFactor, ft, opts);  
            ModelParams.(CurrentTrainingCondition).(CurrentComparment).Bootstrap.K.Vector(it) = cf.K * ConcentrationScalingFactor;   % in microM      
        end
                
        CurrentStructure         = ModelParams.(CurrentTrainingCondition).(CurrentComparment).Bootstrap.K;
        Vector                   = CurrentStructure.Vector;        
        CurrentStructure.Mean    = mean(Vector);
        CurrentStructure.Median  = median(Vector);
        CurrentStructure.STD     = std(Vector,1);
        ModelParams.(CurrentTrainingCondition).(CurrentComparment).Bootstrap.K = CurrentStructure;                
    end
end    
toc
warning('on','curvefit:fit:noStartPoint')

%% Calculate K for Egl4 mutants- based on one available Threshold
TrainingConditions = {'CX17255','Egl4','Egl4_ad450'};
ModelParams = Find_K_FromSingleT_WithBootstrap(ModelParams, TrainingConditions, CompartmentFieldNames);

save('H:\BootstrapThresholdForStepResponse_Combined.mat','Data','ModelParams','-v7.3')

disp([datestr(now),'  --  Done boostrapping. Number of iterations = ',num2str(NumOfIterations)])


return

function [Data, ModelParams] = CalculateThresholdFor50PercentResponse_CalculateK (Data, RelevantStrainIndices)

StrainNames           = fieldnames(Data); 
CompartmentFieldNames = {'Soma','Dendrite','Axon'};

% Fit and Fit Options. Concentration units = Dilution!!
ConcentrationScalingFactor = 11.16*1e6;  % Concentration units in microM
ft_T_and_n         = fittype('T^n/(T^n+C^n) *100', 'dependent',{'y'},'independent',{'C'},'coefficients',{'T','n'});
opts               = fitoptions( 'Method', 'NonlinearLeastSquares','robust','Bisquare' );
opts.Display       = 'Off';
opts.MaxFunEvals   = 10000;
opts.MaxIter       = 10000;
opts.TolFun        = 1e-14;
opts.TolX          = 1e-14;
opts.DiffMinChange = 1e-10;
tic
for s_ind = 1:length(RelevantStrainIndices)
    CurrentStrainIndex   = RelevantStrainIndices(s_ind);
    CurrentStrainName    = StrainNames{CurrentStrainIndex}    
    InitialConcentration = Data.(CurrentStrainName).PulsesInfo.InitialConcentration;
    Data.(CurrentStrainName).Stats.InitialConcentration = InitialConcentration;        % Units of microM
    InitialConcentration_Dilution = InitialConcentration /ConcentrationScalingFactor;  % Units of Dilution
    
    for c_ind = 1:length(CompartmentFieldNames)        
        CurrentComparment     = CompartmentFieldNames{c_ind};
        CurrentField          = [CurrentComparment,'_Responsiveness']; 
        CurrentResponsiveness = Data.(CurrentStrainName).Stats.(CurrentField).Mean;
        FinalConcentrations   = Data.(CurrentStrainName).Stats.([CurrentComparment,'_Information']).FinalConcentrations;    
        FinalConcentrations_Dilution = FinalConcentrations / ConcentrationScalingFactor;  % Units of Dilution
        [CurrentThreshold, Current_n, Tguess] = CalculateThreshold_MM(double(FinalConcentrations_Dilution), double(CurrentResponsiveness), double(InitialConcentration_Dilution), ft_T_and_n, opts);  % using all data, without bootstraping
      
        Data.(CurrentStrainName).Stats.([CurrentComparment,'_InitialThresholdGuess']) = Tguess           * ConcentrationScalingFactor; % Units of microM         
        Data.(CurrentStrainName).Stats.([CurrentComparment,'_Threshold'])             = CurrentThreshold * ConcentrationScalingFactor; % Units of microM
        Data.(CurrentStrainName).Stats.([CurrentComparment,'_n'])                     = Current_n;            
    end
end
toc

%% Generate ModelParams
ListOfStrains            = {'str2_CX17256_FromD7_Fast','str2_CX17256_FromD6_Fast','str2_CX17256_FromD5_Fast','str2_CX17256_FromD4_Fast'};
CurrentTrainingCondition = 'Naive';
ModelParams.(CurrentTrainingCondition) = GenerateThresholdStructure(Data, ListOfStrains, CompartmentFieldNames);
ListOfStrains            = {'str2_CX17256_Desensitized_FromD6_Fast','str2_CX17256_Desensitized_FromD5_Fast','str2_CX17256_Desensitized_FromD4_Fast','str2_CX17256_Desensitized_FromD3_Fast'};
CurrentTrainingCondition = 'Desensitized';
ModelParams.(CurrentTrainingCondition) = GenerateThresholdStructure(Data, ListOfStrains, CompartmentFieldNames);

%% Calculate K using MM or using exponential function
%%%%%%% Initialization of exponential fitting function: T = K*(1-exp(-x/K)) %%%%%
% Concentration units = Dilution!!
ft    = fittype('K*(1-exp(-x/K))','dependent',{'y'},'independent',{'x'},'coefficients',{'K'});
opts  = fitoptions( 'Method', 'NonlinearLeastSquares','robust','Bisquare' );
opts.MaxFunEvals = 10000;
opts.MaxIter     = 10000;
opts.TolFun      = 1e-17;
opts.TolX        = 1e-17;
opts.DiffMinChange = 1e-9;
opts.DiffMaxChange = 1e-3;
% opts.StartPoint  = [1e-6 1 ]; 
opts.Upper       = 1e-3 ; 
opts.Lower       = 1e-9 ; 

warning('off','curvefit:fit:noStartPoint')
TrainingConditions = {'Naive','Desensitized'};
for t_ind = 1:length(TrainingConditions)
    CurrentTrainingCondition = TrainingConditions{t_ind};
    for c_ind = 1:length(CompartmentFieldNames)
        CurrentComparment = CompartmentFieldNames{c_ind};    
        cf = fit(ModelParams.(CurrentTrainingCondition).(CurrentComparment).InitialConcetrations'/ ConcentrationScalingFactor, ...
                 ModelParams.(CurrentTrainingCondition).(CurrentComparment).Threshold'/ ConcentrationScalingFactor, ...
                 ft, opts);  
        ModelParams.(CurrentTrainingCondition).(CurrentComparment).K = cf.K * ConcentrationScalingFactor;  % Units of microM       
    end
end    
warning('on','curvefit:fit:noStartPoint')
toc
return

function ModelParams = Find_K_FromSingleT_WithBootstrap(ModelParams, TrainingConditions, CompartmentFieldNames)
% (T/K) = 1-exp(-C/K)
% Calculations in units of microM
ConcentrationScalingFactor = 11.16*1e6;  % Concentration units in microM
K_vec = ([0.01:0.005:0.095   0.1:0.01:1.98 2:0.02:2.98 3:0.05:3.95 4:0.1:7  7.5:0.5:10])*1e-6 * ConcentrationScalingFactor;

for t_ind = 1:length(TrainingConditions)
    CurrentTrainingCondition = TrainingConditions{t_ind};
    for c_ind = 1:length(CompartmentFieldNames)
        CurrentComparment = CompartmentFieldNames{c_ind};    
        C    = ModelParams.(CurrentTrainingCondition).(CurrentComparment).InitialConcetrations; 
        T    = ModelParams.(CurrentTrainingCondition).(CurrentComparment).Threshold; 
        LHS  = T./K_vec;
        RHS  = 1-exp(-C./K_vec);
        BestK_Index = find(LHS<=RHS,1,'first');
        K           = K_vec(BestK_Index);
        
%         figure; plot(K_vec,LHS,'b.-'); hold on; plot(K_vec,RHS,'r.-'); plot(K,T/K,'ko','markerfacecolor','k'); ylim([0 1])
        
        ModelParams.(CurrentTrainingCondition).(CurrentComparment).K = K;
        
        % For the bootstrapped T:
        Tvec      = ModelParams.(CurrentTrainingCondition).(CurrentComparment).Bootstrap.Threshold.Mat; 
        BestK_vec = ones(1,length(Tvec))*NaN;
        for T_ind = 1:length(Tvec)
            T    = Tvec(T_ind);
            LHS  = T./K_vec;
            RHS  = 1-exp(-C./K_vec);
            BestK_Index = find(LHS<=RHS,1,'first');
            if isempty(BestK_Index)                
                continue
            end
            K           = K_vec(BestK_Index);
            % assign to output
            BestK_vec(T_ind) = K;
              
        end
        ModelParams.(CurrentTrainingCondition).(CurrentComparment).Bootstrap.K.Vector = BestK_vec;
        ModelParams.(CurrentTrainingCondition).(CurrentComparment).Bootstrap.K.Mean   = nanmean(BestK_vec);
        ModelParams.(CurrentTrainingCondition).(CurrentComparment).Bootstrap.K.Median = nanmedian(BestK_vec);
        ModelParams.(CurrentTrainingCondition).(CurrentComparment).Bootstrap.K.STD    = nanstd(BestK_vec,1);        
    end
end

return

function StructureOut = GenerateThresholdStructure(Data, ListOfStrains, CompartmentFieldNames)

for c_ind = 1:length(CompartmentFieldNames)
    CurrentComparment = CompartmentFieldNames{c_ind};
    Threshold_InitialGuess_vec = zeros(1,length(ListOfStrains))*NaN;
    Threshold_vec              = zeros(1,length(ListOfStrains))*NaN;
    n_vec                      = zeros(1,length(ListOfStrains))*NaN;
    InitialConcetrations       = zeros(1,length(ListOfStrains))*NaN;
    for s_ind = 1:length(ListOfStrains)
        CurrentStrain = ListOfStrains{s_ind};
        Threshold_InitialGuess_vec(s_ind) = Data.(CurrentStrain).Stats.([CurrentComparment,'_InitialThresholdGuess']);
        Threshold_vec(s_ind)              = Data.(CurrentStrain).Stats.([CurrentComparment,'_Threshold']);
        n_vec(s_ind)                      = Data.(CurrentStrain).Stats.([CurrentComparment,'_n']);
        InitialConcetrations(s_ind)       = Data.(CurrentStrain).Stats.InitialConcentration;
    end
    StructureOut.(CurrentComparment).ListOfStrains          = ListOfStrains;                
    StructureOut.(CurrentComparment).InitialConcetrations   = InitialConcetrations;                            
    StructureOut.(CurrentComparment).Threshold_InitialGuess = Threshold_InitialGuess_vec;                            
    StructureOut.(CurrentComparment).Threshold              = Threshold_vec;                            
    StructureOut.(CurrentComparment).n                      = n_vec;                
end    

return

function StructureOut = GenerateThresholdStructure_WithBootstrap(Data, ListOfStrains, CompartmentFieldNames)

for c_ind = 1:length(CompartmentFieldNames)
    CurrentComparment = CompartmentFieldNames{c_ind};
    Threshold_InitialGuess_vec = zeros(1,length(ListOfStrains))*NaN;
    Threshold_vec              = zeros(1,length(ListOfStrains))*NaN;
    n_vec                      = zeros(1,length(ListOfStrains))*NaN;
    InitialConcetrations       = zeros(1,length(ListOfStrains))*NaN;
    for s_ind = 1:length(ListOfStrains)
        CurrentStrain = ListOfStrains{s_ind};
        Threshold_InitialGuess_vec(s_ind) = Data.(CurrentStrain).Stats.([CurrentComparment,'_InitialThresholdGuess']);
        Threshold_vec(s_ind)              = Data.(CurrentStrain).Stats.([CurrentComparment,'_Threshold']);
        n_vec(s_ind)                      = Data.(CurrentStrain).Stats.([CurrentComparment,'_n']);
        InitialConcetrations(s_ind)       = Data.(CurrentStrain).Stats.InitialConcentration;
    end
    StructureOut.(CurrentComparment).ListOfStrains          = ListOfStrains;                
    StructureOut.(CurrentComparment).InitialConcetrations   = InitialConcetrations;                            
    StructureOut.(CurrentComparment).Threshold_InitialGuess = Threshold_InitialGuess_vec;                            
    StructureOut.(CurrentComparment).Threshold              = Threshold_vec;                            
    StructureOut.(CurrentComparment).n                      = n_vec;          
    
    % Bootstrap        
    LEN = length(Data.(CurrentStrain).Stats.Bootstrap.([CurrentComparment,'_InitialThresholdGuess']).Vector);
    Threshold_InitialGuess_mat = zeros(LEN,length(ListOfStrains))*NaN;
    Threshold_mat              = zeros(LEN,length(ListOfStrains))*NaN;
    n_mat                      = zeros(LEN,length(ListOfStrains))*NaN;
    InitialConcetrations       = zeros(1,length(ListOfStrains))*NaN;
    MeanNumOfRepeats           = zeros(1,length(ListOfStrains))*NaN;
    for s_ind = 1:length(ListOfStrains)
        CurrentStrain = ListOfStrains{s_ind};
        Threshold_InitialGuess_mat(:,s_ind) = Data.(CurrentStrain).Stats.Bootstrap.([CurrentComparment,'_InitialThresholdGuess']).Vector;
        Threshold_mat(:,s_ind)              = Data.(CurrentStrain).Stats.Bootstrap.([CurrentComparment,'_Threshold']).Vector;
        n_mat(:,s_ind)                      = Data.(CurrentStrain).Stats.Bootstrap.([CurrentComparment,'_n']).Vector;
        InitialConcetrations(s_ind)         = Data.(CurrentStrain).Stats.InitialConcentration;
        MeanNumOfRepeats(s_ind)             = Data.(CurrentStrain).Stats.Bootstrap.([CurrentComparment,'_Threshold']).MeanNumOfRepeats;
    end
    StructureOut.(CurrentComparment).InitialConcetrations                    = InitialConcetrations;  
    StructureOut.(CurrentComparment).Bootstrap.MeanNumOfRepeats              = MeanNumOfRepeats;                            
    
    StructureOut.(CurrentComparment).Bootstrap.Threshold_InitialGuess.Mat    = Threshold_InitialGuess_mat;                            
    StructureOut.(CurrentComparment).Bootstrap.Threshold.Mat                 = Threshold_mat;                            
    StructureOut.(CurrentComparment).Bootstrap.n.Mat                         = n_mat;         
    StructureOut.(CurrentComparment).Bootstrap.Threshold_InitialGuess.Mean   = mean(Threshold_InitialGuess_mat,1);                            
    StructureOut.(CurrentComparment).Bootstrap.Threshold.Mean                = mean(Threshold_mat,1);                            
    StructureOut.(CurrentComparment).Bootstrap.n.Mean                        = mean(n_mat,1);          
    StructureOut.(CurrentComparment).Bootstrap.Threshold_InitialGuess.Median = median(Threshold_InitialGuess_mat,1);                            
    StructureOut.(CurrentComparment).Bootstrap.Threshold.Median              = median(Threshold_mat,1);                            
    StructureOut.(CurrentComparment).Bootstrap.n.Median                      = median(n_mat,1);          
    StructureOut.(CurrentComparment).Bootstrap.Threshold_InitialGuess.STD    = std(Threshold_InitialGuess_mat,1,1);                            
    StructureOut.(CurrentComparment).Bootstrap.Threshold.STD                 = std(Threshold_mat,1,1);                            
    StructureOut.(CurrentComparment).Bootstrap.n.STD                         = std(n_mat,1,1);          
    StructureOut.(CurrentComparment).Bootstrap.Threshold_InitialGuess.SEM    = std(Threshold_InitialGuess_mat,1,1) ./ sqrt(MeanNumOfRepeats);                            
    StructureOut.(CurrentComparment).Bootstrap.Threshold.SEM                 = std(Threshold_mat,1,1)./ sqrt(MeanNumOfRepeats);                            
    StructureOut.(CurrentComparment).Bootstrap.n.SEM                         = std(n_mat,1,1)./ sqrt(MeanNumOfRepeats);           
end    

return

function [Threshold, n, Tguess] = CalculateThreshold_MM(UniqueFinalConcentrations, ResponsivenessVec, InitialC, ft_T_and_n, opts)

% Find 50% Value by linear interpolation for fitting starting point  
Tguess          = CalculateThresholdByLinearInterpolation(UniqueFinalConcentrations, ResponsivenessVec, InitialC); 

% Fit  
opts.StartPoint = [Tguess 2];
opts.Lower      = [max([UniqueFinalConcentrations(1)   Tguess/2]) 1 ];
opts.Upper      = [min([UniqueFinalConcentrations(end) Tguess*2]) 30];

% [cf,gof,output] = fit(UniqueFinalConcentrations',ValuesByFinalConcentration,ft_T_and_n, opts);
cf        = fit(UniqueFinalConcentrations',ResponsivenessVec',ft_T_and_n, opts);
Threshold = cf.T;
n         = cf.n;

return

function Threshold = CalculateThresholdByLinearInterpolation(FinalConcentrations, ResponsivenessVec, InitialC)

CurrentUpperRange = 100; % upper range limit = 100% responders

%% Interpolate missing final concentration values
LongConcentrationsVec = ones(1,(length(FinalConcentrations)-1)*100+1);                   
LongValuesVec         = ones(1,(length(FinalConcentrations)-1)*100+1);                   
first_ind = 1;                    
for conc_ind = 1:length(FinalConcentrations)-1
    LongConcentrationsVec(first_ind:first_ind+100) = linspace(FinalConcentrations(conc_ind),FinalConcentrations(conc_ind+1),101) ;
    LongValuesVec(first_ind:first_ind+100)         = linspace(ResponsivenessVec(conc_ind),ResponsivenessVec(conc_ind+1),101) ;
    first_ind = first_ind+100;
end

%% Find 50% responsiveness
% [~, ClosestPointIndex] = min(abs(LongValuesVec-CurrentUpperRange/2));  % Not robust: this is a good strategy only if there is no double-peak   
[~, ClosestPointIndex] = find(LongValuesVec<CurrentUpperRange/2,1,'first');
if isempty(ClosestPointIndex)
    ClosestPointIndex = length(LongValuesVec);
end
ClosestPointValue      = LongValuesVec(ClosestPointIndex);  
HalfRangeConcentration = LongConcentrationsVec(ClosestPointIndex);    
if (abs(ClosestPointValue/CurrentUpperRange- 1/2) > 0.1) 
    % If half range was not found. As an initial guess assume the half range is between the largest final concentration tested and the initial concentration   
    HalfRangeConcentration = (FinalConcentrations(end)+InitialC)/2;
end                        
Threshold = HalfRangeConcentration;

return

% Calculate adaptation times (tau)
function CalculateACTModelAdaptationTimeAndPredictions_WithBootstrap(Data, ModelParams)
% Compute adaptation times with bootstrap
% Compute predictions, including responsiveness and latency

tic
%% WT experiments
ComputeSpecialFits   = true;
CurrentStrain_Slow   = 'str2_CX17256_FromD6_Slow';
disp([datestr(now),'  --  calculating ''ACT_CX17256'''])
[ACT_CX17256, ~, ModelParams.Naive]   = ComputeAdaptationTime_Bootstrap (Data, CurrentStrain_Slow, ModelParams.Naive, ComputeSpecialFits); 

CurrentStrain_Slow    = 'str2_CX17256_ToBuffer_Slow';
CurrentModelParams    = ModelParams.Naive;
ReferenceACTStructure = ACT_CX17256;      % computed from CurrentStrain_Slow= 'str2_CX17256_FromD6_Slow';
ReferenceACTfitFields = {'FitBasedOnResponsivenessAndLatency_UseDataSTD','FitBasedOnOnePulseResponsivenessAndLatency','FitBasedOnTwoRepeatPerCondition',...
                         'FitBasedOnAllResponsiveness','FitBasedOnOnePulseResponsiveness'};
ACT_VaryInitialC = CalculateACTpredictionsWithGivenModelParams (Data, CurrentStrain_Slow, CurrentModelParams, ReferenceACTStructure, ReferenceACTfitFields);  

%% egl-4 experiments and their control
ComputeSpecialFits   = false;
CurrentStrain_Slow = 'str2_CX17255_AllSlow';
disp([datestr(now),'  --  calculating ''ACT_CX17255_FastCX17255'''])
[ACT_CX17255, ~, ModelParams.CX17255] = ComputeAdaptationTime_Bootstrap (Data, CurrentStrain_Slow, ModelParams.CX17255, ComputeSpecialFits);

CurrentStrain_Slow = 'str2_Egl4_AllSlow';
disp([datestr(now),'  --  calculating ''ACT_Egl4'''])
[ACT_Egl4, ~, ModelParams.Egl4]       = ComputeAdaptationTime_Bootstrap (Data, CurrentStrain_Slow, ModelParams.Egl4, ComputeSpecialFits);

CurrentStrain_Slow = 'str2_Egl4_ad450_AllSlow';
disp([datestr(now),'  --  calculating ''ACT_Egl4_ad450'''])
[ACT_Egl4_ad450, ~, ModelParams.Egl4_ad450]    = ComputeAdaptationTime_Bootstrap (Data, CurrentStrain_Slow, ModelParams.Egl4_ad450, ComputeSpecialFits);

%% save
disp([datestr(now),'  --  saving file'])
save('H:\ProcessedData_ACT_ModelParams.mat','Data','ModelParams','ACT_CX17256','ACT_CX17255','ACT_Egl4','ACT_Egl4_ad450','ACT_VaryInitialC','-v7.3');
disp([datestr(now),'  --  File saved'])

toc
return

function [ACT, MeasuredStats, CurrentModelParams_Output] = ComputeAdaptationTime_Bootstrap (Data, CurrentStrain_Slow, CurrentModelParams, ComputeSpecialFits)
% K and Beta=1/AdadaptationTime are parameters used for the adaptive concentration threshold model.   

%% Example for Input arguments examples:
% CurrentStrain_Slow    = 'str2_CX17256_FromD6_Slow';
% CurrentModelParams    = ModelParams.Naive;
% ComputeSpecialFits    = true;
if ~exist('ComputeSpecialFits','var')
    ComputeSpecialFits = false;
end
ComputeSpecialFitsForAllCompartments = false;  

%% Initialize structures
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
MaxNumberOfRepeats    = 40;  % free parameter. It just needs to be higher than the maximal number of repeats
DyeInformation        = Data.(CurrentStrain_Slow).DyeInformation;
PulsesNames           = Data.(CurrentStrain_Slow).PulsesInfo.PulseFromButanone;
NumberOfPulses        = length(PulsesNames);            

AdaptationTimeVec = [0.01 0.05 0.1:0.05:0.95 1:0.02:2 2.2:0.2:12.8 13:0.1:19 19.2:0.2:21 22:1:30 32:2:40 45:5:70 80:10:120 150:50:1e4 1e5]; % For screen
BetaVec           = 1./AdaptationTimeVec;
C                 = DyeInformation.C_Smoothed;    
dCdt              = DyeInformation.dCdt_Smoothed;
DeltaCoverC       = DyeInformation.DeltaCoverC;
d2Cdt2            = DyeInformation.d2Cdt2;
NaN_Mat           = all(isnan(C),3); % true where (p_ind,rep_ind) is empty (only NaNs)
TimeVec           = DyeInformation.TimeVectorInSec;
DyeVectorLength   = length(TimeVec);

% Find Measured stats for each compartment and condition (average over repeats)  
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};   
    MeasuredStats.(CurrentCompartment).Responsiveness = Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Responsiveness']);
    MeasuredStats.(CurrentCompartment).Latency        = Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Latency_InSec']);
    MeasuredStats.(CurrentCompartment).Latency_MA     = Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Latency_MA_InSec']);
    MeasuredStats.(CurrentCompartment).C_AtNeuronActivation.Matrix           = DyeInformation.Stats.C_AtNeuronActivation.(CurrentCompartment);
    MeasuredStats.(CurrentCompartment).dCdt_AtNeuronActivation.Matrix        = DyeInformation.Stats.dCdt_AtNeuronActivation.(CurrentCompartment);
    MeasuredStats.(CurrentCompartment).d2Cdt2_AtNeuronActivation.Matrix      = DyeInformation.Stats.d2Cdt2_AtNeuronActivation.(CurrentCompartment);
    MeasuredStats.(CurrentCompartment).DeltaCoverC_AtNeuronActivation.Matrix = DyeInformation.Stats.DeltaCoverC_AtNeuronActivation.(CurrentCompartment);
end
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};   
    MeasuredStats.(CurrentCompartment).Responsiveness.NumberOfRepeats = sum(~isnan(Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Responsiveness']).Matrix),2)';  
    MeasuredStats.(CurrentCompartment).Latency.NumberOfRepeats        = sum(~isnan(Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Latency_InSec']).Matrix),2)';  
    MeasuredStats.(CurrentCompartment).Latency_MA.NumberOfRepeats     = sum(~isnan(Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Latency_MA_InSec']).Matrix),2)'; 
end

% Initialization
ACT.Information.CompartmentFieldNames = CompartmentFieldNames;
ACT.Information.StrainName            = CurrentStrain_Slow;
ACT.Information.ModelParams           = CurrentModelParams;
ACT.Information.MeasuredStats         = MeasuredStats;
ACT.Information.InitialConcentrationPerPulse = DyeInformation.InitialConcentrationPerPulse;
ACT.Information.FinalConcentrationPerPulse   = DyeInformation.FinalConcentrationPerPulse;
for c_ind = 1:length(CompartmentFieldNames)        
    CurrentCompartment = CompartmentFieldNames{c_ind};
    ACT.(CurrentCompartment).K                                  = CurrentModelParams.(CurrentCompartment).K;
    ACT.(CurrentCompartment).AdaptationTimeVec                  = AdaptationTimeVec;     
    ACT.(CurrentCompartment).StatsPerPulse.FinalConcentration   = DyeInformation.FinalConcentrationPerPulse;     
    ACT.(CurrentCompartment).StatsPerPulse.InitialConcentration = DyeInformation.InitialConcentrationPerPulse;    
    ACT.(CurrentCompartment).StatsPerPulse.NumberOfRepeats      = Data.(CurrentStrain_Slow).Stats.Soma_Information.NumberOfRepeats;
    ACT.(CurrentCompartment).ACT                       = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,DyeVectorLength,'single')*NaN; 
    ACT.(CurrentCompartment).Responsiveness            = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
    ACT.(CurrentCompartment).Latency                   = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
    ACT.(CurrentCompartment).C_at_Activation           = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
    ACT.(CurrentCompartment).dCdt_at_Activation        = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
    ACT.(CurrentCompartment).d2Cdt2_at_Activation      = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
    ACT.(CurrentCompartment).DeltaCoverC_at_Activation = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
end

%% Find adaptive threshold for each compartment and condition
for p_ind = 1:NumberOfPulses         
%     disp([datestr(now), '  --  ',PulsesNames{p_ind}])
    for rep_ind = 1:MaxNumberOfRepeats            
        if NaN_Mat(p_ind,rep_ind)
            % no more traces for this pulse
            break
        end
        Current_C           = double(squeeze(C(p_ind,rep_ind,:)));
        Current_dCdt        = double(squeeze(dCdt(p_ind,rep_ind,:)));
        Current_d2Cdt2      = double(squeeze(d2Cdt2(p_ind,rep_ind,:)));
        Current_DeltaCoverC = double(squeeze(DeltaCoverC(p_ind,rep_ind,:)));
        
        for c_ind = 1:length(CompartmentFieldNames)
            CurrentCompartment = CompartmentFieldNames{c_ind};     
            K = ACT.(CurrentCompartment).K;
        
            for B_ind = 1:length(BetaVec)                
                Beta = BetaVec(B_ind);

                %% Loop over Beta + Compare to actual initiation time
                [AdaptiveThresholdVector, PredictedActivationIndex, PredictedActivationTime] = CalculateAdaptiveThresholdFromDyePattern_ExpSaturation_03 (TimeVec, Current_C, K, Beta);
                ACT.(CurrentCompartment).ACT(B_ind,p_ind,rep_ind,1:DyeVectorLength)          = AdaptiveThresholdVector;

                if ~isempty(PredictedActivationTime)                    
                    ACT.(CurrentCompartment).Responsiveness(B_ind,p_ind,rep_ind)            = 100;       
                    ACT.(CurrentCompartment).Latency(B_ind,p_ind,rep_ind)                   = PredictedActivationTime;    
                    ACT.(CurrentCompartment).C_at_Activation(B_ind,p_ind,rep_ind)           = Current_C(PredictedActivationIndex)   ;    
                    ACT.(CurrentCompartment).dCdt_at_Activation(B_ind,p_ind,rep_ind)        = Current_dCdt(PredictedActivationIndex);    
                    ACT.(CurrentCompartment).d2Cdt2_at_Activation(B_ind,p_ind,rep_ind)      = Current_d2Cdt2(PredictedActivationIndex);    
                    ACT.(CurrentCompartment).DeltaCoverC_at_Activation(B_ind,p_ind,rep_ind) = Current_DeltaCoverC(PredictedActivationIndex);    
                else            
                    ACT.(CurrentCompartment).Responsiveness(B_ind,p_ind,rep_ind)            = 0;   
                end
            end
        end                
    end
end
disp([datestr(now), '  --  ','Finished ACT screening'])           
     
% Find ACT-predicted stats for each compartment and condition (average over repeats)  
FieldsForStats = {'Responsiveness','Latency','C_at_Activation','dCdt_at_Activation','d2Cdt2_at_Activation','DeltaCoverC_at_Activation'};
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};   
    for f_ind = 1:length(FieldsForStats)
        CurrentField = FieldsForStats{f_ind};
        ACT.(CurrentCompartment).StatsPerPulse.(CurrentField).Mean   = nanmean(ACT.(CurrentCompartment).(CurrentField),3);
        ACT.(CurrentCompartment).StatsPerPulse.(CurrentField).Median = nanmedian(ACT.(CurrentCompartment).(CurrentField),3);
        ACT.(CurrentCompartment).StatsPerPulse.(CurrentField).STD    = nanstd(ACT.(CurrentCompartment).(CurrentField),[],3);
        ACT.(CurrentCompartment).StatsPerPulse.(CurrentField).SEM    = ACT.(CurrentCompartment).StatsPerPulse.(CurrentField).STD ./ ...
                                                           sqrt(repmat(ACT.(CurrentCompartment).StatsPerPulse.NumberOfRepeats',length(BetaVec),1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Find Best fit based on different criteria and different training sets (portions of data) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% NO BOOTSTRAP - USE ALL TRACES.
disp([datestr(now),'  --  Adaptation time calculation, without bootstrapping']);
%%%%  Initialization  %%%%
NumberOfIterations = 1; it=1;
MaxNumberOfPulses  = size(MeasuredStats.Soma.Responsiveness.Matrix,1);

InitializationStructure.AdaptationTimeVec          = AdaptationTimeVec;
InitializationStructure.BestFitIndexInVec          = ones(1,NumberOfIterations,'single')*NaN;       
InitializationStructure.LowerLimitIndex            = ones(1,NumberOfIterations,'single')*NaN;       
InitializationStructure.UpperLimitIndex            = ones(1,NumberOfIterations,'single')*NaN;       
InitializationStructure.BestFitAdaptationTime      = ones(1,NumberOfIterations,'single')*NaN;  
InitializationStructure.LowerLimitAdaptationTime   = ones(1,NumberOfIterations,'single')*NaN;  
InitializationStructure.UpperLimitAdaptationTime   = ones(1,NumberOfIterations,'single')*NaN;  
InitializationStructure.ReferencePulseIndex        = ones(1,NumberOfIterations,'single')*NaN;  
InitializationStructure.ReferencePulseIndex2       = ones(1,NumberOfIterations,'single')*NaN;  
InitializationStructure.ReferencePulseName         = cell(1,NumberOfIterations);       
InitializationStructure.MeanResponsiveness         = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       
InitializationStructure.MeanLatency                = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       
InitializationStructure.ReferenceResponsiveness    = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       
InitializationStructure.ReferenceResponsivenessSEM = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       
InitializationStructure.ReferenceLatency           = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       
InitializationStructure.ReferenceLatencySEM        = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       

%%%%  Loop  %%%%
for c_ind = 1:length(CompartmentFieldNames)       
    CurrentCompartment = CompartmentFieldNames{c_ind}; 
    if ComputeSpecialFitsForAllCompartments || c_ind==1            
        CurrentComputeSpecialFits = ComputeSpecialFits;
    else
        CurrentComputeSpecialFits = false;
    end
    
    %% Get Data and model predictions for different adaptation times
    ReferenceResponsivenessMatrix = MeasuredStats.(CurrentCompartment).Responsiveness.Matrix;    
    ReferenceLatencyMatrix        = MeasuredStats.(CurrentCompartment).Latency.Matrix;   
    ResponsivenessNumberOfRepeats = MeasuredStats.(CurrentCompartment).Responsiveness.NumberOfRepeats;  
    LatencyNumberOfRepeats        = MeasuredStats.(CurrentCompartment).Latency.NumberOfRepeats;    
    
    ReferenceResponsivenessMean   = nanmean(ReferenceResponsivenessMatrix,2)';
    ReferenceResponsivenessSEM    = nanstd(ReferenceResponsivenessMatrix,[],2)' ./ sqrt(ResponsivenessNumberOfRepeats);
    ReferenceLatencyMean          = nanmean(ReferenceLatencyMatrix,2)';
    ReferenceLatencyMean(LatencyNumberOfRepeats<=1) = NaN;                   % Do not rely on data that comes from ONE SINGLE data point 
    ReferenceLatencySEM           = nanstd(ReferenceLatencyMatrix,[],2)' ./ sqrt(LatencyNumberOfRepeats);
    ReferenceResponsivenessMat          = repmat(ReferenceResponsivenessMean,length(BetaVec),1);
    ReferenceResponsivenessSEM_Mat      = repmat(ReferenceResponsivenessSEM,length(BetaVec),1);
    ResponsivenessNumberOfRepeats_Mat   = repmat(ResponsivenessNumberOfRepeats,length(BetaVec),1);
    ReferenceLatencyMat                 = repmat(ReferenceLatencyMean,length(BetaVec),1);    
    
    PredictedResponsiveness_Mat    = ACT.(CurrentCompartment).Responsiveness;
    PredictedLatencies_Mat         = ACT.(CurrentCompartment).Latency;        
    PredictedResponsivenessMeanMat = nanmean(PredictedResponsiveness_Mat,3);
    PredictedLatencyMeanMat        = nanmean(PredictedLatencies_Mat,3);

    
    %% Structures initialization
    InitializationStructure.ResponsivenessNumberOfRepeats  = ResponsivenessNumberOfRepeats;
    InitializationStructure.LatencyNumberOfRepeats         = LatencyNumberOfRepeats;
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD = InitializationStructure;  % Use responsiveness and Latency. UseDataSTD  
    if CurrentComputeSpecialFits      
        ACT.(CurrentCompartment).FitBasedOnAllResponsiveness                   = InitializationStructure;  % Use responsiveness 
        ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness              = InitializationStructure;  % Use responsiveness in one pulse  
        ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency    = InitializationStructure;  % Use responsiveness and Latency in one pulse  
    end
       
    %%      FitBasedOnResponsivenessAndLatency_UseDataSTD
    %%%%%%  Compute "cost" = error between prediction and data        
    DistanceBetweenPredictionAndMeasurement_Mat = (PredictedResponsivenessMeanMat - ReferenceResponsivenessMat).^2;              % SE
    IndicesWithNoDistanceAllowed                = (ReferenceResponsivenessSEM_Mat==0)&(ResponsivenessNumberOfRepeats_Mat>1);     % zero standard deviation
    DistanceBetweenPredictionAndMeasurement_Mat(IndicesWithNoDistanceAllowed & (DistanceBetweenPredictionAndMeasurement_Mat>0)) = inf;    
    DistanceBetweenPredictionAndMeasurement     = sqrt(nanmean(DistanceBetweenPredictionAndMeasurement_Mat,2));           % RMSE = RootMeanSquareErrors
    Cost_Responsiveness                         = DistanceBetweenPredictionAndMeasurement; 

    DistanceBetweenPredictionAndMeasurement_Mat = (PredictedLatencyMeanMat- ReferenceLatencyMat).^2;               % SE
    DistanceBetweenPredictionAndMeasurement     = sqrt(nanmean(DistanceBetweenPredictionAndMeasurement_Mat,2));   % RMSE = RootMeanSquareErrors
    Cost_Latency                                = DistanceBetweenPredictionAndMeasurement; 

%     Cost_Total = Cost_Responsiveness + Cost_Latency;  
%     figure; plot(Cost_Responsiveness,'b'); hold on;plot(Cost_Latency,'g');plot(Cost_Total,'r'); %plot(BestFitIndex,Cost_Total(BestFitIndex),'ko')

    % Do not allow "infinite error cost" of either responsiveness or latency   
    Cost_Responsiveness(~isfinite(Cost_Latency))        = NaN;
    Cost_Latency(~isfinite(Cost_Responsiveness))        = NaN;            
    Cost_Responsiveness(~isfinite(Cost_Responsiveness)) = NaN;
    Cost_Latency(~isfinite(Cost_Latency))               = NaN;  

    %%%%%%   Normalize to dimensionless parameters in the range of [0 1]   
    if (nanmax(Cost_Responsiveness)-nanmin(Cost_Responsiveness))>0  % Avoid normalizing by zeros
        Cost_Responsiveness = (Cost_Responsiveness - min(Cost_Responsiveness))./(max(Cost_Responsiveness)-min(Cost_Responsiveness));
    end
    if (nanmax(Cost_Latency)-nanmin(Cost_Latency))>0                % Avoid normalizing by zeros
        Cost_Latency        = (Cost_Latency - min(Cost_Latency))./(max(Cost_Latency)-min(Cost_Latency));    
    end
%     Cost_Total          = Cost_Responsiveness + Cost_Latency;  
%     figure; plot(Cost_Responsiveness,'b'); hold on;plot(Cost_Latency,'g');plot(Cost_Total,'r'); %plot(BestFitIndex,Cost_Total(BestFitIndex),'ko')
%     Cost_Responsiveness_all = Cost_Responsiveness;
%     Cost_Latency_all        = Cost_Latency;       
%     Cost_Total_all          = Cost_Responsiveness + Cost_Latency; % For plot

    %%%%%%   Find range of possible good fits. 
    % 1. Between minima for Responsiveness and Latency costs
    % 2. Reponsiveness cost is not out of range

    Indices_Responsiveness = find(Cost_Responsiveness==nanmin(Cost_Responsiveness));
    Indices_Latency        = find(Cost_Latency==nanmin(Cost_Latency));
    Indices_OutOfRange     = false(size(Cost_Responsiveness));

    if Indices_Responsiveness(1) >= Indices_Latency(end)       % Keep only indices between the two minima
        Indices_OutOfRange(1:Indices_Latency(1)-1)            = true;
        Indices_OutOfRange(Indices_Responsiveness(end)+1:end) = true;
    elseif Indices_Responsiveness(end) <= Indices_Latency(1)
        Indices_OutOfRange(1:Indices_Responsiveness(1)-1)     = true;
        Indices_OutOfRange(Indices_Latency(end)+1:end)        = true;
    else
        % Overlap exist. Anything out of the overlapping area is excluded.
        Indices_OutOfRange  = true(size(Cost_Responsiveness));
        Indices_OutOfRange(intersect(Indices_Responsiveness,Indices_Latency))= false;                
    end    

    ResponsivenessCostOutOfRange  = Cost_Responsiveness - nanmin(Cost_Responsiveness) >  0.5 ;  % Huge error (on a scale of [0,1])
    Indices_OutOfRange            = Indices_OutOfRange | ResponsivenessCostOutOfRange;

    Cost_Responsiveness(Indices_OutOfRange)= NaN;
    Cost_Latency(Indices_OutOfRange)       = NaN;    
    Cost_Total = Cost_Responsiveness + Cost_Latency; 

    [~,BestFitIndex]  = min(Cost_Total);      
    IndicesAroundMin  = Cost_Total <= 1.1 * min(Cost_Total);                
    LowerLimitIndex   = find(IndicesAroundMin,1,'first'); 
    UpperLimitIndex   = find(IndicesAroundMin,1,'last'); 

%         figure; plot(Cost_Responsiveness_all,'b--'); hold on; plot(Cost_Latency_all,'g--');plot(Cost_Total_all,'r--'); 
%         plot(Cost_Responsiveness,'b','linewidth',2); hold on; plot(Cost_Latency,'g','linewidth',2);plot(Cost_Total,'r','linewidth',2); plot(BestFitIndex,Cost_Total(BestFitIndex),'ko')
%         title(AdaptationTimeVec(BestFitIndex))
    %%%%%%  Assign to structure
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.BestFitIndexInVec(it)            = BestFitIndex;
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.LowerLimitIndex(it)              = LowerLimitIndex;
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.UpperLimitIndex(it)              = UpperLimitIndex;
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.LowerLimitAdaptationTime(it)     = AdaptationTimeVec(LowerLimitIndex);
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.UpperLimitAdaptationTime(it)     = AdaptationTimeVec(UpperLimitIndex);
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.BestFitAdaptationTime(it)        = AdaptationTimeVec(BestFitIndex);           
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.MeanResponsiveness(it,:)         = PredictedResponsivenessMeanMat(BestFitIndex,:);
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.MeanLatency(it,:)                = PredictedLatencyMeanMat(BestFitIndex,:);
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.ReferenceResponsiveness(it,:)    = ReferenceResponsivenessMean;
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.ReferenceResponsivenessSEM(it,:) = ReferenceResponsivenessSEM;
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.ReferenceLatency(it,:)           = ReferenceLatencyMean;
    ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.ReferenceLatencySEM(it,:)        = ReferenceLatencySEM;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if ~CurrentComputeSpecialFits
        continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    %%   FitBasedOnAllResponsiveness            
    DistanceBetweenPredictionAndMeasurement_Mat = (PredictedResponsivenessMeanMat - ReferenceResponsivenessMat).^2;     % SE            
    DistanceBetweenPredictionAndMeasurement     = nansum(DistanceBetweenPredictionAndMeasurement_Mat,2);

    LogicalVec               = DistanceBetweenPredictionAndMeasurement == min(DistanceBetweenPredictionAndMeasurement);        
    LowerLimitIndex          = find(LogicalVec,1,'first');
    UpperLimitIndex          = find(LogicalVec,1,'last');
    LowerLimitAdaptationTime = AdaptationTimeVec(LowerLimitIndex);
    UpperLimitAdaptationTime = AdaptationTimeVec(UpperLimitIndex);        
    BestFitAdaptationTime    = (LowerLimitAdaptationTime + UpperLimitAdaptationTime)/2;
    BestFitIndex             = find(AdaptationTimeVec >=  BestFitAdaptationTime, 1, 'first');    

    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.BestFitIndexInVec(it)            = BestFitIndex;
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.LowerLimitIndex(it)              = LowerLimitIndex;
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.UpperLimitIndex(it)              = UpperLimitIndex;
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.LowerLimitAdaptationTime(it)     = AdaptationTimeVec(LowerLimitIndex);
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.UpperLimitAdaptationTime(it)     = AdaptationTimeVec(UpperLimitIndex);
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.BestFitAdaptationTime(it)        = AdaptationTimeVec(BestFitIndex);           
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.MeanResponsiveness(it,:)         = PredictedResponsivenessMeanMat(BestFitIndex,:);
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.MeanLatency(it,:)                = PredictedLatencyMeanMat(BestFitIndex,:);
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.ReferenceResponsiveness(it,:)    = ReferenceResponsivenessMean;
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.ReferenceResponsivenessSEM(it,:) = ReferenceResponsivenessSEM;
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.ReferenceLatency(it,:)           = ReferenceLatencyMean;
    ACT.(CurrentCompartment).FitBasedOnAllResponsiveness.ReferenceLatencySEM(it,:)        = ReferenceLatencySEM;
   
    %%  FitBasedOnOnePulseResponsiveness             
    ReferenceResponsivenessVec                                   = ReferenceResponsivenessMean;    
    ReferenceResponsivenessVec(ResponsivenessNumberOfRepeats<=1) = NaN;    % Do not rely on data that comes from ONE SINGLE data point     
    ReferencePulseIndex                     = find(ReferenceResponsivenessVec > 0, 1, 'last'); % last responsive condition
    if ReferenceResponsivenessVec(ReferencePulseIndex)==100 
        % 100% responsiveness and next pulse is 0% responsiveness. In this case TWO pulses needs to be used (include also ReferencePulseIndex+1)   
        ReferencePulseIndex2 = ReferencePulseIndex + 1;
        DistanceBetweenPredictionAndMeasurement = (PredictedResponsivenessMeanMat(:,ReferencePulseIndex)  - ReferenceResponsivenessVec(ReferencePulseIndex)).^2 + ...
                                                  (PredictedResponsivenessMeanMat(:,ReferencePulseIndex2) - ReferenceResponsivenessVec(ReferencePulseIndex2)).^2;
    else
        % At this pulse: 0 < responsiveness < 100. Sufficient for adaptation time calculation
        ReferencePulseIndex2 = NaN;
        DistanceBetweenPredictionAndMeasurement = (PredictedResponsivenessMeanMat(:,ReferencePulseIndex) - ReferenceResponsivenessVec(ReferencePulseIndex)).^2;
    end

    LogicalVec               = DistanceBetweenPredictionAndMeasurement == min(DistanceBetweenPredictionAndMeasurement);        
    LowerLimitIndex          = find(LogicalVec,1,'first');
    UpperLimitIndex          = find(LogicalVec,1,'last');
    LowerLimitAdaptationTime = AdaptationTimeVec(LowerLimitIndex);
    UpperLimitAdaptationTime = AdaptationTimeVec(UpperLimitIndex);        
    BestFitAdaptationTime    = (LowerLimitAdaptationTime + UpperLimitAdaptationTime)/2;
    BestFitIndex             = find(AdaptationTimeVec >=  BestFitAdaptationTime, 1, 'first');   

    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.ReferencePulseIndex(it)          = ReferencePulseIndex;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.ReferencePulseIndex2(it)         = ReferencePulseIndex2;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.ReferencePulseName{it}           = PulsesNames{ReferencePulseIndex};
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.BestFitIndexInVec(it)            = BestFitIndex;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.LowerLimitIndex(it)              = LowerLimitIndex;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.UpperLimitIndex(it)              = UpperLimitIndex;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.LowerLimitAdaptationTime(it)     = AdaptationTimeVec(LowerLimitIndex);
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.UpperLimitAdaptationTime(it)     = AdaptationTimeVec(UpperLimitIndex);
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.BestFitAdaptationTime(it)        = AdaptationTimeVec(BestFitIndex);           
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.MeanResponsiveness(it,:)         = PredictedResponsivenessMeanMat(BestFitIndex,:);
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.MeanLatency(it,:)                = PredictedLatencyMeanMat(BestFitIndex,:);
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.ReferenceResponsiveness(it,:)    = ReferenceResponsivenessMean;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.ReferenceResponsivenessSEM(it,:) = ReferenceResponsivenessSEM;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.ReferenceLatency(it,:)           = ReferenceLatencyMean;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.ReferenceLatencySEM(it,:)        = ReferenceLatencySEM;

    %%  FitBasedOnOnePulseResponsivenessAndLatency             
    ReferenceResponsivenessVec                                   = ReferenceResponsivenessMean;    
    ReferenceResponsivenessVec(ResponsivenessNumberOfRepeats<=1) = NaN;    % Do not rely on data that comes from ONE SINGLE data point     
    ReferenceLatencyVec                                          = ReferenceLatencyMean;    
    ReferenceLatencyVec(LatencyNumberOfRepeats<=1)               = NaN;    % Do not rely on data that comes from ONE SINGLE data point 

    LogicalVector                             = ReferenceResponsivenessVec > 0;
    LogicalVector(isnan(ReferenceLatencyVec)) = false;                
    ReferencePulseIndex                       = find(LogicalVector, 1, 'last'); % last responsive condition with Latency data
    DistanceBetweenPredictionAndMeasurement   = (PredictedResponsivenessMeanMat(:,ReferencePulseIndex) - ReferenceResponsivenessVec(ReferencePulseIndex)).^2;
    % Allow Flexibility Around BestResponsiveness
    LogicalVec_Responsiveness1                  = DistanceBetweenPredictionAndMeasurement <= prctile(DistanceBetweenPredictionAndMeasurement,1);
    LogicalVec_Responsiveness2                  = DistanceBetweenPredictionAndMeasurement <= 1.1*min(DistanceBetweenPredictionAndMeasurement);
    LogicalVec_Responsiveness                   = LogicalVec_Responsiveness1 | LogicalVec_Responsiveness2;

    DistanceBetweenPredictionAndMeasurement   = (PredictedLatencyMeanMat(:,ReferencePulseIndex) - ReferenceLatencyVec(ReferencePulseIndex)).^2;          
    DistanceBetweenPredictionAndMeasurement(~LogicalVec_Responsiveness) = inf;               % First condition:  must best fit reponsiveness

    LogicalVec               = DistanceBetweenPredictionAndMeasurement == min(DistanceBetweenPredictionAndMeasurement);        
    LowerLimitIndex          = find(LogicalVec,1,'first');
    UpperLimitIndex          = find(LogicalVec,1,'last');
    LowerLimitAdaptationTime = AdaptationTimeVec(LowerLimitIndex);
    UpperLimitAdaptationTime = AdaptationTimeVec(UpperLimitIndex);        
    BestFitAdaptationTime    = (LowerLimitAdaptationTime + UpperLimitAdaptationTime)/2;
    BestFitIndex             = find(AdaptationTimeVec >=  BestFitAdaptationTime, 1, 'first');   

    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.ReferencePulseIndex(it)          = ReferencePulseIndex;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.ReferencePulseName{it}           = PulsesNames{ReferencePulseIndex};
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.BestFitIndexInVec(it)            = BestFitIndex;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.LowerLimitIndex(it)              = LowerLimitIndex;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.UpperLimitIndex(it)              = UpperLimitIndex;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.LowerLimitAdaptationTime(it)     = AdaptationTimeVec(LowerLimitIndex);
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.UpperLimitAdaptationTime(it)     = AdaptationTimeVec(UpperLimitIndex);
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.BestFitAdaptationTime(it)        = AdaptationTimeVec(BestFitIndex);           
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.MeanResponsiveness(it,:)         = PredictedResponsivenessMeanMat(BestFitIndex,:);
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.MeanLatency(it,:)                = PredictedLatencyMeanMat(BestFitIndex,:);
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.ReferenceResponsiveness(it,:)    = ReferenceResponsivenessMean;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.ReferenceResponsivenessSEM(it,:) = ReferenceResponsivenessSEM;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.ReferenceLatency(it,:)           = ReferenceLatencyMean;
    ACT.(CurrentCompartment).FitBasedOnOnePulseResponsivenessAndLatency.ReferenceLatencySEM(it,:)        = ReferenceLatencySEM;                  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% BOOTSTRAP OVER TRACES
%%%%  Initialization  %%%%
rng('default')
NumberOfIterations = 2e3;
MaxNumberOfPulses  = size(MeasuredStats.Soma.Responsiveness.Matrix,1);

InitializationStructure.AdaptationTimeVec          = AdaptationTimeVec;
InitializationStructure.BestFitIndexInVec          = ones(1,NumberOfIterations,'single')*NaN;       
InitializationStructure.LowerLimitIndex            = ones(1,NumberOfIterations,'single')*NaN;       
InitializationStructure.UpperLimitIndex            = ones(1,NumberOfIterations,'single')*NaN;       
InitializationStructure.BestFitAdaptationTime      = ones(1,NumberOfIterations,'single')*NaN;  
InitializationStructure.LowerLimitAdaptationTime   = ones(1,NumberOfIterations,'single')*NaN;  
InitializationStructure.UpperLimitAdaptationTime   = ones(1,NumberOfIterations,'single')*NaN;  
InitializationStructure.ReferencePulseIndex        = ones(1,NumberOfIterations,'single')*NaN;  
InitializationStructure.ReferencePulseIndex2       = ones(1,NumberOfIterations,'single')*NaN;  
InitializationStructure.ReferencePulseName         = cell(1,NumberOfIterations);       
InitializationStructure.MeanResponsiveness         = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       
InitializationStructure.MeanLatency                = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       
InitializationStructure.ReferenceResponsiveness    = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       
InitializationStructure.ReferenceResponsivenessSEM = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       
InitializationStructure.ReferenceLatency           = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       
InitializationStructure.ReferenceLatencySEM        = ones(NumberOfIterations,MaxNumberOfPulses,'single')*NaN;       

%%%%  Loop  %%%%
for c_ind = 1:length(CompartmentFieldNames)       
    CurrentCompartment = CompartmentFieldNames{c_ind};         
    if ComputeSpecialFitsForAllCompartments || c_ind==1            
        CurrentComputeSpecialFits = ComputeSpecialFits;
    else
        CurrentComputeSpecialFits = false;
    end

    disp([datestr(now),'  --  Calculating adaptation time for ',CurrentCompartment,' with boostrapping']);
    
    %% Get Data and model predictions for different adaptation times
    ReferenceResponsivenessMatrix = MeasuredStats.(CurrentCompartment).Responsiveness.Matrix;    
    ReferenceLatencyMatrix        = MeasuredStats.(CurrentCompartment).Latency.Matrix;   
    ResponsivenessNumberOfRepeats = MeasuredStats.(CurrentCompartment).Responsiveness.NumberOfRepeats;  
    LatencyNumberOfRepeats        = MeasuredStats.(CurrentCompartment).Latency.NumberOfRepeats;    

    ReferenceCell.Responsiveness.RepeatsValues                = cell(1,NumberOfPulses);
    ReferenceCell.Responsiveness.RepeatsColumnIndicesInMatrix = cell(1,NumberOfPulses);
    ReferenceCell.Responsiveness.NumberOfRepeats              = ones(1,NumberOfPulses)*NaN;
    ReferenceCell.Latency                                     = ReferenceCell.Responsiveness;
    for p_ind = 1:NumberOfPulses
        Vec             = ReferenceResponsivenessMatrix(p_ind,:); 
        IndicesInMatrix = find(~isnan(Vec));
        Vec             = Vec(IndicesInMatrix);
        ReferenceCell.Responsiveness.RepeatsColumnIndicesInMatrix{p_ind}   = IndicesInMatrix;
        ReferenceCell.Responsiveness.RepeatsValues{p_ind}                  = Vec;
        ReferenceCell.Responsiveness.NumberOfRepeats(p_ind)                = length(Vec);
           
        Vec             = ReferenceLatencyMatrix(p_ind,:); 
        IndicesInMatrix = find(~isnan(Vec));
        Vec             = Vec(IndicesInMatrix);
        ReferenceCell.Latency.RepeatsColumnIndicesInMatrix{p_ind}   = IndicesInMatrix;
        ReferenceCell.Latency.RepeatsValues{p_ind}                  = Vec;
        ReferenceCell.Latency.NumberOfRepeats(p_ind)                = length(Vec);           
    end    
    %% Structures initialization
    InitializationStructure.ResponsivenessNumberOfRepeats  = ResponsivenessNumberOfRepeats;
    InitializationStructure.LatencyNumberOfRepeats         = LatencyNumberOfRepeats;
    ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD = InitializationStructure;  % Use responsiveness and Latency. UseDataSTD  
    
    if CurrentComputeSpecialFits        
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness                   = InitializationStructure;  % Use responsiveness 
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness              = InitializationStructure;  % Use responsiveness in one pulse  
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency    = InitializationStructure;  % Use responsiveness and Latency in one pulse  
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnTwoRepeatPerCondition               = InitializationStructure;  % Use responsiveness and Latency in one pulse          
    end
    
    for it = 1:NumberOfIterations 
        %% Compute Reference (data) and predictions (ACT model) for this iteration
        ReferenceResponsivenessMean   = ones(1,NumberOfPulses)*NaN;
        ReferenceResponsivenessSEM    = ones(1,NumberOfPulses)*NaN;
        PredictedResponsiveness_Mat   = ones(size(ACT.(CurrentCompartment).Responsiveness))*NaN;
        ReferenceLatencyMean          = ones(1,NumberOfPulses)*NaN;
        ReferenceLatencySEM           = ones(1,NumberOfPulses)*NaN;
        ReferenceLatencyNumberOfRepeats = ones(1,NumberOfPulses)*NaN;
        PredictedLatencies_Mat        = ones(size(ACT.(CurrentCompartment).Responsiveness))*NaN;
        
        ReferenceResponsiveness_OneRandomRepeat     = ones(1,NumberOfPulses)*NaN; 
        ReferenceLatency_OneRandomRepeat            = ones(1,NumberOfPulses)*NaN; 
        PredictedResponsiveness_Mat_OneRandomRepeat = ones(length(AdaptationTimeVec),NumberOfPulses)*NaN;
        PredictedLatency_Mat_OneRandomRepeat        = ones(length(AdaptationTimeVec),NumberOfPulses)*NaN;
        
        ReferenceResponsiveness_TwoRandomRepeat     = ones(1,NumberOfPulses)*NaN; 
        ReferenceResponsivenessSEM_TwoRandomRepeat  = ones(1,NumberOfPulses)*NaN; 
        ReferenceLatency_TwoRandomRepeat            = ones(1,NumberOfPulses)*NaN; 
        ReferenceLatencySEM_TwoRandomRepeat         = ones(1,NumberOfPulses)*NaN; 
        PredictedResponsiveness_Mat_TwoRandomRepeat = ones(length(AdaptationTimeVec),NumberOfPulses)*NaN;
        PredictedLatency_Mat_TwoRandomRepeat        = ones(length(AdaptationTimeVec),NumberOfPulses)*NaN;
        ResponsivenessNumberOfRepeats_TwoRandomRepeat = ones(1,NumberOfPulses)*NaN; 
        LatencyNumberOfRepeats_TwoRandomRepeat      = ones(1,NumberOfPulses)*NaN; 
        
        for p_ind = 1:NumberOfPulses
            CurrentNumOfRepeats                  = ReferenceCell.Responsiveness.NumberOfRepeats(p_ind);
            if CurrentNumOfRepeats > 0
                RandNumbers                          = randi(CurrentNumOfRepeats,1,CurrentNumOfRepeats);
                CurrentRepeatIndicesInMatrix         = ReferenceCell.Responsiveness.RepeatsColumnIndicesInMatrix{p_ind}(RandNumbers);
                Vec                                  = ReferenceResponsivenessMatrix(p_ind,CurrentRepeatIndicesInMatrix);
                ReferenceResponsivenessMean(p_ind)   = nanmean(Vec);
                ReferenceResponsivenessSEM(p_ind)    = nanstd(Vec)./sqrt(CurrentNumOfRepeats);
                Matrix                               = squeeze(ACT.(CurrentCompartment).Responsiveness(:,p_ind,CurrentRepeatIndicesInMatrix));
                PredictedResponsiveness_Mat(:,p_ind,1:CurrentNumOfRepeats)= Matrix;

                ReferenceResponsiveness_OneRandomRepeat(p_ind)       = Vec(1);
                PredictedResponsiveness_Mat_OneRandomRepeat(:,p_ind) = Matrix(:,1);
                
                Vec_TwoRandom = Vec(1:min([2 CurrentNumOfRepeats]));                
                ReferenceResponsiveness_TwoRandomRepeat(p_ind)       = nanmean(Vec_TwoRandom); % 1 or 2
                ReferenceResponsivenessSEM_TwoRandomRepeat(p_ind)    = nanstd(Vec_TwoRandom)./sqrt(length(Vec_TwoRandom));
                PredictedResponsiveness_Mat_TwoRandomRepeat(:,p_ind) = nanmean(Matrix(:,1:length(Vec_TwoRandom)),2);
                ResponsivenessNumberOfRepeats_TwoRandomRepeat(p_ind) = length(find(~isnan(Vec_TwoRandom))); % number of repeats that are not NaN
            end
            
            CurrentNumOfRepeats                  = ReferenceCell.Latency.NumberOfRepeats(p_ind);
            if CurrentNumOfRepeats > 0
                Vec                                  = ReferenceLatencyMatrix(p_ind,CurrentRepeatIndicesInMatrix); % Vec also includes NaNs
                ReferenceLatencyMean(p_ind)          = nanmean(Vec);
                N                                    = sum(~isnan(Vec)); % N can be different from 'CurrentNumOfRepeats' since it includes repeats with no responsiveness, i.e. with latency == NaN
                ReferenceLatencyNumberOfRepeats(p_ind) = N;
                ReferenceLatencySEM(p_ind)           = nanstd(Vec)./sqrt(N);
                Matrix                               = squeeze(ACT.(CurrentCompartment).Latency(:,p_ind,CurrentRepeatIndicesInMatrix));
                PredictedLatencies_Mat(:,p_ind,1:length(CurrentRepeatIndicesInMatrix)) = Matrix; % May include NaNs
                
                ReferenceLatency_OneRandomRepeat(p_ind)       = Vec(1);
                PredictedLatency_Mat_OneRandomRepeat(:,p_ind) = Matrix(:,1);                
                                
                Vec_TwoRandom = Vec(1:min([2 CurrentNumOfRepeats]));
                ReferenceLatency_TwoRandomRepeat(p_ind)       = nanmean(Vec_TwoRandom); % 1 or 2
                ReferenceLatencySEM_TwoRandomRepeat(p_ind)    = nanstd(Vec_TwoRandom)./sqrt(length(Vec_TwoRandom));
                PredictedLatency_Mat_TwoRandomRepeat(:,p_ind) = nanmean(Matrix(:,1:length(Vec_TwoRandom)),2);
                LatencyNumberOfRepeats_TwoRandomRepeat(p_ind) = length(find(~isnan(Vec_TwoRandom))); % number of repeats that are not NaN
            end
        end        
        ReferenceLatencyMean(ReferenceLatencyNumberOfRepeats<=1)                    = NaN;  % Do not rely on data that comes from ONE SINGLE data point 
        ReferenceLatency_TwoRandomRepeat(LatencyNumberOfRepeats_TwoRandomRepeat<=1) = NaN;  % Do not rely on data that comes from ONE SINGLE data point 
        PredictedResponsivenessMeanMat = nanmean(PredictedResponsiveness_Mat,3);
        PredictedLatencyMeanMat        = nanmean(PredictedLatencies_Mat,3);
        
        ReferenceResponsivenessMat                  = repmat(ReferenceResponsivenessMean,length(BetaVec),1);
        ReferenceResponsivenessSEM_Mat              = repmat(ReferenceResponsivenessSEM,length(BetaVec),1);
        ResponsivenessNumberOfRepeats_Mat           = repmat(ResponsivenessNumberOfRepeats,length(BetaVec),1);
        ReferenceLatencyMat                         = repmat(ReferenceLatencyMean,length(BetaVec),1);               
        
        %%      FitBasedOnResponsivenessAndLatency_UseDataSTD
        %%%%%%  Compute "cost" = error between prediction and data        
        DistanceBetweenPredictionAndMeasurement_Mat = (PredictedResponsivenessMeanMat - ReferenceResponsivenessMat).^2;              % SE
        IndicesWithNoDistanceAllowed                = (ReferenceResponsivenessSEM_Mat==0)&(ResponsivenessNumberOfRepeats_Mat>1);     % zero standard deviation
        DistanceBetweenPredictionAndMeasurement_Mat(IndicesWithNoDistanceAllowed & (DistanceBetweenPredictionAndMeasurement_Mat>0)) = inf;    
        DistanceBetweenPredictionAndMeasurement     = sqrt(nanmean(DistanceBetweenPredictionAndMeasurement_Mat,2));           % RMSE = RootMeanSquareErrors
        Cost_Responsiveness                         = DistanceBetweenPredictionAndMeasurement; 
        if all(~isfinite(Cost_Responsiveness)) % Skip Setting inf on non-allowed predictions
            DistanceBetweenPredictionAndMeasurement_Mat = (PredictedResponsivenessMeanMat - ReferenceResponsivenessMat).^2;       % SE
            DistanceBetweenPredictionAndMeasurement     = sqrt(nanmean(DistanceBetweenPredictionAndMeasurement_Mat,2));           % RMSE = RootMeanSquareErrors
            Cost_Responsiveness                         = DistanceBetweenPredictionAndMeasurement; 
        end
        
        DistanceBetweenPredictionAndMeasurement_Mat = (PredictedLatencyMeanMat- ReferenceLatencyMat).^2;              % SE
        DistanceBetweenPredictionAndMeasurement     = sqrt(nanmean(DistanceBetweenPredictionAndMeasurement_Mat,2));   % RMSE = RootMeanSquareErrors
        Cost_Latency                                = DistanceBetweenPredictionAndMeasurement; 

    %     Cost_Total = Cost_Responsiveness + Cost_Latency;  
    %     figure; plot(Cost_Responsiveness,'b'); hold on;plot(Cost_Latency,'g');plot(Cost_Total,'r'); %plot(BestFitIndex,Cost_Total(BestFitIndex),'ko')

        % Do not allow "infinite error cost" of either responsiveness or latency   
        Cost_Responsiveness(~isfinite(Cost_Latency))        = NaN;
        Cost_Latency(~isfinite(Cost_Responsiveness))        = NaN;            
        Cost_Responsiveness(~isfinite(Cost_Responsiveness)) = NaN;
        Cost_Latency(~isfinite(Cost_Latency))               = NaN;  

        %%%%%%   Normalize to dimensionless parameters in the range of [0 1]   
        if (nanmax(Cost_Responsiveness)-nanmin(Cost_Responsiveness))>0  % Avoid normalizing by zeros
            Cost_Responsiveness = (Cost_Responsiveness - min(Cost_Responsiveness))./(max(Cost_Responsiveness)-min(Cost_Responsiveness));
        end
        if (nanmax(Cost_Latency)-nanmin(Cost_Latency))>0                % Avoid normalizing by zeros
            Cost_Latency        = (Cost_Latency - min(Cost_Latency))./(max(Cost_Latency)-min(Cost_Latency));    
        end
    %     Cost_Total          = Cost_Responsiveness + Cost_Latency;  
    %     figure; plot(Cost_Responsiveness,'b'); hold on;plot(Cost_Latency,'g');plot(Cost_Total,'r'); %plot(BestFitIndex,Cost_Total(BestFitIndex),'ko')
%         Cost_Responsiveness_all = Cost_Responsiveness;
%         Cost_Latency_all        = Cost_Latency;       
%         Cost_Total_all          = Cost_Responsiveness + Cost_Latency; % For plot

        %%%%%%   Find range of possible good fits. 
        % 1. Between minima for Responsiveness and Latency costs
        % 2. Reponsiveness cost is not out of range

        Indices_Responsiveness = find(Cost_Responsiveness==nanmin(Cost_Responsiveness));
        Indices_Latency        = find(Cost_Latency==nanmin(Cost_Latency));
        Indices_OutOfRange     = false(size(Cost_Responsiveness));                  
            
        if Indices_Responsiveness(1) >= Indices_Latency(end)       % Keep only indices between the two minima
            Indices_OutOfRange(1:Indices_Latency(1)-1)            = true;
            Indices_OutOfRange(Indices_Responsiveness(end)+1:end) = true;
        elseif Indices_Responsiveness(end) <= Indices_Latency(1)
            Indices_OutOfRange(1:Indices_Responsiveness(1)-1)     = true;
            Indices_OutOfRange(Indices_Latency(end)+1:end)        = true;
        else
            % Overlap exist. Anything out of the overlapping area is excluded.
            Indices_OutOfRange  = true(size(Cost_Responsiveness));
            Indices_OutOfRange(intersect(Indices_Responsiveness,Indices_Latency))= false;     
            if all(Indices_OutOfRange)
                % e.g. The wierd situation of two minima points in Cost_Responsiveness that are not connected 
                Indices_OutOfRange(Indices_Responsiveness(1))=false; % keep one value
            end
        end    

        ResponsivenessCostOutOfRange  = Cost_Responsiveness - nanmin(Cost_Responsiveness) >  0.5 ;  % Huge error (on a scale of [0,1])
        Indices_OutOfRange            = Indices_OutOfRange | ResponsivenessCostOutOfRange;

        Cost_Responsiveness(Indices_OutOfRange)= NaN;
        Cost_Latency(Indices_OutOfRange)       = NaN;    
        Cost_Total = Cost_Responsiveness + Cost_Latency; 

        [~,BestFitIndex]  = min(Cost_Total);      
        IndicesAroundMin  = Cost_Total <= 1.1 * min(Cost_Total);                
        LowerLimitIndex   = find(IndicesAroundMin,1,'first'); 
        UpperLimitIndex   = find(IndicesAroundMin,1,'last'); 

%         figure; plot(Cost_Responsiveness_all,'b--'); hold on; plot(Cost_Latency_all,'g--');plot(Cost_Total_all,'r--'); 
%         plot(Cost_Responsiveness,'b','linewidth',2); hold on; plot(Cost_Latency,'g','linewidth',2);plot(Cost_Total,'r','linewidth',2); plot(BestFitIndex,Cost_Total(BestFitIndex),'ko')
%         title(AdaptationTimeVec(BestFitIndex))
        %%%%%%  Assign to structure
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.BestFitIndexInVec(it)            = BestFitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.LowerLimitIndex(it)              = LowerLimitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.UpperLimitIndex(it)              = UpperLimitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.LowerLimitAdaptationTime(it)     = AdaptationTimeVec(LowerLimitIndex);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.UpperLimitAdaptationTime(it)     = AdaptationTimeVec(UpperLimitIndex);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.BestFitAdaptationTime(it)        = AdaptationTimeVec(BestFitIndex);           
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.MeanResponsiveness(it,:)         = PredictedResponsivenessMeanMat(BestFitIndex,:);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.MeanLatency(it,:)                = PredictedLatencyMeanMat(BestFitIndex,:);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.ReferenceResponsiveness(it,:)    = ReferenceResponsivenessMean;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.ReferenceResponsivenessSEM(it,:) = ReferenceResponsivenessSEM;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.ReferenceLatency(it,:)           = ReferenceLatencyMean;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.ReferenceLatencySEM(it,:)        = ReferenceLatencySEM;
        
        %% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        if ~CurrentComputeSpecialFits        
            continue
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        %%   FitBasedOnAllResponsiveness
        DistanceBetweenPredictionAndMeasurement_Mat = (PredictedResponsivenessMeanMat - ReferenceResponsivenessMat).^2;     % SE            
        DistanceBetweenPredictionAndMeasurement     = nansum(DistanceBetweenPredictionAndMeasurement_Mat,2);
        
        LogicalVec               = DistanceBetweenPredictionAndMeasurement == min(DistanceBetweenPredictionAndMeasurement);        
        LowerLimitIndex          = find(LogicalVec,1,'first');
        UpperLimitIndex          = find(LogicalVec,1,'last');
        LowerLimitAdaptationTime = AdaptationTimeVec(LowerLimitIndex);
        UpperLimitAdaptationTime = AdaptationTimeVec(UpperLimitIndex);        
        BestFitAdaptationTime    = (LowerLimitAdaptationTime + UpperLimitAdaptationTime)/2;
        BestFitIndex             = find(AdaptationTimeVec >=  BestFitAdaptationTime, 1, 'first');    
        
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.BestFitIndexInVec(it)            = BestFitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.LowerLimitIndex(it)              = LowerLimitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.UpperLimitIndex(it)              = UpperLimitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.LowerLimitAdaptationTime(it)     = AdaptationTimeVec(LowerLimitIndex);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.UpperLimitAdaptationTime(it)     = AdaptationTimeVec(UpperLimitIndex);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.BestFitAdaptationTime(it)        = AdaptationTimeVec(BestFitIndex);           
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.MeanResponsiveness(it,:)         = PredictedResponsivenessMeanMat(BestFitIndex,:);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.MeanLatency(it,:)                = PredictedLatencyMeanMat(BestFitIndex,:);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.ReferenceResponsiveness(it,:)    = ReferenceResponsivenessMean;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.ReferenceResponsivenessSEM(it,:) = ReferenceResponsivenessSEM;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.ReferenceLatency(it,:)           = ReferenceLatencyMean;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnAllResponsiveness.ReferenceLatencySEM(it,:)        = ReferenceLatencySEM;
                
        %%  FitBasedOnOnePulseResponsiveness             
        ReferenceResponsivenessVec                                   = ReferenceResponsivenessMean;    
        ReferenceResponsivenessVec(ResponsivenessNumberOfRepeats<=1) = NaN;    % Do not rely on data that comes from ONE SINGLE data point     
        ReferencePulseIndex                     = find(ReferenceResponsivenessVec > 0, 1, 'last'); % last responsive condition
        if ReferenceResponsivenessVec(ReferencePulseIndex)==100 
            % 100% responsiveness and next pulse is 0% responsiveness. In this case TWO pulses needs to be used (include also ReferencePulseIndex+1)   
            ReferencePulseIndex2 = ReferencePulseIndex + 1;
            DistanceBetweenPredictionAndMeasurement = (PredictedResponsivenessMeanMat(:,ReferencePulseIndex)  - ReferenceResponsivenessVec(ReferencePulseIndex)).^2 + ...
                                                      (PredictedResponsivenessMeanMat(:,ReferencePulseIndex2) - ReferenceResponsivenessVec(ReferencePulseIndex2)).^2;
        else
            % At this pulse: 0 < responsiveness < 100. Sufficient for adaptation time calculation
            ReferencePulseIndex2 = NaN;
            DistanceBetweenPredictionAndMeasurement = (PredictedResponsivenessMeanMat(:,ReferencePulseIndex) - ReferenceResponsivenessVec(ReferencePulseIndex)).^2;
        end

        LogicalVec               = DistanceBetweenPredictionAndMeasurement == min(DistanceBetweenPredictionAndMeasurement);        
        LowerLimitIndex          = find(LogicalVec,1,'first');
        UpperLimitIndex          = find(LogicalVec,1,'last');
        LowerLimitAdaptationTime = AdaptationTimeVec(LowerLimitIndex);
        UpperLimitAdaptationTime = AdaptationTimeVec(UpperLimitIndex);      
        if UpperLimitAdaptationTime > 10*LowerLimitAdaptationTime % Out of bound- numerical error
            UpperLimitIndex          = LowerLimitIndex;
            UpperLimitAdaptationTime = LowerLimitAdaptationTime;
        end
        BestFitAdaptationTime    = (LowerLimitAdaptationTime + UpperLimitAdaptationTime)/2;
        BestFitIndex             = find(AdaptationTimeVec >=  BestFitAdaptationTime, 1, 'first');   

        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.ReferencePulseIndex(it)          = ReferencePulseIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.ReferencePulseIndex2(it)         = ReferencePulseIndex2;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.ReferencePulseName{it}           = PulsesNames{ReferencePulseIndex};
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.BestFitIndexInVec(it)            = BestFitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.LowerLimitIndex(it)              = LowerLimitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.UpperLimitIndex(it)              = UpperLimitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.LowerLimitAdaptationTime(it)     = AdaptationTimeVec(LowerLimitIndex);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.UpperLimitAdaptationTime(it)     = AdaptationTimeVec(UpperLimitIndex);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.BestFitAdaptationTime(it)        = AdaptationTimeVec(BestFitIndex);           
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.MeanResponsiveness(it,:)         = PredictedResponsivenessMeanMat(BestFitIndex,:);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.MeanLatency(it,:)                = PredictedLatencyMeanMat(BestFitIndex,:);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.ReferenceResponsiveness(it,:)    = ReferenceResponsivenessVec;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.ReferenceResponsivenessSEM(it,:) = ReferenceResponsivenessSEM;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.ReferenceLatency(it,:)           = ReferenceLatencyMean;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.ReferenceLatencySEM(it,:)        = ReferenceLatencySEM;
                
        %%  FitBasedOnOnePulseResponsivenessAndLatency             
        ReferenceResponsivenessVec                                   = ReferenceResponsivenessMean;    
        ReferenceResponsivenessVec(ResponsivenessNumberOfRepeats<=1) = NaN;    % Do not rely on data that comes from ONE SINGLE data point     
        ReferenceLatencyVec                                          = ReferenceLatencyMean;    
        
        LogicalVector                             = ReferenceResponsivenessVec > 0;
        LogicalVector(isnan(ReferenceLatencyVec)) = false;                
        ReferencePulseIndex                       = find(LogicalVector, 1, 'last'); % last responsive condition with Latency data
        DistanceBetweenPredictionAndMeasurement   = (PredictedResponsivenessMeanMat(:,ReferencePulseIndex) - ReferenceResponsivenessVec(ReferencePulseIndex)).^2;
        % Allow Flexibility Around BestResponsiveness
        LogicalVec_Responsiveness1                  = DistanceBetweenPredictionAndMeasurement <= prctile(DistanceBetweenPredictionAndMeasurement,1);
        LogicalVec_Responsiveness2                  = DistanceBetweenPredictionAndMeasurement <= 1.1*min(DistanceBetweenPredictionAndMeasurement);
        LogicalVec_Responsiveness                   = LogicalVec_Responsiveness1 | LogicalVec_Responsiveness2;
        
        DistanceBetweenPredictionAndMeasurement   = (PredictedLatencyMeanMat(:,ReferencePulseIndex) - ReferenceLatencyVec(ReferencePulseIndex)).^2;          
        DistanceBetweenPredictionAndMeasurement(~LogicalVec_Responsiveness) = inf;               % First condition:  must best fit reponsiveness
        
        LogicalVec               = DistanceBetweenPredictionAndMeasurement == min(DistanceBetweenPredictionAndMeasurement);        
        LowerLimitIndex          = find(LogicalVec,1,'first');
        UpperLimitIndex          = find(LogicalVec,1,'last');
        LowerLimitAdaptationTime = AdaptationTimeVec(LowerLimitIndex);
        UpperLimitAdaptationTime = AdaptationTimeVec(UpperLimitIndex);        
        BestFitAdaptationTime    = (LowerLimitAdaptationTime + UpperLimitAdaptationTime)/2;
        BestFitIndex             = find(AdaptationTimeVec >=  BestFitAdaptationTime, 1, 'first');   

        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.ReferencePulseIndex(it)          = ReferencePulseIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.ReferencePulseName{it}           = PulsesNames{ReferencePulseIndex};
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.BestFitIndexInVec(it)            = BestFitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.LowerLimitIndex(it)              = LowerLimitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.UpperLimitIndex(it)              = UpperLimitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.LowerLimitAdaptationTime(it)     = AdaptationTimeVec(LowerLimitIndex);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.UpperLimitAdaptationTime(it)     = AdaptationTimeVec(UpperLimitIndex);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.BestFitAdaptationTime(it)        = AdaptationTimeVec(BestFitIndex);           
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.MeanResponsiveness(it,:)         = PredictedResponsivenessMeanMat(BestFitIndex,:);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.MeanLatency(it,:)                = PredictedLatencyMeanMat(BestFitIndex,:);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.ReferenceResponsiveness(it,:)    = ReferenceResponsivenessVec;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.ReferenceResponsivenessSEM(it,:) = ReferenceResponsivenessSEM;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.ReferenceLatency(it,:)           = ReferenceLatencyVec;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsivenessAndLatency.ReferenceLatencySEM(it,:)        = ReferenceLatencySEM;                       
     
        %%  FitBasedOnTwoRepeatPerCondition  (Based on "Simple" Responsiveness and Latency analysis, no flexibility for Responsivenss)         
        ReferenceResponsivenessMat_TwoRandomRepeat = repmat(ReferenceResponsiveness_TwoRandomRepeat,length(BetaVec),1);
        ReferenceLatencyMat_TwoRandomRepeat        = repmat(ReferenceLatency_TwoRandomRepeat,length(BetaVec),1); 

        DistanceBetweenPredictionAndMeasurement_Mat = (PredictedResponsiveness_Mat_TwoRandomRepeat - ReferenceResponsivenessMat_TwoRandomRepeat).^2;              % SE
        DistanceBetweenPredictionAndMeasurement     = sqrt(nanmean(DistanceBetweenPredictionAndMeasurement_Mat,2));           % RMSE = RootMeanSquareErrors        
        LogicalVec_Responsiveness                   = DistanceBetweenPredictionAndMeasurement == min(DistanceBetweenPredictionAndMeasurement);

        DistanceBetweenPredictionAndMeasurement_Mat = (PredictedLatency_Mat_TwoRandomRepeat- ReferenceLatencyMat_TwoRandomRepeat).^2;               % SE
        DistanceBetweenPredictionAndMeasurement     = nanmean(DistanceBetweenPredictionAndMeasurement_Mat,2);
        DistanceBetweenPredictionAndMeasurement(~LogicalVec_Responsiveness) = inf;                   % First condition:  must best fit reponsiveness
        [~,BestFitIndex]                            = min(DistanceBetweenPredictionAndMeasurement);  % Second condition: best fit latency
                
        %%%%%%  Assign to structure
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnTwoRepeatPerCondition.BestFitIndexInVec(it)            = BestFitIndex;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnTwoRepeatPerCondition.BestFitAdaptationTime(it)        = AdaptationTimeVec(BestFitIndex);           
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnTwoRepeatPerCondition.MeanResponsiveness(it,:)         = PredictedResponsiveness_Mat_TwoRandomRepeat(BestFitIndex,:);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnTwoRepeatPerCondition.MeanLatency(it,:)                = PredictedLatency_Mat_TwoRandomRepeat(BestFitIndex,:);
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnTwoRepeatPerCondition.ReferenceResponsiveness(it,:)    = ReferenceResponsiveness_TwoRandomRepeat;
        ACT.(CurrentCompartment).Bootstrap.FitBasedOnTwoRepeatPerCondition.ReferenceLatency(it,:)           = ReferenceLatency_TwoRandomRepeat;                        
               
    end    
end
disp([datestr(now),'  --  Adaptation time calculation is complete']);

for c_ind = 1:length(CompartmentFieldNames)       
    CurrentCompartment = CompartmentFieldNames{c_ind}; 
    NumberOfRepeats.Total.(CurrentCompartment)           = sum(MeasuredStats.Soma.Responsiveness.NumberOfRepeats);
    NumberOfRepeats.TwoPerCondition.(CurrentCompartment) = NumberOfPulses * 2;
    NumberOfRepeats.Percentage.(CurrentCompartment)      = 100* NumberOfRepeats.TwoPerCondition.(CurrentCompartment) / NumberOfRepeats.Total.(CurrentCompartment); 
end

%% Reduce structure size (~1 order of magnitude difference)
for c_ind = 1:length(CompartmentFieldNames)        
    CurrentCompartment = CompartmentFieldNames{c_ind};
    ACT.(CurrentCompartment) = rmfield(ACT.(CurrentCompartment),'ACT');     
end

%% Assign values into ModelParams structure 
CurrentModelParams_Output = CurrentModelParams;
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};
    if ComputeSpecialFitsForAllCompartments || c_ind==1            
        CurrentComputeSpecialFits = ComputeSpecialFits;
    else
        CurrentComputeSpecialFits = false;
    end    
    
    CurrentModelParams_Output.(CurrentCompartment).ExperimentForAdaptationTime = CurrentStrain_Slow;
    % Best Fit using all data
    CurrentModelParams_Output.(CurrentCompartment).Tau_BestFitAllData = ...
                          ACT.(CurrentCompartment).FitBasedOnResponsivenessAndLatency_UseDataSTD.BestFitAdaptationTime;
    Vector =              ACT.(CurrentCompartment).Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.BestFitAdaptationTime;
    CurrentModelParams_Output.(CurrentCompartment).Bootstrap.Tau_BestFitAllData.Vector = Vector;
    CurrentModelParams_Output.(CurrentCompartment).Bootstrap.Tau_BestFitAllData.Mean   = nanmean(Vector);
    CurrentModelParams_Output.(CurrentCompartment).Bootstrap.Tau_BestFitAllData.Median = nanmedian(Vector);
    CurrentModelParams_Output.(CurrentCompartment).Bootstrap.Tau_BestFitAllData.STD    = nanstd(Vector);
                      
    if ~CurrentComputeSpecialFits        
        continue
    end
    CurrentModelParams_Output.(CurrentCompartment).Tau_BestFitOnePulseResponsiveness = ...
                          ACT.(CurrentCompartment).FitBasedOnOnePulseResponsiveness.BestFitAdaptationTime;
    Vector =              ACT.(CurrentCompartment).Bootstrap.FitBasedOnOnePulseResponsiveness.BestFitAdaptationTime;
    CurrentModelParams_Output.(CurrentCompartment).Bootstrap.Tau_BestFitOnePulseResponsiveness.Vector = Vector;
    CurrentModelParams_Output.(CurrentCompartment).Bootstrap.Tau_BestFitOnePulseResponsiveness.Mean   = nanmean(Vector);
    CurrentModelParams_Output.(CurrentCompartment).Bootstrap.Tau_BestFitOnePulseResponsiveness.Median = nanmedian(Vector);
    CurrentModelParams_Output.(CurrentCompartment).Bootstrap.Tau_BestFitOnePulseResponsiveness.STD    = nanstd(Vector);
end

return

function [ACT, MeasuredStats] = CalculateACTpredictionsWithGivenModelParams (Data, CurrentStrain_Slow, CurrentModelParams, ReferenceACTStructure, ReferenceACTfitFields)
% K and Beta are parameters used for the adaptive threshold model.   
% Here we find predictions for sensory responses using on K and Beta that were calculated for a different set of experiments  

%% Example for Input arguments examples:
% CurrentStrain_Slow    = 'str2_CX17256_ToBuffer_Slow';
% CurrentModelParams    = ModelParams.Naive;
% ReferenceACTStructure = ACT_CX17256; % computed for 'str2_CX17256_FromD6_Slow'
% ReferenceACTfitFields = {'FitBasedOnResponsivenessAndLatency_UseDataSTD','FitBasedOnTwoRepeatPerCondition',...
%                          'FitBasedOnAllResponsiveness','FitBasedOnOnePulseResponsiveness'};
% [ACT, MeasuredStats] = CalculateACTpredictionsWithGivenModelParams (Data, CurrentStrain_Slow, CurrentModelParams, ReferenceACTStructure, ReferenceACTfitFields);  

%% Initialize structures
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
MaxNumberOfRepeats    = 40;  % free parameter. It just needs to be higher than the maximal number of repeats

DyeInformation      = Data.(CurrentStrain_Slow).DyeInformation;
PulsesNames         = Data.(CurrentStrain_Slow).PulsesInfo.PulseFromButanone;
NumberOfPulses      = length(PulsesNames);            

AdaptationTimeVec = [0.01 0.05 0.1:0.05:0.95 1:0.02:2 2.2:0.2:12.8 13:0.1:19 19.2:0.2:21 22:1:30 32:2:40 45:5:70 80:10:120 150:50:1e4 1e5]; % For screen
BetaVec           = 1./AdaptationTimeVec;
C                 = DyeInformation.C_Smoothed;    
dCdt              = DyeInformation.dCdt_Smoothed;
DeltaCoverC       = DyeInformation.DeltaCoverC;
d2Cdt2            = DyeInformation.d2Cdt2;
NaN_Mat           = all(isnan(C),3); % true is no vector is empty (only NaNs)
TimeVec           = DyeInformation.TimeVectorInSec;
DyeVectorLength   = length(TimeVec);

% Find Measured stats for each compartment and condition (average over repeats)  
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};   
    MeasuredStats.(CurrentCompartment).Responsiveness = Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Responsiveness']);
    MeasuredStats.(CurrentCompartment).Latency        = Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Latency_InSec']);
    MeasuredStats.(CurrentCompartment).Latency_MA     = Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Latency_MA_InSec']);
    MeasuredStats.(CurrentCompartment).C_AtNeuronActivation.Matrix           = DyeInformation.Stats.C_AtNeuronActivation.(CurrentCompartment);
    MeasuredStats.(CurrentCompartment).dCdt_AtNeuronActivation.Matrix        = DyeInformation.Stats.dCdt_AtNeuronActivation.(CurrentCompartment);
    MeasuredStats.(CurrentCompartment).d2Cdt2_AtNeuronActivation.Matrix      = DyeInformation.Stats.d2Cdt2_AtNeuronActivation.(CurrentCompartment);
    MeasuredStats.(CurrentCompartment).DeltaCoverC_AtNeuronActivation.Matrix = DyeInformation.Stats.DeltaCoverC_AtNeuronActivation.(CurrentCompartment);
end
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};   
    MeasuredStats.(CurrentCompartment).Responsiveness.NumberOfRepeats = sum(~isnan(Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Responsiveness']).Matrix),2)';  
    MeasuredStats.(CurrentCompartment).Latency.NumberOfRepeats        = sum(~isnan(Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Latency_InSec']).Matrix),2)';  
    MeasuredStats.(CurrentCompartment).Latency_MA.NumberOfRepeats     = sum(~isnan(Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Latency_MA_InSec']).Matrix),2)'; 
end

% Initialization
ACT.Information.CompartmentFieldNames = CompartmentFieldNames;
ACT.Information.StrainName            = CurrentStrain_Slow;
ACT.Information.ModelParams           = CurrentModelParams;
ACT.Information.MeasuredStats         = MeasuredStats;
ACT.Information.InitialConcentrationPerPulse = DyeInformation.InitialConcentrationPerPulse;
ACT.Information.FinalConcentrationPerPulse   = DyeInformation.FinalConcentrationPerPulse;
FieldsInFitToCopyFromReference = {'BestFitIndexInVec','LowerLimitIndex','UpperLimitIndex',...
                                  'BestFitAdaptationTime','LowerLimitAdaptationTime','UpperLimitAdaptationTime'};
for c_ind = 1:length(CompartmentFieldNames)        
    CurrentCompartment = CompartmentFieldNames{c_ind};
    ACT.(CurrentCompartment).K  = CurrentModelParams.(CurrentCompartment).K;        

    for fit_ind = 1:length(ReferenceACTfitFields)
        CurrentFit = ReferenceACTfitFields{fit_ind};
        if isfield(ReferenceACTStructure.(CurrentCompartment),CurrentFit)            
            for field_ind = 1:length(FieldsInFitToCopyFromReference)
                CurrentField = FieldsInFitToCopyFromReference{field_ind};            
                ACT.(CurrentCompartment).(CurrentFit).(CurrentField)          = ReferenceACTStructure.(CurrentCompartment).(CurrentFit).(CurrentField);
                ACT.(CurrentCompartment).Bootstrap.(CurrentFit).(CurrentField)= ReferenceACTStructure.(CurrentCompartment).Bootstrap.(CurrentFit).(CurrentField);
            end        
        end
    end
    
    ACT.(CurrentCompartment).AdaptationTimeVec                  = AdaptationTimeVec;     
    ACT.(CurrentCompartment).StatsPerPulse.FinalConcentration   = DyeInformation.FinalConcentrationPerPulse;     
    ACT.(CurrentCompartment).StatsPerPulse.InitialConcentration = DyeInformation.InitialConcentrationPerPulse;    
    ACT.(CurrentCompartment).StatsPerPulse.NumberOfRepeats      = Data.(CurrentStrain_Slow).Stats.Soma_Information.NumberOfRepeats;
    ACT.(CurrentCompartment).ACT                       = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,DyeVectorLength,'single')*NaN; 
    ACT.(CurrentCompartment).Responsiveness            = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
    ACT.(CurrentCompartment).Latency                   = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
    ACT.(CurrentCompartment).C_at_Activation           = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
    ACT.(CurrentCompartment).dCdt_at_Activation        = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
    ACT.(CurrentCompartment).d2Cdt2_at_Activation      = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
    ACT.(CurrentCompartment).DeltaCoverC_at_Activation = zeros(length(BetaVec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN; 
end

%% Find adaptive threshold for each compartment and condition
for p_ind = 1:NumberOfPulses         
%     disp([datestr(now), '  --  ',PulsesNames{p_ind}])
    for rep_ind = 1:MaxNumberOfRepeats            
        if NaN_Mat(p_ind,rep_ind)
            % no more traces for this pulse
            break
        end
        Current_C           = double(squeeze(C(p_ind,rep_ind,:)));
        Current_dCdt        = double(squeeze(dCdt(p_ind,rep_ind,:)));
        Current_d2Cdt2      = double(squeeze(d2Cdt2(p_ind,rep_ind,:)));
        Current_DeltaCoverC = double(squeeze(DeltaCoverC(p_ind,rep_ind,:)));
        
        for c_ind = 1:length(CompartmentFieldNames)
            CurrentCompartment = CompartmentFieldNames{c_ind};     
            K = ACT.(CurrentCompartment).K;
        
            for B_ind = 1:length(BetaVec)                
                Beta = BetaVec(B_ind);

                %% Loop over Beta + Compare to actual initiation time
                [AdaptiveThresholdVector, PredictedActivationIndex, PredictedActivationTime] = CalculateAdaptiveThresholdFromDyePattern_ExpSaturation_03 (TimeVec, Current_C, K, Beta);
                ACT.(CurrentCompartment).ACT(B_ind,p_ind,rep_ind,1:DyeVectorLength)          = AdaptiveThresholdVector;

                if ~isempty(PredictedActivationTime)                    
                    ACT.(CurrentCompartment).Responsiveness(B_ind,p_ind,rep_ind)            = 100;       
                    ACT.(CurrentCompartment).Latency(B_ind,p_ind,rep_ind)                   = PredictedActivationTime;    
                    ACT.(CurrentCompartment).C_at_Activation(B_ind,p_ind,rep_ind)           = Current_C(PredictedActivationIndex)   ;    
                    ACT.(CurrentCompartment).dCdt_at_Activation(B_ind,p_ind,rep_ind)        = Current_dCdt(PredictedActivationIndex);    
                    ACT.(CurrentCompartment).d2Cdt2_at_Activation(B_ind,p_ind,rep_ind)      = Current_d2Cdt2(PredictedActivationIndex);    
                    ACT.(CurrentCompartment).DeltaCoverC_at_Activation(B_ind,p_ind,rep_ind) = Current_DeltaCoverC(PredictedActivationIndex);    
                else            
                    ACT.(CurrentCompartment).Responsiveness(B_ind,p_ind,rep_ind)            = 0;   
                end
            end
        end                
    end
end
disp([datestr(now), '  --  ','Finished ACT screening'])           
     
% Find ACT-predicted stats for each compartment and condition (average over repeats)  
FieldsForStats = {'Responsiveness','Latency','C_at_Activation','dCdt_at_Activation','d2Cdt2_at_Activation','DeltaCoverC_at_Activation'};
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};   
    for f_ind = 1:length(FieldsForStats)
        CurrentField = FieldsForStats{f_ind};
        ACT.(CurrentCompartment).StatsPerPulse.(CurrentField).Mean   = nanmean(ACT.(CurrentCompartment).(CurrentField),3);
        ACT.(CurrentCompartment).StatsPerPulse.(CurrentField).Median = nanmedian(ACT.(CurrentCompartment).(CurrentField),3);
        ACT.(CurrentCompartment).StatsPerPulse.(CurrentField).STD    = nanstd(ACT.(CurrentCompartment).(CurrentField),[],3);
        ACT.(CurrentCompartment).StatsPerPulse.(CurrentField).SEM    = ACT.(CurrentCompartment).StatsPerPulse.(CurrentField).STD ./ ...
                                                           sqrt(repmat(ACT.(CurrentCompartment).StatsPerPulse.NumberOfRepeats',length(BetaVec),1));
    end
end

for c_ind = 1:length(CompartmentFieldNames)       
    CurrentCompartment = CompartmentFieldNames{c_ind}; 
    NumberOfRepeats.Total.(CurrentCompartment)           = sum(MeasuredStats.Soma.Responsiveness.NumberOfRepeats);
    NumberOfRepeats.TwoPerCondition.(CurrentCompartment) = NumberOfPulses * 2;
    NumberOfRepeats.Percentage.(CurrentCompartment)      = 100* NumberOfRepeats.TwoPerCondition.(CurrentCompartment) / NumberOfRepeats.Total.(CurrentCompartment); 
end

%% Reduce structure size (~1 order of magnitude difference)
for c_ind = 1:length(CompartmentFieldNames)        
    CurrentCompartment = CompartmentFieldNames{c_ind};
    ACT.(CurrentCompartment) = rmfield(ACT.(CurrentCompartment),'ACT');     
end

return

% Calculate prediction of alternative (non ACT) models
function [ModelPredictions, MeasuredStats] = ComputeAlternativeModelParams (Data, CurrentStrain_Slow, CurrentStrain_Fast, PreviousMeasuredBestFit)

%% Input arguments examples:
% CurrentStrain_Slow  = 'str2_CX17256_FromD6_Slow';
% CurrentStrain_Fast  = 'str2_CX17256_FromD6_Fast';
% CurrentModelParams  = ModelParams.Naive;

% Initialization
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
MaxNumberOfRepeats    = 40;  % free parameter. It just needs to be higher than the maximal number of repeats
DyeInformation      = Data.(CurrentStrain_Slow).DyeInformation;
PulsesNames         = Data.(CurrentStrain_Slow).PulsesInfo.PulseFromButanone;
NumberOfPulses      = length(PulsesNames);        
C           = DyeInformation.C_Smoothed;    
dCdt        = DyeInformation.dCdt_Smoothed;
DeltaCoverC = DyeInformation.DeltaCoverC;
d2Cdt2      = DyeInformation.d2Cdt2;
NaN_Mat     = all(isnan(C),3); % true is no vector is empty (only NaNs)
TimeVec     = DyeInformation.TimeVectorInSec;

% Find Measured stats for each compartment 
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};   
    MeasuredStats.(CurrentCompartment).Threshold_FromStepResponses = Data.(CurrentStrain_Fast).Stats.([CurrentCompartment,'_Threshold']);
    MeasuredStats.(CurrentCompartment).Responsiveness = Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Responsiveness']);
    MeasuredStats.(CurrentCompartment).Latency        = Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Latency_InSec']);
    MeasuredStats.(CurrentCompartment).Latency_MA     = Data.(CurrentStrain_Slow).Stats.([CurrentCompartment,'_Latency_MA_InSec']);
    MeasuredStats.(CurrentCompartment).C_AtNeuronActivation.Matrix           = DyeInformation.Stats.C_AtNeuronActivation.(CurrentCompartment);
    MeasuredStats.(CurrentCompartment).dCdt_AtNeuronActivation.Matrix        = DyeInformation.Stats.dCdt_AtNeuronActivation.(CurrentCompartment);
    MeasuredStats.(CurrentCompartment).d2Cdt2_AtNeuronActivation.Matrix      = DyeInformation.Stats.d2Cdt2_AtNeuronActivation.(CurrentCompartment);
    MeasuredStats.(CurrentCompartment).DeltaCoverC_AtNeuronActivation.Matrix = DyeInformation.Stats.DeltaCoverC_AtNeuronActivation.(CurrentCompartment);
    MeasuredStats.(CurrentCompartment).dCdt_TimeToPeak.Matrix                = DyeInformation.Stats.dCdt_TimeToPeak;
    MeasuredStats.(CurrentCompartment).d2Cdt2_TimeToPeak.Matrix              = DyeInformation.Stats.d2Cdt2_TimeToPeak;
    MeasuredStats.(CurrentCompartment).DeltaCoverC_TimeToPeak.Matrix         = DyeInformation.Stats.DeltaCoverC_TimeToPeak;
    
    BestFit.(CurrentCompartment).C.Mean_AtNeuronActivation           = nanmean(MeasuredStats.(CurrentCompartment).C_AtNeuronActivation.Matrix(:));
    BestFit.(CurrentCompartment).dCdt.Mean_AtNeuronActivation        = nanmean(MeasuredStats.(CurrentCompartment).dCdt_AtNeuronActivation.Matrix(:));
    BestFit.(CurrentCompartment).d2Cdt2.Mean_AtNeuronActivation      = nanmean(MeasuredStats.(CurrentCompartment).d2Cdt2_AtNeuronActivation.Matrix(:));
    BestFit.(CurrentCompartment).DeltaCoverC.Mean_AtNeuronActivation = nanmean(MeasuredStats.(CurrentCompartment).DeltaCoverC_AtNeuronActivation.Matrix(:));
    BestFit.(CurrentCompartment).C_FitStepResponseThreshold.Mean_AtNeuronActivation  = MeasuredStats.(CurrentCompartment).Threshold_FromStepResponses;
end
if exist('PreviousMeasuredBestFit','var')
%    disp('Forcing BestFit based on previous measurements')
   BestFit = PreviousMeasuredBestFit;    
end

% Initialization
ModelPredictions.Information.CompartmentFieldNames = CompartmentFieldNames;
ModelPredictions.Information.StrainName            = CurrentStrain_Slow;
ModelPredictions.Information.StrainName_Fast       = CurrentStrain_Fast;
ModelPredictions.Information.MeasuredStats         = MeasuredStats;
ModelPredictions.Information.BestFit               = BestFit;
ModelPredictions.Information.InitialConcentrationPerPulse = DyeInformation.InitialConcentrationPerPulse;
ModelPredictions.Information.FinalConcentrationPerPulse   = DyeInformation.FinalConcentrationPerPulse;
ModelPredictions.Information.PulsesInfo                   = Data.(CurrentStrain_Slow).PulsesInfo;

% Find model predictions using BestFit parameters
Model_Fields   = {'C',         'dCdt',         'd2Cdt2','DeltaCoverC','C_FitStepResponseThreshold'};
MatchingFields = {'C_Smoothed','dCdt_Smoothed','d2Cdt2','DeltaCoverC','C_Smoothed'};
TimeZeroIndex  = find(TimeVec==0);
StructureInitialization.IndexInTimeVec  = ones(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;
StructureInitialization.Responsiveness  = ones(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;
StructureInitialization.Latency         = ones(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;
StructureInitialization.C               = ones(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;
StructureInitialization.dCdt            = ones(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;
StructureInitialization.d2Cdt2          = ones(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;
StructureInitialization.DeltaCoverC     = ones(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;

for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};  
    for m_ind = 1:length(Model_Fields)
        CurrentModel            = Model_Fields{m_ind};
        CurrentMatchingField    = MatchingFields{m_ind};
        Mean_AtNeuronActivation = BestFit.(CurrentCompartment).(CurrentModel).Mean_AtNeuronActivation;
        ModelPredictions.(CurrentCompartment).(CurrentModel) = StructureInitialization;
        for p_ind = 1:NumberOfPulses
            for rep_ind = 1:MaxNumberOfRepeats
                if NaN_Mat(p_ind,rep_ind)
                    break
                end
                CurrentVector = squeeze(DyeInformation.(CurrentMatchingField)(p_ind,rep_ind,:));
                index         = find(CurrentVector <= Mean_AtNeuronActivation, 1,'first');
                if isempty(index)
                    ModelPredictions.(CurrentCompartment).(CurrentModel).Responsiveness(p_ind,rep_ind) = 0;
                else
                    ModelPredictions.(CurrentCompartment).(CurrentModel).Responsiveness(p_ind,rep_ind) = 100;
                    if index<TimeZeroIndex
                        index = TimeZeroIndex;
                    end
                    ModelPredictions.(CurrentCompartment).(CurrentModel).IndexInTimeVec(p_ind,rep_ind) = index;
                    ModelPredictions.(CurrentCompartment).(CurrentModel).Latency(p_ind,rep_ind)        = TimeVec(index);
                    ModelPredictions.(CurrentCompartment).(CurrentModel).C(p_ind,rep_ind)              = C(p_ind,rep_ind,index);
                    ModelPredictions.(CurrentCompartment).(CurrentModel).dCdt(p_ind,rep_ind)           = dCdt(p_ind,rep_ind,index);
                    ModelPredictions.(CurrentCompartment).(CurrentModel).d2Cdt2(p_ind,rep_ind)         = d2Cdt2(p_ind,rep_ind,index);
                    ModelPredictions.(CurrentCompartment).(CurrentModel).DeltaCoverC(p_ind,rep_ind)    = DeltaCoverC(p_ind,rep_ind,index);
                end                                       
            end
        end                               
    end    
end

%% Add L-N model to model predictions
if isfield(Data.(CurrentStrain_Slow).(PulsesNames{1}),'Soma_LNmodelPredictedDeltaFoverF')
    % Parameters used for Measured Traces analysis of responsiveness and latency    
    DeflectionPoint_LowDiffThreshold                = 5e-3;    % until you reach this lower threshold
    DynamicRangeThresholdForResponsiveness.Soma     = 0.2;  
    DynamicRangeThresholdForResponsiveness.Dendrite = 0.2;  
    DynamicRangeThresholdForResponsiveness.Axon     = 0.2;     % axon responses tend to be more moderate >> requires lower threhold  
    MinimalChangeRequiredForActivityInitiation      = 0.05;     % 5% deltaF/F        
    
    CurrentModel = 'LNmodel';
    for c_ind = 1:length(CompartmentFieldNames(1:3))
        CurrentCompartment  = CompartmentFieldNames{c_ind}; 
        CurrentLNmodelField = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF'];
        CurrentDynamicRangeThreshold = DynamicRangeThresholdForResponsiveness.(CurrentCompartment);
        
        ModelPredictions.(CurrentCompartment).(CurrentModel) = StructureInitialization;
        
        for p_ind = 1:NumberOfPulses
            CurrentPulse         = PulsesNames{p_ind};
            PredictedActivityMat = Data.(CurrentStrain_Slow).(CurrentPulse).(CurrentLNmodelField);
            DynamicRangeVec      = nanmax(PredictedActivityMat,[],2);
            ResponsivenessVec    = 100* (DynamicRangeVec > CurrentDynamicRangeThreshold);  
            TimeVec              = Data.(CurrentStrain_Slow).(CurrentPulse).Dye_TimeVec;
            TimeZeroIndex        = find(TimeVec==0,1);            
            
            LatencyVec           = ones(1,MaxNumberOfRepeats)*NaN;
            IndexForLatencyVec   = ones(1,MaxNumberOfRepeats)*NaN;
            for rep_ind = 1:MaxNumberOfRepeats
                if NaN_Mat(p_ind,rep_ind)
                    break
                end
                if ResponsivenessVec(rep_ind)==0                     
                    ModelPredictions.(CurrentCompartment).(CurrentModel).Responsiveness(p_ind,rep_ind) = 0;                    
                else
                    ModelPredictions.(CurrentCompartment).(CurrentModel).Responsiveness(p_ind,rep_ind) = 100;
                    CurrentVec     = PredictedActivityMat(rep_ind,:);
                    CurrentDiffVec = [NaN diff(CurrentVec)];
                                        
                    %% Find Latency by deflection point and by threshold value
                    % Deflection point by differential
                    IndexAfterDeflectionPoint = find(CurrentVec>CurrentDynamicRangeThreshold,1,'first');
                    IndexOfDeflectionPoint    = find(CurrentDiffVec(1:IndexAfterDeflectionPoint)<=DeflectionPoint_LowDiffThreshold,1,'last'); 
                    if TimeZeroIndex > IndexOfDeflectionPoint  
                       IndexOfDeflectionPoint = TimeZeroIndex;
                    end 

                    % Deflection point by max value with at least minimal differential            
                    IndexAfterDeflectionPoint = find(CurrentVec>CurrentDynamicRangeThreshold,1,'first');
                    IndexOfHighActivity       = find(CurrentVec(1:IndexAfterDeflectionPoint)<=MinimalChangeRequiredForActivityInitiation,1,'last')+1;
                    if IndexOfHighActivity< TimeZeroIndex
                        % assume error, assign maximum length
                         IndexOfHighActivity = length(CurrentVec);
                    end
  
                    %%%% Deflection point, by activity or by diff(activity)?
                    if isempty(IndexOfDeflectionPoint) || (IndexOfHighActivity < IndexOfDeflectionPoint) 
                        IndexForLatency = IndexOfHighActivity;
                    else
                        IndexForLatency = IndexOfDeflectionPoint;
                    end
                    IndexForLatencyVec(rep_ind) = IndexForLatency;    
                    LatencyVec(rep_ind)         = TimeVec(IndexForLatency);                                                       
                                        
                    ModelPredictions.(CurrentCompartment).(CurrentModel).IndexInTimeVec(p_ind,rep_ind) = IndexForLatency;
                    ModelPredictions.(CurrentCompartment).(CurrentModel).Latency(p_ind,rep_ind)        = TimeVec(IndexForLatency);
                    ModelPredictions.(CurrentCompartment).(CurrentModel).C(p_ind,rep_ind)              = C(p_ind,rep_ind,IndexForLatency);
                    ModelPredictions.(CurrentCompartment).(CurrentModel).dCdt(p_ind,rep_ind)           = dCdt(p_ind,rep_ind,IndexForLatency);
                    ModelPredictions.(CurrentCompartment).(CurrentModel).d2Cdt2(p_ind,rep_ind)         = d2Cdt2(p_ind,rep_ind,IndexForLatency);
                    ModelPredictions.(CurrentCompartment).(CurrentModel).DeltaCoverC(p_ind,rep_ind)    = DeltaCoverC(p_ind,rep_ind,IndexForLatency);
                    
                end                                           
            end
        end                                
    end       
end

if 0
%%  Plot Showing that all models are not consistent with the data, using only MeasuredStats  
%%% Models are consistent only ABOVE the dotted line (Latency <= Peak time), but all data is below it 
% Note: Measured time to peak is similar for every compartment.
Y      = MeasuredStats.(CurrentCompartment).dCdt_TimeToPeak.Matrix;
Y_Mean = nanmean(Y,2);
Y_SEM  = nanstd(Y,[],2) ./ sqrt(sum(~isnan(Y),2));
Z      = MeasuredStats.(CurrentCompartment).d2Cdt2_TimeToPeak.Matrix;
Z_Mean = nanmean(Z,2);
Z_SEM  = nanstd(Z,[],2) ./ sqrt(sum(~isnan(Z),2));
W      = MeasuredStats.(CurrentCompartment).DeltaCoverC_TimeToPeak.Matrix;
W_Mean = nanmean(W,2);
W_SEM  = nanstd(W,[],2) ./ sqrt(sum(~isnan(W),2));
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};  
    figure('name',CurrentCompartment);    
    
    X      = MeasuredStats.(CurrentCompartment).Latency.Matrix;
    X_Mean = nanmean(X,2);
    X_SEM  = nanstd(X,[],2) ./ sqrt(sum(~isnan(X),2));
    
    errorbar(X_Mean,Y_Mean,Y_SEM,Y_SEM,X_SEM,X_SEM, 'color','k','linestyle','none','marker','.');   hold on; %  errorbar(x,y,yneg,ypos,xneg,xpos)
    errorbar(X_Mean,Z_Mean,Z_SEM,Z_SEM,X_SEM,X_SEM, 'color','g','linestyle','none','marker','.');   hold on; %  errorbar(x,y,yneg,ypos,xneg,xpos)
    errorbar(X_Mean,W_Mean,W_SEM,W_SEM,X_SEM,X_SEM, 'color','b','linestyle','none','marker','.');   hold on; %  errorbar(x,y,yneg,ypos,xneg,xpos)            
    
    hold on; plot([0 100],[0 100],'k:'); axis([-3 51 -3 51])
    ylabel('Time to peak [s]')
    xlabel('Latency [s]')
    legend('dC/dt','d^2C/dt^2','\DeltaC/C')
end

if isfield(Data.(CurrentStrain_Slow).(PulsesNames{1}),'Soma_LNmodelPredictedDeltaFoverF')
    % LN model plot
    for c_ind = 1:length(CompartmentFieldNames)
        CurrentCompartment = CompartmentFieldNames{c_ind};  
        figure('name',[CurrentCompartment,' , LN model LAtency predictions']);    

%         Y      = ModelPredictions.(CurrentCompartment).;
%         Y_Mean = nanmean(Y,2);
%         Y_SEM  = nanstd(Y,[],2) ./ sqrt(sum(~isnan(Y),2));

        X      = MeasuredStats.(CurrentCompartment).Latency.Matrix;
        X_Mean = nanmean(X,2);
        X_SEM  = nanstd(X,[],2) ./ sqrt(sum(~isnan(X),2));

        errorbar(X_Mean,Y_Mean,Y_SEM,Y_SEM,X_SEM,X_SEM, 'color','k','linestyle','none','marker','.');   hold on; %  errorbar(x,y,yneg,ypos,xneg,xpos)

        hold on; plot([0 100],[0 100],'k:'); axis([-3 51 -3 51])
        ylabel('Predicted Latency [s]')
        xlabel('Latency [s]')
        legend('LN model')
    end
end

end
return

function Data = ComputeImpulseFunctions(Data)

%  Compute impulse functions for each fast transition pulse (unique one for each strain, condition and compartment)   
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);         

RelevantStrains = { 'str2_CX17256_FromD7_Fast'
                    'str2_CX17256_FromD6_Fast'
                    'str2_CX17256_FromD5_Fast'
                    'str2_CX17256_FromD4_Fast'
                    'str2_CX17256_Desensitized_FromD6_Fast'
                    'str2_CX17256_Desensitized_FromD5_Fast'
                    'str2_CX17256_Desensitized_FromD4_Fast'
                    'str2_CX17256_Desensitized_FromD3_Fast'
                    'str2_CX17255_AllFast'
                    'str2_Egl4_AllFast'
                    'str2_Egl4_ad450_AllFast' };

Ymax = 17.4; n=2.7; Kd= 307;   % Normalization: From DeltaFoverF to calcium traces [nM] (GCaMP5a, Akerboom et al., 2012)
WindowSize = 5;
b = ones(1,WindowSize)/WindowSize;
a = 1;
D = round(mean(grpdelay(b,a)));  

for s_ind = 1:length(RelevantStrains)
    CurrentStrain = RelevantStrains{s_ind};
    PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    NumberOfPulses      = length(PulsesNames);         
    for p_ind = 1:NumberOfPulses
        CurrentPulse  = PulsesNames{p_ind};
        if ~isfield(Data.(CurrentStrain),CurrentPulse)
            continue
        end
        TimeVec       = Data.(CurrentStrain).(CurrentPulse).TimeVec;
        time_zero_ind = find(TimeVec==0);
        DeltaFOverFLength = length(TimeVec);

        for c_ind = 1:NumberOfCompartment
            CurrentCompartment = CompartmentFieldNames{c_ind};
            Current_Activity_Field      = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed'];
            DeltaFOverF_StepResponseMat = Data.(CurrentStrain).(CurrentPulse).(Current_Activity_Field);
        
            % Calcium        
            Calcium_StepResponseMat          = CalculateCalciumFromDeltaFOverF_UsingGCaMPproperties (DeltaFOverF_StepResponseMat,  Kd, Ymax, n); 
            Calcium_StepResponseMat_Negative = CalculateCalciumFromDeltaFOverF_UsingGCaMPproperties (-DeltaFOverF_StepResponseMat, Kd, Ymax, n); 
            NegativeValuesIndices            = DeltaFOverF_StepResponseMat(:)<0;
            Calcium_StepResponseMat(NegativeValuesIndices) = -Calcium_StepResponseMat_Negative(NegativeValuesIndices);            
            Calcium_StepResponseMat = Calcium_StepResponseMat - repmat(nanmean(Calcium_StepResponseMat(:,10:20),2),1,size(Calcium_StepResponseMat,2)); % normalize by one second before pulse initiation
                                   
            Calcium_StepResponseMatFiltered                                            = filter(b, a, Calcium_StepResponseMat, [], 2);
            Calcium_StepResponseMatFiltered                                            = Calcium_StepResponseMatFiltered(:,D+1:end);
            Calcium_StepResponseMatFiltered(:,(DeltaFOverFLength-D):DeltaFOverFLength) = NaN;
            
            % Impulse function                      
            AverageCalciumStepReponse = nanmean(Calcium_StepResponseMatFiltered,1);
            ImpulseFunction           = [NaN diff(AverageCalciumStepReponse)]; 
            ImpulseFunction           = ImpulseFunction(time_zero_ind:end);
            
            Data.(CurrentStrain).(CurrentPulse).([CurrentCompartment,'_Calcium_StepResponse'])         = Calcium_StepResponseMat;                       
            Data.(CurrentStrain).(CurrentPulse).([CurrentCompartment,'_Calcium_StepResponseFiltered']) = Calcium_StepResponseMatFiltered;                       
            Data.(CurrentStrain).(CurrentPulse).([CurrentCompartment,'_ImpulseFunction_NoScaling'])    = ImpulseFunction;                       
        end        
    end    
end

%% Test the influence of GCaMP5a Kernel. It will later be shown that our results are not sensitive to GCaMP5a kinetics.
%  This is expected since GCaMP5a kinetics are much faster (<1sec) than the timescale of activity initiation in response to slow inputs is slow (10-20sec)   
GCaMP_HalfLife  = 0.28;                  % seconds (for GCaMP5a, see Larsch et al. Cell Reports, 2015)    
GCaMP_DecayTime = GCaMP_HalfLife/log(2); % ~0.4 seconds
t_Kernel        = 0:0.1:1.6;             
Kernel          = exp(-t_Kernel./GCaMP_DecayTime); Kernel = Kernel/sum(Kernel);

for s_ind = 1:length(RelevantStrains)
    CurrentStrain = RelevantStrains{s_ind};
    PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    NumberOfPulses      = length(PulsesNames);         
    for p_ind = 1:NumberOfPulses
        CurrentPulse  = PulsesNames{p_ind};
        if ~isfield(Data.(CurrentStrain),CurrentPulse)
            continue
        end
        TimeVec       = Data.(CurrentStrain).(CurrentPulse).TimeVec;
        time_zero_ind = find(TimeVec==0);
        DeltaFOverFLength = length(TimeVec);

        for c_ind = 1:NumberOfCompartment
            CurrentCompartment = CompartmentFieldNames{c_ind};
            Current_Activity_Field      = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed'];
            DeltaFOverF_StepResponseMat = Data.(CurrentStrain).(CurrentPulse).(Current_Activity_Field);
        
            % Calcium BEFORE GCaMP kernel correction        
            Calcium_StepResponseMat          = CalculateCalciumFromDeltaFOverF_UsingGCaMPproperties (DeltaFOverF_StepResponseMat,  Kd, Ymax, n); 
            Calcium_StepResponseMat_Negative = CalculateCalciumFromDeltaFOverF_UsingGCaMPproperties (-DeltaFOverF_StepResponseMat, Kd, Ymax, n); 
            NegativeValuesIndices            = DeltaFOverF_StepResponseMat(:)<0;
            Calcium_StepResponseMat(NegativeValuesIndices) = -Calcium_StepResponseMat_Negative(NegativeValuesIndices);            
            Calcium_StepResponseMat = Calcium_StepResponseMat - repmat(nanmean(Calcium_StepResponseMat(:,10:20),2),1,size(Calcium_StepResponseMat,2)); % normalize by one second before pulse initiation
                        
            Calcium_StepResponseMatFiltered                                            = filter(b, a, Calcium_StepResponseMat, [], 2);
            Calcium_StepResponseMatFiltered                                            = Calcium_StepResponseMatFiltered(:,D+1:end);
            Calcium_StepResponseMatFiltered(:,(DeltaFOverFLength-D):DeltaFOverFLength) = NaN;
            
            AverageCalciumStepReponse = nanmean(Calcium_StepResponseMatFiltered,1);

            % Deconvolve GCaMP kernel before extraction of impulse function   
            AverageCalciumStepReponse_GCaMP_Deconv = deconv(AverageCalciumStepReponse,Kernel); 
            AverageCalciumStepReponse_GCaMP_Deconv(end+1:length(AverageCalciumStepReponse)) = NaN;            
            % figure; plot(AverageCalciumStepReponse); hold on; plot(AverageCalciumStepReponse_GCaMP_Deconv)  
            
            % Impulse function                           
            ImpulseFunction           = [NaN diff(AverageCalciumStepReponse_GCaMP_Deconv)]; 
            ImpulseFunction           = ImpulseFunction(time_zero_ind:end);
            
            Data.(CurrentStrain).(CurrentPulse).([CurrentCompartment,'_ImpulseFunction_NoScaling_GCaMP_Deconv']) = ImpulseFunction;                       
        end        
    end    
end

%% Saturation effect Screen: here we extract impulse function for various saturation levels. It will later be shown that our results are not sensitive to saturation level
%  Since the saturation level (Ymax) cannot be lower than the measured maximal DeltaF/F 
Ymax_vec = Ymax*[1/3 1/2 1 2 3];   % from 3 lower to 3x higher than in vitro measurements

for s_ind = 1:length(RelevantStrains)
    CurrentStrain = RelevantStrains{s_ind};
    PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    NumberOfPulses      = length(PulsesNames);         
    for p_ind = 1:NumberOfPulses
        CurrentPulse  = PulsesNames{p_ind};
        if ~isfield(Data.(CurrentStrain),CurrentPulse)
            continue
        end
        TimeVec       = Data.(CurrentStrain).(CurrentPulse).TimeVec;
        time_zero_ind = find(TimeVec==0);
        DeltaFOverFLength = length(TimeVec);

        for c_ind = 1:NumberOfCompartment
            CurrentCompartment = CompartmentFieldNames{c_ind};
            Current_Activity_Field      = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed'];
            DeltaFOverF_StepResponseMat = Data.(CurrentStrain).(CurrentPulse).(Current_Activity_Field);
            % Initialization
            ImpulseFunctionLength = length(Data.(CurrentStrain).(CurrentPulse).([CurrentCompartment,'_ImpulseFunction_NoScaling']));
            Data.(CurrentStrain).(CurrentPulse).([CurrentCompartment,'_ImpulseFunction_NoScaling_SaturationScreen']) = ones(length(Ymax_vec),ImpulseFunctionLength,'single')*NaN;
            for y_ind = 1:length(Ymax_vec)
                CurrentYmax = Ymax_vec(y_ind);

                % Calcium        
                Calcium_StepResponseMat          = CalculateCalciumFromDeltaFOverF_UsingGCaMPproperties (DeltaFOverF_StepResponseMat,  Kd, CurrentYmax, n); 
                Calcium_StepResponseMat_Negative = CalculateCalciumFromDeltaFOverF_UsingGCaMPproperties (-DeltaFOverF_StepResponseMat, Kd, CurrentYmax, n); 
                NegativeValuesIndices            = DeltaFOverF_StepResponseMat(:)<0;
                Calcium_StepResponseMat(NegativeValuesIndices) = -Calcium_StepResponseMat_Negative(NegativeValuesIndices);            
                Calcium_StepResponseMat = Calcium_StepResponseMat - repmat(nanmean(Calcium_StepResponseMat(:,10:20),2),1,size(Calcium_StepResponseMat,2)); % normalize by one second before pulse initiation

                Calcium_StepResponseMatFiltered                                            = filter(b, a, Calcium_StepResponseMat, [], 2);
                Calcium_StepResponseMatFiltered                                            = Calcium_StepResponseMatFiltered(:,D+1:end);
                Calcium_StepResponseMatFiltered(:,(DeltaFOverFLength-D):DeltaFOverFLength) = NaN;

                % Impulse function                      
                AverageCalciumStepReponse = nanmean(Calcium_StepResponseMatFiltered,1);
                ImpulseFunction           = [NaN diff(AverageCalciumStepReponse)]; 
                ImpulseFunction           = ImpulseFunction(time_zero_ind:end);

                Data.(CurrentStrain).(CurrentPulse).([CurrentCompartment,'_ImpulseFunction_NoScaling_SaturationScreen'])(y_ind,:) = ImpulseFunction;                       
            end
        end        
    end    
end

return

function X = CalculateCalciumFromDeltaFOverF_UsingGCaMPproperties (Y, Kd, Ymax, n)
% Y =  X^n/(X^n + K^n) * Ymax. 
% Y is deltaFoverF,  X is calcium concentration
% Units: 
%    Y and Ymax are uniteless: deltaF/F
%    n is unitless
%    Kd is concentration [nM]
%
% For GCaMP5a:  Ymax = 17.4; n=2.7; Kd= 307;   % figure; X=0:1e3; plot(X, (X.^n./(X.^n + Kd^n) * Ymax) )

X = Kd ./ ((Ymax./Y - 1).^(1/n));

% If Y<<Ymax -->  X = (K/Ymax^(1/n)) *Y.^(1/n)

return

function DeltaFoverF = ComputeDeltaFOverF_FromCalcium_UsingGCaMPproperties(Calcium, Kd, DeltaFoverF_max, n)
% Y =  X^n/(X^n + K^n) * Ymax. 
% Y is deltaFoverF,  X is calcium concentration
% Units: 
%    Y and Ymax are uniteless: deltaF/F
%    n is unitless
%    Kd is concentration [nM]
%
% % For GCaMP5a:  
% DeltaFoverF_max = 17.4; 
% n    = 2.7; 
% Kd   = 307; 

Calcium = double(Calcium);

DeltaFoverF = DeltaFoverF_max * Calcium.^n ./ (Calcium.^n + Kd.^n  );

return

function Data = PredictResponse_LNmodel(Data)

ReferenceForImpulseFunction.StrainName  = 'str2_CX17256_FromD6_Fast'; 
ReferenceForImpulseFunction.PulseName   = 'Butanone_d6_To_Buffer'; 
RelevantStrains = {     'str2_CX17256_FromD7_Fast'
                        'str2_CX17256_FromD6_Fast'
                        'str2_CX17256_FromD5_Fast'
                        'str2_CX17256_FromD4_Fast'
                        'str2_CX17256_FromD6_Slow'
                        'str2_CX17256_ToBuffer_Slow' };   % All str2_CX17256 experiments

Ymax = 17.4; n=2.7; Kd= 307;   % Normalization: From DeltaFoverF to calcium traces [nM] (GCaMP5a, Akerboom et al., 2012)
MaxNumberOfRepeats    = 40;    % free parameter. It just needs to be higher than the maximal number of repeats
ReferenceStrainName   = ReferenceForImpulseFunction.StrainName;
ReferencePulseName    = ReferenceForImpulseFunction.PulseName;
ReferenceInitialConcentration = Data.(ReferenceStrainName).PulsesInfo.InitialConcentration; % Units of dilution

CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);         

% Reference Impulse Function
ImpulseFunctionLength = length(Data.(ReferenceStrainName).(ReferencePulseName).Soma_ImpulseFunction_NoScaling);
ImpulseFunctionMat    = zeros(NumberOfCompartment,ImpulseFunctionLength)*NaN;
for c_ind = 1:NumberOfCompartment
    CurrentCompartment         = CompartmentFieldNames{c_ind};
    ImpulseFunctionMat(c_ind,:)= Data.(ReferenceStrainName).(ReferencePulseName).([CurrentCompartment,'_ImpulseFunction_NoScaling']);
end
ReferenceForImpulseFunction.ReferenceInitialConcentration = ReferenceInitialConcentration;
ReferenceForImpulseFunction.ImpulseFunctionMat            = ImpulseFunctionMat;

for s_ind = 1:length(RelevantStrains)
    CurrentStrain              = RelevantStrains{s_ind};
    PulsesNames                = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;     % units of microM
    NumberOfPulses             = length(PulsesNames);         
    CurrentInitialConcentration = Data.(CurrentStrain).PulsesInfo.InitialConcentration;  % units of microM
    % Extract concentration information
    if isfield(Data.(CurrentStrain),'DyeInformation')   % Slow transition, dye info available
        TimeVec = Data.(CurrentStrain).DyeInformation.TimeVectorInSec;
        TimeZeroIndex     = find(TimeVec==0,1);
        C       = Data.(CurrentStrain).DyeInformation.C_Smoothed;
    else                                                % Fast transition, dye info unavailable- build step concentration inputs     
        TimeVec           = Data.(CurrentStrain).(PulsesNames{1}).TimeVec;
        TimeZeroIndex     = find(TimeVec==0,1);
        C                 = ones(NumberOfPulses,MaxNumberOfRepeats,length(TimeVec),'single')*NaN;
        C(:,:,1:TimeZeroIndex-1)    = CurrentInitialConcentration;
        for p_ind=1:NumberOfPulses
            C(p_ind,:,TimeZeroIndex:end) = FinalConcentrationPerPulse(p_ind);
        end                            
    end   
    
    Data.(CurrentStrain).ReferenceForImpulseFunction_LNmodel            = ReferenceForImpulseFunction;
 
    for p_ind = 1:NumberOfPulses
        CurrentPulse  = PulsesNames{p_ind};
        if ~isfield(Data.(CurrentStrain),CurrentPulse)
            continue
        end
        NumberOfRepeats = size(Data.(CurrentStrain).(CurrentPulse).Soma_DeltaFOverF,1);
                    
        for c_ind = 1:NumberOfCompartment
            CurrentCompartment          = CompartmentFieldNames{c_ind};
            ImpulseFunction             = ImpulseFunctionMat(c_ind,:);
            ImpulseFunction             = ImpulseFunction(1:find(~isnan(ImpulseFunction),1,'last'));
            CurrentFieldName_Input       = [CurrentCompartment,'_LNmodelInput'];
            CurrentFieldName_Calcium     = [CurrentCompartment,'_LNmodelPredictedCalcium'];
            CurrentFieldName_DeltaFoverF = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF'];
            
            Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Input)       = zeros(NumberOfRepeats,length(TimeVec),'single')*NaN;
            Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Calcium)     = zeros(NumberOfRepeats,length(TimeVec),'single')*NaN;
            Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_DeltaFoverF) = zeros(NumberOfRepeats,length(TimeVec),'single')*NaN;
                                            
            for rep_ind = 1:NumberOfRepeats
                CurrentConcentration        = squeeze(C(p_ind,rep_ind,:));                      % decreasing
                CurrentConcentrationFlipped = max(CurrentConcentration) - CurrentConcentration; % increasing
                CurrentConcentrationFlipped(1:TimeZeroIndex-1) = 0;                             % supposed to be zeros anyway
                CurrentInput   = CurrentConcentrationFlipped / ReferenceInitialConcentration;   % Normalized with respect to the amplitude="1" defined for the Impulse Function
                
                % Linear kernel
                PredictedCalciumResponse = conv(CurrentInput,ImpulseFunction);          % using Impulse Function calculated from reference step function 
                PredictedCalciumResponse = PredictedCalciumResponse(1:length(TimeVec));                
                
                % Non-linear transformation: % Normalization From calcium traces [nM] to DeltaFoverF  (GCaMP5a, Akerboom et al., 2012) 
                PredictedResponse_DeltaFoverF = ComputeDeltaFOverF_FromCalcium_UsingGCaMPproperties(PredictedCalciumResponse, Kd, Ymax, n);
                
                Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Input)(rep_ind,:)       = CurrentInput;
                Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Calcium)(rep_ind,:)     = PredictedCalciumResponse;
                Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_DeltaFoverF)(rep_ind,:) = PredictedResponse_DeltaFoverF;
            end            
        end   
    end        
end


%% Test the influence of GCaMP5a Kernel. It will later be shown that our results are not sensitive to GCaMP5a kinetics.
%  This is expected since GCaMP5a kinetics are much faster (<1sec) than the timescale of activity initiation in response to slow inputs is slow (10-20sec)   
GCaMP_HalfLife  = 0.28;                  % seconds (for GCaMP5a, see Larsch et al. Cell Reports, 2015)    
GCaMP_DecayTime = GCaMP_HalfLife/log(2); % ~0.4 seconds
t_Kernel        = 0:0.1:1.6;             
Kernel          = exp(-t_Kernel./GCaMP_DecayTime); Kernel = Kernel/sum(Kernel);

% Reference Impulse Function
ImpulseFunctionMat_GCaMP_Deconv    = zeros(NumberOfCompartment,ImpulseFunctionLength)*NaN;
for c_ind = 1:NumberOfCompartment
    CurrentCompartment         = CompartmentFieldNames{c_ind};
    ImpulseFunctionMat_GCaMP_Deconv(c_ind,:)= Data.(ReferenceStrainName).(ReferencePulseName).([CurrentCompartment,'_ImpulseFunction_NoScaling_GCaMP_Deconv']);
end
ReferenceForImpulseFunction_GCaMP_Deconv.ReferenceInitialConcentration = ReferenceInitialConcentration;
ReferenceForImpulseFunction_GCaMP_Deconv.ImpulseFunctionMat            = ImpulseFunctionMat_GCaMP_Deconv;

for s_ind = 1:length(RelevantStrains)
    CurrentStrain              = RelevantStrains{s_ind};
    PulsesNames                = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;     % units of microM
    NumberOfPulses             = length(PulsesNames);         
    CurrentInitialConcentration = Data.(CurrentStrain).PulsesInfo.InitialConcentration;  % units of microM
    % Extract concentration information
    if isfield(Data.(CurrentStrain),'DyeInformation')   % Slow transition, dye info available
        TimeVec = Data.(CurrentStrain).DyeInformation.TimeVectorInSec;
        TimeZeroIndex     = find(TimeVec==0,1);
        C       = Data.(CurrentStrain).DyeInformation.C_Smoothed;
    else                                                % Fast transition, dye info unavailable- build step concentration inputs     
        TimeVec           = Data.(CurrentStrain).(PulsesNames{1}).TimeVec;
        TimeZeroIndex     = find(TimeVec==0,1);
        C                 = ones(NumberOfPulses,MaxNumberOfRepeats,length(TimeVec),'single')*NaN;
        C(:,:,1:TimeZeroIndex-1)    = CurrentInitialConcentration;
        for p_ind=1:NumberOfPulses
            C(p_ind,:,TimeZeroIndex:end) = FinalConcentrationPerPulse(p_ind);
        end                            
    end   
    
    Data.(CurrentStrain).ReferenceForImpulseFunction_LNmodel_GCaMP_Deconv  = ReferenceForImpulseFunction_GCaMP_Deconv;
 
    for p_ind = 1:NumberOfPulses
        CurrentPulse  = PulsesNames{p_ind};
        if ~isfield(Data.(CurrentStrain),CurrentPulse)
            continue
        end
        NumberOfRepeats = size(Data.(CurrentStrain).(CurrentPulse).Soma_DeltaFOverF,1);
                    
        for c_ind = 1:NumberOfCompartment
            CurrentCompartment          = CompartmentFieldNames{c_ind};
            ImpulseFunction             = ImpulseFunctionMat_GCaMP_Deconv(c_ind,:);
            ImpulseFunction             = ImpulseFunction(1:find(~isnan(ImpulseFunction),1,'last'));
            CurrentFieldName_Calcium     = [CurrentCompartment,'_LNmodelPredictedCalcium_GCaMP_Deconv'];
            CurrentFieldName_DeltaFoverF = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF_GCaMP_Deconv'];
            
            Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Calcium)     = zeros(NumberOfRepeats,length(TimeVec),'single')*NaN;
            Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_DeltaFoverF) = zeros(NumberOfRepeats,length(TimeVec),'single')*NaN;
                                            
            for rep_ind = 1:NumberOfRepeats
                CurrentConcentration        = squeeze(C(p_ind,rep_ind,:));                      % decreasing
                CurrentConcentrationFlipped = max(CurrentConcentration) - CurrentConcentration; % increasing
                CurrentConcentrationFlipped(1:TimeZeroIndex-1) = 0;                             % supposed to be zeros anyway
                CurrentInput   = CurrentConcentrationFlipped / ReferenceInitialConcentration;   % Normalized with respect to the amplitude="1" defined for the Impulse Function
                
                % Linear kernel
                PredictedCalciumResponse = conv(CurrentInput,ImpulseFunction);          % using Impulse Function calculated from reference step function 
                PredictedCalciumResponse = PredictedCalciumResponse(1:length(TimeVec));                
                
                % GCaMP kernel
                PredictedCalciumResponse_GCaMP_Conv = conv(PredictedCalciumResponse,Kernel);       % correction for (fast) GCaMP kinetics 
                PredictedCalciumResponse_GCaMP_Conv = PredictedCalciumResponse_GCaMP_Conv(1:length(TimeVec));                

                % Non-linear transformation: % Normalization From calcium traces [nM] to DeltaFoverF  (GCaMP5a, Akerboom et al., 2012) 
                PredictedResponse_DeltaFoverF = ComputeDeltaFOverF_FromCalcium_UsingGCaMPproperties(PredictedCalciumResponse_GCaMP_Conv, Kd, Ymax, n);
                
                Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Calcium)(rep_ind,:)     = PredictedCalciumResponse;
                Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_DeltaFoverF)(rep_ind,:) = PredictedResponse_DeltaFoverF;
            end            
        end   
    end        
end


%% Here we extract LN model predictions for various saturation levels. It will later be shown that our results are not sensitive to saturation level
Ymax_vec = Ymax*[1/3 1/2 1 2 3];   % from 3 lower to 3x higher than in vitro measurements

% Reference Impulse Function
ImpulseFunctionMat_SaturationScreen    = zeros(length(Ymax_vec),NumberOfCompartment,ImpulseFunctionLength,'single')*NaN;
for c_ind = 1:NumberOfCompartment
    CurrentCompartment         = CompartmentFieldNames{c_ind};
    ImpulseFunctionMat_SaturationScreen(:,c_ind,:)= Data.(ReferenceStrainName).(ReferencePulseName).([CurrentCompartment,'_ImpulseFunction_NoScaling_SaturationScreen']);
end

ReferenceForImpulseFunction_SaturationScreen.ReferenceInitialConcentration = ReferenceInitialConcentration;
ReferenceForImpulseFunction_SaturationScreen.ImpulseFunctionMat            = ImpulseFunctionMat_SaturationScreen;

for s_ind = 1:length(RelevantStrains)
    CurrentStrain              = RelevantStrains{s_ind};
    PulsesNames                = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;     % units of microM
    NumberOfPulses             = length(PulsesNames);         
    CurrentInitialConcentration = Data.(CurrentStrain).PulsesInfo.InitialConcentration;  % units of microM
    % Extract concentration information
    if isfield(Data.(CurrentStrain),'DyeInformation')   % Slow transition, dye info available
        TimeVec = Data.(CurrentStrain).DyeInformation.TimeVectorInSec;
        TimeZeroIndex     = find(TimeVec==0,1);
        C       = Data.(CurrentStrain).DyeInformation.C_Smoothed;
    else                                                % Fast transition, dye info unavailable- build step concentration inputs     
        TimeVec           = Data.(CurrentStrain).(PulsesNames{1}).TimeVec;
        TimeZeroIndex     = find(TimeVec==0,1);
        C                 = ones(NumberOfPulses,MaxNumberOfRepeats,length(TimeVec),'single')*NaN;
        C(:,:,1:TimeZeroIndex-1)    = CurrentInitialConcentration;
        for p_ind=1:NumberOfPulses
            C(p_ind,:,TimeZeroIndex:end) = FinalConcentrationPerPulse(p_ind);
        end                            
    end   
    
    Data.(CurrentStrain).ReferenceForImpulseFunction_LNmodel_SaturationScreen = ReferenceForImpulseFunction_SaturationScreen;
 
    for p_ind = 1:NumberOfPulses
        CurrentPulse  = PulsesNames{p_ind};
        if ~isfield(Data.(CurrentStrain),CurrentPulse)
            continue
        end
        NumberOfRepeats = size(Data.(CurrentStrain).(CurrentPulse).Soma_DeltaFOverF,1);
                    
        for c_ind = 1:NumberOfCompartment
            CurrentCompartment          = CompartmentFieldNames{c_ind};
            CurrentFieldName_Calcium     = [CurrentCompartment,'_LNmodelPredictedCalcium_SaturationScreen'];
            CurrentFieldName_DeltaFoverF = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF_SaturationScreen'];
            
            Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Calcium)     = zeros(length(Ymax_vec),NumberOfRepeats,length(TimeVec),'single')*NaN;
            Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_DeltaFoverF) = zeros(length(Ymax_vec),NumberOfRepeats,length(TimeVec),'single')*NaN;
            
            for y_ind=1:length(Ymax_vec) 
                CurrentYmax = Ymax_vec(y_ind);
                ImpulseFunction             = squeeze(ImpulseFunctionMat_SaturationScreen(y_ind,c_ind,:));
                ImpulseFunction             = ImpulseFunction(1:find(~isnan(ImpulseFunction),1,'last'));
                for rep_ind = 1:NumberOfRepeats
                    CurrentConcentration        = squeeze(C(p_ind,rep_ind,:));                      % decreasing
                    CurrentConcentrationFlipped = max(CurrentConcentration) - CurrentConcentration; % increasing
                    CurrentConcentrationFlipped(1:TimeZeroIndex-1) = 0;                             % supposed to be zeros anyway
                    CurrentInput   = CurrentConcentrationFlipped / ReferenceInitialConcentration;   % Normalized with respect to the amplitude="1" defined for the Impulse Function

                    % Linear kernel
                    PredictedCalciumResponse = conv(CurrentInput,ImpulseFunction);          % using Impulse Function calculated from reference step function 
                    PredictedCalciumResponse = PredictedCalciumResponse(1:length(TimeVec));                

                    % Non-linear transformation: % Normalization From calcium traces [nM] to DeltaFoverF  (GCaMP5a, Akerboom et al., 2012) 
                    PredictedResponse_DeltaFoverF = ComputeDeltaFOverF_FromCalcium_UsingGCaMPproperties(PredictedCalciumResponse, Kd, CurrentYmax, n);

                    Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Calcium)(y_ind,rep_ind,:)     = PredictedCalciumResponse;
                    Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_DeltaFoverF)(y_ind,rep_ind,:) = PredictedResponse_DeltaFoverF;
                end     
            end
        end   
    end        
end

return

function Data = OptimizeRectifierForLNmodel(Data)

CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% FOR EACH MODEL:
% -- Compute linear kernel (impulse function) for different (b) 
% -- Compute responses to fast and slow transitions for different (b) 
% -- Plot Predictions against data
% -- quantify max(deltaF/F), Responsiveness and Latency -- compare to data    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% MODEL 1: Using a rectifier AFTER linear kernel
%
%                          Linear Kernel            Rectifier with threshold        GCaMP static non-linearity 
% INPUT(concentration)   ----------------->     x      ----------------->       y       ----------------->        OUTPUT (DeltaF/F)  
%
% 
% DEFINITION: Rectifying linear unit (input x and output y) with a threshold (b):
%     y = x-b at x>b
%     y = 0   at x<b
%
% NOTE!!  The impulse function is calculated using dx/dt upon step inputs, but dx/dt = dy/dt regrdless of b (up to discontinuity at t=0) 
%
% Strategy:
% -- In Responsive Fast transitions: DeltaF/F(t>0)>0  -->   y(t>0)>0 --> x(t>0)=y+b  --> x can be calculated for any OUTPUT, and the INPUT is known 
% -- Calculate linear kernel (impulse function) for different b   
% -- Find (b) that best fit the data
%
% Bottom line: adding significant rectification can properly increase response latency for slow inputs, but that also    
%              inaccurately inhibits responsiveness to both fast and slow inputs.   
% Generate plot:   Responsiveness accuracy vs. b_vec, and latency accuracy vs. b_vec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% MODEL 2: Using a rectifier PRIOR to linear kernel
%
%                     Rectifier with threshold            Linear Kernel           GCaMP static non-linearity 
% INPUT(concentration)   ----------------->     x      ----------------->       y       ----------------->        OUTPUT (DeltaF/F)  
%
% DEFINITION: Rectifying linear unit (input x and output y) with a threshold (b):
%     x = INPUT-b at INPUT>b
%     x = 0       at INPUT<b
%
% Strategy:
% -- OUTPUT is known --> Calculate y 
% -- INPUT is known  --> Calculate x for different (b) 
% -- Calculate linear kernel (impulse function) for different b   
% -- Find (b) that best fit the data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%
DeltaFoverF_max       = 17.4;   % (GCaMP5a, Akerboom et al., 2012)
n                     = 2.7; 
Kd                    = 307; 
RelevantStrains       = {'str2_CX17256_FromD6_Fast','str2_CX17256_FromD6_Slow'};
MaxNumberOfRepeats    = 40;         % free parameter. It just needs to be higher than the maximal number of repeats

% Parameters used for Measured Traces analysis of responsiveness and latency    
DeflectionPoint_LowDiffThreshold                = 5e-3;    % until you reach this lower threshold
DynamicRangeThresholdForResponsiveness.Soma     = 0.2;  
DynamicRangeThresholdForResponsiveness.Dendrite = 0.2;  
DynamicRangeThresholdForResponsiveness.Axon     = 0.2;     % axon responses tend to be more moderate >> requires lower threhold  
MinimalChangeRequiredForActivityInitiation      = 0.05;    % 5% deltaF/F        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Model 1
b_vec     = 0:100; 
ModelName = 'LNmodelWithRec_1';
for c_ind = 1:NumberOfCompartment
    CurrentCompartment     = CompartmentFieldNames{c_ind};    
    CurrentFieldName_Model = [CurrentCompartment,'_PredictedDeltaFoverF_',ModelName];
    
    % Impulse function based on step inputs. For model 1: (dx/dt = dy/dt)                        
    ImpulseFunction = Data.str2_CX17256_FromD6_Fast.Butanone_d6_To_Buffer.([CurrentCompartment,'_ImpulseFunction_NoScaling']);
    ImpulseFunction = ImpulseFunction(1:find(~isnan(ImpulseFunction),1,'last'));    
    
    for s_ind = 1:length(RelevantStrains)
        CurrentStrain  = RelevantStrains{s_ind};  
        PulsesNames    = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        NumberOfPulses = length(PulsesNames);     
        if s_ind == 1
            TimeVec    = Data.(CurrentStrain).Butanone_d6_To_Buffer.TimeVec;
        else   % (s_ind==2)
            TimeVec    = Data.(CurrentStrain).Butanone_d6_To_Buffer.Dye_TimeVec;
        end                

        for p_ind = 1:NumberOfPulses        
            CurrentPulse                = PulsesNames{p_ind};                    
            CurrentFieldName_Input      = [CurrentCompartment,'_LNmodelInput'];
            CurrentInputMat             = Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Input);
            NumberOfRepeats             = size(CurrentInputMat,1);                 
            PredictedCalciumResponseMat = ones(length(b_vec),NumberOfRepeats,length(TimeVec),'single')*NaN; %
            Data.(CurrentStrain).(CurrentPulse).(['b_vec_for_',ModelName]) = b_vec;
            
            for rep_ind = 1:NumberOfRepeats     
                CurrentInput = CurrentInputMat(rep_ind,:);
                %% Model 1: Linear Kernel between Input to x
                Predicted_x  = conv(CurrentInput,ImpulseFunction);     % That's the x. using Impulse Function calculated from reference step function  
                Predicted_x  = Predicted_x(1:length(CurrentInput));    % That's the x  
                for b_ind = 1:length(b_vec)    % b = rectifier threshold
                    b = b_vec(b_ind);        
                    % Predicted Calcium Response to slow inputs  
                    PredictedCalciumResponse                     = Predicted_x - b;  % That's the y    
                    PredictedCalciumResponse(Predicted_x<b)      = 0 ;               % That's the y          
                    PredictedCalciumResponseMat(b_ind,rep_ind,:) = PredictedCalciumResponse;                                                          
                end                
            end
            PredictedDeltaFoverFMat = single(ComputeDeltaFOverF_FromCalcium_UsingGCaMPproperties(PredictedCalciumResponseMat, Kd, DeltaFoverF_max, n));            
            Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Model) = PredictedDeltaFoverFMat;                          
        end        
    end    
end   % Model 1 specific

ModelName = 'LNmodelWithRec_1';
%%%% Model Stats %%%%
for c_ind = 1:NumberOfCompartment
    CurrentCompartment           = CompartmentFieldNames{c_ind};    
    CurrentFieldName_Model      = [CurrentCompartment,'_PredictedDeltaFoverF_',ModelName];    
    CurrentDynamicRangeThreshold = DynamicRangeThresholdForResponsiveness.(CurrentCompartment); % For stats analysis
    
    for s_ind = 1:length(RelevantStrains)
        CurrentStrain  = RelevantStrains{s_ind};  
        PulsesNames    = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        NumberOfPulses = length(PulsesNames);     
        if s_ind == 1
            TimeVec    = Data.(CurrentStrain).Butanone_d6_To_Buffer.TimeVec;
        else   % (s_ind==2)
            TimeVec    = Data.(CurrentStrain).Butanone_d6_To_Buffer.Dye_TimeVec;
        end
        TimeZeroIndex = find(TimeVec==0,1);
                
        StructureInitialization.Matrix = ones(length(b_vec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;
        
        Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_DynamicRange'])   = StructureInitialization;
        Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Responsiveness']) = StructureInitialization;
        Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_LatencyIndex'])   = StructureInitialization;
        Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Latency'])        = StructureInitialization;

        for p_ind = 1:NumberOfPulses        
            CurrentPulse                = PulsesNames{p_ind};                    
            CurrentFieldName_Input      = [CurrentCompartment,'_LNmodelInput'];
            CurrentInputMat             = Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Input);
            NumberOfRepeats             = size(CurrentInputMat,1);                 
            PredictedDeltaFoverFMat     = Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Model); %
                        
            DynamicRange   = nanmax(PredictedDeltaFoverFMat,[],3);
            Responsiveness = 100* (DynamicRange > CurrentDynamicRangeThreshold); 
            LatencyIndex   = DynamicRange*NaN; % initialization
            Latency        = DynamicRange*NaN; % initialization

            for rep_ind = 1:NumberOfRepeats     
                for b_ind = 1:length(b_vec)    % b = rectifier threshold
                    CurrentResponsiveness = Responsiveness(b_ind, rep_ind);                                 
                    if CurrentResponsiveness
                       % Compute latency 
                        CurrentVec     = squeeze(PredictedDeltaFoverFMat(b_ind,rep_ind,:))'; 
                        CurrentDiffVec = [NaN diff(CurrentVec)];

                        %% Find Latency by deflection point and by threshold value
                        % Deflection point by differential
                        IndexAfterDeflectionPoint = find(CurrentVec>CurrentDynamicRangeThreshold,1,'first');
                        IndexOfDeflectionPoint    = find(CurrentDiffVec(1:IndexAfterDeflectionPoint)<=DeflectionPoint_LowDiffThreshold,1,'last'); 
                        if TimeZeroIndex > IndexOfDeflectionPoint  
                           IndexOfDeflectionPoint = TimeZeroIndex;
                        end 

                        % Deflection point by max value with at least minimal differential            
                        IndexAfterDeflectionPoint = find(CurrentVec>CurrentDynamicRangeThreshold,1,'first');
                        IndexOfHighActivity       = find(CurrentVec(1:IndexAfterDeflectionPoint)<=MinimalChangeRequiredForActivityInitiation,1,'last')+1;
                        if IndexOfHighActivity< TimeZeroIndex
                            % assume error, assign maximum length
                             IndexOfHighActivity = length(CurrentVec);
                        end

                        %%%% Deflection point, by activity or by diff(activity)?
                        if isempty(IndexOfDeflectionPoint) || (IndexOfHighActivity < IndexOfDeflectionPoint) 
                            IndexForLatency = IndexOfHighActivity;
                        else
                            IndexForLatency = IndexOfDeflectionPoint;
                        end                        
                        LatencyIndex(b_ind,rep_ind)   = IndexForLatency;           
                        Latency(b_ind,rep_ind)        = TimeVec(IndexForLatency);                                 
                    end
                end
            end            
            Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_DynamicRange']).Matrix(:,p_ind,1:NumberOfRepeats)   = DynamicRange;                    
            Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Responsiveness']).Matrix(:,p_ind,1:NumberOfRepeats) = Responsiveness; 
            Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_LatencyIndex']).Matrix(:,p_ind,1:NumberOfRepeats)   = LatencyIndex; 
            Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Latency']).Matrix(:,p_ind,1:NumberOfRepeats)        = Latency;                                  
        end
        
        RelevantFields = {[CurrentCompartment,'_DynamicRange'],[CurrentCompartment,'_Responsiveness'],...
                          [CurrentCompartment,'_LatencyIndex'],[CurrentCompartment,'_Latency']};
        for f_ind = 1:length(RelevantFields)
            CurrentField = RelevantFields{f_ind};
            BigMat          = Data.(CurrentStrain).Stats.(ModelName).(CurrentField).Matrix;
            NumberOfRepeats = sum(~isnan(BigMat),3);
            Data.(CurrentStrain).Stats.(ModelName).(CurrentField).Mean   = nanmean(BigMat,3);
            Data.(CurrentStrain).Stats.(ModelName).(CurrentField).Median = nanmedian(BigMat,3);
            Data.(CurrentStrain).Stats.(ModelName).(CurrentField).STD    = nanstd(BigMat,[],3);
            Data.(CurrentStrain).Stats.(ModelName).(CurrentField).STD    = nanstd(BigMat,[],3) ./ sqrt(NumberOfRepeats);                       
        end          
    end    
end

if 0
% The graphs below shows that:
%    Low  'b' is required for fitting the responsiveness and timing data for FAST odor transitions (b<30)    
%    High 'b' is required for fitting the responsiveness and timing data for SLOW odor transitions (67<b<83)   
%    --> NO OVERLAP --> NO b can fit the data of both Fast and Slow transitions  

%%%% Plot Model Stats %%%%
ModelName          = 'LNmodelWithRec_1';
CurrentCompartment = 'Soma';
b_vec_ticks        = 1:10:101;
b_vec              = Data.str2_CX17256_FromD6_Fast.Butanone_d6_To_Buffer.(['b_vec_for_',ModelName]);

for s_ind = 1:length(RelevantStrains)
    CurrentStrain  = RelevantStrains{s_ind};  
    XtickLabels    = Data.(CurrentStrain).PulsesInfo.FinalConcentration; 
    figure('name',[CurrentStrain,' , Responsiveness'],'position',[ 455   105   990   865]); 
        Matrix = Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Responsiveness']).Mean;
        imagesc(Matrix);
        set(gca,'ytick',b_vec_ticks,'yticklabel',b_vec(b_vec_ticks),'xtick',1:length(XtickLabels),'xticklabel',XtickLabels)
    figure('name',[CurrentStrain,' , Latency'],'position',[ 455   105   990   865]); 
        Matrix = Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Latency']).Mean;
        imagesc_withoutNaNs(1:size(Matrix,2), 1:size(Matrix,1), Matrix, [0 ceil(nanmax(Matrix(:)))], [], [1 1 1]); colorbar;
        set(gca,'ytick',b_vec_ticks,'yticklabel',b_vec(b_vec_ticks),'xtick',1:length(XtickLabels),'xticklabel',XtickLabels)
end
    
%%%% Plot traces - for illustrator show only:  buffer and 2e-7  %%%%    
ModelName   = 'LNmodelWithRec_1';
b_vec_ticks = 1:20:101;
c_ind                   = 1;
b_vec                   = Data.str2_CX17256_FromD6_Fast.Butanone_d6_To_Buffer.(['b_vec_for_',ModelName]);
CurrentCompartment      = CompartmentFieldNames{c_ind};    
CurrentFieldName_Model  = [CurrentCompartment,'_PredictedDeltaFoverF_',ModelName];
CurrentFieldName_Data   = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed'];

for s_ind = 1:length(RelevantStrains)
    CurrentStrain  = RelevantStrains{s_ind};        
    PulsesNames    = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    YlabelNames    = num2cell(Data.(CurrentStrain).PulsesInfo.FinalConcentration);
    NumberOfPulses = length(PulsesNames);     
    if s_ind == 1
        TimeVec          = Data.(CurrentStrain).Butanone_d6_To_Buffer.TimeVec;
        TimeVecPredicted = TimeVec;
    else   % (s_ind==2)
        TimeVec          = Data.(CurrentStrain).Butanone_d6_To_Buffer.TimeVec;
        TimeVecPredicted = Data.(CurrentStrain).Butanone_d6_To_Buffer.Dye_TimeVec;
    end

    figure('name',[CurrentCompartment,' , ',CurrentStrain,'. Measured (k) vs. Model1 (r)']); 
    gcaHandles     = [];    
    CompartmentMax = 0;            

    for p_ind = 1:NumberOfPulses        
        CurrentPulse          = PulsesNames{p_ind};                    
        CurrentMeasuredMat    = Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Data);   % 2D
        CurrentPredictedMat3D = Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Model);  % 3D

        for b_ind = 1:length(b_vec_ticks)    
            b_ind_InVec = b_vec_ticks(b_ind);
            b           = b_vec(b_ind_InVec);
            subplot(NumberOfPulses,length(b_vec_ticks),(p_ind-1)*length(b_vec_ticks)+b_ind);
            CurrentPredictedMat = squeeze(CurrentPredictedMat3D(b_ind_InVec,:,:));             % 2D

            plot(TimeVec         , CurrentMeasuredMat',           '-','color',[0.7 0.7 0.7]); hold on;
            plot(TimeVecPredicted, CurrentPredictedMat',          '-','color',[1 0.7 0.7]);   hold on;
            plot(TimeVec         , nanmean(CurrentMeasuredMat,1), '-','color','k'); hold on;
            plot(TimeVecPredicted, nanmean(CurrentPredictedMat,1),'-','color','r'); hold on;
            xlim([-2 min([TimeVec(end) 70])]);             
            if p_ind == 1
                title(b);
            end                
            if b_ind == 1
                ylabel(YlabelNames{p_ind});
            end
            if p_ind~=NumberOfPulses
                set(gca,'xticklabel',[]);
            end
            if b_ind~=1
                set(gca,'yticklabel',[]);
            end
            gcaHandles     = [gcaHandles gca];
            CompartmentMax = max([CompartmentMax CurrentMeasuredMat(:)' CurrentPredictedMat(:)']);                               
        end                       
    end
    YLIMITS = [-0.2 CompartmentMax*1.1];
    set(gcaHandles,'ylim',YLIMITS)           
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Model 2
b_vec     = (0:99)/100; 
ModelName = 'LNmodelWithRec_2';

for c_ind = 1:NumberOfCompartment
    CurrentCompartment     = CompartmentFieldNames{c_ind};    
    CurrentFieldName_Model = [CurrentCompartment,'_PredictedDeltaFoverF_',ModelName];
    
    % Impulse function based on step input to buffer. For model 2 it needs to be scaled by: (1-b)                      
    ImpulseFunction_Original = Data.str2_CX17256_FromD6_Fast.Butanone_d6_To_Buffer.([CurrentCompartment,'_ImpulseFunction_NoScaling']);
    ImpulseFunction_Original = ImpulseFunction_Original(1:find(~isnan(ImpulseFunction_Original),1,'last'));    
    
    for s_ind = 1:length(RelevantStrains)
        CurrentStrain  = RelevantStrains{s_ind};  
        PulsesNames    = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        NumberOfPulses = length(PulsesNames);     
        if s_ind == 1
            TimeVec    = Data.(CurrentStrain).Butanone_d6_To_Buffer.TimeVec;
        else   % (s_ind==2)
            TimeVec    = Data.(CurrentStrain).Butanone_d6_To_Buffer.Dye_TimeVec;
        end                

        for p_ind = 1:NumberOfPulses        
            CurrentPulse                = PulsesNames{p_ind};                    
            CurrentFieldName_Input      = [CurrentCompartment,'_LNmodelInput'];
            CurrentInputMat             = Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Input);
            NumberOfRepeats             = size(CurrentInputMat,1);                 
            PredictedCalciumResponseMat = ones(length(b_vec),NumberOfRepeats,length(TimeVec),'single')*NaN; %
            Data.(CurrentStrain).(CurrentPulse).(['b_vec_for_',ModelName]) = b_vec;
            
            for rep_ind = 1:NumberOfRepeats     
                CurrentInput = CurrentInputMat(rep_ind,:);                
                
                %% Model 2: Linear Kernel between x to y(calcium)                                
                for b_ind = 1:length(b_vec)    % b = rectifier threshold
                    b                           = b_vec(b_ind);        
                    % x = Input-b at Input>b, therefore the reference step function was between 0 and (1-b) and not between 0 and 1
                    CurrentImpulseFunction      = ImpulseFunction_Original /(1-b);  
                    
                    Predicted_x                 = CurrentInput-b;
                    Predicted_x(CurrentInput<b) = 0 ;          
                    PredictedCalciumResponse    = conv(Predicted_x,CurrentImpulseFunction);            % That's the y. 
                    PredictedCalciumResponse    = PredictedCalciumResponse(1:length(CurrentInput));    % That's the y                                                                               
                    PredictedCalciumResponseMat(b_ind,rep_ind,:) = PredictedCalciumResponse;                                                          
                end
            end
            PredictedDeltaFoverFMat = single(ComputeDeltaFOverF_FromCalcium_UsingGCaMPproperties(PredictedCalciumResponseMat, Kd, DeltaFoverF_max, n));            
            Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Model) = PredictedDeltaFoverFMat;                          
        end        
    end    
end   % Model 2 specific

ModelName = 'LNmodelWithRec_2';
%%%% Model Stats %%%%
for c_ind = 1:NumberOfCompartment
    CurrentCompartment           = CompartmentFieldNames{c_ind};    
    CurrentFieldName_Model      = [CurrentCompartment,'_PredictedDeltaFoverF_',ModelName];    
    CurrentDynamicRangeThreshold = DynamicRangeThresholdForResponsiveness.(CurrentCompartment); % For stats analysis
    
    for s_ind = 1:length(RelevantStrains)
        CurrentStrain  = RelevantStrains{s_ind};  
        PulsesNames    = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        NumberOfPulses = length(PulsesNames);     
        if s_ind == 1
            TimeVec    = Data.(CurrentStrain).Butanone_d6_To_Buffer.TimeVec;
        else   % (s_ind==2)
            TimeVec    = Data.(CurrentStrain).Butanone_d6_To_Buffer.Dye_TimeVec;
        end
        TimeZeroIndex = find(TimeVec==0,1);
                
        StructureInitialization.Matrix = ones(length(b_vec),NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;
        
        Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_DynamicRange'])   = StructureInitialization;
        Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Responsiveness']) = StructureInitialization;
        Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_LatencyIndex'])   = StructureInitialization;
        Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Latency'])        = StructureInitialization;

        for p_ind = 1:NumberOfPulses        
            CurrentPulse                = PulsesNames{p_ind};                    
            CurrentFieldName_Input      = [CurrentCompartment,'_LNmodelInput'];
            CurrentInputMat             = Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Input);
            NumberOfRepeats             = size(CurrentInputMat,1);                 
            PredictedDeltaFoverFMat     = Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Model); %
                        
            DynamicRange   = nanmax(PredictedDeltaFoverFMat,[],3);
            Responsiveness = 100* (DynamicRange > CurrentDynamicRangeThreshold); 
            LatencyIndex   = DynamicRange*NaN; % initialization
            Latency        = DynamicRange*NaN; % initialization

            for rep_ind = 1:NumberOfRepeats     
                for b_ind = 1:length(b_vec)    % b = rectifier threshold
                    CurrentResponsiveness = Responsiveness(b_ind, rep_ind);                                 
                    if CurrentResponsiveness
                       % Compute latency 
                        CurrentVec     = squeeze(PredictedDeltaFoverFMat(b_ind,rep_ind,:))'; 
                        CurrentDiffVec = [NaN diff(CurrentVec)];

                        %% Find Latency by deflection point and by threshold value
                        % Deflection point by differential
                        IndexAfterDeflectionPoint = find(CurrentVec>CurrentDynamicRangeThreshold,1,'first');
                        IndexOfDeflectionPoint    = find(CurrentDiffVec(1:IndexAfterDeflectionPoint)<=DeflectionPoint_LowDiffThreshold,1,'last'); 
                        if TimeZeroIndex > IndexOfDeflectionPoint  
                           IndexOfDeflectionPoint = TimeZeroIndex;
                        end 

                        % Deflection point by max value with at least minimal differential            
                        IndexAfterDeflectionPoint = find(CurrentVec>CurrentDynamicRangeThreshold,1,'first');
                        IndexOfHighActivity       = find(CurrentVec(1:IndexAfterDeflectionPoint)<=MinimalChangeRequiredForActivityInitiation,1,'last')+1;
                        if IndexOfHighActivity< TimeZeroIndex
                            % assume error, assign maximum length
                             IndexOfHighActivity = length(CurrentVec);
                        end

                        %%%% Deflection point, by activity or by diff(activity)?
                        if isempty(IndexOfDeflectionPoint) || (IndexOfHighActivity < IndexOfDeflectionPoint) 
                            IndexForLatency = IndexOfHighActivity;
                        else
                            IndexForLatency = IndexOfDeflectionPoint;
                        end                        
                        LatencyIndex(b_ind,rep_ind)   = IndexForLatency;           
                        Latency(b_ind,rep_ind)        = TimeVec(IndexForLatency);                                 
                    end
                end
            end            
            Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_DynamicRange']).Matrix(:,p_ind,1:NumberOfRepeats)   = DynamicRange;                    
            Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Responsiveness']).Matrix(:,p_ind,1:NumberOfRepeats) = Responsiveness; 
            Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_LatencyIndex']).Matrix(:,p_ind,1:NumberOfRepeats)   = LatencyIndex; 
            Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Latency']).Matrix(:,p_ind,1:NumberOfRepeats)        = Latency; 
                                 
        end
        
        RelevantFields = {[CurrentCompartment,'_DynamicRange'],[CurrentCompartment,'_Responsiveness'],...
                          [CurrentCompartment,'_LatencyIndex'],[CurrentCompartment,'_Latency']};
        for f_ind = 1:length(RelevantFields)
            CurrentField = RelevantFields{f_ind};
            BigMat          = Data.(CurrentStrain).Stats.(ModelName).(CurrentField).Matrix;
            NumberOfRepeats = sum(~isnan(BigMat),3);
            Data.(CurrentStrain).Stats.(ModelName).(CurrentField).Mean   = nanmean(BigMat,3);
            Data.(CurrentStrain).Stats.(ModelName).(CurrentField).Median = nanmedian(BigMat,3);
            Data.(CurrentStrain).Stats.(ModelName).(CurrentField).STD    = nanstd(BigMat,[],3);
            Data.(CurrentStrain).Stats.(ModelName).(CurrentField).STD    = nanstd(BigMat,[],3) ./ sqrt(NumberOfRepeats);                       
        end          
    end    
end

if 0
%%%% Plot Model Stats %%%%
ModelName          = 'LNmodelWithRec_2';
CurrentCompartment = 'Soma';
b_vec_ticks        = 1:20:100;
b_vec              = Data.str2_CX17256_FromD6_Fast.Butanone_d6_To_Buffer.(['b_vec_for_',ModelName]);

for s_ind = 1:length(RelevantStrains)
    CurrentStrain  = RelevantStrains{s_ind};  
    XtickLabels    = Data.(CurrentStrain).PulsesInfo.FinalConcentration; 
    figure('name',[CurrentStrain,' , Responsiveness'],'position',[ 455   105   990   865]); 
        Matrix = Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Responsiveness']).Mean;
        imagesc(Matrix);
        set(gca,'ytick',b_vec_ticks,'yticklabel',b_vec(b_vec_ticks),'xtick',1:length(XtickLabels),'xticklabel',XtickLabels)
    figure('name',[CurrentStrain,' , Latency'],'position',[ 455   105   990   865]); 
        Matrix = Data.(CurrentStrain).Stats.(ModelName).([CurrentCompartment,'_Latency']).Mean;
        imagesc_withoutNaNs(1:size(Matrix,2), 1:size(Matrix,1), Matrix, [0 ceil(nanmax(Matrix(:)))], [], [1 1 1]); colorbar;
        set(gca,'ytick',b_vec_ticks,'yticklabel',b_vec(b_vec_ticks),'xtick',1:length(XtickLabels),'xticklabel',XtickLabels)
end
    
%%%%  Plot traces - for illustrator show only:  buffer and 2e-7  %%%%     
ModelName   = 'LNmodelWithRec_2';
b_vec_ticks = [1:20:81 99];
% b_vec_ticks = [1:20:41 51 61 71 99];
c_ind                   = 1;
b_vec                   = Data.str2_CX17256_FromD6_Fast.Butanone_d6_To_Buffer.(['b_vec_for_',ModelName]);
CurrentCompartment      = CompartmentFieldNames{c_ind};    
CurrentFieldName_Model  = [CurrentCompartment,'_PredictedDeltaFoverF_',ModelName];
CurrentFieldName_Data   = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed'];

for s_ind = 1:length(RelevantStrains)
    CurrentStrain  = RelevantStrains{s_ind};        
    PulsesNames    = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    YlabelNames    = num2cell(Data.(CurrentStrain).PulsesInfo.FinalConcentration);
    NumberOfPulses = length(PulsesNames);     
    if s_ind == 1
        TimeVec          = Data.(CurrentStrain).Butanone_d6_To_Buffer.TimeVec;
        TimeVecPredicted = TimeVec;
    else   % (s_ind==2)
        TimeVec          = Data.(CurrentStrain).Butanone_d6_To_Buffer.TimeVec;
        TimeVecPredicted = Data.(CurrentStrain).Butanone_d6_To_Buffer.Dye_TimeVec;
    end

    figure('name',[CurrentCompartment,' , ',CurrentStrain,'. Measured (k) vs. Model1 (r)']); 
    gcaHandles     = [];    
    CompartmentMax = 0;            

    for p_ind = 1:NumberOfPulses        
        CurrentPulse          = PulsesNames{p_ind};                    
        CurrentMeasuredMat    = Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Data);  % 2D
        CurrentPredictedMat3D = Data.(CurrentStrain).(CurrentPulse).(CurrentFieldName_Model); % 3D

        for b_ind = 1:length(b_vec_ticks)    
            b_ind_InVec = b_vec_ticks(b_ind);
            b           = b_vec(b_ind_InVec);
            subplot(NumberOfPulses,length(b_vec_ticks),(p_ind-1)*length(b_vec_ticks)+b_ind);
            CurrentPredictedMat = squeeze(CurrentPredictedMat3D(b_ind_InVec,:,:));             % 2D

            plot(TimeVec         , CurrentMeasuredMat',           '-','color',[0.7 0.7 0.7]); hold on;
            plot(TimeVecPredicted, CurrentPredictedMat',          '-','color',[1 0.7 0.7]);   hold on;
            plot(TimeVec         , nanmean(CurrentMeasuredMat,1), '-','color','k'); hold on;
            plot(TimeVecPredicted, nanmean(CurrentPredictedMat,1),'-','color','r'); hold on;
            xlim([-2 min([TimeVec(end) 70])]);             
            if p_ind == 1
                title(b);
            end                
            if b_ind == 1
                ylabel(YlabelNames{p_ind});
            end
            if p_ind~=NumberOfPulses
                set(gca,'xticklabel',[]);
            end
            if b_ind~=1
                set(gca,'yticklabel',[]);
            end
            gcaHandles     = [gcaHandles gca];
            CompartmentMax = max([CompartmentMax CurrentMeasuredMat(:)' CurrentPredictedMat(:)']);                               
        end                       
    end
    YLIMITS = [-0.2 CompartmentMax*1.1];
    set(gcaHandles,'ylim',YLIMITS)           
end  

end

return

% Useful function for combining outputs of bootstrap functions after parallel computing
function Combine_BootstrapThresholdForStepResponse

load('H:\BootstrapThresholdForStepResponse_01.mat','Data','ModelParams')
Data_Combined = Data;
FileNames = {'H:\BootstrapThresholdForStepResponse_02.mat','H:\BootstrapThresholdForStepResponse_03.mat',...
             'H:\BootstrapThresholdForStepResponse_04.mat','H:\BootstrapThresholdForStepResponse_05.mat'};
StrainNames           = fieldnames(Data);
RelevantStrainIndices = setdiff(1:length(StrainNames), find(strcmpi(StrainNames,'str2_CX17256_ToBuffer_Slow')));
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
for f=1:length(FileNames)
    FileNames{f}
    load(FileNames{f},'Data');
    CurrentComparment  = 'Soma';
    CurrentStrainName  = 'str2_CX17256_FromD7_Fast';   
    StartIndex         = length(Data_Combined.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_Threshold']).Vector) + 1; 
    EndIndex           = length(Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_Threshold']).Vector) + StartIndex - 1; 
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrainIndex = RelevantStrainIndices(s_ind);
        CurrentStrainName  = StrainNames{CurrentStrainIndex};    

        for c_ind = 1:length(CompartmentFieldNames)        
            CurrentComparment     = CompartmentFieldNames{c_ind};
            Data_Combined.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_InitialThresholdGuess']).Vector(StartIndex:EndIndex)  = ...
                     Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_InitialThresholdGuess']).Vector ;
            Data_Combined.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_Threshold']).Vector(StartIndex:EndIndex)              = ...
                     Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_Threshold']).Vector             ;
            Data_Combined.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_n']).Vector(StartIndex:EndIndex)                      = ....
                     Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_n']).Vector                     ;
        end
    end          
end    
Data = Data_Combined;

for s_ind = 1:length(RelevantStrainIndices)
    CurrentStrainIndex = RelevantStrainIndices(s_ind);
    CurrentStrainName  = StrainNames{CurrentStrainIndex};    
    for c_ind = 1:length(CompartmentFieldNames)        
        CurrentComparment     = CompartmentFieldNames{c_ind};
        RepeatsPerCondition   = Data.(CurrentStrainName).Stats.([CurrentComparment,'_Information']).NumberOfRepeats;               
        MeanNumOfRepeats      = mean(RepeatsPerCondition);
        
        CurrentStructure         = Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_Threshold']);
        Vector                   = CurrentStructure.Vector;        
        CurrentStructure.Mean    = mean(Vector);
        CurrentStructure.Median  = median(Vector);
        CurrentStructure.STD     = std(Vector,1);
        CurrentStructure.SEM     = std(Vector,1)./ sqrt(MeanNumOfRepeats);
        CurrentStructure.MeanNumOfRepeats = MeanNumOfRepeats;
        Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_Threshold']) = CurrentStructure;

        CurrentStructure         = Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_n']);
        Vector                   = CurrentStructure.Vector;        
        CurrentStructure.Mean    = mean(Vector);
        CurrentStructure.Median  = median(Vector);
        CurrentStructure.STD     = std(Vector,1);
        CurrentStructure.SEM     = std(Vector,1)./ sqrt(MeanNumOfRepeats);
        CurrentStructure.MeanNumOfRepeats = MeanNumOfRepeats;                
        Data.(CurrentStrainName).Stats.Bootstrap.([CurrentComparment,'_n']) = CurrentStructure;
    end
end         
save('H:\BootstrapThresholdForStepResponse_Combined.mat','Data','ModelParams','-v7.3')

% load('H:\BootstrapThresholdForStepResponse_Combined.mat','Data')

return


