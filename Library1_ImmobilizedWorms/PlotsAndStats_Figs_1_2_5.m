function PlotsAndStats_Figs_1_2_5
%
%  Written by Sagi Levy, September 2019
%
% This function generates plots and statistics of data from AWC(ON) calcium imaging experiments in immobilized animals   
% Related to figures 1, 2 and 5 and their respective supplementary figures  

%% load processed data (see details in the function 'ProcessDataForPaper'
% load('H:\Data_AWConImaging_Processed.mat',...
load('H:\Data_AWConImaging_Processed.mat',...
             'Data','ModelParams',...
             'ACT_CX17256','ACT_CX17255','ACT_Egl4','ACT_Egl4_ad450','ACT_VaryInitialC',...
             'ModelPredictions_VaryFinalC',  'MeasuredStats_VaryFinalC',...
             'ModelPredictions_VaryInitialC','MeasuredStats_VaryInitialC');  % load processed data

%% load and plot threshold constants (K)
Pvalues = PlotThresholdConstants (Data, ModelParams);                              % Related to Figure 2 and S3

%% load and plot adaptation times (tau)
Pvalues = PlotAdaptationTimes(ACT_CX17256, ACT_CX17255, ACT_Egl4, ACT_Egl4_ad450); % Related to Figure 2, Figure S3 and Figure 5E

%% Compare model performances. Related to Figure 2 and S3
% Here ACT predictions are based on very small fraction of fitting data
CompartmentFieldNames = {'Soma'};
FitField              = 'FitBasedOnOnePulseResponsiveness';               % Fit using only data of responsiveness (w/o latency) that comes from a single experiment type (slow decrease from 11 to 2 microM)
CompareModelPerformances(ModelPredictions_VaryFinalC,   ACT_CX17256,      FitField, CompartmentFieldNames);       
CompareModelPerformances(ModelPredictions_VaryInitialC, ACT_VaryInitialC, FitField, CompartmentFieldNames, true);  

% Here ACT predictions are based on all data
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
FitField              = 'FitBasedOnResponsivenessAndLatency_UseDataSTD'; % Best fit using all data
CompareModelPerformances(ModelPredictions_VaryFinalC,   ACT_CX17256,      FitField, CompartmentFieldNames);       
CompareModelPerformances(ModelPredictions_VaryInitialC, ACT_VaryInitialC, FitField, CompartmentFieldNames, true); 

%% Mismatch between predicted and measured Latency in slow odor decrease experiments for non-ACT models (Figure 1 and S2)
PlotLatencyPredictions_NotACT (ModelPredictions_VaryFinalC, MeasuredStats_VaryFinalC)

%% LN model plots (Figure S2)
RelevantStrainNames = {'str2_CX17256_FromD6_Fast'};
PlotImpulseFunctions(Data, RelevantStrainNames);

% RelevantStrainIndices = [1:4 9 10];                                           % All str2_CX17256 experiments
RelevantStrainIndices   = [2 9];                                                % Fast vs. slow decrease of odor from 11 microM ('str2_CX17256_FromD6_Fast','str2_CX17256_FromD6_Slow') 
PlotPredictedResponses_LNmodel_Traces(Data, RelevantStrainIndices);             % used for Figure S2D
RelevantStrainIndices   = 9;                                                    % Figure S2E. Slow decrease of odor from 11 microM ('str2_CX17256_FromD6_Slow') 
PlotPredictedResponses_LNmodel_Stats(Data,  RelevantStrainIndices, true, false) % all
RelevantStrainIndices   = 1:4;                                                  % Figure S2F. Fast decrease of odor from various initial conditions
PlotPredictedResponses_LNmodel_Stats(Data,  RelevantStrainIndices, false, true) % all

PlotOptimizeRectifierForLNmodel(Data); % Figures S2H and S2I

% The parameter used for GCaMP5a saturation level does not modify our results: LN model cannot predict slow changing stimuli 
RelevantStrainIndices   = [2 9];                                                % Fast vs. slow decrease of odor from 11 microM ('str2_CX17256_FromD6_Fast','str2_CX17256_FromD6_Slow') 
PlotPredictedResponses_LNmodel_Traces_SaturationScreen(Data, RelevantStrainIndices);             % used for Figure S2D
RelevantStrainIndices   = 9;                                                    % Figure S2E. Slow decrease of odor from 11 microM ('str2_CX17256_FromD6_Slow') 
PlotPredictedResponses_LNmodel_Stats_SaturationScreen(Data,  RelevantStrainIndices, true) % all

% Deconvolving GCaMP5a kernel does not modify our results: LN model cannot predict slow changing stimuli 
RelevantStrainIndices   = [2 9];                                                % Fast vs. slow decrease of odor from 11 microM ('str2_CX17256_FromD6_Fast','str2_CX17256_FromD6_Slow') 
PlotPredictedResponses_LNmodel_Traces_GCaMP_Deconv(Data, RelevantStrainIndices);             % used for Figure S2D
RelevantStrainIndices   = 9;                                                    % Figure S2E. Slow decrease of odor from 11 microM ('str2_CX17256_FromD6_Slow') 
PlotPredictedResponses_LNmodel_Stats_GCaMP_Deconv(Data,  RelevantStrainIndices, true, false) % all
RelevantStrainIndices   = 1:4;                                                  % Figure S2F. Fast decrease of odor from various initial conditions
PlotPredictedResponses_LNmodel_Stats_GCaMP_Deconv(Data,  RelevantStrainIndices, false, true) % all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%% Plot raw data and basic data features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%% Plot raw sensory activity data. Figures 1C and 1D
CompartmentsToPlot    = {'Soma','Dendrite','Axon'};  % subplots
Concentrations        = [0    0.1116    0.5580    1.1160    2.2320    3.3480    4.4640    8.9280   11.1600]; % microM
ColorsPerConcetration = {[1 0 0],[0.7 0 0.3],[0.3 0 0.7],[0 0 1],[0 0.5 0.5],[0 1 0],[0.5 0.5 0],[0.7 0.7 0.2],[0 0 0]}; 

 % raw data, (all) single repeats
PlotType = 0;  YLIMITS= [-0.1 6];                            
XLIMITS= [-2 30];  StrainName = 'str2_CX17256_FromD6_Fast';
PlotResponseToStimulus(Data, StrainName, CompartmentsToPlot, Concentrations, ColorsPerConcetration, YLIMITS, XLIMITS, PlotType);
XLIMITS= [-2 70];  StrainName = 'str2_CX17256_FromD6_Slow';
PlotResponseToStimulus(Data, StrainName, CompartmentsToPlot, Concentrations, ColorsPerConcetration, YLIMITS, XLIMITS, PlotType);
% average over repeats and SEM
PlotType = 1;  YLIMITS= [-0.1 3.5];   
XLIMITS= [-2 30];  StrainName = 'str2_CX17256_FromD6_Fast';
PlotResponseToStimulus(Data, StrainName, CompartmentsToPlot, Concentrations, ColorsPerConcetration, YLIMITS, XLIMITS, PlotType);
XLIMITS= [-2 70];  StrainName = 'str2_CX17256_FromD6_Slow';
PlotResponseToStimulus(Data, StrainName, CompartmentsToPlot, Concentrations, ColorsPerConcetration, YLIMITS, XLIMITS, PlotType);

%% Plot dye and dye derivatives dynamics 
RelevantStrainNames = {'str2_CX17256_FromD6_Slow','str2_CX17256_ToBuffer_Slow'};
PlotDyeStats(Data,RelevantStrainNames)     % Plot Dye statistics, including calculation of smoothed C, dC/dt, d2C/dt2, DeltaC/C. Used for Figure 1F 

%% Plot WT responses to fast versus slow decrease in odor concentration, related to Figure 1
RelevantStrainIndices = [2 9];  % RelevantStrainNames = {'str2_CX17256_FromD6_Slow','str2_CX17256_ToBuffer_Slow'};
PlotInfo.Xlimits = [0.4 20];
Data = PlotFastVsSlowResponsiveness(Data, RelevantStrainIndices, PlotInfo); % Figure 1E
PlotInfo.Xlimits             = [0.04 20];
PlotInfo.Ylimits_DeltaFoverF = [-0.1 4];
Data = PlotFastVsSlowResponsiveness(Data, RelevantStrainIndices, PlotInfo); % Figure S2A

%% egl-4 
Pvalue_CompareStrains = PlotEgl4Stats(Data);                        % Statistics: Figure 5C, S5C, S5D
Plot_egl4_Traces(Data)                                              % Traces: Figure 5A, S5A, S5B
[Pvalue, OddsRatio, FDR] = PlotFastVsSlowResponsiveness_egl4(Data); % Responsiveness: Figure 5B 
Pvalues = PlotLatency_egl4(Data);                                   % Latency: Figure 5D, S5E 

%%  Plot Dye variability (Figure S1C)
PlotDyeVariability(Data)

%% Schematics showing models similar responses to step responses, and how to challenge models' predictions, related to Figure 1A 
SchematicsForModelingComparison_Fig1A

%% Number Of Repeats
[NumberOfRepeats, CellToPrint] = FindNumberOfRepeats(Data); 
MaxNumberOfRepeats             = max(NumberOfRepeats(:));
% xlswrite('D:\Paper\Repeats_Test.xlsx',CellToPrint);     % Generate Excel file with number of repeats    

% set(get(0,'children'),'WindowStyle','dock')

return

%% Schematics, Figure 1A
function  SchematicsForModelingComparison_Fig1A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%   Figure 1A for ON neuron (activated at HIGH concentrations)   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% First row of 1A
DeltaT = 0.001;
time   = -2:DeltaT:5;
ZeroIndex = find(time==0);
tao  = 0.05;
C1   = ones(1,length(time));
C2   = ones(1,length(time));
C1(ZeroIndex:end) = 2-(2-1)*exp(-time(ZeroIndex:end)/tao);
C2(ZeroIndex:end) = 3-(3-1)*exp(-time(ZeroIndex:end)/tao);
Derivative1       = [diff(C1) NaN]/DeltaT;
Derivative2       = [diff(C2) NaN]/DeltaT;

Ylimits.C         = [0.85 3*1.1];
Ylimits.dCdt      = [-3 45];
Ylimits.dCdtoverC = [-3 45];
XlimitsZoom.C         = [0.4/3*time(1)/time(end)  0.4];
XlimitsZoom.dCdt      = [0.1/3*time(1)/time(end)  0.1];
XlimitsZoom.dCdtoverC = [0.1/3*time(1)/time(end)  0.1];

figure('position',[680   812   560   166]); 
subplot(1,3,1); plot(time,C2,'r',time,C1,'k'); ylim([Ylimits.C]); xlim(XlimitsZoom.C); title('C')
subplot(1,3,2); plot(time,Derivative2,'r',time,Derivative1,'k'); ylim([Ylimits.dCdt]); xlim(XlimitsZoom.dCdt); title('dC/dt')
subplot(1,3,3); plot(time,Derivative2./C2,'r',time,Derivative1./C1,'k'); ylim([Ylimits.dCdtoverC]); xlim(XlimitsZoom.dCdtoverC); title('(dC/dt)/C');

%% Second row of 1A
tao  = 0.05;
C1                = ones(1,length(time));
C1(ZeroIndex:end) = 2-(2-1)*exp(-time(ZeroIndex:end)/tao);

C2   = ones(1,length(time));
Cf  = 2;
C0  = 1;
tao = 0.8;
t   = time(ZeroIndex:end);
C2(ZeroIndex:end) = C0 * (Cf/C0).^( 1 - exp(-t/tao).*(1 + t/tao) );

Derivative1       = [diff(C1) NaN]/DeltaT;
Derivative2       = [diff(C2) NaN]/DeltaT;

Ylimits.C         = [0.9 2.4];
Ylimits.dCdt      = [-0.05 0.6];
Ylimits.dCdtoverC = [-0.04 0.5];
Xlimits.C         = [-0.5 4];
Xlimits.dCdt      = [-0.5 3];
Xlimits.dCdtoverC = [-0.5 3];

figure('position',[680   812   560   166]); 
subplot(1,3,1); plot(time,C2,'r',time,C1,'k'); ylim([Ylimits.C]); xlim(Xlimits.C); title('C')
subplot(1,3,2); plot(time,Derivative2,'r',time,Derivative1,'k'); ylim([Ylimits.dCdt]); xlim(Xlimits.dCdt); title('dC/dt')
subplot(1,3,3); plot(time,Derivative2./C2,'r',time,Derivative1./C1,'k');  ylim([Ylimits.dCdtoverC]); xlim(Xlimits.dCdtoverC); title('(dC/dt)/C');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%   Figure 1A for OFF neuron (activated at LOW concentrations)   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% First row of 1A
DeltaT = 0.001;
time   = -2:DeltaT:5;
ZeroIndex = find(time==0);
tao  = 0.05;
C1   = 3*ones(1,length(time));
C2   = 3*ones(1,length(time));
C1(ZeroIndex:end) = 2+1*exp(-time(ZeroIndex:end)/tao);
C2(ZeroIndex:end) = 1+2*exp(-time(ZeroIndex:end)/tao);
Derivative1       = [diff(C1) NaN]/DeltaT;
Derivative2       = [diff(C2) NaN]/DeltaT;

Ylimits.C         = [0.8 3*1.07];
Ylimits.dCdt      = [-44 4.5];
Ylimits.dCdtoverC = [-15 1.5];
XlimitsZoom.C         = [0.4/3*time(1)/time(end)  0.4];
XlimitsZoom.dCdt      = [0.1/3*time(1)/time(end)  0.1];
XlimitsZoom.dCdtoverC = [0.1/3*time(1)/time(end)  0.1];

figure('position',[680   812   560   166]); 
subplot(1,3,1); plot(time,C2,'r',time,C1,'k'); ylim([Ylimits.C]); xlim(XlimitsZoom.C); title('C')
subplot(1,3,2); plot(time,Derivative2,'r',time,Derivative1,'k'); ylim([Ylimits.dCdt]); xlim(XlimitsZoom.dCdt); title('dC/dt')
subplot(1,3,3); plot(time,Derivative2./C2,'r',time,Derivative1./C1,'k'); ylim([Ylimits.dCdtoverC]); xlim(XlimitsZoom.dCdtoverC); title('(dC/dt)/C');

%% Second row of 1A
tao  = 0.05;
C1   = 3*ones(1,length(time));
C1(ZeroIndex:end) = 1+2*exp(-time(ZeroIndex:end)/tao);

C2  = 3*ones(1,length(time));
Cf  = 1;
C0  = 3;
tao = 0.8;
t   = time(ZeroIndex:end);
C2(ZeroIndex:end) = C0 * (Cf/C0).^( 1 - exp(-t/tao).*(1 + t/tao) );

Derivative1       = [diff(C1) NaN]/DeltaT;
Derivative2       = [diff(C2) NaN]/DeltaT;

Ylimits.C         = [0.4 3.2];
Ylimits.dCdt      = [-1.6 0.14];
Ylimits.dCdtoverC = [-0.65 0.07];
Xlimits.C         = [-0.5 4];
Xlimits.dCdt      = [-0.5 3];
Xlimits.dCdtoverC = [-0.5 3];

figure('position',[680   812   560   166]); 
subplot(1,3,1); plot(time,C2,'r',time,C1,'k'); ylim([Ylimits.C]); xlim(Xlimits.C); title('C')
subplot(1,3,2); plot(time,Derivative2,'r',time,Derivative1,'k'); ylim([Ylimits.dCdt]); xlim(Xlimits.dCdt); title('dC/dt')
subplot(1,3,3); plot(time,Derivative2./C2,'r',time,Derivative1./C1,'k');  ylim([Ylimits.dCdtoverC]); xlim(Xlimits.dCdtoverC); title('(dC/dt)/C');

return

%% Plots of raw data, and basic data features

function PlotResponseToStimulus(Data, StrainName, CompartmentsToPlot, Concentrations, ColorsPerConcetration, YLIMITS, XLIMITS, PlotType)

RelevantDataFieldName = '_DeltaFOverF';

PulseTypeNames     = Data.(StrainName).PulsesInfo.PulseFromButanone; 
FinalConcentration = Data.(StrainName).PulsesInfo.FinalConcentration; 
NumOfPulses        = length(PulseTypeNames);
TimeVec            = Data.(StrainName).(PulseTypeNames{1}).TimeVec;

% FindColors
ColorPerPulse = cell(1,NumOfPulses);
LEGENDS       = cell(1,NumOfPulses);
for p_ind = 1:NumOfPulses
    [~, index]           = min(abs(FinalConcentration(p_ind) - Concentrations));
    ColorPerPulse{p_ind} = ColorsPerConcetration{index}; 
    LEGENDS{p_ind}       = num2str(FinalConcentration(p_ind));
end

%% plot

if PlotType==1
    figure('name',StrainName,'position',[644    64   520   934]);  
    for c_ind = 1:length(CompartmentsToPlot)
        CurrentCompartment = CompartmentsToPlot{c_ind};            
        subplot(length(CompartmentsToPlot),1,c_ind)
        for p_ind = 1:NumOfPulses
            CurrentPulse = PulseTypeNames{p_ind};                                  
            MAT   = Data.(StrainName).(CurrentPulse).([CurrentCompartment,RelevantDataFieldName]);                   
            MEAN  = nanmean(MAT,1);  
            SEM   = nanstd(MAT,1,1)/sqrt(size(MAT,1));    
            upper = MEAN + SEM;
            lower = MEAN - SEM;     
            fill([TimeVec, TimeVec(end:-1:1)],[upper, lower(end:-1:1)],[0.85 0.85 0.85 ],'edgecolor','none','FaceAlpha',1); hold on;              
        end
        for p_ind = 1:NumOfPulses
            CurrentPulse = PulseTypeNames{p_ind};                                  
            MAT   = Data.(StrainName).(CurrentPulse).([CurrentCompartment,RelevantDataFieldName]);                   
            MEAN  = nanmean(MAT,1);  
            plot(TimeVec,MEAN,'-','linewidth',0.75,'color',ColorPerPulse{p_ind});  hold on;                       
        end        
        ylim(YLIMITS); xlim(XLIMITS); 
        ylabel([CurrentCompartment,' \DeltaF/F']);
    end
    
elseif PlotType==0
    for c_ind = 1:length(CompartmentsToPlot)
        CurrentCompartment = CompartmentsToPlot{c_ind};   
        figure('name',[StrainName,' , ',CurrentCompartment],'position',[644    64   520   934]);  
        for p_ind = 1:NumOfPulses
            CurrentPulse = PulseTypeNames{p_ind};                                  
            MAT          = Data.(StrainName).(CurrentPulse).([CurrentCompartment,RelevantDataFieldName]);                   
            subplot(NumOfPulses,1,p_ind)
            plot(TimeVec,MAT,'-','linewidth',0.3,'color',ColorPerPulse{p_ind});          
            ylim(YLIMITS); xlim(XLIMITS);   
            legend(LEGENDS{p_ind});
            ylabel(['\DeltaF/F']);
        end
    end 
    
end

return

function PlotFastVsSlowResponsiveness(Data, RelevantStrainIndices, PlotInfo)

StrainNames           = fieldnames(Data);
ColorsPerStrain       = {'k',[0.5 0.5 0.5]};
MarkerSizePerStrain   = [6 7];
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);         
FigureString          = [StrainNames{RelevantStrainIndices(1)},' (k) vs. ', StrainNames{RelevantStrainIndices(2)},' (grey)'];
C_curve               = 10.^(-1:0.005:1.5); % in microM

% Responsiveness
figure('name',['Responsiveness,  ',FigureString],'position',[652   296  252  545]); 
for c_ind = 1:NumberOfCompartment
    CurrentCompartment      = CompartmentFieldNames{c_ind};
    subplot(NumberOfCompartment,1,c_ind)
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
        CurrentStats  = Data.(CurrentStrain).Stats;
        
        T             = CurrentStats.([CurrentCompartment,'_Threshold']) ;
        n             = CurrentStats.([CurrentCompartment,'_n']);
        R_curve       = 100*T^n ./ (T^n + C_curve.^n);
        plot(C_curve, R_curve,'color',ColorsPerStrain{s_ind}); hold on;
        
        C             = CurrentStats.([CurrentCompartment,'_Information']).FinalConcentrations;
        errorbar(C, CurrentStats.([CurrentCompartment,'_Responsiveness']).Mean, CurrentStats.([CurrentCompartment,'_Responsiveness']).SEM, 'o',...
                 'color',ColorsPerStrain{s_ind},'markerfacecolor',ColorsPerStrain{s_ind},'markersize',MarkerSizePerStrain(s_ind)); hold on;
    end
    set(gca,'xscale','log','ylim',[-5 105],'box','off') 
    if exist('PlotInfo','var') && isfield(PlotInfo,'Xlimits')
        xlim(PlotInfo.Xlimits)
    end
end    

% DeltaF/F
figure('name',['DeltaF/F,  ',FigureString],'position',[652   296  252  545]); 
for c_ind = 1:NumberOfCompartment
    CurrentCompartment      = CompartmentFieldNames{c_ind};
    subplot(NumberOfCompartment,1,c_ind)
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
        CurrentStats  = Data.(CurrentStrain).Stats;
                
        C             = CurrentStats.([CurrentCompartment,'_Information']).FinalConcentrations;
        C_ForDisplay  = C;
        C_ForDisplay(C_ForDisplay==0) = 0.05;  % FOR LOG SCALE DISPLAY INCLUDING ZERO. DON'T FORGET TO ADD AN AXIS BREAK: '//' 
        
        errorbar(C_ForDisplay, CurrentStats.([CurrentCompartment,'_DeltaFOverF_MAX']).Mean, CurrentStats.([CurrentCompartment,'_DeltaFOverF_MAX']).SEM, 'o',...
                 'color',ColorsPerStrain{s_ind},'markerfacecolor',ColorsPerStrain{s_ind},'markersize',MarkerSizePerStrain(s_ind)); hold on;
    end
    set(gca,'xscale','log','ylim',[-.05 3],'box','off') 
    if exist('PlotInfo','var') && isfield(PlotInfo,'Xlimits')
        xlim(PlotInfo.Xlimits)
    end
    if exist('PlotInfo','var') && isfield(PlotInfo,'Ylimits_DeltaFoverF')
        ylim(PlotInfo.Ylimits_DeltaFoverF)
    end
end    

return

function PlotDyeStats(Data, RelevantStrainNames)
%% Input arguments examples:
% RelevantStrainNames = {'str2_CX17256_FromD6_Slow','str2_CX17256_ToBuffer_Slow'};
warning('off')
MaxNumberOfRepeats = 40;  % free parameter. It just needs to be higher than the maximal number of repeats

for s_ind = 1:length(RelevantStrainNames)
    CurrentStrain  = RelevantStrainNames{s_ind};  
    PulsesNames    = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    NumberOfPulses = length(PulsesNames);    
    DyeInformation = Data.(CurrentStrain).DyeInformation;

    for p_ind = 1:NumberOfPulses
        CurrentPulseName = PulsesNames{p_ind}; 
        DyeTimeVec    = Data.(CurrentStrain).(CurrentPulseName).Dye_TimeVec;    
        
        FigureString = [CurrentStrain,', ',CurrentPulseName];
        row          = 0;
        
        figure('name',FigureString);    
        for rep_ind = 1:MaxNumberOfRepeats
            if all(isnan(DyeInformation.C(p_ind,rep_ind,:)),3)
                % no more traces for this
                break
            end
            row = row + 1;
            if row > 8 % up to 8 rows per figure
                row = 1;
                figure('name',FigureString);    
            end       
                        
            % Plots
            subplot(8,4,(row-1)*4+1); 
            plot(DyeTimeVec,squeeze(DyeInformation.C_Smoothed(p_ind,rep_ind,:)),'r-'); hold on;
%             plot(DyeTimeVec,squeeze(DyeInformation.C(p_ind,rep_ind,:)),'k:'); 
            if row==1, title('C'); end
            ylabel(rep_ind)
            xlim([-10 90])
            
            subplot(8,4,(row-1)*4+2); 
            plot(DyeTimeVec,squeeze(DyeInformation.dCdt_Smoothed(p_ind,rep_ind,:)),'r-'); hold on;
%             plot(DyeTimeVec,squeeze(DyeInformation.dCdt(p_ind,rep_ind,:)),'k:'); 
            if row==1, title('dCdt'); end
            xlim([-10 90])
            
            subplot(8,4,(row-1)*4+3); 
            plot(DyeTimeVec,squeeze(DyeInformation.d2Cdt2(p_ind,rep_ind,:)),'r-'); hold on;
            if row==1, title('d^2C/dt^2'); end
            xlim([-10 90])
            
            subplot(8,4,(row-1)*4+4); 
            plot(DyeTimeVec,squeeze(DyeInformation.DeltaCoverC(p_ind,rep_ind,:)),'r-'); hold on;
            if row==1, title('\DeltaC/C'); end                                                
            xlim([-10 90])
        end
    end   
end
warning('on')

return

function PlotFastVsSlowMaxDeltaFoverF(Data) 

RelevantStrainNames   = {'str2_CX17256_FromD6_Fast','str2_CX17256_FromD6_Slow'};
ColorPerStrain        = {'k',[0.6 0.6 0.6]};
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);   
MaxNumberOfRepeats    = 40;  % free parameter. It just needs to be higher than the maximal number of repeats
  
% Max DelatFoverF vs. concentration for both measured and predicted   
figure('name','Maximal responses to slow (grey) and fast (black) decrease in odor concentration','position', [384   218   367   580]);       
for s_ind = 1:length(RelevantStrainNames)
    CurrentStrain       = RelevantStrainNames{s_ind};
    CurrentColor        = ColorPerStrain{s_ind};
    
    PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;  % units of dilution
    NumberOfPulses      = length(PulsesNames);         

    for c_ind = 1:NumberOfCompartment            
        CurrentCompartment     = CompartmentFieldNames{c_ind};
        MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
%       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
        MaxStructure.Measured.Matrix  = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;            
        MaxStructure.Measured.Mean    = zeros(NumberOfPulses,1,'single')*NaN;            
        MaxStructure.Measured.STD     = zeros(NumberOfPulses,1,'single')*NaN;            
        MaxStructure.Measured.SEM     = zeros(NumberOfPulses,1,'single')*NaN;            
        
        for p_ind = 1:NumberOfPulses
            CurrentPulse  = PulsesNames{p_ind};
            if ~isfield(Data.(CurrentStrain),CurrentPulse)
                continue
            end
            MeasuredMat            = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field);            
            CurrentMaxMeasuredVec  = nanmax(MeasuredMat,[],2);   % Result: max of different repeats of the same pulse            
            MaxStructure.Measured.Matrix(p_ind,1:length(CurrentMaxMeasuredVec))  = CurrentMaxMeasuredVec;
        end 
        Mat                         = MaxStructure.Measured.Matrix;
        MaxStructure.Measured.Mean  = nanmean(Mat,2);    
        MaxStructure.Measured.STD   = nanstd(Mat,[],2);    
        MaxStructure.Measured.SEM   = nanstd(Mat,[],2) ./ sqrt(sum(~isnan(Mat),2));    
        
        subplot(NumberOfCompartment,1,c_ind)
        FinalConcentrationPerPulse_ForDisplay = FinalConcentrationPerPulse;
        FinalConcentrationPerPulse_ForDisplay(FinalConcentrationPerPulse_ForDisplay==0) = 0.01;    % FOR DISPLAY ON LOG SCALE, DON'T FORGET TO ADD A BREAK SIGN "//" on X AXIS !!!!
        
        % Max DelatFoverF vs. concentration for both measured and predicted   
        errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Measured.Mean,  MaxStructure.Measured.SEM, 'color',CurrentColor,'linestyle','none','marker','.');   hold on; 
        set(gca,'xscale','log');
    end        
end

return

function [NumberOfRepeats, CellToPrint] = FindNumberOfRepeats(Data) 

%% Extract Number Of Repeats
StrainNames           = fieldnames(Data);
MaxNumberOfPulseTypes = 20; 
NeuronFieldNames      = {'soma','dendrite','axon'};
NumberOfRepeats       = zeros(length(StrainNames), MaxNumberOfPulseTypes, length(NeuronFieldNames), 'single') * NaN; % Initialization
CellToPrint           = {'Strain Name','Pulse Name','n (Soma)','n (Axon)','n (Dendrite)'};
CellRow               = 1;

for s_count = 1:length(StrainNames)
    CurrentStrain     = StrainNames{s_count};
    CurrentStructure  = Data.(CurrentStrain);
    CurrentPulseNames = CurrentStructure.PulsesInfo.PulseFromButanone; % some of which may be empty ! 
    NumOfPulses       = length(CurrentPulseNames); % some of which may be empty ! 
    
    for p_ind = 1:NumOfPulses
        PulseName = CurrentPulseNames{p_ind};
        if isfield(CurrentStructure,PulseName)
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

%% egl-4
function Pvalue_CompareStrains = PlotEgl4Stats(Data)

%% Compute max(DeltaFoverF)
% egl-4 experiments were done in this order: 
%     Final concentrations of (buffer, 2d7, 4d7)
%     Final concentrations of (buffer, 6d7, 8d7)
%     Final concentrations of (buffer, d6) OR (d6), for fast and slow tansitions respectively. 

Fields_Strings        = {'str2_CX17255_All','str2_Egl4_All','str2_Egl4_ad450_All'}; % WT, egl-4(lf), egl-4(gf)
NewFields             = {'CX17255','egl4','egl4_ad450'}; 
PulsesNames           = {'Butanone_d6_To_Buffer'
                         'Butanone_d6_To_Butanone_2d7'
                         'Butanone_d6_To_Butanone_4d7'
                         'Butanone_d6_To_Butanone_6d7'
                         'Butanone_d6_To_Butanone_8d7'
                         'Butanone_d6_To_Butanone_d6'};
%     Data.str2_Egl4_AllSlow.PulsesInfo.PulseFromButanone; % similar pulses for all strains
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
clear MaxDeltaFOverF_All

for f_ind = 1:length(Fields_Strings)
    CurrentStrainField_Fast = [Fields_Strings{f_ind},'Fast'];
    CurrentStrainField_Slow = [Fields_Strings{f_ind},'Slow'];
    CurrentNewField         = NewFields{f_ind};
    MaxNumberOfRepeats      = length(Data.(CurrentStrainField_Fast).(PulsesNames{1}).Soma_DeltaFOverF_MAX); 
    
    for c_ind = 1:length(CompartmentFieldNames)            
        CurrentCompartment         = CompartmentFieldNames{c_ind};            
        CurrentMaxActivityFieldRaw = [CurrentCompartment,'_DeltaFOverF'];                        
        MaxMatrix.Fast             = ones(MaxNumberOfRepeats,length(PulsesNames),'single')*NaN; 
        MaxMatrix.Slow             = ones(MaxNumberOfRepeats,length(PulsesNames),'single')*NaN;         

        running_ind = 1;
        for p_ind = 1:length(PulsesNames)
            CurrentPulse = PulsesNames{p_ind};               
            MaxVecFast   = nanmax(Data.(CurrentStrainField_Fast).(CurrentPulse).(CurrentMaxActivityFieldRaw),[],2);
            MaxVecSlow   = nanmax(Data.(CurrentStrainField_Slow).(CurrentPulse).(CurrentMaxActivityFieldRaw),[],2);
            
            MaxMatrix.Fast(running_ind:running_ind+length(MaxVecFast)-1,p_ind) = MaxVecFast;
            MaxMatrix.Slow(running_ind:running_ind+length(MaxVecSlow)-1,p_ind) = MaxVecSlow;  
            if p_ind==3 || p_ind==5   % see notes above
                running_ind = find(~isnan(MaxMatrix.Fast(:,p_ind)),1,'last')+1;
            end
        end               
        % assign to structure
        MaxDeltaFOverF_All.(CurrentNewField).(CurrentCompartment) = MaxMatrix;                
    end        
end

%% Plot DeltaFoverF(fast) / DeltaFoverF(slow)
Pvalue_CompareStrains = DeltaFOverF_Stats(MaxDeltaFOverF_All);

return

function Pvalue_CompareStrains = DeltaFOverF_Stats(MaxDeltaFOverF)
NeuronFieldsForMaxDeltaFOverF = {'Soma','Dendrite','Axon'};
StrainFieldsForMaxDeltaFOverF = {'egl4','CX17255','egl4_ad450'};
% ColorPerStrain              = {'r','k','b'};
NumberOfPulses                = size(MaxDeltaFOverF.CX17255.Soma.Fast,2);


%% log(ratio) stats
for n_ind = 1:length(NeuronFieldsForMaxDeltaFOverF)
    CurrentFieldForMaxDeltaFOverF = NeuronFieldsForMaxDeltaFOverF{n_ind};
        
    Fast       = MaxDeltaFOverF.CX17255.(CurrentFieldForMaxDeltaFOverF).Fast;
    Slow       = MaxDeltaFOverF.CX17255.(CurrentFieldForMaxDeltaFOverF).Slow;
    RatioRef   = Fast ./ Slow;     
    Fast       = MaxDeltaFOverF.egl4.(CurrentFieldForMaxDeltaFOverF).Fast;
    Slow       = MaxDeltaFOverF.egl4.(CurrentFieldForMaxDeltaFOverF).Slow;        
    Ratio_egl4 = Fast ./ Slow;
    Fast       = MaxDeltaFOverF.egl4_ad450.(CurrentFieldForMaxDeltaFOverF).Fast;
    Slow       = MaxDeltaFOverF.egl4_ad450.(CurrentFieldForMaxDeltaFOverF).Slow;        
    Ratio_egl4_ad450 = Fast ./ Slow;
     
    % CX17255 versus egl-4 (n478)
    MW_pValue     = ones(1,NumberOfPulses);
    MW_EffectSize = zeros(1,NumberOfPulses);
    for p_ind = 1:NumberOfPulses
        Vec_WT   = RatioRef(:,p_ind);   Vec_WT   = Vec_WT(~isnan(Vec_WT));
        Vec_egl4 = Ratio_egl4(:,p_ind); Vec_egl4 = Vec_egl4(~isnan(Vec_egl4));
        if ~(isempty(Vec_WT)||isempty(Vec_egl4))            
            [MW_pValue(p_ind), ~, ~, MW_EffectSize(p_ind)] = ranksum_inline(Vec_WT, Vec_egl4,'tail','right');
        end
    end
    egl4_SmallDiff.Pvalue.(CurrentFieldForMaxDeltaFOverF)     = MW_pValue;
    egl4_SmallDiff.EffectSize.(CurrentFieldForMaxDeltaFOverF) = MW_EffectSize;

    % CX17255 versus egl-4 (ad450)
    MW_pValue     = ones(1,NumberOfPulses);
    MW_EffectSize = zeros(1,NumberOfPulses);
    for p_ind = 1:NumberOfPulses
        Vec_WT   = RatioRef(:,p_ind);         Vec_WT   = Vec_WT(~isnan(Vec_WT));
        Vec_egl4 = Ratio_egl4_ad450(:,p_ind); Vec_egl4 = Vec_egl4(~isnan(Vec_egl4));
        if ~(isempty(Vec_WT)||isempty(Vec_egl4))            
            [MW_pValue(p_ind), ~, ~, MW_EffectSize(p_ind)] = ranksum_inline(Vec_WT, Vec_egl4,'tail','left');
        end
    end
    egl4_ad450_LargeDiff.Pvalue.(CurrentFieldForMaxDeltaFOverF)     = MW_pValue;
    egl4_ad450_LargeDiff.EffectSize.(CurrentFieldForMaxDeltaFOverF) = MW_EffectSize;
    
end


%% (Fast/Slow ratio), Figures 5C and S5C 
StrainFieldsForMaxDeltaFOverF_BarsDisplay = StrainFieldsForMaxDeltaFOverF;   % egl4, WT, egl4_ad450 
BaseLines = (0:0.2:1)*11.16; % Scaling from dilution to microM
clear NumOfRepeats
for n_ind = 1:length(NeuronFieldsForMaxDeltaFOverF)
    CurrentFieldForMaxDeltaFOverF = NeuronFieldsForMaxDeltaFOverF{n_ind};
    figure('name',[CurrentFieldForMaxDeltaFOverF,' (dFast/dSlow) ratio bar plot. dF=max(deltaF/F)'],'position',[ 593   655   180   306])
    
    BarsMeans = zeros(NumberOfPulses, length(StrainFieldsForMaxDeltaFOverF_BarsDisplay));
    BarsSEMs  = zeros(NumberOfPulses, length(StrainFieldsForMaxDeltaFOverF_BarsDisplay));    
        
    for s_ind = 1:length(StrainFieldsForMaxDeltaFOverF_BarsDisplay)            
        CurrentStrainForMaxDeltaFOverF = StrainFieldsForMaxDeltaFOverF_BarsDisplay{s_ind}; 
        Fast  = MaxDeltaFOverF.(CurrentStrainForMaxDeltaFOverF).(CurrentFieldForMaxDeltaFOverF).Fast;
        Slow  = MaxDeltaFOverF.(CurrentStrainForMaxDeltaFOverF).(CurrentFieldForMaxDeltaFOverF).Slow;             
        for p_ind = 1:NumberOfPulses
            Vec   = Fast(:,p_ind)./Slow(:,p_ind);                       
            % ErrorBar 
            BarsSEMs(p_ind,s_ind)    = nanstd(Vec,[],1)./sqrt(length(find(~isnan(Vec))));
            BarsMeans(p_ind,s_ind)   = nanmean(Vec);
            NumOfRepeats.(CurrentFieldForMaxDeltaFOverF)(p_ind,s_ind)= length(find(~isnan(Vec)));
        end
    end   
%     barwitherr_inline(BarsMeans', BarsSEMs' , Properties); ylim([0 34])      
    errorbar(BaseLines, BarsMeans(:,2), BarsSEMs(:,2),'k.','linewidth',1,'markersize',15,'capsize',10); hold on;
    errorbar(BaseLines, BarsMeans(:,1), BarsSEMs(:,1),'r.','linewidth',1,'markersize',15,'capsize',10); hold on;
    errorbar(BaseLines, BarsMeans(:,3), BarsSEMs(:,3),'b.','linewidth',1,'markersize',15,'capsize',10); hold on;
    set(gca,'xtick',0:5:10,'xticklabel',0:5:10,'ytick',0:5:35,'yticklabel',0:5:30)
    ylim([0 34]); xlim([-1 12]) 
end


%% Djs and significance testing by bootstrapping
%%%%% The same with dot plot    %%%%%
clear NumOfRepeats
bootstrap_iterations = 1e4; 
NumberOfStrains      = length(StrainFieldsForMaxDeltaFOverF_BarsDisplay);
ReferenceStrain      = 2; % WT
tic 
for n_ind = 1:length(NeuronFieldsForMaxDeltaFOverF)
    CurrentFieldForMaxDeltaFOverF = NeuronFieldsForMaxDeltaFOverF{n_ind};
    
    BarsMeans    = zeros(NumberOfPulses, NumberOfStrains);
    BarsSEMs     = zeros(NumberOfPulses, NumberOfStrains);    
    BarsMeans_bs = zeros(bootstrap_iterations, NumberOfPulses, NumberOfStrains);
       
    for s_ind = 1:NumberOfStrains
        CurrentStrainForMaxDeltaFOverF = StrainFieldsForMaxDeltaFOverF_BarsDisplay{s_ind}; 
        Fast  = MaxDeltaFOverF.(CurrentStrainForMaxDeltaFOverF).(CurrentFieldForMaxDeltaFOverF).Fast;
        Slow  = MaxDeltaFOverF.(CurrentStrainForMaxDeltaFOverF).(CurrentFieldForMaxDeltaFOverF).Slow; 
        FastDivideBySlow = Fast ./ Slow;        

        for p_ind = 1:NumberOfPulses
            Values      = FastDivideBySlow(:,p_ind);                       
            Values      = Values(~isnan(Values));    
            CurrentN    = length(Values);
            RandomDraws = datasample(Values,bootstrap_iterations*CurrentN);      % vector
            RandomDraws = reshape(RandomDraws,[CurrentN, bootstrap_iterations]);  % Matrix
            
            BarsMeans(p_ind,s_ind)        = mean(Values);
            BarsSEMs(p_ind,s_ind)         = std(Values,[],1)/sqrt(CurrentN);
            BarsMeans_bs(:,p_ind,s_ind)   = mean(RandomDraws,1);
            NumOfRepeats.(CurrentFieldForMaxDeltaFOverF)(p_ind,s_ind)= CurrentN; 
        end               
    end   
    AverageRepeatsPerPulse = mean(NumOfRepeats.(CurrentFieldForMaxDeltaFOverF),1);
    
    Reference = BarsMeans(:,ReferenceStrain);
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_AllData = ones(1,NumberOfStrains)*NaN;
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs      = ones(bootstrap_iterations, NumberOfStrains,'single')*NaN;
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs_mean = ones(1,NumberOfStrains)*NaN;
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs_std  = ones(1,NumberOfStrains)*NaN;
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs_sem  = ones(1,NumberOfStrains)*NaN;
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_Pvalue_UsingRefbs = ones(1,NumberOfStrains)*NaN;
    
    for s_ind = 1:NumberOfStrains
        CurrentMeans   = BarsMeans(:,s_ind);
        Djs            = JehnsenShannonDivergence_inline (CurrentMeans,Reference);         
        Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_AllData (s_ind) = Djs;        
            
        for it = 1:bootstrap_iterations
            CurrentMeans   = squeeze(BarsMeans_bs(it,:,s_ind))';
            Djs            = JehnsenShannonDivergence_inline (CurrentMeans,Reference)';         
            Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs(it,s_ind)   = Djs;            
        end        
    end
    
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs_mean = mean(Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs,1);
    Current_STD_vec                                           = std(Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs,[],1);
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs_std  = Current_STD_vec;
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs_sem  = Current_STD_vec ./ AverageRepeatsPerPulse; % COMPUTED FROM AVERAGE NUMBER IN A PULSE (more strict), NOT BY NUMBER OF REPEATS IN ALL PULSES
                
    % What are the chances of getting the current Djs given the reference distribution        
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_Pvalue  = ones(1,NumberOfStrains)*NaN;            
    Djs_structure.(CurrentFieldForMaxDeltaFOverF).EffectSize  = ones(1,NumberOfStrains)*NaN;            
    RefVec = Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs(:,ReferenceStrain);  
    for s_ind = 1:NumberOfStrains
        CurrentVec = Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_bs(:,s_ind);  
        [P, ~, ~, MW_EffectSize] = ranksum_inline(RefVec, CurrentVec);               
        if P < 1/bootstrap_iterations
            P = 1/bootstrap_iterations; % Lower limit of Pvalue depends on number of iterations
        end       
        Djs_structure.(CurrentFieldForMaxDeltaFOverF).Djs_Pvalue(s_ind)  = P;            
        Djs_structure.(CurrentFieldForMaxDeltaFOverF).EffectSize(s_ind)  = MW_EffectSize;            
    end
    
end
toc
%% Figure S5D
Colors                 = [1 0 0; 0 0 1; 0 0 0]; % {'egl-4 (n478)','egl-4 (ad450)','Control'}
RepeatCutOffForDisplay = 50;
DataArray              = Djs_structure.Soma.Djs_bs(1:RepeatCutOffForDisplay,[1 3 2]); % egl-4 (LF), egl-4 (GF) and then WT
CurrentMeans           = Djs_structure.Soma.Djs_bs_mean([1 3 2]); % egl-4 (LF), egl-4 (GF) and then WT
figure('position',[680   558   315   420])
UnivarScatter_SL(DataArray,'Label',{'egl-4 (n478)','egl-4 (ad450)','Control'},'MarkerFaceColor',Colors,'MarkerEdgeColor',Colors,'PointSize',20,'PlotOnlyDots',true);
hold on; plot(get(gca,'xlim'),CurrentMeans(end)*ones(1,2),'k:')
ylabel('Distance from control (Djs)')
set(gca,'PlotBoxAspectRatioMode','auto','YLimMode','manual','XLimMode','manual');


%% Compute P-values on Djs
% P-values by comparison of bootstrap distributions 
ReferenceIndex = 2;    % WT
Vec_Reference  = Djs_structure.Soma.Djs_bs(:,ReferenceIndex);
Pvalue_LargerThanReference  = zeros(1,NumberOfStrains)*NaN;
for s_ind = 1:NumberOfStrains 
    X = Vec_Reference;
    Y = Djs_structure.Soma.Djs_bs(:,s_ind);
    [~, Pvalue_Y_LargerThan_X]         = CompareBootstrapVectors(X, Y);
    Pvalue_LargerThanReference(s_ind)  = Pvalue_Y_LargerThan_X;
end
Cutoff_Pvalue =  1/length(Vec_Reference);                                                % Lower Limit due to finite sampling
Pvalue_LargerThanReference(Pvalue_LargerThanReference<Cutoff_Pvalue) = Cutoff_Pvalue;  % Lower Limit due to finite sampling

Pvalue_CompareStrains.ReferenceIndex              = ReferenceIndex;
Pvalue_CompareStrains.ReferenceName               = StrainFieldsForMaxDeltaFOverF_BarsDisplay{ReferenceIndex};
Pvalue_CompareStrains.AllNames                    = StrainFieldsForMaxDeltaFOverF_BarsDisplay;
Pvalue_CompareStrains.Pvalue_LargerThanReference  = Pvalue_LargerThanReference;

% FDR correction for multiple fits
NonReferenceIndex                                 = setdiff(1:NumberOfStrains, ReferenceIndex);
Pvalue_CompareStrains.FDR_LargerThanReference     = ComputeFDR_MultipleSignificanceLevels(Pvalue_LargerThanReference, NonReferenceIndex);

return

function MultiPlots(Data, PulsesNames, PulsesTitles, CurrentStrainNames, Colors, CurrentNeuronField, CurrentTime, CurrentXlim, CurrentYlim, CurrentYlabel, YTicks, DisplayWithFilter, AveragingWindowFrame)

if exist('DisplayWithFilter','var') && DisplayWithFilter &&  exist('AveragingWindowFrame','var') 
    a = 1;
    b = ones(1,AveragingWindowFrame)/AveragingWindowFrame;
    D = round(mean(grpdelay(b,a)));  
else
    DisplayWithFilter = false;
end    

for p_ind = 1:length(PulsesNames)
    CurrentPulse = PulsesNames{p_ind};
    CurrentTitle = PulsesTitles{p_ind};
    for s_ind = 1:length(CurrentStrainNames)
        CurrentStrain = CurrentStrainNames{s_ind};
        CurrentColor  = Colors{s_ind};
        sp_ind        = (s_ind-1)*length(PulsesNames) + p_ind;
        subplot(length(CurrentStrainNames),length(PulsesNames),sp_ind);
        MAT    = Data.(CurrentStrain).(CurrentPulse).(CurrentNeuronField);
        
        if DisplayWithFilter
            VectorLength = size(MAT,2);
            MAT = filter(b,a,MAT,[],2); 
            MAT = MAT(:,D+1:end);
            MAT(:,(VectorLength-D):VectorLength) = NaN;                    
        end        
        plot(CurrentTime, MAT','color',CurrentColor);      hold on; 
        xlim(CurrentXlim); ylim(CurrentYlim); set(gca,'ytick',YTicks)
        if p_ind==1
            ylabel(CurrentYlabel);
        end
        if s_ind==1
            title(CurrentTitle);
        end
    end
end

return

function Plot_egl4_Traces(Data)

PulsesNames     = {'Butanone_d6_To_Buffer','Butanone_d6_To_Butanone_2d7','Butanone_d6_To_Butanone_4d7','Butanone_d6_To_Butanone_6d7','Butanone_d6_To_Butanone_8d7','Butanone_d6_To_Butanone_d6'};
PulsesTitles    = {'0','2','4','6','9','11'}; % microM
Colors          = {'r','k','b'};             % egl-4 (lof), WT, egl-4 (gof)
% Colors          = {'r',[0.6 0.6 0.6],'b'};     % egl-4 (lof), WT, egl-4 (gof)
StrainNamesFast = {'str2_Egl4_AllFast','str2_CX17255_AllFast','str2_Egl4_ad450_AllFast'};
StrainNamesSlow = {'str2_Egl4_AllSlow','str2_CX17255_AllSlow','str2_Egl4_ad450_AllSlow'};
TimeFast        = -2:0.1:30;
TimeSlow        = -2:0.1:90;
Xlim_Fast       = [-2 30];
Xlim_Slow       = [-2 70];
% DisplayWithFilter    = false;
DisplayWithFilter    = true;
AveragingWindowFrame = 5;

%% Soma Raw deltaF/F %%%%%%
CurrentNeuronField = 'Soma_DeltaFOverF';
CurrentYlabel      = '\DeltaF / F_0';
CurrentYlim        = [-0.2 6];
YTicks             = 0:6;
% Step response 
figure('position',get(0,'ScreenSize'),'name','Soma step response. Raw \DeltaF / F_0'); 
MultiPlots(Data, PulsesNames, PulsesTitles, StrainNamesFast, Colors, CurrentNeuronField, TimeFast, Xlim_Fast, CurrentYlim, CurrentYlabel, YTicks, DisplayWithFilter, AveragingWindowFrame)
% Response to graded signal
figure('position',get(0,'ScreenSize'),'name','Soma response to graded signal. Raw \DeltaF / F_0'); 
MultiPlots(Data, PulsesNames, PulsesTitles, StrainNamesSlow, Colors, CurrentNeuronField, TimeSlow, Xlim_Slow, CurrentYlim, CurrentYlabel, YTicks, DisplayWithFilter, AveragingWindowFrame)

return

function [Pvalue_structure, OddsRatio_structure, FDR_structure] = PlotFastVsSlowResponsiveness_egl4(Data)
% Figure 5B

StrainNamesFast       = {'str2_Egl4_AllFast','str2_CX17255_AllFast','str2_Egl4_ad450_AllFast'};
StrainNamesSlow       = {'str2_Egl4_AllSlow','str2_CX17255_AllSlow','str2_Egl4_ad450_AllSlow'};
NewFields             = {'egl4','CX17255','egl4_ad450'}; 
ColorsPerStrain       = {'r','k','b'};
MarkerSize            = 7;
% CompartmentFieldNames = {'Soma','Dendrite','Axon'}; % Use this too plot responses of other compartments 
CompartmentFieldNames = {'Soma'};
NumberOfCompartment   = length(CompartmentFieldNames);         

clear StructForTable Pvalue_structure  OddsRatio_structure  Chi2_structure;

%% Responsiveness plots
for c_ind = 1:NumberOfCompartment
    CurrentCompartment      = CompartmentFieldNames{c_ind};
    figure('name',[CurrentCompartment, ' responsiveness'],'position',[652   296  252  545]); 
    for s_ind = 1:length(StrainNamesFast)
        subplot(length(StrainNamesFast),1,s_ind)
        % Fast
        CurrentStrain = StrainNamesFast{s_ind};        
        CurrentStats  = Data.(CurrentStrain).Stats;                
        C             = CurrentStats.([CurrentCompartment,'_Information']).FinalConcentrations;
        errorbar(C, CurrentStats.([CurrentCompartment,'_Responsiveness']).Mean, CurrentStats.([CurrentCompartment,'_Responsiveness']).SEM, 'o',...
                 'color',ColorsPerStrain{s_ind},'markerfacecolor',ColorsPerStrain{s_ind},'markersize',MarkerSize,'linestyle',':'); hold on;
        % Slow
        CurrentStrain = StrainNamesSlow{s_ind};        
        CurrentStats  = Data.(CurrentStrain).Stats;                
        C             = CurrentStats.([CurrentCompartment,'_Information']).FinalConcentrations;
        errorbar(C, CurrentStats.([CurrentCompartment,'_Responsiveness']).Mean, CurrentStats.([CurrentCompartment,'_Responsiveness']).SEM, 'o',...
                 'color',ColorsPerStrain{s_ind},'markerfacecolor',ColorsPerStrain{s_ind},'markersize',MarkerSize,'linestyle','-'); hold on;
        
        xlabel('Final Concentration');
        ylabel([NewFields{s_ind},' \DeltaF/F'])
        set(gca,'ylim',[-5 105],'xlim',[-1 12],'box','off') 
    end
end    

%% P-values: Compare responsiveness to fast and slow odor concentration decrease
%  Compute Pvalues for the hypothesis that: * slower signals result in less responsiveness *  
for c_ind = 1:NumberOfCompartment
    CurrentCompartment      = CompartmentFieldNames{c_ind};
    for s_ind = 1:length(StrainNamesFast)
        CurrentStrainNewField = NewFields{s_ind};        
        % Fast
        CurrentStrain = StrainNamesFast{s_ind};        
        CurrentStats  = Data.(CurrentStrain).Stats;                
        FastMat = CurrentStats.([CurrentCompartment,'_Responsiveness']).Matrix;        
        % Slow
        CurrentStrain = StrainNamesSlow{s_ind};        
        CurrentStats  = Data.(CurrentStrain).Stats;                
        SlowMat = CurrentStats.([CurrentCompartment,'_Responsiveness']).Matrix;               
        % structure for fisher test tables    
        StructForTable.Response_Fast   = nansum(FastMat==100,2);
        StructForTable.NoResponse_Fast = nansum(FastMat==0,2);
        StructForTable.Response_Slow   = nansum(SlowMat==100,2);
        StructForTable.NoResponse_Slow = nansum(SlowMat==0,2);
        % compute P-values
        NumOfPulses       = length(StructForTable.Response_Fast);
        CurrentPvalues    = zeros(1,NumOfPulses);
        CurrentOddsRatio  = zeros(1,NumOfPulses);
%         CurrentEffectSize = zeros(1,NumOfPulses);
        for pulse_ind = 1:NumOfPulses   
            Column1      = [StructForTable.Response_Fast(pulse_ind)  ; StructForTable.Response_Slow(pulse_ind)];
            Column2      = [StructForTable.NoResponse_Fast(pulse_ind); StructForTable.NoResponse_Slow(pulse_ind)];
            CurrentTable = table(Column1,Column2,'VariableNames',{'Response','NoResponse'},'RowNames',{'Fast','Slow'});
            [h,p,stats]  = fishertest(CurrentTable,'tail','right');
            CurrentPvalues(pulse_ind)   = p;    
            CurrentOddsRatio(pulse_ind) = stats.OddsRatio;    
%             % For Chi2
%             Chi2Mat      = [Column1 Column2];
%             TotalRepeats = sum(Chi2Mat(:));
%             Chi2Mat(Chi2Mat==0)=eps;
%             Expected     = [Chi2Mat(1,:);Chi2Mat(1,:)]; 
%             CurrentChi2  = sum((Chi2Mat(:)-Expected(:)).^2 ./ Expected(:));
%             CurrentEffectSize(pulse_ind) = sqrt(CurrentChi2/TotalRepeats);
        end
        Pvalue_structure.(CurrentStrainNewField).(CurrentCompartment)    = CurrentPvalues;
        OddsRatio_structure.(CurrentStrainNewField).(CurrentCompartment) = CurrentOddsRatio;
%         Chi2_structure.(CurrentStrainNewField).(CurrentCompartment)      = CurrentEffectSize;                  
    end
end    

%% FDR correction for multiple hypotheses
FDRfields = {'FDR005','FDR0005','FDR0001','FDR00005'};
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind}; 
    for s_ind = 1:length(StrainNamesFast)             
        CurrentStrainNewField = NewFields{s_ind};    
        for fdr_ind = 1:length(FDRfields)
            CurrentFDRfield = FDRfields{fdr_ind};
            FDR_structure.(CurrentStrainNewField).(CurrentCompartment).(CurrentFDRfield) = ones(1,NumOfPulses)*NaN;
        end
    end
    PvaluesMat         = [Pvalue_structure.egl4.(CurrentCompartment)', Pvalue_structure.CX17255.(CurrentCompartment)', Pvalue_structure.egl4_ad450.(CurrentCompartment)'];
    for p_ind = 1:NumOfPulses  
        CurrentPvalues = PvaluesMat(p_ind,:);
        CurrentFDR = ComputeFDR_MultipleSignificanceLevels(CurrentPvalues,1:length(CurrentPvalues));
        for fdr_ind = 1:length(FDRfields)
            CurrentFDRfield = FDRfields{fdr_ind};
            FDR_structure.egl4.(CurrentCompartment).(CurrentFDRfield)(p_ind)       = CurrentFDR.(CurrentFDRfield)(1);
            FDR_structure.CX17255.(CurrentCompartment).(CurrentFDRfield)(p_ind)    = CurrentFDR.(CurrentFDRfield)(2);
            FDR_structure.egl4_ad450.(CurrentCompartment).(CurrentFDRfield)(p_ind) = CurrentFDR.(CurrentFDRfield)(3);
        end
    end
end
         
return

function Pvalues = PlotLatency_egl4(Data)
% Figure 5B

StrainNamesSlow       = {'str2_Egl4_AllSlow','str2_CX17255_AllSlow','str2_Egl4_ad450_AllSlow'};
NewFields             = {'egl4','CX17255','egl4_ad450'}; 
% CompartmentFieldNames = {'Soma','Dendrite','Axon'}; % Use this too plot responses of other compartments 
CompartmentFieldNames = {'Soma'};
NumberOfCompartment   = length(CompartmentFieldNames);         

%% Latency plot
Properties.BarProps.EdgeColor      = {'r','k','b'};  % egl-4(n478), WT, egl-4(ad450). Change 'barwitherr_inline' function inputs accordingly  
Properties.BarProps.FaceColor      = {'r','k','b'};
Properties.ErrorBarProps.linewidth = {1,1,1,1,1,1,1};
Properties.ErrorBarProps.color     = {'k','k','k','k','k','k'};
Properties.gcaProps.xticklabel     = {'0','2','4','6','9','11'};   % corresponding to ticks of 1:size(Y,1)   
% Properties.gcaProps.xlim           = [0.5 6.5];
Properties.gcaProps.xlim           = [0.5 2.5];

clear LatencyCell
for c_ind = 1:NumberOfCompartment
    CurrentCompartment      = CompartmentFieldNames{c_ind};
    figure('name',[CurrentCompartment, ' Latency'],'position',[652   296  252  545]); 
    for s_ind = 1:length(StrainNamesSlow)
        CurrentStrain = StrainNamesSlow{s_ind};        
        CurrentStats  = Data.(CurrentStrain).Stats;    
        LatencyCell.Matrix{s_ind} = CurrentStats.([CurrentCompartment,'_Latency_InSec']).Matrix;
        LatencyCell.Mean{s_ind}   = CurrentStats.([CurrentCompartment,'_Latency_InSec']).Mean;
        LatencyCell.SEM{s_ind}    = CurrentStats.([CurrentCompartment,'_Latency_InSec']).SEM;                
    end 
    barwitherr_inline([LatencyCell.Mean{1}'  LatencyCell.Mean{2}'  LatencyCell.Mean{3}']', ...
                  [LatencyCell.SEM{1}'   LatencyCell.SEM{2}'   LatencyCell.SEM{3}']', Properties); ylim([0 27])   % egl-4(n478), WT, egl-4(ad450)

    xlabel('Final Concentration');
    ylabel([NewFields{s_ind},' latency [sec]'])
    set(gca,'box','off')     
end

%% P-values and FDR correction: Compare Latency of different genotypes 
clear Pvalues
ReferenceIndex = 2; % WT
for c_ind = 1:NumberOfCompartment
    CurrentCompartment      = CompartmentFieldNames{c_ind};
    for s_ind = 1:length(StrainNamesSlow)
        CurrentStrain = StrainNamesSlow{s_ind};        
        CurrentStats  = Data.(CurrentStrain).Stats;    
        LatencyCell.Matrix{s_ind} = CurrentStats.([CurrentCompartment,'_Latency_InSec']).Matrix;
    end 
    NumOfPulses      = size(LatencyCell.Matrix{1},1);    
    LatencyMatrix_WT = LatencyCell.Matrix{ReferenceIndex};     
    
    %% Wann Whitney U test: compare mutant to WT latency
    MW_Pvalue_SmallerThanWT = ones(1,NumOfPulses)*NaN;
    MW_Pvalue_LargerThanWT  = ones(1,NumOfPulses)*NaN;
    MW_EffectSize           = zeros(1,NumOfPulses)*NaN;     
    for s_ind = 1:length(StrainNamesSlow)
        CurrentStrain = NewFields{s_ind};     
        for p_ind = 1:NumOfPulses
            Vec_WT     = LatencyMatrix_WT(p_ind,:);          Vec_WT     = Vec_WT(~isnan(Vec_WT));
            Vec_mutant = LatencyCell.Matrix{s_ind}(p_ind,:); Vec_mutant = Vec_mutant(~isnan(Vec_mutant));
            if ~(isempty(Vec_WT)||isempty(Vec_mutant))            
                [MW_Pvalue_SmallerThanWT(p_ind), ~, ~, MW_EffectSize(p_ind)] = ranksum_inline(Vec_WT, Vec_mutant,'tail','right');
                 MW_Pvalue_LargerThanWT(p_ind)                               = ranksum_inline(Vec_WT, Vec_mutant,'tail','left');
            end
        end
        Pvalues.(CurrentStrain).(CurrentCompartment).Pvalue_SmallerThanWT = MW_Pvalue_SmallerThanWT;
        Pvalues.(CurrentStrain).(CurrentCompartment).Pvalue_LargerThanWT  = MW_Pvalue_LargerThanWT;
        Pvalues.(CurrentStrain).(CurrentCompartment).EffectSize           = MW_EffectSize;
    end    
end

FDRfields = {'FDR005','FDR0005','FDR0001','FDR00005'};
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind}; 
    for s_ind = 1:length(NewFields)             
        CurrentStrain = NewFields{s_ind};    
        for fdr_ind = 1:length(FDRfields)
            CurrentFDRfield = FDRfields{fdr_ind};
            Pvalues.(CurrentStrain).(CurrentCompartment).FDR_SmallerThanWT.(CurrentFDRfield) = zeros(1,NumOfPulses);
            Pvalues.(CurrentStrain).(CurrentCompartment).FDR_LargerThanWT.(CurrentFDRfield)  = zeros(1,NumOfPulses);
        end
    end
    % FDR correction for smaller Than WT
    PvaluesMat         = [Pvalues.egl4.(CurrentCompartment).Pvalue_SmallerThanWT', ...
                          Pvalues.CX17255.(CurrentCompartment).Pvalue_SmallerThanWT', ...
                          Pvalues.egl4_ad450.(CurrentCompartment).Pvalue_SmallerThanWT'];                      
    for p_ind = 1:NumOfPulses  
        CurrentPvalues  = PvaluesMat(p_ind,:);
        RelevantIndices = setdiff(1:length(CurrentPvalues),ReferenceIndex);
        CurrentFDR      = ComputeFDR_MultipleSignificanceLevels(CurrentPvalues,RelevantIndices);
        for fdr_ind = 1:length(FDRfields)
            CurrentFDRfield = FDRfields{fdr_ind};
            Pvalues.egl4.(CurrentCompartment).FDR_SmallerThanWT.(CurrentFDRfield)(p_ind)       = CurrentFDR.(CurrentFDRfield)(1);
            Pvalues.egl4_ad450.(CurrentCompartment).FDR_SmallerThanWT.(CurrentFDRfield)(p_ind) = CurrentFDR.(CurrentFDRfield)(2);
        end
    end
    % FDR correction for larger Than WT
    PvaluesMat         = [Pvalues.egl4.(CurrentCompartment).Pvalue_LargerThanWT', ...
                          Pvalues.CX17255.(CurrentCompartment).Pvalue_LargerThanWT', ...
                          Pvalues.egl4_ad450.(CurrentCompartment).Pvalue_LargerThanWT'];                      
    for p_ind = 1:NumOfPulses  
        CurrentPvalues  = PvaluesMat(p_ind,:);
        RelevantIndices = setdiff(1:length(CurrentPvalues),ReferenceIndex);
        CurrentFDR      = ComputeFDR_MultipleSignificanceLevels(CurrentPvalues,RelevantIndices);
        for fdr_ind = 1:length(FDRfields)
            CurrentFDRfield = FDRfields{fdr_ind};
            Pvalues.egl4.(CurrentCompartment).FDR_LargerThanWT.(CurrentFDRfield)(p_ind)       = CurrentFDR.(CurrentFDRfield)(1);
            Pvalues.egl4_ad450.(CurrentCompartment).FDR_LargerThanWT.(CurrentFDRfield)(p_ind) = CurrentFDR.(CurrentFDRfield)(2);
        end
    end    
end
       
%% Wiskers (median and percentiles)
Colors  = [1 0 0;       0 0 0;          0 0 1];       % {'egl-4 (n478)','WT','egl-4 (ad450)'}
clear LatencyCell
for c_ind = 1:NumberOfCompartment
    CurrentCompartment      = CompartmentFieldNames{c_ind};
    for s_ind = 1:length(StrainNamesSlow)
        CurrentStrain = StrainNamesSlow{s_ind};        
        CurrentStats  = Data.(CurrentStrain).Stats;    
        LatencyCell.(CurrentCompartment).Matrix{s_ind} = CurrentStats.([CurrentCompartment,'_Latency_InSec']).Matrix;
    end 
end
for c_ind = 1:NumberOfCompartment
    CurrentCompartment             = CompartmentFieldNames{c_ind};  
    figure('name',[CurrentCompartment, ' latency']);     
    
    p_ind = 1; subplot(1,2,p_ind); 
    DataArray       = [LatencyCell.(CurrentCompartment).Matrix{1}(p_ind,:)', ...
                       LatencyCell.(CurrentCompartment).Matrix{2}(p_ind,:)', ...
                       LatencyCell.(CurrentCompartment).Matrix{3}(p_ind,:)']; % egl-4 (LF), WT, egl-4 (GF) 
    h_b = boxplot(DataArray); 
    for i=1:3, set(h_b(:,i),'Color',Colors(i,:),'MarkerEdgeColor',Colors(i,:),'LineStyle','-'); end
    ylabel('Latency [sec]')
    ylim([-2 50])

    p_ind = 2; subplot(1,2,p_ind); 
    DataArray       = [LatencyCell.(CurrentCompartment).Matrix{1}(p_ind,:)',...
                       LatencyCell.(CurrentCompartment).Matrix{2}(p_ind,:)',...
                       LatencyCell.(CurrentCompartment).Matrix{3}(p_ind,:)']; % egl-4 (LF), WT, egl-4 (GF) 
    h_b = boxplot(DataArray); 
    for i=1:3, set(h_b(:,i),'Color',Colors(i,:),'MarkerEdgeColor',Colors(i,:),'LineStyle','-'); end
    ylabel('Latency [sec]')    
    ylim([-2 50])
end


return

%% Plots of models parameters and performances 
function Pvalues = PlotThresholdConstants (Data, ModelParams)

TrainingConditionNames = {'Naive','Desensitized'};
CompartmentFieldNames  = {'Soma','Dendrite','Axon'};
C = [1e-10:1e-10:9.9e-9  1e-8:1e-9:9.9e-8  1e-7:1e-8:9.9e-7  1e-6:1e-7:9.9e-6  1e-5:1e-6:9.9e-5 1e-4:1e-5:9.9e-4 1e-3:1e-4:5e-3]*1e6; % Concentration in microM

%% Saturation constant in different compartments and strains
NumOfCompartments = length(CompartmentFieldNames);
NumberOfIteration = length(ModelParams.Naive.Soma.Bootstrap.K.Vector);
DataArrayCombined = [];
Y_AllDataCombined = [];
for t_ind = 1:length(TrainingConditionNames) 
    CurrentTrainingCondition = TrainingConditionNames{t_ind};        
    figure('name',[CurrentTrainingCondition , 'Saturation Constant in different compartments'],'position',[680   558   243   420]);
    Y_AllData  = ones(1,NumOfCompartments,'single')*NaN;
    DataArray  = ones(NumberOfIteration,NumOfCompartments,'single')*NaN;

    for c_ind = 1:NumOfCompartments
        CurrentNeuronField = CompartmentFieldNames{c_ind}; 
        Y_AllData(c_ind)   = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).K;                 
        DataArray(:,c_ind) = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).Bootstrap.K.Vector;     
    end
    DataArrayCombined = [DataArrayCombined DataArray]; % Naive and then desensitized data in one matrix 
    Y_AllDataCombined = [Y_AllDataCombined Y_AllData]; % Naive and then desensitized data in one matrix 
    
    boxplot_SL(DataArray,'BoxColor',[0.7 0.7 0.7],'BoxWidthValue',0.7,'WhiskerDisplay','ON','WhiskerPercentiles',[5 95],...
                         'MedianColor','r','MedianWidth',4); hold on;
    plot(1:NumOfCompartments, Y_AllData, 'ko','markerfacecolor','k'); hold on
    xlim([0 NumOfCompartments+1]); %ylim([0 30])

    % P-values by comparison of bootstrap distributions 
    ReferenceIndex = 1;
    Vec_Reference  = DataArray(:,ReferenceIndex);
    Pvalue_SmallerThanReference = zeros(1,size(DataArray,2))*NaN;
    Pvalue_LargerThanReference  = zeros(1,size(DataArray,2))*NaN;
    for c_ind = 1:NumOfCompartments 
        X = Vec_Reference;
        Y = DataArray(:,c_ind);
        [Pvalue_X_LargerThan_Y, Pvalue_Y_LargerThan_X] = CompareBootstrapVectors(X, Y);
        Pvalue_SmallerThanReference(c_ind) = Pvalue_X_LargerThan_Y;
        Pvalue_LargerThanReference(c_ind)  = Pvalue_Y_LargerThan_X;
    end
    Cutoff_Pvalue =  1/size(DataArray,1);                                                    % Lower Limit due to finite sampling
    Pvalue_SmallerThanReference(Pvalue_SmallerThanReference<Cutoff_Pvalue) = Cutoff_Pvalue;  % Lower Limit due to finite sampling
    Pvalue_LargerThanReference(Pvalue_LargerThanReference<Cutoff_Pvalue)   = Cutoff_Pvalue;  % Lower Limit due to finite sampling

    Pvalue_CompareCompartments.(CurrentTrainingCondition).ReferenceIndex              = ReferenceIndex;
    Pvalue_CompareCompartments.(CurrentTrainingCondition).ReferenceName               = CompartmentFieldNames{ReferenceIndex};
    Pvalue_CompareCompartments.(CurrentTrainingCondition).AllNames                    = CompartmentFieldNames;
    Pvalue_CompareCompartments.(CurrentTrainingCondition).Pvalue_SmallerThanReference = Pvalue_SmallerThanReference;
    Pvalue_CompareCompartments.(CurrentTrainingCondition).Pvalue_LargerThanReference  = Pvalue_LargerThanReference;

    % FDR correction for multiple fits
    NonReferenceIndex                                                              = setdiff(1:size(DataArray,2), ReferenceIndex);
    Pvalue_CompareCompartments.(CurrentTrainingCondition).FDR_SmallerThanReference = ComputeFDR_MultipleSignificanceLevels(Pvalue_SmallerThanReference, NonReferenceIndex);
    Pvalue_CompareCompartments.(CurrentTrainingCondition).FDR_LargerThanReference  = ComputeFDR_MultipleSignificanceLevels(Pvalue_LargerThanReference, NonReferenceIndex);

end

% Compare desensitized to naive
figure('name','Saturation Constant in different compartments and strains','position',[668   511   352   420]);
boxplot_SL(DataArrayCombined,'BoxColor',[0.7 0.7 0.7],'BoxWidthValue',0.7,'WhiskerDisplay','ON','WhiskerPercentiles',[5 95],...
                             'MedianColor','r','MedianWidth',2); hold on;
plot(1:length(Y_AllDataCombined), Y_AllDataCombined, 'ko','markerfacecolor','k'); hold on
xlim([0 length(Y_AllDataCombined)+1]); ylim([2 5e3])
set(gca,'yscale','log')

% P-values by comparison of bootstrap distributions 
Pvalue_SmallerThanReference = zeros(1,size(DataArray,2))*NaN;
Pvalue_LargerThanReference  = zeros(1,size(DataArray,2))*NaN;
for c_ind = 1:NumOfCompartments 
    X = DataArrayCombined(:,c_ind);
    Y = DataArrayCombined(:,c_ind+NumOfCompartments);
    [Pvalue_X_LargerThan_Y, Pvalue_Y_LargerThan_X] = CompareBootstrapVectors(X, Y);
    Pvalue_SmallerThanReference(c_ind) = Pvalue_X_LargerThan_Y;
    Pvalue_LargerThanReference(c_ind)  = Pvalue_Y_LargerThan_X;
end
Cutoff_Pvalue =  1/size(DataArray,1);                                                    % Lower Limit due to finite sampling
Pvalue_SmallerThanReference(Pvalue_SmallerThanReference<Cutoff_Pvalue) = Cutoff_Pvalue;  % Lower Limit due to finite sampling
Pvalue_LargerThanReference(Pvalue_LargerThanReference<Cutoff_Pvalue)   = Cutoff_Pvalue;  % Lower Limit due to finite sampling

Pvalue_CompareDesensitizedToNaive.Pvalue_SmallerThanReference = Pvalue_SmallerThanReference;
Pvalue_CompareDesensitizedToNaive.Pvalue_LargerThanReference  = Pvalue_LargerThanReference;

% FDR correction for multiple fits
NonReferenceIndex                                          = 1:length(Pvalue_SmallerThanReference);
Pvalue_CompareDesensitizedToNaive.FDR_SmallerThanReference = ComputeFDR_MultipleSignificanceLevels(Pvalue_SmallerThanReference, NonReferenceIndex);
Pvalue_CompareDesensitizedToNaive.FDR_LargerThanReference  = ComputeFDR_MultipleSignificanceLevels(Pvalue_LargerThanReference, NonReferenceIndex);

%% Assign P-values to structure
Pvalues.CompareCompartments        = Pvalue_CompareCompartments;
Pvalues.CompareDesensitizedToNaive = Pvalue_CompareDesensitizedToNaive;

%% Compare threshold to model-based threshold Curves in black, dots in grey (Figure 2)
%  For supplementary: display all compartments, For figure 2: display only soma   
for t_ind = 1:length(TrainingConditionNames)    
    CurrentTrainingCondition = TrainingConditionNames{t_ind};        
    figure('name',[CurrentTrainingCondition , '.   dots= data, curve= model. NOTE!!! BREAK X AXIS!!! '],'position',get(0,'Screensize'));
    if t_ind==1
        PlotAxis          = [0.7e-3 1e4 -10 110];    
        ZeroForLogDisplay = 1e-3;
    elseif t_ind==2
        PlotAxis          = [0.7e-2 1e5 -10 110];    
        ZeroForLogDisplay = 1e-2;
    end    

    for n_ind = 1:length(CompartmentFieldNames)   % 3 columns in subplot for neuron areas
        CurrentNeuronField = CompartmentFieldNames{n_ind}; 
        K                       = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).K; 
        InitialConcetration_vec = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).InitialConcetrations; 
        n_vec                   = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).n; 
        Threshold_ByModel_vec   = K*(1-exp(-InitialConcetration_vec/K));
        ListOfStrains           = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).ListOfStrains; 

        for dil_ind = 1:length(InitialConcetration_vec)
            n                    = n_vec(dil_ind);
            Tmodel               = Threshold_ByModel_vec(dil_ind);
            Responsiveness_Model = (Tmodel^n)./(Tmodel^n+C.^n) *100;
            
            C_Data             = Data.(ListOfStrains{dil_ind}).Stats.([CurrentNeuronField,'_Information']).FinalConcentrations; 
            Responsiveness     = Data.(ListOfStrains{dil_ind}).Stats.([CurrentNeuronField,'_Responsiveness']).Mean;
            Responsiveness_SEM = Data.(ListOfStrains{dil_ind}).Stats.([CurrentNeuronField,'_Responsiveness']).SEM; 
            
            if isempty(C_Data)
                continue
            end
            FinalConc       = C_Data;
            FinalConc_Model = C;
            
            FinalConcForLogDisplay                            = FinalConc; 
            FinalConcForLogDisplay(FinalConcForLogDisplay==0) = ZeroForLogDisplay;    %  NOTE!!! BREAK X AXIS!!!

            sp_ind = n_ind + (dil_ind-1)*length(CompartmentFieldNames);
            subplot(length(InitialConcetration_vec),length(CompartmentFieldNames),sp_ind);

            plot(FinalConc_Model ,Responsiveness_Model,'-','color','k','linewidth',0.5); hold on;
            errorbar(FinalConcForLogDisplay,Responsiveness,Responsiveness_SEM,'.','color',[0.6 0.6 0.6],'linewidth',1, 'markerfacecolor',[0.6 0.6 0.6],'markersize',12,'capsize',0);
            set(gca,'xscale','log');
            hold on;

            axis(PlotAxis);
        end
    end 
    set(get(gcf,'Children'),'YTickMode','Manual','XTickMode','Manual','Ytick',0:50:100,'Xtick',10.^(-3:5),'box','off','XMinorTick','off');
    set(get(gcf,'Children'),'Xtick',10.^(-3:5),'Xticklabel',{'',-2,'',0,'',2,'',4,''});
    set(gcf,'position',[113   440   521   509]);
end

%% Plot actual vs. model-derived thresholds        
figure('name','Model vs. data, Thresholds','position', [695   719   868   255]);
ColorPerTrainingCondition     = {'k','r'};
MarkerPerTrainingCondition    = {'o','o'};
FillColorPerTrainingCondition = {[0.5 0.5 0.5],[1 0.5 0.5]};
for n_ind = 1:length(CompartmentFieldNames)   % 3 columns in subplot for neuron areas        
    CurrentNeuronField = CompartmentFieldNames{n_ind};   
    subplot(1,length(CompartmentFieldNames),n_ind);
    
    for t_ind = 1:length(TrainingConditionNames)    
        CurrentTrainingCondition = TrainingConditionNames{t_ind};             
        CurrentColor     = ColorPerTrainingCondition{t_ind};
        CurrentFillColor = FillColorPerTrainingCondition{t_ind};
        CurrentMarker    = MarkerPerTrainingCondition{t_ind};

        K                       = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).K; 
        InitialConcetration_vec = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).InitialConcetrations; 
        T_vec                   = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).Threshold; 
        K_std                   = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).Bootstrap.K.STD; 
        T_vec_std               = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).Bootstrap.Threshold.STD; 

        T_Error_High = T_vec_std;
        T_Error_Low  = T_vec_std;
        
        C_Curve           =  exp(log(InitialConcetration_vec(1)):0.01:log(InitialConcetration_vec(end)));
        Tmodel_curve      =  K *(1-exp(-C_Curve/K));
        C_Fill            =  exp(log(InitialConcetration_vec(1)):0.01:log(InitialConcetration_vec(end)));
        Tmodel_Fill       = (K+K_std) *(1-exp(-C_Fill/(K+K_std)));
        Tmodel_Fill2      = (K-K_std) *(1-exp(-C_Fill/(K-K_std)));
        C_Fill            = [C_Fill C_Fill(end:-1:1)];
        Tmodel_Fill       = [Tmodel_Fill Tmodel_Fill2(end:-1:1)];        
       
        fill(C_Fill,Tmodel_Fill,CurrentFillColor,'edgecolor','none'); hold on;
        plot(C_Curve,Tmodel_curve,'color',CurrentColor); hold on;
        errorbar(InitialConcetration_vec, T_vec,      T_Error_Low,      T_Error_High,     'linestyle','none','marker',CurrentMarker,'color',CurrentColor,'markerfacecolor',CurrentColor); hold on;
        set(gca,'xscale','log','yscale','log','ylim',[0.5 2e4],'xlim',[0.5 2e4],'xtick',10.^(-1:5),'ytick',10.^(-1:5));%,'YTickMode','Manual','XTickMode','Manual')
        if n_ind == 1, ylabel('Threshold'); end               
        xlabel('C_0')                    
    end 
end

%% For Supplementary Figures
% Hill coefficients         
figure('name','Model vs. data, Hill coefficients','position', [695   719   868   255]);
ColorPerTrainingCondition     = {'k','r'};
MarkerPerTrainingCondition    = {'^','o'};
for n_ind = 1:length(CompartmentFieldNames)   % 3 columns in subplot for neuron areas        
    CurrentNeuronField = CompartmentFieldNames{n_ind};   
    subplot(1,length(CompartmentFieldNames),n_ind);
    
    for t_ind = 1:length(TrainingConditionNames)    
        CurrentTrainingCondition = TrainingConditionNames{t_ind};             
        CurrentColor     = ColorPerTrainingCondition{t_ind};
        CurrentMarker    = MarkerPerTrainingCondition{t_ind};

        InitialConcetration_vec = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).InitialConcetrations; 
        n_vec                   = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).n; 
        n_vec_sem               = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).Bootstrap.n.SEM;         
   
        errorbar(InitialConcetration_vec, n_vec, n_vec_sem,'linestyle','none','marker',CurrentMarker,'color',CurrentColor,'markerfacecolor',CurrentColor); hold on;
        set(gca,'xscale','log','ylim',[0 35],'xlim',[0.5 2e4],'xtick',10.^(-1:5),'ytick',0:5:30);%,'YTickMode','Manual','XTickMode','Manual')
        if n_ind == 1, ylabel('Hill coefficient'); end               
        xlabel('C_0')                    
    end 
end

% Dot plot instead of bars for Thresholds and model-derived thresholds BASED ONLY ON LAST CONCETRATION MEASUREMENT         
figure('name','Model (based on T at saturating C) vs. data, Thresholds','position', [695   719   868   255]);
ColorPerTrainingCondition     = {'k','r'};
MarkerPerTrainingCondition    = {'o','o'};
FillColorPerTrainingCondition = {[0.5 0.5 0.5],[1 0.5 0.5]};
for n_ind = 1:length(CompartmentFieldNames)   % 3 columns in subplot for neuron areas        
    CurrentNeuronField = CompartmentFieldNames{n_ind};   
    subplot(1,length(CompartmentFieldNames),n_ind);
    
    for t_ind = 1:length(TrainingConditionNames)    
        CurrentTrainingCondition = TrainingConditionNames{t_ind};             
        CurrentColor     = ColorPerTrainingCondition{t_ind};
        CurrentFillColor = FillColorPerTrainingCondition{t_ind};
        CurrentMarker    = MarkerPerTrainingCondition{t_ind};

        InitialConcetration_vec = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).InitialConcetrations; 
        T_vec                   = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).Threshold; 
        T_vec_std               = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).Bootstrap.Threshold.STD; 
        
        % MODEL BASED ON THRESHOLD MEASUREMENT AT SATURATING INITIAL CONCENTRATION 
        K     = T_vec(end);
        K_std = T_vec_std(end);

        T_Error_High = T_vec_std;
        T_Error_Low  = T_vec_std;
        
        C_Curve           =  exp(log(InitialConcetration_vec(1)):0.01:log(InitialConcetration_vec(end)));
        Tmodel_curve      =  K *(1-exp(-C_Curve/K));
        C_Fill            =  exp(log(InitialConcetration_vec(1)):0.01:log(InitialConcetration_vec(end)));
        Tmodel_Fill       = (K+K_std) *(1-exp(-C_Fill/(K+K_std)));
        Tmodel_Fill2      = (K-K_std) *(1-exp(-C_Fill/(K-K_std)));
        C_Fill            = [C_Fill C_Fill(end:-1:1)];
        Tmodel_Fill       = [Tmodel_Fill Tmodel_Fill2(end:-1:1)];        
       
        fill(C_Fill,Tmodel_Fill,CurrentFillColor,'edgecolor','none'); hold on;
        plot(C_Curve,Tmodel_curve,'color',CurrentColor); hold on;
        errorbar(InitialConcetration_vec, T_vec,      T_Error_Low,      T_Error_High,     'linestyle','none','marker',CurrentMarker,'color',CurrentColor,'markerfacecolor',CurrentColor); hold on;
%         errorbar(InitialConcetration_vec, Tmodel_vec, Tmodel_Error_Low, Tmodel_Error_High,'linestyle','none','marker','x','color',CurrentColor); hold on;
        set(gca,'xscale','log','yscale','log','ylim',[0.5 2e4],'xlim',[0.5 2e4],'xtick',10.^(-1:5),'ytick',10.^(-1:5));%,'YTickMode','Manual','XTickMode','Manual')
        if n_ind == 1, ylabel('Threshold'); end
%         title(['K= ',num2str(K,5),' (',num2str(K_std,5),')']);
               
        xlabel('C_0')                    
    end 
end

% Threshold Concentration, Concentration Change and Fold Change: data, all compartments on the same graph     
figure('name','Model vs. data, Thresholds','position', [695   719   868   255]);
ColorPerTrainingCondition     = {'k','r'};
MarkerPerCompartment          = {'x','.','+'};
for n_ind = 1:length(CompartmentFieldNames)    
    CurrentNeuronField = CompartmentFieldNames{n_ind};   
    CurrentMarker      = MarkerPerCompartment{n_ind};
    
    for t_ind = 1:length(TrainingConditionNames)    
        CurrentTrainingCondition = TrainingConditionNames{t_ind};             
        CurrentColor     = ColorPerTrainingCondition{t_ind};

        K                       = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).K; 
        InitialConcetration_vec = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).InitialConcetrations; 
        T_vec                   = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).Threshold; 
        T_vec_std               = ModelParams.(CurrentTrainingCondition).(CurrentNeuronField).Bootstrap.Threshold.STD; 

        T_Error_High = T_vec_std;
        T_Error_Low  = T_vec_std;        
        
        C_Curve           =  exp(log(InitialConcetration_vec(1)):0.01:log(InitialConcetration_vec(end)));
        Tmodel_curve      =  K *(1-exp(-C_Curve/K));
        
        % Concentration change and fold change  
        %   Reminder: VAR(a(X+b))=a^2*VAR(X) --> std(C-T) = std(T), std((C-T)/C)=std(T)/C  
        ConcentrationChange             =  InitialConcetration_vec - T_vec;
        FoldChange                      = (InitialConcetration_vec - T_vec) ./ InitialConcetration_vec;
        ConcentrationChange_std         = T_vec_std;
        FoldChange_std                  = T_vec_std ./ InitialConcetration_vec;        
        
        % Plots                     
        subplot(1,3,1);
        errorbar(InitialConcetration_vec, T_vec, T_Error_Low, T_Error_High, 'linestyle','none','marker',CurrentMarker,'color',CurrentColor,'markerfacecolor',CurrentColor); hold on;
        set(gca,'xscale','log','yscale','log','ylim',[0.5 2e4],'xlim',[0.5 2e4],'xtick',10.^(-1:5),'ytick',10.^(-1:5));%,'YTickMode','Manual','XTickMode','Manual')
        if n_ind == 1, ylabel('Threshold'); end
        xlabel('C_0')      
        
        subplot(1,3,2);
        errorbar(InitialConcetration_vec, ConcentrationChange, ConcentrationChange_std, ConcentrationChange_std,'linestyle','none','marker',CurrentMarker,'color',CurrentColor,'markerfacecolor',CurrentColor); hold on;
        set(gca,'xscale','log','yscale','log','ylim',[1e-2 1e5],'xlim',[0.5 2e4],'xtick',10.^(-1:5),'ytick',10.^(-2:5));%,'YTickMode','Manual','XTickMode','Manual')
        if n_ind == 1, ylabel('Concentration Change'); end
        xlabel('C_0')     
                
        subplot(1,3,3);
        errorbar(InitialConcetration_vec, FoldChange, FoldChange_std, FoldChange_std, 'linestyle','none','marker',CurrentMarker,'color',CurrentColor,'markerfacecolor',CurrentColor); hold on;
        set(gca,'xscale','log','ylim',[0 1],'xlim',[0.5 2e4],'xtick',10.^(-1:5),'ytick',0:0.2:1);%,'YTickMode','Manual','XTickMode','Manual')
        if n_ind == 1, ylabel('Fold Change'); end
        xlabel('C_0')      
    end 
end

return

function Pvalues = PlotAdaptationTimes(ACT_CX17256, ACT_CX17255, ACT_Egl4, ACT_Egl4_ad450)

CompartmentFieldNames = {'Soma','Dendrite','Axon'};   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Calculation of adaptation times using various methods and parts of the data
figure('name','Calculation of adaptation times, various methods and parts of the data','position',[680   558   359   420])
FitNames   = {'FitBasedOnResponsivenessAndLatency_UseDataSTD','FitBasedOnOnePulseResponsivenessAndLatency','FitBasedOnTwoRepeatPerCondition',...
              'FitBasedOnAllResponsiveness','FitBasedOnOnePulseResponsiveness'};
NumberOfFits      = length(FitNames);
NumberOfIteration = length(ACT_CX17256.Soma.Bootstrap.FitBasedOnResponsivenessAndLatency_UseDataSTD.BestFitAdaptationTime);
c_ind = 1; 
CurrentCompartment = CompartmentFieldNames{c_ind};  
ACT        = ACT_CX17256;
Y_AllData  = ones(1,NumberOfFits,'single')*NaN;
DataArray  = ones(NumberOfIteration,NumberOfFits,'single')*NaN;
for f_ind = 1:NumberOfFits
    CurrentFit = FitNames{f_ind};
    if strcmpi(CurrentFit,'FitBasedOnTwoRepeatPerCondition')
        Y_AllData(f_ind) = nanmean(ACT.(CurrentCompartment).Bootstrap.(CurrentFit).BestFitAdaptationTime);                
    else
        Y_AllData(f_ind) = ACT.(CurrentCompartment).(CurrentFit).BestFitAdaptationTime;
    end
    BestFitAdaptationTime_BS = ACT.(CurrentCompartment).Bootstrap.(CurrentFit).BestFitAdaptationTime;
    DataArray(:,f_ind)       = BestFitAdaptationTime_BS;
end
boxplot_SL(DataArray,'BoxColor',[0.7 0.7 0.7],'BoxWidthValue',0.5,'WhiskerDisplay','ON','WhiskerPercentiles',[5 95],...
                     'MedianColor','r','MedianWidth',4); hold on;
plot(1:NumberOfFits, Y_AllData, 'ko','markerfacecolor','k'); hold on
xlim([0 NumberOfFits+1]); ylim([0 30])

% P-values by comparison of bootstrap distributions 
ReferenceIndex = 1;
Vec_Reference  = DataArray(:,ReferenceIndex);
Pvalue_SmallerThanReference = zeros(1,size(DataArray,2))*NaN;
Pvalue_LargerThanReference  = zeros(1,size(DataArray,2))*NaN;
for f_ind = 1:NumberOfFits 
    X = Vec_Reference;
    Y = DataArray(:,f_ind);
    [Pvalue_X_LargerThan_Y, Pvalue_Y_LargerThan_X] = CompareBootstrapVectors(X, Y);
    Pvalue_SmallerThanReference(f_ind) = Pvalue_X_LargerThan_Y;
    Pvalue_LargerThanReference(f_ind)  = Pvalue_Y_LargerThan_X;
end
Cutoff_Pvalue =  1/size(DataArray,1);                                                    % Lower Limit due to finite sampling
Pvalue_SmallerThanReference(Pvalue_SmallerThanReference<Cutoff_Pvalue) = Cutoff_Pvalue;  % Lower Limit due to finite sampling
Pvalue_LargerThanReference(Pvalue_LargerThanReference<Cutoff_Pvalue)   = Cutoff_Pvalue;  % Lower Limit due to finite sampling

Pvalue_CompareFits.ReferenceIndex              = ReferenceIndex;
Pvalue_CompareFits.ReferenceName               = FitNames{ReferenceIndex};
Pvalue_CompareFits.AllNames                    = FitNames;
Pvalue_CompareFits.Pvalue_SmallerThanReference = Pvalue_SmallerThanReference;
Pvalue_CompareFits.Pvalue_LargerThanReference  = Pvalue_LargerThanReference;

% FDR correction for multiple fits
NonReferenceIndex                           = setdiff(1:size(DataArray,2), ReferenceIndex);
Pvalue_CompareFits.FDR_SmallerThanReference = ComputeFDR_MultipleSignificanceLevels(Pvalue_SmallerThanReference, NonReferenceIndex);
Pvalue_CompareFits.FDR_LargerThanReference  = ComputeFDR_MultipleSignificanceLevels(Pvalue_LargerThanReference, NonReferenceIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Adaptation time in different compartments 
figure('name',' Adaptation time in different compartments','position',[680   558   243   420])
CompartmentFieldNames = {'Soma','Dendrite','Axon'};  
NumOfCompartments     = length(CompartmentFieldNames);
CurrentFit            = 'FitBasedOnResponsivenessAndLatency_UseDataSTD';
ACT                   = ACT_CX17256;
Y_AllData  = ones(1,NumOfCompartments,'single')*NaN;
DataArray  = ones(NumberOfIteration,NumOfCompartments,'single')*NaN;
for c_ind = 1:NumOfCompartments
    CurrentCompartment       = CompartmentFieldNames{c_ind};          
    Y_AllData(c_ind)         = ACT.(CurrentCompartment).(CurrentFit).BestFitAdaptationTime;    
    BestFitAdaptationTime_BS = ACT.(CurrentCompartment).Bootstrap.(CurrentFit).BestFitAdaptationTime;
    DataArray(:,c_ind)       = BestFitAdaptationTime_BS;
end
boxplot_SL(DataArray,'BoxColor',[0.7 0.7 0.7],'BoxWidthValue',0.7,'WhiskerDisplay','ON','WhiskerPercentiles',[5 95],...
                     'MedianColor','r','MedianWidth',4); hold on;
plot(1:NumOfCompartments, Y_AllData, 'ko','markerfacecolor','k'); hold on
xlim([0 NumOfCompartments+1]); 
ylim([0 30])

% P-values by comparison of bootstrap distributions 
ReferenceIndex = 1;
Vec_Reference  = DataArray(:,ReferenceIndex);
Pvalue_SmallerThanReference = zeros(1,size(DataArray,2))*NaN;
Pvalue_LargerThanReference  = zeros(1,size(DataArray,2))*NaN;
for f_ind = 1:NumOfCompartments 
    X = Vec_Reference;
    Y = DataArray(:,f_ind);
    [Pvalue_X_LargerThan_Y, Pvalue_Y_LargerThan_X] = CompareBootstrapVectors(X, Y);
    Pvalue_SmallerThanReference(f_ind) = Pvalue_X_LargerThan_Y;
    Pvalue_LargerThanReference(f_ind)  = Pvalue_Y_LargerThan_X;
end
Cutoff_Pvalue =  1/size(DataArray,1);                                                    % Lower Limit due to finite sampling
Pvalue_SmallerThanReference(Pvalue_SmallerThanReference<Cutoff_Pvalue) = Cutoff_Pvalue;  % Lower Limit due to finite sampling
Pvalue_LargerThanReference(Pvalue_LargerThanReference<Cutoff_Pvalue)   = Cutoff_Pvalue;  % Lower Limit due to finite sampling

Pvalue_CompareCompartments.ReferenceIndex              = ReferenceIndex;
Pvalue_CompareCompartments.ReferenceName               = CompartmentFieldNames{ReferenceIndex};
Pvalue_CompareCompartments.AllNames                    = CompartmentFieldNames;
Pvalue_CompareCompartments.Pvalue_SmallerThanReference = Pvalue_SmallerThanReference;
Pvalue_CompareCompartments.Pvalue_LargerThanReference  = Pvalue_LargerThanReference;

% FDR correction for multiple fits
NonReferenceIndex                                    = setdiff(1:size(DataArray,2), ReferenceIndex);
Pvalue_CompareCompartments.FDR_SmallerThanReference = ComputeFDR_MultipleSignificanceLevels(Pvalue_SmallerThanReference, NonReferenceIndex);
Pvalue_CompareCompartments.FDR_LargerThanReference  = ComputeFDR_MultipleSignificanceLevels(Pvalue_LargerThanReference, NonReferenceIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Calculation of relative adaptation times for different strains (relative to CX17255, the WT background of egl4 imaging lines) 
figure('name','Calculation of adaptation times for different strains','position',[680   558   359   420])
ACT_StructureNames = {'ACT_CX17255','ACT_Egl4','ACT_Egl4_ad450'};
NumberOfStrains    = length(ACT_StructureNames);
CurrentCompartment = 'Soma';  
CurrentFit         = 'FitBasedOnResponsivenessAndLatency_UseDataSTD';
Y_AllData  = ones(1,NumberOfStrains,'single')*NaN;
DataArray  = ones(NumberOfIteration,NumberOfStrains,'single')*NaN;
 
for act_ind = 1:NumberOfStrains      
    eval(['ACT=',ACT_StructureNames{act_ind},';']);
    Y_AllData(act_ind)         = ACT.(CurrentCompartment).(CurrentFit).BestFitAdaptationTime;    
    BestFitAdaptationTime_BS   = ACT.(CurrentCompartment).Bootstrap.(CurrentFit).BestFitAdaptationTime;
    DataArray(:,act_ind)       = BestFitAdaptationTime_BS;
end
DataArray = DataArray / Y_AllData(1); % Relative adaptation time
Y_AllData = Y_AllData / Y_AllData(1); % Relative adaptation time
boxplot_SL(DataArray,'BoxColor',[0.7 0.7 0.7],'BoxWidthValue',0.5,'WhiskerDisplay','ON','WhiskerPercentiles',[5 95],...
                     'MedianColor','r','MedianWidth',4); hold on;
plot(1:NumberOfStrains, Y_AllData, 'ko','markerfacecolor','k'); hold on
ylabel('Relative adaptation time')
xlim([0 NumberOfStrains+1]); ylim([0.1 500])
set(gca,'yscale','log')

% P-values by comparison of bootstrap distributions 
ReferenceIndex = 1;
Vec_Reference  = DataArray(:,ReferenceIndex);
Pvalue_SmallerThanReference = zeros(1,size(DataArray,2))*NaN;
Pvalue_LargerThanReference  = zeros(1,size(DataArray,2))*NaN;
for f_ind = 1:NumberOfStrains 
    X = Vec_Reference;
    Y = DataArray(:,f_ind);
    [Pvalue_X_LargerThan_Y, Pvalue_Y_LargerThan_X] = CompareBootstrapVectors(X, Y);
    Pvalue_SmallerThanReference(f_ind) = Pvalue_X_LargerThan_Y;
    Pvalue_LargerThanReference(f_ind)  = Pvalue_Y_LargerThan_X;
end
Cutoff_Pvalue =  1/size(DataArray,1);                                                    % Lower Limit due to finite sampling
Pvalue_SmallerThanReference(Pvalue_SmallerThanReference<Cutoff_Pvalue) = Cutoff_Pvalue;  % Lower Limit due to finite sampling
Pvalue_LargerThanReference(Pvalue_LargerThanReference<Cutoff_Pvalue)   = Cutoff_Pvalue;  % Lower Limit due to finite sampling

Pvalue_CompareStrains.ReferenceIndex              = ReferenceIndex;
Pvalue_CompareStrains.ReferenceName               = ACT_StructureNames{ReferenceIndex};
Pvalue_CompareStrains.AllNames                    = ACT_StructureNames;
Pvalue_CompareStrains.Pvalue_SmallerThanReference = Pvalue_SmallerThanReference;
Pvalue_CompareStrains.Pvalue_LargerThanReference  = Pvalue_LargerThanReference;

% FDR correction for multiple fits
NonReferenceIndex                              = setdiff(1:size(DataArray,2), ReferenceIndex);
Pvalue_CompareStrains.FDR_SmallerThanReference = ComputeFDR_MultipleSignificanceLevels(Pvalue_SmallerThanReference, NonReferenceIndex);
Pvalue_CompareStrains.FDR_LargerThanReference  = ComputeFDR_MultipleSignificanceLevels(Pvalue_LargerThanReference, NonReferenceIndex);


%% Assign to output structure
clear Pvalues
Pvalues.CompareFits         = Pvalue_CompareFits;
Pvalues.CompareCompartments = Pvalue_CompareCompartments;
Pvalues.CompareStrains      = Pvalue_CompareStrains;

return

function PlotLatencyPredictions_NotACT (ModelPredictions, MeasuredStats)
% Figure 1G and S2B and S2G
CompartmentFieldNames = {'Soma','Dendrite','Axon'};

%%  Plot Showing that the absolute concentration, derivative and fold-change models are not consistent with the data, using only MeasuredStats  
%%% Models are consistent only ABOVE the trace line where (Latency <= Peak time), but all data is below it 
% Note: Measured time to peak is similar for every compartment.
CurrentCompartment = 'Soma';
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

% LN model plot
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};  
    figure('name',[CurrentCompartment,' , LN model Latency predictions']);    

    Y      = ModelPredictions.(CurrentCompartment).LNmodel.Latency;
    Y_Mean = nanmean(Y,2);
    Y_SEM  = nanstd(Y,[],2) ./ sqrt(sum(~isnan(Y),2));
    
    X      = MeasuredStats.(CurrentCompartment).Latency.Matrix;
    X_Mean = nanmean(X,2);
    X_SEM  = nanstd(X,[],2) ./ sqrt(sum(~isnan(X),2));

    errorbar(X_Mean,Y_Mean,Y_SEM,Y_SEM,X_SEM,X_SEM, 'color','k','linestyle','none','marker','.');   hold on; %  errorbar(x,y,yneg,ypos,xneg,xpos)

    hold on; plot([0 100],[0 100],'k:'); axis([-3 51 -3 51])
    ylabel('Predicted Latency [s]')
    xlabel('Latency [s]')
    legend('LN model')
end

return

function [GOF_structure, Pvalues_RowModel_BetterThan_ColumnModel] = CompareModelPerformances(ModelPredictions, ACT, FitField, CompartmentFieldNames, PlotForVaryingCo)
%% Compare alternative models and ACT model predictions to measured data
% Initialization
if ~exist('PlotForVaryingCo','var')
    PlotForVaryingCo = false; % if true - plots from various initial conditions to buffer
end
% CompartmentFieldNames = {'Soma','Dendrite','Axon','AnyCompartment'};
Fields         = {'Responsiveness','Latency','C',              'dCdt',              'd2Cdt2',              'DeltaCoverC'}; 
MatchingFields = {'Responsiveness','Latency','C_at_Activation','dCdt_at_Activation','d2Cdt2_at_Activation','DeltaCoverC_at_Activation'}; 
MeasuredStats  = ModelPredictions.Information.MeasuredStats ;

% assign ACT model structure into ModelPredictions structure
for c_ind = 1:length(CompartmentFieldNames)        
    CurrentComparment = CompartmentFieldNames{c_ind};
    if isfield(ACT.(CurrentComparment),FitField)
        BestFitIndexInVec = ACT.(CurrentComparment).(FitField).BestFitIndexInVec;
        for f_ind = 1:length(Fields)
            ModelPredictions.(CurrentComparment).ACT.(Fields{f_ind}) = ...
                 squeeze(ACT.(CurrentComparment).(MatchingFields{f_ind})(BestFitIndexInVec,:,:));
        end   
    end
end

%% Plot models predictions. One figure per compartment, one subplot per comparison parameter, curves for each model. All concentration in microM
%  Related to Figure 2 and S3
ModelsToPlot  = {'C_FitStepResponseThreshold','dCdt','d2Cdt2','DeltaCoverC','ACT'};
ModelsColors  = {[0.5 0.5 0.5],[0.5 0.5 0.5],[0.5 0.5 0.5],[0.5 0.5 0.5],'r'};
ModelsMarkers = {'^','s','x','d','o'};
% ModelsLines   = {':',':',':',':','--'}; % For Matlab display
ModelsLines   = {'-','-','-','-','-'};  % For editing in Illustrator
NumberOfEstimatedParameters = [1 1 1 1 2]; 
if isfield(ModelPredictions.Soma,'LNmodel')
    ModelsToPlot  = {'C_FitStepResponseThreshold','dCdt','d2Cdt2','DeltaCoverC','LNmodel','ACT'};
    ModelsColors  = {[0.5 0.5 0.5],[0.5 0.5 0.5],[0.5 0.5 0.5],[0.5 0.5 0.5],[0.5 0.5 0.5],'r'};
    ModelsMarkers = {'^','s','x','d','+','o'};
%     ModelsLines   = {':',':',':',':',':','--'};  % For Matlab display
    ModelsLines   = {'-','-','-','-','-','-'};   % For editing in Illustrator
    NumberOfEstimatedParameters = [1 1 1 1 5 2]; % 5 dof are estimated for the impulse function, based on Kato et al., excluding 3 additional non-linearity parameters

end
FieldsInData  = {'Responsiveness','Latency','C_AtNeuronActivation','dCdt_AtNeuronActivation','d2Cdt2_AtNeuronActivation','DeltaCoverC_AtNeuronActivation'};
FieldsInModels= {'Responsiveness','Latency','C',                   'dCdt',                   'd2Cdt2',                   'DeltaCoverC'};
YLIMITS       = {[-8 108],         [-3 55],  [0 13],              [-0.7 0.1],               [-0.6 0.1],                  [-0.1 0.01]};
% ScalingIsRequired = [false, false, false, false, false, false];
% MicroM_ScalingFactor = 11.16e6;
C_Axis_For_Display    = ModelPredictions.Information.FinalConcentrationPerPulse;
C_Axis_For_Display(1) = 0.01;      % FOR DISPLAY OF "0" ON LOG SCALE, DON'T FORGET TO ADD A BREAK SIGN "//" on X AXIS !!!!
X_Axis_For_Display    = C_Axis_For_Display;
XLIMITS               = [0.007 20];

if PlotForVaryingCo % Plots from various initial conditions to buffer
    C_Axis_For_Display = ModelPredictions.Information.InitialConcentrationPerPulse;
    X_Axis_For_Display = C_Axis_For_Display;
    XLIMITS            = [0.5 200];
    YLIMITS            = {[-8 108], [-3 65], [-6 120], [-7 0.1], [-4.5 0.1], [-0.1 0.01]};
end

for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};           
    figure('name',CurrentCompartment); 
    
    for f_ind=1:length(FieldsInData)  
        subplot(2,3,f_ind)
        Measured       = MeasuredStats.(CurrentCompartment).(FieldsInData{f_ind}).Matrix;
        Measured_Mean  = nanmean(Measured,2);
        Measured_SEM   = nanstd(Measured,[],2) ./ sqrt(sum(~isnan(Measured),2));            
        errorbar(X_Axis_For_Display, Measured_Mean,  Measured_SEM, 'color','k','linestyle','-','marker','o','markerfacecolor','k','linewidth',2); hold on;
        
        for m_ind = 1:length(ModelsToPlot)
            CurrentModel   = ModelsToPlot{m_ind};
            Predicted      = ModelPredictions.(CurrentCompartment).(CurrentModel).(FieldsInModels{f_ind});
            Predicted_Mean = nanmean(Predicted,2);
                        
            % Show predicted DeltaCOverC even when it is less reliable (at very long times after stimulus due to very low dye signal and dye derivatives at late times)    
            if strcmpi(FieldsInModels{f_ind},'DeltaCoverC')
                NaNIndices = isnan(Predicted_Mean);
                if any(NaNIndices)
                   DeltaCOverC = nanmean(ModelPredictions.(CurrentCompartment).(CurrentModel).dCdt ./ ...
                                         ModelPredictions.(CurrentCompartment).(CurrentModel).C, 2);
                   Predicted_Mean(NaNIndices) = DeltaCOverC(NaNIndices);                                                      
                end
            end
%             Predicted_SEM  = nanstd(Predicted,[],2) ./ sqrt(sum(~isnan(Predicted),2));                           
%             errorbar(X_Axis_For_Display, Predicted_Mean, Predicted_SEM,...
%                      'color',ModelsColors{m_ind},'markerfacecolor',ModelsColors{m_ind},'linestyle',ModelsLines{m_ind},'marker',ModelsMarkers{m_ind})

            plot(X_Axis_For_Display, Predicted_Mean, ...
                     'color',ModelsColors{m_ind},'markerfacecolor',ModelsColors{m_ind},'linestyle',ModelsLines{m_ind},'marker',ModelsMarkers{m_ind})
        end
        set(gca,'xscale','log'); xlim(XLIMITS); ylim(YLIMITS{f_ind}); title(FieldsInModels{f_ind});
    end    
end

%% Calculate goodness of fit using BIC and AIC
StructureInitialization.All_GoodnessOfFit  = zeros(length(ModelsToPlot), length(FieldsInData))*NaN;
StructureInitialization.All_Smin_Div_dof   = zeros(length(ModelsToPlot), length(FieldsInData))*NaN;
StructureInitialization.All_dof            = zeros(length(ModelsToPlot), length(FieldsInData))*NaN;
StructureInitialization.All_BIC            = zeros(length(ModelsToPlot), length(FieldsInData))*NaN;
StructureInitialization.All_AIC            = zeros(length(ModelsToPlot), length(FieldsInData))*NaN;

for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};    
    GOF_structure.(CurrentCompartment) = StructureInitialization;
    
    for f_ind=1:length(FieldsInData)  
        Measured       = MeasuredStats.(CurrentCompartment).(FieldsInData{f_ind}).Matrix;
        % Data variance
        DataVariance                     = nanvar(Measured,0,2);                     % Measured variance at each experimental condition. WITH Bessel correction
        DataVarianceMat                  = repmat(DataVariance,1,size(Measured,2)); 
        DataVarianceMat(isnan(Measured)) = NaN;
        
        for m_ind = 1:length(ModelsToPlot)
            CurrentModel   = ModelsToPlot{m_ind};
            Predicted      = ModelPredictions.(CurrentCompartment).(CurrentModel).(FieldsInModels{f_ind});
            
            %% Compute model errors                
            ResidualsSquare                  = (Measured- Predicted).^2;                       % (Data-ModelFit)^2
            SminMAT                          = ResidualsSquare ./ DataVarianceMat;     % Residuals square relative to the expected (experimentally measured) variance
            SminMAT(ResidualsSquare==0)      = 0;
            Repeats                          = sum(~isnan(SminMAT),2);

            %% Compute statistics F- distribution   (Global - over all conditions)      
            dof                   = sum(Repeats) - NumberOfEstimatedParameters(m_ind); % degrees of freedom
            Smin                  = nansum(SminMAT(:));                         % Smin 
            GoodnessOfFit         = chi2cdf(Smin,dof,'upper');                        
            Smin_Div_dof          = Smin ./ dof;        
            BIC                   = Smin + NumberOfEstimatedParameters(m_ind)*log(sum(Repeats));
            AIC                   = Smin + NumberOfEstimatedParameters(m_ind)*2;

            %% Assign to matrices
            GOF_structure.(CurrentCompartment).All_GoodnessOfFit(m_ind,f_ind)    = GoodnessOfFit;
            GOF_structure.(CurrentCompartment).All_Smin_Div_dof(m_ind,f_ind)     = Smin_Div_dof;
            GOF_structure.(CurrentCompartment).All_dof(m_ind,f_ind)              = dof;
            GOF_structure.(CurrentCompartment).All_BIC(m_ind,f_ind)              = BIC;
            GOF_structure.(CurrentCompartment).All_AIC(m_ind,f_ind)              = AIC;

        end
    end    
end

%% Plot BIC and AIC (Figure 2 and S3 )
MaxForDisplay     = 1e3;                   % CutOff for display. Do not forget to add break sign "//" on Y axis !!
PerformanceFields = {'All_BIC','All_AIC'};
for p_ind = 1:length(PerformanceFields)
    CurrentPerformanceField = PerformanceFields{p_ind};
    figure('name',['Compare models to data using ',CurrentPerformanceField],'position',[280 479 1314 220]); 
    for c_ind = 1:length(CompartmentFieldNames)
        CurrentCompartment = CompartmentFieldNames{c_ind};   
        subplot(1,length(CompartmentFieldNames),c_ind);

        PerformanceMat = GOF_structure.(CurrentCompartment).(CurrentPerformanceField); % rows=models, columns=tested data parameters
        PerformanceMat(PerformanceMat>MaxForDisplay) = MaxForDisplay;

        for m_ind = 1:length(ModelsToPlot)
            plot(1:length(FieldsInData),PerformanceMat(m_ind,:),'color',ModelsColors{m_ind},'markerfacecolor',ModelsColors{m_ind},...
                 'linestyle',ModelsLines{m_ind},'marker',ModelsMarkers{m_ind}); hold on;
        end
        ylim([0 MaxForDisplay*1.1]); xlim([0.5 length(FieldsInData)+0.5]); title(CurrentCompartment);
        set(gca,'xtick',1:length(FieldsInData));
    end
end

%% Compute T statistics and compare to F-distrtibution
NumberForInfinity = 1e10;   % To avoid numeric errors
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};  
    Pvalues_RowModel_BetterThan_ColumnModel.(CurrentCompartment) = zeros(length(FieldsInData),length(ModelsToPlot),length(ModelsToPlot));
    All_dof          = GOF_structure.(CurrentCompartment).All_dof;
    All_Smin_Div_dof = GOF_structure.(CurrentCompartment).All_Smin_Div_dof;
    All_Smin_Div_dof(All_Smin_Div_dof==inf) = NumberForInfinity;
    
    for f_ind = 1:length(FieldsInData)         % Which parameter is tested     
        for row = 1:length(ModelsToPlot)       % Which row model
            dof_row        = All_dof(row,f_ind);
            SminDivdof_row = All_Smin_Div_dof(row,f_ind);
            for col = 1:length(ModelsToPlot)   % Which column model
                dof_col        = All_dof(col,f_ind);
                SminDivdof_col = All_Smin_Div_dof(col,f_ind);
                T              = SminDivdof_row ./ SminDivdof_col; % T stats
                if SminDivdof_col==0 && SminDivdof_row==0
                    T = 1;
                elseif SminDivdof_col==0
                    T = inf;
                elseif SminDivdof_row==0
                    T = 0;                         
                end            
                P = fcdf(T,dof_row,dof_col); % F distribution upper tail
                Pvalues_RowModel_BetterThan_ColumnModel.(CurrentCompartment)(f_ind, row, col) = P;  
            end
        end
    end
end

%% Plots based on T statistics (Figure S3)
for c_ind = 1:length(CompartmentFieldNames)
    CurrentCompartment = CompartmentFieldNames{c_ind};  
    figure('name',[CurrentCompartment, ' , -log(Pvalue) for Row model better than column model. -log(0.05)=3 , -log(10^-4)=9.2'],'position',[134  611  1702 150]);
    for f_ind = 1:length(FieldsInData)  % Which parameter is tested    
        subplot(1,length(FieldsInData),f_ind);
        MAT = squeeze(Pvalues_RowModel_BetterThan_ColumnModel.(CurrentCompartment)(f_ind,:,:));
        imagesc(-log(MAT)); set(gca,'clim',[0 10]); colorbar    
        title(FieldsInData{f_ind})
        set(gca,'xtick',[],'ytick',[])
    end
end

return 

function PlotDyeVariability(Data)
% Figure S1C
Dye  = Data.str2_CX17256_FromD6_Slow.Butanone_d6_To_Buffer.Dye;
Time = Data.str2_CX17256_FromD6_Slow.Butanone_d6_To_Buffer.Dye_TimeVec;
MicroM_ScalingFactor  = 11.16;  % in microM
Concentrations = MicroM_ScalingFactor*(1-Dye);
figure; plot(Time, Concentrations,'-','color',[0.7 0.7 0.7]); hold on; plot(Time,nanmean(Concentrations,1),'k-');

return

%% Display and statistics functions
function boxplot_SL(data, varargin)
%% Example for inputs
% figure; boxplot_SL(Distributions_At1microM','BoxColor',[0.7 0.7 0.7],'BoxWidthValue',0.7,'WhiskerDisplay','ON','WhiskerPercentiles',[5 95],'MedianColor','r','MedianWidth',2); hold on;

stringVars=varargin(1:2:end);
valueVars=varargin(2:2:end);

% Ind = strmatch('BoxType',stringVars);       % 'Filled' or 'Empty'.
% if ~isempty(Ind) 
%     BoxTypeValue = valueVars{Ind};
% else
%     BoxTypeValue = 'Filled';
% end
Ind = strmatch('BoxColor',stringVars);      % color sign or rgb
if ~isempty(Ind) 
    BoxColorValue = valueVars{Ind};
else
    BoxColorValue = 'Empty';
end
Ind = strmatch('WhiskerDisplay',stringVars); % 'ON' or 'OFF'.
if ~isempty(Ind) 
    WhiskerDisplayValue = valueVars{Ind};
else
    WhiskerDisplayValue = 'OFF';
end
Ind = strmatch('WhiskerPercentiles',stringVars); % [lower, upper] percentiles 
if ~isempty(Ind) 
    WhiskerPercentilesValue = valueVars{Ind};
else
    WhiskerPercentilesValue = 'OFF';
end
Ind = strmatch('MedianDisplay',stringVars);     % 'ON' or 'OFF'.
if ~isempty(Ind) 
    MedianDisplayValue = valueVars{Ind};
else
    MedianDisplayValue = 'ON';
end
Ind = strmatch('MedianColor',stringVars);        % color sign or rgb
if ~isempty(Ind) 
    MedianColorValue = valueVars{Ind};
else
    MedianColorValue = 'k';
end
Ind = strmatch('WhiskerColor',stringVars);        % color sign or rgb
if ~isempty(Ind) 
    WhiskerColorValue = valueVars{Ind};
else
    WhiskerColorValue = 'k';
end
Ind = strmatch('Width',stringVars);             % width of median, box, and whisker 
if ~isempty(Ind) 
    WidthValue = valueVars{Ind};
else
    WidthValue = 0.5;
end

NumOfGroups     = size(data,2);
% NumOfDataPoints = size(data,1);
HandlesMatrix = zeros(3, NumOfGroups); % whisker, Box, Median
% Xaxis = 1:NumOfGroups;

if ~strcmpi(WhiskerDisplayValue,'OFF')
    WhiskerMatrix = prctile(data,WhiskerPercentilesValue);
    for x = 1:NumOfGroups
        HandlesMatrix(1,x) = plot([x x],WhiskerMatrix(:,x),'color',WhiskerColorValue); hold on;
    end
end

% Box display    
BoxMatrix = prctile(data,[25 75]);
for x = 1:NumOfGroups
    % Box display
%     HandlesMatrix(2,x) = rectangle('position',[x-WidthValue/2 BoxMatrix(1,x) x+WidthValue/2 BoxMatrix(2,x)],WhiskerMatrix(:,x),...
%                                    'FaceColor',BoxColorValue,'EdgeColor',BoxColorValue);
    HandlesMatrix(2,x) = rectangle('position',[x-WidthValue/2 BoxMatrix(1,x) WidthValue BoxMatrix(2,x)-BoxMatrix(1,x)],...
                                   'FaceColor',BoxColorValue,'EdgeColor',BoxColorValue); hold on;
end

if ~strcmpi(MedianDisplayValue,'OFF')
    MedianVec = median(data,1);
    for x = 1:NumOfGroups
        HandlesMatrix(3,x) = plot([x-WidthValue/2 x+WidthValue/2],MedianVec(x)*ones(1,2),'color',MedianColorValue); hold on;
    end
end

return

function [Pvalue_X_LargerThan_Y, Pvalue_Y_LargerThan_X, DifferenceMatrix] = CompareBootstrapVectors(X, Y)

%% Comparison of bootstrapped distributions
%  Problem:     Two bootstrap vectors with values X=[X1, ... ,Xn] and Y=[Y1, ... , Ym]. 
%  Test Matrix: M(nxm) with  Mij = Xi - Yj
%  Hypothesis:  X is larger than Y.   Pvalue = length(find(M(:)<0)) / (nxm)  
%  Hypothesis:  X is smaller than Y.  Pvalue = length(find(M(:)>0)) / (nxm)  

%% Make sure X is row vector and Y is column vector
SizeX = size(X);
SizeY = size(Y);
if length(SizeX)~=2 || length(SizeY)~=2
    disp('X and Y must be 1D vectors');
    return
end
% X must be a column vector
if SizeX(2)>1       % not column vector
    if SizeX(1)==1  % row vector
        X = X';
    else            % matrix
        disp('X must be a 1D vector');
        return        
    end    
end
% Y must be a row vector
if SizeY(1)>1       % not row vector
    if SizeY(2)==1  % column vector
        Y = Y';
    else            % matrix
        disp('Y must be a 1D vector');
        return        
    end    
end

%% Compute test matrix
X = single(X); 
Y = single(Y);
LenX = length(X);
LenY = length(Y);

MatX = repmat(X,1,LenY);
MatY = repmat(Y,LenX,1);

DifferenceMatrix      = MatX - MatY;
Pvalue_X_LargerThan_Y = length(find(DifferenceMatrix(:)<0))/ (LenX*LenY);
Pvalue_Y_LargerThan_X = length(find(DifferenceMatrix(:)>0))/ (LenX*LenY);

return

function Significance = ComputeFDR(Pvalue_vec, RequiredFDRconfidence)
%% Example:
% RequiredFDRconfidence = 0.05; % required FDR confidence level
% Pvalue_vec            = [0.01 0.6 0.8 0.03 0.2];

% Pvalue_vec has to be a row vector 
if size(Pvalue_vec,1) > size(Pvalue_vec,2)
    Pvalue_vec = Pvalue_vec';
end
%% Computing post-hoc FDR for multiple analysis
m                 = length(Pvalue_vec);
[PvaluesSorted,I] = sort(Pvalue_vec);    
FDRval            = RequiredFDRconfidence.*(1:m)/m;
FDR_Mat           = [PvaluesSorted', (1:m)', FDRval']; 
LastDiscovery     = find(FDR_Mat(:,1)<=FDR_Mat(:,3),1,'last');
AllDiscoveryIndicesInOriginalVector = I(1:LastDiscovery);
Significance      = false(1,m);
Significance(AllDiscoveryIndicesInOriginalVector) = true;               

return

function SignificanceStructure = ComputeFDR_MultipleSignificanceLevels(Pvalue_vec, RelevantIndices)
% RelevantIndices --> should not include P-value of the reference! A value is not compared to itself. 
Pvalue_ToTest                                  = Pvalue_vec(RelevantIndices);
SignificanceStructure.FDR005                   = false(1,length(Pvalue_vec));
SignificanceStructure.FDR0005                  = false(1,length(Pvalue_vec));
SignificanceStructure.FDR0001                  = false(1,length(Pvalue_vec));
SignificanceStructure.FDR00005                 = false(1,length(Pvalue_vec));
SignificanceStructure.FDR005(RelevantIndices)  = ComputeFDR(Pvalue_ToTest, 0.05); 
SignificanceStructure.FDR0005(RelevantIndices) = ComputeFDR(Pvalue_ToTest, 0.005);
SignificanceStructure.FDR0001(RelevantIndices) = ComputeFDR(Pvalue_ToTest, 0.001);
SignificanceStructure.FDR00005(RelevantIndices)= ComputeFDR(Pvalue_ToTest, 0.0005);

return

function Djs = JehnsenShannonDivergence_inline (P,Q)
% P and Q are distribution vectors of the same size
% Djs is a scalar

% avoid numerical problems
P(P==0)=eps;
Q(Q==0)=eps;

M   = (P+Q)/2;
% Djs = nansum( 1/2* (P.*log2(P./M)) +  1/2* (Q.*log2(Q./M)) );  % THIS MEASURE IS GOOD FOR PROBABILITY DITRBUTIONS SINCE IT USES LOG-2 THAT SCALES IT [0-1] 
Djs = nansum( 1/2* (P.*log(P./M)) +  1/2* (Q.*log(Q./M)) );      % This is the general term

return

function [p, h, stats, EffectSize_XsmallerthanY] = ranksum_inline(x,y,varargin)
%RANKSUM Wilcoxon rank sum test for equal medians. Modified from original Matlab function 'ranksum'   
%   P = RANKSUM(X,Y) performs a two-sided rank sum test of the hypothesis
%   that two independent samples, in the vectors X and Y, come from
%   distributions with equal medians, and returns the p-value from the
%   test.  P is the probability of observing the given result, or one more
%   extreme, by chance if the null hypothesis ("medians are equal") is
%   true.  Small values of P cast doubt on the validity of the null
%   hypothesis.  The two sets of data are assumed to come from continuous
%   distributions that are identical except possibly for a location shift,
%   but are otherwise arbitrary.  X and Y can be different lengths.
%   RANKSUM treats NaNs in X or Y as missing values, and removes them.
%   The two-sided p-value is computed by doubling the most significant
%   one-sided value.
%
%   The Wilcoxon rank sum test is equivalent to the Mann-Whitney U test.
%
%   [P,H] = RANKSUM(...) returns the result of the hypothesis test,
%   performed at the 0.05 significance level, in H.  H=0 indicates that
%   the null hypothesis ("medians are equal") cannot be rejected at the 5%
%   level. H=1 indicates that the null hypothesis can be rejected at the
%   5% level.
%
%   [P,H] = RANKSUM(...,'alpha',ALPHA) returns the result of the hypothesis
%   test performed at the significance level ALPHA.
%
%   [P,H] = RANKSUM(...,'method',M) computes the p-value exactly if M is
%   'exact', or uses a normal approximation if M is 'approximate'.  If you
%   omit this argument, RANKSUM uses the exact method for small samples and
%   the approximate method for larger samples.
%
%   [P,H] = RANKSUM(...,'tail',TAIL) performs the test against the
%   alternative hypothesis specified by TAIL:
%       'both'  -- "medians are not equal" (two-tailed test, default)
%       'right' -- "median of X is greater than median of Y" (right-tailed test)
%       'left'  -- "median of X is less than median of Y" (left-tailed test)
%   TAIL must be a single string.
%
%   [P,H,STATS] = RANKSUM(...) returns STATS, a structure with one or two
%   fields.  The field 'ranksum' contains the value of the rank sum
%   statistic for X.  For the 'approximate' method, the field 'zval'
%   contains the value of the normal (Z) statistic.
%
%   See also SIGNTEST, SIGNRANK, KRUSKALWALLIS, TTEST2.

%   References:
%      [1] Hollander, M. and D. A. Wolfe.  Nonparametric Statistical
%          Methods. Wiley, 1973.
%      [2] Gibbons, J.D.  Nonparametric Statistical Inference,
%          2nd ed.  M. Dekker, 1985.


% Check that x and y are vectors
if ~isvector(x) || ~isvector(y)
   error(message('stats:ranksum:InvalidData'));
end

% Remove missing data
x = x( ~isnan(x) );
y = y( ~isnan(y) );
if isempty(x) || isempty(y)
	error(message('stats:signrank:NotEnoughData'));
end

% Determine value for 'alpha' and parse inputs
alpha = 0.05;   % default
if nargin>2 && isnumeric(varargin{1})
   % Grandfathered syntax:  ranksum(x,y,alpha)
   alpha = varargin{1};
   varargin(1) = [];
end
%
oknames = {'alpha' 'method' 'tail'};
dflts   = {alpha   ''   'both'};
[alpha,method,tail] = internal.stats.parseArgs(oknames,dflts,varargin{:});

% Check value of 'alpha'
if ~isscalar(alpha) || ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
   error(message('stats:ranksum:BadAlpha'));
end

% Check value of 'tail'
tail = internal.stats.getParamVal(tail, {'both'  'right'  'left'}, '''tail''');

% Determine and check value for 'method'
nx = numel(x);
ny = numel(y);
ns = min(nx,ny);
if isempty(method)
   if (ns < 10)  &&  ((nx+ny) < 20)
      method = 'exact';
   else
      method = 'approximate';
   end
elseif strcmpi(method, 'oldexact')
	method = 'oldexact';
else   % method not recognized, throw error
   method = internal.stats.getParamVal(method,{'exact' 'approximate'},'''method''');
end

% Determine computational 'technique'
switch method
	case 'approximate'
		technique = 'normal_approximation';
	case 'oldexact'
		technique = 'full_enumeration';
	case 'exact'
		if (nx+ny) < 10
			technique = 'full_enumeration';
		else
			technique = 'network_algorithm';
		end
end

%      %      %      %      %      %      %      %      %      %
% Calculations for Rank Sum Test

x = x(:);   % ensure columns
y = y(:);
if nx <= ny
   smsample = x;
   lgsample = y;
   same_order = true;
else
   smsample = y;
   lgsample = x;
   same_order = false;
end

% Compute the rank sum statistic based on the smaller sample
[ranks, tieadj] = tiedrank([smsample; lgsample]);
srank = ranks(1:ns);
w = sum(srank);


switch technique
	case 'full_enumeration'
		allpos = nchoosek(ranks,ns);   % enumerate all possibilities
		sumranks = sum(allpos,2);
		np = size(sumranks, 1);
		
		switch tail
			case 'both'
				plo = sum( sumranks <= w) / np ;
				phi = sum( sumranks >= w) / np ;
				p_tail = min(plo,phi);
				p = min(2*p_tail, 1);   % 2-sided, p>1 means middle is double-counted
				
			case 'right'
				switch same_order
					case true
						p = sum( sumranks >= w) / np ;
					case false
						p = sum( sumranks <= w) / np;
				end
				
			case 'left'
				switch same_order
					case true
						p = sum( sumranks <= w) / np ;
					case false
						p = sum( sumranks >= w) / np;
				end
				
		end
		
		%     %     %     %     %     %     %      %      %      %
		
	case 'network_algorithm'
		[p_net, pvals] = exactprob(smsample, lgsample, w);
		
		if any(isnan(p_net)) || any(isnan(pvals))
			warning(message('stats:ranksum:NanResult'));
			p = NaN;
			
		else
			switch tail
				case 'both'   % two-tailed test
					p = min(2*p_net, 1);   % p>1 means the middle is double-counted
					
				case 'right'   % right-tail test
					switch same_order
						case true
							p =  pvals(2) + pvals(3);
						case false
							p = pvals(2) + pvals(1);
					end
					
				case 'left'   % left-tail test
					switch same_order
						case true
							p =  pvals(2) + pvals(1);
						case false
							p = pvals(2) + pvals(3);
					end
					
			end   % conditional on 'tail'
		
		end
		
		%     %     %     %     %     %     %      %      %      %
				
	case 'normal_approximation'
		wmean = ns*(nx + ny + 1)/2;
		tiescor = 2 * tieadj / ((nx+ny) * (nx+ny-1));
		wvar  = nx*ny*((nx + ny + 1) - tiescor)/12;
		wc = w - wmean;

		% compute z-value, including continuity correction
		switch tail
			case 'both'
				z = (wc - 0.5 * sign(wc))/sqrt(wvar);
				if ~same_order
					z = -z;
				end
				p = 2*normcdf(-abs(z));
				
			case 'right'
				if same_order
					z = (wc - 0.5)/sqrt(wvar);
				else
					z = -(wc + 0.5)/sqrt(wvar);
				end
				
				p = normcdf(-z);
				
			case 'left'
				if same_order
					z = (wc + 0.5)/sqrt(wvar);
				else
					z = -(wc - 0.5)/sqrt(wvar);
				end
				
				p = normcdf(z);
		end


		if (nargout > 2)   % handle additional outputs
			stats.zval = z;
		end
		
end   % conditional on 'technique'



% Handle additional outputs
if nargout > 1,
   h = (p<=alpha);

   if (nargout > 2)
	   if same_order
		   stats.ranksum = w;
	   else
		   stats.ranksum = sum( ranks(ns+1:end) );
	   end
   end   
end

%% Calculate effect size
try
    [X,Y]                   = meshgrid(single(x),single(y));
    MAT_XsmallerthanY       = Y-X; 
    Fraction_XsmallerthanY  = length(find(MAT_XsmallerthanY(:)>0)) / length(MAT_XsmallerthanY(:));    % another way: 1-2*(length(find(~MAT_XsmallerthanY(:)))/2)/ (length(x)*length(y))
catch
    disp('Large vectors- calculating effect size in loop');
    NumberXsmallerthanY = 0;
    for ind=1:length(x)
        NumberXsmallerthanY = NumberXsmallerthanY + length(find(x(ind)<y));
    end   
    Fraction_XsmallerthanY = NumberXsmallerthanY/(length(x)*length(y));
end

EffectSize_XsmallerthanY = 2*Fraction_XsmallerthanY - 1;   % FractionSmaller- FractionHigher.  [-1 1]
    
return

function [p1, pvals] = exactprob(x,y,w)
%EXACTPROB Exact P-values for Wilcoxon Mann Whitney nonparametric test
%   [P1,PVALS]=EXACTPROB(X,Y,W) computes the p-value P for the test
%   statistic W in a Wilcoxon-Mann-Whitney nonparametric test of the
%   hypothesis that X and Y come from distributions with equal medians.

% Create a contingency table with element (i,j) indicating how many
% times u(j) appears in sample i
u = unique([x(:); y(:)]);
t = zeros(2,length(u));
t(1,:) = histc(x,u)';
t(2,:) = histc(y,u)';

% Compute weights for wmw test
colsum = sum(t,1);
tmp = cumsum(colsum);
wts = [0 tmp(1:end-1)] + .5*(1+diff([0 tmp]));

% Compute p-value using network algorithm for contingency tables
[p1, pvals] = statctexact_SL(t,wts,w);
return

function handles = barwitherr_inline(Y, errors, Properties)
% Input arguments
% Y         = mxn matrix. Designates m bar plots with n time points.
% errorbars = mxn matrix. Optional. Designates errors for the bars in Y.
% Properties = An optional structure 
%   Properties.BarProps.FieldsName       =  FieldsValue;
%   Properties.ErrorBarProps.FieldsName  =  FieldsValue;
%   Properties.gcaProps.FieldsName       =  FieldsValue;
%     
%   Examples for Properties structure values:
%     Properties.BarProps.EdgeColor      = {'r','b','g','c'};
%     Properties.BarProps.FaceColor      = {'r','b','g','c'};
%     Properties.ErrorBarProps.linewidth = {3,3,3,3};
%     Properties.ErrorBarProps.color     = {'r','b','g','c'};
%     Properties.ErrorBarProps.color     = {'k','k','k','k'};
%     Properties.gcaProps.xticklabel     = [0 1 2 3];   % corresponding to ticks of 1:size(Y,1)   
%     Properties.gcaProps.xlim           = [-0.1 3.1];
%     Properties.gcaProps.fontsize       = 20;
%     Properties.xlabelProps.String      = 'Time [hours]';
%     Properties.ylabelProps.String      = 'Learning Index';
%     Properties.titleProps.String       = 'CI';
%     Properties.legendProps.String      =  num2cell(UniqueDilutions);
% 

if ~exist('Properties','var')
    Properties.stam = 0;
end
if exist('errors','var') && ~isempty(errors)
    plot_errorbars = true;
else
    plot_errorbars = false;
end

NumOfBarPlots     = size(Y,1); 
Xaxis             = 1:size(Y,2);

if NumOfBarPlots>1
    InterTimePointShift = 0.3;
    BarToInterbarRatio  = 3;
    barwidth            = (1-InterTimePointShift)/(NumOfBarPlots+(NumOfBarPlots-1)/BarToInterbarRatio);
    Interbarwidth       = barwidth/BarToInterbarRatio;
    
    RelativeBarXshift = (barwidth/2):(barwidth+Interbarwidth):((barwidth/2)+(NumOfBarPlots-1)*(barwidth+Interbarwidth));
    RelativeBarXshift = RelativeBarXshift-mean(RelativeBarXshift);  % around 0   
else
    RelativeBarXshift = 0;
    barwidth          = 0.8;
end

handles.bars      = zeros(1,NumOfBarPlots);
handles.errorbars = zeros(1,NumOfBarPlots);

for bar_num = 1:NumOfBarPlots
    CurrentXaxis          = RelativeBarXshift(bar_num) + Xaxis;
    handles.bars(bar_num) = bar(CurrentXaxis,Y(bar_num,:),'barwidth',barwidth); hold on;
    
    if isfield(Properties,'BarProps')
        FieldNames = fieldnames(Properties.BarProps);
        for f_num = 1:length(FieldNames)
            set(handles.bars(bar_num),FieldNames(f_num),Properties.BarProps.(FieldNames{f_num})(bar_num));
        end
    end    
end

if plot_errorbars
    for bar_num = 1:NumOfBarPlots
        CurrentXaxis               = RelativeBarXshift(bar_num) + Xaxis;       
        handles.errorbars(bar_num) = errorbar(CurrentXaxis, Y(bar_num,:), errors(bar_num,:),'.');
                
        if isfield(Properties,'ErrorBarProps')
            FieldNames = fieldnames(Properties.ErrorBarProps);
            for f_num = 1:length(FieldNames)
                set(handles.errorbars(bar_num),FieldNames(f_num),Properties.ErrorBarProps.(FieldNames{f_num})(bar_num));
            end
        end
    end
end

handles.gca = gca; 
set(gca,'xtick',Xaxis);

if isfield(Properties,'gcaProps')
    FieldNames = fieldnames(Properties.gcaProps);
    for f_num = 1:length(FieldNames)
        set(gca,FieldNames{f_num},Properties.gcaProps.(FieldNames{f_num}));
    end
end
if isfield(Properties,'xlabelProps')
    FieldNames    = fieldnames(Properties.xlabelProps);
    xlabel_handle = xlabel(Properties.xlabelProps.String);
    for f_num = 1:length(FieldNames)
        set(xlabel_handle,FieldNames{f_num},Properties.xlabelProps.(FieldNames{f_num}));
    end
end
if isfield(Properties,'ylabelProps')
    FieldNames    = fieldnames(Properties.ylabelProps);
    ylabel_handle = ylabel(Properties.ylabelProps.String);
    for f_num = 1:length(FieldNames)
        set(ylabel_handle,FieldNames{f_num},Properties.ylabelProps.(FieldNames{f_num}));
    end
end
if isfield(Properties,'titleProps')
    FieldNames   = fieldnames(Properties.titleProps);
    title_handle = title(Properties.titleProps.String);
    for f_num = 1:length(FieldNames)
        set(title_handle,FieldNames{f_num},Properties.titleProps.(FieldNames{f_num}));
    end
end
if isfield(Properties,'legendProps')
    FieldNames    = fieldnames(Properties.legendProps);
    legend_handle = legend(Properties.legendProps.String);
    for f_num = 1:length(FieldNames)
        set(legend_handle,FieldNames{f_num},Properties.legendProps.(FieldNames{f_num}));
    end
end

return

%% Plots related to LN model details
function PlotImpulseFunctions(Data, RelevantStrainNames)
     
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);         

%% Plot impulse functions calculated from step inputs 11 microM --> buffer, for all compartments 
ColorPerCompartment = {'k','b',[1 0.7 0]};
for s_ind = 1:length(RelevantStrainNames)
    CurrentStrain = RelevantStrainNames{s_ind};
    PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    figure('name',[CurrentStrain,' Impulse functions']);
    
    CurrentPulse  = PulsesNames{1};  % to buffer

    for c_ind = 1:NumberOfCompartment
        CurrentColor           = ColorPerCompartment{c_ind};    
        CurrentCompartment     = CompartmentFieldNames{c_ind};
        Current_Activity_Field = [CurrentCompartment,'_ImpulseFunction_NoScaling'];
        ImpulseFunction        = Data.(CurrentStrain).(CurrentPulse).(Current_Activity_Field);
        X                      = (0:length(ImpulseFunction)-1)/10; % 10 Hz
        plot(X, ImpulseFunction,'-','color',CurrentColor); hold on;
        xlim([-0.1 8]);        
    end   
    legend(CompartmentFieldNames)

end

return

function PlotPredictedResponses_LNmodel_Traces(Data, RelevantStrainIndices) 

StrainNames           = fieldnames(Data);
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);    
XLIMITS               = [-2 70];

% Traces
for s_ind = 1:length(RelevantStrainIndices)
    CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
    PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    LegendPerPulseType  = num2cell(Data.(CurrentStrain).PulsesInfo.FinalConcentration);
    NumberOfPulses      = length(PulsesNames);         
    figure('name',[CurrentStrain,' Predicted (LNmodel) vs. measured responses'],'position', get(0,'ScreenSize'));       

    for c_ind = 1:NumberOfCompartment            
        CurrentCompartment     = CompartmentFieldNames{c_ind};
%       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
        MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
        PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF'];
        gcaHandles     = [];    
        CompartmentMax = 0;            
        
        for p_ind = 1:NumberOfPulses
            CurrentPulse  = PulsesNames{p_ind};
            if ~isfield(Data.(CurrentStrain),CurrentPulse)
                continue
            end
            subplot(NumberOfCompartment,NumberOfPulses, p_ind+(c_ind-1)*NumberOfPulses)
            TimeVec      = Data.(CurrentStrain).(CurrentPulse).TimeVec;
            MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
            PredictedMat = Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(:,1:length(TimeVec));             
            
            plot(TimeVec, MeasuredMat',           '-','color',[0.5 0.5 0.5]); hold on;
            plot(TimeVec, PredictedMat',          '-','color',[1 0.5 0.5]);   hold on;
            plot(TimeVec, nanmean(MeasuredMat,1), '-','color','k'); hold on;
            plot(TimeVec, nanmean(PredictedMat,1),'-','color','r'); hold on;
            CurrentXLIMITS = XLIMITS; CurrentXLIMITS(end)=min([CurrentXLIMITS(end) TimeVec(end)]);            
            xlim(CurrentXLIMITS);  
            if p_ind==1
                ylabel([CurrentCompartment,' \DeltaF/F']);
            end
            if c_ind==1
                title(LegendPerPulseType{p_ind});
            end
            gcaHandles     = [gcaHandles gca];
            CompartmentMax = max([CompartmentMax MeasuredMat(:)' PredictedMat(:)']);
        end 
        YLIMITS = [-0.2 CompartmentMax*1.1];
        set(gcaHandles,'ylim',YLIMITS)
    end    
    
end

return

function PlotPredictedResponses_LNmodel_Traces_SaturationScreen(Data, RelevantStrainIndices) 

StrainNames           = fieldnames(Data);
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);    
XLIMITS               = [-2 70];
Ymax                  = 17.4;
RelativeSaturation    = [1/3 1/2 1 2 3];
Ymax_vec              = Ymax*RelativeSaturation; 
% Traces
for s_ind = 1:length(RelevantStrainIndices)
    CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
    PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    LegendPerPulseType  = num2cell(Data.(CurrentStrain).PulsesInfo.FinalConcentration);
    NumberOfPulses      = length(PulsesNames);         

    for c_ind = 1:NumberOfCompartment            
        CurrentCompartment     = CompartmentFieldNames{c_ind};
        figure('name',[CurrentStrain,', ',CurrentCompartment,' Predicted (LNmodel) vs. measured responses'],'position', get(0,'ScreenSize'));       
%       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
        MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
        PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF_SaturationScreen'];
        gcaHandles = [];    
        CurrentMax = 0;            
        
        for p_ind = 1:NumberOfPulses
            CurrentPulse  = PulsesNames{p_ind};
            if ~isfield(Data.(CurrentStrain),CurrentPulse)
                continue
            end
            TimeVec      = Data.(CurrentStrain).(CurrentPulse).TimeVec;
            MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
            for y_ind = 1:length(Ymax_vec)
                CurrentRelativeSaturation  = RelativeSaturation(y_ind);
                PredictedMat = squeeze(Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(y_ind,:,1:length(TimeVec)));              
                subplot(length(Ymax_vec),NumberOfPulses, p_ind+(y_ind-1)*NumberOfPulses)

                plot(TimeVec, MeasuredMat',           '-','color',[0.5 0.5 0.5]); hold on;
                plot(TimeVec, PredictedMat',          '-','color',[1 0.5 0.5]);   hold on;
                plot(TimeVec, nanmean(MeasuredMat,1), '-','color','k'); hold on;
                plot(TimeVec, nanmean(PredictedMat,1),'-','color','r'); hold on;
                CurrentXLIMITS = XLIMITS; CurrentXLIMITS(end)=min([CurrentXLIMITS(end) TimeVec(end)]);            
                xlim(CurrentXLIMITS);  
                if p_ind==1
                    ylabel(['Relative Saturation=',num2str(CurrentRelativeSaturation),char(10),CurrentCompartment,' \DeltaF/F']);
                end
                if y_ind==1
                    title(LegendPerPulseType{p_ind});
                end
                gcaHandles     = [gcaHandles gca];
                CurrentMax = max([CurrentMax MeasuredMat(:)' PredictedMat(:)']);
            end
        end 
        YLIMITS = [-0.2 CurrentMax*1.1];
        set(gcaHandles,'ylim',YLIMITS)
    end       
end

return

function PlotPredictedResponses_LNmodel_Traces_GCaMP_Deconv(Data, RelevantStrainIndices) 

StrainNames           = fieldnames(Data);
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);    
XLIMITS               = [-2 70];

% Traces
for s_ind = 1:length(RelevantStrainIndices)
    CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
    PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
    LegendPerPulseType  = num2cell(Data.(CurrentStrain).PulsesInfo.FinalConcentration);
    NumberOfPulses      = length(PulsesNames);         
    figure('name',[CurrentStrain,' Predicted (LNmodel) vs. measured responses'],'position', get(0,'ScreenSize'));       

    for c_ind = 1:NumberOfCompartment            
        CurrentCompartment     = CompartmentFieldNames{c_ind};
%       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
        MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
        PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF_GCaMP_Deconv'];
        gcaHandles     = [];    
        CompartmentMax = 0;            
        
        for p_ind = 1:NumberOfPulses
            CurrentPulse  = PulsesNames{p_ind};
            if ~isfield(Data.(CurrentStrain),CurrentPulse)
                continue
            end
            subplot(NumberOfCompartment,NumberOfPulses, p_ind+(c_ind-1)*NumberOfPulses)
            TimeVec      = Data.(CurrentStrain).(CurrentPulse).TimeVec;
            MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
            PredictedMat = Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(:,1:length(TimeVec));             
            
            plot(TimeVec, MeasuredMat',           '-','color',[0.5 0.5 0.5]); hold on;
            plot(TimeVec, PredictedMat',          '-','color',[1 0.5 0.5]);   hold on;
            plot(TimeVec, nanmean(MeasuredMat,1), '-','color','k'); hold on;
            plot(TimeVec, nanmean(PredictedMat,1),'-','color','r'); hold on;
            CurrentXLIMITS = XLIMITS; CurrentXLIMITS(end)=min([CurrentXLIMITS(end) TimeVec(end)]);            
            xlim(CurrentXLIMITS);  
            if p_ind==1
                ylabel([CurrentCompartment,' \DeltaF/F']);
            end
            if c_ind==1
                title(LegendPerPulseType{p_ind});
            end
            gcaHandles     = [gcaHandles gca];
            CompartmentMax = max([CompartmentMax MeasuredMat(:)' PredictedMat(:)']);
        end 
        YLIMITS = [-0.2 CompartmentMax*1.1];
        set(gcaHandles,'ylim',YLIMITS)
    end    
    
end

return

function PlotPredictedResponses_LNmodel_Stats(Data, RelevantStrainIndices, PlotCompartmentsComparison, PlotConditionsComparison) 

StrainNames           = fieldnames(Data);
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);    

MaxNumberOfRepeats   = 40;  % free parameter. It just needs to be higher than the maximal number of repeats

if PlotCompartmentsComparison
    % Responsiveness vs. concentration for both measured and predicted  
    DynamicsRangeThresholdForResponsiveness = 0.2; 
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
        PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;  
        NumberOfPulses      = length(PulsesNames);         
        figure('name',[CurrentStrain,' Predicted (LNmodel, red) vs. measured (black) max responses'],'position', [384  486  1164  220]);       

        for c_ind = 1:NumberOfCompartment            
            CurrentCompartment     = CompartmentFieldNames{c_ind};
            MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
    %       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
            PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF'];
            MaxStructure.Measured.Matrix  = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;            
            MaxStructure.Measured.Mean    = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.STD     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.SEM     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Predicted        = MaxStructure.Measured;  

            for p_ind = 1:NumberOfPulses
                CurrentPulse  = PulsesNames{p_ind};
                if ~isfield(Data.(CurrentStrain),CurrentPulse)
                    continue
                end
                MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
                PredictedMat = Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(:,1:size(MeasuredMat,2));   

                CurrentMaxMeasuredVec  = nanmax(MeasuredMat,[],2);   % Result: max of different repeats of the same pulse
                CurrentMaxPredictedVec = nanmax(PredictedMat,[],2);  % Result: max of different repeats of the same pulse

                MaxStructure.Measured.Matrix(p_ind,1:length(CurrentMaxMeasuredVec))  = CurrentMaxMeasuredVec;
                MaxStructure.Predicted.Matrix(p_ind,1:length(CurrentMaxMeasuredVec)) = CurrentMaxPredictedVec;                      
            end 
            Mat                           = MaxStructure.Measured.Matrix;
            MatResponsiveness             = (Mat>DynamicsRangeThresholdForResponsiveness)*100;
            MatResponsiveness(isnan(Mat)) = NaN;
            MaxStructure.Measured.Mean    = nanmean(MatResponsiveness,2);    
            MaxStructure.Measured.STD     = nanstd(MatResponsiveness,[],2);    
            MaxStructure.Measured.SEM     = nanstd(MatResponsiveness,[],2) ./ sqrt(sum(~isnan(MatResponsiveness),2));    
            Mat                           = MaxStructure.Predicted.Matrix;
            MatResponsiveness             = (Mat>DynamicsRangeThresholdForResponsiveness)*100;
            MatResponsiveness(isnan(Mat)) = NaN;
            MaxStructure.Predicted.Mean   = nanmean(MatResponsiveness,2);    
            MaxStructure.Predicted.STD    = nanstd(MatResponsiveness,[],2);    
            MaxStructure.Predicted.SEM    = nanstd(MatResponsiveness,[],2) ./ sqrt(sum(~isnan(MatResponsiveness),2));    

            subplot(1,NumberOfCompartment,c_ind)
            FinalConcentrationPerPulse_ForDisplay = FinalConcentrationPerPulse;
            FinalConcentrationPerPulse_ForDisplay(FinalConcentrationPerPulse_ForDisplay==0) = 0.01;    % FOR DISPLAY ON LOG SCALE, DON'T FORGET TO ADD A BREAK SIGN "//" on X AXIS !!!!

            % Responsiveness vs. concentration for both measured and predicted   
            % Max DelatFoverF vs. concentration for both measured and predicted   
            errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Measured.Mean,  MaxStructure.Measured.SEM, 'ko-');   hold on; 
            plot(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Predicted.Mean,'rx-');   hold on; 
            set(gca,'xscale','log');
            ylim([-5 105])
            title(CurrentCompartment);
            if c_ind==1
                ylabel ('% responders')
            end
        end        
    end

    % Max DelatFoverF vs. concentration for both measured and predicted   
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
        PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;  
        NumberOfPulses      = length(PulsesNames);         
        figure('name',[CurrentStrain,' Predicted (LNmodel, red) vs. measured (black) max responses'],'position', [384  486  1164  220]);       

        for c_ind = 1:NumberOfCompartment            
            CurrentCompartment     = CompartmentFieldNames{c_ind};
            MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
    %       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
            PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF'];
            MaxStructure.Measured.Matrix  = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;            
            MaxStructure.Measured.Mean    = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.STD     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.SEM     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Predicted        = MaxStructure.Measured;  

            for p_ind = 1:NumberOfPulses
                CurrentPulse  = PulsesNames{p_ind};
                if ~isfield(Data.(CurrentStrain),CurrentPulse)
                    continue
                end
                MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
                PredictedMat = Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(:,1:size(MeasuredMat,2));   

                CurrentMaxMeasuredVec  = nanmax(MeasuredMat,[],2);   % Result: max of different repeats of the same pulse
                CurrentMaxPredictedVec = nanmax(PredictedMat,[],2);  % Result: max of different repeats of the same pulse

                MaxStructure.Measured.Matrix(p_ind,1:length(CurrentMaxMeasuredVec))  = CurrentMaxMeasuredVec;
                MaxStructure.Predicted.Matrix(p_ind,1:length(CurrentMaxMeasuredVec)) = CurrentMaxPredictedVec;                      
            end 
            Mat                         = MaxStructure.Measured.Matrix;
            MaxStructure.Measured.Mean  = nanmean(Mat,2);    
            MaxStructure.Measured.STD   = nanstd(Mat,[],2);    
            MaxStructure.Measured.SEM   = nanstd(Mat,[],2) ./ sqrt(sum(~isnan(Mat),2));    
            Mat                         = MaxStructure.Predicted.Matrix;
            MaxStructure.Predicted.Mean = nanmean(Mat,2);    
            MaxStructure.Predicted.STD  = nanstd(Mat,[],2);    
            MaxStructure.Predicted.SEM  = nanstd(Mat,[],2) ./ sqrt(sum(~isnan(Mat),2));    

            subplot(1,NumberOfCompartment,c_ind)
            FinalConcentrationPerPulse_ForDisplay = FinalConcentrationPerPulse;
            FinalConcentrationPerPulse_ForDisplay(FinalConcentrationPerPulse_ForDisplay==0) = 0.01;    % FOR DISPLAY ON LOG SCALE, DON'T FORGET TO ADD A BREAK SIGN "//" on X AXIS !!!!

            % Max DelatFoverF vs. concentration for both measured and predicted   
            errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Measured.Mean,  MaxStructure.Measured.SEM, 'color','k','linestyle','none','marker','.');   hold on; 
            errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Predicted.Mean, MaxStructure.Predicted.SEM,'color','r','linestyle','none','marker','.');   hold on; 
            set(gca,'xscale','log');

        end        
    end
end

if PlotConditionsComparison
    % Responsiveness vs. concentration for both measured and predicted  
    DynamicsRangeThresholdForResponsiveness = 0.2; 
    figure('name','Predicted (LNmodel, red) vs. measured (black) max responses','position', [623   200   429   699]);   
    c_ind = 1;
    CurrentCompartment     = CompartmentFieldNames{c_ind};
    
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
        PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration; 
        InitialConcentration       = Data.(CurrentStrain).PulsesInfo.InitialConcentration; 
        NumberOfPulses      = length(PulsesNames);         
        subplot(length(RelevantStrainIndices),1,s_ind)
          
        MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
%       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
        PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF'];
        MaxStructure.Measured.Matrix  = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;            
        MaxStructure.Measured.Mean    = zeros(NumberOfPulses,1,'single')*NaN;            
        MaxStructure.Measured.STD     = zeros(NumberOfPulses,1,'single')*NaN;            
        MaxStructure.Measured.SEM     = zeros(NumberOfPulses,1,'single')*NaN;            
        MaxStructure.Predicted        = MaxStructure.Measured;  

        for p_ind = 1:NumberOfPulses
            CurrentPulse  = PulsesNames{p_ind};
            if ~isfield(Data.(CurrentStrain),CurrentPulse)
                continue
            end
            MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
            PredictedMat = Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(:,1:size(MeasuredMat,2));   

            CurrentMaxMeasuredVec  = nanmax(MeasuredMat,[],2);   % Result: max of different repeats of the same pulse
            CurrentMaxPredictedVec = nanmax(PredictedMat,[],2);  % Result: max of different repeats of the same pulse

            MaxStructure.Measured.Matrix(p_ind,1:length(CurrentMaxMeasuredVec))  = CurrentMaxMeasuredVec;
            MaxStructure.Predicted.Matrix(p_ind,1:length(CurrentMaxMeasuredVec)) = CurrentMaxPredictedVec;                      
        end 
        Mat                           = MaxStructure.Measured.Matrix;
        MatResponsiveness             = (Mat>DynamicsRangeThresholdForResponsiveness)*100;
        MatResponsiveness(isnan(Mat)) = NaN;
        MaxStructure.Measured.Mean    = nanmean(MatResponsiveness,2);    
        MaxStructure.Measured.STD     = nanstd(MatResponsiveness,[],2);    
        MaxStructure.Measured.SEM     = nanstd(MatResponsiveness,[],2) ./ sqrt(sum(~isnan(MatResponsiveness),2));    
        Mat                           = MaxStructure.Predicted.Matrix;
        MatResponsiveness             = (Mat>DynamicsRangeThresholdForResponsiveness)*100;
        MatResponsiveness(isnan(Mat)) = NaN;
        MaxStructure.Predicted.Mean   = nanmean(MatResponsiveness,2);    
        MaxStructure.Predicted.STD    = nanstd(MatResponsiveness,[],2);    
        MaxStructure.Predicted.SEM    = nanstd(MatResponsiveness,[],2) ./ sqrt(sum(~isnan(MatResponsiveness),2));    

        FinalConcentrationPerPulse_ForDisplay = FinalConcentrationPerPulse;
        FinalConcentrationPerPulse_ForDisplay(FinalConcentrationPerPulse_ForDisplay==0) = 0.01;    % FOR DISPLAY ON LOG SCALE, DON'T FORGET TO ADD A BREAK SIGN "//" on X AXIS !!!!

        % Responsiveness vs. concentration for both measured and predicted   
        errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Measured.Mean,  MaxStructure.Measured.SEM, 'ko-');   hold on; 
        plot(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Predicted.Mean,'rx-');   hold on; 
        set(gca,'xscale','log');
        ylim([-5 105])
        title(InitialConcentration);
        ylabel ('% responders')                        
    end       
end



return

function PlotPredictedResponses_LNmodel_Stats_SaturationScreen(Data, RelevantStrainIndices, PlotCompartmentsComparison) 

StrainNames           = fieldnames(Data);
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);    

MaxNumberOfRepeats   = 40;  % free parameter. It just needs to be higher than the maximal number of repeats
Ymax                 = 17.4;
RelativeSaturation   = [1/3 1/2 1 2 3];
Ymax_vec             = Ymax*RelativeSaturation; 

if PlotCompartmentsComparison
    % Responsiveness vs. concentration for both measured and predicted  
    DynamicsRangeThresholdForResponsiveness = 0.2; 
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
        PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;  
        NumberOfPulses      = length(PulsesNames);         

        for c_ind = 1:NumberOfCompartment            
            CurrentCompartment     = CompartmentFieldNames{c_ind};
            figure('name',[CurrentStrain,', ',CurrentCompartment,' Predicted (LNmodel, red) vs. measured (black) max responses'],'position', [ 384   174   326   749]);       
            
            MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
    %       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
            PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF_SaturationScreen'];
            MaxStructure.Measured.Matrix  = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;            
            MaxStructure.Measured.Mean    = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.STD     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.SEM     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Predicted        = MaxStructure.Measured;  
            
            for y_ind = 1:length(Ymax_vec)
                CurrentRelativeSaturation  = RelativeSaturation(y_ind);
                for p_ind = 1:NumberOfPulses
                    CurrentPulse  = PulsesNames{p_ind};
                    if ~isfield(Data.(CurrentStrain),CurrentPulse)
                        continue
                    end
                    MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
                    PredictedMat = squeeze(Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(y_ind,:,1:size(MeasuredMat,2)));   

                    CurrentMaxMeasuredVec  = nanmax(MeasuredMat,[],2);   % Result: max of different repeats of the same pulse
                    CurrentMaxPredictedVec = nanmax(PredictedMat,[],2);  % Result: max of different repeats of the same pulse

                    MaxStructure.Measured.Matrix(p_ind,1:length(CurrentMaxMeasuredVec))  = CurrentMaxMeasuredVec;
                    MaxStructure.Predicted.Matrix(p_ind,1:length(CurrentMaxMeasuredVec)) = CurrentMaxPredictedVec;                      
                end 
                Mat                           = MaxStructure.Measured.Matrix;
                MatResponsiveness             = (Mat>DynamicsRangeThresholdForResponsiveness)*100;
                MatResponsiveness(isnan(Mat)) = NaN;
                MaxStructure.Measured.Mean    = nanmean(MatResponsiveness,2);    
                MaxStructure.Measured.STD     = nanstd(MatResponsiveness,[],2);    
                MaxStructure.Measured.SEM     = nanstd(MatResponsiveness,[],2) ./ sqrt(sum(~isnan(MatResponsiveness),2));    
                Mat                           = MaxStructure.Predicted.Matrix;
                MatResponsiveness             = (Mat>DynamicsRangeThresholdForResponsiveness)*100;
                MatResponsiveness(isnan(Mat)) = NaN;
                MaxStructure.Predicted.Mean   = nanmean(MatResponsiveness,2);    
                MaxStructure.Predicted.STD    = nanstd(MatResponsiveness,[],2);    
                MaxStructure.Predicted.SEM    = nanstd(MatResponsiveness,[],2) ./ sqrt(sum(~isnan(MatResponsiveness),2));    

                subplot(length(Ymax_vec),1,y_ind)
                FinalConcentrationPerPulse_ForDisplay = FinalConcentrationPerPulse;
                FinalConcentrationPerPulse_ForDisplay(FinalConcentrationPerPulse_ForDisplay==0) = 0.01;    % FOR DISPLAY ON LOG SCALE, DON'T FORGET TO ADD A BREAK SIGN "//" on X AXIS !!!!

                % Responsiveness vs. concentration for both measured and predicted   
                % Max DelatFoverF vs. concentration for both measured and predicted   
                errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Measured.Mean,  MaxStructure.Measured.SEM, 'ko-');   hold on; 
                plot(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Predicted.Mean,'rx-');   hold on; 
                set(gca,'xscale','log');
                ylim([-5 105])
                title(['Relative Saturation=',num2str(CurrentRelativeSaturation)]);
                ylabel('% responders');
            end            
        end        
    end

    % Max DelatFoverF vs. concentration for both measured and predicted   
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
        PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;  
        NumberOfPulses      = length(PulsesNames);         

        for c_ind = 1:NumberOfCompartment            
            CurrentCompartment     = CompartmentFieldNames{c_ind};
            figure('name',[CurrentStrain,', ',CurrentCompartment,' Predicted (LNmodel, red) vs. measured (black) max responses'],'position', [ 384   174   326   749]);       
            
            MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
    %       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
            PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF_SaturationScreen'];
            MaxStructure.Measured.Matrix  = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;            
            MaxStructure.Measured.Mean    = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.STD     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.SEM     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Predicted        = MaxStructure.Measured;  

            for y_ind = 1:length(Ymax_vec)
                CurrentRelativeSaturation  = RelativeSaturation(y_ind);
                for p_ind = 1:NumberOfPulses
                    CurrentPulse  = PulsesNames{p_ind};
                    if ~isfield(Data.(CurrentStrain),CurrentPulse)
                        continue
                    end
                    MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
                    PredictedMat = squeeze(Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(y_ind,:,1:size(MeasuredMat,2)));   

                    CurrentMaxMeasuredVec  = nanmax(MeasuredMat,[],2);   % Result: max of different repeats of the same pulse
                    CurrentMaxPredictedVec = nanmax(PredictedMat,[],2);  % Result: max of different repeats of the same pulse

                    MaxStructure.Measured.Matrix(p_ind,1:length(CurrentMaxMeasuredVec))  = CurrentMaxMeasuredVec;
                    MaxStructure.Predicted.Matrix(p_ind,1:length(CurrentMaxMeasuredVec)) = CurrentMaxPredictedVec;                      
                end 
                Mat                         = MaxStructure.Measured.Matrix;
                MaxStructure.Measured.Mean  = nanmean(Mat,2);    
                MaxStructure.Measured.STD   = nanstd(Mat,[],2);    
                MaxStructure.Measured.SEM   = nanstd(Mat,[],2) ./ sqrt(sum(~isnan(Mat),2));    
                Mat                         = MaxStructure.Predicted.Matrix;
                MaxStructure.Predicted.Mean = nanmean(Mat,2);    
                MaxStructure.Predicted.STD  = nanstd(Mat,[],2);    
                MaxStructure.Predicted.SEM  = nanstd(Mat,[],2) ./ sqrt(sum(~isnan(Mat),2));    

                subplot(length(Ymax_vec),1,y_ind)
                FinalConcentrationPerPulse_ForDisplay = FinalConcentrationPerPulse;
                FinalConcentrationPerPulse_ForDisplay(FinalConcentrationPerPulse_ForDisplay==0) = 0.01;    % FOR DISPLAY ON LOG SCALE, DON'T FORGET TO ADD A BREAK SIGN "//" on X AXIS !!!!

                % Max DelatFoverF vs. concentration for both measured and predicted   
                errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Measured.Mean,  MaxStructure.Measured.SEM, 'color','k','linestyle','none','marker','.');   hold on; 
                errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Predicted.Mean, MaxStructure.Predicted.SEM,'color','r','linestyle','none','marker','.');   hold on; 
                set(gca,'xscale','log');
                ylabel('max(\DeltaF/F)');
                title(['Relative Saturation=',num2str(CurrentRelativeSaturation)]);
            end
        end        
    end
end

return

function PlotPredictedResponses_LNmodel_Stats_GCaMP_Deconv(Data, RelevantStrainIndices, PlotCompartmentsComparison, PlotConditionsComparison) 

StrainNames           = fieldnames(Data);
CompartmentFieldNames = {'Soma','Dendrite','Axon'};
NumberOfCompartment   = length(CompartmentFieldNames);    

MaxNumberOfRepeats   = 40;  % free parameter. It just needs to be higher than the maximal number of repeats

if PlotCompartmentsComparison
    % Responsiveness vs. concentration for both measured and predicted  
    DynamicsRangeThresholdForResponsiveness = 0.2; 
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
        PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;  
        NumberOfPulses      = length(PulsesNames);         
        figure('name',[CurrentStrain,' Predicted (LNmodel, red) vs. measured (black) max responses'],'position', [384  486  1164  220]);       

        for c_ind = 1:NumberOfCompartment            
            CurrentCompartment     = CompartmentFieldNames{c_ind};
            MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
    %       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
            PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF_GCaMP_Deconv'];
            MaxStructure.Measured.Matrix  = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;            
            MaxStructure.Measured.Mean    = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.STD     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.SEM     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Predicted        = MaxStructure.Measured;  

            for p_ind = 1:NumberOfPulses
                CurrentPulse  = PulsesNames{p_ind};
                if ~isfield(Data.(CurrentStrain),CurrentPulse)
                    continue
                end
                MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
                PredictedMat = Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(:,1:size(MeasuredMat,2));   

                CurrentMaxMeasuredVec  = nanmax(MeasuredMat,[],2);   % Result: max of different repeats of the same pulse
                CurrentMaxPredictedVec = nanmax(PredictedMat,[],2);  % Result: max of different repeats of the same pulse

                MaxStructure.Measured.Matrix(p_ind,1:length(CurrentMaxMeasuredVec))  = CurrentMaxMeasuredVec;
                MaxStructure.Predicted.Matrix(p_ind,1:length(CurrentMaxMeasuredVec)) = CurrentMaxPredictedVec;                      
            end 
            Mat                           = MaxStructure.Measured.Matrix;
            MatResponsiveness             = (Mat>DynamicsRangeThresholdForResponsiveness)*100;
            MatResponsiveness(isnan(Mat)) = NaN;
            MaxStructure.Measured.Mean    = nanmean(MatResponsiveness,2);    
            MaxStructure.Measured.STD     = nanstd(MatResponsiveness,[],2);    
            MaxStructure.Measured.SEM     = nanstd(MatResponsiveness,[],2) ./ sqrt(sum(~isnan(MatResponsiveness),2));    
            Mat                           = MaxStructure.Predicted.Matrix;
            MatResponsiveness             = (Mat>DynamicsRangeThresholdForResponsiveness)*100;
            MatResponsiveness(isnan(Mat)) = NaN;
            MaxStructure.Predicted.Mean   = nanmean(MatResponsiveness,2);    
            MaxStructure.Predicted.STD    = nanstd(MatResponsiveness,[],2);    
            MaxStructure.Predicted.SEM    = nanstd(MatResponsiveness,[],2) ./ sqrt(sum(~isnan(MatResponsiveness),2));    

            subplot(1,NumberOfCompartment,c_ind)
            FinalConcentrationPerPulse_ForDisplay = FinalConcentrationPerPulse;
            FinalConcentrationPerPulse_ForDisplay(FinalConcentrationPerPulse_ForDisplay==0) = 0.01;    % FOR DISPLAY ON LOG SCALE, DON'T FORGET TO ADD A BREAK SIGN "//" on X AXIS !!!!

            % Responsiveness vs. concentration for both measured and predicted   
            % Max DelatFoverF vs. concentration for both measured and predicted   
            errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Measured.Mean,  MaxStructure.Measured.SEM, 'ko-');   hold on; 
            plot(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Predicted.Mean,'rx-');   hold on; 
            set(gca,'xscale','log');
            ylim([-5 105])
            title(CurrentCompartment);
            if c_ind==1
                ylabel ('% responders')
            end
        end        
    end

    % Max DelatFoverF vs. concentration for both measured and predicted   
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
        PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration;  
        NumberOfPulses      = length(PulsesNames);         
        figure('name',[CurrentStrain,' Predicted (LNmodel, red) vs. measured (black) max responses'],'position', [384  486  1164  220]);       

        for c_ind = 1:NumberOfCompartment            
            CurrentCompartment     = CompartmentFieldNames{c_ind};
            MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
    %       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
            PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF_GCaMP_Deconv'];
            MaxStructure.Measured.Matrix  = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;            
            MaxStructure.Measured.Mean    = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.STD     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Measured.SEM     = zeros(NumberOfPulses,1,'single')*NaN;            
            MaxStructure.Predicted        = MaxStructure.Measured;  

            for p_ind = 1:NumberOfPulses
                CurrentPulse  = PulsesNames{p_ind};
                if ~isfield(Data.(CurrentStrain),CurrentPulse)
                    continue
                end
                MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
                PredictedMat = Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(:,1:size(MeasuredMat,2));   

                CurrentMaxMeasuredVec  = nanmax(MeasuredMat,[],2);   % Result: max of different repeats of the same pulse
                CurrentMaxPredictedVec = nanmax(PredictedMat,[],2);  % Result: max of different repeats of the same pulse

                MaxStructure.Measured.Matrix(p_ind,1:length(CurrentMaxMeasuredVec))  = CurrentMaxMeasuredVec;
                MaxStructure.Predicted.Matrix(p_ind,1:length(CurrentMaxMeasuredVec)) = CurrentMaxPredictedVec;                      
            end 
            Mat                         = MaxStructure.Measured.Matrix;
            MaxStructure.Measured.Mean  = nanmean(Mat,2);    
            MaxStructure.Measured.STD   = nanstd(Mat,[],2);    
            MaxStructure.Measured.SEM   = nanstd(Mat,[],2) ./ sqrt(sum(~isnan(Mat),2));    
            Mat                         = MaxStructure.Predicted.Matrix;
            MaxStructure.Predicted.Mean = nanmean(Mat,2);    
            MaxStructure.Predicted.STD  = nanstd(Mat,[],2);    
            MaxStructure.Predicted.SEM  = nanstd(Mat,[],2) ./ sqrt(sum(~isnan(Mat),2));    

            subplot(1,NumberOfCompartment,c_ind)
            FinalConcentrationPerPulse_ForDisplay = FinalConcentrationPerPulse;
            FinalConcentrationPerPulse_ForDisplay(FinalConcentrationPerPulse_ForDisplay==0) = 0.01;    % FOR DISPLAY ON LOG SCALE, DON'T FORGET TO ADD A BREAK SIGN "//" on X AXIS !!!!

            % Max DelatFoverF vs. concentration for both measured and predicted   
            errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Measured.Mean,  MaxStructure.Measured.SEM, 'color','k','linestyle','none','marker','.');   hold on; 
            errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Predicted.Mean, MaxStructure.Predicted.SEM,'color','r','linestyle','none','marker','.');   hold on; 
            set(gca,'xscale','log');

        end        
    end
end

if PlotConditionsComparison
    % Responsiveness vs. concentration for both measured and predicted  
    DynamicsRangeThresholdForResponsiveness = 0.2; 
    figure('name','Predicted (LNmodel, red) vs. measured (black) max responses','position', [623   200   429   699]);   
    c_ind = 1;
    CurrentCompartment     = CompartmentFieldNames{c_ind};
    
    for s_ind = 1:length(RelevantStrainIndices)
        CurrentStrain = StrainNames{RelevantStrainIndices(s_ind)};
        PulsesNames         = Data.(CurrentStrain).PulsesInfo.PulseFromButanone;
        FinalConcentrationPerPulse = Data.(CurrentStrain).PulsesInfo.FinalConcentration; 
        InitialConcentration       = Data.(CurrentStrain).PulsesInfo.InitialConcentration; 
        NumberOfPulses      = length(PulsesNames);         
        subplot(length(RelevantStrainIndices),1,s_ind)
          
        MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF_MA_Smoothed']; %             
%       MeasuredActivity_Field  = [CurrentCompartment,'_DeltaFOverF'];           
        PredictedActivity_Field = [CurrentCompartment,'_LNmodelPredictedDeltaFoverF_GCaMP_Deconv'];
        MaxStructure.Measured.Matrix  = zeros(NumberOfPulses,MaxNumberOfRepeats,'single')*NaN;            
        MaxStructure.Measured.Mean    = zeros(NumberOfPulses,1,'single')*NaN;            
        MaxStructure.Measured.STD     = zeros(NumberOfPulses,1,'single')*NaN;            
        MaxStructure.Measured.SEM     = zeros(NumberOfPulses,1,'single')*NaN;            
        MaxStructure.Predicted        = MaxStructure.Measured;  

        for p_ind = 1:NumberOfPulses
            CurrentPulse  = PulsesNames{p_ind};
            if ~isfield(Data.(CurrentStrain),CurrentPulse)
                continue
            end
            MeasuredMat  = Data.(CurrentStrain).(CurrentPulse).(MeasuredActivity_Field); 
            PredictedMat = Data.(CurrentStrain).(CurrentPulse).(PredictedActivity_Field)(:,1:size(MeasuredMat,2));   

            CurrentMaxMeasuredVec  = nanmax(MeasuredMat,[],2);   % Result: max of different repeats of the same pulse
            CurrentMaxPredictedVec = nanmax(PredictedMat,[],2);  % Result: max of different repeats of the same pulse

            MaxStructure.Measured.Matrix(p_ind,1:length(CurrentMaxMeasuredVec))  = CurrentMaxMeasuredVec;
            MaxStructure.Predicted.Matrix(p_ind,1:length(CurrentMaxMeasuredVec)) = CurrentMaxPredictedVec;                      
        end 
        Mat                           = MaxStructure.Measured.Matrix;
        MatResponsiveness             = (Mat>DynamicsRangeThresholdForResponsiveness)*100;
        MatResponsiveness(isnan(Mat)) = NaN;
        MaxStructure.Measured.Mean    = nanmean(MatResponsiveness,2);    
        MaxStructure.Measured.STD     = nanstd(MatResponsiveness,[],2);    
        MaxStructure.Measured.SEM     = nanstd(MatResponsiveness,[],2) ./ sqrt(sum(~isnan(MatResponsiveness),2));    
        Mat                           = MaxStructure.Predicted.Matrix;
        MatResponsiveness             = (Mat>DynamicsRangeThresholdForResponsiveness)*100;
        MatResponsiveness(isnan(Mat)) = NaN;
        MaxStructure.Predicted.Mean   = nanmean(MatResponsiveness,2);    
        MaxStructure.Predicted.STD    = nanstd(MatResponsiveness,[],2);    
        MaxStructure.Predicted.SEM    = nanstd(MatResponsiveness,[],2) ./ sqrt(sum(~isnan(MatResponsiveness),2));    

        FinalConcentrationPerPulse_ForDisplay = FinalConcentrationPerPulse;
        FinalConcentrationPerPulse_ForDisplay(FinalConcentrationPerPulse_ForDisplay==0) = 0.01;    % FOR DISPLAY ON LOG SCALE, DON'T FORGET TO ADD A BREAK SIGN "//" on X AXIS !!!!

        % Responsiveness vs. concentration for both measured and predicted   
        errorbar(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Measured.Mean,  MaxStructure.Measured.SEM, 'ko-');   hold on; 
        plot(FinalConcentrationPerPulse_ForDisplay, MaxStructure.Predicted.Mean,'rx-');   hold on; 
        set(gca,'xscale','log');
        ylim([-5 105])
        title(InitialConcentration);
        ylabel ('% responders')                        
    end       
end



return

function PlotOptimizeRectifierForLNmodel(Data)
% Figures S2H and S2I
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
    
%%%% Plot traces - used for Figure S2 display for transition to buffer and 2e-7  %%%%    
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
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Model 2

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
    
%%%% Plot traces - used for Figure S2 display for transition to buffer and 2e-7  %%%%    
ModelName   = 'LNmodelWithRec_2';
b_vec_ticks = [1:20:81 99];
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


return


