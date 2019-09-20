function NoiseFilteringAnalysis_Fig6
%
%  Written by Sagi Levy, September 2019
%
%% General theoretical considerations
% WGN  = White Gaussian Noise
% Input Signals to explore:
%     S1 = S0
%     S2 = S0 - b*t
%     S3 = S0*exp(-a*t)
% Signals with noise:  
%     Si --> Si + WGN
% Models and relevant noise variance:
%       S < T(absolute)             NoiseVar(absolute)    = var(    S+WGN                - S )        = var(WGN).  
%      dS/dt < T(derivative)        NoiseVar(derivative)  = var(  d(S+WGN)/dt            - dS/dt )    = var(dWGN/dt)  
%     (dS/dt)/S < T(fold-change)    NoiseVar(fold-change) = var( (d(S+WGN)/dt))/(S+WGN) - (dS/dt)/S )   
%       S < T(ACT)                  NoiseVar(ACT)         = var(    S+WGN                - S )        = var(WGN).  
% 
%     NOTES:   For Si=constant -->  NoiseVar(fold-change) = var((d(WGN)/dt))/(S+WGN))

% Model SNR criteria (ability to filter noise by buffering noise-induced "false alarms")
%  General:  SNR(t) = (signal-threshold)^2/ var(noise)
%            When signal, threshold, and/or noise variance are time-dependent, we can compute:  mean(SNR(t)), over relevant time course 
%            When signal and threshold are time-dependent but noise isn't, mean can be computed only on the numerator. 
%            For the derivative model:  signal, var(noise) --> d(signal)/t, var(d(noise)/dt).  All the rest is the same
%            The variance of noise in fold-change is different: it depends on the signal, signal derivative, var(noise) and var(noise derivative) 
%               Therefore, var(fold-change noise) is time dependent if signal is time dependent,  
%               even if the variance(noise) and var(noise derivative) are constant in time. 
%
%  Exact definitions: 
%     SNR(absolute)    = mean(  (S - T(absolute))^2   / NoiseVar(absolute) )  
%                      = mean(  (S - T(absolute))^2 ) / var(WGN)    
%     SNR(derivative)  = mean(  (dS/dt - T(derivative))^2   / NoiseVar(derivative) )  
%                      = mean(  (dS/dt - T(derivative))^2 ) / var(d(WGN)/dt)    
%     SNR(fold-change) = mean( ((dS/dt)/S - T(fold-change))^2 / NoiseVar(fold-change) )  
%                        This needs to be computed by simulating noise.   NOTE: NoiseVar(fold-change) is time-dependent!!  
%     SNR(ACT)         = mean(  (S - T(ACT))^2   / NoiseVar(ACT) )  
%                      = mean(  (S - T(ACT))^2 ) / var(WGN)    
%  Notes:
%     1. For constant signals:   
%        SNR(derivative ) = T(derivative)^2  / NoiseVar(derivative) 
%        SNR(fold-change) = T(fold-change)^2 / NoiseVar(fold-change) 
%                           and for constant signal: NoiseVar(fold-change)= var((dn(t)/dt)/(S(t)+n(t))    
%     2. Averaging for mean(SNR(t))is done only over the time in which the relevant variable is larger than its threshold,   
%        For absolute model, for example, only for time in which S(t)> T(absolute).

%% Output .mat files:    
% 'NoiseAnalysis_Figure6.mat'    SNR and speed calculated for ACT model and for the absolute, derivative and fold-change models at various conditions   
%                                For details on output variables see comments below.   
NoiseAnalysis_part1_FileName  = 'D:\Models\NoiseAnalysis_part1.mat';     % Intermediate file for part 1        
NoiseAnalysis_part2_FileName  = 'D:\Models\NoiseAnalysis_part2.mat';     % Intermediate file for part 2     
NoiseAnalysis_FileName        = 'D:\Models\NoiseAnalysis_Figure6.mat';   % All parts, all variables 

%% Initialization
%%%  Free Parameters  %%%
SampleRate_TimeVec  = 5000; % Hz
TotalTime           = 60;  % seconds
Time                = single(0:1/SampleRate_TimeVec:TotalTime);
length_Time         = length(Time);

% Thresholds information - use previously extracted parameter information  
Tabsolute   = 5.5;   % microM
Tderivative = -0.28; % microM/sec
Tfoldchange = -0.06; % 1/sec
K           = 5.5;   % microM
Tao         = 17;    % seconds

%%%  Vectors For Screen  %%%
% MinimalWormSampleRate = 5;                   % Hz. Sample rate cannot be slower than our observed response times - around 200 msec
SampleRate_vec = [5000 500 50 5];              % 1/DeltaT
S0_vec         = [1 10 1e2 1e3 1e4];           % in microM.     For constant, linear decay and exponential decay 
NoiseSTD_vec   = 10.^(-5:0.05:6);              % in microM
b_vec          = [0.01 0.1 1 10 1e2 1e3 1e4];  % in microM/sec. For linear decay
a_vec          = [0.01 0.02 0.05 0.1 1 10];    % in 1/sec.      For exponential decay

rng('default')
NoiseVec_STD1 = randn(1,length_Time,'single'); 
iterations    = 1e2;
NoiseMat_STD1 = randn(iterations,length_Time,'single'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Part 1: Calculate SNR given different inputs and noise levels for various simple models (without low-pass filters) %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  Output structure: SNR.(InputType).(ModelType)  , where ModelType = Absolute, Derivative, FoldChange, ACT 
%         SNR.ConstantInput.(ModelType) = 3-D matrix (Input amplitude, Sample rate, noise level)
%         SNR.LinearInput.(ModelType)   = 4-D matrix (Input amplitude, Sample rate, decay rate(b_vec), noise level)
%         SNR.ExponentInput.(ModelType) = 4-D matrix (Input amplitude, Sample rate, decay rate(a_vec), noise level)
%  Structure is saved in 'NoiseAnalysis_part1_FileName.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Input 1: S1 = S0  %%%
disp([datestr(now),'   --   Computing SNR for constant signal']);
SNR.ConstantInput.Absolute   = NaN* ones(length(S0_vec), length(SampleRate_vec), length(NoiseSTD_vec),'single');
SNR.ConstantInput.Derivative = NaN* ones(length(S0_vec), length(SampleRate_vec), length(NoiseSTD_vec),'single');
SNR.ConstantInput.FoldChange = NaN* ones(length(S0_vec), length(SampleRate_vec), length(NoiseSTD_vec),'single');
SNR.ConstantInput.ACT        = NaN* ones(length(S0_vec), length(SampleRate_vec), length(NoiseSTD_vec),'single');    
for S0_ind = 1:length(S0_vec)
    S0 = S0_vec(S0_ind);    
%         ACT_ForConstantSignal = K*(1-exp(-S0/K));    
    for SR_ind = 1:length(SampleRate_vec)
        CurrentSampleRate    = SampleRate_vec(SR_ind);
        CurrentDeltaT        = 1/CurrentSampleRate;
        IntervalsInTimeVec   = round(SampleRate_TimeVec / CurrentSampleRate);
        CurrentNoiseVec_STD1 = NoiseVec_STD1(:,1:IntervalsInTimeVec:end);
        CurrentTime          = Time(1:IntervalsInTimeVec:end);

        % Correct noise to prevent negative concentration values  (S+noise >= 0) 
        CurrentNoiseMat                      = NoiseSTD_vec' * CurrentNoiseVec_STD1;      
        CurrentNoiseMat(CurrentNoiseMat<=-S0)= -0.9999*S0;
        CurrentNoiseDerivativeMat            = [diff(CurrentNoiseMat,[],2), NaN*ones(length(NoiseSTD_vec),1)] / CurrentDeltaT;             
        CurrentNoiseSTD                      = nanstd(CurrentNoiseMat,[],2);
        CurrentNoiseDerivativeSTD            = nanstd(CurrentNoiseDerivativeMat,[],2);                                

        %%% Absolute
        DistanceFromThreshold                         = S0 - Tabsolute;
        DistanceFromThreshold(DistanceFromThreshold<0)= NaN; % Look only on false alarms (passing the threshold when not supposed to)          
        SNR_PerNoiseSTD_Absolute                      = (DistanceFromThreshold)^2 ./ CurrentNoiseSTD.^2;
        SNR.ConstantInput.Absolute(S0_ind,SR_ind,:)   = SNR_PerNoiseSTD_Absolute;

        %%% ACT                                   
        InputMatWithNoise = S0+CurrentNoiseMat; InputMatWithNoise(:,1)=S0;
        for n_ind = 1:length(CurrentNoiseSTD)
            InputWithNoise = InputMatWithNoise(n_ind,:);
            CurrentACT = CalculateAdaptiveThreshold_inline (CurrentTime, InputWithNoise, K, 1/Tao); % full vector
            DistanceFromThreshold_OverTime                  = S0 - CurrentACT;     % Signal = Input(without noise) - Threshold
            DistanceFromThreshold_OverTime(DistanceFromThreshold_OverTime<0)= 0;   % make sure signal is zero whenever a "miss" is identified          
            AverageDistanceSquareFromThreshold              = nanmean(DistanceFromThreshold_OverTime.^2);
            SNR_PerNoiseSTD_ACT                             = AverageDistanceSquareFromThreshold ./ CurrentNoiseSTD(n_ind).^2;  % Use non-filtered noise for ACT model                            
            SNR.ConstantInput.ACT(S0_ind,SR_ind,n_ind)      = SNR_PerNoiseSTD_ACT;        
        end

        %%% Derivative
        DistanceFromThreshold                         = 0 - Tderivative;
        DistanceFromThreshold(DistanceFromThreshold<0)= NaN; % Look only on false alarms (passing the threshold when not supposed to)          
        SNR_PerNoiseSTD_Derivative                    = (DistanceFromThreshold)^2 ./ CurrentNoiseDerivativeSTD.^2;
        SNR.ConstantInput.Derivative(S0_ind,SR_ind,:) = SNR_PerNoiseSTD_Derivative;        

        %%% Fold-Change
        DistanceFromThreshold                         = 0 - Tfoldchange;
        DistanceFromThreshold(DistanceFromThreshold<0)= NaN; % Look only on false alarms (passing the threshold when not supposed to)          
        % compute  NoiseVar(fold-change), and for Si=constant -->  var((d(WGN)/dt))/(Si+WGN))
        if ~isnan(DistanceFromThreshold)
            FoldChangeNoiseMat  = CurrentNoiseDerivativeMat ./ (S0+CurrentNoiseMat);
            NoiseVar_FoldChange = nanvar(FoldChangeNoiseMat,[],2);                
            SNR_PerNoiseSTD_FoldChange                     = (DistanceFromThreshold)^2 ./ NoiseVar_FoldChange;
            SNR.ConstantInput.FoldChange(S0_ind,SR_ind,:)  = SNR_PerNoiseSTD_FoldChange;                                   
       end                                          
    end        
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Input 2: S2 = S0 - b*t  %%%
disp([datestr(now),'   --   Computing mean(SNR(t)) for linear decaying signal']);    
SNR.LinearInput.Absolute   = NaN* ones(length(S0_vec), length(SampleRate_vec), length(b_vec), length(NoiseSTD_vec),'single');
SNR.LinearInput.Derivative = NaN* ones(length(S0_vec), length(SampleRate_vec), length(b_vec), length(NoiseSTD_vec),'single');
SNR.LinearInput.FoldChange = NaN* ones(length(S0_vec), length(SampleRate_vec), length(b_vec), length(NoiseSTD_vec),'single');
SNR.LinearInput.ACT        = NaN* ones(length(S0_vec), length(SampleRate_vec), length(b_vec), length(NoiseSTD_vec),'single');
for S0_ind = 1:length(S0_vec)
    S0 = S0_vec(S0_ind);
    disp([datestr(now),'   --   Computing mean(SNR(t)) for linear decaying signal from ',num2str(S0),' microM']);
    for SR_ind = 1:length(SampleRate_vec)
        CurrentSampleRate    = SampleRate_vec(SR_ind);
        CurrentDeltaT        = 1/CurrentSampleRate;
        IntervalsInTimeVec   = round(SampleRate_TimeVec / CurrentSampleRate);
        CurrentTime          = Time(1:IntervalsInTimeVec:end);

        CurrentNoiseMatForRepeats_STD1 = NoiseMat_STD1(:,1:IntervalsInTimeVec:end);
        CurrentNoiseVec_STD1           = NoiseVec_STD1(:,1:IntervalsInTimeVec:end);            
        CurrentNoiseMat_BeforeNonNegativeCorrection  = NoiseSTD_vec' * CurrentNoiseVec_STD1;               

        for b_ind = 1:length(b_vec)
            Current_b      = b_vec(b_ind);
            Input          = S0 - Current_b*CurrentTime;
            Input(Input<0) = NaN;                
            InputMat       = repmat(Input,length(NoiseSTD_vec),1);
            InputMat2      = repmat(Input,iterations,1);

            % Correct noise to prevent negative concentration values  (S+noise >= 0) 
            CurrentNoiseMat                            = CurrentNoiseMat_BeforeNonNegativeCorrection;      
            CurrentNoiseMat(CurrentNoiseMat<=-InputMat)= -0.9999*InputMat(CurrentNoiseMat<=-InputMat);
            CurrentNoiseDerivativeMat     = [diff(CurrentNoiseMat,[],2), NaN*ones(length(NoiseSTD_vec),1)] / CurrentDeltaT;             
            CurrentNoiseSTD               = nanstd(CurrentNoiseMat,[],2);
            CurrentNoiseDerivativeSTD     = nanstd(CurrentNoiseDerivativeMat,[],2);           

            %%% Absolute
            DistanceFromThreshold_OverTime                  = Input - Tabsolute;   % Signal = Input - Threshold
            DistanceFromThreshold_OverTime(DistanceFromThreshold_OverTime<0)= NaN; % Look only on false alarms (passing the threshold when not supposed to)          
            AverageDistanceSquareFromThreshold              = nanmean(DistanceFromThreshold_OverTime.^2);
            SNR_PerNoiseSTD_Absolute                        = AverageDistanceSquareFromThreshold ./ CurrentNoiseSTD.^2;
            SNR.LinearInput.Absolute(S0_ind,SR_ind,b_ind,:) = SNR_PerNoiseSTD_Absolute;    

            %%% ACT                                   
            CurrentACTWithoutNoise = CalculateAdaptiveThreshold_ExpSaturation_BeforeActivation (CurrentTime, Input, K, 1/Tao); % NaNs from activation time
            DistanceFromThreshold_OverTime_WithoutNoise = Input - CurrentACTWithoutNoise;  % Signal = Input - Threshold
            FirstNaN          = find(DistanceFromThreshold_OverTime_WithoutNoise<0,1,'first');  % Look only on false alarms (passing the threshold when not supposed to)  
            InputMatWithNoise = InputMat+CurrentNoiseMat; InputMatWithNoise(:,1)=InputMat(:,1);
            for n_ind = 1:length(CurrentNoiseSTD)
                InputWithNoise = InputMatWithNoise(n_ind,:);
                CurrentACT = CalculateAdaptiveThreshold_inline (CurrentTime, InputWithNoise, K, 1/Tao); % full vector
                DistanceFromThreshold_OverTime                                  = Input - CurrentACT;      % Signal = Input(without noise) - Threshold
                DistanceFromThreshold_OverTime(DistanceFromThreshold_OverTime<0)= 0;    % make sure signal is zero whenever a "miss" is identified          
                DistanceFromThreshold_OverTime(FirstNaN:end)                    = NaN;  % Look only on false alarms (passing the threshold when not supposed to)           
                AverageDistanceSquareFromThreshold              = nanmean(DistanceFromThreshold_OverTime.^2);
                SNR_PerNoiseSTD_ACT                             = AverageDistanceSquareFromThreshold ./ CurrentNoiseSTD(n_ind).^2;  % Use non-filtered noise for ACT model                            
                SNR.LinearInput.ACT(S0_ind,SR_ind,b_ind,n_ind)  = SNR_PerNoiseSTD_ACT;        
            end                

            %%% Derivative
            DistanceFromThreshold                           = - Current_b - Tderivative;  % Signal = Input - Threshold
            if DistanceFromThreshold > 0      % Look only on false alarms (passing the threshold when not supposed to)  
                SNR_PerNoiseSTD_Derivative                        = DistanceFromThreshold^2 ./ CurrentNoiseDerivativeSTD.^2;
                SNR.LinearInput.Derivative(S0_ind,SR_ind,b_ind,:) = SNR_PerNoiseSTD_Derivative;    
            end

            %%% Fold-Change
            FoldChangeInputWithoutNoise                     = -Current_b./Input;
            DistanceFromThreshold_OverTime                  = FoldChangeInputWithoutNoise - Tfoldchange;  % Signal = Input - Threshold
            FirstForNaN                                     = find(DistanceFromThreshold_OverTime<0,1,'first');
            DistanceFromThreshold_OverTime(FirstForNaN:end) = NaN; % Look only on false alarms (passing the threshold when not supposed to)  
            DistanceFromThreshold_OverTime_Squared          = DistanceFromThreshold_OverTime.^2; 

            if isempty(FirstForNaN) || FirstForNaN > 1
                % For each noise STD, compute NoiseVar(fold-change) at each time point, extract SNR(t), compute and store mean(SNR(t))      
                FoldChangeInputWithoutNoiseMat    = repmat(FoldChangeInputWithoutNoise,iterations,1);
                SNR_PerNoiseSTD_FoldChange        = NaN*ones(1,length(CurrentNoiseSTD),'single');    

                for n_ind = 1:length(NoiseSTD_vec)
                    CurrentNoiseMatForRepeats = CurrentNoiseMatForRepeats_STD1 * CurrentNoiseSTD(n_ind);                                        
                    % Correct noise to prevent negative concentration values  (S+noise >= 0) 
                    CurrentNoiseMatForRepeats(CurrentNoiseMatForRepeats<=-InputMat2) = -0.9999*InputMat2(CurrentNoiseMatForRepeats<=-InputMat2);
                    CurrentNoiseDerivativeMatForRepeats     = [diff(CurrentNoiseMatForRepeats,[],2), NaN*ones(iterations,1)] / CurrentDeltaT;  
                    CurrentNoiseMatForRepeats(:,FirstForNaN:end)           = NaN;
                    CurrentNoiseDerivativeMatForRepeats(:,FirstForNaN:end) = NaN;                       

    %                 NoiseVariance = var( FoldChangeInputWithNoise - FoldChangeInputWithoutNoise );
                    FoldChangeInputWithNoise  = (-Current_b + CurrentNoiseDerivativeMatForRepeats)./(InputMat2+CurrentNoiseMatForRepeats);             
                    WithMinusWithoutNoiseMat  = FoldChangeInputWithNoise - FoldChangeInputWithoutNoiseMat;
                    VarianceInTime            = var(WithMinusWithoutNoiseMat,[],1);  % ERGODIC assumption: run variance over repeats at a given time 't' gives variance(t)
                    SNR_PerNoiseSTD_FoldChange(n_ind) = nanmean(DistanceFromThreshold_OverTime_Squared ./ VarianceInTime); % mean(SNR(t))
                end                 
                SNR.LinearInput.FoldChange(S0_ind,SR_ind,b_ind,:) = SNR_PerNoiseSTD_FoldChange;                             
            end
        end
    end        
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Input 3: S3 = S0*exp(-a*t)  %%%  
SNR.ExponentInput.Absolute   = NaN* ones(length(S0_vec), length(SampleRate_vec), length(a_vec), length(NoiseSTD_vec),'single');
SNR.ExponentInput.Derivative = NaN* ones(length(S0_vec), length(SampleRate_vec), length(a_vec), length(NoiseSTD_vec),'single');
SNR.ExponentInput.FoldChange = NaN* ones(length(S0_vec), length(SampleRate_vec), length(a_vec), length(NoiseSTD_vec),'single');
SNR.ExponentInput.ACT        = NaN* ones(length(S0_vec), length(SampleRate_vec), length(a_vec), length(NoiseSTD_vec),'single');
for S0_ind = 1:length(S0_vec)
    S0 = S0_vec(S0_ind);   
    disp([datestr(now),'   --   Computing mean(SNR(t)) for exponential decaying signal from ',num2str(S0),' microM']);
    for SR_ind = 1:length(SampleRate_vec)
        CurrentSampleRate    = SampleRate_vec(SR_ind);
        CurrentDeltaT        = 1/CurrentSampleRate;
        IntervalsInTimeVec   = round(SampleRate_TimeVec / CurrentSampleRate);
        CurrentTime          = Time(1:IntervalsInTimeVec:end);                       

        CurrentNoiseMatForRepeats_STD1 = NoiseMat_STD1(:,1:IntervalsInTimeVec:end);
        CurrentNoiseVec_STD1           = NoiseVec_STD1(:,1:IntervalsInTimeVec:end);            
        CurrentNoiseMat_BeforeNonNegativeCorrection  = NoiseSTD_vec' * CurrentNoiseVec_STD1;              

        for a_ind = 1:length(a_vec)
            Current_a      = a_vec(a_ind);
            Input          = S0 * exp(-Current_a*CurrentTime);                                
            InputMat       = repmat(Input,length(NoiseSTD_vec),1);
            InputMat2      = repmat(Input,iterations,1);

            % Correct noise to prevent negative concentration values  (S+noise >= 0) 
            CurrentNoiseMat                            = CurrentNoiseMat_BeforeNonNegativeCorrection;      
            CurrentNoiseMat(CurrentNoiseMat<=-InputMat)= -0.9999*InputMat(CurrentNoiseMat<=-InputMat);
            CurrentNoiseDerivativeMat     = [diff(CurrentNoiseMat,[],2), NaN*ones(length(NoiseSTD_vec),1)] / CurrentDeltaT;             
            CurrentNoiseSTD               = nanstd(CurrentNoiseMat,[],2);
            CurrentNoiseDerivativeSTD     = nanstd(CurrentNoiseDerivativeMat,[],2);         

            %%% Absolute
            DistanceFromThreshold_OverTime                    = Input - Tabsolute;   % Signal = Input - Threshold
            DistanceFromThreshold_OverTime(DistanceFromThreshold_OverTime<0)= NaN; % Look only on false alarms (passing the threshold when not supposed to)          
            AverageDistanceSquareFromThreshold                = nanmean(DistanceFromThreshold_OverTime.^2);
            SNR_PerNoiseSTD_Absolute                          = AverageDistanceSquareFromThreshold ./ CurrentNoiseSTD.^2;
            SNR.ExponentInput.Absolute(S0_ind,SR_ind,a_ind,:) = SNR_PerNoiseSTD_Absolute;    

            %%% ACT                                   
            CurrentACTWithoutNoise = CalculateAdaptiveThreshold_ExpSaturation_BeforeActivation (CurrentTime, Input, K, 1/Tao); % NaNs from activation time
            DistanceFromThreshold_OverTime_WithoutNoise = Input - CurrentACTWithoutNoise;  % Signal = Input - Threshold
            FirstNaN          = find(DistanceFromThreshold_OverTime_WithoutNoise<0,1,'first');  % Look only on false alarms (passing the threshold when not supposed to)  
            InputMatWithNoise = InputMat+CurrentNoiseMat; InputMatWithNoise(:,1)=InputMat(:,1);
            for n_ind = 1:length(CurrentNoiseSTD)
                InputWithNoise = InputMatWithNoise(n_ind,:);
                CurrentACT = CalculateAdaptiveThreshold_inline (CurrentTime, InputWithNoise, K, 1/Tao); % full vector
                DistanceFromThreshold_OverTime                                  = Input - CurrentACT;      % Signal = Input(without noise) - Threshold
                DistanceFromThreshold_OverTime(DistanceFromThreshold_OverTime<0)= 0;    % make sure signal is zero whenever a "miss" is identified          
                DistanceFromThreshold_OverTime(FirstNaN:end)                    = NaN;  % Look only on false alarms (passing the threshold when not supposed to)           
                AverageDistanceSquareFromThreshold              = nanmean(DistanceFromThreshold_OverTime.^2);
                SNR_PerNoiseSTD_ACT                             = AverageDistanceSquareFromThreshold ./ CurrentNoiseSTD(n_ind).^2;  % Use non-filtered noise for ACT model                            
                SNR.ExponentInput.ACT(S0_ind,SR_ind,a_ind,n_ind)= SNR_PerNoiseSTD_ACT;        
            end

            %%% Derivative
            DistanceFromThreshold_OverTime                      = - Current_a*Input - Tderivative;  % Signal = Input - Threshold
            FirstForNaN                                         = find(DistanceFromThreshold_OverTime<0,1,'first');
            DistanceFromThreshold_OverTime(FirstForNaN:end)     = NaN;                 % Look only on false alarms (passing the threshold when not supposed to)                      
            AverageDistanceSquareFromThreshold                  = nanmean(DistanceFromThreshold_OverTime.^2);            
            SNR_PerNoiseSTD_Derivative                          = AverageDistanceSquareFromThreshold ./ CurrentNoiseDerivativeSTD.^2;
            SNR.ExponentInput.Derivative(S0_ind,SR_ind,a_ind,:) = SNR_PerNoiseSTD_Derivative;    

            %%% Fold-Change
            FoldChangeInputWithoutNoise       = -Current_a;                        
            DistanceFromThreshold             = FoldChangeInputWithoutNoise - Tfoldchange;  % Signal = Input - Threshold                    
            DistanceFromThresholdSquared      = DistanceFromThreshold^2;                               

            if DistanceFromThreshold > 0  % Look only on false alarms (passing the threshold when not supposed to) 
                % For each noise STD, compute NoiseVar(fold-change) at each time point, extract SNR(t), compute and store mean(SNR(t))      
                FoldChangeInputWithoutNoiseMat    = repmat(FoldChangeInputWithoutNoise,iterations,length(Input));
                SNR_PerNoiseSTD_FoldChange        = NaN*ones(1,length(NoiseSTD_vec),'single');      

                for n_ind = 1:length(NoiseSTD_vec)
                    CurrentNoiseMatForRepeats = CurrentNoiseMatForRepeats_STD1 * CurrentNoiseSTD(n_ind);                                        
                    % Correct noise to prevent negative concentration values  (S+noise >= 0) 
                    CurrentNoiseMatForRepeats(CurrentNoiseMatForRepeats<=-InputMat2) = -0.9999*InputMat2(CurrentNoiseMatForRepeats<=-InputMat2);
                    CurrentNoiseDerivativeMatForRepeats     = [diff(CurrentNoiseMatForRepeats,[],2), NaN*ones(iterations,1)] / CurrentDeltaT;  

    %                 NoiseVariance = var( FoldChangeInputWithNoise - FoldChangeInputWithoutNoise );
                    FoldChangeInputWithNoise  = (-Current_a*InputMat2 + CurrentNoiseDerivativeMatForRepeats)./(InputMat2+CurrentNoiseMatForRepeats);  
                    WithMinusWithoutNoiseMat  = FoldChangeInputWithNoise - FoldChangeInputWithoutNoiseMat;
                    VarianceInTime            = var(WithMinusWithoutNoiseMat,[],1);  % ERGODIC assumption: run variance over repeats at a given time 't' gives variance(t)
                    SNR_PerNoiseSTD_FoldChange(n_ind) = nanmean(DistanceFromThresholdSquared ./ VarianceInTime); % mean(SNR(t))
                end                              
                SNR.ExponentInput.FoldChange(S0_ind,SR_ind,a_ind,:) = SNR_PerNoiseSTD_FoldChange;                            
           end                                  
        end                                
    end        
end
save(NoiseAnalysis_part1_FileName,'*','-v7.3');

% load(NoiseAnalysis_part2_FileName);     % Intermediate file for part 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Part 2: Calculate the effect of integration time on SNR in the derivative and fold change models (low pass filter) and in the ACT model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%  Input structure:  SNR.(InputType).(ModelType)  from Part 1
%  Output structure: SNR.(InputType).(ModelType)  , where ModelType = Absolute, Derivative, FoldChange, ACT
%     Same as part 1:
%         SNR.ConstantInput.(ModelType) = 3-D matrix (Input amplitude, Sample rate, noise level)
%         SNR.LinearInput.(ModelType)   = 4-D matrix (Input amplitude, Sample rate, decay rate(b_vec), noise level)
%         SNR.ExponentInput.(ModelType) = 4-D matrix (Input amplitude, Sample rate, decay rate(a_vec), noise level)
%     New for part 2:
%         SNR.ConstantInputWithLPFOnInput.(ModelType)      = 4-D matrix for Derivative or fold change models (Input amplitude, Sample rate, noise level, Integration time(Tao_vec) )
%                                                            The low-pass filter is applied on input before differentiation.  
%         SNR.ConstantInputWithLPFOnVariable.(ModelType)   = 4-D matrix for Derivative or fold change models (Input amplitude, Sample rate, noise level, Integration time(Tao_vec) )
%                                                            The low-pass filter is applied after differentiation, on the relevant variable.  
%         SNR.ConstantInput_VaryTao.ACT                    = 4-D matrix for ACT model                        (Input amplitude, Sample rate, noise level, Integration time(Tao_vec) )
%  Structure is saved in 'NoiseAnalysis_part2_FileName.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%%%  Free Parameters  %%%
SampleRate_TimeVec2 = 500; % Hz   % We will skip 5kHz due to long calculation times 
TotalTime2          = 200;  % seconds
Time2               = single(0:1/SampleRate_TimeVec2:TotalTime2);
length_Time2        = length(Time2);

% Thresholds information - use previously extracted parameter information  
Tderivative = -0.28; % microM/sec
Tfoldchange = -0.06; % 1/sec
K_ACT       = 5.5;   % microM

%%%  Vectors For Screen  %%%
% MinimalWormSampleRate = 5;          % Hz. Sample rate cannot be slower than our observed response times - around 200 msec
SampleRate_vec = [5000 500 50 5];     % 1/DeltaT      % Here we skip 5kHz due to long running time 
S0_vec         = [1 10 1e2 1e3 1e4];  % in microM.      For constant, linear decay and exponential decay 
NoiseSTD_vec   = 10.^(-5:0.05:6);     % in microM
Tao_vec        = [0.01:0.002:0.03 0.035:0.005:0.05 0.06:0.01:0.12 0.14:0.02:0.3 0.35:0.05:0.5 0.6:0.1:1.2 1.4:0.2:3 3.5:0.5:5 6:1:12 14:2:30 35:5:50]; % IntegrationTime_vec % sec

rng('default')
NoiseVec_STD1 = randn(1,length_Time2,'single'); 

SNR.ConstantInputWithLPFOnInput.Derivative    = NaN* ones(length(S0_vec), length(SampleRate_vec), length(NoiseSTD_vec), length(Tao_vec),'single');
SNR.ConstantInputWithLPFOnInput.FoldChange    = NaN* ones(length(S0_vec), length(SampleRate_vec), length(NoiseSTD_vec), length(Tao_vec),'single');
SNR.ConstantInputWithLPFOnVariable.Derivative = NaN* ones(length(S0_vec), length(SampleRate_vec), length(NoiseSTD_vec), length(Tao_vec),'single');
SNR.ConstantInputWithLPFOnVariable.FoldChange = NaN* ones(length(S0_vec), length(SampleRate_vec), length(NoiseSTD_vec), length(Tao_vec),'single');
SNR.ConstantInput_VaryTao.ACT                 = NaN* ones(length(S0_vec), length(SampleRate_vec), length(NoiseSTD_vec), length(Tao_vec),'single');
tic
a = 1;  % for filter
for Tao_ind = 1:length(Tao_vec)  % integration times
    Current_tao = Tao_vec(Tao_ind);
    disp([datestr(now),'   --   Computing SNR for constant signal with LPF of ',num2str(Current_tao),' seconds']);

    for S0_ind = 1:length(S0_vec)
        S0 = S0_vec(S0_ind);    
        for SR_ind = 1:length(SampleRate_vec)
            if SR_ind==1  % We will skip 5kHz due to long running times 
                continue
            end
            CurrentSampleRate    = SampleRate_vec(SR_ind);
            CurrentDeltaT        = 1/CurrentSampleRate;
            IntervalsInTimeVec   = round(SampleRate_TimeVec2 / CurrentSampleRate);
            CurrentNoiseVec_STD1 = NoiseVec_STD1(:,1:IntervalsInTimeVec:end);     % noise (unfiltered) (STD1)
            CurrentTime          = Time2(1:IntervalsInTimeVec:end);
            
            % Correct noise to prevent negative concentration values  (S+noise >= 0) 
            CurrentNoiseMat                      = NoiseSTD_vec' * CurrentNoiseVec_STD1;      % noise (unfiltered) (All STDs)  
            CurrentNoiseMat(CurrentNoiseMat<=-S0)= -0.9999*S0;
            CurrentNoiseDerivativeMat            = [diff(CurrentNoiseMat,[],2), NaN*ones(length(NoiseSTD_vec),1)] / CurrentDeltaT;   % noise (unfiltered) derivative (All STDs)        
            CurrentNoiseSTD                      = nanstd(CurrentNoiseMat,[],2);

            %%% ACT                           
            InputMat      = S0+CurrentNoiseMat; InputMat(:,1)=S0;
            for n_ind = 1:length(CurrentNoiseSTD)
                Input = InputMat(n_ind,:);
                CurrentACT = CalculateAdaptiveThreshold_inline (CurrentTime, Input, K_ACT, 1/Current_tao); % ACT, full vector
                DistanceFromThreshold_OverTime                  = S0 - CurrentACT;  % Signal = Input(without noise) - Threshold
                DistanceFromThreshold_OverTime(DistanceFromThreshold_OverTime<0)= 0;   % make sure signal is zero whenever a "miss" is identified          
                AverageDistanceSquareFromThreshold              = nanmean(DistanceFromThreshold_OverTime.^2);
                SNR_PerNoiseSTD_ACT                             = AverageDistanceSquareFromThreshold ./ CurrentNoiseSTD(n_ind).^2;  % Use non-filtered noise for ACT model
                SNR.ConstantInput_VaryTao.ACT(S0_ind,SR_ind,n_ind,Tao_ind) = SNR_PerNoiseSTD_ACT;    
            end
                                   
            % Find filter parameters
            CurrentTimeScale   = Current_tao * CurrentSampleRate;      % Relevant timescale for filtering at the current sampling rate         
            CurrentWindowSize  = max([round(CurrentTimeScale * 3) 2]); % at least 2 points
            Current_b          = exp(-(0:CurrentWindowSize-1)/CurrentTimeScale);  Current_b=Current_b/sum(Current_b);                           
            
            %%% filtering noise %%%
            CurrentNoiseMat_AfterLPF                          = filter(Current_b,a,CurrentNoiseMat,[],2);
            CurrentNoiseMat_AfterLPF(:,1:CurrentWindowSize-1) = NaN;                           
            
            %%% filtering noise derivative %%%
            CurrentNoiseDerivativeMat_LPFonDerivative                          = filter(Current_b,a,CurrentNoiseDerivativeMat,[],2);
            CurrentNoiseDerivativeMat_LPFonDerivative(:,1:CurrentWindowSize-1) = NaN;                           
            NoiseDerivativeSTDvec_LPFonDerivative                              = nanstd(CurrentNoiseDerivativeMat_LPFonDerivative,[],2);    % LPF(Derivative(noise))  
            
            %%% derivative of filtered noise (As control, mathematically identical to the definition above)              
            CurrentNoiseDerivativeMat_DerivativeOnLPF = [diff(CurrentNoiseMat_AfterLPF,[],2), NaN*ones(length(NoiseSTD_vec),1)] / CurrentDeltaT; 
            NoiseDerivativeSTDvec_DerivativeOnLPF     = nanstd(CurrentNoiseDerivativeMat_DerivativeOnLPF,[],2);
                        
            %% Derivative SNR
            %%% Derivative. CurrentNoise = Derivative(LPF(noise))
            DistanceFromThreshold                         = 0 - Tderivative;
            DistanceFromThreshold(DistanceFromThreshold<0)= NaN; % Look only on false alarms (passing the threshold when not supposed to)          
            SNR_PerNoiseSTD_Derivative                    = (DistanceFromThreshold)^2 ./ NoiseDerivativeSTDvec_DerivativeOnLPF.^2;
            SNR.ConstantInputWithLPFOnInput.Derivative(S0_ind,SR_ind,:,Tao_ind) = SNR_PerNoiseSTD_Derivative;        
                        
            %%% Derivative. CurrentNoise = LPF(derivative(noise)). That's a Control, should be mathematically identical to the one above   
            SNR_PerNoiseSTD_Derivative                     = (DistanceFromThreshold)^2 ./ NoiseDerivativeSTDvec_LPFonDerivative.^2;
            SNR.ConstantInputWithLPFOnVariable.Derivative(S0_ind,SR_ind,:,Tao_ind) = SNR_PerNoiseSTD_Derivative;                                         
                        
            %% Fold-change SNR.   For S=constant: NoiseVar(fold-change) = var((d(WGN)/dt))/(S+WGN))
            DistanceFromThreshold                         = 0 - Tfoldchange;
            DistanceFromThreshold(DistanceFromThreshold<0)= NaN;              % Look only on false alarms (passing the threshold when not supposed to)     
                         
            if ~isnan(DistanceFromThreshold) % skip these time-limiting calculations if they are not necessary
                %%% Fold-Change with LPF on original signal,    CurrentNoise = Derivative(LPF(noise)) / (S0+ LPF(noise))                                     
                FoldChangeNoiseMat         = CurrentNoiseDerivativeMat_DerivativeOnLPF ./ (S0+CurrentNoiseMat_AfterLPF);
                NoiseVar_FoldChange        = nanvar(FoldChangeNoiseMat,[],2);                
                SNR_PerNoiseSTD_FoldChange = (DistanceFromThreshold)^2 ./ NoiseVar_FoldChange;                       
                SNR.ConstantInputWithLPFOnInput.FoldChange(S0_ind,SR_ind,:,Tao_ind) = SNR_PerNoiseSTD_FoldChange;                          

                %%% Fold-Change with LPF on signal fold-change, CurrentNoise = LPF( Derivative(noise) / (S0+ noise) ) 
                FoldChangeNoiseMat_BeforeFilter             = CurrentNoiseDerivativeMat ./ (S0+CurrentNoiseMat);
                FoldChangeNoiseMat                          = filter(Current_b,a,FoldChangeNoiseMat_BeforeFilter,[],2);                
                FoldChangeNoiseMat(:,1:CurrentWindowSize-1) = NaN;                            
                NoiseVar_FoldChange                         = nanvar(FoldChangeNoiseMat,[],2);                
                SNR_PerNoiseSTD_FoldChange                  = (DistanceFromThreshold)^2 ./ NoiseVar_FoldChange;                       
                SNR.ConstantInputWithLPFOnVariable.FoldChange(S0_ind,SR_ind,:,Tao_ind) = SNR_PerNoiseSTD_FoldChange;                                  
            end            
        end        
    end
    
end
save(NoiseAnalysis_part2_FileName,'*','-v7.3');
% load(NoiseAnalysis_part2_FileName);    % Intermediate file for part 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Part 3: Calculate the effect of integration time on response speed in the derivative and fold change models (low pass filter) and in the ACT model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%  How does integration effect the responses to exponential decaying inputs?  Input 3: S3 = S0*exp(-at)  
%  Output structures:
%    Latency, DelayTime, Speed    >>>> Outputs of screen over 75 integraton times (Tao_vec),   6 decay speeds (a_vec1), 5 initial input amplitudes (S0_vec) 
%    Latency2, DelayTime2, Speed2 >>>> Outputs of screen over 4  integraton times (Tao_vec2), 91 decay speeds (a_vec2), 1 initial input amplitudes (S0_vec2) 
%    Second fields:
%        Derivative_LPFOnInput       >> Low-pass derivative model  with integration applied on input before differentiation  
%        FoldChange_LPFOnInput       >> Low-pass fold-change model with integration applied on input before differentiation 
%        Derivative_LPFOnVariable    >> Low-pass derivative model  with integration applied on derivative 
%        FoldChange_LPFOnVariable    >> Low-pass fold-change model with integration applied on fold-change 
%        ACT                         >> ACT model with different integration time (neuron adaptation time)  
%        ReferenceLatency_Derivative >> Derivative model without input integration (no low-pass filter) 
%        ReferenceLatency_FoldChange >> Fold-change model without input integration (no low-pass filter) 
%        ReferenceLatency_ACTmodel   >> ACT model without input integration (threshold is not adapting) 
%  Structure is saved in 'NoiseAnalysis_part2_FileName.mat'
%  Each screen will be done in 3 steps:
%       Step 1. Find response latency as a function of integration time
%       Step 2. Define minimal latency: derivative and fold-change at tao = 0, ACT model at tao = inf   
%       Step 3. ResponseSpeed = 1/DelayTime, where DelayTime = latency - minimal latency.    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Initialization
%%%  Free Parameters  %%%
SampleRate_TimeVec3 = 5000; % Hz   % We will skip 5kHz due to long running times 
TotalTime3          = 150;  % seconds
Time3               = single(-TotalTime3:1/SampleRate_TimeVec3:TotalTime3);
ZeroIndex           = TotalTime3*SampleRate_TimeVec3 + 1;
Time3_FromZero      = Time3(ZeroIndex:end); 

% Thresholds information - use previously extracted parameter information  
Tderivative = -0.28; % microM/sec
Tfoldchange = -0.06; % 1/sec
K_ACT       = 5.5;   % microM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First Loop, over 75 integraton times, 6 decay speeds, 5 initial amplitudes
%%%  Vectors For Screen  %%%
S0_vec         = [1 10 1e2 1e3 1e4];    % in microM.     For constant, linear decay and exponential decay 
a_vec1         = [0.3 0.5 1 5 10 100];  % in 1/sec.      For exponential decay

%%% Initialization
Latency.Derivative_LPFOnInput    = zeros(length(Tao_vec),length(a_vec1),length(S0_vec),'single')*NaN; % Rows: 2 Types of integration calculation, Columns: integration times
Latency.Derivative_LPFOnVariable = zeros(length(Tao_vec),length(a_vec1),length(S0_vec),'single')*NaN; % Rows: 2 Types of integration calculation, Columns: integration times
Latency.FoldChange_LPFOnInput    = zeros(length(Tao_vec),length(a_vec1),length(S0_vec),'single')*NaN;
Latency.FoldChange_LPFOnVariable = zeros(length(Tao_vec),length(a_vec1),length(S0_vec),'single')*NaN;
Latency.ACT                      = zeros(length(Tao_vec),length(a_vec1),length(S0_vec),'single')*NaN;
Latency.ReferenceLatency_ACTmodel   =              zeros(length(a_vec1),length(S0_vec),'single')*NaN;
Latency.ReferenceLatency_Derivative =              zeros(length(a_vec1),length(S0_vec),'single')*NaN;
Latency.ReferenceLatency_FoldChange =              zeros(length(a_vec1),length(S0_vec),'single')*NaN;

%%% loop
a = 1;   % For filter
for Tao_ind = 1:length(Tao_vec)
    Tao = Tao_vec(Tao_ind);
    disp([datestr(now),'   --   Computing reponse latency for exponentially decaying signal with LPF of ',num2str(Tao),' seconds']);

    Tao        = Tao_vec(Tao_ind);          % in seconds
    TimeScale  = Tao*SampleRate_TimeVec3;   % in frames
    windowSize = round(TimeScale*3);
    b          = exp(-(0:windowSize-1)/TimeScale);  b=b/sum(b); 

    for a_ind = 1:length(a_vec1)
        Current_a                = a_vec1(a_ind);
        %%% Input
        InputAmp1                = ones(1,length(Time3),'single');       % S0=1
        InputAmp1(ZeroIndex:end) = exp(-Current_a*Time3(ZeroIndex:end));  % S0=1
        InputAmp1_FromZero       = InputAmp1(ZeroIndex:end);
        %%% Derivative
        InputAmp1_Derivative          = [NaN, diff(InputAmp1)]*SampleRate_TimeVec3;
        %%% filter(input)
        InputAmp1_AfterLPF          = filter(b,a,InputAmp1);   
        InputAmp1_AfterLPF(1:windowSize-1) = NaN;
        InputAmp1_AfterLPF_FromZero = InputAmp1_AfterLPF(ZeroIndex:end); 
        %%% filter(derivative(input))
        InputDerivativeAmp1_LPFonDerivative          = filter(b,a,InputAmp1_Derivative);       
        InputDerivativeAmp1_LPFonDerivative_FromZero = InputDerivativeAmp1_LPFonDerivative(ZeroIndex:end); 
        %%% derivative(filter(input))
        InputDerivativeAmp1_DerivativeOnLPF          = [NaN, diff(InputAmp1_AfterLPF)]*SampleRate_TimeVec3;
        InputDerivativeAmp1_DerivativeOnLPF_FromZero = InputDerivativeAmp1_DerivativeOnLPF(ZeroIndex:end);         
        %%% FoldChange(derivative(filter(input) / filter(input))
        InputFoldChange_LPFonInput_FromZero     = InputDerivativeAmp1_DerivativeOnLPF_FromZero ./ InputAmp1_AfterLPF_FromZero;
        %%% filter(FoldChange)        
        InputFoldChange                         = zeros(1,length(Time3),'single');      
        InputFoldChange(ZeroIndex:end)          = -Current_a;     
        InputFoldChange_LPFonVariable           = filter(b,a,InputFoldChange);  
        InputFoldChange_LPFonVariable_FromZero  = InputFoldChange_LPFonVariable(ZeroIndex:end);   
        
        %%% Latency for fold-change model (independent of S0 for exponential decaying signal)
        CurrentTime = Time3_FromZero(find(InputFoldChange_LPFonInput_FromZero <= Tfoldchange, 1, 'first'));
        if ~isempty(CurrentTime)
            Latency.FoldChange_LPFOnInput(Tao_ind,a_ind,:) = CurrentTime;
        end
        CurrentTime = Time3_FromZero(find(InputFoldChange_LPFonVariable_FromZero <= Tfoldchange, 1, 'first'));
        if ~isempty(CurrentTime)
            Latency.FoldChange_LPFOnVariable(Tao_ind,a_ind,:) = CurrentTime;
        end                
        
        %%% Latency for derivative and ACT models (S0 dependent)
        for S0_ind = 1:length(S0_vec)
            S0 = S0_vec(S0_ind);            
            Input_FromZero                             = S0 * InputAmp1_FromZero;   % Not filtered. For ACT
            InputDerivative_LPFonDerivative_FromZero   = S0 * InputDerivativeAmp1_LPFonDerivative_FromZero;    % Filtered. For derivative, option 1         
            InputDerivative_DerivativeOnLPF_FromZero   = S0 * InputDerivativeAmp1_DerivativeOnLPF_FromZero;    % Filtered. For derivative, option 2, should be identical to option 1    
                      
            %%% Latency for derivative model            
            CurrentTime = Time3_FromZero(find(InputDerivative_DerivativeOnLPF_FromZero <= Tderivative, 1, 'first'));
            if ~isempty(CurrentTime)
                Latency.Derivative_LPFOnInput(Tao_ind,a_ind,S0_ind) = CurrentTime;      % derivative(filter(input))
            end
            CurrentTime = Time3_FromZero(find(InputDerivative_LPFonDerivative_FromZero <= Tderivative, 1, 'first'));
            if ~isempty(CurrentTime)
                Latency.Derivative_LPFOnVariable(Tao_ind,a_ind,S0_ind) = CurrentTime;   % filter(derivative(input)) 
            end                                  

            %%%  Latency for ACT model
            Beta = 1/Tao;                        
            [~, ~, CurrentTime] = CalculateAdaptiveThreshold_ExpSaturation_BeforeActivation(Time3_FromZero, Input_FromZero, K_ACT, Beta);
            if ~isnan(CurrentTime)
                Latency.ACT(Tao_ind,a_ind,S0_ind) = CurrentTime;
            end            
        end
    end   
end

%%% Reference latencies 
for a_ind = 1:length(a_vec1)
    Current_a                = a_vec1(a_ind);
    %%% Input
    InputAmp1                = zeros(1,length(Time3),'single');       % S0=1
    InputAmp1(ZeroIndex:end) = exp(-Current_a*Time3(ZeroIndex:end));  % S0=1
    InputAmp1_FromZero       = InputAmp1(ZeroIndex:end);
    %%% Derivative
    InputAmp1_Derivative_FromZero = -Current_a * InputAmp1_FromZero;
    InputFoldChange_FromZero      = -Current_a * ones(1,length(Time3_FromZero),'single');      % Not Filtered. For fold-change
    
    %%% Latency for fold-change model (independent of S0 for exponential decaying signal)
    CurrentTime = Time3_FromZero(find(InputFoldChange_FromZero <= Tfoldchange, 1, 'first'));
    if ~isempty(CurrentTime)
        Latency.ReferenceLatency_FoldChange(a_ind,:) = CurrentTime;
    end
    
    %%% Latency for derivative and ACT models (S0 dependent)
    for S0_ind = 1:length(S0_vec)
        S0 = S0_vec(S0_ind);            
        Input_FromZero           = S0 * InputAmp1_FromZero;              % Not filtered. For ACT
        InputDerivative_FromZero = S0 * InputAmp1_Derivative_FromZero;   % Not filtered. For Derivatice

        %%% Latency for derivative model            
        CurrentTime = Time3_FromZero(find(InputDerivative_FromZero <= Tderivative, 1, 'first'));
        if ~isempty(CurrentTime)
            Latency.ReferenceLatency_Derivative(a_ind,S0_ind) = CurrentTime;  
        end

        %%% Latency for ACT model           
        ConstantThreshold = K*(1-exp(-S0/K));  % no adaptation
        CurrentTime       = Time3_FromZero(find(Input_FromZero <= ConstantThreshold, 1, 'first'));
        if ~isempty(CurrentTime)
            Latency.ReferenceLatency_ACTmodel(a_ind,S0_ind) = CurrentTime;  
        end                                                    
    end
end   

%%% Compute DelayTimes and Speed structures
FieldsInStructure          = {'Derivative_LPFOnInput',      'Derivative_LPFOnVariable',   'FoldChange_LPFOnInput',      'FoldChange_LPFOnVariable',   'ACT'};
ReferenceFieldsInStructure = {'ReferenceLatency_Derivative','ReferenceLatency_Derivative','ReferenceLatency_FoldChange','ReferenceLatency_FoldChange','ReferenceLatency_ACTmodel'};
for f_ind = 1:length(FieldsInStructure)
    CurrentField             = FieldsInStructure{f_ind};
    CurrentReferenceField    = ReferenceFieldsInStructure{f_ind};
    DelayTime.(CurrentField) = Latency.(CurrentField)*NaN;         % Initialization 
    for Tao_ind = 1:length(Tao_vec)
        DelayTime.(CurrentField)(Tao_ind,:,:) = squeeze(Latency.(CurrentField)(Tao_ind,:,:)) - Latency.(CurrentReferenceField);        
    end    
    Speed.(CurrentField) = 1 ./ DelayTime.(CurrentField);       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second Loop, over 4 integraton times, 91 decay speeds, 1 initial amplitude
%%%  Vectors For Screen  %%%
S0_vec2         = [10];    % in microM.     For constant, linear decay and exponential decay 
a_vec2          = [0.3:0.02:0.78 0.8:0.05:2 2.1:0.1:4 4.2:0.2:7 7.5:0.5:10];  % in 1/sec.      For exponential decay
Tao_vec2        = [0.1 1 10 50];         % in sec

%%% Initialization
Latency2.Derivative_LPFOnInput    = zeros(length(Tao_vec2),length(a_vec2),length(S0_vec2),'single')*NaN; % Rows: 2 Types of integration calculation, Columns: integration times
Latency2.Derivative_LPFOnVariable = zeros(length(Tao_vec2),length(a_vec2),length(S0_vec2),'single')*NaN; % Rows: 2 Types of integration calculation, Columns: integration times
Latency2.FoldChange_LPFOnInput    = zeros(length(Tao_vec2),length(a_vec2),length(S0_vec2),'single')*NaN;
Latency2.FoldChange_LPFOnVariable = zeros(length(Tao_vec2),length(a_vec2),length(S0_vec2),'single')*NaN;
Latency2.ACT                      = zeros(length(Tao_vec2),length(a_vec2),length(S0_vec2),'single')*NaN;
Latency2.ReferenceLatency_ACTmodel   =               zeros(length(a_vec2),length(S0_vec2),'single')*NaN;
Latency2.ReferenceLatency_Derivative =               zeros(length(a_vec2),length(S0_vec2),'single')*NaN;
Latency2.ReferenceLatency_FoldChange =               zeros(length(a_vec2),length(S0_vec2),'single')*NaN;

%%% loop
a = 1;   % For filter
for Tao_ind = 1:length(Tao_vec2)
    Tao = Tao_vec2(Tao_ind);
    disp([char(10),char(10),char(10),datestr(now),'   --   Computing reponse latency for exponentially decaying signal with LPF of ',num2str(Tao),' seconds']);

    Tao        = Tao_vec2(Tao_ind);         % in seconds
    TimeScale  = Tao*SampleRate_TimeVec3;   % in frames
    windowSize = round(TimeScale*3);
    b          = exp(-(0:windowSize-1)/TimeScale);  b=b/sum(b); 
    % figure; plot((0:windowSize-1)/SampleRate_TimeVec3, b); title(['filter of ',num2str(TimeScale),' frames'])    

    for a_ind = 1:length(a_vec2)
        Current_a                = a_vec2(a_ind);
        disp([datestr(now),'   --   a_ind = ',num2str(a_ind),' / ',num2str(length(a_vec2))]);
        %%% Input
        InputAmp1                = ones(1,length(Time3),'single');       % S0=1
        InputAmp1(ZeroIndex:end) = exp(-Current_a*Time3(ZeroIndex:end));  % S0=1
        InputAmp1_FromZero       = InputAmp1(ZeroIndex:end);
        %%% Derivative
        InputAmp1_Derivative          = [NaN, diff(InputAmp1)]*SampleRate_TimeVec3;
        %%% filter(input)
        InputAmp1_AfterLPF          = filter(b,a,InputAmp1);   
        InputAmp1_AfterLPF(1:windowSize-1) = NaN;
        InputAmp1_AfterLPF_FromZero = InputAmp1_AfterLPF(ZeroIndex:end); 
        %%% filter(derivative(input))
        InputDerivativeAmp1_LPFonDerivative          = filter(b,a,InputAmp1_Derivative);       
        InputDerivativeAmp1_LPFonDerivative_FromZero = InputDerivativeAmp1_LPFonDerivative(ZeroIndex:end); 
        %%% derivative(filter(input))
        InputDerivativeAmp1_DerivativeOnLPF          = [NaN, diff(InputAmp1_AfterLPF)]*SampleRate_TimeVec3;
        InputDerivativeAmp1_DerivativeOnLPF_FromZero = InputDerivativeAmp1_DerivativeOnLPF(ZeroIndex:end);         
        %%% FoldChange(derivative(filter(input) / filter(input))
        InputFoldChange_LPFonInput_FromZero     = InputDerivativeAmp1_DerivativeOnLPF_FromZero ./ InputAmp1_AfterLPF_FromZero;
        %%% filter(FoldChange)        
        InputFoldChange                         = zeros(1,length(Time3),'single');      
        InputFoldChange(ZeroIndex:end)          = -Current_a;     
        InputFoldChange_LPFonVariable           = filter(b,a,InputFoldChange);  
        InputFoldChange_LPFonVariable_FromZero  = InputFoldChange_LPFonVariable(ZeroIndex:end);   
        
        %%% Latency for fold-change model (independent of S0 for exponential decaying signal)
        CurrentTime = Time3_FromZero(find(InputFoldChange_LPFonInput_FromZero <= Tfoldchange, 1, 'first'));
        if ~isempty(CurrentTime)
            Latency2.FoldChange_LPFOnInput(Tao_ind,a_ind,:) = CurrentTime;
        end
        CurrentTime = Time3_FromZero(find(InputFoldChange_LPFonVariable_FromZero <= Tfoldchange, 1, 'first'));
        if ~isempty(CurrentTime)
            Latency2.FoldChange_LPFOnVariable(Tao_ind,a_ind,:) = CurrentTime;
        end                
        
        %%% Latency for derivative and ACT models (S0 dependent)
        for S0_ind = 1:length(S0_vec2)
            S0 = S0_vec2(S0_ind);            
            Input_FromZero                             = S0 * InputAmp1_FromZero;   % Not filtered. For ACT
            InputDerivative_LPFonDerivative_FromZero   = S0 * InputDerivativeAmp1_LPFonDerivative_FromZero;    % Filtered. For derivative, option 1         
            InputDerivative_DerivativeOnLPF_FromZero   = S0 * InputDerivativeAmp1_DerivativeOnLPF_FromZero;    % Filtered. For derivative, option 2, should be identical to option 1    
                      
            %%% Latency for derivative model            
            CurrentTime = Time3_FromZero(find(InputDerivative_DerivativeOnLPF_FromZero <= Tderivative, 1, 'first'));
            if ~isempty(CurrentTime)
                Latency2.Derivative_LPFOnInput(Tao_ind,a_ind,S0_ind) = CurrentTime;      % derivative(filter(input))
            end
            CurrentTime = Time3_FromZero(find(InputDerivative_LPFonDerivative_FromZero <= Tderivative, 1, 'first'));
            if ~isempty(CurrentTime)
                Latency2.Derivative_LPFOnVariable(Tao_ind,a_ind,S0_ind) = CurrentTime;   % filter(derivative(input)) 
            end                                  

            %%%  Latency for ACT model
            Beta = 1/Tao;                        
            [~, ~, CurrentTime] = CalculateAdaptiveThreshold_ExpSaturation_BeforeActivation(Time3_FromZero, Input_FromZero, K_ACT, Beta);
            if ~isnan(CurrentTime)
                Latency2.ACT(Tao_ind,a_ind,S0_ind) = CurrentTime;
            end            
        end
    end   
end
toc

% Reference latencies (running time < 1 sec)
for a_ind = 1:length(a_vec2)
    Current_a                = a_vec2(a_ind);
    %%% Input
    InputAmp1                = zeros(1,length(Time3),'single');       % S0=1
    InputAmp1(ZeroIndex:end) = exp(-Current_a*Time3(ZeroIndex:end));  % S0=1
    InputAmp1_FromZero       = InputAmp1(ZeroIndex:end);
    %%% Derivative
    InputAmp1_Derivative_FromZero = -Current_a * InputAmp1_FromZero;
    InputFoldChange_FromZero      = -Current_a * ones(1,length(Time3_FromZero),'single');      % Not Filtered. For fold-change
    
    %%% Latency for fold-change model (independent of S0 for exponential decaying signal)
    CurrentTime = Time3_FromZero(find(InputFoldChange_FromZero <= Tfoldchange, 1, 'first'));
    if ~isempty(CurrentTime)
        Latency2.ReferenceLatency_FoldChange(a_ind,:) = CurrentTime;
    end
    
    %%% Latency for derivative and ACT models (S0 dependent)
    for S0_ind = 1:length(S0_vec2)
        S0 = S0_vec2(S0_ind);            
        Input_FromZero           = S0 * InputAmp1_FromZero;              % Not filtered. For ACT
        InputDerivative_FromZero = S0 * InputAmp1_Derivative_FromZero;   % Not filtered. For Derivatice

        %%% Latency for derivative model            
        CurrentTime = Time3_FromZero(find(InputDerivative_FromZero <= Tderivative, 1, 'first'));
        if ~isempty(CurrentTime)
            Latency2.ReferenceLatency_Derivative(a_ind,S0_ind) = CurrentTime;  
        end

        %%% Latency for ACT model           
        ConstantThreshold = K*(1-exp(-S0/K));  % no adaptation
        CurrentTime       = Time3_FromZero(find(Input_FromZero <= ConstantThreshold, 1, 'first'));
        if ~isempty(CurrentTime)
            Latency2.ReferenceLatency_ACTmodel(a_ind,S0_ind) = CurrentTime;  
        end                                                    
    end
end   

% Compute DelayTimes and Speed structures
FieldsInStructure          = {'Derivative_LPFOnInput',      'Derivative_LPFOnVariable',   'FoldChange_LPFOnInput',      'FoldChange_LPFOnVariable',   'ACT'};
ReferenceFieldsInStructure = {'ReferenceLatency_Derivative','ReferenceLatency_Derivative','ReferenceLatency_FoldChange','ReferenceLatency_FoldChange','ReferenceLatency_ACTmodel'};
for f_ind = 1:length(FieldsInStructure)
    CurrentField             = FieldsInStructure{f_ind};
    CurrentReferenceField    = ReferenceFieldsInStructure{f_ind};
    DelayTime2.(CurrentField) = Latency2.(CurrentField)*NaN;         % Initialization 
    for Tao_ind = 1:length(Tao_vec2)
        DelayTime2.(CurrentField)(Tao_ind,:) = squeeze(Latency2.(CurrentField)(Tao_ind,:)) - Latency2.(CurrentReferenceField)';        
    end    
    Speed2.(CurrentField) = 1 ./ DelayTime2.(CurrentField);       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Save all structures (from all parts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
save(NoiseAnalysis_FileName,'*','-v7.3');

% load(NoiseAnalysis_FileName);
return


function [AdaptiveThresholdVector, PredictedActivationIndex, PredictedActivationTime] = ...
             CalculateAdaptiveThreshold_inline (TimeVector, DyeVector, K, Beta)

TimeDiffPerStep = diff(TimeVector([1 2])); % seconds
C0              = DyeVector(1);
TimeZeroIndex   = find(TimeVector>=0,1,'first');
VectorLength    = length(TimeVector);
if VectorLength ~= length(DyeVector)
    disp('Non-matching vector lengths. Aborting...')
    return
elseif size(TimeVector,1) == size(DyeVector,2) % row versus column vector -- fix it
    DyeVector = DyeVector';
end

AdaptiveThresholdVector = K*(1-exp(-C0/K))*ones(size(TimeVector,1),size(TimeVector,2));
for ind = (TimeZeroIndex+1):VectorLength
    CurrentDyeValue              = DyeVector(ind);
    AdaptiveThresholdVector(ind) = (AdaptiveThresholdVector(ind-1)*exp(-Beta*TimeDiffPerStep) + K*(1-exp(-CurrentDyeValue/K))*(1- exp(-Beta*TimeDiffPerStep)));
end 
PredictedActivationIndex = find(AdaptiveThresholdVector>DyeVector,1,'first'); % Threhold curve cross dye curve
PredictedActivationTime  = TimeVector(PredictedActivationIndex);

return



