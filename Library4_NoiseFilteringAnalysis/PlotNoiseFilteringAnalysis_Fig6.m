function  PlotNoiseFilteringAnalysis_Fig6
%
%  Written by Sagi Levy, September 2019
%

NoiseAnalysis_FileName = 'D:\Models\NoiseAnalysis_Figure6.mat';
load(NoiseAnalysis_FileName);

% SNR display limits
LimitForDisplay.Upper = 10;     
LimitForDisplay.Lower = 1/10;  

% LinesPerModel = {'k-','r-','k--'};
LinesPerModel = {'k-','r-','b-'};
S0_vec_Titles = {'10^0','10^1','10^2','10^3','10^4'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Plots for definitions and examples 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% SNR definition: Figure 6A
%%%  Free Parameters  %%%
SampleRate_TimeVec_Example  = 5; % Hz
TotalTime_Example    = 220;  % seconds
CurrentTime          = single(0:1/SampleRate_TimeVec_Example:TotalTime_Example);
S0                   = 10; % microM
CurrentDeltaT        = 1/SampleRate_TimeVec_Example;

%%% Thresholds information - use previously extracted parameter information  %%%  
Tderivative = -0.28; % microM/sec
Tfoldchange = -0.06; % 1/sec
K           = 5.5;   % microM
Tao         = 17;    % seconds

%%% High noise values %%%
NoiseSTD = 1;  % microM 
rng('default')
NoiseVec                = randn(1,length(CurrentTime),'single')*NoiseSTD; 
NoiseVec(NoiseVec<=-S0) = -0.9999*S0;  % Correct noise to prevent negative concentration values  (S+noise >= 0) 
NoiseDerivativeVec      = [diff(NoiseVec,[],2), NaN] / CurrentDeltaT;             
    
InputWithNoise   = S0+NoiseVec; InputWithNoise(1)=S0;
CurrentACT       = CalculateAdaptiveThreshold_inline (CurrentTime, InputWithNoise, K, 1/Tao); % ACT, full vector
FoldChangeNoise  = NoiseDerivativeVec ./ (S0 + NoiseVec); % compute  NoiseVar(fold-change), and for Si=constant -->  var((d(WGN)/dt))/(Si+WGN))

figure('name','SNR definitions- high noise values','position',[212   144   292   822]);
subplot(3,1,1)
    plot(CurrentTime, S0*ones(1,length(CurrentTime)),'k-'); hold on;
    plot(CurrentTime, CurrentACT,'b-'); hold on;
    plot(CurrentTime, S0 + NoiseVec,'g-'); hold on; xlim([0 30]);
    ylim([0 15]); ylabel('C')
subplot(3,1,2)
    plot(CurrentTime, zeros(1,length(CurrentTime)),'k-'); hold on;
    plot(CurrentTime, Tderivative*ones(1,length(CurrentTime)),'b-'); hold on;
    plot(CurrentTime, NoiseDerivativeVec,'g-'); hold on; xlim([0 30]);
    ylim([-30 30]); ylabel('dC/dt')
subplot(3,1,3)
    plot(CurrentTime, zeros(1,length(CurrentTime)),'k-'); hold on;
    plot(CurrentTime, Tfoldchange*ones(1,length(CurrentTime)),'b-'); hold on;
    plot(CurrentTime, FoldChangeNoise,'g-'); hold on; xlim([0 30]);
    ylim([-4 4]);    ylabel('(dC/dt)/C')
    xlabel('Time [sec]')
          
%%% Low noise values
NoiseSTD = 0.01;  % microM 
rng('default')
NoiseVec                  = randn(1,length(CurrentTime),'single')*NoiseSTD; 
NoiseVec(NoiseVec<=-S0)   = -0.9999*S0;    % Correct noise to prevent negative concentration values  (S+noise >= 0) 
NoiseDerivativeVec        = [diff(NoiseVec,[],2), NaN] / CurrentDeltaT;             

InputWithNoise  = S0+NoiseVec; InputWithNoise(1)=S0;
CurrentACT      = CalculateAdaptiveThreshold_inline (CurrentTime, InputWithNoise, K, 1/Tao); % ACT, full vector
FoldChangeNoise = NoiseDerivativeVec ./ (S0 + NoiseVec); % compute  NoiseVar(fold-change), and for Si=constant -->  var((d(WGN)/dt))/(Si+WGN))

figure('name','SNR definitions- low noise values','position',[212   144   292   822]);
subplot(3,1,1)
    plot(CurrentTime, S0*ones(1,length(CurrentTime)),'k-'); hold on;
    plot(CurrentTime, CurrentACT,'b-'); hold on;
    plot(CurrentTime, S0 + NoiseVec,'g-'); hold on; xlim([0 30]);
    ylim([0 15]); ylabel('C')
subplot(3,1,2)
    plot(CurrentTime, zeros(1,length(CurrentTime)),'k-'); hold on;
    plot(CurrentTime, Tderivative*ones(1,length(CurrentTime)),'b-'); hold on;
    plot(CurrentTime, NoiseDerivativeVec,'g-'); hold on; xlim([0 30]);
    ylim([-0.6 0.3]); ylabel('dC/dt')
subplot(3,1,3)
    plot(CurrentTime, zeros(1,length(CurrentTime)),'k-'); hold on;
    plot(CurrentTime, Tfoldchange*ones(1,length(CurrentTime)),'b-'); hold on;
    plot(CurrentTime, FoldChangeNoise,'g-'); hold on; xlim([0 30]);
    ylim([-0.12 0.06]); ylabel('(dC/dt)/C')

%% Speed definition, Figure 6D    
%%%  Free Parameters  %%%
SampleRate_TimeVec_Example  = 50; % Hz
TotalTime_Example    = 220;  % seconds
ZeroTimeInSeconds    = 160;  % Seconds, i.e. time vector [-160,60] seconds   
ZeroIndex            = ZeroTimeInSeconds * SampleRate_TimeVec_Example + 1;
CurrentTime          = single(0:1/SampleRate_TimeVec_Example:TotalTime_Example) - ZeroTimeInSeconds;

%%%  Thresholds information - use previously extracted parameter information  %%%  
Tderivative = -0.28; % microM/sec
Tfoldchange = -0.06; % 1/sec
K           = 5.5;   % microM
Tao         = 20;    % seconds

% Input information
S0        = 10; % microM
Current_a = 0.3;

%%% LPF     
CurrentTimeScale   = Tao * SampleRate_TimeVec_Example;      % Relevant timescale for filtering at the current sampling rate         
CurrentWindowSize  = max([round(CurrentTimeScale * 3) 2]); % at least 2 points
Current_b          = exp(-(0:CurrentWindowSize-1)/CurrentTimeScale);  Current_b=Current_b/sum(Current_b);               

%%% Input
Input                = S0*ones(1,length(CurrentTime),'single');         
Input(ZeroIndex:end) = S0*exp(-Current_a*CurrentTime(ZeroIndex:end)); 
%%% Derivative
Input_Derivative     = [NaN, diff(Input)]*SampleRate_TimeVec_Example;
%%% filter(derivative(input))
InputDerivative_LPFonDerivative          = filter(Current_b,a,Input_Derivative);       
%%% filter(FoldChange)        
InputFoldChange                         = zeros(1,length(CurrentTime),'single');      
InputFoldChange(ZeroIndex:end)          = -Current_a;     
InputFoldChange_LPFonVariable           = filter(Current_b,a,InputFoldChange);  

Beta = 1/Tao;                        
CurrentACT = CalculateAdaptiveThreshold_inline(CurrentTime, Input, K, Beta);
NonAdaptiveThreshold = K * (1-exp(-S0/K)); 

figure('name','Speed definitions- a = 0.3, S0=10','position',[212   144   292   822]);
subplot(3,1,1)
    plot(CurrentTime, Input,'k--'); hold on;
    plot(CurrentTime, CurrentACT,'b-'); hold on;
    plot(CurrentTime, NonAdaptiveThreshold * ones(1,length(CurrentTime),'single'),'b--'); hold on; 
    xlim([-1 10]); ylim([0 12]); ylabel('C')
subplot(3,1,2)
    plot(CurrentTime, Input_Derivative,'k--'); hold on;
    plot(CurrentTime, InputDerivative_LPFonDerivative,'k-'); hold on;
    plot(CurrentTime, Tderivative*ones(1,length(CurrentTime)),'b-'); hold on;
    xlim([-1 10]); ylim([-3.5 0.5]); ylabel('dC/dt')
subplot(3,1,3)
    plot(CurrentTime, InputFoldChange,'k--'); hold on;
    plot(CurrentTime, InputFoldChange_LPFonVariable,'k-'); hold on; 
    plot(CurrentTime, Tfoldchange*ones(1,length(CurrentTime)),'b-'); hold on;
    xlim([-1 10]); ylim([-0.4 0.05]); ylabel('(dC/dt)/C')      
    xlabel('Time [sec]')
      
%% SNR definition for each Input type: (1) S=S0.  (2) S = S0-b*t.  (3) S = S0*exp(-a*t).  Related to Figure S6A       
%%%  Free Parameters  %%%
SampleRate_TimeVec_Example  = 5; % Hz
TotalTime_Example    = 220;  % seconds
CurrentTime          = single(0:1/SampleRate_TimeVec_Example:TotalTime_Example);
S0       = 10;   % microM
NoiseSTD = 0.5;  % microM 

%%%  Thresholds information - use previously extracted parameter information  %%% 
K           = 5.5;   % microM
Tao         = 17;    % seconds

rng('default')
NoiseVec_Raw = randn(1,length(CurrentTime),'single')*NoiseSTD; 

%%%  Input 1: S = S0    
NoiseVec                  = NoiseVec_Raw; 
NoiseVec(NoiseVec<=-S0)   = -0.9999*S0;   % Correct noise to prevent negative concentration values  (S+noise >= 0) 

InputWithNoise = S0+NoiseVec; InputWithNoise(1)=S0;
CurrentACT = CalculateAdaptiveThreshold_inline (CurrentTime, InputWithNoise, K, 1/Tao); % ACT, full vector

Input1          = S0*ones(1,length(CurrentTime));
ACT1            = CurrentACT;
InputWithNoise1 = S0 + NoiseVec;

%%%  Input 2: S = S0-b*t    
b_ind     = 2;       % Only 0.1 microM/sec   
Current_b = b_vec(b_ind);       
Input          = S0 - Current_b*CurrentTime;
Input(Input<0) = NaN;                      

NoiseVec                  = NoiseVec_Raw;                        
NoiseVec(NoiseVec<=-Input)= -0.9999*Input(NoiseVec<=-Input);    % Correct noise to prevent negative concentration values  (S+noise >= 0)         
InputWithNoise = Input+NoiseVec; InputWithNoise(1)=S0;
CurrentACT = CalculateAdaptiveThreshold_inline (CurrentTime, InputWithNoise, K, 1/Tao); % ACT, full vector

Input2          = Input;
ACT2            = CurrentACT;
InputWithNoise2 = Input + NoiseVec;

%%%  Input 3: S = S0*exp(-a*t)   
a_ind          = 2;  % 0.5 [1/sec] 
Current_a      = a_vec(a_ind);   
Input          = S0 * exp( -Current_a*CurrentTime);
Input(Input<0) = NaN;                      

NoiseVec                  = NoiseVec_Raw;                        
NoiseVec(NoiseVec<=-Input)= -0.9999*Input(NoiseVec<=-Input);          % Correct noise to prevent negative concentration values  (S+noise >= 0) 
 
InputWithNoise = Input+NoiseVec; InputWithNoise(1)=S0;
CurrentACT     = CalculateAdaptiveThreshold_inline (CurrentTime, InputWithNoise, K, 1/Tao); % ACT, full vector

Input3          = Input;
ACT3            = CurrentACT;
InputWithNoise3 = Input + NoiseVec;

%%% Figure S6A   
figure('name','SNR definitions, 3 models','position',[287   581   822   190]);
subplot(1,3,1)
    plot(CurrentTime, Input1,'k-'); hold on;
    plot(CurrentTime, ACT1,'b-'); hold on;
    plot(CurrentTime, InputWithNoise1,'g-'); hold on; xlim([0 100]);
    ylim([0 13]); xlim([0 100]); set(gca,'xtick',[0 50 100],'ytick',[0 5 10]);
    ylabel('C')
    xlabel('Time [sec]')
subplot(1,3,2)
    plot(CurrentTime, Input2,'k-'); hold on;
    plot(CurrentTime, ACT2,'b-'); hold on;
    plot(CurrentTime, InputWithNoise2,'g-'); hold on; xlim([0 100]);
    ylim([0 13]); xlim([0 100]); set(gca,'xtick',[0 50 100],'ytick',[0 5 10])
    xlabel('Time [sec]')
subplot(1,3,3)
    plot(CurrentTime, Input3,'k-'); hold on;
    plot(CurrentTime, ACT3,'b-'); hold on;
    plot(CurrentTime, InputWithNoise3,'g-'); hold on; xlim([0 100]);
    ylim([0 13]); xlim([0 100]); set(gca,'xtick',[0 50 100],'ytick',[0 5 10])
    xlabel('Time [sec]')
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Effect of sampling rate on noise filtering, Figure S6B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Input 1: S1 = S0
ModelNames = {'Derivative','ACT','FoldChange'};
S0_ind     = 2;  % 10 microM
rows       = 1;
columns    = 3;
figure('name','Constant Input. ACT(r), Der(k), FC(b)','position',[212  747   1020  219]);     
sp_ind = 0;    
for model_ind = 1:length(ModelNames)    
    sp_ind = sp_ind+1;            
    subplot(rows, columns, sp_ind)                
    CurrentModel = ModelNames{model_ind};
    for SR_ind = 1:length(SampleRate_vec)       
        Y = squeeze(SNR.ConstantInput.(CurrentModel)(S0_ind,SR_ind,:));
        Y(Y>LimitForDisplay.Upper)=LimitForDisplay.Upper; % For display
        Y(Y<LimitForDisplay.Lower)=LimitForDisplay.Lower; % For display            
        loglog(NoiseSTD_vec, Y,LinesPerModel{model_ind}); hold on;
    end
    xlim([5e-6 2e6]);
    ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]);
    set(gca,'xtick',10.^(-5:6),'xticklabel',{'','10^{-4}','','10^{-2}','','10^{0}','','10^2','','10^4','','10^6'},'ytick',10.^(-4:4));
    set(gca,'box','off','xminortick','off','yminortick','off');
   
    ylabel('SNR');    
    xlabel('\sigma_{noise} [\muM]');                     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Effect of signal and noise amplitudes on noise filtering, Figures 6B and S6C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Figure 6B
%%% Input 1: S1 = S0
ModelNames    = {'Derivative','ACT','FoldChange'};
figure('name','Constant Input. ACT(r), Der(k), FC(b)');     
SR_ind = length(SampleRate_vec);     % only 5Hz   
S0_ind = 2; % 1 microM                  
for model_ind = 1:length(ModelNames)    
    CurrentModel = ModelNames{model_ind};
    Y = squeeze(SNR.ConstantInput.(CurrentModel)(S0_ind,SR_ind,:));
    Y(Y>LimitForDisplay.Upper)=LimitForDisplay.Upper; % For display
    Y(Y<LimitForDisplay.Lower)=LimitForDisplay.Lower; % For display            
    loglog(NoiseSTD_vec, Y,LinesPerModel{model_ind}); hold on;
end
xlim([1e-3 1e3]);
ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]);
set(gca,'xtick',10.^(-5:6),'xticklabel',{'','10^{-4}','','10^{-2}','','10^{0}','','10^2','','10^4','','10^6'},'ytick',10.^(-4:4));
set(gca,'box','off','xminortick','off','yminortick','off');
ylabel(['average SNR',char(10),'S_0 = ',S0_vec_Titles{S0_ind},' \muM']);
xlabel('\sigma_{noise} [\muM]');

%% Figure S6C
%%% Input 1: S1 = S0
ModelNames    = {'Derivative','ACT','FoldChange'};
columns       = 1;
rows          = length(S0_vec);
figure('name','Constant Input. ACT(r), Der(k), FC(b)','position',[212   144   292   822]);     
sp_ind = 0;    
SR_ind = length(SampleRate_vec);     % only 5Hz   
for S0_ind = 1:length(S0_vec)                  
    sp_ind = sp_ind+1;            
    subplot(rows, columns, sp_ind)                
    for model_ind = 1:length(ModelNames)    
        CurrentModel = ModelNames{model_ind};
        Y = squeeze(SNR.ConstantInput.(CurrentModel)(S0_ind,SR_ind,:));
        Y(Y>LimitForDisplay.Upper)=LimitForDisplay.Upper; % For display
        Y(Y<LimitForDisplay.Lower)=LimitForDisplay.Lower; % For display            
        loglog(NoiseSTD_vec, Y,LinesPerModel{model_ind}); hold on;
    end
    xlim([5e-6 2e6]);
    ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]);
    set(gca,'xtick',10.^(-5:6),'xticklabel',{'','10^{-4}','','10^{-2}','','10^{0}','','10^2','','10^4','','10^6'},'ytick',10.^(-4:4));
    set(gca,'box','off','xminortick','off','yminortick','off');
    ylabel(['average SNR',char(10),'S_0 = ',S0_vec_Titles{S0_ind},' \muM']);
    
    if S0_ind==length(S0_vec)
        xlabel('\sigma_{noise} [\muM]');
    end                               
end

%%% Input 2: S2 = S0 - b*t
ModelNames    = {'Derivative','ACT','FoldChange'};
rows       = length(S0_vec);
columns    = 1;  
b_ind     = 2;       % Only 0.1 microM/sec   
Current_b = b_vec(b_ind);    
figure('name',['Linear decay, slope= ',num2str(Current_b),'microM/sec. ACT(r), Der(k), FC(b)'],'position',[212   144   292   822]);     
sp_ind = 0;    
SR_ind = length(SampleRate_vec);        
for S0_ind = 1:length(S0_vec)                  
    sp_ind = sp_ind+1;            
    subplot(rows, columns, sp_ind)                
    for model_ind = 1:length(ModelNames)    
        CurrentModel = ModelNames{model_ind};
        Y = squeeze(SNR.LinearInput.(CurrentModel)(S0_ind,SR_ind,b_ind,:));
        Y(Y>LimitForDisplay.Upper)=LimitForDisplay.Upper; % For display
        Y(Y<LimitForDisplay.Lower)=LimitForDisplay.Lower; % For display            
        loglog(NoiseSTD_vec, Y,LinesPerModel{model_ind}); hold on;
    end
    xlim([5e-6 2e6]);
    ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]);
    set(gca,'xtick',10.^(-5:6),'xticklabel',{'','10^{-4}','','10^{-2}','','10^{0}','','10^2','','10^4','','10^6'},'ytick',10.^(-4:4));
    set(gca,'box','off','xminortick','off','yminortick','off');
    
    ylabel(['average SNR',char(10),'Input @ ',S0_vec_Titles{S0_ind},' \muM']);
    
    if S0_ind==length(S0_vec)
        xlabel('\sigma_{noise} [\muM]');
    end                               
end

%%% Input 3: S3 = S0*exp(-a*t)
ModelNames    = {'Derivative','ACT','FoldChange'};
rows       = length(S0_vec); 
columns    = 1; 
a_ind      = 2;  % Only 0.02 [1/sec] 
Current_a = a_vec(a_ind);    
figure('name',['Exponential decay, time constant= ',num2str(Current_a),' 1/sec. ACT(r), Der(k), FC(b)'],'position',[212   144   292   822]);     
sp_ind = 0;    
SR_ind = length(SampleRate_vec);        
for S0_ind = 1:length(S0_vec)                  
    sp_ind = sp_ind+1;            
    subplot(rows, columns, sp_ind)                
    for model_ind = 1:length(ModelNames)    
        CurrentModel = ModelNames{model_ind};
        Y = squeeze(SNR.ExponentInput.(CurrentModel)(S0_ind,SR_ind,a_ind,:));
        Y(Y>LimitForDisplay.Upper)=LimitForDisplay.Upper; % For display
        Y(Y<LimitForDisplay.Lower)=LimitForDisplay.Lower; % For display            
        loglog(NoiseSTD_vec, Y,LinesPerModel{model_ind}); hold on;
    end
    xlim([5e-6 2e6]);
    ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]);
    set(gca,'xtick',10.^(-5:6),'xticklabel',{'','10^{-4}','','10^{-2}','','10^{0}','','10^2','','10^4','','10^6'},'ytick',10.^(-4:4));
    set(gca,'box','off','xminortick','off','yminortick','off');    
    ylabel(['average SNR',char(10),'Input @ ',S0_vec_Titles{S0_ind},' \muM']);

    if S0_ind==length(SampleRate_vec)
        xlabel('\sigma_{noise} [\muM]');             
    end                               
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Effect of Integration on noise filtering of constant signal, Figures 6C, S6D and S6E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Tao_ind_ToDisplay = [20 40 60 75]; % 0.1, 1, 10, 50

%%%  Figure 6C
ModelNames    = {'Derivative','ACT','FoldChange'};
Tao_ind_ToDisplay = [20 40 60 75]; % 0.1, 1, 10, 50
columns       = length(ModelNames);
rows          = 1;   
SR_ind        = length(SampleRate_vec);   % Only 5 Hz  
figure('name',['Constant Signal + noise with LPF on variable, sample rate= ',num2str(SampleRate_vec(SR_ind)),' Hz. ACT(r), Der(k), FC(b)'],...
       'position',[ 417   546   946   302]);  
S0_ind = 2; % 10 microM
           
for model_ind = 1:length(ModelNames)               
    CurrentModel = ModelNames{model_ind};                             
    subplot(rows, columns, model_ind)                
    for Tao_ind = Tao_ind_ToDisplay             
%         Current_Tao = Tao_vec(Tao_ind); 
        if strcmpi(CurrentModel,'ACT')
            Y = squeeze(SNR.ConstantInput_VaryTao.(CurrentModel)(S0_ind,SR_ind,:,Tao_ind));
        else
            Y = squeeze(SNR.ConstantInputWithLPFOnVariable.(CurrentModel)(S0_ind,SR_ind,:,Tao_ind));
        end
        Y(Y>LimitForDisplay.Upper)=LimitForDisplay.Upper; % For display
        Y(Y<LimitForDisplay.Lower)=LimitForDisplay.Lower; % For display            
        loglog(NoiseSTD_vec, Y,LinesPerModel{model_ind}); hold on;
    end
    xlim([1e-3 1e3]);
    ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]);
    set(gca,'xtick',10.^(-5:6),'xticklabel',{'','10^{-4}','','10^{-2}','','10^{0}','','10^2','','10^4','','10^6'},'ytick',10.^(-4:4));
    set(gca,'box','off','xminortick','off','yminortick','off');
    xlabel('\sigma_{noise} [\muM]');            

    if model_ind==1
        ylabel(['SNR',char(10),'Input @ ',S0_vec_Titles{S0_ind},' \muM']);            
    end                               
end

%%%  Figure S6D. Columns are integration times , rows are input amplitudes, X-axis for noise  
ModelNames    = {'Derivative','ACT','FoldChange'};
Tao_ind_ToDisplay = [20 60];       % 0.1, 10
columns       = length(Tao_ind_ToDisplay);
rows          = length(S0_vec);   
SR_ind        = length(SampleRate_vec);   % Only 5 Hz  
figure('name',['Constant Signal + noise with LPF on variable, sample rate= ',num2str(SampleRate_vec(SR_ind)),' Hz. ACT(r), Der(k), FC(b)'],...
       'position',[212   144   584   822]);     
sp_ind = 0;         
for S0_ind = 1:length(S0_vec)              
    for Tao_ind = Tao_ind_ToDisplay             
        Current_Tao = Tao_vec(Tao_ind);                  
        sp_ind = sp_ind+1;            
        subplot(rows, columns, sp_ind)                
        for model_ind = 1:length(ModelNames)    
            CurrentModel = ModelNames{model_ind};
            if strcmpi(CurrentModel,'ACT')
                Y = squeeze(SNR.ConstantInput_VaryTao.(CurrentModel)(S0_ind,SR_ind,:,Tao_ind));
            else
                Y = squeeze(SNR.ConstantInputWithLPFOnVariable.(CurrentModel)(S0_ind,SR_ind,:,Tao_ind));
            end
            Y(Y>LimitForDisplay.Upper)=LimitForDisplay.Upper; % For display
            Y(Y<LimitForDisplay.Lower)=LimitForDisplay.Lower; % For display            
            loglog(NoiseSTD_vec, Y,LinesPerModel{model_ind}); hold on;
        end
        xlim([5e-6 2e6]);
        ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]);
        set(gca,'xtick',10.^(-5:6),'xticklabel',{'','10^{-4}','','10^{-2}','','10^{0}','','10^2','','10^4','','10^6'},'ytick',10.^(-4:4));
        set(gca,'box','off','xminortick','off','yminortick','off');
        if S0_ind==1
            title(['\tau = ',num2str(Current_Tao),' s']);
        elseif S0_ind==length(S0_vec)
            xlabel('\sigma_{noise} [\muM]');            
        end
        if Tao_ind==Tao_ind_ToDisplay(1)
            ylabel(['SNR',char(10),'Input @ ',S0_vec_Titles{S0_ind},' \muM']);            
        end                               
    end
end

%%% Figure S6E. Rows are signal amplitudes (with input/noise values are kept constant), X-axis for integration times  
ModelNames    = {'Derivative','ACT','FoldChange'};
Noise_ind_ToDisplay = 61:20:221;
columns       = 1;
rows          = length(S0_vec);  
SR_ind        = length(SampleRate_vec);   % Only 5 Hz 
figure('name',['Constant Signal + noise with LPF on variable, sample rate = ',num2str(SampleRate_vec(SR_ind)),' Hz. ACT(r), Der(k), FC(b)'],...
    'position',[212   144   195   822]);     
for S0_ind = 1:length(S0_vec)   
    S0         = S0_vec(S0_ind);
    for n_ind = Noise_ind_ToDisplay                       
        Current_Noise = NoiseSTD_vec(n_ind);  
        if (Current_Noise < S0/10) || (Current_Noise >= S0)
            continue
        end
        sp_ind = S0_ind;         
        subplot(rows, columns, sp_ind)                
        for model_ind = 1:length(ModelNames)    
            CurrentModel = ModelNames{model_ind};
            if strcmpi(CurrentModel,'ACT')
                Y = squeeze(SNR.ConstantInput_VaryTao.(CurrentModel)(S0_ind,SR_ind,n_ind,:));
            else
                Y = squeeze(SNR.ConstantInputWithLPFOnVariable.(CurrentModel)(S0_ind,SR_ind,n_ind,:));
%                     Y = squeeze(SNR.ConstantInputWithLPFOnInput.(CurrentModel)(S0_ind,SR_ind,n_ind,:));
            end
            Y(Y>LimitForDisplay.Upper)=LimitForDisplay.Upper; % For display
            Y(Y<LimitForDisplay.Lower)=LimitForDisplay.Lower; % For display            
            loglog(Tao_vec, Y,LinesPerModel{model_ind}); hold on;
        end
        xlim([5e-3 1e2]);
        ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]);
        set(gca,'xtick',10.^(-3:2),'ytick',10.^(-4:4));        
        set(gca,'box','off','xminortick','off','yminortick','off');    
        if S0_ind==length(S0_vec)
            xlabel('Integration time [sec]');
        end
%         ylabel(['SNR',char(10),'Input= ',S0_vec_Titles{S0_ind},', \sigma_{noise}= ',num2str(Current_Noise)]);                                                      
        ylabel(['SNR',char(10),'I= ',S0_vec_Titles{S0_ind},', \sigma= ',num2str(Current_Noise)]);                                                      
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Effect of Integration on speed for exponential decay signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%  Figure S6F. Columns are input decay rates times , rows are input amplitudes, X-axis for integration times    
ModelNames    = {'Derivative_LPFOnVariable','ACT','FoldChange_LPFOnVariable'};
a_ind_ToShow  = [1 3 5];
columns       = length(a_ind_ToShow);
rows          = length(S0_vec);  

figure('name','Effect Of integration time on Speed. ACT(r), Der(k), FC(b)','position',[212   144   585   822]); 
column=0;
for a_ind = a_ind_ToShow        
    Current_a = a_vec1(a_ind); 
    column=column+1;
    for S0_ind = 1:length(S0_vec)       
%         S0         = S0_vec(S0_ind);        
        sp_ind = (S0_ind-1)*length(a_ind_ToShow) + column;           
        subplot(rows, columns, sp_ind)                
        for model_ind = 1:length(ModelNames)    
            CurrentModel = ModelNames{model_ind};
            Y            = squeeze(Speed.(CurrentModel)(:,a_ind,S0_ind));
        
            Y(Y>LimitForDisplay.Upper)   = LimitForDisplay.Upper;   % For display
            Y(Y<LimitForDisplay.Lower)   = LimitForDisplay.Lower;   % For display   
            Y(isnan(Y))                  = LimitForDisplay.Lower;   % NaNs are "no activation" errors, Speed --> 0
            loglog(Tao_vec, Y,LinesPerModel{model_ind}); hold on;
        end
        xlim([5e-3 1e2]);
        ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]);
        set(gca,'xtick',10.^(-3:2),'ytick',10.^(-2:2));
        set(gca,'box','off','xminortick','off','yminortick','off');    
        if column==1
            ylabel(['Speed [1/s]',char(10),'Input= ',S0_vec_Titles{S0_ind}]);
        end 
        if S0_ind==length(S0_vec)
            xlabel('Integration time [sec]');
        elseif S0_ind==1
            title(['a = ',num2str(Current_a)]);                                        
        end
    end
end

%%% Figure 6E. Speed performance (Y axis) for different input decay rates (X axis) ans integration times (curves)
figure('name','Speed dependence on integration times (0.1, 1, 10, 50 sec) for various input decay rates (X axis)','position',[ 672   235   296   636]);
subplot(3,1,1); 
    Y = squeeze(Speed2.ACT');
    Y(Y>LimitForDisplay.Upper) = LimitForDisplay.Upper;
    Y(Y<LimitForDisplay.Lower) = LimitForDisplay.Lower;    
    Y(isnan(Y))                = LimitForDisplay.Lower;    
    loglog(a_vec2, Y, 'r-'); ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]); 
    ylabel('Speed [1/sec]')
    title('ACT'); 
    xlim([a_vec2(1) a_vec2(end)]); 
subplot(3,1,2); 
    Y = squeeze(Speed2.Derivative_LPFOnVariable');
    Y(Y>LimitForDisplay.Upper) = LimitForDisplay.Upper;
    Y(Y<LimitForDisplay.Lower) = LimitForDisplay.Lower;    
    Y(isnan(Y))                = LimitForDisplay.Lower;    
    loglog(a_vec2, Y, 'k-'); ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]); 
    ylabel('Speed [1/sec]')
    title('Derivative'); 
    xlim([a_vec2(1) a_vec2(end)]); 
subplot(3,1,3); 
    Y = squeeze(Speed2.FoldChange_LPFOnVariable');
    Y(Y>LimitForDisplay.Upper) = LimitForDisplay.Upper;
    Y(Y<LimitForDisplay.Lower) = LimitForDisplay.Lower;    
    Y(isnan(Y))                = LimitForDisplay.Lower;    
    loglog(a_vec2, Y, 'b-'); ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]); 
    ylabel('Speed [1/sec]')
    title('Fold Change'); 
    xlim([a_vec2(1) a_vec2(end)]); 
    xlabel('Input decay rate [1/sec]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Tradeoff between noise filtering capacity and speed. Figure 6F 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
S0         = 10;  % microM
Noise      = 1;   % microM
a          = 0.3; % 1/sec 
SampleRate = 5;   % 1/sec 
S0_ind = find(S0_vec==S0,1);
n_ind  = find(NoiseSTD_vec==Noise,1);
% a_ind  = find(a_vec==a,1);
a_ind  = find(a_vec1==a,1);
SR_ind = find(SampleRate_vec==SampleRate,1);
for i=1
figure('name','Tradeoff figure','position',[274  546  1391  253]);
subplot(1,3,1); 
    yyaxis right
    Y = squeeze(SNR.ConstantInputWithLPFOnVariable.Derivative(S0_ind,SR_ind,n_ind,:));
    Y(Y>LimitForDisplay.Upper) = LimitForDisplay.Upper;
    Y(Y<LimitForDisplay.Lower) = LimitForDisplay.Lower;    
    loglog(Tao_vec, Y); ylim([LimitForDisplay.Lower/1.9 LimitForDisplay.Upper*1.9]); 
    ylabel('SNR')
    yyaxis left
    Y = squeeze(Speed.Derivative_LPFOnVariable(:,a_ind,S0_ind));
    Y(Y>LimitForDisplay.Upper) = LimitForDisplay.Upper;
    Y(Y<LimitForDisplay.Lower) = LimitForDisplay.Lower;    
    Y(isnan(Y))                = LimitForDisplay.Lower;    
    loglog(Tao_vec, Y); ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]); 
    ylabel('Speed [1/sec]')
    title('Derivative'); 
    xlim([0.1 Tao_vec(end)]); 
subplot(1,3,2); 
    yyaxis right
    Y = squeeze(SNR.ConstantInputWithLPFOnVariable.FoldChange(S0_ind,SR_ind,n_ind,:));
    Y(Y>LimitForDisplay.Upper) = LimitForDisplay.Upper;
    Y(Y<LimitForDisplay.Lower) = LimitForDisplay.Lower;    
    loglog(Tao_vec, Y); ylim([LimitForDisplay.Lower/1.9 LimitForDisplay.Upper*1.9]); 
    ylabel('SNR')
    yyaxis left
    Y = squeeze(Speed.FoldChange_LPFOnVariable(:,a_ind,S0_ind));
    Y(Y>LimitForDisplay.Upper) = LimitForDisplay.Upper;
    Y(Y<LimitForDisplay.Lower) = LimitForDisplay.Lower;    
    Y(isnan(Y))                = LimitForDisplay.Lower;    
    loglog(Tao_vec, Y); ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]); 
    ylabel('Speed [1/sec]')
    title('Fold Change'); 
    xlim([0.1 Tao_vec(end)]); 
subplot(1,3,3); 
    yyaxis right
    Y = squeeze(SNR.ConstantInput_VaryTao.ACT(S0_ind,SR_ind,n_ind,:));
    Y(Y>LimitForDisplay.Upper) = LimitForDisplay.Upper;
    Y(Y<LimitForDisplay.Lower) = LimitForDisplay.Lower;    
    loglog(Tao_vec, Y); ylim([LimitForDisplay.Lower/1.9 LimitForDisplay.Upper*1.9]); 
    ylabel('SNR')
    yyaxis left
    Y = squeeze(Speed.ACT(:,a_ind,S0_ind));
    Y(Y>LimitForDisplay.Upper) = LimitForDisplay.Upper;
    Y(Y<LimitForDisplay.Lower) = LimitForDisplay.Lower;    
    Y(isnan(Y))                = LimitForDisplay.Lower;    
    loglog(Tao_vec, Y); ylim([LimitForDisplay.Lower/2 LimitForDisplay.Upper*2]); 
    ylabel('Speed [1/sec]')
    title('ACT'); 
    xlim([0.1 Tao_vec(end)]); 
end


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







