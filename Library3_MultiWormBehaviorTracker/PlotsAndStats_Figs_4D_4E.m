function PlotsAndStats_Figs_4D_4E
% 
% Figures 4D and 4E, Behavioral responses to fast and slow decrease in odor concentration 
% 
% September 2019, Sagi Levy
%  
%%  load processed data
load('BehavioralResponses_Figs_4D_4E.mat','ExpInfo',...
            'ChangeInAversiveBehaviorProbability','ChangeInAversiveBehaviorProbability_Dynamics',...
            'ChangeInSpeed','ChangeInSpeed_Dynamics');

FrameRate          = ExpInfo.FrameRate; 
FlowDelay          = ExpInfo.FlowDelay; 
PatternsOnset_Fast = ExpInfo.PatternsOnset_Fast; 
StimuliTitles      = ExpInfo.StimuliTitles;      % Corresponds to columns 1:4 in data matrices 
TimeVecInSec       = ExpInfo.TimeVecInSec;       % Corresponds to X-axis in data matrices 
ConditionsTitles   = ExpInfo.ConditionsTitles;   % Corresponds to rows 1:4 in data matrices 

%%  Plots and Stats
IncludeAnimalsAfterTreatment    = false;

%%% Plot change in aversive behavior probability
[~, ttest_Pvalue_FastHigher, ~] = ComputeSlowAndFastSwitchesStats(ChangeInAversiveBehaviorProbability);
YlimitsDynamics                 = [0.25 0.75];
PlotComparisonBetweenSlowAndFastSwitches(ChangeInAversiveBehaviorProbability, ChangeInAversiveBehaviorProbability_Dynamics, TimeVecInSec, StimuliTitles, 'P(after)-P(before)', IncludeAnimalsAfterTreatment, YlimitsDynamics)
                      
%%% Plot change in speed
[~, ~, ttest_Pvalue_FastSlower] = ComputeSlowAndFastSwitchesStats(ChangeInSpeed);
YlimitsDynamics                 = [110 180];
PlotComparisonBetweenSlowAndFastSwitches(ChangeInSpeed, ChangeInSpeed_Dynamics, TimeVecInSec, StimuliTitles, '(S(after)-S(before))/S(before)', IncludeAnimalsAfterTreatment, YlimitsDynamics)

%%% Display Fast vs. slow dye flow
figure('position',[680   725   396   253]); 
colors = {'k','k','r','r'};
for ar=[1 3]%  Naive animals arenas
    CurrentDelayInSec    = FlowDelay.DelayBetweenValveAndArenas(ar);
    CurrentDelayInFrames = round(CurrentDelayInSec*FrameRate);
    CurrentDye           = PatternsOnset_Fast(ar,:);
    CurrentTime          = ((0:length(CurrentDye)-1) - CurrentDelayInFrames)/FrameRate;
    plot(CurrentTime, 1-CurrentDye,'color',colors{ar},'linewidth',0.75); hold on    
end
xlim([0 15]); ylim([-0.02 1.05]); ylabel('dye concentration'); xlabel('Time [sec]')

return

function BehaviorCode
%  My Codes:  
%    Low Level
%             0 'Out of bounds'
%             1 'Forward'
%             2 'Curve'
%             3 'Pause'
%             4 'Reverse'
%             5 'Omega'
%             6 'SharpTurn'
%
%    High Level
%             0 'Out of bounds'
%             1 'Forward'
%             2 'Curve'
%             3 'Pause'
%             4 'Reversal'
%             5 'Omega-Pause'
%             6 'SharpTurn'
%             7 'Pirouette. Initial reversal'
%             8 'Pirouette. After reversal'



return

function PlotComparisonBetweenSlowAndFastSwitches(StructureIN, AllRawValues, TimeVecInSec, StimuliTitles, Ytitle, IncludeAnimalsAfterTreatment, YLimits4Dynamics)
         
PlotDynamics = false;

Properties.BarProps.EdgeColor      = {'k','k'};
Properties.BarProps.FaceColor      = {[0.8 0.8 0.8],[0.2 0.2 0.2]};
Properties.ErrorBarProps.linewidth = {3,3,3,3};
Properties.ErrorBarProps.color     = {'k','k'};
Properties.ErrorBarProps.color     = {'k','k'};
Properties.gcaProps.xticklabel     = StimuliTitles;   % corresponding to ticks of 1:size(Y,1)   
Properties.gcaProps.xlim           = [0.2 length(StimuliTitles)+1-0.2];
Properties.xlabelProps.String      = 'stimuli';
Properties.ylabelProps.String      = 'Aversive behavior probability';

LEN = size(StructureIN.MEANS,2);

% Compare slow to fast switches:
figure('name','Naive: Compare slow to fast switches','position',[151   267   504   372]);
barwitherr_SL([StructureIN.MEANS(3,:); StructureIN.MEANS(1,:)], [StructureIN.SEMs(1,:); StructureIN.SEMs(3,:)], Properties);
ylabel(Ytitle)
xlabel('Concentration after switch (starts at 10^-^6)')
legend('Slow switch','Fast switch')
if IncludeAnimalsAfterTreatment
    figure('name','AfterTreatment: Compare slow to fast switches','position',[151   267   504   372]);
    barwitherr_SL([StructureIN.MEANS(4,:); StructureIN.MEANS(2,:)], [StructureIN.SEMs(1,:); StructureIN.SEMs(3,:)], Properties);
    ylabel(Ytitle)
    xlabel('Concentration after switch (starts at 10^-^6)')
    legend('Slow switch','Fast switch')
end

if PlotDynamics
    %% Plot behavioral response dynamics
    figure('name','Naive: Compare slow to fast switches dynamics','position', [151   267   268   612]);
    for pulse_ind=1:LEN
        subplot(LEN,1,pulse_ind);
        % Fast 
        ConditionsIndex = 1;
        Vec = squeeze(AllRawValues.MeanDynamics(ConditionsIndex,pulse_ind,:));
        SEM = squeeze(AllRawValues.SEMDynamics(ConditionsIndex,pulse_ind,:));
        upper = Vec+SEM;
        lower = Vec-SEM;
        fill([TimeVecInSec, TimeVecInSec(end:-1:1)],[upper; lower(end:-1:1)]',[0.85 0.85 0.85 ],'edgecolor','none','FaceAlpha',1); hold on;              
        plot(TimeVecInSec, Vec,'k-','linewidth',0.75); hold on;    
        % Slow 
        ConditionsIndex = 3;
        Vec = squeeze(AllRawValues.MeanDynamics(ConditionsIndex,pulse_ind,:));
        SEM = squeeze(AllRawValues.SEMDynamics(ConditionsIndex,pulse_ind,:));
        upper = Vec+SEM;
        lower = Vec-SEM;
        fill([TimeVecInSec, TimeVecInSec(end:-1:1)],[upper; lower(end:-1:1)]',[1 0.8 0.8],'edgecolor','none','FaceAlpha',1); hold on;              
        plot(TimeVecInSec, Vec,'r-','linewidth',0.75); hold on;
        xlim(TimeVecInSec([1 end]));   
        if exist('YLimits4Dynamics','var')
            ylim(YLimits4Dynamics)
        end
    end
    if IncludeAnimalsAfterTreatment
        % Compare slow to fast switches dynamics:
        figure('name','AfterTreatment: Compare slow to fast switches dynamics','position', [151   267   268   612]);
        for pulse_ind=1:LEN
            subplot(LEN,1,pulse_ind);
            % Fast 
            ConditionsIndex = 2;
            Vec = squeeze(AllRawValues.MeanDynamics(ConditionsIndex,pulse_ind,:));
            SEM = squeeze(AllRawValues.SEMDynamics(ConditionsIndex,pulse_ind,:));
            upper = Vec+SEM;
            lower = Vec-SEM;
            fill([TimeVecInSec, TimeVecInSec(end:-1:1)],[upper; lower(end:-1:1)]',[0.85 0.85 0.85 ],'edgecolor','none','FaceAlpha',1); hold on;              
            plot(TimeVecInSec, Vec,'k-','linewidth',0.75); hold on;    
            % Slow 
            ConditionsIndex = 4;
            Vec = squeeze(AllRawValues.MeanDynamics(ConditionsIndex,pulse_ind,:));
            SEM = squeeze(AllRawValues.SEMDynamics(ConditionsIndex,pulse_ind,:));
            upper = Vec+SEM;
            lower = Vec-SEM;
            fill([TimeVecInSec, TimeVecInSec(end:-1:1)],[upper; lower(end:-1:1)]',[1 0.8 0.8],'edgecolor','none','FaceAlpha',1); hold on;              
            plot(TimeVecInSec, Vec,'r-','linewidth',0.75); hold on;
            xlim(TimeVecInSec([1 end]));   
            if exist('YLimits4Dynamics','var')
                ylim(YLimits4Dynamics)
            end
        end
    end      

    %% Relative Scale: Compare slow to fast switches dynamics:
    ZeroIndex          = find(TimeVecInSec==0,1);
    figure('name','Relative Scale. Naive: Compare slow to fast switches dynamics','position', [151   267   268   612]);
    for pulse_ind=1:LEN
        subplot(LEN,1,pulse_ind);
        % Fast 
        ConditionsIndex = 1;
        Vec = squeeze(AllRawValues.MeanDynamics(ConditionsIndex,pulse_ind,:));
        Vec = Vec-Vec(ZeroIndex);
        SEM = squeeze(AllRawValues.SEMDynamics(ConditionsIndex,pulse_ind,:));
        upper = Vec+SEM;
        lower = Vec-SEM;
        fill([TimeVecInSec, TimeVecInSec(end:-1:1)],[upper; lower(end:-1:1)]',[0.85 0.85 0.85 ],'edgecolor','none','FaceAlpha',1); hold on;              
        plot(TimeVecInSec, Vec,'k-','linewidth',0.75); hold on;    
        % Slow 
        ConditionsIndex = 3;
        Vec = squeeze(AllRawValues.MeanDynamics(ConditionsIndex,pulse_ind,:));
        Vec = Vec-Vec(ZeroIndex);
        SEM = squeeze(AllRawValues.SEMDynamics(ConditionsIndex,pulse_ind,:));
        upper = Vec+SEM;
        lower = Vec-SEM;
        fill([TimeVecInSec, TimeVecInSec(end:-1:1)],[upper; lower(end:-1:1)]',[1 0.8 0.8],'edgecolor','none','FaceAlpha',1); hold on;              
        plot(TimeVecInSec, Vec,'r-','linewidth',0.75); hold on;
        xlim(TimeVecInSec([1 end]));   
        if exist('YLimits4Dynamics','var')
            Yscale = diff(YLimits4Dynamics)/2;
            ylim([-Yscale Yscale])
        end
    end
    if IncludeAnimalsAfterTreatment
        % Compare slow to fast switches dynamics:
        figure('name','AfterTreatment: Compare slow to fast switches dynamics','position', [151   267   268   612]);
        for pulse_ind=1:LEN
            subplot(LEN,1,pulse_ind);
            % Fast 
            ConditionsIndex = 2;
            Vec = squeeze(AllRawValues.MeanDynamics(ConditionsIndex,pulse_ind,:));
            Vec = Vec-Vec(ZeroIndex);
            SEM = squeeze(AllRawValues.SEMDynamics(ConditionsIndex,pulse_ind,:));
            upper = Vec+SEM;
            lower = Vec-SEM;
            fill([TimeVecInSec, TimeVecInSec(end:-1:1)],[upper; lower(end:-1:1)]',[0.85 0.85 0.85 ],'edgecolor','none','FaceAlpha',1); hold on;              
            plot(TimeVecInSec, Vec,'k-','linewidth',0.75); hold on;    
            % Slow 
            ConditionsIndex = 4;
            Vec = squeeze(AllRawValues.MeanDynamics(ConditionsIndex,pulse_ind,:));
            Vec = Vec-Vec(ZeroIndex);
            SEM = squeeze(AllRawValues.SEMDynamics(ConditionsIndex,pulse_ind,:));
            upper = Vec+SEM;
            lower = Vec-SEM;
            fill([TimeVecInSec, TimeVecInSec(end:-1:1)],[upper; lower(end:-1:1)]',[1 0.8 0.8],'edgecolor','none','FaceAlpha',1); hold on;              
            plot(TimeVecInSec, Vec,'r-','linewidth',0.75); hold on;
            xlim(TimeVecInSec([1 end]));   
            if exist('YLimits4Dynamics','var')
                Yscale = diff(YLimits4Dynamics)/2;
                ylim([-Yscale Yscale])
            end
        end
    end      
end
return

function [ttest_Pvalue, ttest_Pvalue_FastHigher, ttest_Pvalue_FastSlower] = ComputeSlowAndFastSwitchesStats(StructureIN)

%% Initialization
AnalyzeAnimalsAfterTreatment = false;
LEN                          = size(StructureIN.MEANS,2);
ttest_Pvalue.Naive           = zeros(1,LEN);
if AnalyzeAnimalsAfterTreatment
    ttest_Pvalue.AfterTreatment  = zeros(1,LEN);
end
ttest_Pvalue_FastHigher      = ttest_Pvalue;
ttest_Pvalue_FastSlower      = ttest_Pvalue;

%%  Paired ttest: Compare slow to fast switches
%%% Both tails
FastIndex = 1; 
SlowIndex = 3;
for condition_ind = 1:LEN
    [h,ttest_Pvalue.Naive(condition_ind)] = ttest( StructureIN.MeanValues{FastIndex,condition_ind}, StructureIN.MeanValues{SlowIndex,condition_ind} );
end
if AnalyzeAnimalsAfterTreatment
    FastIndex = 2; 
    SlowIndex = 4;
    for condition_ind = 1:LEN
        [h,ttest_Pvalue.AfterTreatment(condition_ind)] = ttest( StructureIN.MeanValues{FastIndex,condition_ind}, StructureIN.MeanValues{SlowIndex,condition_ind} );
    end
end

%%% One tail: Fast > Slow
FastIndex = 1; 
SlowIndex = 3;
for condition_ind = 1:LEN
    [h,ttest_Pvalue_FastHigher.Naive(condition_ind)] = ttest( StructureIN.MeanValues{FastIndex,condition_ind}, StructureIN.MeanValues{SlowIndex,condition_ind} ,'Tail','right');
end
if AnalyzeAnimalsAfterTreatment
    FastIndex = 2; 
    SlowIndex = 4;
    for condition_ind = 1:LEN
        [h,ttest_Pvalue_FastHigher.AfterTreatment(condition_ind)] = ttest( StructureIN.MeanValues{FastIndex,condition_ind}, StructureIN.MeanValues{SlowIndex,condition_ind},'Tail','right' );
    end
end

%%% One tail: Fast < Slow
FastIndex = 1; 
SlowIndex = 3;
for condition_ind = 1:LEN
    [h,ttest_Pvalue_FastSlower.Naive(condition_ind)] = ttest( StructureIN.MeanValues{FastIndex,condition_ind}, StructureIN.MeanValues{SlowIndex,condition_ind} ,'Tail','left');
end
if AnalyzeAnimalsAfterTreatment
    FastIndex = 2; 
    SlowIndex = 4;
    for condition_ind = 1:LEN
        [h,ttest_Pvalue_FastSlower.AfterTreatment(condition_ind)] = ttest( StructureIN.MeanValues{FastIndex,condition_ind}, StructureIN.MeanValues{SlowIndex,condition_ind},'Tail','left' );
    end
end

return


