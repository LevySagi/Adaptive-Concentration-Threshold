function AttractionTests = Information_Exp1
% 20170501  

AttractionTests.TrackFile              = '';
AttractionTests.FlowDelay.Device       = 'SB1_Setup2';
%                                    buffer-buffer   High-buffer       High-Low          High-intermediate-Low-buffer             Dye test             
AttractionTests.DelayBetween_2valve_and_8valve_in_seconds = 5;
AttractionTests.positions_valve1    = [ 1 1 1 2     2 2 2 2 2 2     2 2 2 2 2 2  2    2 4 4 2  2 4 4 2  2 4 4 2  2 4 4 2   1 6 6 6 6 6 6 6 6 1 1 ];  
AttractionTests.directions_valve1   = [ 1 1 1 1     1 1 1 1 1 1     1 1 1 1 1 1  1    1 1 1 0  1 1 1 0  1 1 1 0  1 1 1 0   1 0 1 1 1 1 1 1 1 1 1 ]; 
AttractionTests.positions_valve2    = [ 1 1 1 1     1 1 1 1 1 1     4 4 4 4 4 4  1    3 3 1 1  3 3 1 1  3 3 1 1  3 3 1 1   1 1 1 1 1 1 1 1 1 7 1 ];  
AttractionTests.directions_valve2   = [ 1 1 1 1     1 1 1 1 1 1     1 1 1 1 1 1  0    1 1 0 1  1 1 0 1  1 1 0 1  1 1 0 1   1 1 1 1 1 1 1 1 1 0 1 ];
AttractionTests.positions_FastValve = [ 0 1 0 1     0 1 0 1 0 1     0 1 0 1 0 1  1    0 1 0 1  0 1 0 1  0 1 0 1  0 1 0 1   0 0 1 0 1 0 1 0 1 1 1 ];
AttractionTests.StimulusTimes       = [ 2 2 2 2     2 2 2 2 2 5     2 2 2 2 2 2  5    2 2 2 2  2 2 2 2  2 2 2 2  2 2 2 2   1 7 2 2 2 2 2 2 2 7 7 ];

AttractionTests.Concentration          =  [0 1e-6 3e-7 5e-8 0 0 0 0];  % valve position 1:8. Butanone concentration in units of dilution (no dilution = 11.16M) 
AttractionTests.Concentration_Strings  =  {'Buffer','10^-^6Bu','3*10^-^7Bu','5*10^-^8Bu','Buffer','Dye/Buffer','Buffer/Dye','Buffer'} ;  % valve position 1:8. Butanone concentration in units of dilution 
AttractionTests.ExpType                = '';
AttractionTests.Keywords               = {'Butanone','SB1','N2','Dye','Setup2','Pulses'};

% arenas 1,2 = near inlets with fast flow, 3,4 = near outlets with slow flow  
AttractionTests.ArenaNumbers                  = 1:4;  
AttractionTests.ArenaLocationRelativeToInlets = {'Inlet','Inlet','Outlet','Outlet'};
AttractionTests.ArenaLocationRightLeft        = {'Left','Right','Left','Right'};  % Right is the far side from the sash, left is the close side
AttractionTests.FlowDelay.Device              = 'SB1';
AttractionTests.Genotypes                     = {'N2','N2','N2','N2'};  
AttractionTests.TrainingCondition             = {'AfterTreatment','Naive','Naive','AfterTreatment'}; 
AttractionTests.AssayTimesInHours             = 0; 

return
