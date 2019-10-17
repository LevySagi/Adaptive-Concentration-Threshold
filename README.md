# Adaptive-Concentration-Threshold

General information

•	This repository codes used for data analysis, plots and statistics in the manuscript: “An adaptive-threshold mechanism for odor sensation and animal navigation” by Sagi Levy and Cori Bargmann which was published in Neuron. This manuscript should be cited for any use of the attached code or data. 
•	All codes are written in Matlab (MathWorks) and detailed commentary and instructions are available within each function. Codes are arranged in five libraries as discussed below. 
•	Raw Data repository: most raw data are stored in small-size ‘.mat’ files and are available in this repository. All data, including large-size raw and processed data, can be downloaded from Mendeley Data, http://dx.doi.org/10.17632/bgvk49bjcx.1 
•	Input raw mat files are loaded within the functions. Please make sure that files are present in the relevant directory or Matlab path. Alternatively, please update correct file paths within the functions.



Code organization and instructions

Library 1: Imaging of AWC(ON) calcium activity in animals immobilized in microfluidic devices, Related to figures 1, 2, 5, S1-S3 and S5

•	Basic image processing: extraction of raw fluorescence in different AWCON compartments of immobilized worms. Run functions in the order shown above:
  •	ImmobilizedWorm_SomaTracker.m
  •	ImmobilizedWorm_ProcessTracker.m
  •	ImmobilizedWorm_DyeTracker.m
  •	Notes: run functions in the order shown above. See inline comments for details.

•	Data processing 
  •	Name: 	 ProcessData_Figs_1_2_5.m
  •	Input:	 Data_AWConImaging_Raw.mat 
  •	Output:  Data_AWConImaging_Processed.mat 
  •	Notes:	 The function includes slow inline bootstrap functions to estimate data standard deviations. For faster performance please use parallel computing for bootstrap functions (see inline commentary). Alternatively, to avoid slow running time of bootstrap functions, files with intermediate data processing are also available. See inline comments for more details.

•	Plots and statistics 
  •	Name: 	PlotsAndStats_Figs_1_2_5.m
  •	Input:	Data_AWConImaging_Processed.mat
 

 
Library 2: AWC(ON) calcium and behavioral responses in animals freely moving in microfluidic devices, Related to figures 3, 4 and S4

•	General notes
  •	To run the code, users must define a directory named 'Matlab_Temporary_Files' in 'C:\'. 
  •	In all functions, concentration is shown in units of dilution (undiluted butanone = 11.16 M)
  •	Run functions in the order shown below
  •	3 important structures are used in all functions:
    •	File	   Information about movie files, experiment information, and parameters for analysis
    •	Tracks   Animal data for segmented tracks
    •	Data	   Animal data, organized in (nxm) matrices, with ‘n’ the animal number and ‘m’ is the frame number.

•	Functions List
1)    Name:                 [File, All_Files] = Tracker_NeuronalImagingGradientDevice_v02
      External functions:   (a)	FragmentTracker_SL_ImagingSetup_v09(Fragment), Image processing in a single movie segments.
                            (b)	SegmentationSettings = SegmentationSettings_ImagingSetup_v02, Includes behavior segmentation parameters
    	Description:         	Basic image processing and extraction of animal features 
    	Output:     	        Directory named 'trackfile parts' containing mat files, each file contains image processing of one movie fragment  

2)	  Name: 			          File = StitchTracks_ExtractNeuronPosition_02(File)
    	Description:         	concatenate data from fragments, Stitch tracks, correct midline and head/tail positions, extract neuron position and raw fluorescence value   
    	Output:     	        This function saves the following ‘mat’ files: 
                          	(a)	'MovieName_StitchedPatternMatrices_Arena1.mat’, with the variable 'PatternMatrix'
                          	(b)	'MovieName_StitchedAndManuallyCorrected_Arena1.mat', with all variables, except of 'PatternMatrix' and neuron fluorescence data
                          	(c)	'MovieName_NeuronsDataBeforeManualCorrection_Arena1.mat', with all variables, except of 'PatternMatrix'    

3)    Name:                 AnimalAndNeuronTrackerGUI_02(File, BehaviorDataMatrices, SwitchPathDefinitionInTracks)
    	Inputs:              	(a)	BehaviorDataMatrices  = [], The GUI contains an option to check BehaviorDataMatrices from part 4, see inline commentary. 
                          	(b)	SwitchPathDefinitionInTracks = false, Default=false, when the GUI runs on Tracks before behavior segmentation (part 4).  ‘true’ only if the GUI runs on Tracks after behavior segmentation (part 4)
    	External functions:  	Tracks = ExtractNeuralActivityFromTracks_AfterGUI_02(Tracks, File, plotme).
    	Description:         	Optional features include: Manual inspection of tracking results. Check of flaged frames. Automated or manual correction of tracks. Generate movies. Note that while manual inspection is recommended, it is not required, especially if movies have high quality. 
                          	Required part (Save and Exit button): runs the external function: ExtractNeuralActivityFromTracks_AfterGUI_02  
    	Output:     	        This function saves the file:  'MovieName_AfterNeuronPositionManualCorrection.mat'     

4)    Name:                 [Tracks, File, TracksStats, Data] = SegmentBehavior_ImagingInGradient_v03 (File, Tracks, TimeOfGradientStartValveSwitch, TimeOfGradientEndValveSwitch)
    	Inputs:              	'TimeOfGradientStartValveSwitch' and 'TimeOfGradientEndValveSwitch' are the elapsed time in which gradient valves were turned on and off, respectively.
    	External functions:   AnalyzeLocomotionProperties_ImagingSetup_v2.m   >> Locomotion 
                            AnalyzeMorphologyProperties_ImagingSetup_v3.m   >> Morphology 
                            SegmentBehavior_SingleTrack_ImagingSetup_v3.m   >> Behavior 
    	Description:         	Analyzes animal locomotion, morphology and behavior  
    	Output:     	        This function saves the file 'MovieName_DataMatrices_Arena1.mat’ that contains 'Data' structure that stores all data so far. Note that only low level behavior is used in the paper, and that this code is simplified later (see part 6: 'StitchAndProcessGradientImagingExperiments_v02')

5)	  Name:                 [File, Data] = CalculateActivityAndStimulusFeatures_04(File, Data, ExpInfo, plotme);
    	Inputs:               ‘File’ and ‘Data’ are the outputs of the previous section. ‘ExpInfo’ contains information about gradient concentrations in the current experiment. This function also loads ‘ManualInspection_Activity.mat’ that contains data from manual inspection of activity patterns for single frame resolution.
    	Description:         	Filter and cleanup activity and stimulus data. Define gradient locomotion features. Extract adaptive concentration threshold.
    	Output:               This function saves the file: 'MovieName_DataMatrices_WithGradient.mat' that contains Data structure that stores all data so far

6)    Name:                 Data = StitchAndProcessGradientImagingExperiments_v02;
    	Inputs:  	            This function loads data of all previous experiments. Experiment names and paths are defined in the function first rows. It also loads ‘ManualInspection_GradientLocomotion.mat’ with data from manual inspection of gradient locomotion patterns for single frame resolution.  
    	Description:         	Concatenate data from all experiments. Additional processing, including analysis of gradient conditions, navigation decisions and behavior code simplification.  
    	Output:     	        This function saves the file 'GradientImagingMetaAnalysis.mat' that contains raw and processed data from all experiments
    
7)    Name:                 PlotsAndStats_Figs_3_4(Data) 
    	Inputs:               ‘Data’ structure from previous section, stored in ‘GradientImagingMetaAnalysis.mat’, and colormap for behavior matrices stored at ‘ColormapLowLevelBehavior.mat’    
    	Description:          Sort data, plot and calculate statistics related to all panels in figure 3 and 4A-C.



Library 3: behavioral responses in animals freely moving in high-throughput microfluidic devices, Related to figures 4D, 4E, S4I and S4J

•	General notes
  •	To run the code, users must define a directory named 'Matlab_Temporary_Files' in 'C:\'. 
  •	Run functions in the order shown below
  •	3 important structures are used in all functions:
    •	File    Information about movie files, experiment information, and parameters for analysis
    •	Tracks  Animal data for segmented tracks
    •	Data    Animal data, organized in (nxm) matrices, with ‘n’ the animal number and ‘m’ is the frame number.

•	Functions List

1)	  Name: 			        [File, All_Files] = MultiWormTracker_HighThroughputDevices_SL
    	External functions: (a)	SegmentationSettings_ScreeningSetup_SL_v02, Includes behavior segmentation parameters and matches high-throughput microfluidic device behaviors
                        	(b)	FragmentTracker_SL_ScreeningSetup_v10, Image processing in a single movie segment.
                        	(c)	StitchFiles_MultipleSessions_v02, Concatenates all animal tracks in a single arena  
                        	(d)	SegmentBehavior_HighThroughputDevice_02, Analyzes animal locomotion, morphology and behavior, Stitch tracks and generates ‘Data’ Structures with data for all animals in a single arena. This function uses three additional external functions that analyzes dynamic features of a single animal:
                            		I.	AnalyzeLocomotionProperties_SL_HighThroughputDevice_02, locomotion analysis
                            		II.	AnalyzeMorphologyProperties_HighThroughputDevice_02, morphology analysis
                            		III.	SegmentBehavior_SingleTrack_HighThroughputDevice_02, behavior analysis
                        	(e)	FragmentTrackDye_02, Tracks dye pattern in one movie fragment
    	Description:        Image processing, extraction of animal features, stitching animal tracks, segmentation of animal behavior and dye tracking.  
    	Outputs:     	      (a)	MovieName_DataMatrices_ArenaX.mat, Contains data for all animals in arena X (with X=1-4)
                        	(b)	MovieName_DyePatterns.mat, Contains information about the measured microfluidic flow parameters
                        	Note that additional mat files are saved after each processing stage.

2)	  Name: 			        CorrectFlowDelaysForPulseDevices_02(File, ExperimentInformationFunctionName)
    	Inputs:             (a)	‘ExperimentInformationFunctionName’ is a string with the name of the external function storing the relevant experiment-specific information, including valve switching protocol and arena-specific animals’ identity. Examples: ‘Information_Exp1’, ‘Information_Exp2’
                        	(b)	The function loads the mat files from the previous section:  ‘MovieName_DataMatrices_ArenaX.mat’ (with X=1-4) and ‘MovieName_DyePatterns.mat’ 
    	External functions: ‘Information_Exp1.m’ and ‘‘Information_Exp2.m’ with experiment-specific information
    	Description:        Loads experiment-specific information and corrects for flow delays, namely the time from valve switching until the new odor reaches the animal nose. This way all animals are aligned with stimulus initiation regardless of their arena identity and location within the arena.
    	Outputs:            This function saves for each arena (X) files with delay corrected data:  
                        	(a)	'MovieName_DataMatrices_ArenaDelayCorrectedX.mat' with delay corrections based only on arena identity
                        	(b)	‘MovieName_DataMatrices_AllDelaysCorrectedX.mat’ with delay corrections based on arena identity and the animal nose location within the arena

3)	  Name:               ProcessData_Figs_4D_4E 
    	Inputs:          		This function loads the following mat files
                        	(a)	‘MovieName_DataMatrices_AllDelaysCorrected1-4.mat’ containing all data after correction of flow delays.
                        	(b)	‘MovieName_DyePatterns.mat’ with Dye dynamics
    	Description:        Extract behavioral responses to fast and slow decrease in odor concentration, analysis specific to Figures 4D and 4E. 
    	Output:     	      This function saves the file:  'BehavioralResponses_Figs_4D_4E.mat'     

4)	  Name:               PlotsAndStats_Figs_4D_4E 
    	Inputs:   	        Loads processed data from:  ‘BehavioralResponses_Figs_4D_4E.mat’  
    	Description:        Plot and calculate statistics of behavioral responses to fast and slow decrease in odor concentration, related to Figures 4D and 4E. 



Library 4: Theoretical analysis of noise filtering capacity, Related to Figure 6 and S6

• Theoretical analysis 
  • Name:     NoiseFilteringAnalysis_Fig6.m
  • Output:   NoiseAnalysis_Fig6.mat  
  • Notes:    The function includes SNR and speed calculations for each model over various noise and input conditions. Files with intermediate analysis are also available (see inline commentary).
• Plots 
  • Name:     PlotNoiseFilteringAnalysis_Fig6.m
  • Input:    NoiseAnalysis_Fig6.mat 

Library 5: Zebrafish OT calcium responses to looming visual stimuli, Related to Figure 7 and S7

• Data processing 
  • Name:     ProcessZebrafishOT_LoomingResponse_Fig7.m
  • Input:    Data_ZebrafishImaging_Raw.mat
  • Output:   Data_ZebrafishImaging_Processed.mat 	
• Plots 
  • Name:     PlotZebrafishOT_LoomingResponse_Fig7.m
  • Input:    Data_ZebrafishImaging_Processed.mat 	



                         

