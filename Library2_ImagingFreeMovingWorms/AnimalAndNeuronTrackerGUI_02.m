function AnimalAndNeuronTrackerGUI_02(File, BehaviorDataMatrices, SwitchPathDefinitionInTracks)
% This function use the external function:   Tracks = ExtractNeuralActivityFromTracks_AfterGUI_02(Tracks, NumOfPixelsInNeuron, plotme); 
% SwitchPathDefinitionInTracks = false; % If using the Tracks structure saved before behavior segmentation  
% SwitchPathDefinitionInTracks = true;  % If using the Tracks structure saved after behavior segmentation (X-Y definition switch) 

%% Initialization and load data 
% Initialization
ArenaID                      = 1; 
SpecificFrames               = [];
if ~exist('BehaviorDataMatrices','var')
    BehaviorDataMatrices = [];
    MovieOptions.ShowBehaviorFrom_Data_Matrices = false;
else
    MovieOptions.ShowBehaviorFrom_Data_Matrices = true;
end

% Load Tracks data
FileName                     = [File.TrackFile(1:end-4),'_NeuronsDataBeforeManualCorrection_Arena',num2str(ArenaID),'.mat'];  
load(FileName,'File','Tracks','NeuronDisplayMatrix','NeuronCoordinatesMatrix','HeadIsSqueezedMAT','background','VignettingPattern','WormArea');

if exist('SwitchPathDefinitionInTracks','var')&& SwitchPathDefinitionInTracks
    for tr=1:length(Tracks)
        Tracks(tr).Path = Tracks(tr).Path(:,[2 1]);
    end
end

%  Compute mask contours
load(File.BackgroundFile,'Mask');
PROPS = regionprops(~Mask,'Area','PixelIdxList');
All_Mask_Countours = [];
for obj_ind = 1:length([PROPS.Area])
    BW                             = false(size(Mask));
    BW(PROPS(obj_ind).PixelIdxList)= true;
    [row, col]                     = find(BW,1,'first');    
    CONTOUR                        = bwtraceboundary(BW, [row, col],'E');
    PROPS(obj_ind).CONTOUR         = CONTOUR;
    All_Mask_Countours = [All_Mask_Countours; CONTOUR; [NaN NaN]];
end
% figure; 
% imshow(VignettingPattern,[]); hold on;
% plot(All_Mask_Countours(:,2),All_Mask_Countours(:,1),'linewidth',2);    

%% Movie options
MovieOptions.All_Mask_Countours = All_Mask_Countours;  % Display arena Mask. Comment line for avoid display  
if exist('ArenaID','var')
    MovieOptions.arena                      = ArenaID;     
end
% MovieOptions.EnlargeArena                 = 200;     
% MovieOptions.ShowEndOfTracks.FramesBefore = 20;        % If interested only in the collision
% MovieOptions.ShowEndOfTracks.FramesAfter  = 20;        % If interested only in the collision
if exist('SpecificFrames','var')
    if ~isempty(SpecificFrames)
        MovieOptions.ShowSpecificFrames = SpecificFrames;
    end
end

MovieOptions.TracksLine                   = false;      % true/false/Num = show/(don't show)/(Show last Num of frames) of track line preceeding to each frame.                                                         
MovieOptions.TracksText.TrackNum          = true;      % Track text diplay settings
MovieOptions.TracksText.Behavior          = false;            % Possible values: 'LowLevel'/'HighLevel'/false.  default: false
% MovieOptions.TracksText.ColorByID       = 0;     % No, same colors for all worms. different colors based on features 
MovieOptions.TracksText.ColorByID         = 2;     % Color by track index
% MovieOptions.TracksText.ColorByID       = 3;     % Color by Chemotaxis index
% MovieOptions.TracksText.ColorByID       = 4;     % Color by behavior code

MovieOptions.ShowDetectedPixels           = {'Midline','NeuronCoordinates','Perimeter'};                 % A cell with the name of fields to display
MovieOptions.ShowDetectedPixelsColors     = {'g','r','m'};                                               % Color for display (for ColorByID=0)
% MovieOptions.ShowDetectedPixels         = {'Perimeter','Midline','Head','Tail','NeuronCoordinates'}; % A cell with the name of fields to display
% MovieOptions.ShowDetectedPixelsColors   = {'b','g','r','m','y'};                                     % Color for display (for ColorByID=0)

MovieOptions.ConstantScale   = true;        % If true, the function will look for the best scaling and scale each frame similarly. If false --> imshow(Mov,[]) 
MovieOptions.ConstantScale   = 2;           % If numeric, the function will look for the best scaling and scale each frame similarly but will use only 'ConstantScale' of the range[%]. 

% MovieOptions.ShowValueOfThisField         = 'CTX_AngularVelocity_during_Runs';          % If exist and not empty then the value of this field is shown instead of CTX_position              
% MovieOptions.ShowValueOfThisField         = 'Speed_during_Runs';                        % If exist and not empty then the value of this field is shown instead of CTX_position              
MovieOptions.ShowCTXposition              = false;                                        % If true, 'BehaviorDataMatrices' must be added as input to this function              

%%% For behavior segmentation display, if needed
BehaviorStructure = [];
if isfield(MovieOptions,'TracksText')
    if isfield(MovieOptions.TracksText,'Behavior')
        if MovieOptions.TracksText.Behavior
            fieldname                             = MovieOptions.TracksText.Behavior;
            BehaviorStructure.fieldname_inTracks  = ['BehaviorCode_',MovieOptions.TracksText.Behavior];
            BehaviorStructure.BehaviorCodeNumbers = Tracks(1).BehaviorCode_Structure.(fieldname).BehaviorCodeNumbers; 
            BehaviorStructure.BehaviorCodeName    = Tracks(1).BehaviorCode_Structure.(fieldname).BehaviorCodeName;    
            if strcmpi(fieldname,'LowLevel')       %  for:    'Out of bounds'    'Forward'    'Curve'    'Pause'    'Reverse'    'Omega'    'SharpTurn'
                BehaviorStructure.BehaviorCodeColors  = {'y','b','b','k','r','g','r'};    
            elseif strcmpi(fieldname,'HighLevel')   %  for:    'Out of bounds', 'Forward', 'Curve', 'Pause', 'Reversal', 'Omega-Pause', 'SharpTurn', 'Pirouette. Initial reversal', 'Pirouette. After reversal'  
                BehaviorStructure.BehaviorCodeColors  = {'y','b','b','k','r','g','r','m','m'}; 
            end
        end
    end
end            


%% Run GUI
RunGUI(Tracks, File, MovieOptions, BehaviorDataMatrices, BehaviorStructure, ...
       NeuronDisplayMatrix, NeuronCoordinatesMatrix, HeadIsSqueezedMAT, background, VignettingPattern, WormArea); 

   
return

function RunGUI (Tracks, File, MovieOptions, BehaviorDataMatrices, BehaviorStructure, ...
    NeuronDisplayMatrices, NeuronCoordinatesMatrix, HeadIsSqueezedMAT, background, VignettingPattern, WormArea)
% Display a movie with Tracks 
% Information about the Tracks and the movie is given in the Tracks and File input variables  
% NeuronDisplayMatrices   >> (tr_ind,frame,:,:) 
% NeuronCoordinatesMatrix >> (tr_ind,frame,2) 
% MovieOptions-  A stucture with infomation about what and how to display, with the following fields:   
%    arena:              Zoom on one specific arena (Default = no Zoom).       
%                        This field is NOT used if specific tracks are shown (arena zoom will be set by the tracks) 
%    EnlargeArena:       Enlarge Arena zoom by the specified number of pixels (Default = no enlargement).       
%    ShowSpecificTracks: A vector with Tracks of interest to display. Frames will be set by Tracks (Default = Show all tracks within the frames of interest).   
%                        Other Tracks within the same region of interest will be shown as well. 
%    ShowEndOfTracks:    A structure that, if exists, specifiy how many frames to display before the end of the Track (Default= show all the track)
%                        This is convenient to see worms before and after collisions  
%         FramesBefore   Number of frames before the first track ends (before first collision).  
%         FramesAfter    Number of frames after the start of the last detected track (after last collision).            
%    ShowSpecificFrames: A vector of Frames of interest to display. All relevant tracks will be shown(Default = All frames in the movie). 
%    TracksText:         Structure defining what fields to display. Values are true or false (Default = show only track number). 
%         TrackNum       if exists and true, show the track number (#)
%         Behavior      Worm behavior as calculated in the segmentation functions.   
% MovieOptions.TracksLine

%% Movie format
MoviePath      = File.MoviePath;
MovieFileNames = File.MovieFileNames;
if     strcmpi(File.VideoFormat,'tiff')
    VideoFormat   = 1;                                  % read new file every loop iteration 
elseif strcmpi(File.VideoFormat,'multi-tiff')
    VideoFormat   = 2;
    FileFullName  = [MoviePath,'\',MovieFileNames{1}];  % one file
elseif strcmpi(File.VideoFormat,'avi')
    VideoFormat   = 3;
    FileFullName  = [MoviePath,'\',MovieFileNames{1}];  % one file
    MovieObj      = VideoReader(FileFullName);          % one object
elseif strcmpi(File.VideoFormat,'multi-avi')
    VideoFormat   = 4;
    PreviousMovieObjFile = 0;                           % read new file object if and when it is needed 
end

MovieOptions.MinMaxNeuronDisplay = prctile(NeuronDisplayMatrices(:),[0.1 99.9]);

%% Find the frames of each track
NumOfTracks     = length(Tracks);
NumOfFrames     = File.NumberOfFrames;
FirstFrames     = zeros(1,NumOfTracks);
LastFrames      = zeros(1,NumOfTracks);
ALLFRAMES       = false(NumOfFrames,NumOfTracks);
ConvertMovieFrameToIndexWithinTracks = zeros(NumOfTracks,NumOfFrames,'single')*NaN;
for tr=1:NumOfTracks
    CurrentFrames               = Tracks(tr).Frames;
    TotalNumOfFramesInTrack     = length(CurrentFrames);
    FirstFrames(tr)             = CurrentFrames(1);
    LastFrames(tr)              = CurrentFrames(end);
    ALLFRAMES(CurrentFrames,tr) = 1;
    ConvertMovieFrameToIndexWithinTracks(tr,CurrentFrames) = 1:TotalNumOfFramesInTrack;
end

%% Set MovieOptions parameters   
ShowTracksText = false;
if isfield(MovieOptions,'TracksText')
    field_names    = fieldnames(MovieOptions.TracksText);
    field_names    = field_names(~strcmpi(field_names,'ColorByID'));
    for f_ind = 1:length(field_names)
        if islogical(MovieOptions.TracksText.(field_names{f_ind}))
            ShowTracksText = ShowTracksText || MovieOptions.TracksText.(field_names{f_ind});  % at least one text field = true
        else
            ShowTracksText = true;
        end
    end          
end
MovieOptions.ShowTracksText = ShowTracksText;  

% Zoom on a specific arena
if ~isfield(MovieOptions,'arena')          
    MovieOptions.arena = NaN;    
end       
% This field exists if user is interested only in Frames. All relevant Tracks will be shown.  
if isfield(MovieOptions,'ShowSpecificFrames')          
    Frames2Show = MovieOptions.ShowSpecificFrames; 
else
    All_Frames = unique([Tracks.Frames]);
    Frames2Show = All_Frames;    
end
MovieOptions.Frames2Show = Frames2Show;    
% The relevant arena Zoom area: 
if ~isnan(MovieOptions.arena)
    ArenaBox =  File.TrackBoxAxis(MovieOptions.arena,:);
    if isfield(MovieOptions,'EnlargeArena')
        Enlarge = MovieOptions.EnlargeArena;
        ArenaBox = [ArenaBox(1)-Enlarge, ArenaBox(2)+Enlarge, ArenaBox(3)-Enlarge, ArenaBox(4)+Enlarge];
    end
    MovieOptions.ArenaBox = ArenaBox;
end    
if ~isfield(MovieOptions,'TracksLine')  
    MovieOptions.TracksLine = false;
end      
if isfield(MovieOptions,'ShowDetectedPixels')  
    MovieOptions.ShowWormFeatures   = true;
    MovieOptions.WormFeatures       = MovieOptions.ShowDetectedPixels;        % A cell with the name of ields to display
    MovieOptions.WormFeaturesCOLORS = MovieOptions.ShowDetectedPixelsColors;  % Color for display
else
    MovieOptions.ShowWormFeatures = false;
end    
if ~isfield(MovieOptions,'ConstantScale')
    MovieOptions.ConstantScale = false;    
end
ConstantScale = MovieOptions.ConstantScale;
MovieOptions.BestScale = [];
if ConstantScale
    Max = 0;
    Min = inf;
    for Frame = Frames2Show(round(1:end/20:end))
        switch VideoFormat
        case 2                             % A single multiple-tiff file
            Mov           = imread(FileFullName, Frame);     
        case 1                             % A sequence of tiff files
            FileFullName  = [MoviePath,'\',MovieFileNames{Frame}]; 
            Mov           = imread(FileFullName);                          
        case 3                             % A single avi movie file
            Mov           = read(MovieObj, Frame);     
            Mov           = Mov(:,:,1);                %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!
        case 4                             % A sequence of avi movie files
            CurrentMovieObjFile = File.MultiAviFrameConversion.MovieFileNumber(Frame);
            if (PreviousMovieObjFile ~= CurrentMovieObjFile)        % create a new file object if necessary
                FileFullName        = [MoviePath,'\',MovieFileNames{CurrentMovieObjFile}]; 
                PreviousMovieObjFile = CurrentMovieObjFile;
                MovieObj             = VideoReader(FileFullName);
            end
            CurrentFrameNumberInFile = File.MultiAviFrameConversion.FrameNumberInFile(Frame);
            Mov                      = read(MovieObj, CurrentFrameNumberInFile);     
            Mov                      = Mov(:,:,1);          %%%% assuming 3 channel movie !!!!!!!!!!!!!!!!!!                                
        end 

        Max = max([Max, max(Mov(:))]);
        Min = min([Min, min(Mov(:))]);
    end
    BestScale = [Min Max];    
    if isnumeric(ConstantScale)   
        if (ConstantScale<=100)&&(ConstantScale>=0)   % Scale based on percentage given by the user
            disp(['re-setting scale (',num2str(ConstantScale),'%)']);
            BestScale = prctile(double(Mov(:)),[ConstantScale 100]);
%            BestScale = prctile(double(Mov(:)),[ConstantScale 100-ConstantScale]);
        end 
    end
    MovieOptions.BestScale = BestScale;
end
ScaleMin = BestScale(1);
ScaleMax = BestScale(2);
          
if exist('BehaviorStructure','var')
    if ~isempty(BehaviorStructure)
        BehaviorFieldName = BehaviorStructure.fieldname_inTracks;
        MovieOptions.BehaviorFieldName = BehaviorFieldName;
    end
end    
Add_CTXdata_And_SmoothPosition = false;
if isfield(MovieOptions,'ShowCTXposition')   
    if MovieOptions.ShowCTXposition
        if ~isempty(Data_Matrices)
            Add_CTXdata_And_SmoothPosition = true;
        end
    end
end    
if isfield(MovieOptions,'ShowValueOfThisField')   
    if ~isempty(MovieOptions.ShowValueOfThisField)
        if ~isempty(Data_Matrices)
            Add_CTXdata_And_SmoothPosition = true;
            MovieOptions.FramesField = 'FrameNumber';
            MovieOptions.X_field     = 'Coordinates_X_Smoothed';
            MovieOptions.Y_field     = 'Coordinates_Y_Smoothed';
            MovieOptions.CTX_field   = MovieOptions.ShowValueOfThisField;                
            disp(['Showing values of ', MovieOptions.CTX_field]); 
         end
    end
end    
MovieOptions.Add_CTXdata_And_SmoothPosition = Add_CTXdata_And_SmoothPosition;

if MovieOptions.ShowBehaviorFrom_Data_Matrices
    if ~isempty(Data_Matrices)
        if ~isfield(MovieOptions,'CTX_field') % if it was not previously defined in 'MovieOptions.ShowValueOfThisField' then show 'CTX_position'
            if isfield(Data_Matrices,'FrameNumber_Smoothed2')
                MovieOptions.FramesField = 'FrameNumber_Smoothed2';
                MovieOptions.X_field     = 'X_Smoothed2';
                MovieOptions.Y_field     = 'Y_Smoothed2';
                MovieOptions.CTX_field   = 'CTX_position_Smoothed2';
            elseif isfield(Data_Matrices,'FrameNumber')
                disp('Showing ''raw'' CTX data --> The shown CTX  value are not yet corrected to the butanone odor side at that experiment')
                MovieOptions.FramesField = 'FrameNumber';
                MovieOptions.X_field     = 'Coordinates_X_Smoothed';
                MovieOptions.Y_field     = 'Coordinates_Y_Smoothed';
                MovieOptions.CTX_field   = 'CTX_position';
            end    
        end
    end                      
end
if MovieOptions.TracksText.ColorByID
    COLORS_FOR_TRACKS   = {'w','g','c','m','b','y',[0.8 0.5 0.2],[0.8 0.5 1]}; % it will start from the green
    if MovieOptions.TracksText.ColorByID==2
        MovieOptions.TracksColorsIndices = mod(1:length(Tracks),length(COLORS_FOR_TRACKS))+1;  
    elseif MovieOptions.TracksText.ColorByID==1
        MovieOptions.TracksColorsIndices = mod(1:length(unique([Tracks.ID])),length(COLORS_FOR_TRACKS))+1;  
    end    
end

%% Figure and data initialization
Fig_handle=figure('position',get(0,'screensize'),'name','GUI for manual inspection of neuron segmentation');  
% data structure associated with the GUI
data.NeuronDisplayMatrices                      = NeuronDisplayMatrices;
data.NeuronDisplayMatrices_Corrected            = NeuronDisplayMatrices;
data.NeuronDisplayMatrices_ManuallyCorrected    = NeuronDisplayMatrices;
data.NeuronCoordinatesMatrix                    = NeuronCoordinatesMatrix;
data.NeuronCoordinatesMatrix_Corrected          = NeuronCoordinatesMatrix;
data.NeuronCoordinatesMatrix_ManuallyCorrected  = NeuronCoordinatesMatrix;
data.ManuallyAssignedAsCorrect                  = false(NumOfTracks,NumOfFrames); % Will be considered correct for correcting other frames
data.ManuallyAssignedAsDisregard                = false(NumOfTracks,NumOfFrames); % NaN
data.ALLFRAMES                                  = ALLFRAMES;
data.ConvertMovieFrameToIndexWithinTracks       = ConvertMovieFrameToIndexWithinTracks;
data.FramesForNeuronPositionReCalculation       = false(NumOfTracks,NumOfFrames); 

%% Setting GUI control buttons
Panel_handle = uipanel('Parent',Fig_handle,'Title','Settings','Fontweight','bold','FontSize',12,'Position',[.77 .1 .2 .85]);
HANDLES.Fig_handle   = Fig_handle;
HANDLES.Panel_handle = Panel_handle;
% 'position' --> [left bottom width height]
left                  = 15;
text_width            = 250;
text_width_long       = 300;
text_height           = 20;  % also for 'edit'
edit_width            = 50;
interval_between_rows = 30; 
Frame                 = Frames2Show(1);

bottom  = 800;
HANDLES.h_movie_path_text    = uicontrol('Parent',Panel_handle,'style','text','String','Movie path:',  'HorizontalAlignment','left','Position',[left             bottom    text_width   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_movie_path_text    = uicontrol('Parent',Panel_handle,'style','text','String',File.MoviePath, 'HorizontalAlignment','left','Position',[left  bottom    text_width_long   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_movie_res_text    = uicontrol('Parent',Panel_handle,'style','text','String','Movie frame rate [Hz]','HorizontalAlignment','left','Position',[left             bottom    text_width   text_height]);
HANDLES.h_movie_res_text    = uicontrol('Parent',Panel_handle,'style','text','String',num2str(File.FrameRate),     'HorizontalAlignment','left','Position',[left+text_width  bottom    edit_width   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_ZoomOnWorm_text     = uicontrol('Parent',Panel_handle,'style','text','String','Zoom on Worm?', 'HorizontalAlignment','left','Position',[left              bottom    text_width      text_height]);
HANDLES.h_ZoomOnWorm_ONOFF    = uicontrol('Parent',Panel_handle,'style','popup','String','ON|OFF','Value',2,'HorizontalAlignment','left','Position',[left+text_width/2   bottom    edit_width+10   text_height]);
TracksNumCell = cell(1,NumOfTracks); for tr=1:NumOfTracks, TracksNumCell{tr} = num2str(tr); end
HANDLES.h_ZoomOnWorm_tracknum = uicontrol('Parent',Panel_handle,'style','popup','String',TracksNumCell,'HorizontalAlignment','left','Position',[left+text_width   bottom    edit_width+10   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_CurrentFrame_text   = uicontrol('Parent',Panel_handle,'style','text','String','Frame', 'HorizontalAlignment','left','Position',[left              bottom    text_width      text_height]);
HANDLES.h_CurrentFrame        = uicontrol('Parent',Panel_handle,'style','edit','String',num2str(Frame),'Value',Frame,'HorizontalAlignment','left','Position',[left+text_width   bottom    edit_width+10   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_FrameInterval_text   = uicontrol('Parent',Panel_handle,'style','text','String','Frame Interval', 'HorizontalAlignment','left','Position',[left              bottom    text_width      text_height]);
HANDLES.h_FrameInterval        = uicontrol('Parent',Panel_handle,'style','edit','String',num2str(10),'Value',10,'HorizontalAlignment','left','Position',[left+text_width   bottom    edit_width+10   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_TimeBwFrames_text   = uicontrol('Parent',Panel_handle,'style','text','String','Time between frames (sec)', 'HorizontalAlignment','left','Position',[left              bottom    text_width      text_height]);
HANDLES.h_TimeBwFrames        = uicontrol('Parent',Panel_handle,'style','edit','String',num2str(0.1),'Value',0.1,'HorizontalAlignment','left','Position',[left+text_width   bottom    edit_width+10   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_StopMovie           = uicontrol('Parent',Panel_handle,'style','pushbutton','String','Stop','FontWeight','bold','Value',2,'BackgroundColor','r',...
                                          'HorizontalAlignment','left','Position',[left+edit_width*1.5    bottom    edit_width   text_height]);
HANDLES.h_ShowSingleFrame     = uicontrol('Parent',Panel_handle,'style','pushbutton','String','Show single frame','FontWeight','bold','Value',2,...
                                                'Callback',{@ShowSingleFrame,Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions},...
                                                'HorizontalAlignment','left','Position',[left+edit_width*4    bottom    edit_width*3   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_PlayMovie           = uicontrol('Parent',Panel_handle,'style','pushbutton','String','Play','FontWeight','bold','Value',2,...
                                          'Callback',{@PlayMovie,Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions},'BackgroundColor','g',...
                                          'HorizontalAlignment','left','Position',[left+edit_width*1.5    bottom    edit_width   text_height]);
HANDLES.h_ShowNextFrame       = uicontrol('Parent',Panel_handle,'style','pushbutton','String','Show next frame','FontWeight','bold','Value',2,...
                                                'Callback',{@ShowNextFrame,Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions},...
                                                'HorizontalAlignment','left','Position',[left+edit_width*4    bottom    edit_width*3   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_PlayMovieBackwards  = uicontrol('Parent',Panel_handle,'style','pushbutton','String','Play back','FontWeight','bold','Value',2,...
                                          'Callback',{@PlayMovieBackwards,Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions},'BackgroundColor','g',...
                                          'HorizontalAlignment','left','Position',[left+edit_width*1.5    bottom    edit_width*2   text_height]);
HANDLES.h_ShowPreviousFrame     = uicontrol('Parent',Panel_handle,'style','pushbutton','String','Show previous frame','FontWeight','bold','Value',2,...
                                                'Callback',{@ShowPreviousFrame,Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions},...
                                                'HorizontalAlignment','left','Position',[left+edit_width*4    bottom    edit_width*3   text_height]);
bottom  = bottom-interval_between_rows;
bottom  = bottom-interval_between_rows;
HANDLES.h_ScaleRange_text   = uicontrol('Parent',Panel_handle,'style','text','String','Scale range (min,max)', 'HorizontalAlignment','left','Position',          [left                bottom    text_width      text_height]);
HANDLES.h_ScaleRange_Min    = uicontrol('Parent',Panel_handle,'style','edit','String',num2str(ScaleMin),'Value',ScaleMin,'HorizontalAlignment','left','Position',[left+text_width*2/3 bottom    edit_width+10   text_height]);
HANDLES.h_ScaleRange_Max    = uicontrol('Parent',Panel_handle,'style','edit','String',num2str(ScaleMax),'Value',ScaleMax,'HorizontalAlignment','left','Position',[left+text_width     bottom    edit_width+10   text_height]);

bottom  = bottom-interval_between_rows;
bottom  = bottom-interval_between_rows;
HANDLES.h_SaveSessionBackup   = uicontrol('Parent',Panel_handle,'style','pushbutton','String','Save session backup','FontWeight','bold','Value',2,...
                                               'Callback',{@SaveSessionBackup,Fig_handle, File},'BackgroundColor',[1 1 1],...
                                               'HorizontalAlignment','left','Position',[left    bottom    text_width*0.6   text_height]);
HANDLES.h_LoadPreviousSession = uicontrol('Parent',Panel_handle,'style','pushbutton','String','Load previous session','FontWeight','bold','Value',2,...
                                               'Callback',{@LoadPreviousSession,Fig_handle, File},'BackgroundColor',[1 1 1],...
                                               'HorizontalAlignment','left','Position',[left+text_width*0.8    bottom    text_width*0.6   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_UpdateAndSaveTracksFile   = uicontrol('Parent',Panel_handle,'style','pushbutton','String','Finish session: update and save Tracks file','FontWeight','bold','Value',2,...
           'Callback',{@UpdateAndSaveTracksFile,Fig_handle, File, Tracks, BehaviorDataMatrices, background, VignettingPattern, WormArea},...
           'BackgroundColor',[0.7 0.7 1],'HorizontalAlignment','left','Position',[left+30    bottom    text_width   text_height]);

       
bottom  = bottom-interval_between_rows;
bottom  = bottom-interval_between_rows;
bottom  = bottom-interval_between_rows;
bottom  = bottom-interval_between_rows;
HANDLES.h_CreateMovie   = uicontrol('Parent',Panel_handle,'style','pushbutton','String','Create Movie','FontWeight','bold','Value',2,...
           'Callback',{@CreateMovie,Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions},...
           'BackgroundColor',[1 1 1],'HorizontalAlignment','left','Position',[left+30    bottom    text_width   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_MovieLastFrame_text   = uicontrol('Parent',Panel_handle,'style','text','String','Movie last frame', 'HorizontalAlignment','left','Position',[left              bottom    text_width/2      text_height]);
HANDLES.h_MovieLastFrame        = uicontrol('Parent',Panel_handle,'style','edit','String',num2str(10),'Value',10,'HorizontalAlignment','left','Position',[left+text_width/2   bottom    edit_width+10   text_height]);
bottom  = bottom-interval_between_rows;
HANDLES.h_MovieDisplayRate_text   = uicontrol('Parent',Panel_handle,'style','text','String','Movie display rate', 'HorizontalAlignment','left','Position',[left              bottom    text_width/2      text_height]);
HANDLES.h_MovieDisplayRate        = uicontrol('Parent',Panel_handle,'style','edit','String',num2str(10),'Value',10,'HorizontalAlignment','left','Position',[left+text_width/2   bottom    edit_width+10   text_height]);
bottom  = bottom-interval_between_rows;
MovieName = [File.MoviePath,'\',File.MovieName,'_Movie1'];
HANDLES.h_MovieName_text   = uicontrol('Parent',Panel_handle,'style','text','String','Movie name', 'HorizontalAlignment','left','Position',[left              bottom    text_width/2      text_height]);
HANDLES.h_MovieName        = uicontrol('Parent',Panel_handle,'style','edit','String',MovieName,'Value',0.1,'HorizontalAlignment','left','Position',[left+text_width/2   bottom    text_width   text_height]);

       
%% Display first frame and identify possible segmentation inaccuracies for manual inspection      
data.HANDLES = HANDLES;
guidata(Fig_handle,data);   % first assignment of structure to the GUI
DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions);
FindPossibleErrorSegmentation(Fig_handle, File,Tracks,HeadIsSqueezedMAT);
data = guidata(Fig_handle); % get data structure from GUI

%% Buttons for manual correction
figure('position',get(0,'ScreenSize')); 
Subplot_Positions = cell(1,NumOfTracks);
for t_ind=1:length(Tracks)
    subplot(length(Tracks),3,3*t_ind);
    position = get(gca,'position'); position(1) = position(1)-0.32;
    Subplot_Positions{t_ind}=position;
end
close

panel_left   = Subplot_Positions{1}(1)+Subplot_Positions{1}(3)*0.85;
panel_height = Subplot_Positions{1}(4)*1.1;
left2                  = 5;
text_width2            = 200;
for tr=1:length(Tracks)
    if length(Tracks)==3
        bottom_WithinPanel = 195;
        interval_between_rows2 = 21;
        interval_between_rows3 = 15;
    else
        bottom_WithinPanel = 300;
        interval_between_rows2 = 31;
        interval_between_rows3 = 22;
    end        
    panel_bottom    = Subplot_Positions{tr}(2);
    Panel_handle_SP(tr) = uipanel('Parent',Fig_handle,'Title',['Segmentation of track ',num2str(tr)],'Fontweight','bold','FontSize',12,...
        'Position',[panel_left panel_bottom .215 panel_height]);    
    % 
    h_PlayMovieForCorrection(tr) = uicontrol('Parent',Panel_handle_SP(tr),'style','pushbutton','String','Play next possible error','FontWeight','bold','Value',2,...
                                               'Callback',{@PlayMovieForCorrection,Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions, tr},...
                                               'HorizontalAlignment','left','Position',[left2    bottom_WithinPanel    text_width2*0.8    text_height]);
    if ~isempty(data.SegmentsToCorrect(tr).InitialFrame)
        FrameIntervalText = ['Relevant frame interval= [',num2str(data.SegmentsToCorrect(tr).InitialFrame(1)),' - ',num2str(data.SegmentsToCorrect(tr).EndFrame(1)),']'];
    else
        FrameIntervalText = 'No flagged frames';
    end
    h_CorrectionFrameInterval_text(tr)= uicontrol('Parent',Panel_handle_SP(tr),'style','text','String',FrameIntervalText, ...
                                        'ForegroundColor','k','HorizontalAlignment','left','Position',[left2+text_width2 bottom_WithinPanel  text_width2  text_height]);    
    bottom_WithinPanel  = bottom_WithinPanel-interval_between_rows2;      
    h_ConfirmProperSegmentation(tr) = uicontrol('Parent',Panel_handle_SP(tr),'style','pushbutton','BackgroundColor','g',...
                                               'String','Confirm proper segmentation within frame interval','FontWeight','bold','Value',2,...
                                               'Callback',{@ConfirmProperSegmentation,Fig_handle, tr},...
                                               'HorizontalAlignment','left','Position',[left2+25    bottom_WithinPanel    text_width*1.2    text_height]);        
    bottom_WithinPanel  = bottom_WithinPanel-interval_between_rows2-12;  
    
    % Input new coordinates from zoomed image or frame image  
                              uicontrol('Parent',Panel_handle_SP(tr),'style','text','String','Manual position correction','FontWeight','bold', ...
                                        'ForegroundColor','k','HorizontalAlignment','left','Position',[left2 bottom_WithinPanel  text_width2  text_height]);    
    bottom_WithinPanel  = bottom_WithinPanel-interval_between_rows3+5;    
    uicontrol('Parent',Panel_handle_SP(tr),'style','text','String','Incorrect. Input new coordinates:', ...
                                        'ForegroundColor','r','HorizontalAlignment','left','Position',[left2+text_width2 bottom_WithinPanel  text_width2  text_height]);    
    bottom_WithinPanel  = bottom_WithinPanel-interval_between_rows3;    
    h_DefineFrameAsCorrect(tr)   = uicontrol('Parent',Panel_handle_SP(tr),'style','pushbutton','String','Correct','FontWeight','bold','Value',2,...
                                               'Callback',{@DefineFrameAsCorrect,Fig_handle, tr},'BackgroundColor','g',...
                                               'HorizontalAlignment','left','Position',[left2  bottom_WithinPanel  text_width2/2  text_height]);
    h_DefineFrameAsDisregard(tr) = uicontrol('Parent',Panel_handle_SP(tr),'style','pushbutton','String','Disregard','FontWeight','bold','Value',2,...
                                               'Callback',{@DefineFrameAsDisregard,Fig_handle, tr, File, Tracks, BehaviorDataMatrices, MovieOptions},'BackgroundColor',[0.7 0.7 0.7],...
                                               'HorizontalAlignment','left','Position',[left2+text_width2/2  bottom_WithinPanel  text_width2/2  text_height]);                                          
    h_InputNeuronCoordinatesFromZoomedImage(tr) = ...
     uicontrol('Parent',Panel_handle_SP(tr),'style','pushbutton','String','using zoomed image','FontWeight','bold','Value',2,...
               'Callback',{@InputNeuronCoordinatesFromZoomedImage,Fig_handle, tr, File, Tracks, BehaviorDataMatrices, MovieOptions,background, VignettingPattern, WormArea},...
               'BackgroundColor',[1 0.3 0.3],'HorizontalAlignment','left','Position',[left2+text_width2  bottom_WithinPanel  text_width2  text_height]);                                                                                 
    bottom_WithinPanel  = bottom_WithinPanel-interval_between_rows2;    
    h_InputNeuronCoordinatesFromFrameImage(tr) = ...
     uicontrol('Parent',Panel_handle_SP(tr),'style','pushbutton','String','using full image','FontWeight','bold','Value',2,...
               'Callback',{@InputNeuronCoordinatesFromFrameImage,Fig_handle, tr, File, Tracks, BehaviorDataMatrices, MovieOptions,background, VignettingPattern, WormArea},...
               'BackgroundColor',[1 0.3 0.3],'HorizontalAlignment','left','Position',[left2+text_width2  bottom_WithinPanel  text_width2  text_height]);                                                                                 
                                     
    % correction frame interval
    bottom_WithinPanel = bottom_WithinPanel-interval_between_rows2-15;    
       uicontrol('Parent',Panel_handle_SP(tr),'style','text','String','Correct position within frame interval:','FontWeight','bold','HorizontalAlignment','left','Position',[left2              bottom_WithinPanel    text_width      text_height]);
    bottom_WithinPanel = bottom_WithinPanel-interval_between_rows3;    
                              uicontrol('Parent',Panel_handle_SP(tr),'style','text','String','First frame', 'HorizontalAlignment','left','Position',[left2              bottom_WithinPanel    text_width2/2      text_height]);
    h_Correct_FirstFrame(tr)= uicontrol('Parent',Panel_handle_SP(tr),'style','edit','String','0','Value',0,'HorizontalAlignment','left','Position',[left2+text_width2/3   bottom_WithinPanel    edit_width+10   text_height]);
    bottom_WithinPanel = bottom_WithinPanel-interval_between_rows3;    
    	                      uicontrol('Parent',Panel_handle_SP(tr),'style','text','String','Last frame', 'HorizontalAlignment','left','Position',[left2              bottom_WithinPanel    text_width2/2      text_height]);
    h_Correct_LastFrame(tr) = uicontrol('Parent',Panel_handle_SP(tr),'style','edit','String','0','Value',0,'HorizontalAlignment','left','Position',[left2+text_width2/3   bottom_WithinPanel    edit_width+10   text_height]);
    h_ClearNonConfirmedCorrections(tr) = ...
                              uicontrol('Parent',Panel_handle_SP(tr),'style','pushbutton','String','Clear non-confirmed corrections','FontWeight','bold','Value',2,...
                                               'Callback',{@ClearNonConfirmedCorrections,Fig_handle, tr},'BackgroundColor',[0.7 0.7 0.7],...
                                               'HorizontalAlignment','left','Position',[left2+text_width2*0.8  bottom_WithinPanel+interval_between_rows3/2  text_width2  text_height]);
    
    bottom_WithinPanel = bottom_WithinPanel-interval_between_rows2;    
    h_RecalculateNeuronPositionWithinSegment(tr) = ...
                              uicontrol('Parent',Panel_handle_SP(tr),'style','pushbutton','String','Re-calculate neuron position','FontWeight','bold','Value',2,...
                                               'Callback',{@RecalculateNeuronPositionWithinSegment,Fig_handle, tr, File, background, VignettingPattern, WormArea},'BackgroundColor',[1 1 1],...
                                               'HorizontalAlignment','left','Position',[left2  bottom_WithinPanel  text_width2  text_height]);
    h_DisregardNeuronPositionWithinSegment(tr) = ...
                              uicontrol('Parent',Panel_handle_SP(tr),'style','pushbutton','String','Desregard neuron position','FontWeight','bold','Value',2,...
                                               'Callback',{@DisregardNeuronPositionWithinSegment,Fig_handle, tr}, 'ForegroundColor','w','BackgroundColor',[0 0 0],...
                                               'HorizontalAlignment','left','Position',[left2+text_width2  bottom_WithinPanel  text_width2  text_height]);
                                     
    bottom_WithinPanel = bottom_WithinPanel-interval_between_rows2;    
    h_MaxPixelsPerFrame(tr)= uicontrol('Parent',Panel_handle_SP(tr),'style','edit','String','2','Value',2,'HorizontalAlignment','left','Position',[left2+text_width2/3   bottom_WithinPanel    edit_width+10   text_height]);

    h_PlayCurrentMovieInterval(tr) = ...
                              uicontrol('Parent',Panel_handle_SP(tr),'style','pushbutton','String','Play current interval','FontWeight','bold','Value',2,...
                                               'Callback',{@PlayCurrentMovieInterval,Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions, tr},'BackgroundColor',[0 1 0],...
                                               'HorizontalAlignment','left','Position',[left2+text_width2  bottom_WithinPanel  text_width2  text_height]);
        
end

%% Add HANDLES to data structure
HANDLES.TracksPanelHandles       = Panel_handle_SP;
HANDLES.h_DefineFrameAsCorrect   = h_DefineFrameAsCorrect;
HANDLES.h_DefineFrameAsDisregard = h_DefineFrameAsDisregard;
HANDLES.h_InputNeuronCoordinatesFromZoomedImage = h_InputNeuronCoordinatesFromZoomedImage;
HANDLES.h_InputNeuronCoordinatesFromFrameImage  = h_InputNeuronCoordinatesFromFrameImage;
HANDLES.h_Correct_FirstFrame                    = h_Correct_FirstFrame;
HANDLES.h_Correct_LastFrame                     = h_Correct_LastFrame;
HANDLES.h_RecalculateNeuronPositionWithinSegment= h_RecalculateNeuronPositionWithinSegment;
HANDLES.h_DisregardNeuronPositionWithinSegment  = h_DisregardNeuronPositionWithinSegment;
HANDLES.h_PlayMovieForCorrection                = h_PlayMovieForCorrection;
HANDLES.h_CorrectionFrameInterval_text          = h_CorrectionFrameInterval_text;
HANDLES.h_MaxPixelsPerFrame                     = h_MaxPixelsPerFrame;
HANDLES.h_PlayCurrentMovieInterval              = h_PlayCurrentMovieInterval;

data.HANDLES = HANDLES;
guidata(Fig_handle,data); % assign structure to the GUI


return

function [SegmentsToCorrect, FramesToDefineAsNans] = FindPossibleErrorSegmentation(Fig_handle, File,Tracks,HeadIsSqueezedMAT)

FrameRate               = File.FrameRate;
IntervalLengthThreshold = FrameRate;
NumOfTracks             = length(Tracks);

FramesForManualCorrection = HeadIsSqueezedMAT;
FramesToDefineAsNans      = HeadIsSqueezedMAT;

for tr=1:NumOfTracks
    FramesWhereWormIsOOB                               = Tracks(tr).Frames(Tracks(tr).OutOfBounds);
    FramesForManualCorrection(tr,FramesWhereWormIsOOB) = false;    
end
FramesToDefineAsNans(FramesForManualCorrection) = false;

% sum(FramesForManualCorrection,2)
% sum(FramesToDefineAsNans,2)
% figure; plot(FramesForManualCorrection(1,:),'b.'); hold on; plot(FramesToDefineAsNans(1,:),'ro'); ylim([0.9 1.1])

%% 
% clear SegmentsToCorrect  SegmentsToExclude
for tr=1:NumOfTracks
    CurrentFrames                      = find(FramesForManualCorrection(tr,:));
    if ~isempty(CurrentFrames)
        IndicesWherePossibleErrorInitiated = [1, find(diff(CurrentFrames)>=IntervalLengthThreshold)+1];
        IndicesWherePossibleErrorEnded     = [find(diff(CurrentFrames)>=IntervalLengthThreshold) length(CurrentFrames)];    
        FramesWherePossibleErrorInitiated  = CurrentFrames(IndicesWherePossibleErrorInitiated);
        FramesWherePossibleErrorEnded      = CurrentFrames(IndicesWherePossibleErrorEnded);    
    %     NumOfFramesInEachSegment           = IndicesWherePossibleErrorEnded-IndicesWherePossibleErrorInitiated+1;    
        NumOfFramesInEachSegment           = FramesWherePossibleErrorEnded-FramesWherePossibleErrorInitiated+1;    
        IndicesOfSegmentsToCorrect         = NumOfFramesInEachSegment>=FrameRate/3;   

        % SegmentsToCorrect = frame intervals where head is squeezed, worm is not OOB, and segment length is significantly long  
        SegmentsToCorrect(tr).InitialFrame         = FramesWherePossibleErrorInitiated(IndicesOfSegmentsToCorrect);
        SegmentsToCorrect(tr).EndFrame             = FramesWherePossibleErrorEnded(IndicesOfSegmentsToCorrect);
        SegmentsToCorrect(tr).NumOfFramesInEachSegment = NumOfFramesInEachSegment(IndicesOfSegmentsToCorrect);

        % SegmentsToExclude = frame intervals where head is squeezed, worm is not OOB, but segment length is short and not worth manual inspection  
        SegmentsToExclude(tr).InitialFrame         = FramesWherePossibleErrorInitiated(~IndicesOfSegmentsToCorrect);
        SegmentsToExclude(tr).EndFrame             = FramesWherePossibleErrorEnded(~IndicesOfSegmentsToCorrect);
        SegmentsToExclude(tr).NumOfFramesInEachSegment = NumOfFramesInEachSegment(~IndicesOfSegmentsToCorrect);

        for correction_ind = 1: length(SegmentsToExclude(tr).InitialFrame)
            FramesToDefineAsNans(tr,(SegmentsToExclude(tr).InitialFrame(correction_ind)):(SegmentsToExclude(tr).EndFrame(correction_ind)))      = true;
            FramesForManualCorrection(tr,(SegmentsToExclude(tr).InitialFrame(correction_ind)):(SegmentsToExclude(tr).EndFrame(correction_ind))) = false;
        end        
    else
        SegmentsToCorrect(tr).InitialFrame         = [];
        SegmentsToCorrect(tr).EndFrame             = [];
        SegmentsToCorrect(tr).NumOfFramesInEachSegment = [];        
    end        
    SegmentsAlreadyCorrected(tr).InitialFrame         = [];
    SegmentsAlreadyCorrected(tr).EndFrame             = [];
    SegmentsAlreadyCorrected(tr).NumOfFramesInEachSegment = [];        
end

% figure; plot(FramesForManualCorrection(3,:),'b.');  hold on;
% plot(FramesWherePossibleErrorInitiated,ones(size(FramesWherePossibleErrorInitiated)),'ro'); ylim([0.9 1.1])
% plot(FramesWherePossibleErrorEnded,ones(size(FramesWherePossibleErrorInitiated)),'go'); ylim([0.9 1.1])
% 
% sum(FramesForManualCorrection,2)
% sum(FramesToDefineAsNans,2)
% figure; plot(FramesForManualCorrection(1,:),'b.'); hold on; plot(FramesToDefineAsNans(1,:),'ro'); ylim([0.9 1.1])
% figure; plot(FramesForManualCorrection(2,:),'b.'); hold on; plot(FramesToDefineAsNans(2,:),'ro'); ylim([0.9 1.1])
% figure; plot(FramesForManualCorrection(3,:),'b.'); hold on; plot(FramesToDefineAsNans(3,:),'ro'); ylim([0.9 1.1])

%% Assign to data structure:
data                    = guidata(Fig_handle); 
for tr=1:NumOfTracks
    data.NeuronDisplayMatrices_ManuallyCorrected(tr,FramesToDefineAsNans(tr,:),:,:) = NaN;
    data.NeuronDisplayMatrices_Corrected(tr,FramesToDefineAsNans(tr,:),:,:)         = NaN;
    data.NeuronDisplayMatrices(tr,FramesToDefineAsNans(tr,:),:,:)                   = NaN;
    data.NeuronCoordinatesMatrix(tr,FramesToDefineAsNans(tr,:),:)                   = NaN;
    data.NeuronCoordinatesMatrix_Corrected(tr,FramesToDefineAsNans(tr,:),:)         = NaN;
    data.NeuronCoordinatesMatrix_ManuallyCorrected(tr,FramesToDefineAsNans(tr,:),:) = NaN;
    data.ManuallyAssignedAsDisregard(tr,FramesToDefineAsNans(tr,:),:)               = true;
end
data.SegmentsToCorrect         = SegmentsToCorrect;
data.SegmentsAlreadyCorrected  = SegmentsAlreadyCorrected;
data.FramesForManualCorrection = FramesForManualCorrection;
guidata(Fig_handle, data);

return

function Plot_WormPixels_on_Frame(CurrentTrack, Frame_Num, PixelsField, COLOR, ZoomOnWorm )

FrameIndex = find( CurrentTrack.Frames == Frame_Num, 1, 'first'); 
if ZoomOnWorm 
    MARKERSIZE  = 6;        
else   
    MARKERSIZE = 2.5;
end

if ~exist('COLOR','var')
    COLOR   = 'r';
end
LINESTYLE   = 'none'; % for PixelIdxList and HeadTail_PixelIdxList
MARKER      = '.';    % for PixelIdxList

if strcmpi(PixelsField,'Midline') 
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        hold on;
        plot(CurrentTrack.Midline.Y_coordinates_short{FrameIndex},CurrentTrack.Midline.X_coordinates_short{FrameIndex},'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);   
    end
elseif strcmpi(PixelsField,'Perimeter')
    MARKERSIZE = 4;
%     LINESTYLE   = '-';    
%     MARKER      = 'none'; 
    hold on;
    plot(single(CurrentTrack.WormPerimeter.Ycoordinate{FrameIndex}),single(CurrentTrack.WormPerimeter.Xcoordinate{FrameIndex}),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);    
elseif strcmpi(PixelsField,'HeadTail')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = '*';  
        if ZoomOnWorm 
            MARKERSIZE  = 15;
        else            
            MARKERSIZE  = 8; 
        end
        hold on;
        HeadTail_Matrix = squeeze(CurrentTrack.HeadTail(FrameIndex,:,:));
        plot(HeadTail_Matrix(:,2),HeadTail_Matrix(:,1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);           
    end
elseif strcmpi(PixelsField,'Head')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = '*';  
        if ZoomOnWorm 
            MARKERSIZE  = 15;
        else
            MARKERSIZE  = 10;
        end
        hold on;
        Head_Vec = CurrentTrack.Head(FrameIndex,:);
        plot(Head_Vec(2),Head_Vec(1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);         
    end
elseif strcmpi(PixelsField,'Tail')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = 'o';  
        if ZoomOnWorm 
            MARKERSIZE  = 15;
        else   
            MARKERSIZE  = 10; 
        end
        hold on;
        Tail_Vec = CurrentTrack.Tail(FrameIndex,:);
        plot(Tail_Vec(2),Tail_Vec(1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);     
    end
elseif strcmpi(PixelsField,'NeuronCoordinates_Estimations')
    if CurrentTrack.Midline.FlagTrueIfReliable(FrameIndex)
        MARKER      = 'o';  
        if ZoomOnWorm 
            MARKERSIZE  = 15;
        else
            MARKERSIZE  = 10;
        end
        hold on;
        NeuronCoordinates_Estimations = CurrentTrack.NeuronCoordinates_Estimations(FrameIndex,:);
        plot(NeuronCoordinates_Estimations(2),NeuronCoordinates_Estimations(1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);         
    end
elseif strcmpi(PixelsField,'NeuronCoordinates')
    
    MARKER      = 's';  
    if ZoomOnWorm 
        MARKERSIZE  = 15;
    else   
        MARKERSIZE  = 10; 
    end
    hold on;
    NeuronCoordinates_Vec = CurrentTrack.NeuronCoordinates(FrameIndex,:);
    plot(NeuronCoordinates_Vec(2),NeuronCoordinates_Vec(1),'color',COLOR,'marker',MARKER,'linestyle',LINESTYLE,'markersize',MARKERSIZE);     
          
end

return

function DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions)
COLORS_FOR_TRACKS   = {'w','g','c','m','b','y',[0.8 0.5 0.2],[0.8 0.5 1]}; % it will start from the green

data                    = guidata(Fig_handle); 
HANDLES                 = data.HANDLES;
ALLFRAMES               = data.ALLFRAMES;
NeuronCoordinatesMatrix = data.NeuronCoordinatesMatrix_Corrected;
NeuronDisplayMatrices   = data.NeuronDisplayMatrices_Corrected;
ZoomOnWorm              = HANDLES.h_ZoomOnWorm_ONOFF.Value==1;
ZoomOnWorm_tracknum     = HANDLES.h_ZoomOnWorm_tracknum.Value;
MovieOptions.BestScale(1) = str2num(HANDLES.h_ScaleRange_Min.String);
MovieOptions.BestScale(2) = str2num(HANDLES.h_ScaleRange_Max.String);

figure(Fig_handle);        
if isfield(HANDLES,'FrameImage')
    FrameImageHandle  = HANDLES.FrameImageHandle;
    ZoomImagesHandles = HANDLES.ZoomImagesHandles;
else
    subplot(1,2,1);
    position = get(gca,'position'); position(1) = position(1)-0.08;
    set(gca,'position',position)
    FrameImageHandle = gca;
    for t_ind = 1:length(Tracks)
        subplot(length(Tracks),3,3*t_ind);
        position = get(gca,'position'); position(1) = position(1)-0.32;
        set(gca,'position',position)
        ZoomImagesHandles(t_ind)= gca;
    end
    HANDLES.FrameImageHandle  = FrameImageHandle;
    HANDLES.ZoomImagesHandles = ZoomImagesHandles;
    data.HANDLES = HANDLES;
    guidata(Fig_handle, data);     
end
   
ShowTracksText   = MovieOptions.ShowTracksText;
ColorByID        = MovieOptions.TracksText.ColorByID;
arena            = MovieOptions.arena;   
TracksLine       = MovieOptions.TracksLine;
ShowWormFeatures = MovieOptions.ShowWormFeatures;
Add_CTXdata_And_SmoothPosition = MovieOptions.Add_CTXdata_And_SmoothPosition;    

% Get Frame:    assume Video Format == 1 !!!!!!!!!!!   -->    A sequence of tiff files  
FileFullName  = [File.MoviePath,'\',File.MovieFileNames{Frame}]; 
Mov           = imread(FileFullName);                            
axes(FrameImageHandle);
imshow(Mov,MovieOptions.BestScale,'initialmagnification',600);    

if isfield(MovieOptions,'All_Mask_Countours')  % Display Mask coundaries
    hold on; plot(MovieOptions.All_Mask_Countours(:,2),...
                  MovieOptions.All_Mask_Countours(:,1), '-', 'color', [0.25 0.25 0.25], 'linewidth',2)
end

if ~ Add_CTXdata_And_SmoothPosition  
    if ~isnan(arena)        
        axis(MovieOptions.ArenaBox); 
        title(['Arena ',num2str(arena),', Frame ',num2str(Frame)]);
    else
        title(['Frame ',num2str(Frame)]);
    end
else
    if ~isnan(arena)        
        axis(MovieOptions.ArenaBox); 
    end
end

if ZoomOnWorm         
    FrameIndexInTrack = find([Tracks(ZoomOnWorm_tracknum).Frames]==Frame,1);
    if ~isempty(FrameIndexInTrack)   
        HalfWormMaxSize  = File.VariablesInformation.MaxWorm_size/0.05 * File.PixelSize / 2;      % File.VariablesInformation.MaxWorm_size ~ 1mm*0.05 mm
        WormCentroid = Tracks(ZoomOnWorm_tracknum).Path(FrameIndexInTrack,:);
        Yaxis(1)   = max(1,           WormCentroid(1)-HalfWormMaxSize*0.75); 
        Yaxis(2)   = min(size(Mov,1),WormCentroid(1)+HalfWormMaxSize*0.75); 
        Xaxis(1)   = max(1,           WormCentroid(2)-HalfWormMaxSize*0.75); 
        Xaxis(2)   = min(size(Mov,2),WormCentroid(2)+HalfWormMaxSize*0.75); 
        Xaxis = round(Xaxis);
        Yaxis = round(Yaxis);
        axis([Yaxis Xaxis]);
    end
end
    
if Add_CTXdata_And_SmoothPosition            
    RelevantTracks = 1:length(Tracks);            
else        
    RelevantTracks = find(ALLFRAMES(Frame,:));
end
        
hold on;
All_CTX_position = [];
for t_ind = 1:length(RelevantTracks);        
    axes(FrameImageHandle)
    t_num       = RelevantTracks(t_ind);
    IndexInTracks = find(Tracks(t_num).Frames==Frame,1);
    if isempty(IndexInTracks)
        Coordinates = [];
    else
        Coordinates = double(Tracks(t_num).Path(IndexInTracks,:));        
        Coordinates = Coordinates([2 1]);   % Flip X-Y for imshow
    end
        
    % Add Track line of previous frames
    if TracksLine
        RelevantPathFrames = find(Tracks(t_num).Frames<=Frame,1e3,'last');                      % Show all previous track path            
        if ~islogical(TracksLine) && (TracksLine < length(RelevantPathFrames)) % Show only the specified number of previous frames)
            RelevantPathFrames = RelevantPathFrames((end-TracksLine+1):end);       
        end
        Path = Tracks(t_num).Path(RelevantPathFrames,:);           
        line(Path(:,2), Path(:,1),'linewidth',1,'color','k','linestyle','-');
    end
    % ColumnIndex in Data_Matrices
    if Add_CTXdata_And_SmoothPosition  || MovieOptions.ShowBehaviorFrom_Data_Matrices 
         ColumnIndex  = find(BehaviorDataMatrices.(MovieOptions.FramesField)(t_num,:) == Frame,1);   
        if isempty(ColumnIndex) 
         ColumnIndex  = find(BehaviorDataMatrices.(MovieOptions.FramesField)(t_num,:) == (Frame-1),1);   
        end     
    end

    % Defining display text for each track
    if ShowTracksText && ~isempty(Coordinates)
        text_for_Movie = '';
        if MovieOptions.TracksText.TrackNum
            text_for_Movie = [text_for_Movie,'(#',num2str(t_num),')'];
        end        
        if MovieOptions.TracksText.Behavior
            if MovieOptions.ShowBehaviorFrom_Data_Matrices &&  ~isempty(ColumnIndex)
                CurrentBehCode = BehaviorDataMatrices.(MovieOptions.BehaviorFieldName)(t_num,ColumnIndex);
            else
                CurrentBehCode = Tracks(t_num).(MovieOptions.BehaviorFieldName)(IndexInTracks);
            end
%                 CurrentBehColor=  BehaviorStructure.BehaviorCodeColors{BehaviorStructure.BehaviorCodeNumbers==CurrentBehCode};
            CurrentBehName =  BehaviorStructure.BehaviorCodeName{BehaviorStructure.BehaviorCodeNumbers==CurrentBehCode};
            if isempty(text_for_Movie)
                text_for_Movie = CurrentBehName;
            else                    
                text_for_Movie = [text_for_Movie,char(10),CurrentBehName];
            end

        end 
    else
%             text_for_Movie = num2str(t_num);
        text_for_Movie = '';
    end
        
    if Add_CTXdata_And_SmoothPosition  &&  ~isempty(ColumnIndex)
        Current_Coordinates_X_Smoothed = double(BehaviorDataMatrices.(MovieOptions.X_field)(t_num,ColumnIndex));
        Current_Coordinates_Y_Smoothed = double(BehaviorDataMatrices.(MovieOptions.Y_field)(t_num,ColumnIndex));
        Current_CTX_position           = BehaviorDataMatrices.(MovieOptions.CTX_field)(t_num,ColumnIndex);
        All_CTX_position = [All_CTX_position, Current_CTX_position];

        CTXtext  = num2str(Current_CTX_position,2);
    %   From G --> R through brown
    %   R        = abs((Current_CTX_position+1)/2); 
    %   G        = abs((Current_CTX_position-1)/2); 
    %   B        = 0; 
    %   From G --> R through black
            
        if Current_CTX_position>0
            R = Current_CTX_position;
            G = 0;
        else
            G = -Current_CTX_position;
            R = 0;
        end
        B        = 0; 
        CTXcolor = [R G B];            
    else
        CTXtext  = '';
        CTXcolor = 'b';                               
    end   

    CurrentColor = 'b';
    if ~isempty(Coordinates)
        if ColorByID ==2      % Color set by track number
            CurrentColor = COLORS_FOR_TRACKS{MovieOptions.TracksColorsIndices(t_num)};
        elseif ColorByID ==1  % Color set by worm ID
            Worm_ID = Tracks(t_num).ID;
            CurrentColor = COLORS_FOR_TRACKS{MovieOptions.TracksColorsIndices(Worm_ID)};
        elseif ColorByID ==3  % Color set by Chemotaxis index
            CurrentColor = CTXcolor;
        elseif ColorByID ==4  % Color set by Behavior code
            if MovieOptions.ShowBehaviorFrom_Data_Matrices &&  ~isempty(ColumnIndex)
                CurrentBehCode =  BehaviorDataMatrices.(MovieOptions.BehaviorFieldName)(t_num,ColumnIndex);
            else
                CurrentBehCode =  Tracks(t_num).(MovieOptions.BehaviorFieldName)(IndexInTracks);
            end
            CurrentBehColor=  BehaviorStructure.BehaviorCodeColors{BehaviorStructure.BehaviorCodeNumbers==CurrentBehCode};
            CurrentColor   = CurrentBehColor;
            CTXcolor       = CurrentColor;
        end
    else
        if ColorByID ==4
            CTXcolor = 'y';
        end                
    end  
        
    if Add_CTXdata_And_SmoothPosition && ~isempty(CTXtext)
        text(Current_Coordinates_Y_Smoothed-30,Current_Coordinates_X_Smoothed-30, CTXtext,'FontWeight','bold','color',CTXcolor);             % Track count
        plot(Current_Coordinates_Y_Smoothed   ,Current_Coordinates_X_Smoothed,'o','markersize',6,'color',CTXcolor,'markerfacecolor',CTXcolor); 
    end
    if ~isempty(Coordinates)
        text(Coordinates(2)+10,Coordinates(1),text_for_Movie,'FontWeight','bold','color',CurrentColor);             % Track count
        plot(Coordinates(2)   ,Coordinates(1),'.','markersize',12,'color',CurrentColor); 
    end
        
    %% Plot detected pixels
    if  ShowWormFeatures  && ~isempty(Coordinates)
        for f_num = 1:length(MovieOptions.WormFeatures)
            CurrentField = MovieOptions.WormFeatures{f_num};        % name of field to display
            if ~strcmpi(CurrentField,'NeuronCoordinates')
                CurrentColor = MovieOptions.WormFeaturesCOLORS{f_num};  % color
                if ColorByID ==2      % Color set by track number
                    CurrentColor = COLORS_FOR_TRACKS{MovieOptions.TracksColorsIndices(t_num)};
                elseif ColorByID ==1  % Color set by worm ID
                    Worm_ID = Tracks(t_num).ID;
                    CurrentColor = COLORS_FOR_TRACKS{MovieOptions.TracksColorsIndices(Worm_ID)};
                elseif ColorByID ==3  % Color set by Chemotaxis index
                    CurrentColor = CTXcolor;                    
                elseif ColorByID ==4  % Color set by Behavior code
                    CurrentColor = CurrentBehColor;                    
                end                         
                Plot_WormPixels_on_Frame(Tracks(t_num), Frame, CurrentField, CurrentColor, ZoomOnWorm)              
            end
        end
    end  
end

%% Display neuron matrices
if ZoomOnWorm, MARKERSIZE  = 15; else  MARKERSIZE  = 10; end
for t_ind = 1:length(Tracks)    
    CurrentNeuronCoordinates = squeeze(NeuronCoordinatesMatrix(t_ind,Frame,:));
    CurrentMAT               = squeeze(NeuronDisplayMatrices(t_ind,Frame,:,:));            
    if (ColorByID ==2)||(ColorByID ==1)      % Color set by track number
        CurrentColor = COLORS_FOR_TRACKS{MovieOptions.TracksColorsIndices(t_ind)};                
    else
        CurrentColor = 'b';
    end
    axes(ZoomImagesHandles(t_ind)); hold off;
    imshow(CurrentMAT,MovieOptions.MinMaxNeuronDisplay,'initialmagnification',1000);
    title(['track ',num2str(t_ind)],'color',CurrentColor);

    if ~isnan(CurrentNeuronCoordinates(1)) 
       hold on; rectangle('position',[10 10 5 5],'edgecolor',CurrentColor); hold off; 
       if ~isempty(find(strcmpi(MovieOptions.WormFeatures,'NeuronCoordinates'),1))
            axes(FrameImageHandle); hold on;
            plot(CurrentNeuronCoordinates(2),CurrentNeuronCoordinates(1),'color',CurrentColor,'marker','s','linestyle','none','markersize',MARKERSIZE);     
        end
    end
end  
%%    
if Add_CTXdata_And_SmoothPosition  
    if ~isnan(arena)                       
        title(['Arena ',num2str(arena),', Frame ',num2str(Frame),', CTX ',num2str(mean(All_CTX_position),2)]);
    else
        title(['Frame ',num2str(Frame),', CTX ',num2str(mean(All_CTX_position),2)]);
    end
end   

axes(FrameImageHandle); 
hold off;


return

function CreateMovie(stam, stam2, Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions)
data                    = guidata(Fig_handle); 
HANDLES                 = data.HANDLES;
StartFrame              = str2num(HANDLES.h_CurrentFrame.String);
LastFrame               = str2num(HANDLES.h_MovieLastFrame.String);
PauseTime               = str2num(HANDLES.h_TimeBwFrames.String);
NumOfFrames             = File.NumberOfFrames;
LastFrame               = min([LastFrame NumOfFrames]); 
MovieName               = HANDLES.h_MovieName.String; 
MovieDisplayRate        = str2num(HANDLES.h_MovieDisplayRate.String);

writerObj           = VideoWriter([MovieName,'.avi'],'Uncompressed AVI'); 
writerObj.FrameRate = MovieDisplayRate;
open(writerObj);

Frame = StartFrame;
while Frame <= LastFrame
    DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions);
    pause(PauseTime);
    if HANDLES.h_StopMovie.Value==0
        break
    end
    FrameInterval = str2num(HANDLES.h_FrameInterval.String);
    Frame         = str2num(HANDLES.h_CurrentFrame.String)+FrameInterval;
    set(HANDLES.h_CurrentFrame,'Value',Frame,'String',num2str(Frame));
    PauseTime     = str2num(HANDLES.h_TimeBwFrames.String);
    set(HANDLES.h_TimeBwFrames,'Value',PauseTime,'String',num2str(PauseTime));
    
    movieframe = getframe(Fig_handle);
    writeVideo(writerObj,movieframe);
end

set(HANDLES.h_StopMovie,'Value',2);

close(writerObj); 

return

function PlayMovie(stam, stam2, Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions)
data                    = guidata(Fig_handle); 
HANDLES                 = data.HANDLES;
StartFrame              = str2num(HANDLES.h_CurrentFrame.String);
PauseTime               = str2num(HANDLES.h_TimeBwFrames.String);
NumOfFrames             = File.NumberOfFrames;

Frame = StartFrame;
while Frame <= NumOfFrames
    DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions);
    pause(PauseTime);
    if HANDLES.h_StopMovie.Value==0
        break
    end
    FrameInterval = str2num(HANDLES.h_FrameInterval.String);
    Frame         = str2num(HANDLES.h_CurrentFrame.String)+FrameInterval;
    set(HANDLES.h_CurrentFrame,'Value',Frame,'String',num2str(Frame));
    PauseTime     = str2num(HANDLES.h_TimeBwFrames.String);
    set(HANDLES.h_TimeBwFrames,'Value',PauseTime,'String',num2str(PauseTime));
end

set(HANDLES.h_StopMovie,'Value',2);

return

function PlayMovieBackwards(stam, stam2, Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions)
data                    = guidata(Fig_handle); 
HANDLES                 = data.HANDLES;
StartFrame              = str2num(HANDLES.h_CurrentFrame.String);
PauseTime               = str2num(HANDLES.h_TimeBwFrames.String);
FrameInterval           = str2num(HANDLES.h_FrameInterval.String);

Frame = StartFrame;
while Frame >= FrameInterval+1
    DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions);
    pause(PauseTime);
    if HANDLES.h_StopMovie.Value==0
        break
    end
    FrameInterval = str2num(HANDLES.h_FrameInterval.String);
    Frame         = str2num(HANDLES.h_CurrentFrame.String)-FrameInterval;
    set(HANDLES.h_CurrentFrame,'Value',Frame,'String',num2str(Frame));
    PauseTime     = str2num(HANDLES.h_TimeBwFrames.String);
    set(HANDLES.h_TimeBwFrames,'Value',PauseTime,'String',num2str(PauseTime));
end

set(HANDLES.h_StopMovie,'Value',2);

return

function ShowSingleFrame(stam, stam2, Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions)
data     = guidata(Fig_handle); 
HANDLES  = data.HANDLES;
Frame    = str2num(HANDLES.h_CurrentFrame.String);  
DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions);
return

function ShowNextFrame(stam, stam2, Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions)
data     = guidata(Fig_handle); 
HANDLES  = data.HANDLES;
Frame    = min([str2num(HANDLES.h_CurrentFrame.String)+1 File.NumberOfFrames]);  
set(HANDLES.h_CurrentFrame,'Value',Frame,'String',num2str(Frame));
DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions);
return

function ShowPreviousFrame(stam, stam2, Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions)
data     = guidata(Fig_handle); 
HANDLES  = data.HANDLES;
Frame    = max([str2num(HANDLES.h_CurrentFrame.String)-1  1]);  
set(HANDLES.h_CurrentFrame,'Value',Frame,'String',num2str(Frame));
DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions);
return

function DefineFrameAsCorrect(stam, stam2, Fig_handle, tr)
data      = guidata(Fig_handle); 
HANDLES   = data.HANDLES;
Frame     = str2num(HANDLES.h_CurrentFrame.String);

data.FramesForNeuronPositionReCalculation(tr,Frame) = true;
guidata(Fig_handle, data);

CorrectionFirstFrame = str2num(HANDLES.h_Correct_FirstFrame(tr).String);
if (CorrectionFirstFrame>Frame)||(CorrectionFirstFrame==0)
    set(HANDLES.h_Correct_FirstFrame(tr),'Value',Frame,'String',num2str(Frame));
end
CorrectionLastFrame = str2num(HANDLES.h_Correct_LastFrame(tr).String);
if (CorrectionLastFrame<Frame)||(CorrectionLastFrame==0)
    set(HANDLES.h_Correct_LastFrame(tr),'Value',Frame,'String',num2str(Frame));
end

return

function DefineFrameAsDisregard(stam, stam2, Fig_handle, tr, File, Tracks, BehaviorDataMatrices, MovieOptions)
data      = guidata(Fig_handle); 
HANDLES   = data.HANDLES;
Frame     = str2num(HANDLES.h_CurrentFrame.String);

% NewCoordinates_IN_FULL_FRAME                                = [NaN NaN];
data.NeuronCoordinatesMatrix_ManuallyCorrected(tr,Frame,:)  = NaN;
data.NeuronCoordinatesMatrix_Corrected(tr,Frame,:)          = NaN;
data.ManuallyAssignedAsDisregard(tr,Frame)                  = true;
data.FramesForNeuronPositionReCalculation(tr,Frame)         = false;

data.NeuronDisplayMatrices_Corrected(tr,Frame,:,:)         = NaN;
data.NeuronDisplayMatrices_ManuallyCorrected(tr,Frame,:,:) = NaN;
guidata(Fig_handle, data);

DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions)

return

function InputNeuronCoordinatesFromZoomedImage(stam, stam2, Fig_handle, tr, File, Tracks, BehaviorDataMatrices, MovieOptions, background, VignettingPattern, WormArea)
data      = guidata(Fig_handle); 
HANDLES   = data.HANDLES;
Frame     = str2num(HANDLES.h_CurrentFrame.String);

ZoomImageHandle = HANDLES.ZoomImagesHandles(tr);
axes(ZoomImageHandle);
DisplaySquareSize = size(data.NeuronDisplayMatrices_Corrected,3);
[x,y] = ginput(1);

CurrentCoordinates_IN_FULL_FRAME = squeeze(data.NeuronCoordinatesMatrix_Corrected(tr,Frame,:))';
NewCoordinates_IN_FULL_FRAME     = round(CurrentCoordinates_IN_FULL_FRAME+[y,x]-(DisplaySquareSize/2)*ones(1,2));

data.NeuronCoordinatesMatrix_ManuallyCorrected(tr,Frame,:)  = NewCoordinates_IN_FULL_FRAME;
data.NeuronCoordinatesMatrix_Corrected(tr,Frame,:)          = NewCoordinates_IN_FULL_FRAME;
data.ManuallyAssignedAsCorrect(tr,Frame)                    = true;
data.FramesForNeuronPositionReCalculation(tr,Frame)         = true;

CorrectionFirstFrame = str2num(HANDLES.h_Correct_FirstFrame(tr).String);
if (CorrectionFirstFrame>Frame)||(CorrectionFirstFrame==0)
    set(HANDLES.h_Correct_FirstFrame(tr),'Value',Frame,'String',num2str(Frame));
end
CorrectionLastFrame = str2num(HANDLES.h_Correct_LastFrame(tr).String);
if (CorrectionLastFrame<Frame)||(CorrectionLastFrame==0)
    set(HANDLES.h_Correct_LastFrame(tr),'Value',Frame,'String',num2str(Frame));
end

IndexInTracks = data.ConvertMovieFrameToIndexWithinTracks(tr,Frame);
NormFrame     = ReadAndNormalizeFrame(tr,Frame, File, background, VignettingPattern, WormArea, IndexInTracks);
DisplayMatrix = Recalculate_DisplayMatrix(NormFrame, NewCoordinates_IN_FULL_FRAME, DisplaySquareSize);

data.NeuronDisplayMatrices_Corrected(tr,Frame,:,:)         = DisplayMatrix;
data.NeuronDisplayMatrices_ManuallyCorrected(tr,Frame,:,:) = DisplayMatrix;
guidata(Fig_handle, data);

DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions)

return

function InputNeuronCoordinatesFromFrameImage(stam, stam2, Fig_handle, tr, File, Tracks, BehaviorDataMatrices, MovieOptions, background, VignettingPattern, WormArea)
data      = guidata(Fig_handle); 
HANDLES   = data.HANDLES;
Frame     = str2num(HANDLES.h_CurrentFrame.String);

FrameImageHandle = HANDLES.FrameImageHandle;
axes(FrameImageHandle);
DisplaySquareSize = size(data.NeuronDisplayMatrices_Corrected,3);
[x,y] = ginput(1);

NewCoordinates_IN_FULL_FRAME     = round([y,x]);

data.NeuronCoordinatesMatrix_ManuallyCorrected(tr,Frame,:)  = NewCoordinates_IN_FULL_FRAME;
data.NeuronCoordinatesMatrix_Corrected(tr,Frame,:)          = NewCoordinates_IN_FULL_FRAME;
data.ManuallyAssignedAsCorrect(tr,Frame)                    = true;
data.FramesForNeuronPositionReCalculation(tr,Frame)         = true;

CorrectionFirstFrame = str2num(HANDLES.h_Correct_FirstFrame(tr).String);
if (CorrectionFirstFrame>Frame)||(CorrectionFirstFrame==0)
    set(HANDLES.h_Correct_FirstFrame(tr),'Value',Frame,'String',num2str(Frame));
end
CorrectionLastFrame = str2num(HANDLES.h_Correct_LastFrame(tr).String);
if (CorrectionLastFrame<Frame)||(CorrectionLastFrame==0)
    set(HANDLES.h_Correct_LastFrame(tr),'Value',Frame,'String',num2str(Frame));
end

IndexInTracks = data.ConvertMovieFrameToIndexWithinTracks(tr,Frame);

NormFrame     = ReadAndNormalizeFrame(tr,Frame, File, background, VignettingPattern, WormArea, IndexInTracks);
DisplayMatrix = Recalculate_DisplayMatrix(NormFrame, NewCoordinates_IN_FULL_FRAME, DisplaySquareSize);

data.NeuronDisplayMatrices_Corrected(tr,Frame,:,:)         = DisplayMatrix;
data.NeuronDisplayMatrices_ManuallyCorrected(tr,Frame,:,:) = DisplayMatrix;
guidata(Fig_handle, data);

DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions)

return

function NormFrame = ReadAndNormalizeFrame(tr,Frame, File, background, VignettingPattern, WormArea, IndexInTracks)

Mov = imread([File.MoviePath,'\',File.MovieFileNames{Frame}]);                          
Mov = imsubtract(Mov,background);
%% Keep Only A single Worm Object and correct for vignetting patterns
NormFrame = zeros(File.FrameSize,'single')*NaN;
if ~isnan(IndexInTracks)
    NormFrame(WormArea(tr).LinearIndices{IndexInTracks}) = Mov(WormArea(tr).LinearIndices{IndexInTracks});
end  
NormFrame = NormFrame./VignettingPattern;
         
return

function DisplayMatrix = Recalculate_DisplayMatrix(NormFrame, NeuronCoordinates, DisplaySquareSize)

DisplayMatrix           = zeros(DisplaySquareSize,DisplaySquareSize,'single')*NaN;
MaximalDistanceInPixels = ones(2)*(DisplaySquareSize-1)/2;
MidPixelInDisplay       = (DisplaySquareSize+1)/2;

while ((NeuronCoordinates(1)-MaximalDistanceInPixels(1,1))<1)
    MaximalDistanceInPixels(1,1) = MaximalDistanceInPixels(1,1)-1;
end
while ((NeuronCoordinates(2)-MaximalDistanceInPixels(2,1))<1)
    MaximalDistanceInPixels(2,1) = MaximalDistanceInPixels(2,1)-1;
end
while ((NeuronCoordinates(1)+MaximalDistanceInPixels(1,2))> size(NormFrame,1))
    MaximalDistanceInPixels(1,2) = MaximalDistanceInPixels(1,2)-1;
end
while ((NeuronCoordinates(2)+MaximalDistanceInPixels(2,2))> size(NormFrame,2))
    MaximalDistanceInPixels(2,2) = MaximalDistanceInPixels(2,2)-1;
end

MAT = NormFrame((NeuronCoordinates(1)-MaximalDistanceInPixels(1,1)):(NeuronCoordinates(1)+MaximalDistanceInPixels(1,2)), ...
                (NeuronCoordinates(2)-MaximalDistanceInPixels(2,1)):(NeuronCoordinates(2)+MaximalDistanceInPixels(2,2)));

DisplayMatrix((MidPixelInDisplay-MaximalDistanceInPixels(1,1)):(MidPixelInDisplay+MaximalDistanceInPixels(1,2)), ...
              (MidPixelInDisplay-MaximalDistanceInPixels(2,1)):(MidPixelInDisplay+MaximalDistanceInPixels(2,2))) = MAT;
         
return

function PlayMovieForCorrection(stam, stam2, Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions, tr)
data                    = guidata(Fig_handle); 
if isempty(data.SegmentsToCorrect(tr).InitialFrame) % no flagged frames for correction
    return
end

HANDLES    = data.HANDLES;
StartFrame = data.SegmentsToCorrect(tr).InitialFrame(1);
EndFrame   = min([data.SegmentsToCorrect(tr).EndFrame(1) File.NumberOfFrames]);
Frame      = StartFrame;
PauseTime  = str2num(HANDLES.h_TimeBwFrames.String);

set(HANDLES.h_CurrentFrame,'Value',Frame,'String',num2str(Frame));

while Frame <= EndFrame
    DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions);
    pause(PauseTime);
    if HANDLES.h_StopMovie.Value==0
        break
    end
    FrameInterval = str2num(HANDLES.h_FrameInterval.String);
    Frame         = str2num(HANDLES.h_CurrentFrame.String)+FrameInterval;
    set(HANDLES.h_CurrentFrame,'Value',Frame,'String',num2str(Frame));
    PauseTime     = str2num(HANDLES.h_TimeBwFrames.String);
    set(HANDLES.h_TimeBwFrames,'Value',PauseTime,'String',num2str(PauseTime));
end

set(HANDLES.h_StopMovie,'Value',2);

return

function PlayCurrentMovieInterval(stam, stam2, Fig_handle, File, Tracks, BehaviorDataMatrices, MovieOptions, tr)
data       = guidata(Fig_handle); 
HANDLES    = data.HANDLES;
StartFrame = str2num(HANDLES.h_Correct_FirstFrame(tr).String);
EndFrame   = str2num(HANDLES.h_Correct_LastFrame(tr).String);

Frame      = StartFrame;
PauseTime  = str2num(HANDLES.h_TimeBwFrames.String);

set(HANDLES.h_CurrentFrame,'Value',Frame,'String',num2str(Frame));

while Frame <= EndFrame
    DisplayFrame(Fig_handle, Frame, File, Tracks, BehaviorDataMatrices, MovieOptions);
    pause(PauseTime);
    if HANDLES.h_StopMovie.Value==0
        break
    end
    FrameInterval = str2num(HANDLES.h_FrameInterval.String);
    Frame         = str2num(HANDLES.h_CurrentFrame.String)+FrameInterval;
    set(HANDLES.h_CurrentFrame,'Value',Frame,'String',num2str(Frame));
    PauseTime     = str2num(HANDLES.h_TimeBwFrames.String);
    set(HANDLES.h_TimeBwFrames,'Value',PauseTime,'String',num2str(PauseTime));
end

set(HANDLES.h_StopMovie,'Value',2);

return

function ConfirmProperSegmentation(stam, stam2, Fig_handle, tr)
data                     = guidata(Fig_handle); 
HANDLES                  = data.HANDLES;
SegmentsToCorrect        = data.SegmentsToCorrect;
SegmentsAlreadyCorrected = data.SegmentsAlreadyCorrected;

SegmentsAlreadyCorrected(tr).InitialFrame             = [SegmentsAlreadyCorrected(tr).InitialFrame,           SegmentsToCorrect(tr).InitialFrame(1)];
SegmentsAlreadyCorrected(tr).EndFrame                 = [SegmentsAlreadyCorrected(tr).EndFrame, SegmentsToCorrect(tr).EndFrame(1)];
SegmentsAlreadyCorrected(tr).NumOfFramesInEachSegment = [SegmentsAlreadyCorrected(tr).NumOfFramesInEachSegment, SegmentsToCorrect(tr).NumOfFramesInEachSegment(1)];

StartFrame = data.SegmentsToCorrect(tr).InitialFrame(1);
EndFrame   = data.SegmentsToCorrect(tr).EndFrame(1);

data.ManuallyAssignedAsCorrect(tr,StartFrame:EndFrame) = true;
data.FramesForNeuronPositionReCalculation(tr,:)        = false;

NumberOfSegmentsLeft = length(SegmentsToCorrect(tr).InitialFrame)-1;    
if NumberOfSegmentsLeft>=1;
    SegmentsToCorrect(tr).InitialFrame             = SegmentsToCorrect(tr).InitialFrame(2:end);
    SegmentsToCorrect(tr).EndFrame                 = SegmentsToCorrect(tr).EndFrame(2:end);
    SegmentsToCorrect(tr).NumOfFramesInEachSegment = SegmentsToCorrect(tr).NumOfFramesInEachSegment(2:end);
    
    NextStartFrame    = SegmentsToCorrect(tr).InitialFrame(1);
    NextEndFrame      = SegmentsToCorrect(tr).EndFrame(1);
    FrameIntervalText = ['Relevant frame interval= [',num2str(NextStartFrame),' - ',num2str(NextEndFrame),']'];
else
    SegmentsToCorrect(tr).InitialFrame             = [];
    SegmentsToCorrect(tr).EndFrame                 = [];
    SegmentsToCorrect(tr).NumOfFramesInEachSegment = [];
    FrameIntervalText                              = 'No flagged frames';
end
set(HANDLES.h_CorrectionFrameInterval_text(tr),'String',FrameIntervalText);
set(HANDLES.h_Correct_FirstFrame(tr),'Value',0,'String','0');
set(HANDLES.h_Correct_LastFrame(tr),'Value',0,'String','0');
    
data.SegmentsToCorrect        = SegmentsToCorrect;
data.SegmentsAlreadyCorrected = SegmentsAlreadyCorrected;       
guidata(Fig_handle, data);

return

function RecalculateNeuronPositionWithinSegment(stam, stam2, Fig_handle, tr, File, background, VignettingPattern, WormArea)

MinimumPosition                                           = 3;
MaximumPosition_vec                                       = File.FrameSize-3;

data                 = guidata(Fig_handle); 
HANDLES              = data.HANDLES;
FirstCorrectionFrame = str2num(HANDLES.h_Correct_FirstFrame(tr).String);
LastCorrectionFrame  = str2num(HANDLES.h_Correct_LastFrame(tr).String);
MaxPixelsPerFrame    = str2num(HANDLES.h_MaxPixelsPerFrame(tr).String);  % 2 is good for "omega" behaviour  

FramesToAnalyze      = FirstCorrectionFrame:LastCorrectionFrame;
DisplaySquareSize    = size(data.NeuronDisplayMatrices_Corrected,3);

FramesFlaggedAsCorrectForCalculation = data.FramesForNeuronPositionReCalculation(tr,:);
% TrueIfFlaggedAsCorrect = data.ManuallyAssignedAsCorrect(tr,FramesToAnalyze);

for frame_ind = 1:length(FramesToAnalyze)
    Frame = FramesToAnalyze(frame_ind);
    if FramesFlaggedAsCorrectForCalculation(Frame) % if this frame doesn't need correction 
        Position = squeeze(data.NeuronCoordinatesMatrix_Corrected(tr,Frame,:))';
    else                                        % if this frame needs correction
        IndexInTracks = data.ConvertMovieFrameToIndexWithinTracks(tr,Frame);
        NormFrame     = ReadAndNormalizeFrame(tr,Frame, File, background, VignettingPattern, WormArea, IndexInTracks);

        [Position, Distance3] = FindLocalMaximum (NormFrame,PreviousPosition, MaxPixelsPerFrame); % 2 pixels
        Position(Position<MinimumPosition)= MinimumPosition;  
        Position(1)                        = min([Position(1)  MaximumPosition_vec(1)]);
        Position(2)                        = min([Position(2)  MaximumPosition_vec(2)]);          
        data.NeuronCoordinatesMatrix_Corrected(tr,Frame,:) = Position;
        data.NeuronDisplayMatrices_Corrected(tr,Frame,:,:) = Recalculate_DisplayMatrix(NormFrame, Position, DisplaySquareSize);
    end
    PreviousPosition = Position;    
end   

guidata(Fig_handle, data);
disp('position re-calculation is finished');

return

function DisregardNeuronPositionWithinSegment(stam, stam2, Fig_handle, tr)
% Give NaN values to neuron coordinates and values within the frame interval (h_Correct_FirstFrame,h_Correct_LastFrame) then move to the next segment.  
% NOTE: Frames within the current segment that are not within the above frame interval are maintained!!    
data                     = guidata(Fig_handle); 
HANDLES                  = data.HANDLES;
FirstCorrectionFrame = str2num(HANDLES.h_Correct_FirstFrame(tr).String);
LastCorrectionFrame  = str2num(HANDLES.h_Correct_LastFrame(tr).String);

FramesToAnalyze        = FirstCorrectionFrame:LastCorrectionFrame;
if isempty(FramesToAnalyze)||(FirstCorrectionFrame==0)||(LastCorrectionFrame==0)
    return
end
data.NeuronCoordinatesMatrix_Corrected(tr,FramesToAnalyze,:) = NaN;
data.NeuronDisplayMatrices_Corrected(tr,FramesToAnalyze,:,:) = NaN;
data.ManuallyAssignedAsDisregard(tr,FramesToAnalyze)         = true;
data.FramesForNeuronPositionReCalculation(tr,:)              = false;  
guidata(Fig_handle, data);

return

function [Position, Distance] = FindLocalMaximum (Mov,NeuronCoordinates_Estimations, MaximalDistanceInPixels)
MaximalDistanceInPixels1 = MaximalDistanceInPixels;
MaximalDistanceInPixels2 = MaximalDistanceInPixels;

while ((NeuronCoordinates_Estimations(1)-MaximalDistanceInPixels1)<1)||((NeuronCoordinates_Estimations(1)+MaximalDistanceInPixels1)>size(Mov,1))
    MaximalDistanceInPixels1 = MaximalDistanceInPixels1-1;
end
while ((NeuronCoordinates_Estimations(2)-MaximalDistanceInPixels2)<1)||((NeuronCoordinates_Estimations(2)+MaximalDistanceInPixels2)>size(Mov,2))
    MaximalDistanceInPixels2 = MaximalDistanceInPixels2-1;
end
    
MAT = Mov((NeuronCoordinates_Estimations(1)-MaximalDistanceInPixels1):(NeuronCoordinates_Estimations(1)+MaximalDistanceInPixels1), ...
          (NeuronCoordinates_Estimations(2)-MaximalDistanceInPixels2):(NeuronCoordinates_Estimations(2)+MaximalDistanceInPixels2));
MATcenter = [MaximalDistanceInPixels1  MaximalDistanceInPixels2]+1;
if all(isnan(MAT(:)))
    I = MATcenter(1);
    J = MATcenter(2);
else
    [M,I]     = max(MAT(:));
    [I,J]     = ind2sub(size(MAT),I);
end

Distance = sqrt(sum( (MATcenter - [I,J]).^2 ));
Position = round(NeuronCoordinates_Estimations) + [I,J] - MATcenter;

return

function ClearNonConfirmedCorrections(stam, stam2, Fig_handle, tr)
data      = guidata(Fig_handle); 
HANDLES   = data.HANDLES;

data.FramesForNeuronPositionReCalculation(tr,:) = false;
guidata(Fig_handle, data);

set(HANDLES.h_Correct_FirstFrame(tr),'Value',0,'String','0');
set(HANDLES.h_Correct_LastFrame(tr),'Value',0,'String','0');

return

function SaveSessionBackup(stam, stam2, Fig_handle,File)
data                  = guidata(Fig_handle); 
data_without_handles  = rmfield(data,'HANDLES'); 
DateString            = datestr(now,'yymmdd_HHMM');
FileName              = [File.MoviePath,'\',File.MovieName,'_NeuronPositionCorrection_SessionBackup_',DateString,'.mat'];
save(FileName,'data_without_handles','-v7.3');
disp(['Backup session was saved at: ',FileName]);
return

function LoadPreviousSession(stam, stam2, Fig_handle,File)

data    = guidata(Fig_handle); 
HANDLES = data.HANDLES;

DialogTitle         = 'Select a session backup file (NeuronPositionCorrection_SessionBackup_yymmdd_HHMM'; 
DefaultFileName     = [File.MoviePath,'\',File.MovieName,'_NeuronPositionCorrection_SessionBackup_*.mat'];
[FileName,PathName] = uigetfile('*.mat',DialogTitle, DefaultFileName);
if ~isempty(find(FileName==0,1))
    return
end
BackupFileName      = [PathName,FileName];
load(BackupFileName,'data_without_handles')
data                = data_without_handles;
data.HANDLES        = HANDLES;
guidata(Fig_handle,data);

UpdateTrackSegmentationFields(Fig_handle);
disp(['Backup session was loaded from: ',BackupFileName]);

return

function UpdateTrackSegmentationFields(Fig_handle)

data              = guidata(Fig_handle); 
HANDLES           = data.HANDLES;
NumOfTracks       = size(data.FramesForManualCorrection,1);
SegmentsToCorrect = data.SegmentsToCorrect;
for tr=1:NumOfTracks
    ClearNonConfirmedCorrections([], [], Fig_handle, tr);
    
    NumberOfSegmentsLeft = length(SegmentsToCorrect(tr).InitialFrame);    
    if NumberOfSegmentsLeft>=1
        NextStartFrame    = SegmentsToCorrect(tr).InitialFrame(1);
        NextEndFrame      = SegmentsToCorrect(tr).EndFrame(1);
        FrameIntervalText = ['Relevant frame interval= [',num2str(NextStartFrame),' - ',num2str(NextEndFrame),']'];
    else
        FrameIntervalText                              = 'No flagged frames';
    end
    set(HANDLES.h_CorrectionFrameInterval_text(tr),'String',FrameIntervalText);
end

return

function UpdateAndSaveTracksFile(stam, stam2, Fig_handle,File, Tracks, BehaviorDataMatrices, background, VignettingPattern, WormArea)
data     = guidata(Fig_handle); 
data     = rmfield(data,'HANDLES'); 
FileName = [File.MoviePath,'\',File.MovieName,'_AfterNeuronPositionManualCorrection.mat'];
for tr=1:length(Tracks)
    Tracks(tr).Neuron.CoordinatesMatrix = squeeze(data.NeuronCoordinatesMatrix_Corrected(tr,:,:));
    Tracks(tr).Neuron.DisplayMatrix     = squeeze(data.NeuronDisplayMatrices_Corrected(tr,:,:,:));
end

plotme              = true;
Tracks              = ExtractNeuralActivityFromTracks_AfterGUI_02(Tracks, File, plotme); %% EXTERNAL FUNCTION
tic
save(FileName,'data','File','Tracks','BehaviorDataMatrices','background','VignettingPattern','WormArea','-v7.3');
toc
disp(['Tracks were corrected based on manual correction and data was saved at: ',FileName]);
return

















