function handles = AlignAndCrop(handles)

% Help for the Align and Crop module:
% Category: Pre-processing
% 
% This module was written for a specific purpose: to take a larger
% image of the outlines of cells, with some space around the outside,
% and align it with a real image of cells (or of nuclei).
% Fortunately, it can work for almost any set of images for which one
% set is larger than the other and rescaling must be done (and, of
% course, there are features shared which can be used to align the
% images).  Just enter the larger images as the traced ones, and the
% smaller ones as the real images.
% 
% This module determines the best alignment of two images.  It expects
% the first to be a scanned, traced transparency of 8 1/2 by 11
% inches, and the second should be a real image from a microscope.
% The addition of colored registration marks is optional, although
% very helpful in alignment.  The module stretches the real image to
% be at the same resolution as the traced image, then creates a black
% background for it of the same size, then aligns the images, and
% finally crops them all to be the same size and taking up the same
% area as on the original real image. In the future, this may be
% developed to include a way of finding the optimal scaling of the
% images to each other; however, no guarantees.
% 
% Note that as long as the input real images of nuclei and of cells
% are the same dimensions, their output files will have the same
% dimensions as well.
%
% SAVING IMAGES: The aligned, cropped traced images and the aligned,
% cropped real images produced by this module can be easily saved
% using the Save Images module, using the names you assign. If you
% want to save other intermediate images, alter the code for this
% module to save those images to the handles structure (see the
% SaveImages module help) and then use the Save Images module.
%
% See also ALIGN.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
% 
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
% 
% Authors:
%   Anne Carpenter <carpenter@wi.mit.edu>
%   Thouis Jones   <thouis@csail.mit.edu>
%   In Han Kang    <inthek@mit.edu>
%
% $Revision$

% PROGRAMMING NOTE
% HELP:
% The first unbroken block of lines will be extracted as help by
% CellProfiler's 'Help for this analysis module' button as well as Matlab's
% built in 'help' and 'doc' functions at the command line. It will also be
% used to automatically generate a manual page for the module. An example
% image demonstrating the function of the module can also be saved in tif
% format, using the same name as the module, and it will automatically be
% included in the manual page as well.  Follow the convention of: purpose
% of the module, description of the variables and acceptable range for
% each, how it works (technical description), info on which images can be 
% saved, and See also CAPITALLETTEROTHERMODULES. The license/author
% information should be separated from the help lines with a blank line so
% that it does not show up in the help displays.  Do not change the
% programming notes in any modules! These are standard across all modules
% for maintenance purposes, so anything module-specific should be kept
% separate.

% PROGRAMMING NOTE
% DRAWNOW:
% The 'drawnow' function allows figure windows to be updated and
% buttons to be pushed (like the pause, cancel, help, and view
% buttons).  The 'drawnow' function is sprinkled throughout the code
% so there are plenty of breaks where the figure windows/buttons can
% be interacted with.  This does theoretically slow the computation
% somewhat, so it might be reasonable to remove most of these lines
% when running jobs on a cluster where speed is important.
drawnow

%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%
drawnow

% PROGRAMMING NOTE
% VARIABLE BOXES AND TEXT: 
% The '%textVAR' lines contain the variable descriptions which are
% displayed in the CellProfiler main window next to each variable box.
% This text will wrap appropriately so it can be as long as desired.
% The '%defaultVAR' lines contain the default values which are
% displayed in the variable boxes when the user loads the module.
% The line of code after the textVAR and defaultVAR extracts the value
% that the user has entered from the handles structure and saves it as
% a variable in the workspace of this module with a descriptive
% name. The syntax is important for the %textVAR and %defaultVAR
% lines: be sure there is a space before and after the equals sign and
% also that the capitalization is as shown. 
% CellProfiler uses VariableRevisionNumbers to help programmers notify
% users when something significant has changed about the variables.
% For example, if you have switched the position of two variables,
% loading a pipeline made with the old version of the module will not
% behave as expected when using the new version of the module, because
% the settings (variables) will be mixed up. The line should use this
% syntax, with a two digit number for the VariableRevisionNumber:
% '%%%VariableRevisionNumber = 01'  If the module does not have this
% line, the VariableRevisionNumber is assumed to be 00.  This number
% need only be incremented when a change made to the modules will affect
% a user's previously saved settings. There is a revision number at
% the end of the license info at the top of the m-file for revisions
% that do not affect the user's previously saved settings files.

%%% Reads the current module number, because this is needed to find 
%%% the variable values that the user entered.
CurrentModule = handles.Current.CurrentModuleNumber;
CurrentModuleNum = str2double(CurrentModule);

%infotypeVAR01 = imagegroup
%textVAR01 = What did you call the traced images?
TracedImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%infotypeVAR02 = imagegroup
%textVAR02 = What did you call the real images?
%defaultVAR02 = OrigReal
RealImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%infotypeVAR03 = imagegroup indep
%textVAR03 = What do you want to call the aligned, cropped traced images?
%defaultVAR03 = ACTraced
FinishedTracedImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%infotypeVAR04 = imagegroup indep
%textVAR04 = What do you want to call the aligned, cropped real images?
%defaultVAR04 = ACReal
FinishedRealImageName = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Enter the printed size of the real image in inches as "height,width" (no quotes).
%choiceVAR05 = 11,8.5
PrintedImageSize = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu custom

%textVAR06 = Enter the page orientation of the traced images (portrait or landscape)
%choiceVAR06 = portrait
%choiceVAR06 = landscape
Orientation = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Separates the entered dimensions into two variables
PrintedImageSizeNumerical = str2double(PrintedImageSize);
PrintedHeight = PrintedImageSizeNumerical(1);
PrintedWidth = PrintedImageSizeNumerical(2);

%%% Checks that the page orientation is a valid entry
if strcmp(upper(Orientation),'PORTRAIT') ~= 1
    if strcmp(upper(Orientation),'LANDSCAPE') ~= 1
        error('Image processing was canceled because you have not entered a valid response for the page orientaion in the AlignAndCrop module.  Type either landscape or portrait and try again.')
    end
end

%%% Checks whether the image to be analyzed exists in the handles structure.
if isfield(handles.Pipeline, TracedImageName) == 0
    %%% If the image is not there, an error message is produced.  The error
    %%% is not displayed: The error function halts the current function and
    %%% returns control to the calling function (the analyze all images
    %%% button callback.)  That callback recognizes that an error was
    %%% produced because of its try/catch loop and breaks out of the image
    %%% analysis loop without attempting further modules.
    error(['Image processing was canceled because the AlignAndCrop module could not find the input image.  It was supposed to be named ', TracedImageName, ' but an image with that name does not exist.  Perhaps there is a typo in the name.'])
end
%%% Reads the image.
TracedImage = handles.Pipeline.(TracedImageName);
% figure, imshow(TracedImage), title('TracedImage')

%%% Checks whether the image to be analyzed exists in the handles structure.
if isfield(handles, RealImageName) == 0
    %%% If the image is not there, an error message is produced.  The error
    %%% is not displayed: The error function halts the current function and
    %%% returns control to the calling function (the analyze all images
    %%% button callback.)  That callback recognizes that an error was
    %%% produced because of its try/catch loop and breaks out of the image
    %%% analysis loop without attempting further modules.
    error(['Image processing was canceled because the AlignAndCrop module could not find the input image.  It was supposed to be named ', RealImageName, ' but an image with that name does not exist.  Perhaps there is a typo in the name.'])
end
%%% Reads the image.
RealImage = handles.Pipeline.(RealImageName);
% figure, imshow(RealImage), title('RealImage')

%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%
drawnow

% PROGRAMMING NOTE
% TO TEMPORARILY SHOW IMAGES DURING DEBUGGING: 
% figure, imshow(BlurredImage, []), title('BlurredImage') 
% TO TEMPORARILY SAVE IMAGES DURING DEBUGGING: 
% imwrite(BlurredImage, FileName, FileFormat);
% Note that you may have to alter the format of the image before
% saving.  If the image is not saved correctly, for example, try
% adding the uint8 command:
% imwrite(uint8(BlurredImage), FileName, FileFormat);
% To routinely save images produced by this module, see the help in
% the SaveImages module.

%%% creates a two-dimensional, b&w Traced Image to work with during this section
ExTracedImage = TracedImage(:,:,1);
PreTracedImage = im2bw(ExTracedImage,.5);
%%% creates a 2D real image to work with for the size
ExRealImage = RealImage(:,:,1);
PreRealImage = im2bw(ExRealImage,.5);
%%% finds the size of each of the images
[TracedY,TracedX] = size(PreTracedImage);
[RealY,RealX] = size(PreRealImage);
    %CPwarndlg(['The dimensions are as follows: TY=',num2str(TracedY),' TX=',num2str(TracedX),' RY=',num2str(RealY),' RX=',num2str(RealX)],'Alert!')
%%% checks that for both dimensions, the traced image is bigger than the
%%% real one
drawnow
if TracedX < RealX
    error('Image processing was canceled in the AlignAndCrop module because the real image is wider than the traced image. Make sure that the resolution is the same for both images.')
elseif TracedY < RealY
    error('Image processing was canceled in the AlignAndCrop module because the real image is taller than the traced image.  Make sure that the resolution is the same for both images.')
end

%%% Finds the resolution along both axes for both sets of images
RealResolution = RealX / PrintedWidth;
RealResolution2 = RealY / PrintedHeight;
if strcmp(upper(Orientation),'PORTRAIT') == 1
    TracedResolution = TracedX /8.5;
    TracedResolution2 = TracedY /11;
        %CPwarndlg('I think that you wrote portrait')
else TracedResolution = TracedX / (11);
    TracedResolution2 = TracedY / (8.5);
        %CPwarndlg('I think that you wrote landscape')
end
%%% activate this line to pop up a dialog box with the resolution values in it
    %CPwarndlg(['RR is ',num2str(RealResolution),',RR2 is ',num2str(RealResolution2),', TR is ',num2str(TracedResolution),', TR2 is ',num2str(TracedResolution2)],'Alert!')
drawnow

if RealResolution ~= RealResolution2
    error('Oops.  Due to some fundamental problem in coding of the AlignAndCrop module or in you, the resolution is different depending upon which axis you calculate it.  Re-measure your printed image: you may be trying to be too accurate.')
end
if TracedResolution ~= TracedResolution2
    error('Oops.  Due to some fundamental problem in coding of the AlignAndCrop or in you, the resolution is different depending upon which axis you calculate it.  Most likely, you scanned the image not perfectly at 8.5 by 11 and caused major problems.')
end
drawnow
%%% Resizes the real image to be at the same resolution as the traced image
ResizingCoefficient = TracedResolution / RealResolution;
    %CPwarndlg(['Resizing Coefficient = ',num2str(ResizingCoefficient)],'NFO 4 U')
ResizedRealImage = imresize(RealImage,ResizingCoefficient,'bicubic');
%%% finds the dimensions of the resized real image
[RRealX,RRealY] = size(ResizedRealImage);
    %CPwarndlg(['RRealX = ',num2str(RRealX),', RRealY = ',num2str(RRealY)],'Yeahhhhh, About That...')
%%% finds the difference in dimensions to create a margin value
XDifference = TracedX - RRealX;
YDifference = TracedY - RRealY;
XMargin = XDifference / 2;
YMargin = YDifference / 2;

drawnow
%%% creates a matrix of zeros the size of the traced image
RealImageBlankSlate = zeros(TracedY,TracedX);
%%% makes sure images are the right size
NewImageSizeX = TracedX - 2*(XMargin); 
NewImageSizeY = TracedY - 2*(YMargin);
if isequal(NewImageSizeX,RRealX) == 0
    error('The value (TracedX - 2(XMargin)) is not the same as the value (RRealX) in the AlignAndCrop module')
elseif isequal(NewImageSizeY,RRealY) == 0
    error('The value (TracedY - 2(YMargin)) is not the same as the value (RRealY) in the AlignAndCrop')
end
drawnow
%%% rounds each number up, then subtracts the original from that and tests
%%% equality with zero to see if the numbers are not whole, which then pops
%%% up messages but does not abort. then it pastes the resized real image
%%% over the blank image
RXM = round(XMargin); 
RYM = round(YMargin);
XMD = RXM - XMargin; 
YMD = RYM - YMargin;
if isequal(XMD,0) == 0
    if isequal(YMD,0) == 0
        CPwarndlg(['Warning: neither margin value is an integer, and therefore problems may appear in the very near future. XM = ',num2str(XMargin),', YM = ',num2str(YMargin)],'Feeling Unwhole')
        YMargin1 = ceil(YMargin); 
        YMargin2 = floor(YMargin);
    else CPwarndlg('Warning: the XMargin number is not an integer, and therefore may cause problems in the very near future. Fortunately, the YMargin value is an integer','Feeling Unwhole')
    end
    XMargin1 = ceil(XMargin); 
    XMargin2 = floor(XMargin);
elseif isequal(YMD,0) == 0
    CPwarndlg('Warning: the YMargin number is not an integer, and therefore may cause problems in the very near future. Fortunately, the XMargin value is an integer','Feeling Unwhole')
    YMargin1 = ceil(YMargin); 
    YMargin2 = floor(YMargin);
else XMargin1 = XMargin; 
    XMargin2 = XMargin;
    YMargin1 = YMargin; 
    YMargin2 = YMargin;
end
YBegin = YMargin2 + 1;
YEnd = TracedY - YMargin1;
XBegin = XMargin2 + 1;
XEnd = TracedX - XMargin1;
RealImageBlankSlate(YBegin:YEnd,XBegin:XEnd) = ResizedRealImage;
ExpandedRealImage = RealImageBlankSlate;

drawnow
%%% makes a 2D traced image to work with
TwoDTracedImage = rgb2gray(TracedImage);
%%% runs the alignment subfunctions on the two images, working with just
%%% one layer but applying movements to the whole image
InvertedTraced = imcomplement(TwoDTracedImage);
[sx,sy] = autoalign(InvertedTraced,ExpandedRealImage);
    %CPwarndlg('The autoalign step has completed','Notice:')
AlignedTracedImage = subim(TwoDTracedImage, sx, sy);
AlignedRealImage = subim(ExpandedRealImage, -sx, -sy);
    %CPwarndlg('The subim steps have completed','Notice')
% Results = ['(Traced vs. Real: X ',num2str(sx),', Y ',num2str(sy),')'];
    %CPwarndlg(['All image processing has completed. Results are ',Results],'Notice:')
%%% Checks that the size of aligned images is the same
if isequal(size(AlignedTracedImage),size(AlignedRealImage)) == 0
    error('After the alignment step was completed in the AlignAndCrop module, the two images were different sizes for some reason.  This is not good.')
end
%%% finds the end and begin points for all dimensions of the new images and
%%% creates the cropping rectangle matrix
YBeginARI = YBegin + sy;
YEndARI = YEnd + sy;
XBeginARI = XBegin + sx;
XEndARI = XEnd + sx;
CropRectARI = [XBeginARI YBeginARI (XEndARI - XBeginARI) (YEndARI - YBeginARI)];
%%% crops both aligned images
CroppedAlignedTracedImage = imcrop(AlignedTracedImage,CropRectARI);
CroppedAlignedRealImage = imcrop(AlignedRealImage,CropRectARI);

%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

% PROGRAMMING NOTE
% DISPLAYING RESULTS:
% Some calculations produce images that are used only for display or
% for saving to the hard drive, and are not used by downstream
% modules. To speed processing, these calculations are omitted if the
% figure window is closed and the user does not want to save the
% images.

%%% Determines the figure number to display in.
fieldname = ['FigureNumberForModule',CurrentModule];
ThisModuleFigureNumber = handles.Current.(fieldname);
if any(findobj == ThisModuleFigureNumber) == 1;
% PROGRAMMING NOTE
% DRAWNOW BEFORE FIGURE COMMAND:
% The "drawnow" function executes any pending figure window-related
% commands.  In general, Matlab does not update figure windows until
% breaks between image analysis modules, or when a few select commands
% are used. "figure" and "drawnow" are two of the commands that allow
% Matlab to pause and carry out any pending figure window- related
% commands (like zooming, or pressing timer pause or cancel buttons or
% pressing a help button.)  If the drawnow command is not used
% immediately prior to the figure(ThisModuleFigureNumber) line, then
% immediately after the figure line executes, the other commands that
% have been waiting are executed in the other windows.  Then, when
% Matlab returns to this module and goes to the subplot line, the
% figure which is active is not necessarily the correct one. This
% results in strange things like the subplots appearing in the timer
% window or in the wrong figure window, or in help dialog boxes.
drawnow
    CPfigure(handles,ThisModuleFigureNumber);
    subplot(2,2,1); imagesc(TracedImage); colormap(gray);
    title(['Traced Input, Image Set # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,2,2); imagesc(RealImage); title('Real Input Image');
    subplot(2,2,3); imagesc(CroppedAlignedTracedImage); colormap(gray); title('Cropped & Aligned Traced Image');
    subplot(2,2,4); imagesc(CroppedAlignedRealImage);title('Cropped & Aligned Real Image');
    CPFixAspectRatio(RealImage);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

% PROGRAMMING NOTE
% HANDLES STRUCTURE:
%       In CellProfiler (and Matlab in general), each independent
% function (module) has its own workspace and is not able to 'see'
% variables produced by other modules. For data or images to be shared
% from one module to the next, they must be saved to what is called
% the 'handles structure'. This is a variable, whose class is
% 'structure', and whose name is handles. The contents of the handles
% structure are printed out at the command line of Matlab using the
% Tech Diagnosis button. The only variables present in the main
% handles structure are handles to figures and gui elements.
% Everything else should be saved in one of the following
% substructures:
%
% handles.Settings:
%       Everything in handles.Settings is stored when the user uses
% the Save pipeline button, and these data are loaded into
% CellProfiler when the user uses the Load pipeline button. This
% substructure contains all necessary information to re-create a
% pipeline, including which modules were used (including variable
% revision numbers), their setting (variables), and the pixel size.
%   Fields currently in handles.Settings: PixelSize, ModuleNames,
% VariableValues, NumbersOfVariables, VariableRevisionNumbers.
%
% handles.Pipeline:
%       This substructure is deleted at the beginning of the
% analysis run (see 'Which substructures are deleted prior to an
% analysis run?' below). handles.Pipeline is for storing data which
% must be retrieved by other modules. This data can be overwritten as
% each image set is processed, or it can be generated once and then
% retrieved during every subsequent image set's processing, or it can
% be saved for each image set by saving it according to which image
% set is being analyzed, depending on how it will be used by other
% modules. Any module which produces or passes on an image needs to
% also pass along the original filename of the image, named after the
% new image name, so that if the SaveImages module attempts to save
% the resulting image, it can be named by appending text to the
% original file name.
%   Example fields in handles.Pipeline: FileListOrigBlue,
% PathnameOrigBlue, FilenameOrigBlue, OrigBlue (which contains the actual image).
%
% handles.Current:
%       This substructure contains information needed for the main
% CellProfiler window display and for the various modules to
% function. It does not contain any module-specific data (which is in
% handles.Pipeline).
%   Example fields in handles.Current: NumberOfModules,
% StartupDirectory, DefaultOutputDirectory, DefaultImageDirectory,
% FilenamesInImageDir, CellProfilerPathname, ImageToolHelp,
% DataToolHelp, FigureNumberForModule01, NumberOfImageSets,
% SetBeingAnalyzed, TimeStarted, CurrentModuleNumber.
%
% handles.Preferences: 
%       Everything in handles.Preferences is stored in the file
% CellProfilerPreferences.mat when the user uses the Set Preferences
% button. These preferences are loaded upon launching CellProfiler.
% The PixelSize, DefaultImageDirectory, and DefaultOutputDirectory
% fields can be changed for the current session by the user using edit
% boxes in the main CellProfiler window, which changes their values in
% handles.Current. Therefore, handles.Current is most likely where you
% should retrieve this information if needed within a module.
%   Fields currently in handles.Preferences: PixelSize, FontSize,
% DefaultModuleDirectory, DefaultOutputDirectory,
% DefaultImageDirectory.
%
% handles.Measurements
%      Data extracted from input images are stored in the
% handles.Measurements substructure for exporting or further analysis.
% This substructure is deleted at the beginning of the analysis run
% (see 'Which substructures are deleted prior to an analysis run?'
% below). The Measurements structure is organized in two levels. At
% the first level, directly under handles.Measurements, there are
% substructures (fields) containing measurements of different objects.
% The names of the objects are specified by the user in the Identify
% modules (e.g. 'Cells', 'Nuclei', 'Colonies').  In addition to these
% object fields is a field called 'Image' which contains information
% relating to entire images, such as filenames, thresholds and
% measurements derived from an entire image. That is, the Image field
% contains any features where there is one value for the entire image.
% As an example, the first level might contain the fields
% handles.Measurements.Image, handles.Measurements.Cells and
% handles.Measurements.Nuclei.
%      In the second level, the measurements are stored in matrices 
% with dimension [#objects x #features]. Each measurement module
% writes its own block; for example, the MeasureAreaShape module
% writes shape measurements of 'Cells' in
% handles.Measurements.Cells.AreaShape. An associated cell array of
% dimension [1 x #features] with suffix 'Features' contains the names
% or descriptions of the measurements. The export data tools, e.g.
% ExportData, triggers on this 'Features' suffix. Measurements or data
% that do not follow the convention described above, or that should
% not be exported via the conventional export tools, can thereby be
% stored in the handles.Measurements structure by leaving out the
% '....Features' field. This data will then be invisible to the
% existing export tools.
%      Following is an example where we have measured the area and
% perimeter of 3 cells in the first image and 4 cells in the second
% image. The first column contains the Area measurements and the
% second column contains the Perimeter measurements.  Each row
% contains measurements for a different cell:
% handles.Measurements.Cells.AreaShapeFeatures = {'Area' 'Perimeter'}
% handles.Measurements.Cells.AreaShape{1} = 	40		20
%                                               100		55
%                                              	200		87
% handles.Measurements.Cells.AreaShape{2} = 	130		100
%                                               90		45
%                                               100		67
%                                               45		22
%
% Which substructures are deleted prior to an analysis run?
%       Anything stored in handles.Measurements or handles.Pipeline
% will be deleted at the beginning of the analysis run, whereas
% anything stored in handles.Settings, handles.Preferences, and
% handles.Current will be retained from one analysis to the next. It
% is important to think about which of these data should be deleted at
% the end of an analysis run because of the way Matlab saves
% variables: For example, a user might process 12 image sets of nuclei
% which results in a set of 12 measurements ("ImageTotalNucArea")
% stored in handles.Measurements. In addition, a processed image of
% nuclei from the last image set is left in the handles structure
% ("SegmNucImg"). Now, if the user uses a different algorithm which
% happens to have the same measurement output name "ImageTotalNucArea"
% to analyze 4 image sets, the 4 measurements will overwrite the first
% 4 measurements of the previous analysis, but the remaining 8
% measurements will still be present. So, the user will end up with 12
% measurements from the 4 sets. Another potential problem is that if,
% in the second analysis run, the user runs only a module which
% depends on the output "SegmNucImg" but does not run a module that
% produces an image by that name, the module will run just fine: it
% will just repeatedly use the processed image of nuclei leftover from
% the last image set, which was left in handles.Pipeline.

%%% The adjusted image is saved to the handles structure so it can be used
%%% by subsequent modules.
handles.Pipeline.(FinishedTracedImageName) = CroppedAlignedTracedImage;
handles.Pipeline.(FinishedRealImageName) = CroppedAlignedRealImage;

%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%

%%% Written by Thouis R. Jones in the Align module

function [shiftx, shifty] = autoalign(in1, in2)
%%% Aligns two images using mutual-information and hill-climbing.
best = mutualinf(in1, in2);
bestx = 0;
besty = 0;
%%% Checks which one-pixel move is best.
for dx=-1:1,
    for dy=-1:1,
        cur = mutualinf(subim(in1, dx, dy), subim(in2, -dx, -dy));
        if (cur > best),
            best = cur;
            bestx = dx;
            besty = dy;
        end
    end
end
if (bestx == 0) && (besty == 0),
    shiftx = 0;
    shifty = 0;
    return;
end
%%% Remembers the lastd direction we moved.
lastdx = bestx;
lastdy = besty;
%%% Loops until things stop improving.
while true,
    [nextx, nexty, newbest] = one_step(in1, in2, bestx, besty, lastdx, lastdy, best);
    if (nextx == 0) && (nexty == 0),
        shiftx = bestx;
        shifty = besty;
        return;
    else
        bestx = bestx + nextx;
        besty = besty + nexty;
        best = newbest;
    end
end

function [nx, ny, nb] = one_step(in1, in2, bx, by, ldx, ldy, best)
%%% Finds the best one pixel move, but only in the same direction(s) we
%%% moved last time (no sense repeating evaluations)
nb = best;
for dx=-1:1,
    for dy=-1:1,
        if (dx == ldx) || (dy == ldy),
            cur = mutualinf(subim(in1, bx+dx, by+dy), subim(in2, -(bx+dx), -(by+dy)));
            if (cur > nb),
                nb = cur;
                nx = dx;
                ny = dy;
            end
        end
    end
end
if (best == nb),
    %%% no change, so quit searching
    nx = 0;
    ny = 0;
end

function sub = subim(im, dx, dy)
%%% Subimage with positive or negative offsets
if (dx > 0),
    sub = im(:,dx+1:end);
else
    sub = im(:,1:end+dx);
end
if (dy > 0),
    sub = sub(dy+1:end,:);
else
    sub = sub(1:end+dy,:);
end

function H = entropy(X)
%%% Entropy of samples X
S = imhist(X,256);
%%% if S is probability distribution function N is 1
N=sum(sum(S));
if ((N>0) && (min(S(:))>=0))
    Snz=nonzeros(S);
    H=log2(N)-sum(Snz.*log2(Snz))/N;
else
    H=0;
end

function H = entropy2(X,Y)
%%% joint entropy of paired samples X and Y
%%% Makes sure images are binned to 256 graylevels
X = double(im2uint8(X));
Y = double(im2uint8(Y));
%%% Creates a combination image of X and Y
XY = 256*X + Y;
S = histc(XY(:),0:(256*256-1));
%%% If S is probability distribution function N is 1
N=sum(sum(S));          
if ((N>0) && (min(S(:))>=0))
    Snz=nonzeros(S);
    H=log2(N)-sum(Snz.*log2(Snz))/N;
else
    H=0;
end

function I = mutualinf(X, Y)
%%% Mutual information of images X and Y
I = entropy(X) + entropy(Y) - entropy2(X,Y);