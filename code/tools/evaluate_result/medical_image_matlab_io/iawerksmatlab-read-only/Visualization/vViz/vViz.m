function varargout = vViz(varargin)
% VVIZ M-file for vViz.fig
%      VVIZ, by itself, creates a new VVIZ or raises the existing
%      singleton*.
%
%      H = VVIZ returns the handle to a new VVIZ or the handle to
%      the existing singleton*.
%
%      VVIZ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VVIZ.M with the given input arguments.
%
%      VVIZ('Property','Value',...) creates a new VVIZ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vViz_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vViz_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vViz

% Last Modified by GUIDE v2.5 19-Jan-2007 23:10:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vViz_OpeningFcn, ...
                   'gui_OutputFcn',  @vViz_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before vViz is made visible.
function vViz_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vViz (see VARARGIN)

% Choose default command line output for vViz
handles.output = hObject;

handles.defaultDim = [50 50 50];
handles.dim = handles.defaultDim;

% load default image
handles.im = createDefaultImage(handles.dim);
handles.imOrigin = [0 0 0];
handles.imSpacing= [1 1 1];

% initialize intensity range and windowing
handles.imageWindowMin = min(handles.im(:));
handles.imageWindowMax = max(handles.im(:));
handles.imageMin = handles.imageWindowMin;
handles.imageMax = handles.imageWindowMax;
updateImageWindowSlider(handles);

% initially there is no field
handles.haveField = false;
handles.vField = [];
handles.detJ = [];
handles.fieldStart = [1 1 1];
handles.fieldSize = [3 0 0 0];

handles.fieldWindowMin = 0;
handles.fieldWindowMax = 0;
handles.fieldMin = 0;
handles.fieldMax = 0+eps;
updateFieldWindowSlider(handles);

handles.orientationStrs = {'sagittal','coronal','axial'}
handles.orientation = 3;
handles.slices = round(handles.dim/2);
%handles.slice = round(handles.dim(3)/2);
handles.minSlice = 1;
handles.maxSlice = handles.dim(3);

handles.majorGrid = false;
handles.minorGrid = false;

set(handles.versionTag,'String','v. 0.6');
set(handles.axialToggle,'Value',1);
set(handles.sagittalToggle,'Value',0);
set(handles.coronalToggle,'Value',0);

set(handles.view2DToggle,'Value',1);
set(handles.view3DToggle,'Value',0);

set(handles.vFieldToggle,'Value',1);
set(handles.hFieldToggle,'Value',0);

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using vViz.
if strcmp(get(hObject,'Visible'),'off')
    updateDisplay(handles);
end

updateImageWindowSlider(handles);
updateSliceSlider(handles);

% UIWAIT makes vViz wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = vViz_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function updateDisplay(handles)
axes(handles.axes1);

fprintf('image orientation: %s\n',handles.orientationStrs{handles.orientation});
fprintf('slice: %d\n',handles.slices);
fprintf('image range: %g--%g\n',handles.imageMin,handles.imageMax);
fprintf('image window: %g--%g\n',handles.imageWindowMin,handles.imageWindowMax);
fprintf('field range: %g--%g\n',handles.fieldMin,handles.fieldMax);
fprintf('field window: %g--%g\n',handles.fieldWindowMin,handles.fieldWindowMax);

if get(handles.view2DToggle,'Value')
  % get image slice
  fprintf('Loading image slice...\n');
  si = getImageSlice(handles);
  
  % window image slice
  if (handles.imageWindowMax > handles.imageWindowMin)
    si = 1 + (si-handles.imageWindowMin) * ...
      255/(handles.imageWindowMax-handles.imageWindowMin);
  end
  si(si<1) = 1;
  si(si>256) = 256;
  
  % get true color copy of image
  cmaps = get(handles.imageColormapPopup,'String');
  cmap = colormap([cmaps{get(handles.imageColormapPopup, 'Value')} '(256)']);
  iTC = reshape(cmap(round(si),:),[size(si) 3]);
  
  fieldViewStrs = get(handles.fieldViewPopup,'String');
  fieldView = fieldViewStrs{get(handles.fieldViewPopup,'Value')};
  if handles.haveField & ~strcmpi(fieldView,'Quiver 2D')
    fprintf('Loading field slice...\n');
    
    % get field slice
    sf = getFieldSliceImage(handles,getFieldSlice(handles));
    
    % window field slice
    sf = 1 + (sf-handles.fieldWindowMin) * ...
      255/(handles.fieldWindowMax-handles.fieldWindowMin);
    sf(sf<1) = 1;
    sf(sf>256) = 256;
    
    % get true color copy of field slice
    cmap = colormap([cmaps{get(handles.fieldColormapPopup, 'Value')}...
      '(256)']);
    fTC = reshape(cmap(round(sf),:),[size(sf) 3]);
    
    % blend field and image slice
    alpha = get(handles.imageFieldBlendSlider,'Value');
    iTC = (1-alpha)*iTC + alpha * fTC;
  end
  
  % display the image
  cla;
  image(permute(iTC,[2 1 3]));
  
  if handles.haveField & strcmpi(fieldView,'Quiver 2D')
    sv = getFieldSlice(handles);
    step = getQuiverSpacing(handles);
    scale = getQuiverScaling(handles);
    fprintf('Size: %d x %d\n',size(sv,3),size(sv,2));
    [y,x] = meshgrid(1:step:size(sv,3),1:step:size(sv,2));
    hold on;
    color = getQuiverColor(handles);
    
    switch handles.orientationStrs{handles.orientation}
      case 'axial'
        quiver(x,y,squeeze(sv(1,1:step:end,1:step:end)),...
          squeeze(sv(2,1:step:end,1:step:end)),scale,color);
      case 'sagittal'
        quiver(x,y,squeeze(sv(2,1:step:end,1:step:end)),...
          squeeze(sv(3,1:step:end,1:step:end)),scale,color);
      case 'coronal'
        quiver(x,y,squeeze(sv(1,1:step:end,1:step:end)),...
          squeeze(sv(3,1:step:end,1:step:end)),scale,color);
    end
    hold off;
  end
  
  %axis xy;
  
  % label axes
  switch handles.orientationStrs{handles.orientation}
    case 'axial'
      xlabel('X_1');
      ylabel('X_2');
    case 'sagittal'
      xlabel('X_2');
      ylabel('X_3');
    case 'coronal'
      xlabel('X_1');
      ylabel('X_3');
  end
elseif get(handles.view3DToggle,'Value') & handles.haveField
  % 3D quiver view
  fprintf('3D Quiver View...\n');
  sv = getFieldSlice(handles);
  step = getQuiverSpacing(handles);
  scale = getQuiverScaling(handles);
  [y,x] = meshgrid(1:step:size(sv,3),1:step:size(sv,2));
  color = getQuiverColor(handles);
  switch handles.orientationStrs{handles.orientation}
    case 'axial'
      plane = handles.slices(3)*ones(size(sv,2),size(sv,3));
      quiver3(x,y,plane(1:step:end,1:step:end),...
        squeeze(sv(1,1:step:end,1:step:end)),...
        squeeze(sv(2,1:step:end,1:step:end)),...
        squeeze(sv(3,1:step:end,1:step:end)),scale,color);
    case 'sagittal'
      plane = handles.slices(1)*ones(size(sv,2),size(sv,3));
      quiver3(plane(1:step:end,1:step:end),x,y,...
        squeeze(sv(1,1:step:end,1:step:end)),...
        squeeze(sv(2,1:step:end,1:step:end)),...
        squeeze(sv(3,1:step:end,1:step:end)),scale,color);
    case 'coronal'
      plane = handles.slices(2)*ones(size(sv,2),size(sv,3));
      quiver3(x,plane(1:step:end,1:step:end),y,...
        squeeze(sv(1,1:step:end,1:step:end)),...
        squeeze(sv(2,1:step:end,1:step:end)),...
        squeeze(sv(3,1:step:end,1:step:end)),scale,color);
  end
  axis([0 handles.dim(1)+1 0 handles.dim(2)+1 0 handles.dim(3)+1]);
  xlabel('X_1');
  ylabel('X_2');
  zlabel('X_3');
end

function si = getFieldSliceImage(handles,s)
strs = get(handles.fieldViewPopup,'String');
si = s(1,:,:);
switch strs{get(handles.fieldViewPopup,'Value')}
  case 'L_2 Norm'
    si = squeeze(sqrt(sum(s .* s)));
  case 'L_1 Norm'
    si = squeeze(sum(abs(s)));
  case 'L_inf Norm'
    si = squeeze(max(abs(s)));
  case 'V_1'
    si = squeeze(s(1,:,:));
  case 'V_2'
    si = squeeze(s(2,:,:));
  case 'V_3'
    si = squeeze(s(3,:,:));
  case 'det(Jacobian)'
    si = s;
end

function s = getFieldSlice(handles)
strs = get(handles.fieldViewPopup,'String');
% special case, s is just a slice of the detJacobian 'image'
if strcmpi(strs{get(handles.fieldViewPopup,'Value')},'det(Jacobian)')
  switch handles.orientationStrs{handles.orientation}
    case 'axial'
    s = ones(handles.dim(1),handles.dim(2));
    s(handles.fieldStart(1)+1:handles.fieldStart(1)+handles.fieldSize(2),...
      handles.fieldStart(2)+1:handles.fieldStart(2)+handles.fieldSize(3))...
      = handles.detJ(:,:,handles.slices(3));
    case 'sagittal'
      s = ones(handles.dim(2),handles.dim(3));
      s(handles.fieldStart(2)+1:handles.fieldStart(2)+handles.fieldSize(3),...
        handles.fieldStart(3)+1:handles.fieldStart(3)+handles.fieldSize(4))...
        = handles.detJ(handles.slices(1),:,:);
    case 'coronal'
      s = ones(handles.dim(1),handles.dim(3));
      s(handles.fieldStart(1)+1:handles.fieldStart(1)+handles.fieldSize(2),...
        handles.fieldStart(3)+1:handles.fieldStart(3)+handles.fieldSize(4))...
        = handles.detJ(:,handles.slices(2),:);
  end
else
  % get 2D slice of field (vectors are still 3D)
  switch handles.orientationStrs{handles.orientation}
    case 'axial'
      s = zeros(3,handles.dim(1),handles.dim(2));
      fieldSlice = handles.slices(3) - handles.fieldStart(3);
      if fieldSlice > 0 & fieldSlice <= handles.fieldSize(4)
        s(:,...
          handles.fieldStart(1)+1:handles.fieldStart(1)+handles.fieldSize(2),...
          handles.fieldStart(2)+1:handles.fieldStart(2)+handles.fieldSize(3))...
          = squeeze(handles.vField(:,:,:,fieldSlice));
      end
      if (get(handles.hFieldToggle,'Value'))
        s(3,:,:) = s(3,:,:) + handles.slices(3);
        s(1:2,:,:) = s(1:2,:,:) + eyeHField([size(s,2) size(s,3)]);
      end
    case 'sagittal'
      s = zeros(3,handles.dim(2),handles.dim(3));
      fieldSlice = handles.slices(1) - handles.fieldStart(1);
      if fieldSlice > 0 & fieldSlice <= handles.fieldSize(2)
        s(:,...
          handles.fieldStart(2)+1:handles.fieldStart(2)+handles.fieldSize(3),...
          handles.fieldStart(3)+1:handles.fieldStart(3)+handles.fieldSize(4))...
          = squeeze(handles.vField(:,fieldSlice,:,:));
      end
      if (get(handles.hFieldToggle,'Value'))
        s(1,:,:) = s(1,:,:) + handles.slices(1);
        s(2:3,:,:) = s(2:3,:,:) + eyeHField([size(s,2) size(s,3)]);
      end
    case 'coronal'
      s = zeros(3,handles.dim(1),handles.dim(3));
      fieldSlice = handles.slices(2) - handles.fieldStart(2);
      if fieldSlice > 0 & fieldSlice <= handles.fieldSize(3)
        s(:,...
          handles.fieldStart(1)+1:handles.fieldStart(1)+handles.fieldSize(2),...
          handles.fieldStart(3)+1:handles.fieldStart(3)+handles.fieldSize(4))...
          = squeeze(handles.vField(:,:,fieldSlice,:));
      end
      if (get(handles.hFieldToggle,'Value'))
        s(2,:,:) = s(2,:,:) + handles.slice(2);
        s([1 3],:,:) = s([1 3],:,:) + eyeHField([size(s,2) size(s,3)]);
      end
  end
end




function s = getImageSlice(handles)
switch handles.orientationStrs{handles.orientation}
 case 'axial'
  s = squeeze(handles.im(:,:,handles.slices(3)));
 case 'sagittal'
  s = squeeze(handles.im(handles.slices(1),:,:));
 case 'coronal'
  s = squeeze(handles.im(:,handles.slices(2),:));
end

function I = createDefaultImage(dim)
[Y,X,Z]=meshgrid(1:dim(2),1:dim(1),1:dim(3));
I = (X*(1/2)+Y*(1/3)+Z*(1/6)-1);
I = I * 255/max(I(:));

% --- Executes during object creation, after setting all properties.
function imageMaxSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageMaxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function imageLevelSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageLevelSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function imageMinSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageMinSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function fieldMaxSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldMaxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function fieldLevelSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldLevelSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function fieldMinSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldMinSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function imageWindowTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageWindowTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function imageWindowPresets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageWindowPresets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', {'Max Range','Custom','Bone','Prostate','Prostate2','Prostate3','Lung','Gas'});
set(hObject,'Value',1);

% --- Executes during object creation, after setting all properties.
function imageFieldBlendSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageFieldBlendSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function fieldWindowTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldWindowTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function sliceSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function sliceNumberTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceNumberTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axialToggle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axialToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function coronalToggle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coronalToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function sagittalToggle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sagittalToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function hFieldToggle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hFieldToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function vFieldToggle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vFieldToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function view2DToggle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view2DToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function view3DToggle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view3DToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function snapPhotoButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snapPhotoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function centerImageButton_CreateFcn(hObject, eventdata, handles)
%

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function fieldWindowPresets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldWindowPresets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', {'Max Range','Custom','Zero-Max'});
set(hObject,'Value',1);

% --- Executes during object creation, after setting all properties.
function imageColormapPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageColormapPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', {'autumn','bone','colorcube','cool','copper','flag','gray','hot','hsv','jet','lines','pink','prism','spring','summer','white','winter'});
set(hObject,'Value',2);

% --- Executes on selection change in imageColormapPopup.
function imageColormapPopup_Callback(hObject, eventdata, handles)
% hObject    handle to imageColormapPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns imageColormapPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imageColormapPopup
updateDisplay(handles);

% --- Executes during object creation, after setting all properties.
function fieldColormapPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldColormapPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', {'autumn','bone','colorcube','cool','copper','flag','gray','hot','hsv','jet','lines','pink','prism','spring','summer','white','winter'});
set(hObject,'Value',10);

% --- Executes on selection change in fieldColormapPopup.
function fieldColormapPopup_Callback(hObject, eventdata, handles)
% hObject    handle to fieldColormapPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fieldColormapPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fieldColormapPopup
updateDisplay(handles);


% --- Executes on slider movement.
function imageFieldBlendSlider_Callback(hObject, eventdata, handles)
% hObject    handle to imageFieldBlendSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
updateDisplay(handles);

function sliceSlider_Callback(hObject, eventdata, handles)
% round to nearest slice
val = round(get(hObject,'Value'));
set(hObject,'Value',val);
handles.slices(handles.orientation) = val;
guidata(hObject,handles);
updateDisplay(handles);
updateSliceTag(handles);

% --- Executes during object creation, after setting all properties.
function fieldViewPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldViewPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', {'L_2 Norm','L_1 Norm','L_inf Norm','det(Jacobian)','V_1','V_2','V_3','Quiver 2D'});
set(hObject,'Value',1);

% --- Executes on selection change in fieldViewPopup.
function fieldViewPopup_Callback(hObject, eventdata, handles)
% hObject    handle to fieldViewPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fieldViewPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fieldViewPopup
updateDisplay(handles);

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)



% --------------------------------------------------------------------
function loadImageMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,pathname] = uigetfile('*.mhd');
if ~isequal(file, 0)
    % load the image
    [handles.im,handles.imOrigin,handles.imSpacing]= loadMETA([pathname file]);

    % unload field if this is a new size
    if any(handles.dim~=size(handles.im))
      handles.vField = [];
      handles.detJ = [];
      handles.haveField = false;
      handles.dim = size(handles.im);
    end

    handles.imageMin=min(handles.im(:));
    handles.imageMax=max(handles.im(:));
    handles.imageWindowMin = handles.imageMin;
    handles.imageWindowMax = handles.imageMax;
    set(handles.imageWindowPresets,'Value',1);
    updateImageWindowSlider(handles);

    % update slice slider settings
    handles.maxSlice = handles.dim(handles.orientation);
    idx = handles.slices > handles.dim;
    handles.slices(idx) = handles.dim(idx)/2;

    % update application data
    guidata(hObject,handles);
    
    % update the image
    updateSliceSlider(handles);
    updateFieldWindowSlider(handles);
    updateDisplay(handles);
end

% --------------------------------------------------------------------
function loadDefaultImageMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.im = createDefaultImage(handles.dim);
handles.imageMin=min(handles.im(:));
handles.imageMax=max(handles.im(:));
handles.imageWindowMin = handles.imageMin;
handles.imageWindowMax = handles.imageMax;
set(handles.imageWindowPresets,'Value',1);
updateImageWindowSlider(handles);
guidata(hObject,handles);
updateSliceSlider(handles);
updateDisplay(handles);

% --------------------------------------------------------------------
function loadBlankImageMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.im(:) = 0;
handles.imageMin=0;
handles.imageMax=0;
handles.imageWindowMin = 0;
handles.imageWindowMax = 0;
set(handles.imageWindowPresets,'Value',1);
updateImageWindowSlider(handles);
guidata(hObject,handles);
updateSliceSlider(handles);
updateDisplay(handles);

% --------------------------------------------------------------------
function loadOneImageMenu_Callback(hObject, eventdata, handles)
% hObject    handle to loadOneImageMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.im(:) = 1;
handles.imageMin=1;
handles.imageMax=1;
handles.imageWindowMin = 1;
handles.imageWindowMax = 1;
set(handles.imageWindowPresets,'Value',1);
updateImageWindowSlider(handles);
guidata(hObject,handles);
updateSliceSlider(handles);
updateDisplay(handles);



% --- Executes on button press in axialToggle.
function axialToggle_Callback(hObject, eventdata, handles)
% hObject    handle to axialToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axialToggle
if get(hObject,'Value')
    % update look of buttons
    set(hObject,'Background',[1 0 0]);
    set(handles.coronalToggle,'Value',0);
    set(handles.sagittalToggle,'Value',0);
    set(handles.coronalToggle,'Background',[0.65 0.65 0]);
    set(handles.sagittalToggle,'Background',[0.65 0.65 0]);    

    % update saved state
    handles.orientation = 3;
    handles.maxSlice = handles.dim(3);
    guidata(hObject,handles);
    
    % update slider and image
    updateSliceSlider(handles);
    updateDisplay(handles);
else
    set(hObject,'Value',1);
end

% --- Executes on button press in sagittalToggle.
function sagittalToggle_Callback(hObject, eventdata, handles)
% hObject    handle to sagittalToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sagittalToggle
if get(hObject,'Value')
    % update look of buttons
    set(hObject,'Background',[1 0 0]);
    set(handles.coronalToggle,'Value',0);
    set(handles.axialToggle,'Value',0);
    set(handles.coronalToggle,'Background',[0.65 0.65 0]);
    set(handles.axialToggle,'Background',[0.65 0.65 0]);    

    % update saved state
    handles.orientation = 1;
    handles.maxSlice = handles.dim(1);
    guidata(hObject,handles);
    
    % update slider and image
    updateSliceSlider(handles);
    updateDisplay(handles);
else
    set(hObject,'Value',1);
end

% --- Executes on button press in coronalToggle.
function coronalToggle_Callback(hObject, eventdata, handles)
% hObject    handle to coronalToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coronalToggle
if get(hObject,'Value')
    % update look of buttons
    set(hObject,'Background',[1 0 0]);
    set(handles.axialToggle,'Value',0);
    set(handles.sagittalToggle,'Value',0);
    set(handles.axialToggle,'Background',[0.65 0.65 0]);
    set(handles.sagittalToggle,'Background',[0.65 0.65 0]);    

    % update saved state
    handles.orientation = 2;
    handles.maxSlice = handles.dim(2);
    guidata(hObject,handles);
    
    % update slider and image
    updateSliceSlider(handles);
    updateDisplay(handles);
else
    set(hObject,'Value',1);
end

function view2DToggle_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    % update look of buttons
    set(hObject,'Background',[1 0 0]);
    set(handles.view3DToggle,'Value',0);
    set(handles.view3DToggle,'Background',[0.65 0.65 0]);
    updateDisplay(handles);
else
    set(hObject,'Value',1);
end

function view3DToggle_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    % update look of buttons
    set(hObject,'Background',[1 0 0]);
    set(handles.view2DToggle,'Value',0);
    set(handles.view2DToggle,'Background',[0.65 0.65 0]);
    updateDisplay(handles);
else
    set(hObject,'Value',1);
end

function hFieldToggle_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    % update look of buttons
    set(hObject,'Background',[1 0 0]);
    set(handles.vFieldToggle,'Value',0);
    set(handles.vFieldToggle,'Background',[0.65 0.65 0]);
    updateDisplay(handles);
else
    set(hObject,'Value',1);
end

function vFieldToggle_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    % update look of buttons
    set(hObject,'Background',[1 0 0]);
    set(handles.hFieldToggle,'Value',0);
    set(handles.hFieldToggle,'Background',[0.65 0.65 0]);
    updateDisplay(handles);
else
    set(hObject,'Value',1);
end

function updateSliceTag(handles)
slice = handles.slices(handles.orientation);
tagStr = sprintf('%03d\n%03d',slice,handles.maxSlice);
set(handles.sliceNumberTag,'String',tagStr);

function updateImageWindowSlider(handles)
set(handles.imageMinSlider,'Min',handles.imageMin);
set(handles.imageMinSlider,'Max',handles.imageMax+eps);
set(handles.imageMinSlider,'Value',handles.imageWindowMin);
set(handles.imageMaxSlider,'Min',handles.imageMin);
set(handles.imageMaxSlider,'Max',handles.imageMax+eps);
set(handles.imageMaxSlider,'Value',handles.imageWindowMax);

set(handles.imageLevelSlider,'Min',handles.imageMin);
set(handles.imageLevelSlider,'Max',handles.imageMax+eps);
set(handles.imageLevelSlider,'Value',(handles.imageWindowMax+handles.imageWindowMin)/2);
tagStr = sprintf('%6.2f:%6.2f\n%6.2f:%6.2f',handles.imageWindowMin,handles.imageWindowMax,handles.imageMin,handles.imageMax);
set(handles.imageWindowTag,'String',tagStr);

function updateFieldWindowSlider(handles)
set(handles.fieldMinSlider,'Min',handles.fieldMin);
set(handles.fieldMinSlider,'Max',handles.fieldMax+eps);
set(handles.fieldMinSlider,'Value',handles.fieldWindowMin);
set(handles.fieldMaxSlider,'Min',handles.fieldMin);
set(handles.fieldMaxSlider,'Max',handles.fieldMax+eps);
set(handles.fieldMaxSlider,'Value',handles.fieldWindowMax);

set(handles.fieldLevelSlider,'Min',handles.fieldMin);
set(handles.fieldLevelSlider,'Max',handles.fieldMax+eps);
set(handles.fieldLevelSlider,'Value',(handles.fieldWindowMax+handles.fieldWindowMin)/2);
tagStr = sprintf('%6.2f:%6.2f\n%6.2f:%6.2f',handles.fieldWindowMin,handles.fieldWindowMax,handles.fieldMin,handles.fieldMax);
set(handles.fieldWindowTag,'String',tagStr);

function updateSliceSlider(handles)
set(handles.sliceSlider,'Min',handles.minSlice);
set(handles.sliceSlider,'Max',handles.maxSlice);
slice = handles.slices(handles.orientation);
set(handles.sliceSlider,'Value',slice);
set(handles.sliceSlider,'SliderStep',...
    [1/(handles.maxSlice-handles.minSlice), ...
     10/(handles.maxSlice-handles.minSlice)]);
updateSliceTag(handles);

% --- Executes on slider movement.
function imageMaxSlider_Callback(hObject, eventdata, handles)
% hObject    handle to imageMaxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
minVal = get(handles.imageMinSlider,'Value');
if val < minVal
    set(hObject,'Value',minVal+eps);
    val = minVal+eps;
end
handles.imageWindowMax = val;
guidata(hObject,handles);
updateImageWindowSlider(handles);
set(handles.imageWindowPresets,'Value',2);
updateDisplay(handles);

% --- Executes on slider movement.
function imageLevelSlider_Callback(hObject, eventdata, handles)
% hObject    handle to imageLevelSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
wdif = (handles.imageWindowMax-handles.imageWindowMin)/2;

handles.imageWindowMax = val+wdif;
handles.imageWindowMin = val-wdif;

if handles.imageWindowMax > handles.imageMax
    handles.imageWindowMax = handles.imageMax;
    handles.imageWindowMin = handles.imageMax - 2*wdif;
elseif handles.imageWindowMin < handles.imageMin
    handles.imageWindowMin = handles.imageMin;
    handles.imageWindowMax = handles.imageMin + 2*wdif;
end
updateImageWindowSlider(handles);
guidata(hObject,handles);
set(handles.imageWindowPresets,'Value',2);
updateDisplay(handles);


% --- Executes on slider movement.
function imageMinSlider_Callback(hObject, eventdata, handles)
% hObject    handle to imageMinSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
maxVal = get(handles.imageMaxSlider,'Value');
if val > maxVal
    set(hObject,'Value',maxVal-eps);
    val = maxVal-eps;
end
handles.imageWindowMin = val;
guidata(hObject,handles);
updateImageWindowSlider(handles);
set(handles.imageWindowPresets,'Value',2);
updateDisplay(handles);


% --- Executes on selection change in imageWindowPresets.
function imageWindowPresets_Callback(hObject, eventdata, handles)
% hObject    handle to imageWindowPresets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns imageWindowPresets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imageWindowPresets
strs = get(hObject,'String');
str  = strs{get(hObject,'Value')};
switch str
 case 'Max Range'
  handles.imageWindowMin = handles.imageMin;
  handles.imageWindowMax = handles.imageMax;
 case 'Bone'
  handles.imageWindowMin = 1150;
  handles.imageWindowMax = 1360;  
 case 'Prostate'
  handles.imageWindowMin = 900;
  handles.imageWindowMax = 1100;      
 case 'Prostate2'
  handles.imageWindowMin = 800;
  handles.imageWindowMax = 1200;
 case 'Prostate3'
  handles.imageWindowMin = 650;
  handles.imageWindowMax = 1250;
 case 'Lung'
  handles.imageWindowMin = 0;
  handles.imageWindowMax = 1456;
 case 'Gas'
  handles.imageWindowMin = 750;
  handles.imageWindowMax = 751;
end
guidata(hObject,handles);
updateImageWindowSlider(handles);
updateDisplay(handles);

% --------------------------------------------------------------------
function unloadFieldMenu_Callback(hObject, eventdata, handles)
% hObject    handle to loadOneImageMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.vField = [];
handles.detJ = [];
handles.fieldSize = [];
handles.fieldStart = [];
handles.haveField = false;
updateFieldWindowSlider(handles);
guidata(hObject,handles);
updateDisplay(handles);

% --------------------------------------------------------------------
function loadVFieldMenu_Callback(hObject, eventdata, handles)
% hObject    handle to loadOneImageMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, pathname] = uigetfile('*.mhd');
if ~isequal(file, 0)
    % load the field
    handles.vField = loadMETA([pathname file]);
    handles.haveField = true;
    handles.detJ = jacobianDetVField(handles.vField);    
    handles.fieldStart = [0 0 0];
    handles.fieldSize  = size(handles.vField);

    % load the default image if this is a new dimension
    fieldDim = size(handles.vField);
    if any(handles.dim~=fieldDim(2:end))
      handles.dim = fieldDim(2:end);
      handles.im = createDefaultImage(handles.dim);
      handles.imageMin=min(handles.im(:));
      handles.imageMax=max(handles.im(:));
      handles.imageWindowMin = handles.imageMin;
      handles.imageWindowMax = handles.imageMax;
      set(handles.imageWindowPresets,'Value',1);
      updateImageWindowSlider(handles);      
    end

    L1 = sum(abs(handles.vField));
    maxL1 = max(L1(:));
    
    handles.fieldMin=-maxL1;
    handles.fieldMax=maxL1;
    handles.fieldWindowMin = handles.fieldMin;
    handles.fieldWindowMax = handles.fieldMax;
    set(handles.fieldWindowPresets,'Value',1);
    updateFieldWindowSlider(handles);
    
    % update slice slider settings
    handles.maxSlice = handles.dim(handles.orientation);
    idx = handles.slices > handles.dim;
    handles.slices(idx) = handles.dim(idx)/2;

    % update application data
    guidata(hObject,handles);

    % update the image
    updateSliceSlider(handles);
    updateDisplay(handles);
end


% --------------------------------------------------------------------
function loadHFieldMenu_Callback(hObject, eventdata, handles)
% hObject    handle to loadOneImageMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, pathname] = uigetfile('*.mhd');
if ~isequal(file, 0)
    % load the image
    field = loadMETA([pathname file])+1;
    hSize = size(field);
    eyeH = eyeHField(hSize(2:end));
    handles.vField = field-eyeH;
    handles.haveField = true;
    handles.detJ = jacobianDetVField(handles.vField);
    handles.fieldStart = [0 0 0];
    handles.fieldSize  = size(handles.vField);

    % load the default image if this is a new dimension
    fieldDim = size(handles.vField);
    if any(handles.dim~=fieldDim(2:end))
      handles.dim = fieldDim(2:end);
      handles.im = createDefaultImage(handles.dim);
      handles.imageMin=min(handles.im(:));
      handles.imageMax=max(handles.im(:));
      handles.imageWindowMin = handles.imageMin;
      handles.imageWindowMax = handles.imageMax;
      set(handles.imageWindowPresets,'Value',1);
      updateImageWindowSlider(handles);      
    end
    
    L1 = sum(abs(handles.vField));
    maxL1 = max(L1(:));
    
    handles.fieldMin=-maxL1;
    handles.fieldMax=maxL1;
    handles.fieldWindowMin = handles.fieldMin;
    handles.fieldWindowMax = handles.fieldMax;
    set(handles.fieldWindowPresets,'Value',1);
    updateFieldWindowSlider(handles);

    % update slice slider settings
    handles.maxSlice = handles.dim(handles.orientation);
    idx = handles.slices > handles.dim;
    handles.slices(idx) = handles.dim(idx)/2;
    
    % update application data
    guidata(hObject,handles);

    % update the image
    updateSliceSlider(handles);
    updateDisplay(handles);
end

% --------------------------------------------------------------------
function loadHFieldROIMenu_Callback(hObject, eventdata, handles)
[file, pathname] = uigetfile('*.mhd');
if ~isequal(file, 0)
    % load the field
    [f,fOrigin,fSpacing] = loadMETA([pathname file]);
    f = f+1;
    handles.fieldSize = size(f);
    handles.fieldStart = round((fOrigin-handles.imOrigin)./handles.imSpacing);
    fprintf('roiStart: %d\n',handles.fieldStart);

    % set as current field
    eyeH = eyeHField(handles.fieldSize(2:end));
    handles.vField = f-eyeH;
    handles.vField(isnan(handles.vField))=0;
    handles.haveField = true;
    handles.detJ = jacobianDetVField(handles.vField);

    L1 = sum(abs(handles.vField));
    maxL1 = max(L1(:));
    
    handles.fieldMin=-maxL1;
    handles.fieldMax=maxL1;
    handles.fieldWindowMin = handles.fieldMin;
    handles.fieldWindowMax = handles.fieldMax;
    set(handles.fieldWindowPresets,'Value',1);
    updateFieldWindowSlider(handles);

    % update application data
    guidata(hObject,handles);

    % update the image
    updateDisplay(handles);
end

function loadVFieldROIMenu_Callback(hObject, eventdata, handles)
[file, pathname] = uigetfile('*.mhd');
if ~isequal(file, 0)
    % load the field
    [f,fOrigin,fSpacing] = loadMETA([pathname file]);
    handles.fieldSize = size(f);

    % find index of field into image
    handles.fieldStart = round((fOrigin-handles.imOrigin)./handles.imSpacing);
    fprintf('roiStart: %d\n',handles.fieldStart);

    % set as current field
    handles.vField = f;
    handles.vField(isnan(handles.vField))=0;
    handles.haveField = true;
    handles.detJ = jacobianDetVField(handles.vField);

    L1 = sum(abs(handles.vField));
    maxL1 = max(L1(:));
    
    handles.fieldMin=-maxL1;
    handles.fieldMax=maxL1;
    handles.fieldWindowMin = handles.fieldMin;
    handles.fieldWindowMax = handles.fieldMax;
    set(handles.fieldWindowPresets,'Value',1);
    updateFieldWindowSlider(handles);

    % update application data
    guidata(hObject,handles);

    % update the image
    updateDisplay(handles);
end

function loadIDField_Callback(hObject, eventdata, handles)
handles.vField = eyeHField(handles.dim);%zeros([3 handles.dim]);
handles.haveField = true;
handles.detJ = [];%ones(size(handles.im));
handles.fieldStart = [0 0 0];
handles.fieldSize  = size(handles.vField);

L1 = sum(abs(handles.vField));
maxL1 = max(L1(:));
%maxL1 = 0;

handles.fieldMin=-maxL1;
handles.fieldMax=maxL1;
handles.fieldWindowMin = handles.fieldMin;
handles.fieldWindowMax = handles.fieldMax;
set(handles.fieldWindowPresets,'Value',1);
updateFieldWindowSlider(handles);

% update application data
guidata(hObject,handles);

% update the image
updateDisplay(handles);

function loadImageGradientField_Callback(hObject, eventdata, handles)
handles.vField = imageGradient(handles.im);
handles.haveField = true;
handles.detJ = jacobianDetVField(handles.vField);
handles.fieldStart = [0 0 0];
handles.fieldSize  = size(handles.vField);

L1 = sum(abs(handles.vField));
maxL1 = max(L1(:));

handles.fieldMin=-maxL1;
handles.fieldMax=maxL1;
handles.fieldWindowMin = handles.fieldMin;
handles.fieldWindowMax = handles.fieldMax;
set(handles.fieldWindowPresets,'Value',1);
updateFieldWindowSlider(handles);

% update application data
guidata(hObject,handles);

% update the image
updateDisplay(handles);

% --- Executes on selection change in fieldWindowPresets.
function fieldWindowPresets_Callback(hObject, eventdata, handles)
% hObject    handle to fieldWindowPresets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fieldWindowPresets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fieldWindowPresets
strs = get(hObject,'String');
str  = strs{get(hObject,'Value')};
switch str
 case 'Max Range'
  handles.fieldWindowMin = handles.fieldMin;
  handles.fieldWindowMax = handles.fieldMax;
 case 'Zero-Max'
  handles.fieldWindowMin = 0;
  handles.fieldWindowMax = handles.fieldMax;  
end
guidata(hObject,handles);
updateFieldWindowSlider(handles);
updateDisplay(handles);


% --- Executes on slider movement.
function fieldMaxSlider_Callback(hObject, eventdata, handles)
% hObject    handle to fieldMaxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
minVal = get(handles.fieldMinSlider,'Value');
if val < minVal
    set(hObject,'Value',minVal+eps);
    val = minVal+eps;
end
handles.fieldWindowMax = val;
guidata(hObject,handles);
updateFieldWindowSlider(handles);
set(handles.fieldWindowPresets,'Value',2);
updateDisplay(handles);

% --- Executes on slider movement.
function fieldLevelSlider_Callback(hObject, eventdata, handles)
% hObject    handle to fieldLevelSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
wdif = (handles.fieldWindowMax-handles.fieldWindowMin)/2;

handles.fieldWindowMax = val+wdif;
handles.fieldWindowMin = val-wdif;

if handles.fieldWindowMax > handles.fieldMax
    handles.fieldWindowMax = handles.fieldMax;
    handles.fieldWindowMin = handles.fieldMax - 2*wdif;
elseif handles.fieldWindowMin < handles.fieldMin
    handles.fieldWindowMin = handles.fieldMin;
    handles.fieldWindowMax = handles.fieldMin + 2*wdif;
end
updateFieldWindowSlider(handles);
guidata(hObject,handles);
set(handles.fieldWindowPresets,'Value',2);
updateDisplay(handles);

% --- Executes on slider movement.
function fieldMinSlider_Callback(hObject, eventdata, handles)
% hObject    handle to fieldMinSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
maxVal = get(handles.fieldMaxSlider,'Value');
if val > maxVal
    set(hObject,'Value',maxVal-eps);
    val = maxVal-eps;
end
handles.fieldWindowMin = val;
guidata(hObject,handles);
updateFieldWindowSlider(handles);
set(handles.fieldWindowPresets,'Value',2);
updateDisplay(handles);

function snapPhotoButton_Callback(hObject, eventdata, handles)
% hObject    handle to snapPhotoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path,fidx] = uiputfile('*.png','Save image as...');
if ~isequal(file, 0)
  fprintf('Saving image: %s%s\n',path,file);
  % save current axes as a png image.
  print('-dpng',[path file]);
end

% --- Executes during object creation, after setting all properties.
function quiverSpacingPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quiverSpacingPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', {'1','2','3','4','5','6','7','8','9','10','15','20'});
set(hObject,'Value',5);

% --- Executes on selection change in quiverSpacingPopup.
function quiverSpacingPopup_Callback(hObject, eventdata, handles)
% hObject    handle to quiverSpacingPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns quiverSpacingPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from quiverSpacingPopup
updateDisplay(handles);

% --- Executes during object creation, after setting all properties.
function quiverColorPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quiverColorPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', {'blue','green','red','cyan','magenta','yellow','black'});
set(hObject,'Value',1);

% --- Executes on selection change in quiverColorPopup.
function quiverColorPopup_Callback(hObject, eventdata, handles)
% hObject    handle to quiverColorPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns quiverColorPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from quiverColorPopup
updateDisplay(handles);



function c = getQuiverSpacing(handles)
strs = get(handles.quiverSpacingPopup,'String');
c = str2num(strs{get(handles.quiverSpacingPopup,'Value')});

function c = getQuiverScaling(handles)
c = get(handles.quiverScalingSlider,'Value');

function c = getQuiverColor(handles)
strs = get(handles.quiverColorPopup,'String');
color = strs{get(handles.quiverColorPopup,'Value')};
switch color
    case 'blue'
        c = 'b';    
    case 'green'
        c = 'g';
    case 'red'
        c = 'r';
    case 'cyan' 
        c = 'c';
    case 'magenta'
        c = 'm'
    case 'yellow'
        c = 'y';
    case 'black'
        c = 'k';
end


% --- Executes during object creation, after setting all properties.
function quiverScalingSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quiverScalingSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject,'Min',0.01);
set(hObject,'Max',5);
set(hObject,'Value',1);

% --- Executes on slider movement.
function quiverScalingSlider_Callback(hObject, eventdata, handles)
% hObject    handle to quiverScalingSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.quiverScalingTag,'String',sprintf('Quiver Scaling: %0.2g',get(hObject,'Value')));
updateDisplay(handles);


% --- Executes on button press in centerImageButton.
function centerImageButton_Callback(hObject, eventdata, handles)
% hObject    handle to centerImageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.slices = handles.dim/2;
guidata(hObject, handles);
updateSliceSlider(handles);
updateDisplay(handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


