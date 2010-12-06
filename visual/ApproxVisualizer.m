function varargout = ApproxVisualizer(varargin)
% APPROXVISUALIZER M-file for ApproxVisualizer.fig
%      APPROXVISUALIZER, by itself, creates a new APPROXVISUALIZER or raises the existing
%      singleton*.
%
%      H = APPROXVISUALIZER returns the handle to a new APPROXVISUALIZER or the handle to
%      the existing singleton*.
%
%      APPROXVISUALIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APPROXVISUALIZER.M with the given input arguments.
%
%      APPROXVISUALIZER('Property','Value',...) creates a new APPROXVISUALIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ApproxVisualizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ApproxVisualizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ApproxVisualizer

% Last Modified by GUIDE v2.5 12-Oct-2010 14:44:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ApproxVisualizer_OpeningFcn, ...
    'gui_OutputFcn',  @ApproxVisualizer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ApproxVisualizer is made visible.
function ApproxVisualizer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ApproxVisualizer (see VARARGIN)

% Choose default command line output for ApproxVisualizer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ApproxVisualizer wait for user response (see UIRESUME)
% uiwait(handles.main);

%%%%%%%%%%%% Custom code %%%%%%%%%%%%%
if isempty(varargin)
    disp('Pass a reduced model to ApproxVisualizer.');
    return;
end
if ~isa(varargin{1},'models.ReducedModel')
    error('ApproxVisualizer takes only models.ReducedModel subclasses as input.');
end

% Store data
r = varargin{1};
m = r.FullModel;
setappdata(handles.main,'r',r);
setappdata(handles.main,'atd',m.Data.ApproxTrainData);

% Initialize GUI
modelToGUI(r, handles);

% Set initial values
if isempty(r.ParamSamples)
    setappdata(handles.main,'pidx',[]);
else
    setappdata(handles.main,'pidx',1);
end
if m.System.InputCount == 0
    setappdata(handles.main,'inidx',[]);
else
    setappdata(handles.main,'inidx',1);
end
setappdata(handles.main,'fidx',1);
setappdata(handles.main,'factor',1);
setappdata(handles.main,'mu',[]);

% Plot default visualization
plotCurrent(handles);


function modelToGUI(r, h)
m = r.FullModel;
% Disable controls specific to model settings
if m.System.ParamCount == 0
    set(h.pslide,'Enable','off');
else
    pc = size(r.ParamSamples,2);
    set(h.pslide,'Max',pc,'Min',1,'Value',1);
    set(h.pslide,'SliderStep',[1/(pc-1) 10/(pc-1)]);
end
if m.System.InputCount == 0
    set(h.islide,'Enable','off');
else
    set(h.islide,'Max',m.System.InputCount,'Min',1,'Value',1);
    set(h.islide,'SliderStep',[1/m.System.InputCount 10/m.System.InputCount]);
end
fdims = size(m.Data.ApproxTrainData,1)-3;
set(h.fslide,'Max',fdims,'Value',1);
set(h.fslide,'SliderStep',[1/(fdims-1) 10/(fdims-1)]);
set(h.txtModel,'String',r.Name);

set(h.pnlParams,'Units','pixels');
pos = get(h.pnlParams,'Position');
set(h.pnlParams,'Units','normalized');
dist = 22;%px
pcnt = 0;
% Create model parameter slides
for pidx = 1:m.ParamCount
    p = m.Params(pidx);
    if p.HasRange
        % Top location
        top = pos(4)-pcnt*dist-40;
        % Create labels
        label = uicontrol('Tag',['runtime_lbl' num2str(pidx)],'Style','text',...
                'Parent',h.pnlParams,'HorizontalAlignment','left');
        set(label,'String',p.Name,'Units','pixels', 'Position',[10 top 50 14]);
        set(label,'Units','normalized');
        label = uicontrol('Tag',['runtime_lbll' num2str(pidx)],'Style','text',...
                'Parent',h.pnlParams,'HorizontalAlignment','right');
        set(label,'String',sprintf('%2.3f',p.MinVal),'Units','pixels',...
            'Position',[80 top 30 14]);
        set(label,'Units','normalized');
        label = uicontrol('Tag',['runtime_lblr' num2str(pidx)],'Style','text',...
                'Parent',h.pnlParams,'HorizontalAlignment','left');
        set(label,'String',sprintf('%2.3f',p.MaxVal),'Units','pixels',...
            'Position',[260 top 30 14]);
        set(label,'Units','normalized');
        % Create slider
        tag = ['runtime_pslide' num2str(pidx)];
        ctrl = uicontrol('Tag',tag,'Parent',...
            h.pnlParams,'Style','slider','UserData',pidx);
        h.(tag) = ctrl;
        % Position
        set(ctrl,'Units','pixels','Position',[120 top 130 16]);
        set(ctrl,'Units','normalized');
        % Range etc
        set(ctrl,'Min',p.MinVal,'Max',p.MaxVal,'Value',p.MinVal);
        set(ctrl,'SliderStep',[0.01 0.1]);
        % Set callback & string
        set(ctrl,'Callback',@(hO,e)(updateUserParam(guidata(hO))));
       
        % increase position counter
        pcnt = pcnt + 1;
    end
end

function updateUserParam(h)
r = getappdata(h.main,'r');
pc = r.FullModel.ParamCount;
mu = zeros(pc,1);
slides = findobj(h.pnlParams,'Style','slider');
for n=1:pc
    if r.FullModel.Params(n).HasRange
        mu(n) = get(slides(pc-n+1),'Value');
    else
        mu(n) = r.FullModel.Params(n).MinVal;
    end
end
setappdata(h.main,'mu',mu);
set(h.txtUParam,'String',['[' sprintf('%2.3f ',mu) ']']);
plotCurrent(h);

function plotCurrent(h)

% Get approximation training data (atd)
r = getappdata(h.main,'r');
inidx = getappdata(h.main,'inidx');

if get(h.rbS,'Value') == 1
    d = r.FullModel.Data;
    sn = d.ApproxTrainData;
    % Select required subset of training data
    pidx = getappdata(h.main,'pidx');
    
    if ~isempty(pidx)
        psel = sn(1,:) == pidx;
    else
        psel = true(1,size(sn,2));
    end
    if ~isempty(inidx)
        isel = sn(2,:) == inidx;
    else
        isel = true(1,size(sn,2));
    end
    sel = logical(psel .* isel);
    sn = sn(:,sel);
    xi = sn(4:end,:);
    ti = sn(3,:);

    % Compute relation to non-extended atd base
    factor = getappdata(h.main,'factor');
    nz = factor*size(sn,2);
    oldidx = 1:factor:nz;
    XI(:,oldidx) = xi;
    TI(oldidx) = ti;
    mui(oldidx) = sn(1,:);
    fx(:,oldidx) = d.ApproxfValues(:,sel);
    if factor > 1
        di = (xi(:,2:end)-xi(:,1:end-1))/factor;
        di(:,end+1) = 0;
        dti = (ti(2:end)-ti(1:end-1))/factor;
        dti(end+1) = 0;

        for fac = 1:factor-1
            newidx = oldidx+fac;
            % xi
            XI(:,newidx) = xi+di*fac;
            % ti
            TI(newidx) = ti+dti*fac;
            % mui
            mui(newidx) = sn(1,:);

            % Evaluate original function at middle points
            for idx=newidx
                fx(:,idx) = r.FullModel.System.f.evaluate(XI(:,idx),TI(idx),d.getParams(mui(idx)));
            end
        end
    end
    MUI = d.getParams(mui);    
else
    mu = getappdata(h.main,'mu');
    [TI,XI] = r.FullModel.computeTrajectory(mu,inidx);
    n = length(TI);
    MUI = repmat(mu,1,n);
    
    fx = zeros(size(XI,1),n);
    for idx=1:n
        fx(:,idx) = r.FullModel.System.f.evaluate(XI(:,idx),TI(idx),MUI(:,idx));
    end
end

afx = r.FullModel.Approx.evaluate(XI,TI,MUI);
setappdata(h.main,'fx',fx);
setappdata(h.main,'afx',afx);
setappdata(h.main,'ti',TI);

plotCurrentFxi(h)

function plotCurrentFxi(h)
% Get current f
fidx = getappdata(h.main,'fidx');
fx = getappdata(h.main,'fx');
afx = getappdata(h.main,'afx');
ti = getappdata(h.main,'ti');
fxi = fx(fidx,:);
afxi = afx(fidx,:);

err = sqrt(sum((fxi-afxi).^2));

plot(ti,fxi,'r',ti,afxi,'b');

% Plot extra markers for sampled data
if get(h.rbS,'Value') == 1
    factor = getappdata(h.main,'factor');
    nz = size(ti,2);
    oldidx = 1:factor:nz;
    newidx = setdiff(1:nz,oldidx);
    hold on;
    plot(ti(oldidx),fxi(oldidx),'rs',ti(newidx),fxi(newidx),'blackx');
    hold off;
end
axis tight;
title(sprintf('f,f'' at dimension %d. L2-Error:%f',fidx,err));
xlabel('t');
ylabel('f(x)');

currentToGUI(getappdata(h.main,'r'),h);

function currentToGUI(r,h)
pidx = getappdata(h.main,'pidx');
if ~isempty(pidx)
    %set(h.pslide,'Value',pidx);
    set(h.txtp,'String',['[' sprintf('%2.3f ',r.ParamSamples(:,pidx)) ']']);
end
inidx = getappdata(h.main,'inidx');
if ~isempty(inidx)
    %set(h.islide,'Value',inidx);
    set(h.txti,'String',sprintf('%d',inidx));
end
set(h.txtf,'String',sprintf('%d',getappdata(h.main,'fidx')));
set(h.txtg,'String',sprintf('%d',getappdata(h.main,'factor')));


% --- Outputs from this function are returned to the command line.
function varargout = ApproxVisualizer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function pslide_Callback(hObject, eventdata, handles)
% hObject    handle to pslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider
setappdata(handles.main,'pidx',round(get(hObject,'Value')));
plotCurrent(handles);

% --- Executes during object creation, after setting all properties.
function pslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function islide_Callback(hObject, eventdata, handles)
% hObject    handle to islide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
setappdata(handles.main,'inidx',round(get(hObject,'Value')));
plotCurrent(handles);

% --- Executes during object creation, after setting all properties.
function islide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to islide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function fslide_Callback(hObject, eventdata, handles)
% hObject    handle to fslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
setappdata(handles.main,'fidx',round(get(hObject,'Value')));
plotCurrentFxi(handles);

% --- Executes during object creation, after setting all properties.
function fslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function factorslide_Callback(hObject, eventdata, handles)
% hObject    handle to factorslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
setappdata(handles.main,'factor',round(get(hObject,'Value')));
plotCurrent(handles);


% --- Executes during object creation, after setting all properties.
function factorslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to factorslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in rbS.
function rbS_Callback(hObject, eventdata, handles)
set(handles.rbU,'Value',0);
set(handles.pnlParams,'Visible','off');
set(handles.pslide,'Enable','on');
set(handles.factorslide,'Enable','on');
plotCurrent(handles);

% --- Executes on button press in rbU.
function rbU_Callback(hObject, eventdata, handles)
set(handles.rbS,'Value',0);
set(handles.pnlParams,'Visible','on');
set(handles.pslide,'Enable','off');
set(handles.factorslide,'Enable','off');
updateUserParam(handles);
