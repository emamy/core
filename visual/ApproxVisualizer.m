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

% Last Modified by GUIDE v2.5 05-May-2011 15:32:05

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
if isempty(m.Approx)
    error('The Approx field of the reduced models full model must not be empty.');
end
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
    setappdata(handles.main,'inidx',m.TrainingInputs(1));
end
setappdata(handles.main,'fidx',1);
setappdata(handles.main,'mu',[]);

% Plot default visualization
plotCurrent(handles);


function modelToGUI(r, h)
m = r.FullModel;
% Disable controls specific to model settings
pc = size(r.ParamSamples,2);
if m.System.ParamCount == 0 || pc <= 1
    set(h.pslide,'Enable','off');
else
    set(h.pslide,'Max',pc,'Min',1,'Value',1);
    set(h.pslide,'SliderStep',[1/(pc-1) 10/(pc-1)]);
end
if m.System.InputCount <= 1 
    set(h.islideF,'Enable','off');
else
    set(h.islideF,'Max',m.System.InputCount,'Min',1,'Value',1);
    set(h.islideF,'SliderStep',[1/m.System.InputCount 10/m.System.InputCount]);
end
if m.TrainingInputCount <= 1 
    set(h.islideS,'Enable','off');
else
    set(h.islideS,'Max',m.TrainingInputCount,'Min',1,'Value',1);
    set(h.islideS,'SliderStep',[1/m.TrainingInputCount 10/m.TrainingInputCount]);
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
for pidx = 1:m.System.ParamCount
    p = m.System.Params(pidx);
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
pc = r.FullModel.System.ParamCount;
if pc > 0
    mu = zeros(pc,1);
    slides = findobj(h.pnlParams,'Style','slider');
    for n=1:pc
        if r.FullModel.System.Params(n).HasRange
            mu(n) = get(slides(pc-n+1),'Value');
        else
            mu(n) = r.FullModel.System.Params(n).MinVal;
        end
    end
    setappdata(h.main,'mu',mu);
    currentToGUI(r,h);
end

function plotCurrent(h)

% Get approximation training data (atd)
r = getappdata(h.main,'r');
inidx = getappdata(h.main,'inidx');

if get(h.rbS,'Value') == 1
    d = r.FullModel.Data;
    
    if get(h.chkFullData,'Value') == 1
        sn = d.TrainingData;
        if ~isempty(r.V) && ~isempty(r.W)
            sn(4:end,:) = r.V*(r.W'*sn(4:end,:));
        end        
    else
        sn = d.ApproxTrainData;
    end
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
    
    if get(h.chkFullData,'Value') == 1
        fx = getappdata(h.main,'fxall');
    else
        fx = d.ApproxfValues;
    end
    sn = sn(:,sel);
    fx = fx(:,sel);
    
    xi = sn(4:end,:);
    ti = sn(3,:);
    mui = d.getParams(sn(1,:));
    
    % Find the currently used center vectors of the expansion, if set
    a = r.FullModel.Approx;
    if isa(a, 'dscomponents.AKernelCoreFun')
        cen = a.Centers.xi;
    elseif isa(a, 'approx.TPWLApprox')
        cen = a.xi;
    else
        cen = [];
    end
    centers = unique(general.Utils.findVecInMatrix(xi,cen));
    if ~isempty(centers) && centers(1) == 0
        centers(1) = [];
    end
    setappdata(h.main,'centers',centers);

else
    mu = getappdata(h.main,'mu');
    [ti,xi] = r.FullModel.computeTrajectory(mu,inidx);
    n = length(ti);
    mui = repmat(mu,1,n);
    
    fx = r.FullModel.System.f.evaluate(xi,ti,mui);
end

afx = r.FullModel.Approx.evaluate(xi,ti,mui);
setappdata(h.main,'fx',fx);
setappdata(h.main,'afx',afx);
setappdata(h.main,'ti',ti);
setappdata(h.main,'xi',xi);
setappdata(h.main,'mui',mui);

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
errinf = max(abs(fxi-afxi));

plot(ti,fxi,'r',ti,afxi,'b');
legend('Orig. Fcn','Approx');

% Plot extra markers for sampled data
if get(h.rbS,'Value') == 1
    nz = size(ti,2);
    idx = 1:nz;
    if get(h.chkCenters,'Value')
        centers = getappdata(h.main,'centers');
        if ~isempty(centers)
            hold on;
            plot(ti(idx(centers)),afxi(idx(centers)),'o','MarkerEdgeColor','black',...
                'MarkerFaceColor','r','MarkerSize',4);
            hold off;
        end
    end
end
axis tight;
title(sprintf('f,f'' at dimension %d. Linf-Err:%5.3e, L2-Err:%5.3e',fidx,errinf,err));
xlabel('t');
ylabel('f(x)');

currentToGUI(getappdata(h.main,'r'),h);

function currentToGUI(r,h)
pidx = getappdata(h.main,'pidx');
if ~isempty(pidx)
    if get(h.rbS,'Value') == 1
        set(h.txtPSampled,'String',['Param: [' sprintf('%2.3f ',r.ParamSamples(:,pidx)) ']']);
    else
        set(h.txtPFull,'String',['Param: [' sprintf('%2.3f ',getappdata(h.main,'mu')) ']']);
    end
end
inidx = getappdata(h.main,'inidx');
if ~isempty(inidx)
    r = getappdata(h.main,'r');
    if get(h.rbS,'Value') == 1
        % Might have more inputs than have been sampled.
        if inidx > r.FullModel.TrainingInputCount
            inidx = r.FullModel.TrainingInputs(end);
            setappdata(h.main,'inidx',inidx);
            set(h.txtISampled,'Value',inidx);
        end
        set(h.txtISampled,'String',sprintf('Input: %d/%d',inidx,max(r.FullModel.TrainingInputs)));
    else
        set(h.txtIFull,'String',sprintf('Input: %d/%d',inidx,r.System.InputCount));
    end
end
set(h.txtf,'String',sprintf('%d',getappdata(h.main,'fidx')));

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
function pslide_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to pslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function islideS_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to islideS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
idx = round(get(hObject,'Value'));
r = getappdata(handles.main,'r');
idx = r.FullModel.TrainingInputs(idx);
setappdata(handles.main,'inidx',idx);
plotCurrent(handles);

% --- Executes during object creation, after setting all properties.
function islideS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to islideS (see GCBO)
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

% --- Executes on button press in rbS.
function rbS_Callback(hObject, eventdata, h)
set(h.rbU,'Value',0);
set(h.pnlParams,'Visible','off');
r = getappdata(h.main,'r');
if size(r.ParamSamples,2) > 1
    set(h.pslide,'Enable','on');
end
if r.FullModel.TrainingInputCount > 1
    set(h.islideS,'Enable','on');
end
set(h.chkGlobal,'Visible','on');
set([h.islideF h.txtIFull h.txtPFull],'Enable','off');
set([h.txtISampled h.txtPSampled h.chkCenters h.chkFullData],'Enable','on');
plotCurrent(h);

function rbU_Callback(hObject, eventdata, h)
set(h.rbS,'Value',0);
set(h.pnlParams,'Visible','on');
set(h.chkGlobal,'Visible','off','Value',0);
set([h.txtIFull h.txtPFull],'Enable','on');
set([h.pslide h.txtISampled h.txtPSampled h.islideS h.chkCenters h.chkFullData],'Enable','off');
r = getappdata(h.main,'r');
if r.System.InputCount > 1
    set(h.islideF,'Enable','on');
end
updateUserParam(h);
plotCurrent(h);

function chkFullData_Callback(hObject, eventdata, handles)
if get(hObject,'Value') == 1 && isempty(getappdata(handles.main,'fxall'))
    r = getappdata(handles.main,'r');
    d = r.FullModel.Data;
    sn = d.TrainingData;
    atdidx = r.FullModel.Approx.TrainDataSelector.LastUsed;
    
    fxall = zeros(size(sn)-[3 0]);
    fxall(:,atdidx) = d.ApproxfValues;
    [tocomp, pos] = setdiff(1:size(sn,2),atdidx);
    if ~isempty(tocomp)
        fxall(:,pos) = r.FullModel.System.f.evaluate(sn(4:end,tocomp),sn(3,tocomp),...
            d.getParams(sn(1,tocomp)));
    end
    setappdata(handles.main,'fxall',fxall);
end
plotCurrent(handles)

function chkCenters_Callback(hObject, eventdata, handles) %#ok<*INUSL>
plotCurrentFxi(handles)

function btnShowErr_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowErr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = handles;
if get(h.rbS,'Value') == 1
    if get(h.chkFullData,'Value') == 1
        showErrors(h, 'Full training data');
    else
        showErrors(h, 'Approx training data');
    end
else
    showErrors(h, 'Current trajectory');
end

function showErrors(h, name)
r = getappdata(h.main,'r');
m = r.FullModel;

if get(h.rbS,'Value') == 1 && get(h.chkGlobal,'Value') == 1
    name = [name ' [global]'];
    d = m.Data;
    if get(h.chkFullData,'Value')
        f = getappdata(h.main,'fxall');
        sn = d.TrainingData;
    else
        f = d.ApproxfValues;
        sn = d.ApproxTrainData;
    end
    xi = sn(4:end,:);
    ti = sn(3,:);
    mui = d.getParams(sn(1,:));
else
    ti = getappdata(h.main,'ti');
    xi = getappdata(h.main,'xi');
    mui = getappdata(h.main,'mui');
    f = getappdata(h.main,'fx');
end

fa = m.Approx.evaluate(xi,ti,mui);
if ~isempty(m.SpaceReducer)
    fap = r.V*r.System.f.evaluate(r.W'*xi,ti,mui);
else
    fap = r.System.f.evaluate(xi,ti,mui);
end

if get(h.rbLinf,'Value')
    errfun = @getLInftyErr;
    fun = 'Linf';
elseif get(h.rbL2,'Value')
    errfun = @getL2Err;
    fun = 'L2';
end
figure;
[v,i,e] = errfun(f,fa);
plotit(e,1)
title(sprintf('%s\nFull vs full approx %s-errors\nMax: %5.3e, mean: %5.3e',name,fun,v,mean(e)));

[v,i,e] = errfun(f,fap);
plotit(e,2)
title(sprintf('%s\nFull vs projected approx %s-errors\nMax: %5.3e, mean: %5.3e',name,fun,v,mean(e)));

[v,i,e] = errfun(fa,fap);
plotit(e,3)
title(sprintf('%s\nFull approx vs projected approx %s-errors\nMax: %5.3e, mean: %5.3e',name,fun,v,mean(e)));
    
function plotit(e,nr)
subplot(1,3,nr);
%plot(e,'r.','MarkerSize',.5);
plot(e,'r');
axis tight;

function [val,idx,errs] = getLInftyErr(a,b)
% computes the 'L^\infty'-approximation error over the
% training set for the current approximation
errs = max(abs(a-b),[],1);
[val, idx] = max(errs);

function [val,idx,errs] = getL2Err(a,b)
% computes the 'L^\infty'-approximation error over the
% training set for the current approximation
errs = sqrt(sum((a-b).^2));
[val, idx] = max(errs);


% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
general.Utils.saveAxes(handles.image, '');


% --- Executes on slider movement.
function islideF_Callback(hObject, eventdata, handles)
% hObject    handle to islideF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
setappdata(handles.main,'inidx',round(get(hObject,'Value')));
plotCurrent(handles);

% --- Executes during object creation, after setting all properties.
function islideF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to islideF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
