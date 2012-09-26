function varargout = DEIMEstimatorAnalyzer(varargin)
% DEIMESTIMATORANALYZER MATLAB code for DEIMEstimatorAnalyzer.fig
%      DEIMESTIMATORANALYZER, by itself, creates a new DEIMESTIMATORANALYZER or raises the existing
%      singleton*.
%
%      H = DEIMESTIMATORANALYZER returns the handle to a new DEIMESTIMATORANALYZER or the handle to
%      the existing singleton*.
%
%      DEIMESTIMATORANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEIMESTIMATORANALYZER.M with the given input arguments.
%
%      DEIMESTIMATORANALYZER('Property','Value',...) creates a new DEIMESTIMATORANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DEIMEstimatorAnalyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DEIMEstimatorAnalyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DEIMEstimatorAnalyzer

% Last Modified by GUIDE v2.5 25-Sep-2012 16:33:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DEIMEstimatorAnalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @DEIMEstimatorAnalyzer_OutputFcn, ...
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


% --- Executes just before DEIMEstimatorAnalyzer is made visible.
function DEIMEstimatorAnalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DEIMEstimatorAnalyzer (see VARARGIN)

% Choose default command line output for DEIMEstimatorAnalyzer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

r = varargin{1};
% Care condition: dont simulate very expensive by default!
e = r.ErrorEstimator;
if e.JacSimTransMaxSize > 0 && e.JacSimTransSize == 0
    e.JacSimTransSize = 1;
end
setappdata(handles.main,'r',r);
createModelParamSliders(handles,r);
modelToGUI(handles, r);

extra = '';
sr = r.FullModel.SpaceReducer;
if sr.IncludeTrajectoryFxiData
    extra = 'InclFxiData:yes';
end
if sr.IncludeFiniteDifferences
    extra = [extra 'InclFinDiff:yes'];
end
n = sprintf('Model "%s", %d-dim, reduced %d-dim, %s',...
    r.Name,r.FullModel.Dimension,size(r.V,2),extra);
set(hObject,'Name',n);
set(handles.rbgAlpha,'SelectionChangeFcn',@(s,e)rbgAlpha_SelectionChanged(handles,s,e));
set(handles.rbgBeta,'SelectionChangeFcn',@(s,e)rbgBeta_SelectionChanged(handles,s,e));
set(handles.rbg_Param,'SelectionChangeFcn',@(s,e)rbgParam_SelectionChanged(handles,s,e));

setappdata(handles.main,'mu',r.ParamSamples(:,1));
in = [];
if r.System.InputCount > 0
    in = 1;
end
setappdata(handles.main,'inputidx',in);

reSimulate(handles);

function createModelParamSliders(h,r)
% Setup training sample slider
set(h.slTrainParam,'Min',1,'Max',size(r.ParamSamples,2),'Value',1);
% Create model parameter slides
parent = h.rbg_Param;
set(parent,'Units','pixels');
pos = get(parent,'Position');
set(parent,'Units','normalized');
dist = 22;%px
pcnt = 0;
pdata = struct;
pdata.sl = zeros(r.System.ParamCount,1);
pdata.lbl = zeros(r.System.ParamCount,1);
for pidx = 1:r.System.ParamCount
    p = r.System.Params(pidx);
    if p.HasRange
        % Top location
        top = pos(4)-pcnt*dist-40-2*dist;
        % Create labels
        label = uicontrol('Tag',sprintf('lblP%d',pidx),'Style','text',...
                'Parent',parent,'HorizontalAlignment','left');
        set(label,'String',p.Name,'Units','pixels', 'Position',[10 top 100 14]);
        set(label,'Units','normalized');
        pdata.lbl(pidx) = label;
        % Labels for min/max values
        label = uicontrol('Tag',['runtime_lbll' num2str(pidx)],'Style','text',...
                'Parent',parent,'HorizontalAlignment','right');
        set(label,'String',sprintf('%g',p.MinVal),'Units','pixels',...
            'Position',[100 top 30 14]);
        set(label,'Units','normalized');
        label = uicontrol('Tag',['runtime_lblr' num2str(pidx)],'Style','text',...
                'Parent',parent,'HorizontalAlignment','left');
        set(label,'String',sprintf('%g',p.MaxVal),'Units','pixels',...
            'Position',[320 top 30 14]);
        set(label,'Units','normalized');
        % Create slider
        ctrl = uicontrol('Tag',sprintf('slP%d',pidx),'Parent',...
            parent,'Style','slider','UserData',pidx);
        % Position
        set(ctrl,'Units','pixels','Position',[140 top 170 16]);
        set(ctrl,'Units','normalized');
        % Range etc
        set(ctrl,'Min',p.MinVal,'Max',p.MaxVal,'Value',p.MinVal);
        set(ctrl,'SliderStep',[0.01 0.1]);
        % Set callback & string
        set(ctrl,'Callback',@(~,~)(updateUserParam(h)));
        pdata.sl(pidx) = ctrl;
       
        % increase position counter
        pcnt = pcnt + 1;
    end
end
setappdata(h.main,'pdata',pdata);

function updateUserParam(h)
r = getappdata(h.main,'r');
pdata = getappdata(h.main,'pdata');
pc = r.System.ParamCount;
if pc > 0
    set(h.rbCustomParam,'Value',1);
    mu = zeros(pc,1);
    for pidx = 1:pc
        if r.FullModel.System.Params(pidx).HasRange
            mu(pidx) = get(pdata.sl(pidx),'Value');
        else
            mu(pidx) = r.FullModel.System.Params(pidx).MinVal;
        end
        s = sprintf('%s: %g',r.System.Params(pidx).Name,mu(pidx));
        set(pdata.lbl(pidx),'String',s);
    end
    setappdata(h.main,'mu',mu);
end
reSimulate(h);

function modelToGUI(h, r)
    
    set(h.slM,'Min',1,'Max',r.System.f.MaxOrder,'Value', r.System.f.Order(1));
    set(h.lblM,'String',sprintf('%d',r.System.f.Order(1)));
    
    % m' estimates
    e = r.ErrorEstimator;
    if e.UseTrueDEIMErr
        set(h.rbTrueDEIM,'Value',1);
    else
        set(h.rbMD,'Value',1);
    end
    mv = r.System.f.MaxOrder-r.System.f.Order(1);
    if mv ~= 0, en='on'; else en='off'; end
    set(h.slMD,'Min',0,'Max',max(mv,1),'Enable',en,'Value', r.System.f.Order(2));
    set(h.lblMD,'String',sprintf('%d',r.System.f.Order(2)));
    
    % mj JacMDEIM
    set(h.slJM,'Min',1,'Max',e.JacMatDEIMMaxOrder,'Value', e.JacMatDEIMOrder);
    set(h.lblJM,'String',sprintf('%d',e.JacMatDEIMOrder));
    % simtrans k
    set(h.slK,'Min',0,'Max',e.JacSimTransMaxSize,'Value', e.JacSimTransSize);
    set(h.lblK,'String',sprintf('%d',e.JacSimTransSize));
    if e.UseTrueLogLipConst
        set(h.rbtrueloclip,'Value',1);
    elseif e.UseJacobianLogLipConst
        set(h.rbjacloclip,'Value',1);
    elseif e.UseFullJacobian
        set(h.rbfulljac,'Value',1);
    else
        set(h.rbjmst,'Value',1);
    end
    
function reSimulate(h)
drawnow;
fprintf('Simulating... ');
r = getappdata(h.main,'r');
mu = getappdata(h.main,'mu');
in = getappdata(h.main,'inputidx');
[~, y, tf, x] = r.FullModel.simulate(mu,in);
[t, yr, tr, xr] = r.simulate(mu,in);
err = Norm.L2(y-yr);
yno = Norm.L2(y);
setappdata(h.main,'x',x);
setappdata(h.main,'xr',xr);
setappdata(h.main,'t',t);
setappdata(h.main,'y',y);
setappdata(h.main,'yr',yr);
tools.LogPlot.cleverPlot(h.axerr,t,err,'b');
tools.LogPlot.cleverPlot(h.axrel,t,err./yno,'b');
if r.ErrorEstimator.Enabled
    hold(h.axerr,'on');
    tools.LogPlot.cleverPlot(h.axerr,t,r.ErrorEstimator.OutputError,'r');
    hold(h.axerr,'off');
    % Relative error
    hold(h.axrel,'on');
    tools.LogPlot.cleverPlot(h.axrel,t,r.ErrorEstimator.OutputError./yno,'r');
    hold(h.axrel,'off');
    % Effectivities
    tools.LogPlot.cleverPlot(h.axeff,t,r.ErrorEstimator.OutputError./err,'g');
    eest = r.ErrorEstimator.OutputError(end);
else
    cla(h.axeff);
    eest = -1;
end
str = sprintf('Times: Full %gs, Red.: %gs, Speedup: %g\n',tf,tr,tf/tr);
str = [str sprintf('Errors(T=%g): True: %g, Est.: %g,\nRel: %g, Est.Rel.: %g, Eff.:%g',...
    t(end),err(end),eest,err(end)/yno(end),eest/yno(end),eest/err(end))];
set(h.lblRes,'String',str);
tools.LogPlot.cleverPlot(h.axhlp,r.Times(2:end),r.ODESolver.LastAlpha,'m');
tools.LogPlot.cleverPlot(h.axhlp2,r.Times(2:end),r.ODESolver.LastBeta,'c');
axis(h.axerr,'tight');
if r.ErrorEstimator.Enabled
    axis(h.axeff,'tight');
    axis(h.axhlp,'tight');
    axis(h.axhlp2,'tight');
end
fprintf('done.\n');

function getDEIMErrorsOnTraj(h)
r = getappdata(h.main,'r');
xr = getappdata(h.main,'xr');
mu = getappdata(h.main,'mu');
[ef,er,fxno] = testing.DEIM.getApproxErrorFullRed(r, xr, r.Times, mu, r.V);
pm = tools.PlotManager;
pm.LeaveOpen = true;
ax = pm.nextPlot('deimerr','True and projected approximation error on current trajectory','time','error');
if max(ef)/min(ef) < 100 || max(er)/min(er) < 100
    plotfun = @plot;
else
    plotfun = @semilogy;
end
plotfun(ax,r.Times,ef,'b',r.Times,er,'r');
legend('Full','Projected');
pm.done;

% --- Outputs from this function are returned to the command line.
function varargout = DEIMEstimatorAnalyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function rbgAlpha_SelectionChanged(h, source, ~)
rb = get(source,'SelectedObject');
r = getappdata(h.main,'r');
e = r.ErrorEstimator;
switch get(rb,'UserData')
    case 1
        e.UseTrueDEIMErr = true;
        e.Enabled = true;
    case 2
        e.UseTrueDEIMErr = false;
end
reSimulate(h);

function rbgBeta_SelectionChanged(h,source,eventdata)
rb = get(source,'SelectedObject');
r = getappdata(h.main,'r');
e = r.ErrorEstimator;
e.UseTrueLogLipConst = false;
e.UseJacobianLogLipConst = false;
e.UseFullJacobian = false;
switch get(rb,'UserData')
    case 1
        e.UseTrueLogLipConst = true;
    case 2
        e.UseJacobianLogLipConst = true;
    case 3
        e.UseFullJacobian = true;
end
reSimulate(h);

function rbgParam_SelectionChanged(h, source, ~)
rb = get(source,'SelectedObject');
switch get(rb,'UserData')
    case 1
        slTrainParam_Callback(h.slTrainParam, [], h);
    case 2
        updateUserParam(h);
end

% --- Executes on slider movement.
function slMD_Callback(hObject, eventdata, handles)
r = getappdata(handles.main,'r');
md = round(get(hObject,'Value'));
if r.System.f.Order(2) ~= md
    fprintf('Setting m''=%d... ',md);
    r.System.f.Order = [r.System.f.Order(1) md];
    r.ErrorEstimator.Enabled = md > 0;
    set(handles.lblMD,'String',sprintf('%d',md));
    set(handles.rbMD,'Value',1);
    reSimulate(handles);
end

% --- Executes during object creation, after setting all properties.
function slMD_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slJM_Callback(hObject, eventdata, handles)
r = getappdata(handles.main,'r');
jm = round(get(hObject,'Value'));
if r.ErrorEstimator.JacMatDEIMOrder ~= jm
    fprintf('Setting m_j=%d... ',jm);
    r.ErrorEstimator.JacMatDEIMOrder = jm;
    set(handles.lblJM,'String',sprintf('%d',jm));
    reSimulate(handles);
end

% --- Executes during object creation, after setting all properties.
function slJM_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slK_Callback(hObject, eventdata, handles)
r = getappdata(handles.main,'r');
k = round(get(hObject,'Value'));
if r.ErrorEstimator.JacSimTransSize ~= k
    fprintf('Setting k=%d... ',k);
    r.ErrorEstimator.JacSimTransSize = k;
    set(handles.lblK,'String',sprintf('%d',k));
    reSimulate(handles);
end

% --- Executes during object creation, after setting all properties.
function slK_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slM_Callback(hObject, eventdata, handles)
r = getappdata(handles.main,'r');
m = round(get(hObject,'Value'));
if r.System.f.Order(1) ~= m
    fprintf('Setting m=%d... ',m);
    r.System.f.Order = [m min(r.System.f.Order(2),r.System.f.MaxOrder-m)];
    set(handles.lblM,'String',sprintf('%d',m));

    % Update range for m' values
    mv = r.System.f.MaxOrder-m;
    if mv ~= 0, en='on'; else en='off'; end
    set(handles.slMD,'Min',0,'Max',max(mv,1),'Enable',en,'Value', r.System.f.Order(2));
    set(handles.lblMD,'String',sprintf('%d',r.System.f.Order(2)));
    reSimulate(handles);
end

% --- Executes during object creation, after setting all properties.
function slM_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnDEIMErr.
function btnDEIMErr_Callback(hObject, eventdata, handles)
getDEIMErrorsOnTraj(handles);


% --- Executes on button press in btnRedSummary.
function btnRedSummary_Callback(hObject, eventdata, h)
r = getappdata(h.main,'r');
ma = tools.ModelAnalyzer(r);
pm = ma.plotReductionOverview;
pm.setFigureNames(sprintf('Reduction summary of "%s"',get(h.main,'Name')));


% --- Executes on slider movement.
function slTrainParam_Callback(hObject, ~, h)
% hObject    handle to slTrainParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
r = getappdata(h.main,'r');
pidx = round(get(hObject,'Value'));
mu = r.ParamSamples(:,pidx);
if ~isequal(getappdata(h.main,'mu'),mu)
    setappdata(h.main,'mu',mu);
    set(h.rbTrainParam,'Value',1);
    pdata = getappdata(h.main,'pdata');
    for pidx = 1:length(mu)
        set(pdata.sl(pidx),'Value',mu(pidx));
        s = sprintf('%s: %g',r.System.Params(pidx).Name,mu(pidx));
        set(pdata.lbl(pidx),'String',s);
    end
    drawnow;
    reSimulate(h);
end

% --- Executes during object creation, after setting all properties.
function slTrainParam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slTrainParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnVisResult.
function btnVisResult_Callback(hObject, eventdata, h)
% hObject    handle to btnVisResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
r = getappdata(h.main,'r');
m = r.FullModel;
t = getappdata(h.main,'t');
if get(h.chkAbsErr,'Value')
    r.plot(t,getappdata(h.main,'yr'));
    set(gcf,'Name','Results of reduced simulation');
    m.plot(t,abs(getappdata(h.main,'y')-getappdata(h.main,'yr')));
    set(gcf,'Name','Absolute errors');
else
    m.plot(t,getappdata(h.main,'y'));
    set(gcf,'Name','Results of full simulation');
    r.plot(t,getappdata(h.main,'yr'));
    set(gcf,'Name','Results of reduced simulation');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
r = getappdata(h.main,'r');
m = r.FullModel;
t = getappdata(h.main,'t');
v = Norm.L2(getappdata(h.main,'x')-r.V*getappdata(h.main,'xr'));
pm = tools.PlotManager;
pm.LeaveOpen = true;
h = pm.nextPlot('','Trajectory subspace projection L2(state) error','time','error');
tools.LogPlot.cleverPlot(h,t,v,'r','LineWidth',2);
pm.done;
