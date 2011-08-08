function varargout = KernelExp2D(varargin)
% KERNELEXP2D (M-file for GUI KernelExp2D.fig)
%
% User interface to visualize kernel expansions using two input dimensions and one putput dimension.
%
% Parameters:
% kexp: The kernel expansion, @type kernels.KernelExpansion
% xi: The space/state centers `x_i`
% fxi: The function evaluations `f(x_i)`
% lbl: [optional] A struct with fields 'x' and 'fx', containing char cell arrays with the names for
% the different dimensions. If not given, the labels `x_1,\ldots,x_n` and `f(x_1),\ldots,f(x_n)` are
% used.
% custcallb: A struct array of custom function callbacks for specific input/output dimension
% choices. Each struct must have a field 'fcn', which contains a function handle to a function `g`
% taking two arguments `x_1,x_2` and yield scalar output `g(x_1,x_2)`. The field 'sel' is required
% to denote at which dimension choices those functions should be plotted additionally. It contains a
% `1\times 3` vector `[i,j,k]` selecting the `[x_i,x_j,f_k(x)]` dimensions. Additionally, the fields
% 'col' and 'alpha' can be used to define the color as in LineSpec and the alpha value,
% respectively.
%
% @author Daniel Wirtz @date 2011-07-26
%
% @change{0,5,dw,2011-08-08} 
% - Made plotting work also if a dimension is of constant values. If such a
% dimension is chosen as one of the base dimensions, a warning is issued and no plot generated.
% - Fixed display of training points if all distances are zero (exchanged "<" by "<=").
%
% @change{0,5,dw,2011-08-02} Extracted the comparison polynomials from PN7_Nils cooperation to more
% general external callbacks and created a small help description.

%      KERNELEXP2D, by itself, creates a new KERNELEXP2D or raises the existing
%      singleton*.
%
%      H = KERNELEXP2D returns the handle to a new KERNELEXP2D or the handle to
%      the existing singleton*.
%
%      KERNELEXP2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KERNELEXP2D.M with the given input arguments.
%
%      KERNELEXP2D('Property','Value',...) creates a new KERNELEXP2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KernelExp2D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KernelExp2D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help KernelExp2D

% Last Modified by GUIDE v2.5 29-Jul-2011 17:25:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KernelExp2D_OpeningFcn, ...
                   'gui_OutputFcn',  @KernelExp2D_OutputFcn, ...
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


% --- Executes just before KernelExp2D is made visible.
function KernelExp2D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KernelExp2D (see VARARGIN)

% Choose default command line output for KernelExp2D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes KernelExp2D wait for user response (see UIRESUME)
% uiwait(handles.main);

%% Custom Code
h = handles;
if length(varargin) < 3
    close(h.main);
    error('This component needs at minimum three arguments: The kernel expansion, the source data xi and the fxi values.');
end
kexp = varargin{1};
xi = varargin{2};
fxi = varargin{3};
custcallb = []; lbl = [];
if length(varargin) == 5
    lbl = varargin{4};
    custcallb = varargin{5};
elseif length(varargin) == 4
    v4 = varargin{4};
    if isfield(v4,'x')
        lbl = v4;
    else
        custcallb = v4;
    end
end
if isempty(lbl)
    % Automatically assign "useful" names to the dimensions
    lbl = struct;
    lbl.x = {};
    for i=1:size(xi,1)
        lbl.x{end+1} = ['x_' num2str(i)];
    end
    lbl.fx = {};
    for i=1:size(fxi,1)
        lbl.fx{end+1} = ['f(x_' num2str(i) ')'];
    end
end
if ~isempty(custcallb) && (~isfield(custcallb,'fcn') || ~isfield(custcallb,'sel'))
    error('If a custom callback struct is passed, it must contain at least the fields ''fcn'' and ''sel''.');
    help KernelExp2D;
end

setappdata(h.main,'kexp',kexp);
setappdata(h.main,'xi',xi);
setappdata(h.main,'fxi',fxi);
setappdata(h.main,'cb',custcallb);

%% Setup default values
conf.dims = size(xi,1);
conf.outdims = size(fxi,1);
conf.lbl = lbl;
% Selected dimensions
conf.d1 = 1; conf.d2 = 2;
conf.dout = 1;
% Refinement factor
conf.ref = 1;
% The base x to use, starting with minimum of all values
conf.basex = min(xi,[],2);
% The percentage of nearest training points to plot along with the approximation
conf.tperc = 3;
setappdata(h.main,'conf',conf)

rotate3d(h.ax,'on');
view(h.ax,-53,48);

%% Fill dropdowns & GUI
set(h.dim1,'String',lbl.x,'Value',conf.d1);
set(h.dim2,'String',lbl.x,'Value',conf.d2);
set(h.dout,'String',lbl.fx,'Value',conf.dout);
set(h.slRefine,'Value',conf.ref);
set(h.slPerc,'Value',conf.tperc);

updateBaseSliders(h, conf);

plotCurrent(h, conf);

%% Plots the current settings
function plotCurrent(h, c)

xi = getappdata(h.main,'xi');
fxi = getappdata(h.main,'fxi');
kexp = getappdata(h.main,'kexp');
nums = round(size(xi,2)*c.ref);

xsel = [c.d1 c.d2];
x = xi(xsel,:);
xr = linspacev(min(x,[],2), max(x,[],2), nums);
x = unique([x xr]','rows')';

x1 = x(1,:);
x2 = x(2,:);

% Only plot if current selected dimensions are non-empty (constant value)
if any(x1(1) ~= x1) && any(x2(1) ~= x2)
    
    % Get triangulation
    tr = delaunay(x1,x2);

    xf = repmat(c.basex,1,size(x,2));
    xf(xsel,:) = x;

    fx = kexp.evaluate(xf);
    fx = fx(c.dout,:);

    %% Plots
    hf = figure;
    cla(h.ax);
    hold(h.ax,'on');
    trisurf(tr, x1, x2, fx,'Parent',h.ax);
    shading(h.ax,'interp');
    % Add original points (select g% subset)
    d = xi - repmat(c.basex,1,size(xi,2));
    d(xsel,:) = [];
    d = sqrt(sum(d.^2,1));
    md = min(d); Md = max(d);
    sel = d <= md + (Md-md)*(c.tperc/100);
    plot3(xi(c.d1,sel), xi(c.d2,sel), fxi(c.dout,sel),'black.','MarkerSize',5,'Parent',h.ax);

    %% Evaluate any given callbacks for the given dimension
    cb = getappdata(h.main,'cb');
    for k=1:length(cb)
        cbk = cb(k);
        if all(cbk.sel([1 2]) == xsel) && c.dout == cbk.sel(3)
            fxp = cbk.fcn(x1,x2);
            col = 'r';
            if isfield(cbk,'col')
                col = cbk.col;
            end
            alph = .7;
            if isfield(cbk,'alpha')
                alph = cbk.alpha;
            end
            hp = trisurf(tr,x1,x2,fxp,'FaceColor',col,'EdgeColor','none','Parent',h.ax);
            alpha(hp,alph);
        end
    end

    xlabel(h.ax,c.lbl.x{c.d1});
    ylabel(h.ax,c.lbl.x{c.d2});
    zlabel(h.ax,c.lbl.fx{c.dout});
    title(h.ax,sprintf('Plot of f(x)=%s against x=%s and y=%s',c.lbl.fx{c.dout},c.lbl.x{c.d1},c.lbl.x{c.d2}));
    hold(h.ax,'off');
    close(hf);
else
    warning('KerMor:KernelExp2D','One of the currently selected dimensions is empty, i.e. the values are constant. Not plotting.');
end

function res = linspacev(x,y,n)
    n = double(n);
    r = 0:n-2;
    res = [repmat(x,1,length(r)) + repmat(r,size(x,1),1) .* repmat((y-x),1,length(r))/(floor(n)-1) y];

function updateBaseSliders(h, c)

xi = getappdata(h.main,'xi');

% Remove old ones if exist
ctrls = getappdata(h.main, 'ctrls');
if ~isempty(ctrls)
    delete(ctrls.labels);
    delete(ctrls.slides);
    clear ctrls;
end

dims = 1:c.dims;
dims([c.d1 c.d2]) = []; 
cnt = 1; 
dist = 22;%px
pos = get(h.pnlBS,'Position');
ctrls.labels = [];
ctrls.slides = [];
for dim = dims
    top = pos(4)-20-cnt*dist;
    name = c.lbl.x{dim};
    m = min(xi(dim,:));
    M = max(xi(dim,:));
    
    % Create label
    label = uicontrol('Tag',['runtime_lbl_' name],'Style','text',...
        'Parent',h.pnlBS,'HorizontalAlignment','left');
    set(label,'String',[name ': '],'Units','pixels','Position',[10 top 40 14]);
    ctrls.labels(end+1) = label;
    
    label = uicontrol('Tag',['runtime_lbl_' name],'Style','text',...
        'Parent',h.pnlBS,'HorizontalAlignment','left',...
        'TooltipString',['Value of ' name ' for current plot']);
    set(label,'String',num2str(c.basex(dim)),'Units','pixels','Position',[50 top 60 14]);
    ctrls.labels(end+1) = label;
    valuelabel = label;
    
    if m ~= M
        label = uicontrol('Tag',['runtime_lbl_min' name],'Style','text',...
            'Parent',h.pnlBS,'HorizontalAlignment','right',...
            'TooltipString',['Minimum value of ' name]);
        set(label,'String',num2str(m),'Units','pixels','Position',[110 top 60 14]);
        ctrls.labels(end+1) = label;

        label = uicontrol('Tag',['runtime_lbl_max' name],'Style','text',...
            'Parent',h.pnlBS,'HorizontalAlignment','left',...
            'TooltipString',['Maximum value of ' name]);
        set(label,'String',num2str(M),'Units','pixels','Position',[370 top 60 14]);
        ctrls.labels(end+1) = label;

        %% Create slider
        ctrl = uicontrol('Tag',['runtime_slide_' name],'Parent',...
            h.pnlBS,'HorizontalAlignment','left');
        set(ctrl,'Style','slider', 'Value', c.basex(dim), 'Min', m, 'Max', M);

        % Position
        set(ctrl,'Units','pixels','Position',[170 top 200 16],'UserData',dim);
        % Set callback & string
        set(ctrl,'Callback',@(h,e)(baseSliderChanged(h,guidata(h),valuelabel)));

        cnt = cnt + 1;
        ctrls.slides(end+1) = ctrl;
    end
end
setappdata(h.main, 'ctrls', ctrls);

function baseSliderChanged(hObj, h, label)
% Get dimension that slider stands for
dim = get(hObj,'UserData');
conf = getappdata(h.main,'conf');
conf.basex(dim) = get(hObj,'Value');
set(label,'String',num2str(conf.basex(dim)));
setappdata(h.main,'conf',conf);

plotCurrent(h, conf);

% --- Outputs from this function are returned to the command line.
function varargout = KernelExp2D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;


% --- Executes on selection change in dim1.
function dim1_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
conf = getappdata(handles.main,'conf');
if get(hObject,'Value') ~= conf.d2
    conf.d1 = get(hObject,'Value');
    setappdata(handles.main,'conf',conf);

    updateBaseSliders(handles, conf);
    plotCurrent(handles, conf);
end

function dim1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in dim2.
function dim2_Callback(hObject, eventdata, handles)
conf = getappdata(handles.main,'conf');
if get(hObject,'Value') ~= conf.d1
    conf.d2 = get(hObject,'Value');
    setappdata(handles.main,'conf',conf);

    updateBaseSliders(handles, conf);
    plotCurrent(handles, conf);
end

% --- Executes during object creation, after setting all properties.
function dim2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slRefine_Callback(hObject, eventdata, handles)
% hObject    handle to slRefine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
c = getappdata(handles.main,'conf');
c.ref = get(hObject,'Value');
set(handles.lblRef,'String',num2str(c.ref));
setappdata(handles.main,'conf',c);
plotCurrent(handles, c);

% --- Executes during object creation, after setting all properties.
function slRefine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slRefine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in dout.
function dout_Callback(hObject, eventdata, handles)
% hObject    handle to dout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns dout contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dout
conf = getappdata(handles.main,'conf');
conf.dout = get(hObject,'Value');
setappdata(handles.main,'conf',conf);
plotCurrent(handles, conf);

% --- Executes during object creation, after setting all properties.
function dout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slPerc_Callback(hObject, eventdata, handles)
conf = getappdata(handles.main,'conf');
conf.tperc = get(hObject,'Value');
set(handles.lblPerc,'String',[num2str(conf.tperc) '%']);
setappdata(handles.main,'conf',conf);
plotCurrent(handles, conf);

% --- Executes during object creation, after setting all properties.
function slPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
