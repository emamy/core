function varargout = FunVis2D(varargin)
% User interface to visualize highdimensional functions using two input dimensions and one output dimension.
%
% This visualization tool can be called with KerMor model instances as well
% as standalone with any function handle or class instance (that has an
% evaluate method).
% For every case, it is assumed to have a function `f(x,t,\mu)`, where the
% last two arguments may each or both be left out; depending on the ranges
% struct (explained below) the according evaluations will be performed.
% When passing KerMor models those ranges are computed automatically; if
% one wants to plot arbitrary functions please see the description of the
% different calling interfaces below.
%
% 'FunVis2D(fullmodel)': Plot the full system's nonlinearity using the
% bounding box of the sampled trajectory data as plotting limits.
%
% 'FunVis2D(reducedmodel)': Plot the reduced model's nonlinearity (i.e. the
% approximation) using the approximation training data's bounding box as
% plotting limits. Additionally, the full model's function can be plotted
% against it or the error between them.
%
% 'FunVis2D(fun_handle, ranges)': Plots the function given by the function
% handle 'fun_handle' in the range given by 'ranges'. Ranges must be a
% struct with at minimum the field 'xrange'. Depending on the structure of
% the passed function handle, it must also contain a field 'trange' or/and
% 'murange' if any of the second arguments `t` or/and `\mu` are expected.
% So, if you leave out 'trange' but pass 'murange', the function will be
% expected to have an interface `f(x,\mu)`.
%
% 'FunVis2D(class_instance, ranges)': Expects the first argument to be a
% class with an 'evaluate' method, corresponding to the evaluation of the
% function to plot. The rest behaves like the call 'FunVis2D(fun_handle,
% ranges)'.
% @note This case also covers the passing of a kernel expansion as function
% argument, using a specified range.
%
% 'FunVis2D(kexp, ranges, trainingdata)': The same as if calling
% 'FunVis2D(class_instance, ranges)', but with two additional plotting
% possibilities: 
% - The first argument is expected to be a 'kernels.KernelExpansion'
% subclass, and the centers of it can additionally be drawn in the plot.
% - Also puts the training data points into the plot using the training
% data 'xi,ti,mui' points in conjunction with the evaluations 'fxi' of the
% full system's nonlinearity at the respective points.
%
% 'FunVis2D(kexp, trainingdata)': The same as if calling
% 'FunVis2D(kexp, ranges, trainingdata)' with the plotting ranges computed
% automatically from the bounding box of the trainingdata.
%
% 'FunVis2D(*, *, *, fun2, lbl)': Any of the above calls can receive two
% more arguments. Those must always be the 4th and 5th arguments and
% any arguments not given before must be passed as '[]'.
% - 'fun2': A second function (either handle or class, see above) to
% compare the main function to. Enables to either add the second function
% to the plot or plot the error.
% - 'lbl': Custom labels for the `x [,t [, mu]]` and `f(x[,t[,\mu]])`
% dimensions. Must be a struct with the fields 'x' and 'fx', each a char cell of length
% equal to the input and output dimensions. If more than just `x` is used
% as function argument, the order is `x`, `t` and `\mu` for the input
% dimensions; the length of the 'x' cell must always equal the sum of all
% arguments dimensions.
%
% @note The 4th argument can also be '[]', if no second function is given
% but custom labels are required.
%
% @author Daniel Wirtz @date 2011-07-26
%
% @new{0,6,dw,2011-11-15} 
% - Centers and approximation training data are now plotted completely if the only non-constant
% dimensions are the ones who are currently plotted (and the corresponding slider is `>0`).
% - Extended the center plotting with lines to the values of the full function on the training
% data if an data.ApproxTrainData class is passed.
%
% @change{0,5,dw,2011-10-24} New interface for simplified calls and some
% improvements and additions (error plots, centers and approx train data
% plots)
%
% @change{0,5,dw,2011-09-15} Fixed display of less dimensions if not the
% full space in the slider-panel is needed
%
% @change{0,5,dw,2011-09-12} Changed the base-x panel so that many
% dimensions can be displayed. A scrollbar callback dynamically creates the
% currently visible sliders. Added ParamTimeKernel expansion support (new
% calling interface with struct for those)
%
% @change{0,5,dw,2011-08-08} 
% - Made plotting work also if a dimension is of constant values. If such a
% dimension is chosen as one of the base dimensions, a warning is issued and no plot generated.
% - Fixed display of training points if all distances are zero (exchanged "<" by "<=").
%
% @change{0,5,dw,2011-08-02} Extracted the comparison polynomials from PN7_Nils cooperation to more
% general external callbacks and created a small help description.
%
% @todo 
% - fix slider display for parameters (not all shown correctly)
% - black dots for centers are not at value of f for center!
%
%
% @ingroup g_visual

%      FunVis2D, by itself, creates a new FunVis2D or raises the existing
%      singleton*.
%
%      H = FunVis2D returns the handle to a new FunVis2D or the handle to
%      the existing singleton*.
%
%      FunVis2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FunVis2D.M with the given input arguments.
%
%      FunVis2D('Property','Value',...) creates a new FunVis2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FunVis2D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FunVis2D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FunVis2D

% Last Modified by GUIDE v2.5 17-Jan-2012 11:03:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FunVis2D_OpeningFcn, ...
                   'gui_OutputFcn',  @FunVis2D_OutputFcn, ...
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


% --- Executes just before FunVis2D is made visible.
function FunVis2D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FunVis2D (see VARARGIN)

% Choose default command line output for FunVis2D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FunVis2D wait for user response (see UIRESUME)
% uiwait(handles.main);

%% Custom Code
h = handles;
ranges = struct('xrange',[],'trange',[],'murange',[]);
conf = struct;
conf.td = [];
conf.fun2 = [];
if isa(varargin{1},'models.BaseFullModel')
    m = varargin{1};
    fun = m.System.f;
    [x,X] = m.Data.getBoundingBox;
    ranges.xrange = [x X];
    ranges.trange = [0 m.T];
    if m.System.ParamCount > 0
        if isempty(m.Data.ParamSamples)
            stop(h,'No param samples available yet. Have you run model.off1_generateParamSamples?');
        end
        [mu, MU] = general.Utils.getBoundingBox(m.Data.ParamSamples);
        ranges.murange = [mu MU];
    end
elseif isa(varargin{1},'models.ReducedModel')
    r = varargin{1};
    fun = r.System.f;
    conf.td = r.FullModel.Data.ApproxTrainData;
    ranges = rangesFromATD(conf.td);
    conf.fun2 = r.FullModel.System.f;
else
    %% Handle first function argument
    fun = varargin{1};
    if isa(fun,'function_handle')
        % Wrap in core fun so that evaluate works
        fun = dscomponents.PointerCoreFun(fun);
    elseif ~ismethod(fun,'evaluate')
        stop(h,'If the first argument is not a function pointer it must be a class with an "evaluate" method.');
    end

    %% Second argument
    conf.td = [];
    sec = varargin{2};
    if isfield(sec,'xrange')
        % Use explicit ranges
        ranges.xrange = sec.xrange;
        if isfield(sec,'trange')
            ranges.trange = sec.trange;
        end
        if isfield(sec,'murange')
            ranges.murange = sec.murange;
        end
        if length(varargin) > 2 && isa(varargin{3},'data.AModelData')
            conf.td = varargin{3}.ApproxTrainData;
        end
    elseif isa(sec,'data.ApproxTrainData')
        % Use the bounds of the training data as ranges
        ranges = rangesFromATD(sec);
        conf.td = sec;
    else
        stop(h,'The second argument has to be a ranges struct or a struct with the fields xi,fxi (+ti,mui if given).');
    end
end
conf.fun = fun;
conf.ranges = ranges;

%% Build global bounding box
conf.hastime = false;
conf.hasparams = false;
box = ranges.xrange;
if ~isempty(ranges.trange)
    conf.hastime = true;
    conf.timeoff = size(ranges.xrange,1)+1;
    box = [box; ranges.trange];
end
if ~isempty(ranges.murange)
    conf.hasparams = true;
    conf.paroff = size(ranges.xrange,1)+1;
    if conf.hastime
        conf.paroff = conf.paroff+1;
    end
    box = [box; ranges.murange];
end
conf.box = box;
conf.xdim = size(ranges.xrange,1);
conf.dims = size(box,1);
% Detect constant fields
conf.const = find(max(box,[],2) == min(box,[],2));

%% Fourth argument (anyone who wants to pass second functions or labels
%% needs to fill the call with [] fields)
fun2 = []; lbl = [];
if length(varargin) > 3
    fun2 = varargin{4};
    if isa(fun2,'function_handle')
        fun2 = dscomponents.PointerCoreFun(fun2,true);
    elseif ~isempty(fun2) && ~ismethod(fun2,'evaluate')
        stop(h,'Any second function-class must have an ''evaluate''-method');
    end
    if length(varargin) > 4
        lbl = varargin{5};
        if ~isfield(lbl,'x') || ~isfield(lbl,'fx')
            stop(h, 'Custom labels must have the struct fields x and fx');
        end
    end
end
if ~isempty(fun2)
    conf.fun2 = fun2;
end

%% Use generic labels if none are given
if isempty(lbl)
    % Automatically assign "useful" names to the dimensions
    lbl = struct;
    lbl.x = {};
    for i=setdiff(1:conf.xdim,conf.const)
        lbl.x{end+1} = ['x_{' num2str(i) '}'];
    end
    args{1} = zeros(conf.xdim,1);
    if conf.hastime && all(conf.const ~= conf.timeoff)
        lbl.x{end+1} = 't';
        args{2} = 0;
    end
    if conf.hasparams
        munum = 1:size(conf.ranges.murange,1);
        [diff, idx] = setdiff(conf.paroff+munum-1,conf.const);
        for i=1:length(diff)
            lbl.x{end+1} = ['mu_{' num2str(munum(idx(i))) '}'];
        end
        args{3} = zeros(munum,1);
    end
    lbl.fx = {};
    if ~isempty(conf.td)
        fdim = size(conf.td.fxi,1);
    else
        % Autodetect output size by evaluating with all zero arguments
        fdim = size(conf.fun.evaluate(args{:}),1);
    end
    for i=1:fdim
        lbl.fx{end+1} = ['f(x_{' num2str(i) '})'];
    end
end
conf.lbl = lbl;

%% Check for kernel expansions
conf.iske = false;
conf.isptke = false;
if isa(fun,'kernels.KernelExpansion')
    conf.iske = true;
    set(h.lblCenters,'Visible','on');
    set(h.lblNumCenters,'Visible','on');
    set(h.slCenters,'Visible','on');
    conf.cperc = 3;
    set(h.slCenters,'Value',conf.cperc);
    
    % Evaluate function also on centers
    cent = conf.fun.Centers;
    if isa(fun,'kernels.ParamTimeKernelExpansion')
        conf.isptke = true;
        conf.centerfx = fun.evaluate(cent.xi,cent.ti,cent.mui);
    else
        conf.centerfx = fun.evaluate(cent.xi,[],[]);
    end
end
%% Visibilities for second function, if given
if ~isempty(conf.fun2)
    set(handles.lbl2nd,'Visible','on');
    set(handles.rbAdd,'Visible','on');
    set(handles.rbErr,'Visible','on');
end

%% Setup default values
[~, conf.idxmap] = setdiff(1:conf.dims,conf.const);
% Selected dimensions
conf.d1 = 1; conf.d2 = 2;
conf.dout = 1;
% Refinement factor
conf.gridpts = 60;
% The base x to use, starting with the center of all values
conf.basex = (box(:,1)+box(:,2))/2;

% Use training data display if any is given
if ~isempty(conf.td)
    % The percentage of nearest training points to plot along with the approximation
    conf.tperc = 3;
    set(h.slPerc,'Value',conf.tperc);
else
    % Hide training point selection slider if no td is given.
    set(h.slPerc,'Visible','off');
    set(h.lblTP,'Visible','off');
    set(h.lblPerc,'Visible','off');
end

setappdata(h.main,'conf',conf)

rotate3d(h.ax,'on');
view(h.ax,-53,48);

%% Fill dropdowns & GUI
set(h.dim1,'String',lbl.x,'Value',conf.d1);
set(h.dim2,'String',lbl.x,'Value',conf.d2);
set(h.dout,'String',lbl.fx,'Value',conf.dout);
set(h.slRefine,'Value',conf.gridpts);

setDimSliders(h);

updateATDPoints(h,conf);
conf = updateCenterPoints(h, conf);
newMesh(h, conf);

function ranges = rangesFromATD(atd)
    ranges = struct('xrange',[],'trange',[],'murange',[]);
    [x, X] = general.Utils.getBoundingBox(atd.xi);
    ranges.xrange = [x X];
    if isprop(atd,'ti') || isfield(atd,'ti')
        ranges.trange = [min(atd.ti) max(atd.ti)];
    end
    if isprop(atd,'mui') || isfield(atd,'mui')
        [mu, MU] = general.Utils.getBoundingBox(atd.mui);
        ranges.murange = [mu MU];
    end
    
function stop(h, errmsg, varargin)
    close(h.main);
    error(errmsg,varargin{:});

function newMesh(h,c)
%% Creates a new mesh and evaluates the function on it.
% Extract the effective dimension indices (const indices are not contained
% in d1,d2 indices)
xsel = c.idxmap([c.d1 c.d2]);
x1 = linspace(c.box(xsel(1),1), c.box(xsel(1),2), c.gridpts);
x2 = linspace(c.box(xsel(2),1), c.box(xsel(2),2), c.gridpts);

[X1,X2] = meshgrid(x1,x2);
setappdata(h.main,'X1',X1);
setappdata(h.main,'X2',X2);

updateFX(h, c);

function updateFX(h,c)
% Updates the fx-values
X1 = getappdata(h.main,'X1');
X2 = getappdata(h.main,'X2');
x = [X1(:)'; X2(:)'];
xf = repmat(c.basex,1,size(x,2));
xsel = c.idxmap([c.d1 c.d2]);
xf(xsel,:) = x;

% Evaluate fcn according to given values
fx2 = [];
if c.hastime && c.hasparams
    fx = c.fun.evaluate(xf(1:c.timeoff-1,:),...
        xf(c.timeoff,:),...
        xf(c.paroff:end,:));
    if ~isempty(c.fun2)
        fx2 = c.fun2.evaluate(xf(1:c.timeoff-1,:),...
        xf(c.timeoff,:),...
        xf(c.paroff:end,:));
    end
elseif c.hastime
    fx = c.fun.evaluate(xf(1:c.timeoff-1,:),...
        xf(c.timeoff,:),[]);
    if ~isempty(c.fun2)
        fx2 = c.fun2.evaluate(xf(1:c.timeoff-1,:),...
        xf(c.timeoff,:),[]);
    end
elseif c.hasparams
    fx = c.fun.evaluate(xf(1:c.paroff-1,:),[],...
        xf(c.paroff:end,:));
    if ~isempty(c.fun2)
        fx2 = c.fun2.evaluate(xf(1:c.paroff-1,:),[],...
        xf(c.paroff:end,:));
    end
else
    fx = c.fun.evaluate(xf,[],[]);
    if ~isempty(c.fun2)
        fx2 = c.fun2.evaluate(xf,[],[]);
    end
end
% Compute some errors
if ~isempty(c.fun2)
    n = (1/numel(fx));
    c.err.l2 = sqrt(n*sum((fx-fx2).^2,2));
    c.err.rel_l2_1 = c.err.l2 ./ sqrt(n*sum(fx.^2,2));
    c.err.rel_l2_2 = c.err.l2 ./ sqrt(n*sum(fx2.^2,2));
    c.err.linf = max(abs(fx-fx2),[],2);
    setappdata(h.main,'conf',c);
end
setappdata(h.main,'fx',fx);
setappdata(h.main,'fx2',fx2);

plotCurrent(h,c);

function c = updateCenterPoints(h, c)
%% Updates the currently displayed center points
if c.iske
    cent = c.fun.Centers;
    C = [cent.xi];
    if c.isptke
        C = [C; cent.ti; cent.mui];
    end
    % Determine the % closest centers
    d = C - repmat(c.basex,1,size(C,2));
    xsel = c.idxmap([c.d1 c.d2]);
    d(xsel,:) = [];
    d = sqrt(sum(d.^2,1));
    md = min(d); Md = max(d);
    if md == Md && c.cperc > 0
        sel = true(size(d));
    else
        sel = d < md + (Md-md)*(c.cperc/100)*1.001;
    end
    set(h.lblNumCenters,'String',sprintf('%d/%d',sum(sel),size(C,2)));
    C = C(:,sel);
    c.curCenters = C;
    c.curCenterSel = sel;
    xf = repmat(c.basex,1,size(C,2));
    xsel = c.idxmap([c.d1 c.d2]);
    xf(xsel,:) = C(xsel,:);
    if c.isptke
        c.curCenterFx = c.fun.evaluate(xf(1:c.timeoff-1,:),...
        xf(c.timeoff,:),...
        xf(c.paroff:end,:));
    else
        c.curCenterFx = c.fun.evaluate(xf,[],[]);
    end
    % compute error if second function is given
    if ~isempty(c.fun2) && get(h.rbErr,'Value') == 1
        c.curCenterFx = c.curCenterFx - c.fun2.evaluate(xf,[],[]);
        if c.hastime || c.hasparams
            warning('a:b','Evaluation is WRONG with params and time set. TODO.');
        end
    end
    setappdata(h.main,'conf',c);
end

function updateATDPoints(h, c)
%% Updates the currently displayed approxtraindata points
if ~isempty(c.td)
    xsel = c.idxmap([c.d1 c.d2]);
    C = [c.td.xi];
    if c.isptke
        C = [C; c.td.ti; c.td.mui];
    end
    % Determine the % closest centers
    d = C - repmat(c.basex,1,size(C,2));
    d(xsel,:) = [];
    d = sqrt(sum(d.^2,1));
    md = min(d); Md = max(d);
    if md == Md && c.tperc > 0
        sel = true(size(d));
    else
        sel = d < md + (Md-md)*(c.tperc/100)*1.001;
    end
    set(h.lblPerc,'String',sprintf('%d/%d',sum(sel),size(C,2)));
    setappdata(h.main,'selATDPoints',sel);
end

%% Plots the current settings

function plotCurrent(h, c)

fx = getappdata(h.main,'fx');
fx2 = getappdata(h.main,'fx2');
X1 = getappdata(h.main,'X1');
X2 = getappdata(h.main,'X2');
    
fx = reshape(fx(c.dout,:),size(X1,1),[]);
cap = sprintf('Plot of %s against %s and %s',c.lbl.fx{c.dout},c.lbl.x{c.d1},c.lbl.x{c.d2});

txt = '';
if ~isempty(fx2) && get(h.rbErr,'Value') == 1
    fx2 = reshape(fx2(c.dout,:),size(X1,1),[]);
    if ~any(isnan(fx2))
        cap = sprintf('Error f_1-f_2 at output %s against %s and %s',c.lbl.fx{c.dout},c.lbl.x{c.d1},c.lbl.x{c.d2});
        %fx = abs(fx-fx2);
        fx = fx-fx2;
        txt = sprintf('L2:%.2e, Linf:%.2e\nrL21:%.2e, rL22:%.2e',c.err.l2(c.dout),c.err.linf(c.dout),...
            c.err.rel_l2_1(c.dout),c.err.rel_l2_2(c.dout));
    else
        fprintf('Warning, current output %s contains NaNs. Not plotting error.\n',c.lbl.fx{c.dout});
    end
end
set(h.lblErr,'String',txt);

mi = min(fx(:));
Ma = max(fx(:));
if mi == 0 && Ma == 0
    mi = -eps; Ma = eps;
elseif mi ~= 0 && abs((mi-Ma) / mi) < 1e-14
    s = sign(mi);
    mi = (1-.001*s)*mi; Ma=(1+.001*s)*Ma;
end

%% Plots
cla(h.ax);
hold(h.ax,'on');
s1 = surf(h.ax,X1,X2,fx);

if ~isempty(fx2) && get(h.rbAdd,'Value') == 1
    fx2 = reshape(fx2(c.dout,:),size(X1,1),[]);
    if ~any(isnan(fx2))
        mi = min(mi,min(fx2(:)));
        Ma = max(Ma,max(fx2(:)));
        cap = sprintf('Plot of %s against %s and %s, both functions',c.lbl.fx{c.dout},c.lbl.x{c.d1},c.lbl.x{c.d2});
        s2 = surf(h.ax,X1,X2,fx2);
        alpha(s2,.7);
        s1 = [s1; s2];
    end
end
if get(h.rbErr,'Value') == 1
    s3 = surf(h.ax,X1,X2,zeros(size(X1)));
    alpha(s3,.7);
    s1 = [s1; s3];
end
mode = 'none';
if get(h.chkGrid,'Value') == 1
    mode = 'interp';
end
set(s1,'EdgeColor',mode);

xsel = c.idxmap([c.d1 c.d2]);
%% Plot center points if desired
if c.iske
    C = c.curCenters;
    plot3(h.ax,C(xsel(1),:),C(xsel(2),:),c.curCenterFx(c.dout,:),'black.','MarkerSize',15);
    % Also plot the centers at their original value
    orig = c.centerfx(c.dout,c.curCenterSel);
    plot3(h.ax,C(xsel(1),:),C(xsel(2),:),orig,'blackx','MarkerSize',15);
    % Plot a connecting line
    plot3(h.ax,[C(xsel(1),:); C(xsel(1),:)],[C(xsel(2),:); C(xsel(2),:)],[c.curCenterFx(c.dout,:); orig],'black');
end

%% Add training data points to plot
% With kernel expansions: option to include centers into mesh
if ~isempty(c.td) && get(h.rbErr,'Value') == 0
    C = [c.td.xi];
    if c.isptke
        C = [C; c.td.ti; c.td.mui];
    end
    C = C(xsel,:);
    sel = getappdata(h.main,'selATDPoints');
    hlpfx = c.td.fxi(c.dout,sel);
    if min(hlpfx) < mi
        mi = min(hlpfx);
    end
    if max(hlpfx) > Ma
        Ma = max(hlpfx);
    end
    plot3(h.ax,C(1,sel),C(2,sel),hlpfx,'red.','MarkerSize',15);
end
axis(h.ax,[X1(1,1) X1(1,end) X2(1,1) X2(end,1) mi Ma]);

xlabel(h.ax,c.lbl.x{c.d1});
ylabel(h.ax,c.lbl.x{c.d2});
zlabel(h.ax,c.lbl.fx{c.dout});
title(h.ax,cap);
hold(h.ax,'off');
grid(h.ax,'on');

function setDimSliders(h)
val = 1-get(h.panelslide,'Value');
c = getappdata(h.main,'conf');
% % Remove old ones if exist
ctrls = getappdata(h.main, 'ctrls');
if ~isempty(ctrls)
    delete(ctrls.labels);
    delete(ctrls.slides);
    clear ctrls;
end

pos = get(h.pnlBS,'Position');
pnlh = pos(4);

dims = c.idxmap;
dims([c.d1 c.d2]) = []; 
dist = 6;%px
height = 16; %px
eh = dist + height; % total element height
dispsliders = min(length(dims), floor(pnlh / eh)-1);
totalheight = (length(dims)-dispsliders) * eh;
if totalheight > pnlh
    curoffset = totalheight * val;
    firstidx = floor(curoffset / eh);
else
    set(h.panelslide,'Visible','off');
    firstidx = 0;
end

ctrls.labels = [];
ctrls.slides = [];
for idx = 1:dispsliders
    dim = dims(firstidx+idx);
    top = pnlh-15-idx*eh;
    name = c.lbl.x{dim};
    % Create labels
    m = c.box(dim,1);
    M = c.box(dim,2);
    tooltip = sprintf('Range [%1.5e - %1.5e]',m,M);
    % Dimension label
    label = uicontrol('Tag',['runtime_lbl_' name],'Style','text',...
        'Parent',h.pnlBS,'HorizontalAlignment','left');
    set(label,'String',[name ': '],'Units','pixels','Position',[10 top 40 14],...
        'TooltipString',tooltip);
    ctrls.labels(end+1) = label;
    
    % Value label
    label = uicontrol('Tag',['runtime_lbl_val_' name],'Style','text',...
        'Parent',h.pnlBS,'HorizontalAlignment','left',...
        'TooltipString',['Value of ' name ' for current plot']);
    set(label,'String',sprintf('%2.4e',c.basex(dim)),'Units','pixels','Position',[60 top 90 14],...
        'TooltipString',tooltip);
    ctrls.labels(end+1) = label;
    valuelabel = label;
    
%     % Min value label
%     label = uicontrol('Tag',['runtime_lbl_min' name],'Style','text',...
%         'Parent',h.pnlBS,'HorizontalAlignment','right',...
%         'TooltipString',['Minimum value of ' name]);
%     set(label,'String',sprintf('%1.1e',m),'Units','pixels','Position',[160 top 50 14]);
%     ctrls.labels(end+1) = label;
        
    % Create slider
    ctrl = uicontrol('Tag',['runtime_slide_' name],'Parent',...
        h.pnlBS,'HorizontalAlignment','left');
    set(ctrl,'Style','slider', 'Value', c.basex(dim),...
        'Min', m, 'Max', M);
    set(ctrl,'Units','pixels','Position',[160 top 250 16],'UserData',dim,...
        'TooltipString',tooltip);
    set(ctrl,'Callback',@(h,e)(baseSliderChanged(h,guidata(h),valuelabel)));
    ctrls.slides(end+1) = ctrl;
    
    % Max value label
%     label = uicontrol('Tag',['runtime_lbl_max' name],'Style','text',...
%         'Parent',h.pnlBS,'HorizontalAlignment','left',...
%         'TooltipString',['Maximum value of ' name]);
%     set(label,'String',sprintf('%1.1e',M),'Units','pixels','Position',[350 top 50 14]);
%     ctrls.labels(end+1) = label;
end
setappdata(h.main, 'ctrls', ctrls);

function baseSliderChanged(hObj, h, label)
% Get dimension that slider stands for
dim = get(hObj,'UserData');
conf = getappdata(h.main,'conf');
conf.basex(dim) = get(hObj,'Value');
set(label,'String',sprintf('%2.4e',conf.basex(dim)));
setappdata(h.main,'conf',conf);

updateATDPoints(h,conf);
conf = updateCenterPoints(h,conf);
updateFX(h, conf);

function varargout = FunVis2D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.main;

function dim1_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
conf = getappdata(handles.main,'conf');
if get(hObject,'Value') ~= conf.d2
    old = conf.d1;
    conf.d1 = get(hObject,'Value');
    setappdata(handles.main,'conf',conf);

    setDimSliders(handles);
    newMesh(handles, conf);
end

function dim1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dim2_Callback(hObject, eventdata, handles)
conf = getappdata(handles.main,'conf');
if get(hObject,'Value') ~= conf.d1
    old = conf.d2;
    conf.d2 = get(hObject,'Value');
    setappdata(handles.main,'conf',conf);

    setDimSliders(handles);
    newMesh(handles, conf);
end

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
c.gridpts = round(get(hObject,'Value'));
set(handles.lblRef,'String',num2str(c.gridpts));
setappdata(handles.main,'conf',c);
newMesh(handles, c);

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
set(handles.lblTP,'String',sprintf('Select %2.2f%% nearest training points',conf.tperc));

setappdata(handles.main,'conf',conf);
updateATDPoints(handles, conf);
plotCurrent(handles,conf);

% --- Executes during object creation, after setting all properties.
function slPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function panelslide_Callback(hObject, eventdata, h)
% hObject    handle to panelslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
setDimSliders(h);

function panelslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panelslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function uipushtool2_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
general.Utils.saveAxes(handles.ax);

function slCenters_Callback(hObject, eventdata, handles)
% hObject    handle to slCenters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
c = getappdata(handles.main,'conf');
c.cperc = get(hObject,'Value');
set(handles.lblCenters,'String',sprintf('Include %2.2f%% nearest centers',c.cperc));
setappdata(handles.main,'conf',c);
c = updateCenterPoints(handles, c);
plotCurrent(handles,c);

function slCenters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slCenters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function rbAdd_Callback(hObject, eventdata, handles)
set(handles.rbErr,'Value',0);
c = updateCenterPoints(handles, getappdata(handles.main,'conf'));
plotCurrent(handles,c);

function rbErr_Callback(hObject, eventdata, handles)
set(handles.rbAdd,'Value',0);
c = updateCenterPoints(handles, getappdata(handles.main,'conf'));
plotCurrent(handles,c);

% --- Executes on button press in chkGrid.
function chkGrid_Callback(hObject, eventdata, handles)
% hObject    handle to chkGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkGrid
plotCurrent(handles,getappdata(handles.main,'conf'));
