classdef FunVis2DHandler < handle
% FunVis2DHandler: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2013-02-15
%
% @new{0,7,dw,2013-02-15} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=private)
        handles;
        Figure;
    end
    
    methods
        function this = FunVis2DHandler(handleFunVis2DFigure)
            this.Figure = handleFunVis2DFigure;
            this.handles = guidata(handleFunVis2DFigure);
        end
        
        function setPlotFunction(this, value)
            set(this.handles.chkPlotFun,'Value',value);
            this.updatePlot;
        end
        
        function setDisplayedCentersPerc(this, value)
            ch = this.handles.slCenters;
            set(ch,'Value',value);
            FunVis2D('slCenters_Callback',ch,[],this.handles);
        end
        
        function setDisplayedTrainingPointsPerc(this, value)
            ch = this.handles.slPerc;
            set(ch,'Value',value);
            FunVis2D('slPerc_Callback',ch,[],this.handles);
        end
        
        function setGridElems(this, num)
            ch = this.handles.slRefine;
            set(ch,'Value',num);
            FunVis2D('slRefine_Callback',ch,[],this.handles);
        end
        
        function setInputDim1(this, value)
            l = this.handles.dim1;
            set(l,'Value',value);
            FunVis2D('dim1_Callback',l,[],this.handles);
        end
        
        function setInputDim2(this, value)
            l = this.handles.dim2;
            set(l,'Value',value);
            FunVis2D('dim2_Callback',l,[],this.handles);
        end
        
        function setOutputDim(this, value)
            l = this.handles.dout;
            set(l,'Value',value);
            FunVis2D('dout_Callback',l,[],this.handles);
        end
        
        function setFreeDim(this, dim, value)
            % Get app data and update
            h = this.handles;
            conf = getappdata(h.main,'conf');
            if value < conf.box(dim,1)
                value = conf.box(dim,1);
                fprintf(2,'Minimum value of %g allowed.\n',conf.box(dim,1));
            elseif value > conf.box(dim,2)
                value = conf.box(dim,2);
                fprintf(2,'Maximum value of %g allowed.\n',conf.box(dim,2));
            end
            
            conf.basex(dim) = value;
            setappdata(h.main,'conf',conf);

            % Update GUI
            ch = get(this.Figure,'Children');
            sh = findobj(ch,'UserData',dim);
            if ~isempty(sh)
                set(sh,'Value',value);
                tag = get(sh,'Tag');
                pos = strfind(tag,'_');
                lbltag = sprintf('runtime_lbl_val_%s',tag(pos(end)+1:end));
                lbl = findobj(ch,'Tag',lbltag);
                set(lbl,'String',sprintf('%2.4e',conf.basex(dim)));
            end
            
            % Call other updating methods
            FunVis2D('updateATDPoints',h,conf);
            conf = FunVis2D('updateCenterPoints',h,conf);
            FunVis2D('updateFX', h, conf);
        end
        
        function setPlotCenterLines(this, value)
            set(this.handles.chkPlotCenterLines,'Value',value);
            this.updatePlot;
        end
        
        function copyAxesTo(this, ax)
            % Copy children
            copyobj(get(this.handles.ax,'Children'), ax);
            % Set same view
            set(ax,'View',get(this.handles.ax,'View'));
            zlabel(ax, this.getOutputDimLabel);
        end
        
        function setView(this, az, el)
            set(this.handles.ax,'View',[az el]);
        end
        
        function lbl = getOutputDimLabel(this)
            sv = get(this.handles.dout,'String');
            lbl = sv(get(this.handles.dout,'Value'));
        end
    end
    methods(Access=private)
        function h = getHandle(name)
            h = findobj(this.fh,'Tag',name);
        end
        
        function updatePlot(this)
            FunVis2D('plotCurrent',this.handles,getappdata(this.handles.main,'conf'));
        end
    end
    
end