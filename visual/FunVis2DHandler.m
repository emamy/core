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
        
        function copyAxesTo(this, ax)
            % Copy children
            copyobj(get(this.handles.ax,'Children'), ax);
            % Set same view
            set(ax,'View',get(this.handles.ax,'View'));
        end
        
        function setView(this, az, el)
            set(this.handles.ax,'View',[az el]);
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