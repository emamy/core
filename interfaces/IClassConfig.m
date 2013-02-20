classdef IClassConfig < KerMorObject
% IClassConfig: Abstract interface for a set of configurations that can be applied to a given
% algorithm
%
% See also: kernels.RBFConfig general.regression.EpsSVRConfig
%
% @author Daniel Wirtz @date 2012-11-22
%
% @new{0,7,dw,2012-11-22} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties
        % The index of the configuration with the best results.
        %
        % Set inside the algorithm this configuration set is used for.
        %
        % @propclass{verbose} Helps identifying the effectively used configuration.
        %
        % @type integer @default []
        vBestConfigIndex = [];
    end
    
    methods(Sealed)
        function lbl = getAxisLabels(this, nrs)
            if nargin < 2
                nrs = 1:this.getNumConfigurations;
            end
            lbl = arrayfun(@convert,nrs,'Unif',false);
            
            function o = convert(e)
                tmp = this.getConfigurationString(e, true);
                o = Utils.implode(tmp,sprintf('/'));
            end
        end
        
        function t = getValueRanges(this)
            t = PrintTable;
            t.HasHeader = true;
            t.HasRowHeader = true;
            t.addRow('Location','Min','Max');
            this.collectRanges(t,{this.getClassName});
        end
    end
    
%     methods
%         function copy = clone(this, copy)
%             copy = clone@KerMorObject(this, copy);
%             copy.vBestConfigIndex = this.vBestConfigIndex;
%         end
%     end
    
    methods(Access=protected)
        function idx = getPartIndices(this, partNr, totalParts)
            if partNr > totalParts
                error('partNr <= totalParts required.');
            end
            n = this.getNumConfigurations;
            ptsize = ceil(n/totalParts);
            idx = ((partNr-1)*ptsize+1):min((partNr*ptsize),n);
        end
        
        function addRange(~, ptable, proppath, minval, maxval)
            head = Utils.implode(proppath,'.');
            ptable.addRow(head,minval,maxval);
        end
    end
    
    methods(Abstract)
        % Returns the number of configurations that can be applied
        %
        % Return values:
        % n: The number of configurations @type integer
        n = getNumConfigurations(this);
        
        % Returns the number of configurations that can be applied
        %
        % Parameters:
        % nr: The configuration number @type integer
        % object: The class object for which to apply the configuration @type handle
        applyConfiguration(this, nr, object);
        
        % Returns the number of configurations that can be applied
        %
        % Parameters:
        % nr: The configuration number @type integer
        % asCell: Flag to indicate that each setting that can be done should be placed in a
        % cell of a cell array. @type logical
        %
        % Return values:
        % str:  @type integer
        str = getConfigurationString(this, nr, asCell);
        
        % Returns a string of the changed properties by this IClassConfig instance
        %
        % Return values:
        % str: The string @type char
        str = getConfiguredPropertiesString(this);
        
        % Returns a sub-part of this configuration as a new instance.
        %
        % Use the helper method getPartIndices to obtain the correct indices of the
        % configurations that belong to a specific part.
        %
        % Parameters:
        % partNr: The part number @type integer
        % totalParts: The total number of parts @type integer
        %
        % Return values:
        % conf: A copy containing the configurations of the specified part @type IClassConfig
        conf = getSubPart(this, partNr, totalParts);
    end
    
    methods(Abstract, Access=protected)
        collectRanges(this, ptable, proppath);
    end
    
    methods(Static)
        function test_ClassConfigPlots
            
            pm = PlotManager(false,2,2);
            pm.LeaveOpen = true;
            runTest(kernels.config.RBFConfig('G',.4:.01:.6));
            runTest(kernels.config.GaussConfig('G',1:10));
            runTest(kernels.config.WendlandConfig('G',1:5,'S',(1:5)/2,'D',2));
            
            pm.done;
            
            function runTest(c)
                fprintf('%s: %s',c.getClassName,c.getConfiguredPropertiesString);
                nc = c.getNumConfigurations;
                x = 1:nc;
                fx = ones(size(x));
                h = pm.nextPlot(c.getClassName,c.getClassName);
                plot(h,x,fx);
                set(h,'XTick',1:nc,'XTickLabel',c.getAxisLabels);
                disp(c.getAxisLabels);
            end            
        end
    end
    
end