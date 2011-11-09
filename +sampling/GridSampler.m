classdef GridSampler < sampling.BaseSampler
% GridSampler: Samples the models params on a grid, using the
% ModelParam.Desired field
%
% @author Daniel Wirtz @date 2010-04-01
%
% @change{0,5,dw,2011-11-09}
% - Fixed logspace grid sampling when any of the range values was zero.
% - New default value for Spacing is 'lin'.
%
% @change{0,5,dw,2011-10-19} Added a linear or logarithmic spacing method.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
%
% @todo maybe move the Spacing property into ModelParam, to enable custom
% spacing for each parameter
    
    properties(SetObservable)
        % The spacing method.
        %
        % Available options are 'log' and 'lin' for logarithmic/linear grid
        % sampling.
        %
        % @propclass{important} Different model/parameters require
        % differently scaled sampling
        %
        % @type char @default lin
        Spacing = 'lin';
    end
    
    methods
        function this = GridSampler
            this = this@sampling.BaseSampler;
            this.registerProps('Spacing');
        end
        
        function samples = performSampling(this, model)
            % Uses the given model and generates a training set by creating
            % a regular grid in joint time/parameter space
            % @ingroup s_grid
            sys = model.System;
            ranges = cell(sys.ParamCount,1);
            
            % Create linearly spaced parameter value ranges
            for pidx=1:sys.ParamCount
                % Cater for single parameters 
                if (sys.Params(pidx).MinVal == sys.Params(pidx).MaxVal)
                    ranges{pidx} = sys.Params(pidx).MinVal;
                else
                    if strcmp(this.Spacing,'lin')
                        ranges{pidx} = linspace(sys.Params(pidx).MinVal,...
                        sys.Params(pidx).MaxVal,...
                        sys.Params(pidx).Desired);
                    else
                        % Avoid NaNs due to log10(0)
                        if sys.Params(pidx).MinVal == 0
                            m = -16; % take double precision zero
                        else
                            m = log10(sys.Params(pidx).MinVal);
                        end
                        if sys.Params(pidx).MaxVal == 0
                            M = -16;
                        else
                            M = log10(sys.Params(pidx).MaxVal);
                        end
                        ranges{pidx} = logspace(m, M, sys.Params(pidx).Desired);
                    end
                end
            end
            
            % Return first range if only one available.
            if sys.ParamCount == 1
                samples = ranges{1};
            else
                samples = general.Utils.createCombinations(ranges);
            end
        end
    end
    
end

