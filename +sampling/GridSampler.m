classdef GridSampler < sampling.BaseSampler
% GridSampler: Samples the models params on a grid, using the
% ModelParam.Desired field
%
% @author Daniel Wirtz @date 2010-04-01
%
% @change{0,7,dw,2013-09-05} Moved the Spacing property to
% data.ModelParam
%
% @change{0,5,dw,2011-11-09}
% - Fixed logspace grid sampling when any of the range values was zero.
% - New default value for Spacing is 'lin'.
%
% @change{0,5,dw,2011-10-19} Added a linear or logarithmic spacing method.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    methods
       
        function samples = performSampling(this, params)
            % Uses the given model and generates a training set by creating
            % a regular grid in joint time/parameter space
            
            nparams = length(params);
            ranges = cell(nparams,1);
            
            % Create linearly spaced parameter value ranges
            for pidx=1:nparams
                % Cater for single parameters 
                if (params(pidx).MinVal == params(pidx).MaxVal)
                    ranges{pidx} = params(pidx).MinVal;
                else
                    if strcmp(params(pidx).Spacing,'lin')
                        ranges{pidx} = linspace(params(pidx).MinVal,...
                        params(pidx).MaxVal,...
                        params(pidx).Desired);
                    else
                        % Avoid NaNs due to log10(0)
                        if params(pidx).MinVal == 0
                            m = -16; % take double precision zero
                        else
                            m = log10(params(pidx).MinVal);
                        end
                        if params(pidx).MaxVal == 0
                            M = -16;
                        else
                            M = log10(params(pidx).MaxVal);
                        end
                        ranges{pidx} = logspace(m, M, params(pidx).Desired);
                    end
                end
            end
            
            % Return first range if only one available.
            if nparams == 1
                samples = ranges{1};
            else
                samples = Utils.createCombinations(ranges);
            end
            
            % Restrict samples to given domain if set
            if ~isempty(this.Domain)
                samples = this.Domain.filter(samples);
            end
        end
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj)
            % Introduced as the "Spacing" property has moved to
            % data.ModelParam
            if ~isa(obj,'sampling.GridSampler')
                obj = sampling.GridSampler;
            end
            obj = loadobj@sampling.BaseSampler(obj);
        end
    end
    
end

