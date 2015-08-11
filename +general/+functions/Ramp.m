classdef Ramp < general.functions.AFunGen
    
    properties(Access=private)
        ramptime;
        max;
        starttime;
    end
    
    methods
        function this = Ramp(ramptime, max, starttime)
            if nargin < 3
                starttime = 0;
                if nargin < 2
                    max = 1;
                    if nargin < 1
                        ramptime = 1;
                    end
                end
            end
            this.ramptime = ramptime;
            this.max = max;
            this.starttime = starttime;            
        end
        
        function [fhandle, df] = getFunction(this)
            rt = this.ramptime;
            if rt <= 0
                fhandle = @(t)0;
            else
                maxval = this.max;
                start = this.starttime;
                fhandle = @(t)(t >= start) .* (maxval * (((t-start)<rt).*(t-start)/rt + (t>=rt+start)));
            end
            df = [];
        end
        
        function plot(this, varargin)
            varargin = [{'R', [this.starttime*.75 this.starttime+this.ramptime*1.33]} varargin];
            plot@general.functions.AFunGen(this, varargin{:});
        end
        
        function str = getConfigStr(this)
            str = sprintf('Start: %g, Ramptime: %g, Max:%g',this.starttime,this.ramptime,this.max);
        end
    end
    
    methods(Static)
        function res = test_Ramp
            f = general.functions.Ramp;
            f.plot;
            f = general.functions.Ramp(10,2,4);
            fh = f.getFunction;
            res = fh(2) == 0 && fh(4) == 0 && fh(9) == 1 && fh(14) == 2;
            f.plot;
        end
    end
    
end

