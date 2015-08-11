classdef ConstantUntil < general.functions.AFunGen
    % A constant function of value 1 util the given time.
    % Uses a ramp to go down from 1 to 0 for continuity. The percent
    % parameter specifies the percentage of the overall nonzero time that
    % is used for the 1-0 transition.
    
    properties(Access=private)
        stoptime;
        rampperc;
    end
    
    methods
        function this = ConstantUntil(time, percent)
            if nargin < 2
                percent = .005;
                if nargin < 1
                    time = 1;
                end
            end
            this.stoptime = time;
            this.rampperc = percent;
        end
        
        function [fhandle, df] = getFunction(this)
            st = this.stoptime;
            p = this.rampperc;
            rampstart = st*(1-p);
                          % 1 for t < rampstart  
            fhandle = @(t)(t < rampstart) + ...
                (t >= rampstart & t < st).*(1-(t-st*(1-p))/(st*p)); 
                % 1 to 0 over ramptime, zero after
            df = @(t)0;
        end
        
        function str = getConfigStr(this)
            str = sprintf('Endtime: %g, Ramp percent: %g',this.stoptime,this.rampperc);
        end
        
        function plot(this, varargin)
            plot@general.functions.AFunGen(this, 'R', [0 this.stoptime*1.5], varargin{:});
        end
    end
    
    methods(Static)
        function res = test_ConstantUntil
            f = general.functions.ConstantUntil;
            f.plot;
            f = general.functions.ConstantUntil(5,.1);
            fh = f.getFunction;
            res = fh(4.5) == 1 && fh(4.75) == .5 && fh(5) == 0;
            f.plot;
        end
    end
    
end

