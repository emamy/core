classdef Sinus < general.functions.AFunGen
    
    properties(Access=private)
        freq;
        opts;
    end
    
    methods
        function this = Sinus(varargin)
            i = inputParser;
            i.addParamValue('Frequency',1,@(v)isscalar(v));
            i.addParamValue('fOff',0,@(v)isscalar(v));
            i.addParamValue('vOff',0,@(v)isscalar(v));
            i.addParamValue('Amplitude',1,@(v)isscalar(v));
            i.parse(varargin{:});
            this.opts = i.Results;
        end
        
        function [fhandle, dfhandle] = getFunction(this)
            o = this.opts;
            fhandle = @(t)o.Amplitude*sin(t/1000*o.Frequency*2*pi+o.fOff)+o.vOff;
            dfhandle = @(t)o.Amplitude*cos(t/1000*o.Frequency*2*pi+o.fOff);
        end
        
        function str = getConfigStr(this)
            o = this.opts;
            str = sprintf('Frequency: %g [Hz], Amplitude: %g [-], ArgOffset: %g, ValueOffset: %g',...
                o.Frequency,o.Amplitude,o.fOff,o.vOff);
        end
    end
    
    methods(Static)
        function res = test_Sinus
            f = general.functions.Sinus;
            f.plot;
            f = general.functions.Sinus('Frequency',50,...
                'fOff',2,'vOff',4,'Amplitude',2);
            f.plot;
            res = true;
        end
    end
    
end

