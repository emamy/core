classdef CompWiseInt < approx.BaseCompWiseKernelApprox
    %COMPWISEINT Component wise kernel interpolation
    %
    
    %| @docupdate 
    
    properties(Access=private)
        % The single dimension kernel interpolation algorithm
        %
        % @type general.interpolation.KernelInterpol
        %
        % See also: general.interpolation.KernelInterpol
        KI;
    end
    
    methods(Sealed)
        function copy = clone(this)
            copy = approx.CompWiseInt;
            copy = clone@approx.BaseCompWiseKernelApprox(this, copy);
        end
    end
    
    methods(Access=protected, Sealed)
        
        function prepareApproximationGeneration(this, K)
            ki = general.interpolation.KernelInterpol;
            ki.UseLU = true;
            ki.K = K;
            this.KI = ki;
        end
        
        function [ai, b, svidx] = calcComponentApproximation(this, fxi)
            [ai,b] = this.KI.interpolate(fxi);
            svidx = [];
        end
    end
end

