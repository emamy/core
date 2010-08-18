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
    
    methods(Access=protected, Sealed)
        
        function prepareApproximationGeneration(this, K)
            this.KI = general.interpolation.KernelInterpol;
            % Assign Kernel matrix
            this.KI.K = K;
        end
        
        function [ai, b, svidx] = calcComponentApproximation(this, fxi)
            ai = this.KI.interpolate(fxi);
            b = 0;
            svidx = [];
        end
        
        function copy = clone(this)
            copy = approx.CompWiseInt;
            copy = clone@approx.BaseCompWiseKernelApprox(this, copy);
        end
    end
end

