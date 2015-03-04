classdef Krylov < spacereduction.BaseSpaceReducer
    %KRYLOV Krylov Subspace generation
    %   @todo Implementation
    
    properties
        mu0;
    end
    
    methods(Access=protected)
        function [V,W] = generateReducedSpaceImpl(this, model)
            
            x0 = model.System.x0(this.mu0);
            
            Ax = model.System.f.evaluate(x0,0,this.mu0);
            
        end
    end
    
end

