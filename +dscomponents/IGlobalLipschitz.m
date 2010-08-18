classdef IGlobalLipschitz < handle
    %IGLOBALLIPSCHITZ Interface for all functions that have a global
    % lipschitz constant for the state/spatial part.
    
    methods(Abstract)
        % For some error estimators, a global Lipschitz constant estimation
        % is performed. This function has to yield the global lipschitz
        % constant for the spatial/state part of a CoreFunction `f(x,t,mu)`
        %
        % @todo: what about cases where no lipschitz constants are
        % available? somehow need to tell the error estimator etc.
        c = getGlobalLipschitz(this, t, mu);
    end
    
end

