classdef Base < handle
% Base: 
%
%
%
% @author Daniel Wirtz @date 2011-07-04
%
% @new{0,5,dw,2011-07-04} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(SetAccess=private, GetAccess=protected)
        M1;
    end
    
    methods
        function this = Base(rm)
            % Obtain the correct snapshots
            % Standard case: the approx function is a kernel expansion. it
            % can also be that the system's core function is already a
            % kernel expansion
            fm = rm.FullModel;
            if ~isempty(fm.Approx)
                % Get full d x N coeff matrix of approx function
                Ma = fm.Approx.Expansion.Ma;
            else
                % Get full d x N coeff matrix of core function
                Ma = fm.System.f.Expansion.Ma;
            end
            
            % Perform any offline computations/preparations
            % Only prepare matrices if projection is used
            if ~isempty(rm.V)
                % Compute projection part matrices, without creating a
                % d x d matrix (too big!)
                M = Ma - rm.V*(rm.W'*Ma);
                hlp = M'*(rm.G*M);
                % Check if matrix needs to be made symmetric
                if any(any(abs(hlp-hlp') > 1e-5))
                    hlp = (hlp + hlp')/2;
                    warning('KerMor:errorest','M1 matrix not sufficiently symmetric, updating (M+M'')/2');
                end
                this.M1 = hlp;
                
                this.inputOfflineComputations(rm, M);
                clear M;
            else
                % No projection means no projection error!
                n = size(Ma,2);
                this.M1 = zeros(n,n);
            end
        end
    end

    methods(Abstract)
        % Computes the `\alpha(x, t,\mu)` value of the error estimator.
        %
        % Template method.
        %
        % Parameters:
        % phi: The kernel vector `\K(x,x_i)` @type colvec
        % ut: The evaluation of the current input `u(t)` at time `t`
        % @type double
        % t: The current time `t` @type double
        % mu: The current parameter `\mu` @type colvec
        a = getAlpha(this, phi, ut, t, mu);
        
        % Performs the offline stage for the error estimators regarding the inputs.
        %
        % Template method.
        %
        % Parameters:
        % model: The full model @type models.BaseFullModel
        % M: The projected coefficient matrix `M_{\alpha} - VW^tM_{\alpha})` @type matrix
        inputOfflineComputations(this, model, M);
    end
    
end