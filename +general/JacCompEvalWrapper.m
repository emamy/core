classdef JacCompEvalWrapper < dscomponents.ACompEvalCoreFun
% JacCompEvalWrapper: Wraps the evaluation of a ACompEvalCoreFun's jacobian
% into a vectorized function.
%
% This class is used to enable matrix DEIM approximation.
%
% @author Daniel Wirtz @date 2012-06-01
%
% @change{0,6,dw,2012-07-16} Added fast implementation of evaluate for already sparse
% jacobians.
%
% @change{0,6,dw,2012-06-08} Improved the speed of the evaluate-function by
% separate code for sparse/full matrices.
%
% @new{0,6,dw,2012-06-01} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        f;
        trafo = {};
    end
    
    properties(SetAccess=private)
        fullxDim;
    end
    
    methods
        function this = JacCompEvalWrapper(f)
          % Creates a new jacobian matrix wrapper for linear indexing,
          % which implements dscomponents.ACompEvalCoreFun
          
          % Call superclass constructor
          this = this@dscomponents.ACompEvalCoreFun;
          
          if nargin == 1
              % Checks
              if ~isa(f,'dscomponents.ACompEvalCoreFun')
                  error('The system''s nonlinearity must be a component-wise evaluable function.');
              end
              % Original full system function
              this.f = f;

              % ACoreFun settings
              this.MultiArgumentEvaluations = true;
              this.xDim = f.xDim;
              % output dim is size of jacobian
              if ~isempty(f.JSparsityPattern)
                  fdim = length(find(f.JSparsityPattern));
              else
                  fdim = f.fDim * f.xDim;
              end
              this.fDim = fdim;
              % So sparsity pattern for the wrapper
              this.JSparsityPattern = [];

              % Extra: Store original full dimension, as f gets possibly
              % projected
              this.fullxDim = f.xDim;
          end
        end
        
        function fx = evaluate(this, x, t, mu)
            % Evaluates the jacobian of the wrapped function and returns a
            % vector corresponding to the natural linear indexing of
            % matlab.
            %
            % Return values:
            % fx: The linearly indexed jacobian of `f`, so `fx =
            % (\d{f_1}{x_1}, \d{f_2}{x_1}, \ldots, \d{f_n}{x_1},
            % \d{f_2}{x_1}, \ldots)`. @type colvec<double>
            
            if ~isempty(this.f.JSparsityPattern)
                nonzero = find(this.f.JSparsityPattern);
            else
                nonzero = 1:this.f.fDim*this.f.xDim;
            end
            fx = zeros(length(nonzero),size(x,2));
            for i=1:size(x,2)
                J = this.f.getStateJacobian(x(:,i),t(i),mu(:,i));
                fx(:,i) = J(nonzero);
            end
        end
        
        function evaluateCoreFun(varargin)
            % Nothing to do here as the evaluate function is overridden
            % directly.
            error('Boo. Do not call me. (directly implemented evaluate)');
        end
        
        function setPointSet(this, nr, pts)
            % Overrides the setPointSet method of
            % dscomponents.ACompEvalCoreFun to allow transformation of the
            % indices from linear to matrix components and derivatives.
            
            % Transform indices to jacobian matrix indices
            % (the jac wrapper only "works" on nonzero jacobian entries)
            if ~isempty(this.f.JSparsityPattern)
                jidx = find(this.f.JSparsityPattern);
                pts = jidx(pts)';
            end
            
            % Ensure row vectors
            if size(pts,1) > 1
                pts = reshape(pts,1,[]);
            end
            if length(unique(pts)) ~= length(pts)
                error('Points have to be unique.');
            end
            
            % Get matrix indexing of the desired points
            % ind2sub direct replacement!
            i = rem(pts-1, this.f.fDim)+1;
            j = (pts-i)/this.f.fDim+1;
            % Only use the i-th effective components of f, the j's
            % correspond to the derivatives of the i-th components with
            % respect to the j-th variable.
            ui = unique(i);
            
            deriv = {};
            pos = [];
            for k = 1:length(ui)
                didx = find(i == ui(k));
                deriv{k} = j(didx); %#ok<*AGROW>
                pos = [pos didx];
            end
            % Can directly set point sets as f instance is clone of
            % original function
            this.f.setPointSet(nr+2, ui, deriv);
            
            l = length(pos);
            this.trafo{nr} = sparse(pos, 1:l, ones(l,1),l,l);
        end
        
        function fx = evaluateComponentSet(this, nr, x, t, mu)
            dfx = this.f.evaluateJacobianSet(nr+2, x, t, mu);
            fx = this.trafo{nr}*dfx;
        end
        
        function copy = clone(this)
            copy = general.JacCompEvalWrapper;
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            copy.f = this.f.clone;
            copy.trafo = this.trafo;
            copy.fullxDim = this.fullxDim;
        end
        
        function proj = project(this, V, W)
            proj = project@dscomponents.ACompEvalCoreFun(this, V, W, this.clone);
            % Project the wrapped function, too (ignores the extra copy.f
            % created in clone as project creates a clone anyways)
            proj.f = this.f.project(V, W);
            proj.xDim = proj.f.xDim;
            % vec-op leads to product of jacobian size dimension
            proj.fDim = proj.f.fDim*proj.f.xDim;
        end
    end
    
    methods(Access=protected)
        function varargout = evaluateComponents(varargin)
            % nothing to do here as this is a wrapper
        end
    end
end