classdef JacCompEvalWrapper < dscomponents.ACompEvalCoreFun
% JacCompEvalWrapper: 
%
%
%
% @author Daniel Wirtz @date 2012-06-01
%
% @new{0,6,dw,2012-06-01} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        f;
        trafo = {};
    end
    
    properties(SetAccess=private)
        fullXDim;
    end
    
    methods
        function this = JacCompEvalWrapper(f)
          % Creates a new jacobian matrix wrapper for linear indexing,
          % which implements dscomponents.ACompEvalCoreFun
          % Checks
          if ~isa(f,'dscomponents.ACompEvalCoreFun')
              error('The system''s nonlinearity must be a component wise evaluable function.');
          end
          
          % Call superclass constructor
          this = this@dscomponents.ACompEvalCoreFun;
          
          % Get copy of original full system function
          this.f = f.clone;
          
          % ACoreFun settings
          this.MultiArgumentEvaluations = true;
          this.XDim = f.XDim;
          this.JSparsityPattern = [];
          
          % Extra: Store original full dimension, as f gets possibly
          % projected
          this.fullXDim = f.XDim;
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
            fx = this.f.getStateJacobian(x(:,1),t(1),mu(:,1));
            fx = fx(:);
            if size(x,2) > 1
                for i=2:size(x,2)
                    J = this.f.getStateJacobian(x(:,i),t(i),mu(:,i));
                    fx(:,i) = J(:);
                end
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
            
            % Ensure row vectors
            if size(pts,1) > 1
                pts = reshape(pts,1,[]);
            end
            if length(unique(pts)) ~= length(pts)
                error('Points have to be unique.');
            end
            
            % Get matrix indexing of the desired points 
            [i, j] = ind2sub([this.fullXDim this.fullXDim], pts);
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
            this.f.setPointSet(nr, ui, deriv);
            
            l = length(pos);
            this.trafo{nr} = sparse(pos, 1:l, ones(l,1),l,l);
        end
        
        function fx = evaluateComponentSet(this, nr, x, t, mu)
            dfx = this.f.evaluateJacobianSet(nr, x, t, mu);
            fx = this.trafo{nr}*dfx;
        end
        
        function varargout = evaluateComponents(varargin)
                % nothing to do here as this is a wrapper
        end
        
        function copy = clone(this)
            copy = general.JacCompEvalWrapper(this.f);
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            copy.f = this.f;
            copy.trafo = this.trafo;
            copy.fullXDim = this.fullXDim;
        end
        
        function proj = project(this, V, W)
            proj = this.clone;
            proj = project@dscomponents.ACompEvalCoreFun(this, V, W, proj);
            proj.f = this.f.project(V,W);
        end
    end
end