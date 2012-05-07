classdef ACompEvalCoreFun < dscomponents.ACoreFun
% ACompEvalCoreFun: A normal CoreFun which supports single-component
% evaluation.
%
% This extra abilities are used within the approx.DEIM approximation.
%
% @author Daniel Wirtz @date 2012-03-26
%
% @new{0,6,dw,2012-03-26} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Dependent)
        PointSets;
    end
    
    properties(Access=private)
        pts = {};
        
        % The x-component selection matrices (precomputed on setting PointSet/AltPointSet)
        S = {};
        
        % Notational vectors for EI-point function argument dependencies
        jrow = {};
        jend = {};
        jself = {};
    end

    methods
        function fx = evaluateComponentSet(this, nr, x, t, mu)
            fx = this.evaluateComponents(this.PointSets{nr},...
                this.jend{nr}, this.jrow{nr}, this.jself{nr}, this.S{nr}*x, t, mu);
        end
        
        function setPointSet(this, nr, pts)
            this.pts{nr} = pts;
            
            jr = []; js = logical.empty(0,1);
            je = zeros(1,length(pts));
            J = this.JSparsityPattern';
            for i=1:length(pts)
                [inew, ~] = find(J(:,pts(i)));
                jr = [jr inew'];%#ok
                js = [js inew' == pts(i)];%#ok
                je(i) = length(jr);
            end
            this.jrow{nr} = jr;
            this.jself{nr} = js;
            this.jend{nr} = je;
            
            % Compose x-entry selection matrix
            len = je(end);
            sel = jr(1:len);
            if ~isempty(this.V)
                this.S{nr} = this.V(sel,:);
            else
                this.S{nr} = sparse(1:len,sel,ones(len,1),len,this.XDim);
            end
        end
        
        function target = project(this, V, W, target)
            if nargin < 4
                target = this.clone;
            end
            target = project@dscomponents.ACoreFun(this, V, W, target);
            for i=1:length(this.S)
                target.S{i} = this.S{i}*V;
            end
        end
        
        function copy = clone(this, copy)
            copy = clone@dscomponents.ACoreFun(this, copy);
            copy.pts = this.pts;
            copy.S = this.S;
            copy.jrow = this.jrow;
            copy.jend = this.jend;
            copy.jself = this.jself;
        end
        
        function psets = get.PointSets(this)
            psets = this.pts;
        end
        
        function res = test_IComponentEvaluableMatch(this, dim, pdim)
            % Tests if the local implementation of
            x = rand(dim,1);
            mu = rand(pdim,1);
            t = rand;
            fx = this.evaluateCoreFun(x, t, mu);
            oldpts = this.PointSets{1};
            this.setPointSet(1,1:dim);
            fxc = this.evaluateComponentSet(1, x, t, mu);
            this.setPointSet(1,oldpts);
            d = norm(fx-fxc);
            fprintf('Norm difference: %e\n',d);
            res = d < sqrt(eps);
        end
    end
    
    methods(Abstract)
        evaluateComponents(this, pts, ends, idx, self, x, t, mu);
    end
end