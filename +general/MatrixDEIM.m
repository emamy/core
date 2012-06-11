classdef MatrixDEIM < general.DEIM
% MatrixDEIM: 
%
%
%
% @author Daniel Wirtz @date 2012-06-04
%
% @new{0,6,dw,2012-06-04} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The number of rows of the output matrices (information needed for
        % reshape command of internal column-wise DEIM)
        NumRows;
    end
    
    properties(SetAccess=private)
        % The custom (partial) similarity transform matrix `S` applied to
        % the approximated matrices via `S^tM(x,t,\mu)S`.
        %
        % The matrix has to be orthonormal if set.
        %
        % @type matrix<double> @default []
        ST = [];
    end
    
    properties(Access=private)
        effNumRows;
    end
    
    methods
        function M = evaluate(this, x, t, mu)
            fx = evaluate@general.DEIM(this, x, t, mu);
            M = reshape(fx, this.effNumRows, []);
        end
        
        function setSimilarityTransform(this, ST)
            this.ST = ST;
            U = this.U_nonproj;
            if isempty(ST) || isequal(ST,1)
                this.effNumRows = this.NumRows;
                unew = U;
            else
                this.effNumRows = size(ST,2);
                unew = zeros(size(this.ST,2)^2, size(U,2));
                for i=1:size(U,2)
                    umat = reshape(U(:,i),size(this.ST,1),[]);
                    umat = this.ST'*(umat*this.ST);
                    unew(:,i) = umat(:);
                end    
            end
            % Overwrite the current U with processed U_nonproj.
            this.U = unew;
        end
        
        function target = project(this, V, ~)
            target = this.clone;
            target = project@general.DEIM(this, V, [], target);
        end
        
        function copy = clone(this)
            copy = general.MatrixDEIM;
            copy = clone@general.DEIM(this, copy);
            copy.ST = this.ST;
            copy.NumRows = this.NumRows;
            copy.effNumRows = this.effNumRows;
        end
    end
    
    methods(Access=protected)
        function updateOrderData(this)
            % @todo make better split-up of updateOrderData (w.r.t. error
            % matrices)
            updateOrderData@general.DEIM(this);
            
            % Apply partial similarity transform (AFTER original DEIM order
            % update)
            % Checks for having ST set are done there.
            this.setSimilarityTransform(this.ST);
            
            % For now: unset the auxiliary DEIM error estimation matrices
            % as they are not used (cannot be used?) for matrix DEIMs yet.
            this.M1 = [];
            this.M2 = [];
%             this.Uerr1 = [];
%             this.Uerr2 = [];
        end
    end
end