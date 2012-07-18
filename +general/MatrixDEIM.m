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
        % The custom (partial) similarity transform matrix `Q_k` applied to
        % the approximated matrices via `Q_k^tM(x,t,\mu)Q_k`.
        %
        % The matrix has to be orthonormal if set.
        %
        % @type matrix<double> @default []
        Qk = [];
    end
    
    properties(Access=private)
        effNumRows;
    end
    
    methods
        function M = evaluate(this, x, t, mu)
            fx = evaluate@general.DEIM(this, x, t, mu);
            M = reshape(fx, this.effNumRows, []);
        end
        
        function setSimilarityTransform(this, Qk)
            this.Qk = Qk;
            U = this.U_nonproj;
            if isempty(Qk) || isequal(Qk,1)
                this.effNumRows = this.NumRows;
                unew = U;
            else
                this.effNumRows = size(Qk,2);
                unew = zeros(size(this.Qk,2)^2, size(U,2));
                for i=1:size(U,2)
                    umat = reshape(U(:,i),size(this.Qk,1),[]);
                    umat = this.Qk'*(umat*this.Qk);
                    unew(:,i) = umat(:);
                end    
            end
            % Overwrite the current U with processed U_nonproj.
            this.U = unew;
        end
        
        function target = project(this, V, ~)
            target = project@general.DEIM(this, V, [], this.clone);
        end
        
        function copy = clone(this)
            copy = general.MatrixDEIM;
            copy = clone@general.DEIM(this, copy);
            copy.Qk = this.Qk;
            copy.NumRows = this.NumRows;
            copy.effNumRows = this.effNumRows;
        end
    end
    
    methods(Access=protected)
        function updateOrderData(this)
            % @todo make better split-up of updateOrderData (w.r.t. error
            % matrices)
            updateOrderData@general.DEIM(this);
            
            % "special" for jaccompevalwrappers:
            % they only evaluate to nonzero entries, which is why the DEIM approximation
            % U, U_nonproj is of size length(J(J ~= 0)) instead of length(J(:))
            if isa(this.f,'general.JacCompEvalWrapper') && ~isempty(this.f.f.JSparsityPattern)
                fsp = this.f.f.JSparsityPattern;
                this.U_nonproj = general.MatUtils.toSparse(this.U_nonproj, find(fsp), numel(fsp)); %#ok<FNDSB>
            end
            
            % Apply partial similarity transform (AFTER original DEIM order
            % update)
            % Checks for having Qk set are done there.
            this.setSimilarityTransform(this.Qk);
            
            % For now: unset the auxiliary DEIM error estimation matrices
            % as they are not used (cannot be used?) for matrix DEIMs yet.
            this.M1 = [];
            this.M2 = [];
%             this.Uerr1 = [];
%             this.Uerr2 = [];
        end
    end
end