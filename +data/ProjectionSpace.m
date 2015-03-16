classdef ProjectionSpace < KerMorObject
% ProjectionSpace: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2015-03-16
%
% @new{0,7,dw,2015-03-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=private)
        V = [];
        W = [];
        
        % The dimensions of the full state space which are associated with
        % this subspace
        Dimensions;
        
        Size;
    end
    
    methods
        function this = ProjectionSpace(V, W, dims)
            this.V = V;
            if isempty(W)
                this.W = V;
            end
            this.Dimensions = dims;
            this.Size = size(V,2);
        end
        
        function set.V(this, value)
            if ~isa(value,'data.FileMatrix') && (~isa(value, 'double') || ~ismatrix(value))
                error('value must be a valid matrix of type double or a data.FileMatrix');
            end
            this.V = value;
        end
        
        function set.W(this, value)
            if ~isa(value,'data.FileMatrix') && (~isa(value, 'double') || ~ismatrix(value))
                error('value must be a valid matrix of type double or a data.FileMatrix');
            end
            this.W = value;
        end
        
        function delete(this)
            this.V = [];
            this.W = [];
        end
    end
    
end