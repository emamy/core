classdef RotationDecorator < spacereduction.BaseSpaceReducer
% RotationDecorator: Decorator for any other space reducer which rotates the resulting matrix.
%
% This space-reduction class has been added in order to be able to obtain slightly worse subspace
% computations as would be given by the underlying space reduction method.
%
% @author Daniel Wirtz @date 2011-05-23
%
% @new{0,4,dw,2011-05-23} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        % `\in [0, 2\pi]`
        %
        % @default .05
        %
        % @propclass{data} Not for real reduction use but rather experiments with controlled
        % projection error.
        Degree = .05;
        
        % Number of dimensions to rotate by @ref Degree
        %
        % @propclass{data} Not for real reduction use but rather experiments with controlled
        % projection error.
        %
        % See also: Degree
        Dims = 5;
    end
    
    properties(Access=private)
        % The inner space reducer
        sp;
    end
    
    methods
        function this = RotationDecorator(s)
            % Creates a new rotation decorator
            % Parameters:
            % s: Subclass of spacereduction.BaseSpaceReducer
            if isempty(s)
                error('Rotation decorator must be given an underlying space reduction class.');
            end
            if ~strcmp(s.ReducableDims,':')
                error('Rotation decorator not implemented for partial reduction (ReducableDims != :)');
            end
            this.sp = s;
            this.registerProps('Degree','Dims');
        end
        
        function [V,W] = generateReducedSpaceImpl(this, model)
            % Computes the subspace given by the underlying subspace reduction class but then
            % rotates for Dims times between two randomly chosen axis
            %
            
            % Call subclass reduction
            [V, W] = this.sp.generateReducedSpace(model);
            
            %rnd = RandStream('mt19937ar','Seed',2564);
            
            n = size(V,1); %#ok<*PROP>
            if this.Dims >= n/2
                warning('KerMor:spacereduction:RotationDecorator','Too large Dims property (%d) for full space dimension %d, setting to %d',this.Dims,n,floor(n/2));
                this.Dims = floor(n/2);                
            end
            R = spdiags(ones(n,1),0,n,n);
            for idx = 1:this.Dims
                Q = spdiags(ones(n,1),0,n,n);
                %idx1 = rnd.randi(n);
                idx1 = 2*(idx-1)+1;
                %idx2 = rnd.randi(n);
                idx2 = idx1+1;
                
                Q(idx1,idx1) = cos(this.Degree);
                Q(idx1,idx2) = -sin(this.Degree);
                Q(idx2,idx1) = sin(this.Degree);
                Q(idx2,idx2) = cos(this.Degree);
                R = R*Q;
            end
            V = R*V;
            W = R*W;
        end
    end
    
end