classdef RectGrid3D < handle
% Rect3D: Rectangular three-dimensional grid
%
% Given the number of elements in each direction, this class contains
% information (=indices) about the different geometric objects like sides, edges etc.
%
% @author Daniel Wirtz @date 2011-10-06
%
% @new{0,5,dw,2011-10-06} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess = private)
        % The indices of the rectangular interior points
        Inner;
        
        % A struct with fields F, Ba, L, R, T, Bo denoting the indices of
        % the interior points of the respective rectangle sides.
        Sides;
        
        % A struct with the fields
        Corners;
        
        % A struct with the fields 
        Edges;
        
        % A 1x3 row vector containing the dimensions in x,y and z
        % direction.
        Dims;
    end
    
    properties(Dependent)
        % The total number of points in the grid.
        Points;
        
        % A matrix of size dim1 x dim2 x dim3 containing the linear indices
        % of the entry at each entry.
        IndexMatrix;
        
        % The front side interior points indices
        F;
        
        % The back side interior point indices
        Ba;
        
        % The left side interior point indices
        L;
        
        % The right side interior point indices
        R;
        
        % The top side interior point indices
        T;
        
        % The bottom side interior point indices
        Bo;
    end
    
    methods
        function this = RectGrid3D(dim1, dim2, dim3)
            this.Dims = [dim1 dim2 dim3];
            m = dim1;
            k = dim2;
            L = this.IndexMatrix;
            
            %% Inner points
            % Find inner points
            this.Inner = find(convn(ones(m,k,dim3),cat(3,[0 0 0; 0 1 0; 0 0 0],[0 1 0;1 0 1; 0 1 0],[0 0 0; 0 1 0; 0 0 0]),'same') == 6);
            
            this.Sides = struct;
            % Front boundary neumann
            this.Sides.F = reshape(L(2:end-1,2:end-1,1),[],1);
            % Back boundary neumann
            this.Sides.Ba = reshape(L(2:end-1,2:end-1,end),[],1);
            % Left boundary neumann
            this.Sides.L = reshape(L(2:end-1,1,2:end-1),[],1);
            % Right boundary neumann
            this.Sides.R = reshape(L(2:end-1,end,2:end-1),[],1);
            % Top boundary neumann
            this.Sides.T = reshape(L(1,2:end-1,2:end-1),[],1);
            % Bottom boundary neumann
            this.Sides.Bo = reshape(L(end,2:end-1,2:end-1),[],1);

            % Corner points
            this.Corners = struct;
            %% Front side corners
            % Front top
            this.Corners.FT = reshape(L(1,2:end-1,1),[],1);
            % Front bottom
            this.Corners.FBo = reshape(L(end,2:end-1,1),[],1);
            % Front left
            this.Corners.FL = reshape(L(2:end-1,1,1),[],1);
            % Front right
            this.Corners.FR = reshape(L(2:end-1,end,1),[],1);
            
            %% Back side corners
            % Back top
            this.Corners.BaT = reshape(L(1,2:end-1,end),[],1);
            % Back bottom
            this.Corners.BaBo = reshape(L(end,2:end-1,end),[],1);
            % Back left
            this.Corners.BaL = reshape(L(2:end-1,1,end),[],1);
            % Back right
            this.Corners.BaR = reshape(L(2:end-1,end,end),[],1);
            
            %% front-back connecting corners
            % top left
            this.Corners.TL = reshape(L(1,1,2:end-1),[],1);
            % top right
            this.Corners.TR = reshape(L(1,end,2:end-1),[],1);
            % bottom left
            this.Corners.BoL = reshape(L(end,1,2:end-1),[],1);
            % bottom right
            this.Corners.BoR = reshape(L(end,end,2:end-1),[],1);

            %% Edge points
            this.Edges = struct;
            % Front Top left
            this.Edges.FTL = 1;
            % Front Top right
            this.Edges.FTR = L(1,end,1);
            % Front bottom left
            this.Edges.FBoL = L(end,1,1);
            % Front bottom right
            this.Edges.FBoR = L(end,end,1);
            % Rear Top left
            this.Edges.BaTL = L(1,1,end);
            % Rear Top right
            this.Edges.BaTR = L(1,end,end);
            % Rear bottom left
            this.Edges.BaBoL = L(end,1,end);
            % Rear bottom right
            this.Edges.BaBoR = L(end,end,end);
            
            clear m k L;
        end
    end
    
    %% Getter
    methods
        function value = get.Points(this)
            value = prod(this.Dims);
        end
        
        function value = get.IndexMatrix(this)            
            value = zeros(this.Dims);
            value(:) = 1:this.Points;
        end
        
        function value = get.F(this)
            value = this.Sides.F;
        end
        
        function value = get.Ba(this)
            value = this.Sides.BA;
        end
        
        function value = get.T(this)
            value = this.Sides.T;
        end
        
        function value = get.Bo(this)
            value = this.Sides.Bo;
        end
        
        function value = get.L(this)
            value = this.Sides.L;
        end
        
        function value = get.R(this)
            value = this.Sides.R;
        end
    end
    
end