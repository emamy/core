classdef MatUtils
    %MATUTILS Matrix utility functions
    
    methods(Static)
        function [A, idxmat] = laplacemat(h, d1, d2)
            %Computes a 2D diffusion sparse matrix with zero neuman
            %boundary conditions.
            %
            % Arguments:
            % h - discretization spatial stepwidth
            % d1 - region dimension 1
            % d2 - region dimension 2
            %
            % @Daniel Wirtz, 11.03.2010
            
            m = d1*d2;
            idxmat = zeros(d1,d2);
            idxmat(:) = 1:m;
            
            inner = find(conv2(ones(size(idxmat)),[0 1 0;1 0 1; 0 1 0],'same')==4);
            
            i = []; j = []; s =[];
            %% Inner points
            % N E S W
            createStencil([ -1 d1 1 -d1],[-4 1 1 1 1], inner);
            
            %% Points that are on a boundaries
            % Left boundary neumann
            createStencil([-1 d1 1] ,[-3 1 1 1],(2:d1-1)');
            % Right boundary neumann
            createStencil([-1 -d1 1] ,[-3 1 1 1],(m-d1+2:m-1)');
            % Top boundary neumann
            createStencil([-d1 1 d1] ,[-3 1 1 1],(d1+1:d1:m-2*d1+1)');
            % Bottom boundary neumann
            createStencil([-d1 -1 d1] ,[-3 1 1 1],(2*d1:d1:m-d1)');
            
            % Edge points
            % Top left
            createStencil([1 d1] ,[-2 1 1],1);
            % Top right
            createStencil([1 -d1] ,[-2 1 1],m-d1+1);
            % bottom left
            createStencil([-1 d1] ,[-2 1 1],d1);
            % bottom right
            createStencil([-1 -d1] ,[-2 1 1],m);
            %% Compile matrix
            A = (1/h^2)*sparse(i,j,s,m,m);
            
            % Some comment on how createStencil works - ABOVE
            function createStencil ( stencil , weights , points )
                % Some comment on how createStencil works - BELOW
                %
                % Parameters:
                % stencil - Testdescription
                %
                % weights - bla bla
                
                i = [i; idxmat(points)];
                j = [j; idxmat(points)];
                s = [s; weights(1)* ones(size(points))];
                idx=2;
                for offset = stencil
                    i = [i; idxmat(points)];
                    j = [j; idxmat(points+offset)];
                    s = [s; weights(idx)*ones(size(points))];
                    idx=idx+1;
                end
            end
        end
        
        
    end
    
end

