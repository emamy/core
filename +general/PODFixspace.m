classdef PODFixspace < general.POD & general.Orthonormalizer
    %PODFIXSPACE Summary of this class goes here
    %   Detailed explanation goes here
      
    methods
        function onvec = computePODFixspace(this, vec, fixvec)
            
            % Orthonormalize fixed vectors
            fixvec = this.orthonormalize(fixvec);
            
            % Subtract projection onto fixed space fixvec
            vec = vec - fixvec * (fixvec' * this.G * vec);
            
            %Gr = vec'*this.G*vec;
            %Gr = 0.5*(Gr+Gr');
            %onvec = this.computePOD(Gr);
            
            onvec = this.computePOD(vec);
            
            onvec = this.orthonormalize(onvec);
        end
    end
    
end

