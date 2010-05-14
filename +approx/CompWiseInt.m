classdef CompWiseInt < approx.BaseKernelApprox
    %COMPWISEINT Component wise kernel interpolation
    
    properties(Access=private)
        adata;
        xi;
    end
    
    methods
        function this = CompWiseInt
            this.CustomProjection = true;
        end
    end
    
    methods(Access=protected)
        
        function fx = evaluate_approximation(this, x)
            % Evaluates the approximated function at point x
            K = this.evaluateKernel(x,this.xi);
            fdims = size(this.adata,1);
            fx = zeros(fdims,size(x,2));
            for fdim = 1:fdims
                fx(fdim,:) = K * this.adata(fdim,:)';
            end
        end
        
        function gen_approximation_data(this, xi, fxi)
            % Computes the approximation according to the concrete
            % approximation strategy. 
            
            in = general.interpolation.KernelInterpol;
            in.K = this.evaluateKernel(xi);
            
            % Make hermetian (no rounding errors)
            %ls.K = .5*(ls.K'+ls.K);
            try
                wh = waitbar(0,'Initializing component-wise kernel LS');
                fdims = size(fxi,1);
                this.adata = zeros(fdims, size(xi,2));
                for fdim = 1:fdims
                    waitbar(fdim/fdims,wh,sprintf('Performing interpolation for dimension %d/%d ... %2.0f %%',fdim,fdims,(fdim/fdims)*100));
                    fx = fxi(fdim,:);
                    this.adata(fdim,:) = in.interpolate(fx);
                end
                close(wh);
            catch ME
                close(wh);
                rethrow(ME);
            end
            this.xi = xi;
        end
        
        function copy = customProject(this, V)
            copy = this.clone;
            
            copy.adata = V'*this.adata;
            
            if this.RotationInvariantKernel
                % Extract system part and project into V space
                [x,t,mu] = this.splitTripleVect(this.xi);
                x = V' * x;
                copy.xi = this.compileTripleVect(x,t,mu);
            end
        end
        
         function target = clone(this)
            % Makes a copy of this instance.
            %
            % See also: ICloneable
            
            target = approx.CompWiseInt;
            % Call superclass clone
            target = clone@approx.BaseKernelApprox(this, target);
            
            target.adata = this.adata;
            target.xi = this.xi;
        end
    end
    
    methods(Static)
        
    end
    
end

