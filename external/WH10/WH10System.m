classdef WH10System < models.BaseDynSystem & dscomponents.ParamTimeKernelCoreFun
    % Numerical experiments class for Paper WH10
    %
    % Current version works with KerMor 0.4
    
    properties(SetObservable) 
        % SV's
        %
        % @propclass{experimental}
        svNum;
    end
    
    methods
        
        function this = WH10System(model)
            
            this = this@models.BaseDynSystem(model);
            
            this.registerProps('svNum');
            
            %% System settings
            dims = model.dim;

            this.f = this;
            this.MaxTimestep = [];
            
            % Sample bases
            this.svNum = 20;
            space = linspace(0,50,this.svNum);
            this.Centers.xi = repmat(space,dims,1);
            this.Centers.ti = [];
            this.Centers.mui = [];
                        
            %% BaseCompWiseKernelApprox settings
            % Choose the kernel with such that each kernel vanishes (<
            % kerneleps) after kernelsupport support vectors. 
            kernelsupport = 2;
            kerneleps = 0.00001;
            
            d = sqrt(dims)*space(2);
            %gamma = -((kernelsupport*d)^2)/log(kerneleps) preprint
            gamma = -((kernelsupport*d))/log(kerneleps)
            this.Kernel = kernels.GaussKernel(gamma);
            this.TimeKernel = kernels.NoKernel;
            this.ParamKernel = kernels.NoKernel;
        end
    end    
end

