classdef RandomModelEstimatorDemo < EstimatorDemo
    % Demo class for the error estimators.
    % Creates a random model using a kernel expansion that can be
    % configured with the provided properties.
    %
    % @author Daniel Wirtz @date 25.11.2010
    %
    % @change{0,2,dw,2011-03-22} Fixed the constructor to only take an
    % optional dimension argument.
    
    properties
        % System dimension
        Dims = 500;
        
        % Number of centers used in kernel expansion
        NumCenters = 10;
        
        % Strictly positive kernel expansion?
        PositiveExpansion = false;
        
        % Uniform expansion or comp-wise separate?
        UniformExpansion = false;
    end
    
    methods
        
        function this = RandomModelEstimatorDemo(dims)
            % Creates a new estimator demo.
            %
            % Parameters:
            % dims: A scalar defining the dimensions of the default test
            % model. Calling the constructor with no arguments causes the
            % Demo to use the default demo model with the default dimension
            % size (500).
            
            %% Model settings
            fm = models.BaseFullModel;
            fm.Name = 'Estimator Demo Model';
            
            fm.T = 1;
            fm.dt = .025;
            
            fm.Approx = [];
            fm.Sampler = [];
            
            %this.ODESolver = solvers.ode.MLWrapper(@ode45);
            fm.ODESolver = solvers.ode.ExplEuler(fm.dt);
            %fm.ODESolver = solvers.ode.Heun(fm.dt);
            
            %% Core function
            cf = dscomponents.ParamTimeKernelCoreFun;
            cf.TimeKernel = kernels.NoKernel;
            cf.ParamKernel = kernels.NoKernel;
            cf.Centers.ti = [];
            cf.Centers.mui = [];
            
            %% System settings
            fm.System = models.BaseDynSystem(fm);
            fm.System.f = cf;
            this.Model = fm;
            
            if nargin == 1
                this.Dims = dims;
            else
                this.newCoeffs;
            end
        end
        
        function setup(this)
            this.Model.System.f.Kernel = kernels.GaussKernel(15);
            x0 = rand(this.Dims,1);
            if this.PositiveExpansion
                base = linspace(0, 40, this.NumCenters);
                this.Model.System.x0 = dscomponents.ConstInitialValue(x0);
            else
                base = linspace(-20, 20, this.NumCenters);
                this.Model.System.x0 = dscomponents.ConstInitialValue(x0-.5);
            end
            this.Model.System.f.Centers.xi = repmat(base,this.Dims,1);
            
            if  this.UniformExpansion
                V = ones(this.Dims,1)*sqrt(1/this.Dims);
                s = spacereduction.ManualReduction(V,V);
            else
                s = spacereduction.PODReducer;
                s.Mode = 'abs';
                s.Value = 1;
            end
            this.Model.SpaceReducer = s;
            
            %% Generation
            this.setModel(this.Model);
        end
        
        function newCoeffs(this)
            % Function coefficients
            offset = .5;
            if this.PositiveExpansion
                offset = 0;
            end
            % Create coefficients
            if this.UniformExpansion
                ai = (rand(1,this.NumCenters)-offset);
                this.Model.System.f.Ma = repmat(ai,this.Dims,1);
            else
                this.Model.System.f.Ma = (rand(this.Dims,this.NumCenters)-offset);
            end
            this.setup;
        end
        
        function set.Dims(this, value)
            this.Dims = value;
            this.newCoeffs; %#ok
            %this.setup; %#ok
        end
        
        function set.NumCenters(this, value)
            this.NumCenters = value;
            this.newCoeffs; %#ok
            %this.setup; %#ok
        end
        
        function set.PositiveExpansion(this, value)
            this.PositiveExpansion = value;
            %this.setup; %#ok
            this.newCoeffs; %#ok
        end
        
        function set.UniformExpansion(this, value)
            this.UniformExpansion = value;
            this.newCoeffs; %#ok
            %this.setup; %#ok
        end
        
    end
end

