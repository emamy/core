classdef BaseModel < handle
    %MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % The actual dynamical system used in the model.
        System = models.BaseDynSystem;
        
        % The name of the Model
        Name = 'Base Model';
        
        % The verbose output level at simulations
        Verbose = 0;
        
        % The final timestep up to which to simulate
        T = 1;
        
        % The desired time-stepsize for simulations
        dt = .1;
        
        % The solver to use for the ODE
        ODESolver = @ode45;
    end
    
    properties(Dependent)
        % Evaluation points of the model as increasing array
        Times;
    end
    
    methods
        
        function [t,y] = simulate(this, mu, inputidx)
            % Simulates the system and produces the system's output.
            %
            % Parameters:
            % mu: The concrete mu parameter sample to simulate for.
            % inputidx: The index of the input function to use.
            %
            % Both parameters are optional. (Which to provide will be
            % determined by the actual system anyways)
            %
            % Return values:
            % t: The times at which the model was evaluated
            % y: Depending on the existance of an output converter, this
            %    either returns the full trajectory or the processed output
            %    at times t.
            
            if nargin < 3
                inputidx = [];
                if nargin < 2
                    mu = [];
                end
            end
            
            [t,y] = this.computeTrajectory(mu, inputidx);
            if ~isempty(this.System.C)
                if this.System.C.TimeDependent
                    % Evaluate the output conversion at each time t
                    for idx=1:length(t)
                        y(:,idx) = this.System.C.evaluate(t(idx),mu)*y(:,idx);
                    end
                else
                    % otherwise it's a constant matrix so multiplication
                    % can be preformed much faster.
                    y = this.System.C.evaluate(0,mu)*y;
                end
            else
                warning('KerMor:NoOutputConversion',['No system output'...
                    'conversion set. Forgot to set property C? Result'...
                    'equals state variables.']);
            end
        end
        
        function value = get.Times(this)
            value = 0:this.dt:this.T;
        end
    end
    
    methods(Access=protected,Sealed)
        
        function [t,x] = computeTrajectory(this, mu, inputidx)
            % Computes a solution/trajectory for the given mu and inputidx.
            %
            % Parameters:
            % mu: The concrete mu parameter sample to simulate for.
            % inputidx: The index of the input function to use.
            %
            % Leave parameters empty if the system does not have the
            % according features (it will be ignored anyway)
            %
            % Return values:
            % t: The times at which the model was evaluated
            % x: Depending on the existance of an output converter, this
            %    either returns the full trajectory or the processed output
            %    at times t.
                        
            % Set ODE options
            opts = [];
            if ~isempty(this.System.MaxTimestep)
                opts = odeset('MaxStep',this.System.MaxTimestep, 'InitialStep',.5*this.System.MaxTimestep);
            end
            
            % Get target ODE function
            odefun = this.System.getODEFun(mu, inputidx);
            
            % Get initial x
            x0 = this.System.x0(mu);
            
            % Solve ODE
            [t,x] = this.ODESolver(odefun, this.Times, x0, opts);
            % Transpose result to match inner data structure
            x = x';
        end
        
        function checkType(this, obj, type)
            % Object typechecker.
            % Checks if a given object is of the specified type and throws
            % an error if not.
            % Convenience method.
            if ~isa(obj, type)
                error(['Wrong type ''' class(obj) ''' for this property. Has to be a ' type]);
            end
        end
    end
    
    methods(Abstract)
        % Abstract simulation method. Computes a solution/trajectory for
        % the given mu and inputidx. Leave parameters empty if the system
        % does not have the according features (it will be ignored anyway)
        %[t,x] = simulate(mu, inputidx);
    end
    
    
    
end

