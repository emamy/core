classdef BaseApprox < dscomponents.ICoreFun & ICloneable
    %BASEAPPROX Summary of this class goes here
    %
    %   In the approximation context, the triple of x,t and mu is put
    %   together to a column vector when evaluated. This central class
    %   manages the combination/splitting process, no other subclass should
    %   (have to) know about the actual combination choice.
    %   See  evaluate and splitTripleVect for details.
    %   
    % @Daniel Wirtz, 11.03.2010
    
    properties(Access=private)
        % The system's full dimension for readonly-use in subclasses.
        xDim;
        % If yet projected, the projection matrix
        % Used if no custom projection of the approximation is possible
        V = [];
    end
    
    properties(Access=protected)
        % Protected Boolean flag.
        % Set to true in subclasses to tell the BaseApprox class that the
        % projection of the approximated function is performed in the
        % subclass.
        % For \c true a setting like `f = f^r(z)` is assumed, for \c false
        % simply `f = f^r(z)` is computed.
        %
        % Defaults to false.
        CustomProjection = false;
    end
    
    methods
        
        function projApprox = project(this, V)
            % Update x (=system) Dimension (needed for triple vector
            % decomposition)
            
            if this.CustomProjection
                % Call template method to ensure subclasses perform
                % projection-specific computations.
                projApprox = this.customProject(V);
            else
                % Call last subclasses' clone (will call this local clone,
                % too)
                projApprox = this.clone;
                % Set V data
                projApprox.V = V;
            end
            projApprox.xDim = size(V,2);
        end
        
        function fx = evaluate(this, x, t, mu)
            % Evaluates the f-approximation. Depending on a possible
            % projection and the CustomProjection-property the function
            % either calls the inner evaluation directly which assumes 
            % `f = f^r(z)` or projects the reduced state variable z into
            % the original space and evaluates the function there, so via
            % `f = V'f(Vz)`
            
            % NOTE: For speed reasons this method does not use
            % compileTripleVect.
            if this.CustomProjection || isempty(this.V)
                fx = this.evaluate_approximation([t; x; mu]);
            else
                fx = this.V'*this.evaluate_approximation([t; this.V*x; mu]);
            end
        end
        
        function approximateCoreFun(this, model)
            ps = model.Data.ParamSamples;
            xi = model.Data.PlainSnapshotArray;
            numPSamples = size(model.Data.Snapshots,3);
            numInSamples = size(model.Data.Snapshots,4);
            
            % Repeat times vector as often as different param/input combinations exist
            times = repmat(model.Times,1,numPSamples*numInSamples);
            
            % Repeat ParamSamples for each time, then later the result for
            % each input. Leave [] if no params are given.
            mu = [];
            if (model.Data.SampleCount > 0)
                idx = repmat(...
                    reshape(...
                    repmat(1:size(ps,2),length(model.Times),1),...
                    1,[]),...
                    1,numInSamples...
                    );
                mu = model.Data.ParamSamples(:,idx);
            end
                
            % Store system dimension for later use
            this.xDim = size(xi,1);
            
            % Compile whole snapshot vector
            data_base = [times; xi; mu];
            
            % Call template method
            this.gen_approximation_data(data_base, model.Data.fValues(:,:));
        end
        
        function set.CustomProjection(this, value)
            if ~islogical(value)
                error('Property "CustomProjection" must be logical/boolean.');
            end
            this.CustomProjection = value;
        end
    end
    
    methods(Access=protected)
        function [x,t,mu] = splitTripleVect(this, vect)
            % Splits up the previously combined triple vectors into its
            % components.
            %
            % See also: compileTripleVect
            t = vect(1,:);
            x = vect(2:this.xDim+1,:);
            mu = [];
            % Only extract mu if parameters are set
            if size(vect,1) > this.xDim+1;
                mu = vect(this.xDim+2:end,:);
            end
        end
        
        function vect = compileTripleVect(this, x, t, mu)%#ok
            % Forms a triple vector from the input arguments x,t,mu.
            % For each the column count must be identical.
            %
            % IMPORTANT: Any changes here must be reflected also within the
            % local "evaluate" function! (speed reasons)
            % 
            % See also: splitTripleVect
            vect = [t; x; mu];
        end
        
        function target = clone(this, target)
            % Creates a copy of this object into the target subclass
            % instance.
            %
            % See also: ICloneable
            if nargin == 1
                error('Cloning a BaseApprox class needs a subclass to fill.');
            end
            target.xDim = this.xDim;
            target.V = this.V;
            target.CustomProjection = this.CustomProjection;
        end
    end   
    
    methods(Abstract, Access=protected)
        % Evaluates the approximated function at point x
        %
        % Template method.
        fx = evaluate_approximation(x);
        
        % Computes the approximation according to the concrete
        % approximation strategy. 
        %
        % Template method.
        gen_approximation_data(xi,fxi);
        
        % Implement this method to perform any approximation data-structure
        % specific operations. Needs to return a NEW INSTANCE of the class
        % in order to keep old settings. Use the clone method to create
        % copies.
        %
        % See also: IProjectable ICloneable BaseApprox/clone
        projApprox = customProject(V);
    end
    
end

