classdef BaseApprox < dscomponents.ACoreFun
    %BASEAPPROX Abstract base class for all core function approximations
    %
    %   In the approximation context, the triple of x,t and mu is put
    %   together to a column vector when evaluated. This central class
    %   manages the combination/splitting process, no other subclass should
    %   (have to) know about the actual combination choice.
    %   See  evaluate and splitTripleVect for details.
    %
    % @author Daniel Wirtz @date 11.03.2010
    
    %     properties(SetAccess=private, GetAccess=protected)
    %         % The system's full dimension for readonly-use in subclasses.
    %         xDim;
    %     end
    
%     properties(Access=private)
%         % The system's full dimension for readonly-use in subclasses.
%         xDim;
%         
% %         % The projection matrix used (if projection is used)
% %         %
% %         % Used if no custom projection of the approximation is possible
% %         V = [];
% %         
% %         % The biorthogonal matrix to V (if projection is used)
% %         %
% %         % Used if no custom projection of the approximation is possible
% %         W = [];
%     end
    
%     properties(Access=protected)
%         % Protected Boolean flag.
%         % Set to true in subclasses to tell the BaseApprox class that the
%         % projection of the approximated function is performed in the
%         % subclass.
%         % For \c true a setting like `f = f^r(z)` is assumed, for \c false
%         % simply `f = f^r(z)` is computed.
%         %
%         % Defaults to false.
%         CustomProjection = false;
%     end
    
    methods
        
%         function projApprox = project(this, V, W)
%             % Default projection method for Approximations.
%             % Override in subclasses to specialize behaviour.
%             
%             
%             %             if this.CustomProjection
%             %                 % Call template method to ensure subclasses perform
%             %                 % projection-specific computations.
%             %                 projApprox = this.customProject(V, W);
%             %             else
%             % Call last subclasses' clone (will call this local clone,
%             % too)
%             projApprox = this.clone;
%             % Set V data
%             projApprox.V = V;
%             projApprox.W = W;
%             %end
%             projApprox.xDim = size(V,2);
%         end
        
%         function fx = evaluate(this, x, t, mu)
%             % Evaluates the f-approximation. Depending on a possible
%             % projection and the CustomProjection-property the function
%             % either calls the inner evaluation directly which assumes
%             % `f = f^r(z)` or projects the reduced state variable z into
%             % the original space and evaluates the function there, so via
%             % `f = V'f(Vz)`
%             if isempty(this.V) || isempty(this.W)
%                 fx = this.evaluate_approximation(x, t, mu);
%             else
%                 fx = this.W'*this.evaluate_approximation(this.V*x, t, mu);
%             end
%         end
        
        function approximateCoreFun(this, model)
            ps = model.Data.ParamSamples;
            xi = model.Data.PlainSnapshotArray;
            numPSamples = size(model.Data.Snapshots,3);
            numInSamples = size(model.Data.Snapshots,4);
            
            % Repeat times vector as often as different param/input combinations exist
            ti = repmat(model.Times,1,numPSamples*numInSamples);
            
            % Repeat ParamSamples for each time, then later the result for
            % each input. Leave [] if no params are given.
            mui = [];
            if (model.Data.SampleCount > 0)
                idx = repmat(...
                    reshape(...
                    repmat(1:size(ps,2),length(model.Times),1),...
                    1,[]),...
                    1,numInSamples...
                    );
                mui = model.Data.ParamSamples(:,idx);
            end
            
            % Store system dimension for later use
            %this.xDim = size(xi,1);
            
            % Compile whole snapshot vector
            %             data_base = [times; xi; mu];
            
            % Call template method
            this.gen_approximation_data(xi, ti, mui, model.Data.fValues(:,:));
        end
        
        
    end
    
%     methods(Access=protected)
%         
%         function target = clone(this, target)
%             % Creates a copy of this object into the target subclass
%             % instance.
%             %
%             % See also: ICloneable
%             if nargin == 1
%                 error('Cloning a BaseApprox class needs a subclass to fill.');
%             end
%             target.xDim = this.xDim;
%             target.V = this.V;
%             target.W = this.W;
%             target.CustomProjection = this.CustomProjection;
%         end
%     end
    
    methods(Abstract, Access=protected)
        % Evaluates the approximated function at point x
        %
        % Parameters:
        % x: The state variable column vectors at which the approximation
        % should be evaluated.
        %
        % Return values:
        % fx: The approximated values of the core function `\hat{f}` at
        % `x`.
        %
        % Template method.
        %fx = evaluate_approximation(x, t, mu);
        
        % Computes the approximation according to the concrete
        % approximation strategy.
        %
        % Template method.
        %
        % Parameters:
        % xi: The state variable snapshots as column vectors, each row
        % representing each system dimension.
        % fxi: The core functions' values at the snapshots xi.
        gen_approximation_data(xi, ti, mui, fxi);
        
        % Implement this method to perform any approximation data-structure
        % specific operations. Needs to return a NEW INSTANCE of the class
        % in order to keep old settings. Use the clone method to create
        % copies.
        %
        % Parameters:
        % V: The columnwise-normed projection matrix
        % W: The V-biorthogonal matrix
        %
        % Return values:
        % projApprox: A new handle instance with a projected version of the
        % approximation function.
        %
        % See also: IProjectable ICloneable BaseApprox/clone
        %projApprox = customProject(V,W);
    end
    
    %         function [x,t,mu] = splitTripleVect(this, vect)
    %             % Splits up the previously combined triple vectors into its
    %             % components.
    %             %
    %             % See also: compileTripleVect
    %             t = vect(1,:);
    %
    %             if this.xDim+1 > size(vect,1)
    %                 disp('adg');
    %             end
    %             x = vect(2:this.xDim+1,:);
    %             mu = [];
    %             % Only extract mu if parameters are set
    %             if size(vect,1) > this.xDim+1;
    %                 mu = vect(this.xDim+2:end,:);
    %             end
    %         end
    %
    %         function vect = compileTripleVect(this, x, t, mu)%#ok
    %             % Forms a triple vector from the input arguments x,t,mu.
    %             % For each the column count must be identical.
    %             %
    %             % IMPORTANT: Any changes here must be reflected also within the
    %             % local "evaluate" function! (speed reasons)
    %             %
    %             % See also: splitTripleVect
    %             vect = [t; x; mu];
    %         end
    
    methods(Static)
        function res = testApproxProjections
            a{1} = approx.CompWiseInt;
            a{2} = approx.CompWiseSVR;
            a{3} = approx.CompWiseLS;
            
            b = cell(length(a),0);
            
            ts = testing.testsettings;
            samples = 30;
            x = rand(ts.testdim, samples);
            fx = ts.fnlin(x,repmat(1:samples,ts.testdim,1));
            
            model = models.BaseFullModel;
            model.Data.Snapshots = x;
            r = spacereduction.PODReducer;
            v = r.generateReducedSpace(model);
            
            res = true;
            for idx=1:length(a)
                try
                    app = a{idx};
                    mc = metaclass(app);
                    name = mc.Name;
                    cprintf(testing.MUnit.GreenCol,['Testing ' name '...']);
                    app.gen_approximation_data(x,fx);
                    b{idx} = app.project(v);
                    
                    ifxfull = app.evaluate_approximation(x);
                    ifxred = v * b{idx}.evaluate_approximation(v' * x);
                    figure(idx+20);
                    plot(1:ts.testdim,sum(abs(ifxfull-ifxred),2));
                    title(name);
                    xlabel('dimension');
                    ylabel('error');
                catch ME
                    res = false;
                    disp(getReport(ME));
                end
            end
        end
    end
    
end

