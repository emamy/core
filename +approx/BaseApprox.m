classdef BaseApprox < dscomponents.ACoreFun
    % Abstract base class for all core function approximations
    %
    % @author Daniel Wirtz @date 2010-03-11
    
    methods        
        
        function approximateCoreFun(this, model)
            % @todo the de-compilation of the ApproxTrainData does not
            % necessarily have to take place in this superclass; as long as
            % the data is in the ModelData class one could also move it
            % into that class. though its not necessary yet..
            
            % Load snapshots
            atd = model.Data.ApproxTrainData;
            
            % Compile necessary data
            xi = atd(4:end,:);
            ti = atd(3,:);
            muidx = atd(1,:);
            if all(muidx == 0)
                mui = [];
            else
                mui = model.Data.ParamSamples(:,muidx);
            end
            
            % Call template method
            this.gen_approximation_data(model, xi, ti, mui);
        end
        
        function atd = selectTrainingData(this, modeldata)%#ok
            % Default approx training data generation method.
            %
            % Simply takes ALL the projection training data.
            %
            % Important:
            % Note that the selected training data is projected into the
            % precomputed subspace if spacereduction is performed.
            %
            % Override in subclasses to tailor the selection behaviour to
            % specific approximation functions' needs.
            %
            % See also:
            % models.BaseFullModel.off4_genApproximationTrainData
            
            % Validity checks
            sn = modeldata.ProjTrainData;
            if isempty(sn)
                error('No projection training data available to take approximation training data from.');
            end

            atd = sn;
        end
        
    end
    
    methods(Abstract, Access=protected)
        % Computes the approximation according to the concrete
        % approximation strategy.
        %
        % Template method.
        %
        % Parameters:
        % model: The full model
        % xi: The state variable snapshots as column vectors, each row
        % representing each system dimension.
        % ti: The times at which the state variable snapshots are taken
        % mui: The parameters used to obtain the snapshots. Equals [] if no
        % parameters have been used.
        gen_approximation_data(this, model, xi, ti, mui);
    end
    
    methods(Static)
        function res = test_ApproxProjections
            a{1} = approx.DefaultCompWiseKernelApprox;
            a{1}.CoeffComp = general.interpolation.KernelInterpol;
            a{2} = approx.DefaultCompWiseKernelApprox;
            a{2}.CoeffComp = general.regression.ScalarEpsSVR;
            a{2} = approx.DefaultCompWiseKernelApprox;
            a{2}.CoeffComp = general.regression.KernelLS;
            
            b = cell(length(a),0);
            
            ts = testing.testsettings;
            samples = 50;
            x = rand(ts.testdim, samples);
            fx = ts.fnlin(x,repmat(1:samples,ts.testdim,1));
            
            model = ts.m;
            model.Data.ProjTrainData = [zeros(3,size(x,2)); x];
            model.Data.ApproxTrainData = model.Data.ProjTrainData;
            model.Data.ApproxfValues = fx;
            v = model.SpaceReducer.generateReducedSpace(model);
            
            res = true;
            for idx=1:length(a)
                try
                    app = a{idx};
                    mc = metaclass(app);
                    name = mc.Name;
                    cprintf(testing.MUnit.GreenCol,['Testing ' name '...\n']);
                    app.gen_approximation_data(model, x, [], []);
                    b{idx} = app.project(v,v);
                    
                    ifxfull = app.evaluate(x,[],[]);
                    ifxred = v * b{idx}.evaluate(v' * x,[],[]);
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

