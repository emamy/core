classdef BaseApprox < dscomponents.ACoreFun
    %BASEAPPROX Abstract base class for all core function approximations
    %
    %
    % @author Daniel Wirtz @date 11.03.2010
    
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
        
        function atd = selectTrainingData(this, modeldata)
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
        gen_approximation_data(model, xi, ti, mui);
    end
    
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
            model.Data.ProjTrainData = x;
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
