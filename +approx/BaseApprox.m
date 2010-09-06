classdef BaseApprox < dscomponents.ACoreFun
    %BASEAPPROX Abstract base class for all core function approximations
    %
    %
    % @author Daniel Wirtz @date 11.03.2010
    
    methods
        function approximateCoreFun(this, model)
            
            % Load snapshots
            sn = model.Data.ApproxTrainData;
            
            % Compile necessary data
            xi = sn(4:end,:);
            ti = sn(3,:);
            muidx = sn(1,:);
            if all(muidx == 0)
                mui = [];
            else
                mui = model.Data.ParamSamples(:,muidx);
            end
            
            % Call template method
            this.gen_approximation_data(xi, ti, mui, model.Data.ApproxfValues);
        end
    end
    
    methods(Abstract, Access=protected)
        % Computes the approximation according to the concrete
        % approximation strategy.
        %
        % Template method.
        %
        % Parameters:
        % xi: The state variable snapshots as column vectors, each row
        % representing each system dimension.
        % ti: The times at which the state variable snapshots are taken
        % mui: The parameters used to obtain the snapshots. Equals [] if no
        % parameters have been used.
        % fxi: The core functions' values at the snapshots xi.
        gen_approximation_data(xi, ti, mui, fxi);
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

