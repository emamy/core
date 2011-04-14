classdef BaseApprox < dscomponents.ACoreFun
    % Abstract base class for all core function approximations
    %
    % Simply provides two methods: 
    % - selectTrainingData: Used to select a (sub-)set of the training data
    % - approximateCoreFun: Abstract template method that performs the
    % actual approximation. Possible algorithms may be i.e. component-wise
    % approximation, multidimensional approximation or any other.
    %
    % @author Daniel Wirtz @date 2010-03-11
    %
    % @new{0,3,dw,2011-04-12} New property
    % approx.BaseApprox.TrainDataSelector. Allows to choose different
    % strategies for training data selection.
    %
    % @change{0,3,dw,2011-04-01} Removed the gen_approximation_data
    % template method as it basically only compiled the single xi,ti,mui
    % values which will be re-joined into the
    % dscomponents.AKernelCoreFun.Centers property in all currently
    % implemented approximation strategies. Some similar intermediate
    % method may be reintroduced later.
    
    properties
        % The algorithm that selects the approximation training data.
        %
        % See also: DefaultSelector LinspaceSelector TimeSelector
        % @default approx.selection.TimeSelector
        TrainDataSelector;
    end
    
    methods
        function this = BaseApprox
            this.TrainDataSelector = approx.selection.TimeSelector;
        end
        
        function copy = clone(this, copy)
            copy = clone@dscomponents.ACoreFun(this, copy);
            copy.TrainDataSelector = this.TrainDataSelector;
        end
    end
    
    methods(Abstract)
        % Computes the approximation according to the concrete
        % approximation strategy.
        %
        % Template method.
        %
        % Parameters:
        % model: The full model
        approximateCoreFun(this, model);
    end
    
    methods(Static)
        function res = test_ApproxProjections
            a{1} = approx.DefaultCompWiseKernelApprox;
            a{1}.CoeffComp = general.interpolation.KernelInterpol;
            a{2} = approx.DefaultCompWiseKernelApprox;
            a{2}.CoeffComp = general.regression.ScalarEpsSVR;
            a{2} = approx.DefaultCompWiseKernelApprox;
            a{2}.CoeffComp = general.regression.KernelLS;
            a{3} = approx.AdaptiveCompWiseKernelApprox;
            a{3}.CoeffComp = general.interpolation.KernelInterpol;
            a{3}.MaxExpansionSize = 20;
            
            b = cell(length(a),0);
            
            ts = testing.testsettings;
            samples = 50;
            x = rand(ts.testdim, samples);
            fx = ts.fnlin(x,repmat(1:samples,ts.testdim,1));
            
            model = ts.m;
            model.Data.TrainingData = [zeros(3,size(x,2)); x];
            model.Data.ApproxTrainData = model.Data.TrainingData;
            model.Data.ApproxfValues = fx;
            v = model.SpaceReducer.generateReducedSpace(model);
            
            res = true;
            for idx=1:length(a)
                try
                    app = a{idx};
                    mc = metaclass(app);
                    name = mc.Name;
                    cprintf(testing.MUnit.GreenCol,['Testing ' name '...\n']);
                    app.approximateCoreFun(model);
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

