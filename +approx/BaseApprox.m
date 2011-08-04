classdef BaseApprox < dscomponents.ACoreFun
    % Abstract base class for all core function approximations inside dynamical systems.
    %
    % Simply provides two methods:
    % - selectTrainingData: Used to select a (sub-)set of the training data
    % - approximateData: Abstract template method that performs the
    % actual approximation. Possible algorithms may be i.e. component-wise
    % approximation, multidimensional approximation or any other.
    %
    % @author Daniel Wirtz @date 2010-03-11
    %
    % @change{0,5,dw,2011-07-07} 
    % - Changed the interface for the old approximateCoreFun to
    % approximateSystemFunction that now takes only the full model instance instead of
    % `x_i,t_i,\mu_i,f_{x_i}` values.
    % - The new default value for the approx train data selector is the
    % approx.selection.DefaultSelector.
    %
    % @change{0,4,dw,2011-05-19} Disconnected the Approx classes from taking a BaseModel instance at
    % approx computation. This way external tools can use the approximation algorithms, too.
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
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
    
    properties(SetObservable)
        % The algorithm that selects the approximation training data.
        %
        % @propclass{important} Determines the strategy used to select the approximation training
        % data
        %
        % @default approx.selection.DefaultSelector
        %
        % @type approx.selection.ASelector
        %
        % @todo MOVE TO BASEFULLMODEL
        %
        % See also: DefaultSelector LinspaceSelector TimeSelector
        TrainDataSelector;
    end
    
    methods
        function set.TrainDataSelector(this, value)
            this.checkType(value, 'approx.selection.ASelector');%#ok
            this.TrainDataSelector = value;
        end
        
        function this = BaseApprox
            this.TrainDataSelector = approx.selection.DefaultSelector;
        end
        
        function copy = clone(this, copy)
            copy = clone@dscomponents.ACoreFun(this, copy);
            copy.TrainDataSelector = this.TrainDataSelector.clone;
        end
    end
    
    methods(Abstract)
        % Computes the approximation according to the concrete
        % approximation strategy.
        %
        % Template method.
        %
        % Parameters:
        % model: The current model.
        approximateSystemFunction(this, model);
    end
    
    methods(Static)
        function res = test_ApproxProjections
            a{1} = approx.algorithms.DefaultCompWiseKernelApprox;
            a{1}.CoeffComp = general.interpolation.KernelInterpol;
            a{2} = approx.algorithms.DefaultCompWiseKernelApprox;
            a{2}.CoeffComp = general.regression.ScalarEpsSVR;
            a{2} = approx.algorithms.DefaultCompWiseKernelApprox;
            a{2}.CoeffComp = general.regression.KernelLS;
            a{3} = approx.algorithms.AdaptiveCompWiseKernelApprox;
            a{3}.CoeffComp = general.interpolation.KernelInterpol;
            a{3}.MaxExpansionSize = 20;
            
            b = cell(length(a),0);
            
            ts = testing.testsettings;
            samples = 50;
            kexp = approx.KernelApprox;
            x = rand(ts.testdim, samples);
            t = 1:size(x,2);
            fx = ts.fnlin(x,repmat(1:samples,ts.testdim,1));

            pr = spacereduction.PODReducer;
            pr.Value = 2;
            pr.Mode = 'abs';
            v = pr.computePOD(x);
            
            res = true;
            for idx=1:length(a)
                try
                    app = a{idx};
                    mc = metaclass(app);
                    name = mc.Name;
                    cprintf(testing.MUnit.GreenCol,['Testing ' name '...\n']);
                    app.computeApproximation(kexp, x, t, [], fx);
                    b{idx} = kexp.project(v,v);
                    
                    ifxfull = kexp.evaluate(x,t,[]);
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

