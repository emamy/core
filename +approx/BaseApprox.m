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
    % data.selection.DefaultSelector.
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
        % @default data.selection.DefaultSelector
        %
        % @type data.selection.ASelector
        %
        % @todo MOVE TO BASEFULLMODEL
        %
        % See also: data.selection.DefaultSelector LinspaceSelector TimeSelector
        TrainDataSelector;
    end
    
    methods
        function set.TrainDataSelector(this, value)
            this.checkType(value, 'data.selection.ASelector');
            this.TrainDataSelector = value;
        end
        
        function this = BaseApprox
            this.TrainDataSelector = data.selection.DefaultSelector;
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
            
            co = approx.algorithms.Componentwise;
            co.ExpConfig = kernels.config.ExpansionConfig;
            co.ExpConfig.StateConfig = kernels.config.GaussConfig('G',1);
            
            a{1} = co;
            a{1}.CoeffConfig = general.interpolation.InterpolConfig;
            
            a{2} = co.clone;
            a{2}.CoeffConfig = general.regression.EpsSVRConfig([.1; 1]);
            a{2}.CoeffConfig.Prototype = general.regression.ScalarEpsSVR;
            
            a{3} = co.clone;
            a{3}.CoeffConfig = general.regression.NuSVRConfig([.3; 1]);
            
            a{4} = co.clone;
            a{4}.CoeffConfig = general.regression.KernelLSConfig(1);
            
            a{5} = approx.algorithms.VKOGA;
            a{5}.MaxExpansionSize = 20;
            a{5}.ExpConfig = co.ExpConfig;
            
            b = cell(length(a),0);
            
            ts = testing.testsettings;
            samples = 50;
            x = rand(ts.testdim, samples);
            t = 1:size(x,2);
            atd = data.ApproxTrainData(x, t, []);
            atd.fxi = ts.fnlin(x,repmat(1:samples,ts.testdim,1));

            pr = spacereduction.PODReducer;
            pr.Value = 2;
            pr.Mode = 'abs';
            v = pr.computePOD(x);
            ap = approx.KernelApprox;
            res = true;
            for idx=1:length(a)
                try
                    app = a{idx};
                    if isa(app,'approx.algorithms.Componentwise')
                        name = sprintf('%s with %s',class(app),class(app.CoeffConfig.Prototype));
                    else
                        name = sprintf('%s with MaxExpansionSize=%d',class(app),app.MaxExpansionSize);
                    end
                    
                    cprintf(MUnit.GreenCol,'Testing %s...',name);
                    ap.Expansion = app.computeApproximation(atd);
                    b{idx} = ap.project(v,v);
                    
                    ifxfull = ap.evaluate(x,t,[]);
                    ifxred = v * b{idx}.evaluate(v' * x,[],[]);
                    fprintf('Max L2 error: %g\n',max(Norm.L2(ifxfull-ifxred)));
                catch ME
                    res = false;
                    disp(getReport(ME));
                end
            end
        end
    end
    
end

