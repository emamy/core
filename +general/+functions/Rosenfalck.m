classdef Rosenfalck < general.functions.AFunGen
    
    properties
        basemV = -80;
        A = 96;
        % Factor to scale the rosenfalck-shape along the x-axis
        %
        % The default value has been obtained by minimizing the averaging
        % absolute error agains shorten-shapes for 200 fibre-type
        % parameters.
        %
        % Using Linf-norm, the ideal factor is 2.6.
        %
        % @default 3
        xscale = 3;
    end
    
    methods
        
        function [fhandle, dfhandle] = getFunction(this)
            minV = this.basemV;
            A = this.A;
            factor = this.xscale;
            fhandle = @(x)A*(factor*x).^3.*exp(-factor*x)+minV;
            dfhandle = @(x)A*(factor*x).^2.*exp(-factor*x).*(3-factor*x);
        end
        
        function pm = plot(this, varargin)
            varargin(end+1:end+2) = {'R' 0:.1:15};
            pm = plot@general.functions.AFunGen(this,varargin{:});
        end
        
        function str = getConfigStr(this)
            str = sprintf('MinV=%g',this.basemV);
        end
    end
    
    methods(Static)
        function detectDefaultScaling
            dd = models.emg.Model.DataDir; 
            s = load(fullfile(dd,'ShapeData_v1.mat'));
            Times = s.Times;
            Shapes = s.Shapes;
            mus = s.mus;

            %%
            r = general.functions.Rosenfalck;
            r.basemV = -80;

            % Need only fibre type param here
            [mus,sidx] = sort(mus(1,:));
            %%
            nmu = size(mus,2);
            factors = 2:.1:4;
            z = 0:.05:10;
            nf = length(factors);
            pm = PlotManager(false,3,3);
            pm.LeaveOpen = true;
            for fidx = 1:nf
                r.xscale = factors(fidx);
                rosenfalck = r.getFunction;
                for idx = 1:nmu
                    k = sidx(idx);
                    rshape = rosenfalck(Times{k});
                    diff = (rshape-Shapes{k})';
                    errors(1,k) = Norm.L2(diff);%#ok
                    %errors(1,k) = Norm.Linf(diff);%#ok
                end
                plot(pm.nextPlot,Times{k},Shapes{k},'r',z,rosenfalck(z),'b',Times{k},rshape,'b.');
                fe(fidx) = max(errors(1,:));%#ok
            end
            plot(pm.nextPlot('','absolute error over xscales'),factors,fe);
        end
    end
    
end