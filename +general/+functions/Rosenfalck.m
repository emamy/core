classdef Rosenfalck < general.functions.AFunGen
    
    properties
        % The baseline for no signal
        basemV = -80;
        
        % The time until the signal starts
        start_ms = 1.65;
        
        % Amplitude factor.
        %
        % Fitted by eye-ball norm to the result plots of
        % general.functions.Rosenfalck.detectDefaultScaling
        %
        % Literature / original: 96
        %
        % @default 84
        A = 84;
        
        % Factor to scale the rosenfalck-shape along the x-axis
        %
        % The default value has been obtained by minimizing the averaging
        % absolute error against precomputed shorten-shapes
        % parameters 121:360 of the mus.mat matrix in +models/+emg/data
        %
        % Using Linf-norm, the ideal factor is 2.6.
        %
        % Best factors for precomputed shapes using L2/Linf: 2.7/2.6
        %
        % @default 2.7
        xscale = 2.7;
    end
    
    methods
        
        function [fhandle, dfhandle] = getFunction(this)
            minV = this.basemV;
            A = this.A;
            s = this.start_ms;
            factor = this.xscale;
            fhandle = @(x)A*(factor*max(0,x-s)).^3.*exp(-factor*max(0,x-s))+minV;
            dfhandle = @(x)A*(factor*max(0,x-s)).^2.*exp(-factor*max(0,x-s)).*(3-factor*max(0,x-s));
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
        % Older function for scaling detection of ShapeData_v1 like shapes
%         function detectDefaultScaling
%             dd = models.emg.Model.DataDir; 
%             s = load(fullfile(dd,'ShapeData_v1.mat'));
%             Times = s.Times;
%             Shapes = s.Shapes;
%             mus = s.fibretypes;
%             
%             %% Error function
% %             efun = @Norm.L2;
%             efun = @Norm.Linf;
% 
%             %%
%             r = general.functions.Rosenfalck;
%             r.basemV = -80;
% 
%             % Need only fibre type param here
%             [mus,sidx] = sort(mus);
%             
%             %%
%             nmu = size(mus,2);
%             factors = 3:.1:5;
%             z = 0:.05:10;
%             nf = length(factors);
%             pm = PlotManager(false,3,3);
%             pm.LeaveOpen = true;
%             for fidx = 1:nf
%                 r.xscale = factors(fidx);
%                 rosenfalck = r.getFunction;
%                 for idx = 1:nmu
%                     k = sidx(idx);
%                     rshape = rosenfalck(Times{k});
%                     diff = (rshape-Shapes{k})';
%                     errors(1,k) = efun(diff);%#ok
%                 end
%                 plot(pm.nextPlot,Times{k},Shapes{k},'r',z,rosenfalck(z),'b',Times{k},rshape,'b.');
%                 fe(fidx) = max(errors(1,:));%#ok
%             end
%             plot(pm.nextPlot('','absolute error over xscales'),factors,fe);
%             [minval, min_idx] = min(fe)
%             bestscaling = factors(min_idx)
%             efun
%         end
    end
    
end