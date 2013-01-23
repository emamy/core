classdef AdaptiveCompWiseKernelApprox < approx.algorithms.BaseAdaptiveCWKA
% Adaptive component-wise kernel approximation algorithm
%
% @author Daniel Wirtz @date 2011-03-31
%
% See also: BaseApprox KernelApprox
%
% @change{0,5,dw,2011-11-02}
% - New interface for approximation computation: Passing an data.ApproxTrainData instance now
% instead of 'xi,ti,mui' parameters.
% - Moved common properties for adaptive algorithms to new base class BaseAdaptiveCWKA and
% using it's provided methods in order to have a more compact algorithm representation.
% - Inserted a global optimum search over all configurations, so that if the best error was
% achieved somewhere in progress this expansion is taken.
%
% @change{0,5,dw,2011-09-09} 
% - Fixed setters for MaxRelErr and MaxAbsErrFactor, now setting values to 'this' instead of
% 'kexp'.
% - Added initial coefficient support
% - Using IKernelMatrix now
%
% @change{0,5,dw,2011-07-28} Changed the algorithm part so that it can also work with
% kernels.KernelExpansion instead of only on kernels.ParamTimeKernelExpansion.
%
% @new{0,5,dw,2011-07-07} Moved the old approx.AdaptiveCompWiseKernelApprox class to this class.
%
% @change{0,4,dw,2011-05-31} Added new experimental properties approx.algorithms.BaseKernelApproxAlgorithm.MinGFactor and approx.algorithms.BaseKernelApproxAlgorithm.MaxGFactor.
%
% @change{0,4,dw,2011-05-19} Disconnected the Approx classes from taking a BaseModel instance at
% approx computation. This way external tools can use the approximation algorithms, too.
%
% @change{0,3,sa,2011-04-21} Implemented Setters for all the properties
% other than NumGammas and ValidationPercent
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% @change{0,3,dw,2011-04-14}
% - Implemented some setters
% - New property approx.algorithms.AdaptiveCompWiseKernelApprox.ValidationPercent enabling a validation set
% to check for best gammas
%
% @change{0,3,dw,2011-04-06} Now works with models that dont have any
% parameters.
%
% @new{0,3,dw,2011-04-01} Added this class.
%
% @todo Think about suitable stopping condition (relative error change?)

    methods    
%         function this = AdaptiveCompWiseKernelApprox
%             this = this@approx.algorithms.BaseAdaptiveCWKA;%#ok
%         end
                        
        function copy = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            copy = approx.algorithms.AdaptiveCompWiseKernelApprox;
            
            copy = clone@approx.algorithms.BaseAdaptiveCWKA(this, copy);
        end
    end
    
    methods(Access=protected, Sealed)
        function startAdaptiveExtension(this, kexp, atd)
            % Performs adaptive approximation generation.
            %
            % Parameters:
            % kexp: The kernel expansion. @type kernels.KernelExpansion
            % atd: The approximation training data instance @type data.ApproxTrainData
            
            ec = this.ExpConfig;
            nc = ec.getNumConfigurations;
            
            this.err = zeros(nc,this.MaxExpansionSize);
            this.relerr = this.err;
            this.expsizes = zeros(nc,1);

            xi = atd.xi.toMemoryMatrix;
            fxi = atd.fxi.toMemoryMatrix;
            fxinorm = this.ErrorFun(fxi);
            fxinorm(fxinorm == 0) = 1;
            N = size(xi,2);
            kexp.clear;
            
            minerr = Inf;
            bestNV = [];
            bestc = [];
            
            %% Run loop for all desired configurations
            pi = tools.ProcessIndicator('Starting AdaptiveCompWiseKernelApprox for %d kernel configurations',nc,false,nc);
            for cidx = 1:nc
                cnt = 1;
                this.initExpansion(kexp, atd);
                
                while true 
                    
                    %% Determine maximum error over training data
                    [val, maxidx, errs] = this.getError(kexp, atd);
                    rel = val / (norm(atd.fxi(maxidx))+eps);
                    this.MaxErrors(cidx,cnt) = val;
                    if val < globminerr
                        globbestg = bestg;
                        globbestMa = bestMa;
                        globbestcnt = cnt;
                        globminerr = val;
                    end

                    %% Verbose stuff
                    if KerMor.App.Verbose > 2
                        doPlots;
                    end

                    %% Stopping condition
                    if this.checkStop(cnt, rel, val)
                        break;
                    end

                    % Add maxidx to list of used centers
                    used(end+1) = maxidx;%#ok

                    %% Extend centers
                    this.extendExpansion(kexp, atd, maxidx);

                    %% Compute new approximation
                    dt = []; dmu = [];
                    if atd.hasTime
                        dt = nt.getMinNN;
                    end
                    if atd.hasParams
                        dmu = atd.muiDia/this.NumGammas;
                        if ~isinf(np.getMinNN)
                            dmu = np.getMinNN;
                        end
                    end
                    dists = this.getDists(atd, nx.getMinNN, dt, dmu);
                
                minerr = Inf;
                for gidx = 1:size(dists,2)
                    g = this.setDistKernelConfig(kexp, dists(:,gidx));
                    %% Compute coefficients
                    % Call coeffcomp preparation method and pass kernel matrix
                    K = data.MemoryKernelMatrix(kexp.getKernelMatrix);
                    this.CoeffComp.init(K, kexp);
                    
                    % Call protected method
                    this.computeCoeffs(kexp, atd.fxi(:,used), [kexp.Ma zeros(size(atd.fxi,1),1)]);
                    
                    % get error on training data
                    val = this.getError(kexp, atd);
                    impro = (val / minerr) * 100;
                    
                    if val < minerr
                        minerr = val;
                        bestg = g;
                        bestMa = kexp.Ma;
                        if KerMor.App.Verbose > 2
                            fprintf(' b: %.5e, %3.2f%%',val,impro);
                        end
                    else
                        if KerMor.App.Verbose > 2
                            fprintf(' w: %.5e, %3.2f%%',val,impro);
                        end
                    end
                    
                    if KerMor.App.Verbose > 2
                        %fprintf(' ||Ma||:%.5e, vd-err:%.5e, prod:%.5e\n',sum(sqrt(sum(kexp.Ma.^2,1))),vdval,val*vdval);
                        %fprintf(' ||Ma||:%.5e\n',sum(sqrt(sum(kexp.Ma.^2,1))));
                    end
                end
                if KerMor.App.Verbose > 2
                    fprintf('\n');
                end
                
                %% Assign best values
                this.setKernelConfig(kexp, bestg);
                kexp.Ma = bestMa;
                
                if KerMor.App.Verbose > 1
                    if this.pte
                        fprintf('-- It: %d ---- Minerr: %f ----- Best values: System:%f, Time:%f, Param:%f ----------\n',cnt,minerr,bestg(1),bestg(2),bestg(3));
                    else
                        fprintf('-- It: %d ---- Minerr: %f ----- Best value: System:%f ----------\n',cnt,minerr,bestg(1));
                    end 
                end
                
%                 FunVis2D(kexp, atd);
%                 pause;
                
                cnt = cnt+1;
            end
            
            if KerMor.App.Verbose > 0
                fprintf('Expansion data from iteration %d.\n',globbestcnt);
            end
            this.setKernelConfig(kexp, globbestg);
            kexp.Ma = globbestMa;
            kexp.Centers.xi = kexp.Centers.xi(:,1:globbestcnt);
            if this.pte && ~isempty(kexp.Centers.ti)
                kexp.Centers.ti = kexp.Centers.ti(:,1:globbestcnt);
            end
            if this.pte && ~isempty(kexp.Centers.mui)
                kexp.Centers.mui = kexp.Centers.mui(:,1:globbestcnt);
            end
            
            if KerMor.App.Verbose > 1
                figure;
                plot(this.MaxErrors,'r');
            end
        
            function doPlots
                fprintf('Max error over training data: %5.20f (relative: %10.20f)\n',val,rel);
                pos = [1 3];
                if KerMor.App.Verbose > 3
                    pos = 1;
                end
                subplot(2,2,pos);
                plot(1:length(errs),errs,'r',used,val,'r*',maxidx,val,'b*');
                subplot(2,2,pos+1);
                hold off;
                plot([atd.Box.xmin, atd.Box.xmax],'black');
                axis tight;
                hold on;
                %plot((BXmin+BXmax)/2,'g');
                if size(kexp.Centers.xi,2) > 1
                    plot(kexp.Centers.xi(:,1:end-1),'.','MarkerSize',2);
                end
                plot(kexp.Centers.xi(:,end),'r*','MarkerSize',3);
                
                % Plot params & time also
                if KerMor.App.Verbose > 3
                    subplot(2,2,3); hold off;
                    plot([atd.Box.tmin, atd.Box.tmin],'black');
                    axis tight;
                    hold on;
                    if size(kexp.Centers.ti,2) > 1
                        plot(kexp.Centers.ti(:,1:end-1),'.','MarkerSize',2);
                    end
                    plot(kexp.Centers.ti(:,end),'r*','MarkerSize',3);
                    subplot(2,2,4); hold off;
                    plot([atd.Box.mumin, atd.Box.mumax],'black');
                    axis tight;
                    hold on;
                    if size(kexp.Centers.mui,2) > 1
                        plot(kexp.Centers.mui(:,1:end-1),'.','MarkerSize',2);
                    end                  
                    pause;
                end
                drawnow;
            end
        end 
    end
end