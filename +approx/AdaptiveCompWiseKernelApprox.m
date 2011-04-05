classdef AdaptiveCompWiseKernelApprox < approx.BaseCompWiseKernelApprox
    % Adaptive component-wise kernel approximation algorithm
    %
    %
    % @author Daniel Wirtz @date 2011-03-31
    %
    % See also: BaseApprox BaseCompWiseKernelApprox
    %
    % @new{0,3,dw,2011-04-01} Added this class.
    
    properties
        % The maximum size of the approximation training data set.
        %
        % @default 5000
        % See also: selectTrainingData
        MaxTrainingSize = 5000;
        
        % The maximum size of the expansion to produce.
        %
        % Equals the maximum number of iterations to perform during
        % adaptive approximation computation as each iteration yields a new
        % center.
        %
        % @default 200
        MaxExpansionSize = 150;
        
        dfact = .05; % 2;
        gameps = .7; %.5;
    end

    methods
        
        function approximateCoreFun(this, model)
            % Performs adaptive approximation generation.
            %
            % @docupdate
            % @todo Think about suitable stopping condition (relative error
            % change?)
            
            %% Checks
            % This algorithm so far works only with Gaussian kernels
            if ~isa(this.SystemKernel,'kernels.GaussKernel') || ...
                    (~isa(this.TimeKernel,'kernels.GaussKernel') && ~isa(this.TimeKernel,'kernels.NoKernel')) || ...
                    (~isa(this.ParamKernel,'kernels.GaussKernel') && ~isa(this.ParamKernel,'kernels.NoKernel'))
                error('Any kernels used have to be Gaussian kernels for this approximation algorithm so far');
            end    
            
            %% Initializations
            atd = model.Data.ApproxTrainData;
            fx = model.Data.ApproxfValues;
            params = model.Data.getParams(atd(1,:));
            nx = general.NNTracker;
            nt = general.NNTracker;
            np = general.NNTracker;
            
            %% Compute bounding boxes & center
            [BXmin, BXmax] = general.Utils.getBoundingBox(atd(4:end,:));
            [BPmin, BPmax] = general.Utils.getBoundingBox(...
                model.Data.getParams(unique(atd(1,:))));
            Btmin = min(atd(3,:)); Btmax = max(atd(3,:));
            thecenter = [BXmin+BXmax; Btmin+Btmax; BPmin + BPmax]/2;
            
            %% Select initial center x0
            % Strategy A: Take the point that is closest to the bounding
            % box center!
            A = repmat(thecenter, 1, size(atd,2));
            B = [atd(3:end,:); params];
            [dummy, inIdx] = min(sum((A-B).^2,1));
            clear A B;
            
            % Strategy B: Select point in indices-middle
            %inIdx = round(size(atd,2)/2);
            
            initial = atd(:,inIdx);
            this.Centers.xi = initial(4:end);
            this.Centers.ti = initial(3);
            this.Centers.mui = params(:,inIdx);
            % Add points to nearest neighbor trackers (for gamma comp)
            nx.addPoint(initial(4:end));
            nt.addPoint(initial(3));
            np.addPoint(this.Centers.mui);
            
            %% Set up initial expansion
            used = inIdx;
            this.Ma = fx(:,inIdx);
            this.off = [];
            
            %% Choose initial gamma
            % Current strategy: Choose Gammas that make Gaussians equal .5
            % over half the shortest bounding boxes dimension!
            def_dx = norm(BXmax - BXmin) / this.dfact;
            this.SystemKernel.setGammaForDistance(def_dx,this.gameps);
            if ~isa(this.TimeKernel,'kernels.NoKernel')
                def_dt = (Btmax - Btmin) / this.dfact;
                this.SystemKernel.setGammaForDistance(def_dt,this.gameps);
            end
            if ~isa(this.ParamKernel,'kernels.NoKernel')
                def_dp = norm(BPmax - BPmin) / this.dfact;
                this.SystemKernel.setGammaForDistance(def_dp,this.gameps);
            end
            
            %% Outer control loop
            cnt = 0;
            while true
                
                %% Determine maximum error over training data
                fhat = this.evaluate(atd(4:end,:),atd(3,:),params);
                errs = sum((model.Data.ApproxfValues - fhat).^2);
                [val, maxidx] = max(errs);
                rel = val / (norm(model.Data.ApproxfValues(maxidx))+eps);
                
                if KerMor.App.Verbose > 1
                    fprintf('Max error over training data: %5.20f (relative: %10.20f)\n',val,rel);
                    pos = [1 3];
                    if KerMor.App.Verbose > 2
                        pos = 1;
                    end
                    subplot(2,2,pos);
                    plot(1:length(errs),errs,'r',used,val,'r*',maxidx,val,'b*');
                    subplot(2,2,pos+1);
                    hold off;
                    plot([BXmin, BXmax],'black');
                    axis tight;
                    hold on;
                    %plot((BXmin+BXmax)/2,'g');
                    if size(this.Centers.xi,2) > 1
                        plot(this.Centers.xi(:,1:end-1),'.','MarkerSize',2);
                    end
                    plot(this.Centers.xi(:,end),'r*','MarkerSize',3);
                    
                    % Plot params & time also
                    if KerMor.App.Verbose > 2
                        subplot(2,2,3); hold off;
                        plot([Btmin, Btmax],'black');
                        axis tight;
                        hold on;
                        if size(this.Centers.ti,2) > 1
                            plot(this.Centers.ti(:,1:end-1),'.','MarkerSize',2);
                        end
                        plot(this.Centers.ti(:,end),'r*','MarkerSize',3);
                        subplot(2,2,4); hold off;
                        plot([BPmin, BPmax],'black');
                        axis tight;
                        hold on;
                        if size(this.Centers.mui,2) > 1
                            plot(this.Centers.mui(:,1:end-1),'.','MarkerSize',2);
                        end
                        plot(this.Centers.mui(:,end),'r*','MarkerSize',3);
                    end
                    
                    pause;
                end
                
                %% Stopping condition
                if cnt == this.MaxExpansionSize || rel < 1e-9
                    %pdata = [1 used size(atd,2)];
                    %plot(pdata,1,'r*');
                    break;
                end
                % Add maxidx to list of used centers
                used(end+1) = maxidx;%#ok
                
                %% Extend centers
                new = atd(:,maxidx);
                this.Centers.xi(:,end+1) = new(4:end);
                this.Centers.ti(end+1) = new(3);
                this.Centers.mui(:,end+1) = params(:,maxidx);
                
                %% Compute new gamma
                % Add points to nearest neighbor trackers (for gamma comp)
                nx.addPoint(new(4:end));
                nt.addPoint(new(3));
                np.addPoint(params(:,maxidx));
                % Update Kernels Gamma values
                d = nx.getMaxNN/this.dfact;
                gs = this.SystemKernel.setGammaForDistance(d,this.gameps);
                if KerMor.App.Verbose > 2
                    fprintf('NN state space distances: %s\n',num2str(nx.NNDists));
                    fprintf('Kernel distances - System:%10f => gamma=%f',d,gs);
                end
%                 if ~isa(this.TimeKernel,'kernels.NoKernel')
%                     if isinf(nt.getMaxNN)
%                         distt = def_dt;
%                     else
%                         distt = nt.getMaxNN;
%                     end
%                     gt = this.SystemKernel.setGammaForDistance(distt/this.dfact,this.gameps);
%                     if KerMor.App.Verbose > 2
%                         fprintf(', Time:%10f => gamma=%10f',distt/this.dfact,gt);
%                     end
%                 end
                if ~isa(this.ParamKernel,'kernels.NoKernel')
                    if isinf(np.getMaxNN)
                        distp = def_dt;
                    else
                        distp = np.getMaxNN;
                    end
                    gp = this.SystemKernel.setGammaForDistance(distp/this.dfact,this.gameps);
                    if KerMor.App.Verbose > 2
                        fprintf(', Param:%10f => gamma=%10f',distp/this.dfact,gp);
                    end
                end
                if KerMor.App.Verbose > 2
                    fprintf('\n');
                end
                                
                %% Compute coefficients
                % Call coeffcomp preparation method and pass kernel matrix
                this.CoeffComp.init(this.getKernelMatrix);

                % Call protected method
                this.computeCoeffs(model.Data.ApproxfValues(:,used));
                
                cnt = cnt+1;
            end
            
                    
%             % Reduce the snapshot array and coeff data to the
%             % really used ones! This means if any snapshot x_n is
%             % not used in any dimension, it is kicked out at this
%             % stage.
%             hlp = sum(this.Ma ~= 0,1);
%             usedidx = find(hlp > 0);
%             if length(usedidx) < n
%                 this.Ma = this.Ma(:,usedidx);
%                 this.Centers.xi = xi(:,usedidx);
%                 if ~isempty(ti)
%                     this.Centers.ti = ti(:,usedidx);
%                 end
%                 if ~isempty(mui)
%                     this.Centers.mui = mui(:,usedidx);
%                 end
%             end
%             
%             % @todo find out when sparse representation is more
%             % efficient!
%             if sum(hlp) / numel(this.Ma) < .5
%                 this.Ma = sparse(this.Ma);
%             end
%             
%             % dont use offset vector if none are given
%             if all(this.off == 0)
%                 this.off = [];
%             end     
        end
                
        function atd = selectTrainingData(this, modeldata)
            % Selects a subset of the projection training data if the size
            % exceeds the value given in the MaxTrainingSize property.
            % If so, MaxTrainingSize samples are taken linearly spaced over
            % the whole training data.
            %
            % Important:
            % Note that the selected training data is projected into the
            % precomputed subspace if spacereduction is performed.
            %
            % Overrides the default method in BaseApprox.
            %
            % See also:
            % models.BaseFullModel.off4_genApproximationTrainData
            
            % Validity checks
            sn = modeldata.TrainingData;
            if isempty(sn)
                error('No projection training data available to take approximation training data from.');
            end
            if (size(sn,2) > this.MaxTrainingSize)
                selection = round(linspace(1,size(sn,2),this.MaxTrainingSize));
                atd = sn(:,selection);
            else
                atd = sn;
            end
        end
        
        function target = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            target = approx.AdaptiveCompWiseKernelApprox;
            
            target = clone@approx.BaseCompWiseKernelApprox(this, target);
            
            % copy local props
            copy.MaxTrainingSize = this.MaxTrainingSize;
            copy.MaxExpansionSize = this.MaxExpansionSize;
            copy.dfact = this.dfact;
            copy.gameps = this.gameps;
        end       
    end   
end


