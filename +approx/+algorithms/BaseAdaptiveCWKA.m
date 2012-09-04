classdef BaseAdaptiveCWKA < approx.algorithms.BaseKernelApproxAlgorithm
% Base class for adaptive component-wise kernel approximation algorithms
%
% @author Daniel Wirtz @date 2011-11-02
%
% @new{0,6,dw,2012-01-26} Added a new approximation stallment detection. The properties
% MinImprovePerc and ImproveRange control at which stage the approximation is to be stopped if
% no sufficient progress in approximation error is made.
%
% @new{0,5,dw,2011-11-02} 
% - Created this class. Collects common properties of the adaptive approx algorithms and
% provides convenience methods for subclasses.
% - New interface for approximation computation: Passing an data.ApproxTrainData instance now
% instead of 'xi,ti,mui' parameters.
%
% See also: BaseApprox KernelApprox approx.algorithms
% BaseKernelApproxAlgorithm AdaptiveCompWiseKernelApprox
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetObservable)
        % The maximum size of the expansion to produce.
        %
        % Equals the maximum number of iterations to perform during
        % adaptive approximation computation as each iteration yields a new
        % center.
        %
        % @propclass{alglimit} 
        % Some text describing the importance of this property.
        %
        % @type integer @default 300
        MaxExpansionSize = 300;
        
        % The number of different Gamma values to try.
        %
        % @propclass{important} 
        %
        % @default 15 @type integer
        NumGammas = 15;
        
        % Percentage `p` of the training data to use as validation data
        %
        % Admissible values are `p\in]0,\frac{1}{2}]`.
        %
        % @propclass{optional} 
        %
        % @default .2
        % @type double
%         ValidationPercent = .2;
        
        % Value for initial Gamma choice.
        %
        % @propclass{experimental} 
        %
        % @type double
        % @default .6
        gameps = .6;
        
        % Stopping condition property. Maximum relative error that may occur
        %
        % @propclass{critical}
        %
        % @default 1e-5 @type double
        MaxRelErr = 1e-5;
        
        % Stopping condition property. 
        %
        % Factor for maximum absolute error that may occur. Factor means that this value will be
        % multiplied by the maximum value of the approximation training data f-Values used to train
        % this approximation. This way, one defines the fraction of the maximum error that is
        % maximally allowed.
        %
        % @propclass{critical}
        %
        % @default 1e-5 @type double
        MaxAbsErrFactor = 1e-5;
        
        % The error functional to use
        % 1 = `L^\infty`-Error (max diff over all vecs & dimensions)
        % 2 = `L^2`-Error (vector-wise L^2, then max)
        %
        % @propclass{experimental} 
        %
        % @default 1 @type integer
        ErrFun = 1;
        
        % 'dfun(this.MinGFactor*bxdia, kexp.MaxGFactor*bxdia);'
        % @propclass{experimental}
        %
        % @default 1 @type double
        MaxGFactor = 1;
        
        % 'dfun(this.MinGFactor*bxdia, kexp.MaxGFactor*bxdia);'
        % @propclass{experimental}
        % @default 0.05 @type double
        MinGFactor = .05;
        
        % Determines how the initial center(s) are chosen.
        %
        % Possible values:
        % - 'center' The training point closest to the geometrical center
        % if the training data's bounding box is used.
        % - 't0' Tries to find a training sample for time zero and uses
        % that. If None are found, the strategy falls back to 'center'. If
        % multiple are found, the one closest to the geometrical center is
        % chosen.
        %
        % @default 'center' @type char
        InitialCenter = 'center';
        
        % The percentage over which the error has to improve compared to the mean error of
        % the ImproveRange earlier steps.
        %
        % Set to empty to disable.
        %
        % @propclass{important} Higher improvement values may terminate the search too early
        % and thus reduce the approximation quality. Too low values might not detect stallment
        % of the approximation process early enough.
        %
        % @type double @default .05
        %
        % See also: ImproveRange
        MinImprovePerc = .05;
        
        % The range over which the error improvement is to be monitored.
        %
        % Set to empty to disable.
        %
        % @propclass{important} A too short monitoring range might stop the algorithm too
        % early, whereas a too long range will detect stallment too late.
        %
        % @type integer @default 15
        %
        % See also: MinImprovePerc
        ImproveRange = 15;
    end
    
    properties(SetAccess=protected)
        % Contains the maximum errors for each iteration/center extension step performed by the last
        % run of this algorithm.
        %
        % @default [] @type rowvec
        MaxErrors = [];
    end
    
    properties(Transient, Access=private)
        effabs;
        lasterrs;
    end
    
    properties(Transient, SetAccess=private, GetAccess=protected)
        % Determines if the kernel expansion is a
        % kernels.ParamTimeKernelExpansion
        % @type logical
        pte;
    end
    
    methods    
        function this = BaseAdaptiveCWKA
            this = this@approx.algorithms.BaseKernelApproxAlgorithm;
            
            % Register default property changed listeners
            this.registerProps('MaxExpansionSize','NumGammas',...
                'gameps','MaxRelErr','MaxAbsErrFactor','ErrFun','MaxGFactor','MinGFactor');
            %'ValidationPercent'
        end
                        
        function copy = clone(this, copy)
            % Clones the instance.
            
            copy = clone@approx.algorithms.BaseKernelApproxAlgorithm(this, copy);
            
            % copy local props
            copy.MaxExpansionSize = this.MaxExpansionSize;
            copy.NumGammas = this.NumGammas;
            copy.gameps = this.gameps;
            copy.MaxRelErr = this.MaxRelErr;
            copy.MaxAbsErrFactor = this.MaxAbsErrFactor;
            copy.MaxErrors = this.MaxErrors;
            copy.ErrFun = this.ErrFun;
            copy.MaxGFactor = this.MaxGFactor;
            copy.MinGFactor = this.MinGFactor;
            copy.effabs = this.effabs;
        end
    end
    
    methods(Access=protected, Sealed)
        function detailedComputeApproximation(this, kexp, atd)
            % Performs adaptive approximation generation.
            %
            % Parameters:
            % kexp: The kernel expansion. @type kernels.KernelExpansion
            % atd: The approximation training data instance @type data.ApproxTrainData
            
            %% Checks
            % This algorithm so far works only with Gaussian kernels
            if ~isa(kexp, 'kernels.KernelExpansion')
                error('Approximation method works only for kernel expansions.');
            elseif ~isa(kexp.Kernel,'kernels.GaussKernel')
                %error('The state kernel has to be a Gaussian for this approximation algorithm so far');
            end
            this.pte = isa(kexp,'kernels.ParamTimeKernelExpansion');
            if this.pte && ((~isa(kexp.TimeKernel,'kernels.GaussKernel') && ~isa(kexp.TimeKernel,'kernels.NoKernel')) || ...
                    (~isa(kexp.ParamKernel,'kernels.GaussKernel') && ~isa(kexp.ParamKernel,'kernels.NoKernel')))
                error('Any kernels used have to be Gaussian kernels for this approximation algorithm so far');
            end
            if size(atd.xi,2) < this.MaxExpansionSize
                warning('BaseAdaptiveCWKA:expansionsize','Only %d training samples but having MaxExpansionSize=%d, changing to %d.',size(atd.xi,2),this.MaxExpansionSize,size(atd.xi,2));
                this.MaxExpansionSize = size(atd.xi,2);
            end
            
            % Stopping condition preps
            this.effabs = this.MaxAbsErrFactor * max(atd.fxi.MaxValue,-atd.fxi.MinValue);
            this.lasterrs = [];
            
            this.detailedAdaptiveApproximation(kexp, atd);
        end
    end
    
    %% Convenience methods for use/override in concrete algorithms
    methods(Access=protected)
        function [c, idx] = getInitialCenter(this, atd)
            % Computes an initial center from the training data.
            %
            % Depending on the InitialCenter property, the following
            % strategy is pursued:
            % - 'center' The training point closest to the geometrical center
            % if the training data's bounding box is used.
            % - 't0' Tries to find a training sample for time zero and uses
            % that. If None are found, the strategy falls back to 'center'. If
            % multiple are found, the one closest to the geometrical center is
            % chosen.
            %
            % @note Might be overridden in subclasses or even exported as
            % strategy pattern if more variants show up.
            %
            % Parameters:
            % atd: The approximation training data @type data.ApproxTrainData
            %
            % Return values:
            % c: The initial center vector `[x; t; \mu]` @type colvec
            % idx: The index of the vector inside the training data @type integer
            %
            % See also: data.ApproxTrainData InitialCenter
            [c, idx] = atd.getClosestToCenter;
            if strcmp(this.InitialCenter,'center')
                return;
            elseif strcmp(this.InitialCenter,'t0')
                if atd.hasTime
                    idx = find(atd.ti == 0);
                    if ~isempty(idx)
                        c = [atd.xi(:,idx); zeros(1,numel(idx))];
                        if atd.hasParams
                            c = [c; atd.mui(:,idx)];
                        end
                        if numel(idx) > 1
                            [~,i] = min(sum((c-repmat(atd.Center,1,numel(idx))).^2,1));
                            c = c(:,i);
                            idx = idx(i);
                        end
                    end
                end
            else
                error('Unknown strategy: %s',this.InitialCenter);
            end
        end
        
        function [val, idx, errs] = getError(this, kexp, atd)
            % Computes the error according to the chosen error function
            % (see ErrFun property) with respect to the current kernel
            % expansion and the `f(x_i)` values in the training data.
            %
            % Parameters:
            % kexp: The kernel expansion @type kernels.KernelExpansion
            % atd: The approximation training data @type data.ApproxTrainData
            %
            % Return values:
            % val: The maximum error `\max ||f(x_i,t_i,\mu_i) -
            % \hat{f}(x_i,t_i,\mu_i)||_{\{2,\infty\}}` @type double
            % idx: The index of the maximum error inside the errs vector @type integer
            % errs: A row vector of the errors for each sample @type rowvec
            %
            % See also: ErrFun
            if this.pte
                afxi = kexp.evaluate(atd.xi, atd.ti, atd.mui);
            else
                afxi = kexp.evaluate(atd.xi);
            end
            if this.ErrFun == 1
                errs = max(abs(atd.fxi-afxi),[],1); % L^\infty error function
            else
                errs = sqrt(sum((atd.fxi-afxi).^2,1)); % L^2 error function
            end
            [val, idx] = max(errs);
        end
        
        function extendExpansion(this, kexp, atd, idx)
            % Extends the kernel expansion 'kexp' by the training data
            % sample in 'atd' given by index 'idx'.
            %
            % Parameters:
            % kexp: The kernel expansion @type kernels.KernelExpansion
            % atd: The approximation training data @type data.ApproxTrainData
            % idx: The index of the center to choose @type integer
            kexp.Centers.xi(:,end+numel(idx)) = atd.xi(:,idx);
            if this.pte
                if atd.hasTime
                    kexp.Centers.ti(end+numel(idx)) = atd.ti(idx);
                end
                if atd.hasParams
                    kexp.Centers.mui(:,end+numel(idx)) = atd.mui(:,idx);
                end
            end
        end
            
        function dists = getDists(this, atd, xm, tm, mum)
            % Computes the distances for the different Gaussian kernel
            % Gamma configurations using the 'atd' data and the algorithms
            % configuration.
            %
            % Parameters:
            % atd: The approximation training data @type data.ApproxTrainData
            % xm: The minimum distance for the state space. @type double @default data.ApproxTrainData.xiDia
            % tm: The minimum time distance. @type double @default data.ApproxTrainData.tiDia
            % mum: The minimum distance for the parameters space. @type double @default data.ApproxTrainData.muiDia
            %
            % Return values:
            % dists: A `3\times n` matrix with `\gamma` values for state,
            % time and parameter kernels (time and parameter if given, but
            % always 2nd and 3rd rows, respectively)
            dfun = @logsp; % gamma distances comp fun (linsp / logsp)
            
            if nargin == 2
                xm = atd.xiDia;
                if atd.hasTime
                    tm = atd.tiDia;
                end
                if atd.hasParams
                    mum = atd.muiDia;
                end
            end
            
            Mf = this.MaxGFactor;
            mf = this.MinGFactor;
            if isscalar(mf)
                mf = ones(3,1)*mf;
            end
            if isscalar(Mf)
                Mf = ones(3,1)*Mf;
            end
            dists = dfun(mf(1)*xm, Mf(1)*atd.xiDia);
            if atd.hasTime
                dists(2,:) = dfun(mf(2)*tm, Mf(2)*atd.tiDia);
            end
            if atd.hasParams
                dists(3,:) = dfun(mf(3)*mum, Mf(3)*atd.muiDia);
            end
            
            function d = linsp(from, to)%#ok
                d = linspace(from,to,this.NumGammas);
            end
            
            function d = logsp(from, to)
                d = logspace(log10(from),log10(to),this.NumGammas);
            end
        end
        
        function g = setDistKernelConfig(this, kexp, value)
            % Sets the configuration of the kernel expansion 'kexp' using
            % the 'dists' values and the algorithm configuration.
            %
            % Parameters:
            % kexp: The kernel expansion @type kernels.KernelExpansion
            % dists: A vector of dimension 3 containing the distances for
            % the state, time and parameter kernels. @type colvec
            %
            % Return values:
            % g: The effecively set `\gamma` values as 3-dim vector for the
            % state, time and parameter kernels.
            %
            % See also: gameps setKernelConfig
            g(1) = kexp.Kernel.setGammaForDistance(value(1),this.gameps);
            if KerMor.App.Verbose > 2
                fprintf('Kernel config: xg:%.5e',g(1));
            end
            if this.pte
                if ~isa(kexp.TimeKernel,'kernels.NoKernel')
                    g(2) = kexp.TimeKernel.setGammaForDistance(value(2),this.gameps);
                    if KerMor.App.Verbose > 2
                        fprintf(', tg:%.5e',g(2));
                    end
                end
                if ~isa(kexp.ParamKernel,'kernels.NoKernel')
                    g(3) = kexp.ParamKernel.setGammaForDistance(value(3),this.gameps);
                    if KerMor.App.Verbose > 2
                        fprintf(', pg=%.5e',g(3));
                    end
                end
            end
            if KerMor.App.Verbose > 2
                fprintf('\n');
            end
        end
        
        function dummy = setPolyKernelDegs(this, kexp, value)
            % Sets the configuration of the kernel expansion 'kexp' using
            % the 'value' values and the algorithm configuration.
            %
            % Parameters:
            % kexp: The kernel expansion @type kernels.KernelExpansion
            % dists: A vector of dimension 3 containing the degrees for
            % the state, time and parameter kernels. @type colvec
            kexp.Kernel.Degree = value(1);
            if this.pte
                if isa(kexp.TimeKernel,'kernels.PolyKernel')
                    kexp.TimeKernel.Degree = value(2);
                end
                if isa(kexp.ParamKernel,'kernels.PolyKernel')
                    kexp.ParamKernel.Degree = value(3);
                end
            end
            dummy = 0;
        end
        
        function setKernelConfig(this, kexp, g)
            % Sets the configuration of the kernel expansion 'kexp' using
            % the 'g' values as `\gamma` values for the state, time and
            % parameter kernel.
            %
            % Parameters:
            % kexp: The kernel expansion @type kernels.KernelExpansion
            % g: A vector of dimension 3 containing the `\gamma` values for
            % the state, time and parameter kernels. @type colvec
            %
            % See also: setDistKernelConfig
            kexp.Kernel.Gamma = g(1);
            if this.pte
                if ~isa(kexp.TimeKernel,'kernels.NoKernel')
                    kexp.TimeKernel.Gamma = g(2);
                end
                if ~isa(kexp.ParamKernel,'kernels.NoKernel')
                    kexp.ParamKernel.Gamma = g(3);
                end
            end
        end
        
        function bool = checkStop(this, cnt, rel, val)
            % Checks the stopping conditions for the adaptive approximation
            % algorithm.
            % Considers maximum expansion size, maximum relative and
            % absolute error (absolute error is computed as
            % 'MaxAbsErrFactor'`\times \max |fxi(:)|`)
            %
            % See also: MaxExpansionSize MaxRelErr MaxAbsErrFactor
            
            % Update lasterrs error improvement record
            this.lasterrs(end+1) = val;
            if numel(this.lasterrs) > this.ImproveRange
                this.lasterrs(1) = [];
            end
            reqimpr = mean(this.lasterrs(1:end-1))*(1-this.MinImprovePerc);
            
            bool = false;
            if cnt == this.MaxExpansionSize
                fprintf('AdaptiveCWKA stopping criteria holds: Max expansion size %d reached.\n',this.MaxExpansionSize);
                bool = true;
            elseif rel < this.MaxRelErr
                fprintf('AdaptiveCWKA stopping criteria holds: Relative error %.7e < %.7e\n',rel,this.MaxRelErr);
                bool = true;
            elseif val < this.effabs
                fprintf('AdaptiveCWKA stopping criteria holds: Absolute error %.7e < %.7e\n',val,this.effabs);
                bool = true;
            elseif numel(this.lasterrs) == this.ImproveRange && reqimpr < val
                fprintf('AdaptiveCWKA stopping criteria holds: Error improvement over mean error of last %d iterations below %2.2f%% percent (required:%e, achieved:%e)\n',...
                    this.ImproveRange, this.MinImprovePerc*100, reqimpr, val);
                bool = true;
                this.lasterrs = [];
            end
        end
    end
    
    methods(Abstract, Access=protected)
        % Runs the actual detailed algorithm.
        %
        % Parameters:
        % kexp: The kernel expansion @type kernels.KernelExpansion
        % atd: The approximation training data @type data.ApproxTrainData
        detailedAdaptiveApproximation(this, kexp, atd);
    end
    
    %% Getter & Setter
    methods
%         function set.ValidationPercent(this, value)
%             if ~isposrealscalar(value) || value > .5
%                 error('The value must be a positive scalar inside the interval ]0,.5[');
%             end
%             this.ValidationPercent = value;
%         end
        
        function set.NumGammas(this,value)
            if ~isposintscalar(value)
                error('Value must be a positive integer.');
            end
            this.NumGammas = value;
        end
        
        function set.MaxExpansionSize(this, value)
            if ~isposintscalar(value)
                error('Value must be a positive integer.');
            end
            this.MaxExpansionSize = value;
        end
                              
        function set.gameps(this, value)
            if ~isposrealscalar(value)
                error('The value must be a positive scalar');
            end
            this.gameps = value;
        end
        
        function set.MaxRelErr(this, value)
            if ~isposrealscalar(value)
                error('The value must be a positive scalar');
            end
            this.MaxRelErr = value;
        end
        
        function set.MaxAbsErrFactor(this, value)
            if ~isposrealscalar(value)
                error('The value must be a positive scalar');
            end
            this.MaxAbsErrFactor = value;
        end
        
        function set.ErrFun(this, value)
            if ~(value == 1 || value == 2)
                error('Value must be either integer 1 or 2.');
            end
            this.ErrFun = value;
        end
        
        function set.ImproveRange(this, value)
            if ~isempty(value) && (round(value) ~= value || value <= 0)
                error('ImproveRange must be a positive natural number.');
            end
            this.ImproveRange = value;
        end
        
        function set.MinImprovePerc(this, value)
            if ~isempty(value) && (value <= 0 || ~isscalar(value))
                error('MinImprovePerc must be a positive double scalar.');
            end
            this.MinImprovePerc = value;
        end
    end
end


