classdef BaseCompWiseKernelApprox < approx.BaseApprox & dscomponents.CompwiseKernelCoreFun
    %Base class for component-wise kernel approximations
    %
    % For each dimension `k` there is a representation
    % ``f_k(x) = \sum\limits_{i=1}^N \alpha_{k,i}\Phi(x,x_i) + b_k``
    % for the approximation. The property cData contains in row `k` all
    % indices `\alpha_k` used in the `k`-th dimension. off contains all
    % `b_k`.
    %
    % The property snData contains all state variable snData that
    % are relevant for the evaluation of the approximated function. (No
    % matter how many originally have been used for the approximation
    % computation!)
    %
    % See also: BaseKernelApprox
       
    methods(Access=protected)
        
        function gen_approximation_data(this, xi, ti, mui, fxi)
            % Computes the approximation according to the concrete
            % approximation strategy.
            % Fills the Ma, off and snData properties of the
            % CompwiseKernelCorefun with data.
            
            this.snData.xi = xi;
            this.snData.ti = ti;
            this.snData.mui = mui;
            n = size(xi,2);
            
            % Call subclass preparation method
            K = this.evaluateAtCenters(xi, ti, mui);
            this.prepareApproximationGeneration(K);
            
            % Make hermetian (no rounding errors)
            %ls.K = .5*(ls.K'+ls.K);
            try
                wh = waitbar(0,'Initializing component-wise kernel approximation');
                fdims = size(fxi,1);
                this.Ma = zeros(fdims, n);
                this.off = zeros(fdims, 1);
                for fdim = 1:fdims
                    waitbar(fdim/fdims+10,wh,sprintf('Computing approximation for dimension %d/%d ... %2.0f %%',fdim,fdims,(fdim/fdims)*100));
                    
                    % Call template method
                    [ai, b, svidx] = this.calcComponentApproximation(fxi(fdim,:));
                    if ~isempty(svidx)
                        this.Ma(fdim,svidx) = ai;
                    else
                        this.Ma(fdim,:) = ai;
                    end
                    this.off(fdim) = b;
                end
                
                waitbar(fdims+5/fdims+10,wh,'Compiling approximation model data...');
                
                % Kick out any too small weights
                %this.Ma(abs(this.Ma) < this.AlphaMinValue) = 0;
                
                % Reduce the snapshot array and coeff data to the
                % really used ones! This means if any snapshot x_n is
                % not used in any dimension, it is kicked out at this
                % stage.
                hlp = sum(this.Ma ~= 0,1);
                usedidx = find(hlp > 0);
                if length(usedidx) < n
                    this.Ma = this.Ma(:,usedidx);
                    this.snData.xi = xi(:,usedidx);
                    if ~isempty(ti)
                        this.snData.ti = ti(:,usedidx);
                    end
                    if ~isempty(mui)
                        this.snData.mui = mui(:,usedidx);
                    end
                end
                
                % @todo find out when sparse representation is more
                % efficient!
                if sum(hlp) / numel(this.Ma) < .5
                    this.Ma = sparse(this.Ma);
                end
                
                % dont use offset vector if none are given
                if all(this.off == 0)
                    this.off = [];
                end
                
                close(wh);
            catch ME
                close(wh);
                rethrow(ME);
            end
            
        end
        
        function target = clone(this, target)
            % Clones the instance.
            %
            % Note:
            % Since we use multiple inheritance we have to call the clone
            % method from both above superclasses. In this case this leads
            % to a double execution of the clone method from
            % dscomponents.ACoreFcn, but this is rather done than ommitted
            % and causing trouble later on.
            target = clone@dscomponents.CompwiseKernelCoreFun(this, target);
            target = clone@approx.BaseApprox(this, target);
            % No local properties to copy.
        end
    end
    
    methods(Abstract, Access=protected)
        % Preparation template method.
        %
        % Is called once before the component-wise approximation
        % computation. Here subclasses may create the conrete
        % single-dimension approximation algorithm implementations
        %
        % Parameters:
        % K: The Kernel matrix created from the snData `x_i`
        %
        % See also: calcComponentApproximation
        prepareApproximationGeneration(this, K);
        
        % Single dimension approximation computation.
        %
        % Here the concrete class performs the approximation calculation
        % for given function evaluation points `fx_i` at the snData
        % `x_i` also used to compute the kernel matrix for
        % prepareApproximationGeneration. Within this method the class
        % properties cData and csMap are to be
        % manipulated/filled with approximation data. cData ist
        % preallocated to dims x size(xi,2).
        %
        % Parameters:
        % fxi: The function values `f(x_i)` as row vector.
        %
        % Return values:
        % ai: The coefficients `\alpha_{k,i}` of `f_k(x)`.
        % b: The offset for `f_k(x)`.
        % svidx: The used support vector indices `i` of `x_i`. Optional,
        % leave empty if all are used.
        %
        % See also: prepareApproximationGeneration
        [ai, b, svidx] = calcComponentApproximation(this, fxi);
    end
    
end

