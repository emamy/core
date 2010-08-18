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
    
%     properties
%         % Minimum `\alpha` coefficient value after projection
%         %
%         % At projection, the coefficients for each component approximation
%         % get mixed leading to an enlarged set of support vectors. Thus,
%         % after combination of the coefficients each new one must be
%         % greater than this threshold to be considered a new one.
%         AlphaMinValue = 1e-5;
%     end
    
%     properties(SetAccess=protected)
%         % The coefficient data for each dimension.
%         Ma;
%     end
%     
%     properties(Access=protected)
%         % The state variable Snapshots used in the approximation.
%         %
%         % This is the union of all snData used within the approximation.
%         snData;
%         
%         % The `b_k` offsets for each dimension.
%         %
%         % Set to empty if no off are used. This property is
%         % preinitialized to [].
%         off = [];
%     end
    
    methods
        
%         function this = BaseCompWiseKernelApprox
%             % Class constructor.
%             %
%             % Sets the CustomProjection property from BaseApprox to true as
%             % component-wise approximations implement their own projection.
%             this.CustomProjection = true;
%         end
        
%         function phi = getPhi(this, zs, t, mu)
%             % Gets the phi(Vz(s),t,mu) value used in error estimations.
%             % @docupdate
%             phi = this.evaluateAtCenters(this.snData.xi, this.snData.ti, this.snData.mui, zs, t, mu);
%         end
        
%         function c = getGlobalLipschitz(this, t, mu)
%             % @todo validate computation
% %             [x,ti,mui] = this.splitTripleVect(this.snData);
%             k = abs(this.TimeKernel.evaluate(this.snData.ti,t).*this.ParamKernel.evaluate(this.snData.mui,mu));
%             c = 0;
%             for idx = 1:size(this.Ma,2)
%                 c = c + k(idx)*norm(this.Ma(:,idx));
%             end
%             c = this.sk.getGlobalLipschitz * c;
%             warning('some:id','not yet implemented/validated correctly!');
%         end
        
    end
    
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
        
%         function fx = evaluate_approximation(this, x, t, mu)
%             d=this.snData;
%             if ~isempty(this.off)
%                 fx = this.Ma * this.evaluateAtCenters(d.xi, d.ti, d.mui, x, t, mu) + repmat(this.off,1,size(x,2));
%             else
%                 fx = this.Ma * this.evaluateAtCenters(d.xi, d.ti, d.mui, x, t, mu);
%             end
%         end
    
%         function copy = customProject(this, V, W)
%             % Performs specific projection computation for component wise
%             % svr.
%             
%             % Create copy to preserve settings in original model
%             copy = this.clone;
%             
%             copy.Ma = W' * this.Ma;
%             if ~isempty(this.off)
%                 copy.off = W' * this.off;
%             end
% 
%             % Project snapshot vectors into reduced space, in case the used
%             % Kernel is rotation invariant (lossless projection, no change
%             % of coefficients necessary)
%             if this.RotationInvariantKernel
%                 % Extract system part and project into V space
% %                 [x,t,mu] = this.splitTripleVect(this.snData);
%                 copy.snData.xi = W' * this.snData.xi;
% %                 copy.snData = this.compileTripleVect(x,t,mu);
%             end
%         end
        
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

