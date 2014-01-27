classdef ImprovedLocalSecantLipschitz < error.lipfun.Base
% ImprovedLocalSecantLipschitz: 
%
%
% @author Daniel Wirtz @date 2011-05-20
%
% @change{0,4,dw,2011-06-07} Moved the ModifiedNewton methods from this class to kernels.BellFunction as
% they are more appropriate at that place.
%
% @new{0,4,dw,2011-05-31} Added new prepareConstants init function from Base.
%
% @new{0,4,dw,2011-05-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(SetObservable)
        % @propclass{optional}
        % @todo make docs, make dependent and re-generate if set to true when not used..
        UseCachedSecants = false;
    end
    
    properties(Access=private, Transient)
        oldrs = [];
    end

    properties(SetAccess=private)
        precomp = [];
    end
    
    methods
        
        function this = ImprovedLocalSecantLipschitz(bellfcn)
            this = this@error.lipfun.Base(bellfcn);
            this.registerProps('UseCachedSecants');
        end
        
        function copy = clone(this)
            copy = error.lipfun.ImprovedLocalSecantLipschitz(this.bellfcn);
            %copy = clone@error.lipfun.Base(this, copy);
            copy.UseCachedSecants = this.UseCachedSecants;
        end
        
        function prepareConstants(this)%#ok
            % nothing to do here!
        end
        
        function ci = evaluate(this, di, C)
            b = this.bellfcn;
            r0 = b.r0;
            if isempty(this.oldrs) || any(isnan(this.oldrs)) || size(this.oldrs,2) ~= size(di,2)
                this.oldrs = r0+.2*sign(r0-di);
            end
            
            % Consider C=Inf case separately for speed reasons
            if isinf(C)
                if this.UseCachedSecants
                    rs = this.CachedModifiedNewton(di);
                else
                    rs = b.ModifiedNewton(this.oldrs,di);
                end
                ci = abs(b.evaluateD1(rs));
                this.oldrs = rs;
            else
                rs = ones(size(di))*r0;
                update = abs(di-r0) < C;
                % Choose suitable starting conditions if no old rs vectors
                % are available
                if any(update)
                    if this.UseCachedSecants
                        rs(update) = this.CachedModifiedNewton(di(update));
                    else
                        rs(update) = b.ModifiedNewton(this.oldrs(update),di(update));
                    end
                    this.oldrs(update) = rs(update);
                end
                left = di + C - rs < 0;
                right = di - C - rs > 0;
                center = ~left & ~right;
                % If C is too small we just take the derivative at this
                % spot
                if C < sqrt(eps)
                    ci(left | right) = abs(b.evaluateD1(di(left | right)));
                else
                    ci(left) = (b.evaluateScalar((di(left))) - b.evaluateScalar((di(left)+C))) / C;
                    ci(right) = (b.evaluateScalar((di(right)-C)) - b.evaluateScalar((di(right)))) / C;
                end
                ci(center) = abs(b.evaluateD1(rs(center)));
            end
        end
        
        function precompMaxSecants(this, maxDistance, num)
            if this.UseCachedSecants
                di = linspace(0,maxDistance,num);
                start = this.bellfcn.r0+.2*sign(this.bellfcn.r0-di);
                if KerMor.App.Verbose > 1
                    fprintf('Precomputing %d max local secants for BellFunction(Gamma=%f,r0=%f) in range [0, %f] ..\n',num,this.Gamma,this.x0,maxDistance);
                end
                rs = this.bellfcn.ModifiedNewton(start,di);
                t = BinTree;
                if KerMor.App.Verbose > 1
                    fprintf('Building tree data structure..\n');
                end
                t.Insert(di,rs);
                this.precomp = t;
            end
        end
    end
    
    methods(Access=private)
        
        function rs = CachedModifiedNewton(this, s)
            if ~isempty(this.precomp)
                t = this.precomp;
                [l,u] = t.FindClosest(s);
                low = s < this.bellfcn.r0;
                rs = zeros(size(y));
                rs(low) = u(low);
                rs(~low) = l(~low);
            else
                error('When using cached secants the method precompMaxSecants has to be run first.');
            end
        end
        
    end
    
    methods(Static)
        
    end
    
end