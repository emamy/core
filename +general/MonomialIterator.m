classdef MonomialIterator < handle
% MonomialIterator: Create list of d-dim monomials
%
% For a given dimension, this class produces an infinite list of monomials which increase in
% degree as all monomials of the last degree have been listed.
% The natural ordering of the monomials is degree first, then lexographically.
% Example for d=3: `0,z,x,y,zz,zy,zx,yy,yx,xx,zzz,zzy,zzx,zyy,zyx\ldots`
%
% Note that the zero-monomial is not included in the list, however, it can be accessed via the
% extra method getNullMonomial.
%
% @author Daniel Wirtz @date 2012-01-20
%
% @new{0,6,dw,2012-01-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetAccess=private)
        Dim;
    end
    
    properties(Access=private)
        incidx;
        deg;
        rnd;
        expidx;
    end
      
    methods
        function this = MonomialIterator(dim)
            % Creates a new monomial iterator.
            %
            % Parameters:
            % dim: The target dimension `d\in\N` @type integer
            if int32(dim) ~= dim || dim <= 0
                error('Dimension d must be a positive integer.');
            end
            this.Dim = int32(dim);
            this.deg = int32(1);
            this.rnd = RandStream('mt19937ar','Seed',2);
            this.expidx = int32(1);
            this.incidx = int32(1);
        end
        
        function alpha = getNullMonomial(this)
            % Get the monomial of degree zero
            %
            % Return values:
            % alpha: The multiindex of size `1\times d` with zero entries, corresponding to the
            % zero degree monomial. @type colvec<double>
            alpha = zeros(this.Dim,1);
        end
        
        function [alpha, deg, exps] = nextMonomial(this)
            % Returns the next monomial in the sequence.
            %
            % Return values:
            % alpha: The multiindex `\alpha\in\Z^d` corresponding to the exponents of the
            % current monomial @type colvec<double>
            % deg: The degree `p = |\alpha|` of the monomial @type int32
            % exps: (Debug) The indices of `(1 \ldots d)` at which an exponent is to be added.
            % Connection: A histogram search of \b exps with bins `(1 \ldots d)` results in
            % `\alpha` @type rowvec<int32>

            % Collect exponents via histogram count for each dimension
            alpha = fliplr(histc(this.expidx, 1:this.Dim))';
            %alpha = histc(this.expidx, 1:this.Dim)';
            deg = this.deg;
            exps = this.expidx;
            
            % Increment exponentials at current incidx
            if this.expidx(this.incidx) == this.Dim && this.incidx > 1
                % Set increase counter one to the left
                this.incidx = this.incidx - 1;
                
                % Einen hoch
                this.expidx(this.incidx) = this.expidx(this.incidx) + 1;
                
                % See 
                if this.expidx(this.incidx) < this.Dim
                    % Daneben alle auf den Wert setzen
                    this.expidx(this.incidx+1:this.deg) = this.expidx(this.incidx);
                    this.incidx = this.deg;
                end
            else
                this.expidx(this.incidx) = this.expidx(this.incidx) + 1;
            end
           
            % Check if current degree is finished
            if this.incidx == 0 || this.expidx(this.incidx) > this.Dim %equiv: all(this.expidx == this.Dim)
                this.deg = this.deg + 1;
                this.expidx = int32(ones(1,this.deg));
                % Increment index is last exponent
                this.incidx = this.deg;
            end
        end
        
        function alpha = getRandMonomial(deg)
            % Returns a random monomial.
            %
            % Parameters:
            % deg: The degree `p = |\alpha|` of the monomial. @type int32
            %
            % Return values:
            % alpha: The multiindex `\alpha\in\Z^d` corresponding to the exponents of the
            % current monomial @type colvec<int32>
            alpha = int32(zeros(this.Dim,1));
            s = 0;
            while s < deg
                idx = this.rnd.randi(this.Dim);
                alpha(idx) = alpha(idx) + this.rnd.randi(deg-s);
                s = sum(alpha);
            end
        end
    end
    
    methods(Static)
        function res = test_MonomialIterator
            % Tests the MonomialIterator for two test cases.
            %
            % Return values:
            % res: True if successful, false else @type logical
            for Dim = 1:3 %#ok<*PROP>
                fprintf('Dim=%d:\n',Dim)
                mi = general.MonomialIterator(Dim);
                deg = 0;
                while deg < 5
                    [alpha, deg, exps] = mi.nextMonomial;
                    fprintf('p=%d, exps=[%s], alpha=[%s]\n',deg, general.Utils.implode(exps,',','%d'),...
                        general.Utils.implode(alpha,',','%d'));
                end
            end
            res = true;
            % test for d=1 and p=1000
            mi = general.MonomialIterator(1);
            for i=1:1000
                [alpha, deg, exps] = mi.nextMonomial; %#ok<*PROP>
                res = res && alpha == i && deg == i && all(exps==1);
            end
            
            % test for d=4 and p=3
            mi = general.MonomialIterator(4);
            cmpp1 = [1; 2; 3; 4];
            cmpp2 = [1 1; 1 2; 1 3; 1 4; 2 2; 2 3; 2 4; 3 3; 3 4; 4 4];
            cmpp3 = [1 1 1; 1 1 2; 1 1 3; 1 1 4; 1 2 2; 1 2 3; 1 2 4; 1 3 3; 1 3 4; 1 4 4; 2 2 2;...
                2 2 3; 2 2 4; 2 3 3; 2 3 4; 2 4 4; 3 3 3; 3 3 4; 3 4 4; 4 4 4];
            for i = 1:size(cmpp1,1)
                [~, deg, exps] = mi.nextMonomial;
                res = res && deg == 1 && all(exps==cmpp1(i,:));
            end
            for i = 1:size(cmpp2,1)
                [~, deg, exps] = mi.nextMonomial;
                res = res && deg == 2 && all(exps==cmpp2(i,:));
            end
            for i = 1:size(cmpp3,1)
                [~, deg, exps] = mi.nextMonomial;
                res = res && deg == 3 && all(exps==cmpp3(i,:));
            end
        end
    end
    
end