classdef LinearSplitOfOne < handle
% LinearSplitOfOne: Computes a sequence of hat functions at equidistant
% nodes from [0,len] to enable an efficient, easy way of division of unity.
%
%
% @author Daniel Wirtz @date 2014-04-11
%
% @new{0,7,dw,2014-04-11} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The variable name in the function strings
        %
        % @type char @default 'x'
        VarName = 'x';
        
        % The number of splits to do. Results in +1 functions, indexed
        % starting at 0.
        %
        % @type integer @default 1
        Splits = 1;
        
        % The length of the interval over which to split unity.
        %
        % @type double @default 1
        IntervalLength = 1;
    end
    
    methods
        
        function this = LinearSplitOfOne(varname, splits, len)
            % Creates a new instance
            %
            % Parameters:
            % varname: The variable name @type char @default 'x'
            % splits: The number of splits to perform @type integer @default 1
            % len: Splitting interval length @type double @default 1
            if nargin > 0
                this.VarName = varname;
                if nargin > 1
                    this.Splits = splits;
                    if nargin > 2
                        this.IntervalLength = len;
                    end
                end
            end
        end
        
        function str = getFunStr(this, idx)
            % Returns the idx-th function string.
            %
            % Pass idx=0 for the function falling from 1 to 0 for arguments
            % starting at 0.
            if idx > this.Splits || idx < 0
                error('Max %d functions (starting at 0) for %d splits',this.Splits+1,this.Splits);
            end
            ilen = this.IntervalLength / this.Splits;
            if idx == 0
                str = sprintf('(1-@@varname@@/%g).*(@@varname@@<%g)',ilen,ilen);
            elseif idx == this.Splits
                str = sprintf('((@@varname@@-%g)/%g).*(@@varname@@>=%g)',(idx-1)*ilen,ilen,(idx-1)*ilen);
            else
                % zero to one
                zto = sprintf('(@@varname@@>=%g).*(@@varname@@<%g).*((@@varname@@-%g)/%g)',...
                    (idx-1)*ilen,idx*ilen,(idx-1)*ilen,ilen);
                % one to zero
                otz = sprintf('(@@varname@@>=%g).*(@@varname@@<%g).*(1-(@@varname@@-%g)/%g)',...
                    idx*ilen,(idx+1)*ilen,idx*ilen,ilen);
                str = [zto '+' otz];
            end
            str = strrep(str,'@@varname@@',this.VarName);
        end
        
        function funs = getAllFunctions(this)
            % Returns all functions in a cell array
            %
            % Return values:
            % funs: The function strings in a cell array @type cell
            funs = cell(this.Splits+1,1);
            for idx = 1:this.Splits+1
                funs{idx} = this.getFunStr(idx-1);
            end
        end
        
        function handle = getAllFunctionsParsed(this)
            % Returns a function handle for a vector function, where each
            % component idx evaluates the (idx-1)-th function of
            % getFunStr(idx)
            funs = this.getAllFunctions;
            str = ['@(' this.VarName ')['];
            for idx = 1:this.Splits+1
                str = [str funs{idx} ';'];%#ok
            end
            str = [str '];'];
            handle = eval(str);
        end
        
        function set.Splits(this, value)
            if value < 1
                error('Minimum of one split necessary.');
            end
            this.Splits = value;
        end
    end
    
    methods(Static)
        function res = test_LinearSplitOfOne
            num = 20;
            s = randi(50,1,num);
            l = rand(1,num);
            res = true;
            for k=1:num
                lso = LinearSplitOfOne('x',s(k),l(k));
                f = lso.getAllFunctionsParsed;
                
                pos = linspace(0,l(k),200);
                for i=1:length(pos)
                    if i < 200
                        res = res && sum(f(pos(i))) == 1;
                    else
                        fprintf('Splits: %d, len: %g, diff at end: %g\n',s(k),l(k),sum(f(pos(i)))-1);
                    end
                end
            end
        end
    end
    
end