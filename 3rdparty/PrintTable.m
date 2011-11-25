classdef PrintTable < handle
% PrintTable: Class that allows table-like output spaced by tabs for multiple rows
%
% Allows to create a table-like output for multiple values and columns.
% A spacer can be used to distinguish different columns (PrintTable.ColSep) and a flag
% (PrintTable.HasHeader) can be set to insert a line of dashes after the first row.
% The table can be printed directly to the console or be written to a file.
%
% The cell content values can be anything that is straightforwardly parseable. You can pass
% char array arguments directly or numeric values; even function handles and classes (handle
% subclasses) can be passed and will automatically be converted to a string representation.
% However, if you need to have a customized string representation, at this stage this must be
% done in your code before calling addRow.
%
% Example:
% t = PrintTable;
% t.addRow('123','456','789');
% t.addRow('1234567','1234567','789');
% t.addRow('1234567','12345678','789');
% t.addRow('12345678','123','789');
% t.addRow('123456789','123','789');
% t.addRow('123456789','12345678910','789');
% t.addRow('adgag',uint8(4),4.6);
% t.addRow(@(mu)mu*4+56*245647869,t,'adgag');
% t.addRow('adgag',4,4.6);
% % Call display
% t.display;
% t.HasHeader = true;
% % or simply state the variable to print
% t
% t.ColSep = ' -@- ';
% t
%
% Printing the table:
% You can also print the table to a file. Any MatLab file handle can be used (any first
% argument for fprintf).
% % run the above example, then
% fid = fopen('mytable.txt','w');
% t.print(fid);
% fclose(fid);
% % this is actually equivalent to calling
% t.saveToFile('mytable.txt');
%
% @note Of course, any editor might have its own setting regarding tab spacing. As the default
% in MatLab and e.g. KWrite is 8 characters, this is what is used here. Change the TabCharLen
% constant to fit to your platform/editor/requirements.
% 
% See also: fprintf
%
% @author Daniel Wirtz @date 2011-11-17
%
% @new{0,6,dw,2011-11-17} Added this class.
%
% This class has originally been developed as part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
%
% Copyright (c) 2011, Daniel Wirtz
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are
% permitted only in compliance with the BSD license, see
% http://www.opensource.org/licenses/bsd-license.php
    
    properties(Constant)
        % Equivalent length of a tab character in single-space characters
        %
        % @default 8 @type integer
        TabCharLen = 8;
    end
    
    properties
        % A char sequence to separate the different columns.
        %
        % @default ' | ' @type char
        ColSep = ' | ';
        
        % Flag that determines if the first row should be used as table header.
        %
        % If true, a separate line with dashes will be inserted after the first printed row.
        %
        % @default false @type logical
        HasHeader = false;
    end
    
    properties(Access=private)
        % The string cell data
        data;
        % Number of min tabs for each column
        ntabs;
    end
    
    methods
        function this = PrintTable
            % Creates a new PrintTable instance.
            this.data = {};
            this.ntabs = [];
        end
        
        function display(this)
            % Overload for the default builtin display method.
            %
            % Calls print with argument 1, i.e. standard output.
            this.print(1);
        end
        
        function print(this, outfile)
            % Prints the current table to a file pointer.
            %
            % Parameters:
            % outfile: The file pointer to print to. Must either be a valid MatLab 'fileID'
            % or can be 1 or 2, for stdout or stderr, respectively. @type integer @default 1
            %
            % See also: fprintf
            if nargin == 1
                outfile = 1;
            end
            tabchar = sprintf('\t');
            for ridx = 1:length(this.data)
                row = this.data{ridx};
                totlen = 1;
                for i = 1:length(row)-1
                    str = row{i};
                    fillstabs = ceil((totlen+length(str))/PrintTable.TabCharLen);
                    tabs = sum(this.ntabs(1:i-1))+this.ntabs(i)-fillstabs+1;
                    fprintf(outfile,'%s%s',[str repmat(tabchar,1,tabs)],this.ColSep);
                    totlen = totlen + this.ntabs(i)*PrintTable.TabCharLen+length(this.ColSep);
                end
                fprintf(outfile,'%s\n',row{end});
                if ridx == 1 && this.HasHeader
                    fprintf(outfile,'%s\n',repmat('_',1,sum(this.ntabs)*PrintTable.TabCharLen));
                end
            end
        end
        
        function saveToFile(this, filename)
            % Prints the current table to a file.
            %
            % Parameters:
            % filename: The file to print to. If the file exists, any contents are discarded.
            % If the file does not exist, an attempt to create a new one is made. @type char
            fid = fopen(filename,'w');
            this.print(fid);
            fclose(fid);
        end
        
        function addRow(this, varargin)
            % Adds a row to the current table.
            %
            % Parameters:
            % varargin: Any number of arguments >= 1, each corresponding to a column of the
            % table. Each argument must be a char array.
            if isempty(varargin)
                error('Not enough input arguments.');
            end
            if isempty(this.data)
                this.data{1} = this.parse(varargin);
                this.ntabs = ones(1,length(varargin));
            else
                if length(this.data{1}) ~= length(varargin)
                    error('Inconsistent row length. Current length: %d, passed: %d',length(this.data{1}),length(varargin));
                end
                this.data{end+1} = this.parse(varargin);
            end
            %% Process tab numbers
            newrow = this.data{end};
            for i=1:length(newrow)
                t = ceil(length(newrow{i})/PrintTable.TabCharLen);
                if t > this.ntabs(i)
                    this.ntabs(i) = t;
                end
            end
        end
        
        function set.ColSep(this, value)
            if ~isempty(value) && ~isa(value,'char')
                error('ColSep must be a char array.');
            end
            this.ColSep = value;
        end
        
        function set.HasHeader(this, value)
            if ~islogical(value) || ~isscalar(value)
                error('HasHeader must be a logical scalar.');
            end
            this.HasHeader = value;
        end
    end
    
    methods(Access=private)
        function parsed = parse(this, data)
            parsed = cell(1,length(data));
            for i=1:length(data)
                el = data{i};
                if isa(el,'char')
                    parsed{i} = el;
                elseif isinteger(el)
                    parsed{i} = sprintf('%d',el);
                elseif isnumeric(el)
                    parsed{i} = sprintf('%e',el);
                elseif isa(el,'function_handle')
                    parsed{i} = func2str(el);
                elseif isa(el,'handle')
                    mc = metaclass(el);
                    parsed{i} = mc.Name;
                else
                    error('Cannot automatically parse an argument of type %s for PrintTable display.',class(el));
                end
            end
        end
    end
    
    methods(Static)
        function test_PrintTable
            t = general.PrintTable;
            t.addRow('123','456','789');
            t.addRow('1234567','1234567','789');
            t.addRow('1234567','12345678','789');
            t.addRow('12345678','123','789');
            t.addRow('123456789','123','789');
            t.addRow('123456789','12345678910','789');
            t.display;
        end
    end
    
end