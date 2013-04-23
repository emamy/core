classdef FileDataCollection < data.FileData
% FileDataCollection: Basic class for storing data given a hashable key value
%
%
% @author Daniel Wirtz @date 2012-09-19
%
% @new{0,6,dw,2012-09-19} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        % The HashMap used to store the indices for each trajectory.
        hm;
    end
    
    methods
        function this = FileDataCollection(data_dir)
            % Creates a new FileDataCollection in a file folder.
            %
            % Parameters:
            % data_dir: Either a data.ModelData instance to infer the storage root from, or a
            % string containing a valid folder. @default A temporary folder within the
            % KerMor.TempDirectory
            if ~usejava('jvm')
                error('FileTrajectoryData cannot be used as java is not enabled.');
            end
            if nargin < 1
                data_dir = fullfile(KerMor.App.TempDirectory,...
                    sprintf('temp_fdc_%s',...
                    IDGenerator.generateID));
            end
            this = this@data.FileData(data_dir);
            this.hm = java.util.HashMap;
            this.clear;
        end
        
        function delete(this)
            % Destructor for FileTrajectoryData
            %
            % Deletes the DataDirectory if no trajectories are stored in it or it has not been
            % saved somewhere.
            if ~this.isSaved || this.hm.size == 0
                this.clear;
            end
            % Superclass delete removes the folder if empty.
            delete@data.FileData(this);
        end
    end
    
    methods(Access=protected)
        function data = getData(this, keydata, varargin)
            % Retrieves data from the collection element of a given key.
            %
            % Parameters:
            % keydata: The unique key data to identify a collection object @type rowvec<double>
            % varargin: Contains the name of the fields to obtain as strings
            %
            % Return values:
            % data: A struct with the fields as specified in varargin. @type struct
            %
            % @todo make varargout instead of data struct
            key = Utils.getHash(keydata);
            data = [];
            if this.hm.containsKey(key)
                try
                    data = load(this.getfile(this.hm.get(key)), varargin{:});
                catch ME
                    warning('KerMor:FileDataCollection',...
                        'Error loading a file from directory %s: %s',...
                        this.DataDirectory,ME.message);
                    data = [];
                end
            else
                % "Backup" function. In case some data is stored and then the
                % FileDataCollection is not saved, the existing trajectory files are recognized
                % and loaded.
                file = fullfile(this.DataDirectory,[key '.mat']);
                if exist(file,'file') == 2
                    data = load(file,varargin{:});
                    this.hm.put(key,[key '.mat']);
                end
            end
        end
        
        function n = getCollectionSize(this)
           n = this.hm.size;
        end
        
        function data = getDataNr(this, nr, varargin)
            % Retrieves data from the nr-st collection element.
            %
            % Parameters:
            % nr: The element number
            % varargin: Contains the name of the fields to obtain as strings
            %
            % Return values:
            % data: A struct with the fields as specified in varargin. @type struct
            %
            % @todo make varargout instead of data struct
            if nr > this.hm.size || nr < 1
                error('Invalid data position number: %d',nr);
            end
            keys = this.hm.keySet.toArray(java_array('java.lang.String',1));
            data = load(this.getfile(this.hm.get(keys(nr))), varargin{:});
        end
        
        function res = hasData(this, keydata)
            res = this.hm.containsKey(Utils.getHash(keydata));
        end
        
        function addData(this, keydata, data)
            key = Utils.getHash(keydata);
            file = [key '.mat'];
            this.hm.put(key,file);
            
            file = fullfile(this.DataDirectory,file);
            try
                fn = fieldnames(data);
                save(file,'-struct','data',fn{:});
            catch ME
                this.hm.remove(key);
                rethrow(ME);
            end
        end
        
        function clear(this)
            ks = this.hm.values.iterator;
            while ks.hasNext
                try
                    file = this.getfile(ks.next);
                catch ME
                    warning('KerMor:data:FileDataCollection',...
                        'Could not delete file "%s": %s\nPlease remove manually.',...
                        file,ME.message);
                end
                delete(file);
            end
            this.hm.clear;
        end
    end
   
    methods(Static, Access=protected)
        function this = loadobj(this, initfrom)
            % Loads a FileTrajectoryData instance.
            %
            % Ensures that the directory associated with this FileTrajectoryData is existent.
            
            created = false;
            if ~isa(this, 'data.FileDataCollection')
                initfrom = this;
                this = data.FileDataCollection(initfrom.DataDirectory);
                created = true;
            end
            if nargin == 2 || created
                this.hm = initfrom.hm;
                this = loadobj@data.FileData(this, initfrom);
            else
               this = loadobj@data.FileData(this);
            end
        end
    end
    
    methods(Static)
        function res = test_FileDataCollection
            % Quickly tests if the directory is persisted(deleted) when (not) saved to a mat
            % file
            res = true;
            thedir = fullfile(KerMor.App.TempDirectory,'testFDC');
            c = data.FileDataCollection(thedir);
            num = 10;
            keys = rand(num,1);
            
            for k=1:num
                s.field = rand(200,300);
                c.addData(keys(k),s);
            end
            clear c;
            
            if exist(thedir,'file') ~= 0
                res = false;
            end
            
            c = data.FileDataCollection(thedir);
            for k=1:num
                s.field = rand(500,500);
                c.addData(keys(k),s);
            end
            save(fullfile(thedir,'tmpsave.mat'),'c');
            clear c;            
            if exist(thedir,'file') ~= 7
                res = false;
            end
            rmdir(thedir,'s');
        end
    end
    
end