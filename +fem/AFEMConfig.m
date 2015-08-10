classdef AFEMConfig < handle
    %AModelConfig
    
    properties(SetAccess=private)
        FEM;
        
        Geometry;
        
        Model;
    end
    
    properties(SetAccess=protected)
        % The coordinate system in which to interpret the applied pressure
        % of neumann boundary conditions.
        %
        % 'local' uses the normals on the faces as given in the reference
        % configuration, i.e. the "true" normals.
        %
        % 'global' uses the global coordinate system of the geometry, i.e.
        % the master element. This can be used to apply forces coming from
        % one fixed direction over a possibly noneven geometry surface.
        %
        % @type char @default 'local'
        NeumannCoordinateSystem = 'local';
        
        Options;
    end
    
    properties
        % Velocity conditions application function
        VelocityBCTimeFun;
    end
    
    properties(Access=private)
        iP;
        optArgs = {};
    end
    
    methods
        function this = AFEMConfig(varargin)
            this.optArgs = varargin;
            
            i = inputParser;
            i.KeepUnmatched = true;
            this.iP = i;
            
            this.addOption('GeoNr',1);
            mc = metaclass(this);
            this.addOption('Tag',mc.Name);
        end
        
        function configureModel(this, model)
            % Overload this method to set model-specific quantities like
            % simulation time etc
            this.Model = model;
        end
        
        function configureModelFinal(this)
        end
        
        function prepareSimulation(this, mu, inputidx)
            % Overload this method to initialize model-specific quantities
            % that are fixed for each simulation
            %
            % Called by override of computeTrajectory in muscle.Model
            
            % do nothing by default
        end
        
        function P = getBoundaryPressure(~, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is megaPascal [MPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            %
            % See also: NeumannCoordinateSystem
            P = [];
        end
        
        function str = getOptionStr(this, withtag)
            if nargin < 2
                withtag = true;
            end
            o = this.Options;
            fieldnames = fields(o);
            strs = {};
            for k=1:length(fieldnames)
                if (withtag || ~strcmp(fieldnames{k},'Tag')) && ~isempty(o.(fieldnames{k}))
                    fmt = '%s-%g';
                    if isa(o.(fieldnames{k}),'char')
                        fmt = '%s-%s';
                    end
                    strs{end+1} = sprintf(fmt,fieldnames{k},o.(fieldnames{k}));%#ok
                end
            end
            str = Utils.implode(strs,'_');
        end
        
        function plotGeometryInfo(this, allnode, elemnr)
            if nargin < 3
                elemnr = 1;
                if nargin < 2
                    allnode = false;
                end
            end
            g = this.FEM.Geometry;
            g.plot(allnode,elemnr);
        end
    end
    
    methods(Access=protected)
       
        function init(this)
            %% Parse the options
            this.iP.parse(this.optArgs{:});
            this.Options = this.iP.Results;
            
            %% Get the geometry
            geo = this.getGeometry;
            if isa(geo,'fem.geometry.Cube27Node')
                this.FEM = fem.HexahedronTriquadratic(geo);
            elseif isa(geo,'fem.geometry.Cube20Node')
                this.FEM = fem.HexahedronSerendipity(geo);
            elseif isa(geo,'fem.geometry.Cube8Node')
                this.FEM = fem.HexahedronTrilinear(geo);
            end
            this.Geometry = geo;
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(~, velo_dir, velo_dir_val)
            % Determines the dirichlet velocities.
            %
            % The unit for the applied quantities is [mm/ms] = [m/s]
            %
            % In the default implementation there are no velocity
            % conditions.
        end
    end
    
    methods(Sealed, Access=protected)
        function addOption(this, name, default, varargin)
            this.iP.addParamValue(name, default, varargin{:});
        end
    end
    
    methods(Abstract, Access=protected)
        displ_dir = setPositionDirichletBC(this, displ_dir);
        
        % Returns the intended geometry for this model config.
        %
        % The options will be set at call time, e.g. "GeoNr" is already set.
        geo = getGeometry(this);
    end
    
    methods(Sealed)
        function [displ_dir, velo_dir, velo_dir_val] = getBC(this)
            N = this.Geometry.NumNodes;
            displ_dir = false(3,N);
            displ_dir = this.setPositionDirichletBC(displ_dir);
            velo_dir = false(3,N);
            velo_dir_val = zeros(3,N);
            [velo_dir, velo_dir_val] = this.setVelocityDirichletBC(velo_dir, velo_dir_val);
            
            if any(any(displ_dir & velo_dir)) 
                error('Cannot impose displacement and velocity dirichlet conditions on same DoF');
            end
        end
        
        function [force, nodeidx, faceswithforce] = getSpatialExternalForces(this)
            fe = this.FEM;
            geo = fe.Geometry;
            ngp = fe.GaussPointsPerElemFace;
            force = zeros(geo.NumNodes * 3,1);
            faceswithforce = false(1,geo.NumFaces);
            
            globalcoord = strcmp(this.NeumannCoordinateSystem,'global');
            for fn = 1:geo.NumFaces
                elemidx = geo.Faces(1,fn);
                faceidx = geo.Faces(2,fn);
                masterfacenodeidx = geo.MasterFaces(faceidx,:);
                % So far: Constant pressure on all gauss points!
                P = this.getBoundaryPressure(elemidx, faceidx);
                if ~isempty(P)
                    faceswithforce(fn) = true;
                    integrand = zeros(3,geo.NodesPerFace);
                    if globalcoord
                        N = geo.FaceNormals(:,faceidx);
                    end
                    for gi = 1:ngp
                        if ~globalcoord
                            N = fe.NormalsOnFaceGP(:,gi,fn);
                        end
                        PN = (P * N) * fe.Ngpface(:,gi,fn)';
                        integrand = integrand + fe.FaceGaussWeights(gi)*PN*fe.face_detjac(fn,gi);
                    end
                    facenodeidx = geo.Elements(elemidx,masterfacenodeidx);
                    facenodeidx = (facenodeidx-1)*3+1;
                    facenodeidx = [facenodeidx; facenodeidx+1; facenodeidx+2];%#ok
                    force(facenodeidx(:)) = force(facenodeidx(:)) + integrand(:);
                end
            end
            % Augment to u,v,w vector
            nodeidx = find(force);
        end
    end
    
end

