classdef POD < KerMorObject
    % POD: Implements proper orthogonal decomposition
    %
    % Mainly relies on the matlab builtin functions svd and svds. This class wraps their
    % functionality to a class integrable with KerMor.
    %
    % @author Daniel Wirtz @date 2010-08-24
    %
    % See also: Mode Value
    %
    % @change{0,6,dw,2012-07-13} POD now also works with data.FileMatrix arguments (only modes
    % 'abs' and 'rel')
    %
    % @change{0,4,sa,2011-05-06} Implemented Setter for UseSVDS flag
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @todo - Create fixed random number stream for reproducable results!
    % - make Value property dependent and check in setter depending on Mode value!
    
    properties(SetObservable)
        % The modus used to generate the reduced space.
        %
        % Possible choices are:
        % - 'sign': All eigenvectors with singular values larger than
        %           'value' percent of the largest eigenvalue are used. Use
        %           0 for all.
        % - 'eps':  All eigenvectors with singular values larger than
        %           'value' are used.
        % - 'rel' : Reduction to 'value' percent of the original space
        %           dimension.
        % - 'abs' : Explicitly specified reduced space dimension. Uses the
        %           first 'value' eigenvectors as space.
        %
        % The fastest modes are 'rel' and 'abs' as for these cases the
        % target dimension can be computed before svd decomposition; for
        % 'sign' and 'eps' all singular values are computed and the
        % dimension is computed using these singular values.
        %
        % @propclass{important} Choices 'sign' and 'eps' can lead to long computation times if the
        % input samples are large as the full decomposition must be computed.
        %
        % @default 'rel'
        %
        % See also: Value
        Mode = 'rel';
        
        % The value associated with the selected Mode
        %
        % @propclass{important} Determines the size of the reduced space. Has different effects
        % depending on the choice of the general.POD.Mode property.
        %
        % @default .3
        %
        % See also: Mode
        Value = .3;
        
        % Flag whether to use the svds routine rather than svd.
        % Applies only for the Modes 'abs' and 'rel'.
        % For both 'sign' and 'eps' modes all singular values have to be
        % computed anyways and svd is faster if all values are needed.
        %
        % IMPORTANT: The results when using svds are NOT REPRODUCABLE as
        % svds uses an randomly initialized Arnoldi algorithm.
        %
        % @propclass{critical} Allows for dramatic reduction of computational costs for the price of
        % less accuracy. See the matlab docs for svds and svd.
        %
        % @default false
        UseSVDS = false;
    end
    
    methods
        
        function this = POD
            this = this@KerMorObject;
            
            this.registerProps('Mode','Value','UseSVDS');
        end
        
        function [podvec, s] = computePOD(this, data, Vexclude)
            % Computes the POD vectors according to the specified settings
            %
            % Depending on the type of data, any given Vinclude will be either projected out of
            % the underlying data (matrix case) or appended to the computed space of target_dim
            % minus Vinclude columns.
            %
            % Parameters:
            % data: The set of column-vectors to perform the POD on. Can be a real matrix or a
            % data.ABlockedData instance. @type matrix<double> or general.ABlockSVD
            % Vexclude: A given set of orthonormal vectors to exclude from the resulting POD.
            % @type matrix<double> @default []
            %
            % Return values:
            % podvec: The POD modes computed. @type matrix<double>
            % s: The singular values of the POD modes. If svds or FileMatrix is used,
            % these are only as many as POD modes, otherwise this vector
            % contains ALL singular values. @type rowvec<double>
            %
            % See also: general.POD#Mode general.POD#Value
            
            %% Dimension checks
            % This is for the "explicit" modes where the target dimension
            % is either a fixed value or a fraction of the full dimension.
            if nargin < 3
                Vexclude = [];
            end
            target_dim=[];
            if strcmpi(this.Mode,'rel')
                target_dim = ceil(min(size(data))*this.Value);
            elseif strcmpi(this.Mode,'abs')
                target_dim = this.Value;
            end
            if ~isempty(target_dim) && target_dim > min(size(data))
                warning('KerMor:general:POD','Targeted reduced space dimension (%d) has to be <= the smallest matrix dimension (%d). Using size %d.',...
                    target_dim,min(size(data)),min(size(data)));
                target_dim = min(size(data));
            end
            
            %% Block SVD case
            if isa(data,'data.ABlockedData')
                if ~any(strcmpi(this.Mode,{'abs','rel'}))
                    error('Cannot use POD on data.ABlockedData with setting "%s" (''abs'',''rel'' are allowed)',this.Mode);
                end
                if KerMor.App.Verbose > 2
                    [n,m] = size(data);
                    fprintf('Starting POD on ABlockedData(%s, %dx%d, %d blocks) with mode ''%s'' and value %g\n',...
                        class(data),n,m,data.getNumBlocks,this.Mode,this.Value);
                end
                [podvec, s] = data.getSVD(target_dim, Vexclude);
                s = diag(s);
            %% Matrix argument case
            else
                % Subtract the space spanned by Vexclude if given
                if ~isempty(Vexclude)
                    data = data - Vexclude*(Vexclude'*data);
                    target_dim = min(target_dim, min(size(data))-size(Vexclude,2));
                end
                
                % Create full matrix out of sparse fxi sets to enable use of
                % svd instead of svds
                nonzero = [];
                if issparse(data)
                    nonzero = find(sum(abs(data),2) ~= 0);
                    % Check if same sparsity pattern holds for each fxi column
                    if ~isempty(nonzero)
                        if KerMor.App.Verbose > 2
                            fprintf('Reducing data dimension from %d to %d as sparse matrix is given.',...
                                size(data,1),length(find(nonzero)));
                        end
                        fulldim = size(data,1);
                        data = full(data(nonzero,:));
                    end
                end
                
                if KerMor.App.Verbose > 2
                    fprintf('Starting POD on %dx%d matrix with mode ''%s'' and value %g (UseSVDS=%d) ... ',...
                        size(data),this.Mode,this.Value,this.UseSVDS);
                end
                %% Reduction for modes 'sign' and 'eps'
                % For these "dynamic" modes the full singular values have
                % to be computed in order to determine how many to use.
                if any(strcmpi(this.Mode,{'sign','eps'}))
                    [u,s] = svd(data,'econ');
                    s = diag(s);
                    if strcmpi(this.Mode,'sign')
                        sig = s >= s(1)*this.Value;
                    elseif strcmpi(this.Mode,'eps')
                        sig = s >= this.Value;
                    end
                    if KerMor.App.Verbose > 2
                        fprintf('selecting %d singular values ... ',length(find(sig)));
                    end
                    % Select wanted subspace
                    u = u(:,sig);

                    err = sum(s(~sig));
                    rerr = 100*err/sum(s);
                else
                    %% Reduction for modes 'abs' and 'rel'
                    % For cases 'abs' or 'rel': fixed target dimension.
                    % So just let svds extract the wanted components!

                    % As tests showed that the svds method is less reliable
                    % computing correct subspaces, an option is added that
                    % explicitly allows to choose if svds is used rather than
                    % svd.
                    if this.UseSVDS || issparse(data)
                        [u, s] = svds(data, target_dim);
                        s = diag(s);
                    else
                        [u,s] = svd(data,'econ');
                        s = diag(s);
                        err = sum(s(target_dim+1:end));
                        rerr = 100*err/sum(s);
                        u(:,target_dim+1:end) = [];
                    end
                end
                if KerMor.App.Verbose > 2
                    fprintf('error %g (%g%% relative)\n',err,rerr);
                end

                % Safety for zero singular values.
                z = find(s(1:size(u,2)) == 0);
                if ~isempty(z)
                    warning('KerMor:POD','%d of %d selected singular values are zero. Assuming output size %d',...
                        length(z),length(s),length(s)-length(z));
                end
                u(:,z) = [];

                % Re-create the sparse matrix structure after POD computation if input argument
                % was a sparse matrix
                if ~isempty(nonzero)
                    u = general.MatUtils.toSparse(u, nonzero, fulldim);
                end
                
                podvec = u;
                % Security checks (only for true matrices so far)
                if any(any(round(podvec'*podvec) ~= eye(size(podvec,2))))
                    warning('KerMor:POD',['Resulting POD vectors not ' ...
                        'sufficiently orthogonal. Orthonormalizing.']);
                    % Sometimes the u vectors generated by svds are not
                    % properly orthonormal
                    o = general.Orthonormalizer;
                    podvec = o.orthonormalize(podvec);
                end
            end
        end
        
        function set.UseSVDS(this, value)
            if ~islogical(value)
                error('Value should be logical');
            end
            this.UseSVDS = value;
        end
        
        function set.Mode(this, value)
            if ~any(strcmpi(value, {'sign','eps','abs','rel'}))
                error(['Unknown POD reduction mode: ''' value '''\nAllowed: sign, eps, abs, rel']);
            end
            this.Mode = lower(value);
        end
        
        function set.Value(this, value)
            if ~isreal(value) || ~isscalar(value)
                error('Value property must be a positive real scalar.');
            end
            this.Value = value;
        end
    end
    
    methods(Static)
        
        function res = test_POD
            pod = general.POD;
            res = true;
            
            % d << N
            vec = rand(50,1000);
            res = res && general.POD.internalPODTest(pod, vec);
            
            % d >> N
            vec = rand(1000,50);
            res = res && general.POD.internalPODTest(pod, vec);
        end
        
    end
    
    methods(Static,Access=private)
        function res = internalPODTest(pod, vec)
            res = true;
            pod.Mode = 'eps';
            pod.Value = 7;
            V = pod.computePOD(vec);
            res = res && isequal(round(V'*V),eye(size(V,2)));
            
            pod.Mode = 'sign';
            pod.Value = .3;
            V = pod.computePOD(vec);
            res = res && isequal(round(V'*V),eye(size(V,2)));
            
            pod.Mode = 'abs';
            pod.Value = 10;
            V = pod.computePOD(vec);
            res = res && isequal(round(V'*V),eye(size(V,2)));
            
            % Vexclude mode
            o = general.Orthonormalizer;
            exclu = o.orthonormalize(rand(size(vec,1),round(pod.Value/2)));
            Vex = pod.computePOD(vec,exclu);
            res = res && isequal(round(Vex'*Vex),eye(size(Vex,2)));
            res = res && norm(exclu'*Vex) < 1e-12;
            
            pod.Mode = 'abs';
            pod.Value = 10;
            Vf = pod.computePOD(data.FileMatrix(vec));
            res = res && isequal(round(Vf'*Vf),eye(size(Vf,2))) && norm(abs(V)-abs(Vf),'fro') < 1e-10;
            
            pod.Mode = 'rel';
            pod.Value = .3;
            V = pod.computePOD(vec);
            res = res && isequal(round(V'*V),eye(size(V,2)));
            
            pod.Mode = 'rel';
            pod.Value = .3;
            Vf = pod.computePOD(data.FileMatrix(vec));
            res = res && isequal(round(Vf'*Vf),eye(size(Vf,2))) && norm(abs(V)-abs(Vf),'fro') < 1e-10;
        end
    end
    
end

