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
        % - 'sign': All eigenvectors with evalues larger than 'value' percent of the
        %         largest eigenvalue are used. Use 0 for all.
        % - 'eps':  All eigenvectors with evalues larger than 'value' are used.
        % - 'rel' : Reduction to 'value' percent of the original space dimension.
        % - 'abs' : Explicitly specified reduced space dimension. Uses the first
        %         'value' eigenvectors as space.
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
        
        function podvec = computePOD(this, vec)
            % Computes the POD vectors according to the specified settings
            %
            % See also: general.POD#Mode general.POD#Value
            
            %% Dimension checks
            % This is for the "explicit" modes where the target dimension
            % is either a fixed value or a fraction of the full dimension.
            target_dim=[];
            if strcmpi(this.Mode,'rel')
                target_dim = ceil(min(size(vec))*this.Value);
            elseif strcmpi(this.Mode,'abs')
                target_dim = this.Value;
            end
            if ~isempty(target_dim) && target_dim > min(size(vec))
                % Yell boo if no reduction achieved!
                error('Targeted reduced space dimension (%d) has to be <= the smallest matrix dimension (%d)!',target_dim,min(size(vec)));
            end
            
            %% Reduction for modes 'sign' and 'eps'
            % For these "dynamic" modes the full singular values have
            % to be computed in order to determine how many to use.
            if any(strcmpi(this.Mode,{'sign','eps'}))
                [u,s,v] = svd(vec,'econ');
                s = diag(s);
                if strcmpi(this.Mode,'sign')
                    sig = s >= s(1)*this.Value;
                elseif strcmpi(this.Mode,'eps')
                    sig = s >= this.Value;
                end
                
                err = sum(s(~sig));
                rerr = 100*err/sum(s);
                fprintf('POD mode ''%s'' with value %2.6f: Selecting %d singular values, error %.4e (%.4e%% relative)\n',this.Mode,this.Value,length(find(sig)),err,rerr);
                % Select wanted subspace
                u = u(:,sig);
                v = v(:,sig);
                s = s(sig);
            else
                %% Reduction for modes 'abs' and 'rel'
                % For cases 'abs' or 'rel': fixed target dimension.
                % So just let svds extract the wanted components!
                
                % As tests showed that the svds method is less reliable
                % computing correct subspaces, an option is added that
                % explicitly allows to choose if svds is used rather than
                % svd.
                if this.UseSVDS
                    fprintf('POD mode ''%s'' with value %2.6f: Warning, SVDS results are non-reproducable (Arnoldi with random seed)\n',this.Mode,this.Value);
                    [u,s,v] = svds(vec, target_dim);
                    s = diag(s);
                else
                    [u,s,v] = svd(vec,'econ');
                    s = diag(s);
                    err = sum(s(target_dim+1:end));
                    rerr = 100*err/sum(s);
                    fprintf('POD mode ''%s'' with value %2.6f: Target dimension is %d singular values, error %.4e (%.4e%% relative)\n',this.Mode,this.Value,target_dim,err,rerr);
                    u(:,target_dim+1:end) = [];
                    v(:,target_dim+1:end) = [];
                    s(target_dim+1:end) = [];
                end
            end
            
            % Safety for zero singular values.
            z = find(s == 0);
            if ~isempty(z)
                warning('general:pod','%d of %d singular values are zero. Assuming output size %d',length(z),length(s),length(s)-length(z));
            end
            
            if size(vec,2) > size(vec,1) % Case N >> d (more columns than dimensions, "normal case")
                % Cut out columns of zero singular values
                u(:,z) = [];
                podvec = u;
            else % Case d >> N ("undersampled")
                % Cut out columns of zero singular values
                s(z) = [];
                v(:,z) = [];
                podvec = vec * v * diag(s.^-1);
            end
            
            if any(isinf(podvec(:)))
                error('Inf values occured in POD vectors.');
            end
            if any(isnan(podvec(:)))
                error('NaN values occured in POD vectors.');
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
            if ~isposrealscalar(value)
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
            vec = rand(50,100);
            res = res && general.POD.internalPODTest(pod, vec);
            
            % d >> N
            vec = rand(100,50);
            res = res && general.POD.internalPODTest(pod, vec);
        end
        
    end
    
    methods(Static,Access=private)
        function res = internalPODTest(pod, vec)
            res = true;
            pod.Mode = 'eps';
            pod.Value = 7;
            fprintf('eps..');
            V = pod.computePOD(vec);
            res = res && isequal(round(V'*V),eye(size(V,2)));
            
            pod.Mode = 'sign';
            pod.Value = .3;
            fprintf('sign..');
            V = pod.computePOD(vec);
            res = res && isequal(round(V'*V),eye(size(V,2)));
            
            pod.Mode = 'rel';
            pod.Value = .3;
            fprintf('rel..');
            V = pod.computePOD(vec);
            res = res && isequal(round(V'*V),eye(size(V,2)));
            
            pod.Mode = 'abs';
            pod.Value = 10;
            fprintf('abs..');
            V = pod.computePOD(vec);
            res = res && isequal(round(V'*V),eye(size(V,2)));
        end
    end
    
end

