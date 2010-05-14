classdef PODReducer < spacereduction.BaseSpaceReducer
    %PODREDUCER Uses POD for reduced space generation.
    %
    % Internally the SVD decomposition of the snapshot array is used.
    % Several modes are supported to enable more specific reduced space
    % selection.
    %
    % @DanielWirtz, 19.03.2010
    
    properties
        % The modus used to generate the reduced space.
        %
        % Possible choices are:
        % 'sign': All eigenvectors with evalues larger than 'value' percent of the
        %         largest eigenvalue are used. Use 0 for all.
        % 'eps':  All eigenvectors with evalues larger than 'value' are used.
        % 'rel' : Reduction to 'value' percent of the original space dimension.
        % 'abs' : Explicitly specified reduced space dimension. Uses the first
        %         'value' eigenvectors as space.
        %
        % The fastest modes are 'rel' and 'abs' as for these cases the
        % target dimension can be computed before svd decomposition; for
        % 'sign' and 'eps' all singular values are computed and the
        % dimension is computed using these singular values.
        %
        % Defaults to 'rel'.
        %
        % See also: Value
        Mode = 'rel';
        
        % The value associated with the selected Mode
        %
        % Defaults to .3
        %
        % See also: Mode
        Value = .3;
    end
    
    methods
        function V = generateReducedSpace(this, model)
            % Implements the abstract method from BaseSpaceReducer
            
            %% Preparation
            % Collapse parameter samples and input dimension.
            % Transpose so that the columns are dimensions and rows measures
            data = model.Data.PlainSnapshotArray;   
            
            %% Dimension checks
            % This is for the "explicit" modes where the target dimension
            % is either a fixed value or a fraction of the full dimension.
            target_dim=[];
            if strcmpi(this.Mode,'rel')
                target_dim = ceil(size(data,1)*this.Value);
            elseif strcmpi(this.Mode,'abs')
                target_dim = this.Value;
            end
            if ~isempty(target_dim) && target_dim >= size(data,1)
                % Yell boo if no reduction achieved!
                error('Reduced space dimension (%d) has to be smaller than the full state space dimension (%d)!',target_dim,size(data,1));
            end
            
            %% Reduction for modes 'sign' and 'eps'
            % For these "dynamic" modes the full singular values have
            % to be computed in order to determine how many to use.
            if any(strcmpi(this.Mode,{'sign','eps'}))    
                [u,s,v] = svd(data,'econ');
                s = diag(s);
                if strcmpi(this.Mode,'sign')
                    sig = s >= s(1)*this.Value;
                elseif strcmpi(this.Mode,'eps')
                    sig = s >= this.Value;
                end
                % Select wanted subspace
                u = u(:,sig);
                v = v(:,sig);
                s = s(sig);
            else
                %% Reduction for modes 'abs' and 'rel'
                % For cases 'abs' or 'rel': fixed target dimension.
                % So just let svds extract the wanted components!
                [u,s,v] = svds(data, target_dim);
                s = diag(s);
            end
            
            % Case N >> d (more samples, "normal case")
            if size(data,2) > size(data,1)
                V = u;
            else
                % Case d >> N ("undersampled")
                V = data * v * diag(s.^-1);
            end
        end
        
        function set.Mode(this, value)
            if ~any(strcmpi(value, {'sign','eps','abs','rel'}))
                error(['Unknown POD reduction mode: ''' value '''\nAllowed: sign, eps, abs, rel']);
            end
            this.Mode = lower(value);
        end
    end
    
    methods(Static)
        
        function res = test_POD
            
            model = models.BaseFullModel;
            pod = spacereduction.PODReducer;
            
            res = true;
            
            % d << N
            model.Data.Snapshots = rand(25,25,10,3);
            res = res && spacereduction.PODReducer.internalPODTest(model,pod);
            
            % d >> N
            model.Data.Snapshots = rand(800,25,5,2);
            res = res && spacereduction.PODReducer.internalPODTest(model,pod);
        end
        
    end
    
    methods(Static,Access=private)
        function res = internalPODTest(model, pod)
            res = true;
            pod.Mode = 'eps';
            pod.Value = 7;
            fprintf('eps..');
            V = pod.generateReducedSpace(model);
            res = res && isequal(round(V'*V),eye(size(V,2)));
            
            pod.Mode = 'sign';
            pod.Value = .3;
            fprintf('sign..');
            V = pod.generateReducedSpace(model);
            res = res && isequal(round(V'*V),eye(size(V,2)));
            
            pod.Mode = 'rel';
            pod.Value = .3;
            fprintf('rel..');
            V = pod.generateReducedSpace(model);
            res = res && isequal(round(V'*V),eye(size(V,2)));
            
            pod.Mode = 'abs';
            pod.Value = 10;
            fprintf('abs..');
            V = pod.generateReducedSpace(model);
            res = res && isequal(round(V'*V),eye(size(V,2)));
        end
    end
    
end

