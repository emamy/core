classdef POD < handle
    %POD Implements proper orthogonal decomposition
    %
    % See also: #Mode #Value
    %
    % @author Daniel Wirtz @date 24.08.2010
    
    properties
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
                % Select wanted subspace
                u = u(:,sig);
                v = v(:,sig);
                s = s(sig);
            else
                %% Reduction for modes 'abs' and 'rel'
                % For cases 'abs' or 'rel': fixed target dimension.
                % So just let svds extract the wanted components!
                [u,s,v] = svds(vec, target_dim);
                s = diag(s);
            end
            
            % Case N >> d (more columns than dimensions, "normal case")
            if size(vec,2) > size(vec,1)
                podvec = u;
            else
                % Case d >> N ("undersampled")
                podvec = vec * v * diag(s.^-1);
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

