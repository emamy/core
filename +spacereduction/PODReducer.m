classdef PODReducer < spacereduction.BaseSpaceReducer
    %PODREDUCER Uses POD for reduced space generation.
    %
    % @DanielWirtz, 11.03.2010
    
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
        % See also: Value
        Mode = 'rel';
        
        % The value associated with the selected Mode
        %
        % See also: Mode
        Value = .3;
    end
    
    methods
        function V = generateReducedSpace(this, model)
            % Implements the abstract method from BaseSpaceReducer
            
            %% Preparation
            % Collapse parameter samples and input dimension
            sn = model.Data.PlainSnapshotArray;
            
            A = sn'*sn;
            
            useeig = false;
            %% Dimension checks
            if strcmpi(this.Mode,'rel')
                target_dim = round(size(sn,1)*this.Value);
                if target_dim > size(A,2)
                    warning('KerMor:PODReducer:mode_rel','Reduced space cant be bigger (%d) than sample size (%d). Using sample size.',target_dim,size(A,2));
                    useeig = true;
                end
            elseif strcmpi(this.Mode,'abs')
                target_dim = this.Value;
                if target_dim > size(A,2)
                    warning('KerMor:PODReducer:mode_abs','''Value'' cant be bigger (%d) than sample size (%d). Setting to sample size.',this.Value,size(A,2));
                    this.Value = size(A,2);
                    useeig = true;
                end
            end
            %% Actual computation
            if any(strcmpi(this.Mode,{'sign','eps'})) || useeig
                % compute full eigenvalues for sign/eps mode or if other modes result
                % in full space, too. (efficiency)
                [ev,ew] = eig(A);
                % Invert
                ev = flipdim(ev,2);
                ew = flipdim(diag(ew),1);
            else
                % Forward verbose setting to eigs fcn
                opts.disp = model.Verbose;
                [ev, ew] = eigs(A,target_dim,'lm',opts);
                ew = diag(ew);
            end
            %% Reduction for modes 'sign' and 'eps'
            if strcmpi(this.Mode,'sign')
                sig = ew >= ew(1)*this.Value;
                % Reduce
                ev = ev(:,sig);
                ew = ew(sig);
            elseif strcmpi(this.Mode,'eps')
                sig = ew >= this.Value;
                % Reduce
                ev = ev(:,sig);
                ew = ew(sig);
            end
            
            % Compute projection matrix
            V = sn * ev * diag(ew.^(-0.5));
        end
        
        function set.Mode(this, value)
            if ~any(strcmpi(value, {'sign','eps','abs','rel'}))
                error(['Unknown POD reduction mode: ''' value '''\nAllowed: sign, eps, abs, rel']);
            end
            this.Mode = lower(value);
        end
    end
    
    methods(Static)
        
        function test_POD
            
            model = models.BaseFullModel;
            model.Data.Snapshots = rand(25,25,10,3);
            
            pod = spacereduction.PODReducer;
            
            pod.Mode = 'eps';
            pod.Value = 1e-2;
            fprintf('eps..');
            pod.generateReducedSpace(model);
            
            pod.Mode = 'sign';
            pod.Value = .3;
            fprintf('sign..');
            pod.generateReducedSpace(model);
            
            pod.Mode = 'rel';
            pod.Value = .3;
            fprintf('rel..');
            pod.generateReducedSpace(model);
            
            pod.Mode = 'abs';
            pod.Value = 10;
            fprintf('abs..');
            pod.generateReducedSpace(model);
        end
    end
    
end

