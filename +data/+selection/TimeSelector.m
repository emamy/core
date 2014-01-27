classdef TimeSelector < data.selection.ASelector
% TimeSelector: Approximation training data selection utilizing time information
%
% This algorithm searches for unique values of times in the training data and determines how many
% samples for each time should be taken according to it's apperance in the training data. Due to
% rounding off numbers usually this process yields less than MaxSize elements, so additionally
% linspaced elements are selected to fill up the MaxSize selection.
%
% @author Daniel Wirtz @date 2011-04-12
%
% @change{0,5,dw,2011-11-09} Re-enabled the use of this selector after adopting to new
% data.ApproxTrainData structure.
%
% @change{0,5,dw,2011-08-04} Disabled the use of this selector, as the new data.ATrajectoryData structure
% does not cater sensefully for this type of approximation training data selection.
%
% @new{0,4,dw,2011-05-06} 
% - Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
% - Also implemented the ICloneable interface.
%
% @change{0,3,sa,2011-04-20} Implemented Setters for the class properties.
%
% @new{0,3,dw,2011-04-12} Added this class for a selection of training samples that utilizes time
% information.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(SetObservable)
        % The (maximum) size of training samples to take
        %
        % @propclass{critical} Determines the (maximum) size of the approximation training data
        % samples to take.
        %
        % @default 10000
        Size = 10000;
        
        % The seed for the per-time-random selection (to enable reproduction of results)
        %
        % @propclass{optional} The seed for the random number generator.
        %
        % @default 1
        Seed = 1;
    end
    
    methods
        function this = TimeSelector
            this = this@data.selection.ASelector;
            this.registerProps('Size','Seed');
        end
        
        function copy = clone(this)
            copy = data.selection.TimeSelector;
            %copy = clone@data.selection.ASelector(this, copy);
            copy.Size = this.Size;
            copy.Seed = this.Seed;
        end
    end
    
    methods(Access=protected,Sealed)
        function [xi, ti, mui] = select(this, model)
            % Performs selection of samples adjusted to the apperances of different times.
            %
            % Parameters:
            % model: The full model with the training data @type models.BaseFullModel
            %
            % Return values:
            % xi: The selected `x_i = x(t_i)` training data @type matrix
            % ti: The selected training times `t_i` @type rowvec
            % mui: The selected parameter samples `\mu_i` with which the states
            % `x_i` have been reached @type matrix
            %
            % @todo Adopt to current KerMor model data structure.
            
            % Build sn array similar to old one
            d = model.Data.TrajectoryData;
            nt = d.getNumTrajectories;
            times = model.Times;
            tlen = length(times);
            
            if nt*tlen > this.Size
                % So many samples for each timestep (fill rest with random samples)
                num = floor(this.Size/tlen);
                s = RandStream.create('mt19937ar','Seed',this.Seed);
                % These are the time step indices that have to be selected overall
                selidx = repmat(s.randperm(tlen),1,num);
                % How many per trajectory?
                % MUST be smaller than tlen, as otherwise selecting more timesteps per
                % trajectory than available would have meant that nt*tlen < this.Size.
                % This means also that repeating the random permutation of 1:tlen will not lead
                % to double selection of a time within the same trajectory.
                protraj = floor(length(selidx)/nt); 
                
                sel = logical.empty;
                for tn = 1:nt
                    % Select nothing by default
                    hlp = false(1,tlen);
                    % select "protraj" many sequentially
                    step = (1:protraj)+(tn-1)*protraj;
                    % set those to true
                    hlp(selidx(step)) = true;
                    % Store in selection array at range for current trajectory
                    sel = [sel hlp]; %#ok<*AGROW>
                end
                
                % Add up to full this.Size selection
                left = this.Size - protraj * nt;
                notsel = find(~sel);
                selsel = s.randperm(length(notsel),left);
                sel(notsel(selsel)) = true;
            else
                sel = true(1,nt*tlen); % Select everything!
            end

            xi = []; ti = []; mui = [];
            for tn = 1:nt
                [x, mu] = d.getTrajectoryNr(tn);
                idx = (1:tlen) + (tn-1)*tlen;
                xi = [xi x(:,sel(idx))];
                newt = times(sel(idx));
                ti = [ti newt];
                if ~isempty(mu)
                    mui = [mui mu(:,ones(1,length(newt)))];
                end
            end
        end
        
%         function [xi, ti, mui] = select(this, model)
%             % Performs selection of samples adjusted to the apperances of different times.
%             %
%             % Parameters:
%             % model: The full model with the training data @type models.BaseFullModel
%             %
%             % Return values:
%             % xi: The selected `x_i = x(t_i)` training data @type matrix
%             % ti: The selected training times `t_i` @type rowvec
%             % mui: The selected parameter samples `\mu_i` with which the states
%             % `x_i` have been reached @type matrix
%             %
%             % @todo Adopt to current KerMor model data structure.
%             
%             % Build sn array similar to old one
%             d = model.Data;
%             tn = d.getNumTrajectories;
%             tlen = length(model.Times);
%             [x, mu] = d.getTrajectoryNr(1);
%             xdim = size(x,1);
%             sn = zeros(xdim+1+size(mu,2),tlen*tn);
%             for i=1:tn
%                 [x, mu] = d.getTrajectoryNr(i);
%                 sn(:,1+(i-1)*tlen:i*tlen) = [model.Times; x; mu(:,ones(1,tlen))];
%             end
%             
%             if (size(sn,2) > this.Size)
%                 times = unique(sn(1,:));
%                 tcnt = length(times);
%                 occ = cell(1,tcnt);
%                 num = zeros(1,tcnt);
%                 
%                 % Collect time apperances
%                 for idx=1:tcnt
%                     occ{idx} = find(sn(1,:) == times(idx));
%                     num(idx) = length(occ{idx});
%                 end
%                 num = round((this.Size/sum(num)) * num);
%                 
%                 % Get random numtimes many samples for each time
%                 selection = [];
%                 s = RandStream.create('mt19937ar','Seed',this.Seed);
%                 for idx=1:tcnt
%                     sel = occ{idx};
%                     selidx = s.randperm(length(sel));
%                     selection(end+1:end+num(idx)) = sel(selidx(1:num(idx)));
%                 end
%                 
%                 % Fill leftover sample places with linspaced leftovers
%                 left = setdiff(1:size(sn,2),selection);
%                 addidx = round(linspace(1,length(left),this.Size-length(selection)));
%                 selection = sort([selection left(addidx)]);
%                 atd = sn(:,selection);
%             else
%                 atd = sn;
%             end
%             
%             % Transform atd back
%             ti = atd(1,:);
%             xi = atd(2:xdim+1,:);
%             mui = atd(xdim+2:end,:);
%         end
    end
    
    methods           
        function set.Size(this, value)
            if ~isposintscalar(value)
                error('The value must be a finite positive integer.');
            end
            this.Size = value;
        end
        
        function set.Seed(this, value)
            if ~isreal(value)
                error('The value must be a real value.');
            end
            this.Seed = value;
        end
    end
end