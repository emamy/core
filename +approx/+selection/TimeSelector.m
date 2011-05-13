classdef TimeSelector < approx.selection.ASelector
% TimeSelector: Approximation training data selection utilizing time information
%
% This algorithm searches for unique values of times in the training data and determines how many
% samples for each time should be taken according to it's apperance in the training data. Due to
% rounding off numbers usually this process yields less than MaxSize elements, so additionally
% linspaced elements are selected to fill up the MaxSize selection.
%
% @author Daniel Wirtz @date 2011-04-12
%
% @new{0,4,dw,2011-05-06} 
% - Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
% - Also implemented the ICloneable interface.
%
% @new{0,3,dw,2011-04-12} Added this class for a selection of training samples that utilizes time
% information.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
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
            this = this@approx.selection.ASelector;
            this.registerProps('Size','Seed');
        end
        
        function copy = clone(this)
            copy = approx.selection.TimeSelector;
            copy = clone@approx.selection.ASelector(this, copy);
            copy.Size = this.Size;
            copy.Seed = this.Seed;
        end
    end
    
    methods(Access=protected,Sealed)
        function atd = select(this, model)
            % Performs selection of samples adjusted to the apperances of different times.
            sn = model.Data.TrainingData;
            if (size(sn,2) > this.Size)
                times = unique(sn(3,:));
                tcnt = length(times);
                occ = cell(1,tcnt);
                num = zeros(1,tcnt);
                
                % Collect time apperances
                for idx=1:tcnt
                    occ{idx} = find(sn(3,:) == times(idx));
                    num(idx) = length(occ{idx});
                end
                num = round((this.Size/sum(num)) * num);
                
                % Get random numtimes many samples for each time
                selection = [];
                s = RandStream.create('mt19937ar','Seed',this.Seed);
                for idx=1:tcnt
                    sel = occ{idx};
                    selidx = s.randperm(length(sel));
                    selection(end+1:end+num(idx)) = sel(selidx(1:num(idx)));
                end
                
                % Fill leftover sample places with linspaced leftovers
                left = setdiff(1:size(sn,2),selection);
                addidx = round(linspace(1,length(left),this.Size-length(selection)));
                selection = sort([selection left(addidx)]);
                atd = sn(:,selection);
                this.LastUsed = selection;
            else
                atd = sn;
                this.LastUsed = 1:size(sn,2);
            end
        end
    end
end