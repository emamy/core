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
% @change{0,3,sa,2011-04-20} Implemented Setters for the class properties
%
% @new{0,3,dw,2011-04-12} Added this class for a selection of training samples that utilizes time
% information.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The (maximum) size of training samples to take
        %
        % @default 10000
        MaxSize = 10000;
        
        % The seed for the per-time-random selection (to enable reproduction of results)
        %
        % @default 1
        Seed = 1;
    end
    
    methods(Access=protected,Sealed)
        function atd = select(this, model)
            % Performs selection of samples adjusted to the apperances of different times.
            sn = model.Data.TrainingData;
            if (size(sn,2) > this.MaxSize)
                times = unique(sn(3,:));
                tcnt = length(times);
                occ = cell(1,tcnt);
                num = zeros(1,tcnt);
                
                % Collect time apperances
                for idx=1:tcnt
                    occ{idx} = find(sn(3,:) == times(idx));
                    num(idx) = length(occ{idx});
                end
                num = round((this.MaxSize/sum(num)) * num);
                
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
                addidx = round(linspace(1,length(left),this.MaxSize-length(selection)));
                selection = sort([selection left(addidx)]);
                atd = sn(:,selection);
                this.LastUsed = selection;
            else
                atd = sn;
                this.LastUsed = 1:size(sn,2);
            end
        end
    end
    
    methods           
        function set.MaxSize(this, value)
            if ~isposintscalar(value)
                error('The value must be a finite positive integer.');
            end
            this.MaxSize = value;
        end
        
        function set.Seed(this, value)
            if ~isposintscalar(value)
                error('The value must be a finite positive integer.');
            end
            this.Seed = value;
        end
    end
end