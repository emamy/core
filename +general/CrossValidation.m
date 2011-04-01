classdef CrossValidation
    %CROSSVALIDATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function crossvalidate(this, data)
            subdata = this.getNextTrainingSet(data);
            while (~isempty(subdata))
                res = this.train_using(subdata);
                this.judge(res);
                subdata = this.getNextTrainingSet(data);
            end
        end
        
    end
    
end

