classdef ReducedModel < models.BaseModel
    %REDUCEDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        
        function this = ReducedModel(fullmodel)
            % Copy common values from the full model
            this.T = fullmodel.T;
            this.dt = fullmodel.dt;
            this.Verbose = fullmodel.Verbose;
            this.ODESolver = fullmodel.ODESolver;
            % Update name ;-)
            this.Name = ['Reduced: ' fullmodel.Name];
        end
        
    end
    
end

