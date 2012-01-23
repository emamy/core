classdef DefaultCompWiseKernelApprox < approx.algorithms.BaseKernelApproxAlgorithm
% Default component-wise kernel approximation algorithm
%
% @author Daniel Wirtz @date 2011-03-31
%
% @change{0,5,dw,2011-11-02} 
% - New interface for approximation computation: Passing an data.ApproxTrainData instance now
% instead of 'xi,ti,mui' parameters.
% - Re-enabled the clone method
%
% @new{0,5,dw,2011-07-07} Moved the old approx.DefaultCompWiseKernelApprox class to this
% class.
%
% @change{0,4,dw,2011-05-20} Removed the ApproxExpansionSize property as this is now determined
% by the approx.TrainDataSelector property
%
% @change{0,3,sa,2011-04-21} Implemented Setter for the class property
%
% @new{0,3,dw,2011-03-31} Added this class to keep old approximation
% generation method.
%
% See also: BaseApprox KernelApprox
       
    methods        
        function copy = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            copy = approx.algorithms.DefaultCompWiseKernelApprox;
            
            copy = clone@approx.algorithms.BaseKernelApproxAlgorithm(this, copy);
        end
    end

    methods(Access=protected) 
        function detailedComputeApproximation(this, kexp, atd)
            % Set AKernelCoreFun centers
            kexp.Centers.xi = atd.xi;
            kexp.Centers.ti = atd.ti;
            kexp.Centers.mui = atd.mui;
            
            % Call coeffcomp preparation method and pass kernel matrix
            this.CoeffComp.init(data.MemoryKernelMatrix(kexp.getKernelMatrix), kexp);
            
            % Call protected method
            this.computeCoeffs(kexp, atd.fxi, []);
        end          
    end
end


