classdef JKerMorExportable < handle
% JKerMorExportable: 
%
%
%
% @author Daniel Wirtz @date 2011-09-29
%
% @new{0,5,dw,2011-09-29} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The package of any functions associated with this model when
        % exported to JKerMor.
        %
        % See also: KerMor.JKerMorSource
        JavaExportPackage = '';
    end
    
    methods(Abstract)
        exportGeometry(this, f);
    end
    
end