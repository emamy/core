classdef TestFunctions
% TestFunctions: Some test functions for nonlinear approximation methods
%
% This class contains some static methods that create function handles to
% be used within test cases for different approximation algorithms in the
% +approx package.
%
% @author Daniel Wirtz @date 2012-04-30
%
% @new{0,6,dw,2012-04-30} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
        
    methods(Static)
%         function f = Nielson
%             f = @(x)x;
%         end
        
        function [f, r] = F3
            % Function "F3" from [RB99]
            f = @(x,~,~)(1.25 + cos(5.4*x(2,:))) ./ (6+6*(3*x(1,:)-1).^2);
            r = [0 1; 0 1];
        end
        
        function [f, r] = F7
            % Function "F10" from [RB99]
            f = @(x,~,~)2*cos(10*x(1,:)).*sin(10*x(2,:)) + sin(10*x(1,:).*x(2,:));
            r = [0 1; 0 1];
        end
        
        function [f, r] = F8
            % Function "F8" from [RB99]
            f = @(x,~,~)exp(-((5-10*x(1,:)).^2)/2) + .75*exp(-((5-10*x(2,:)).^2)/2)...
                +.75*exp(-((5-10*x(1,:)).^2)/2).*exp(-((5-10*x(2,:)).^2)/2);
            r = [0 1; 0 1];
        end
        
        function [f, r] = F10
            % Function "F10" from [RB99]
            f = @(x,~,~)exp(-0.04*sqrt(sum((80*x(1,:)-40).^2 + (90*x(2,:)-45).^2)))...
                *cos(.15*sqrt(sum((80*x(1,:)-40).^2 + (90*x(2,:)-45).^2)));
            r = [0 1; 0 1];
        end
        
        function [f, r] = Franke2D
            % Test function from [FN80]
            f = @(x,~,~).75*exp(-sum((9*x-2).^2)/4) ...
                + .75*exp(-((9*x(1,:)+1).^2)/49 - (9*x(2,:)+1)/10) ...
                + .5*exp(-((9*x(1,:)-7).^2 + (9*x(2,:)-3).^2)/4) ...
                - .2*exp(-(9*x(1,:)-4).^2 - (9*x(2,:)-7).^2);
            r = [0 1; 0 1];
        end
        
        function [f, r] = Franke3D
            % Test function from [LM05, FN80]
            f = @(x,~,~).75*exp(-sum((9*x-2).^2)/4) ...
                + .75*exp(-((9*x(1,:)+1).^2)/49 - sum((9*x(2:3,:)+1)/10)) ...
                + .5*exp(-((9*x(1,:)-7).^2 + (9*x(2,:)-3).^2 + (9*x(3,:)-5).^2)/4) ...
                - .2*exp(-(9*x(1,:)-4).^2 - (9*x(2,:)-7).^2 - (9*x(3,:)-5).^2);
            r = [0 1; 0 1; 0 1];
        end
    end
end