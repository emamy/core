classdef MUnit
    %MUNIT Class Unit Testing Framework for Matlab
    %
    % This class allows to run tests within a OO-Based Matlab program
    % according to the conventions mentioned below:
    %
    % Any test method has to satisfy:
    % - The method's name begins with "TestFunctionPrefix"
    % - Static
    % - No input arguments
    % - One (optional) output argument indicating success or failure of the
    %   test
    % Possible exceptions do not have to be taken care of as they are
    % caught automatically by the system, resulting in a failure of the
    % test.
    %
    % See also: TestFunctionPrefix
    %
    % @author Daniel Wirtz @date 12.03.2010
    
    properties(Constant)
        % The prefix for any function that will be detected in the MUnit
        % testing process.
        %
        % Defaults to 'test_'
        TestFunctionPrefix = 'test_';
        
        % Green color (the default is hardly readable)
        GreenCol = [0 .5 0];
    end
    
    methods(Static)
        function RunClassTests(dir)
            % Runs class tests recursively within the specified path.
            %
            % If no path is given, the current directory (cd) is used.
            %
            % Within the path any static methods from classes beginning
            % with the 'TestFunctionPrefix' are detected and executed.
            %
            % After running all tests a summary is provided.
            
            % Check for no specified path
            if nargin == 0
                dir = cd;
            end
            
            % Check if the current directory path contains packages and
            % extract them if necessary
            curPackage = '';
            pckidx = strfind(dir,'+');
            if ~isempty(pckidx)
                curPackage = dir(pckidx(1)+1:end);
                curPackage = [strrep(strrep(curPackage,'+','.'),'/','') '.'];
            end
            
            % Start recursive run
            [s,f] = testing.MUnit.recursiveRun(dir,curPackage);
            
            % Summary
            fprintf('\n\n All Class Tests finished.\n');
            cprintf(testing.MUnit.GreenCol,'Successful:%d\n',s);
            cprintf('Red','Failed:%d\n',f);
        end
    end
    
    methods(Static, Access=private)
        function [s,f] = recursiveRun(dir,currentPackage)
            w = what(dir);
            
            % Some weird afs stuff goin on!
            w=w(1);
            
            s = 0; f = 0;
            % Descend into any package
            for idx = 1:length(w.packages)
                subdir = fullfile(w.path,['+' w.packages{idx}]);
                %disp(['Descending into ' subdir]);
                [sa,fa] = testing.MUnit.recursiveRun(subdir, [currentPackage w.packages{idx} '.']);
                s = s + sa;
                f = f + fa;
            end
            
            pref = testing.MUnit.TestFunctionPrefix;
            pl = length(pref);
            
            % Run tests in current directory
            for idx = 1:length(w.m)
                [p,n] = fileparts(w.m{idx});
                
                %disp(['Checking ' n '...']);
                
                mc = meta.class.fromName([currentPackage n]);
                % mc is empty if m-file wasnt a class
                if ~isempty(mc)
                    
                    %disp(['Checking ' mc.Name '...']);
                    
                    for midx=1:length(mc.Methods)
                        m = mc.Methods{midx};
                        
                        %disp(['Method: ' mc.Name '.' m.Name]);
                        
                        % check for test function and if the method is
                        % actually declared within the current class
                        % (subclassing would otherwise lead to a repeated
                        % call to the test)
                        if length(m.Name) >= pl+1 ...
                                && strcmpi(m.Name(1:pl),pref) ...
                                && strcmp(m.DefiningClass.Name,mc.Name)
                            % check if the method is static [non-statics
                            % cant be run without instances.. :-)]
                            if m.Static
                                cprintf(testing.MUnit.GreenCol,['Running ' mc.Name ' -> ' m.Name(6:end) '... ']);
                                try
                                    eval(['outargs = nargout(@' mc.Name '.' m.Name ');']);
                                    if outargs > 0
                                        command = ['succ = ' mc.Name '.' m.Name ';'];
                                    else
                                        command = [mc.Name '.' m.Name ';'];
                                    end
                                    eval(command);
                                    if outargs == 0 || succ
                                        cprintf(testing.MUnit.GreenCol,'Success!\n');
                                        s = s+1;
                                    elseif ~succ
                                        cprintf('Red','Failure!\n');
                                        f = f+1;
                                    end
                                catch ME
                                    f = f+1;
                                    cprintf('Red','Failure due to Exception!\n');
                                    disp(getReport(ME));
                                    %rethrow(ME);
                                end
                            else
                                warning('MUnit:NonstaticTest',...
                                    ['Found non-static test function "' m.Name '" in ' mc.Name '\n'...
                                    'Maybe forgot to declare static?']);
                            end
                        end
                    end
                end
            end
        end
    end
end

