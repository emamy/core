classdef SteinwartSVM < dscomponents.ACoreFun
% SteinwartSVM: 
%
%
%
% @author Daniel Wirtz @date 2011-08-09
%
% @new{0,5,dw,2011-08-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        FileFolder = '/home/dwirtz/datastore/kermor/steinwart';
        
        FilePattern = 'bandscheibe';
    end
    
    methods
        function this = SteinwartSVM
            this = this@dscomponents.ACoreFun;
            this.MultiArgumentEvaluations = true;
            this.CustomProjection = false;
        end
        
        function y = evaluateCoreFun(this, x, t, mu)%#ok
            %% Init
            % Get output format
            fstr = [repmat('%f ',1,size(x,1)) '\n'];
            testfile = fullfile(this.FileFolder,'test.nla');
            outfile = fullfile(this.FileFolder,'res.psv');
            
            %% Open file and write values
            f = fopen(testfile,'w+'); 
            fprintf(f,fstr,x);
            fclose(f);
            
            %% Call svm-test
            VERBOSE=1;
%             TYPE='6 7 10';
%             BUILD=$3
            THREADS=4;
            
            y = zeros(size(x));
            % Compute for each dim
            for dim=1:size(x,1)
                bn=fullfile(this.FileFolder,sprintf('%s%d',this.FilePattern,dim));            
            
                cmd = sprintf('~/aghhome/Projects/IngoSVM/bin/svm-test -T %d -V 1 -L 1 -v %d %s.train.psv %s.sol %s none %s',...
                    THREADS,VERBOSE,bn,bn,testfile,outfile);
                system(cmd);
                
                f = fopen(outfile,'r');
                y(dim,:) = fscanf(f, '%f', [1 inf]);
                fclose(f);    
            end
            %delete(outfile);
            %delete(testfile);
        end
    end
    
end