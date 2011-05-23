function trycatch_speedtest(loopsize)
% trycatch_speedtest: 
%
% Demonstrate how slow try-catch blocks are. 
% Copied from http://www.mathworks.com/matlabcentral/newsreader/view_thread/275243.
%
% @author Daniel Wirtz @date 2011-05-18
%
% @new{0,4,dw,2011-05-18} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
% 

t = [];

if nargin < 1
    loopsize = 1000000;
end

tic
for ii = 1:loopsize
       A = ii;   
end
t(end+1) = toc;

tic
try
    for ii = 1:loopsize
           A = ii.^2;   
    end
catch ME
end
t(end+1) = toc;

tic
for ii = 1:loopsize
   try
       A = ii.^2;
   catch ME
   end
end
t(end+1) = toc;

tic
for ii = 1:loopsize
    if ~ii
        % Never here.
%     elseif 0
%         % Never here.
%     elseif mod(ii,2)
%         B = ii;
    else
        C = ii.^2;
    end
end
t(end+1) = toc;

t
reltimes = t./min(t) % Display the relative times of the two loops. 

end

