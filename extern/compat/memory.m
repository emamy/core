function varargout = memory
% memory: Imitates the memory function for all platforms.
%
% Calls the built-in memory function on Windows systems.
%
% See also: memory
%
% @author Daniel Wirtz @date 2012-07-05
if ispc
    [varargout{1:nargout}] = builtin('memory');
else
    userview = struct;
    if isunix
        v = -1;
        cmd = 'cat /proc/meminfo | grep MemFree';
        [s, msg] = system(cmd);
        if s ~= 0
            fprintf(2,'Failed to run system command: %s: %s\n', cmd, msg);
        else
            v = str2double(strtrim(msg(9:end-3)))*1024;
        end
        userview.MaxPossibleArrayBytes = v;
        userview.MemAvailableAllArrays = v;
        
        v = -1;
        cmd = 'ps -e | grep MATLAB';
        [s, msg] = system(cmd);
        msg = strtrim(msg);
        if s ~= 0
            fprintf(2,'Failed to run system command "%s": %s\n', cmd, msg);
        else
            pos = strfind(msg,' ');
            pid = str2double(msg(1:pos(1)-1));
            %pid2 = ppid;
            cmd = sprintf('cat /proc/%d/status | grep VmSize',pid);
            [s, msg] = system(cmd);
            if s ~= 0
                fprintf(2,'Failed to run system command "%s": %s\n', cmd, msg);
            else
                v = str2double(strtrim(msg(8:end-3)))*1024;
            end
        end
        userview.MemUsedMATLAB = v;
        
        % Additional fields
        cmd = 'cat /proc/meminfo | grep MemTotal';
        [s, msg] = system(cmd);
        if s ~= 0
            fprintf(2,'Failed to run system command "%s": %s\n', cmd, msg);
        else
            userview.TotalSystemMemory = str2double(strtrim(msg(10:end-3)))*1024;
        end
    elseif ismac
        error('not implemented for mac yet.');
    end
    if nargout > 0
        varargout{1} = userview;
    else
        disp(userview);
    end
end       
end