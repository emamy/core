classdef SignalAnalyzer < handle
    %SignalAnalyzer Analyze signals
    % Thus far allows to extract the locations of peaks within a discrete
    % signal
    
    properties
        % The minimum signal threshold for a found peak to be considered
        minV = -20;
    end
    
    methods
        
        function peaks = getPeakIdx(this, t, signal)
            peaks = [];
            % Find all locations with negative derivative
            negpos = find(diff(signal(1,:))<0);
            if ~isempty(negpos)
                % Find the positions where the negative derivative starts
                firstpos = [1 find(diff(negpos) > 1)+1];
                % Check that the found locations are peaks and not in the lower noise
                ispeak = signal(1,negpos(firstpos)) > this.minV;
                % Remove unwanted
                firstpos(~ispeak) = [];
                peaks = negpos(firstpos);
                % Check for close consecutive peaks & remove
                tdiff = diff(t(peaks));
                closepeaks = tdiff < median(tdiff)/10;
                peaks(closepeaks) = [];
            end
        end
        
        function [peaks, cov_isi, cov_dr, freq] = analyze(this, t, signal)
            % Analyzes the signal and returns the peak locations,
            % coefficient of variants for the interspike interval/<dr> and
            % the frequency.
            %
            % The time unit equals the time unit of the provided time
            % vector.
            peaks = this.getPeakIdx(t, signal);
            pt = diff(t(peaks));
            cov_isi = std(pt)/mean(pt);
            dr = 1./(pt/1000);
            cov_dr = std(dr)/mean(dr);
            freq = mean(dr);
        end
    end
    
end

