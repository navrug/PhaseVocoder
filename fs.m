% fs.m - synthesis from frequency spectra.
%     [output: waveform] = fs(freqresp: series of spectra)
%
% James Salsman, 1999-2000
%
% Note: sometimes this produces data outside the (1,-1) range, even when 
% run with input from that smaller domain given to af.m.  It is safe to 
% clamp outlying values to the extrema of the corresponding input domain.
% If the resulting signal sounds clipped, it is NOT because of clamping 
% domain outliers of the range (e.g., 1.2, -1.3) to the domain edges (1,-1). 

% To improve the quality, try increasing the length of the fft; if that
% doesn't work, then try increasing the sampling rate.  To improve the speed, 
% rewrite it in C; fftPrep is the only complex vector, and the call to 'sort' 
% is only for the ordering of the dot-product magnitudes, not the sorting.

function [output] = fs(freqresp)

samplingRate = 16000;                  % changed from 8000
frameRate    = 100;                    % should be even, maybe try 160
windowStep   = samplingRate/frameRate; % 16000/100 = 160, shall be even
% windowSize   = samplingRate/40;        % must be even; use 256 for 11k/s rate 
windowSize   = windowStep*2;           % this is better [2002]

[rows, cols] = size(freqresp);         % input dimensions
origFFtSize = rows*2;                  % should be power of two

% Allocate (preextend) all space needed for output array and helpers.
%
output     =  zeros(cols*windowStep + (windowSize - windowStep), 1);
FFtPrep    =  zeros(1, origFFtSize) + 0*i;  % complex ifft input
iFFtReal   =  zeros(1, origFFtSize);        % real ifft output
phases     =  zeros(1, origFFtSize/2);      % resynthesis phase array
newPhases  =  zeros(1, origFFtSize/2);      % uses "-1" for unassigned phases
transInd   =  zeros(1, origFFtSize/2);      % index into bins by magnitude
ignore     =  zeros(1, origFFtSize/2);      % sort values unused but can't hurt

for n = 1:origFFtSize/2,              % Initialize synthesis phases.
    phases(n) = pi*mod(n, 2);         % Each starting phase 180 deg's apart
end;                                  %   from the next one.

% pre-calculate phase transition from one frame to the next for base 2nd bin
%
phaseFactor = 2*pi*samplingRate/(frameRate*origFFtSize);

% resynthesis COLA window should be Portnoff (sinc) because iFFT will 
% shape like the window it was taken from
%
wind = sinc((-windowSize/2:windowSize/2-1)/(windowSize/2));  % pure sinc

for c = 1:cols,                       % iterate over columns
 
    first = (c-1)*windowStep + 1;     % output index
    last = first + windowSize - 1;    %   boundaries
 
    for b = 2:origFFtSize/2,          % iterate over bins 

        % prepare resynthesis
        %
        [re, im] = pol2cart(phases(b), freqresp(b, c)); % polar to cartesian
        FFtPrep(b) = re + i*im;       % rectangular to complex
        FFtPrep(origFFtSize - b+2) = conj(FFtPrep(b)); % the conjugates count

    end;
 
    iFFtReal = real(ifft(FFtPrep));   % resynthesize
 
    output(first:last) = output(first:last) ... % note: column vector (')
      + (iFFtReal( origFFtSize/2 - windowSize/2 : ... 
                   origFFtSize/2 + windowSize/2 - 1  ) .* wind )';  % COLA

    if c < cols,                      % figure out the next set of phases
        
        newPhases = ones(1, origFFtSize/2)*-1;  % "-1" for unassigned phase

        % set transInd to index the least-to-greatest transitioning magnitudes
        %
        [ignore, transInd] = sort(freqresp(:,c)' .* freqresp(:,c+1)');

        % iterate downwards such that indexed bin will be greatest mags first
        %
        for idx = origFFtSize/2:-1:1, % step by -1 from greatest to least
        
            b = transInd(idx);        % get bin from index

            if b-1 > 0,               % adjacent lower neighboring bin
                lower = newPhases(b-1); 
            else
                lower = -1;           % nonexistent = unassigned
            end;

            if b+1 < origFFtSize/2,   % adjacent higher neighboring bin
                higher = newPhases(b+1); 
            else
                higher = -1;          % nonexistent = unassigned
            end;

            if (lower == -1) & (higher == -1),    % both unassigned: propagate
                newPhases(b) = mod(phases(b) + (b-1)*phaseFactor, 2*pi);
            elseif lower == -1,                   % make 180 degs from higher
                newPhases(b) = mod(higher + pi, 2*pi);
            elseif higher == -1,                  % make 180 degs from lower
                newPhases(b) = mod(lower + pi, 2*pi);
            else                                  % both neighbors assigned
                avgAng = (lower + higher)/2;      % Mean will be 180 deg off
                if abs(lower - higher) > pi,      %   when the angles are that 
                   avgAng = avgAng + pi;          %   far apart; if so, adjust
                end;                              %   180 degs from mean angle.
                newPhases(b) = mod(avgAng + pi, 2*pi);
            end;

        end;

        phases = newPhases;   % replace array for next iteration

    end; % if -- skips last iteration

end;