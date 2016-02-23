% af.m - frequency spectrum analysis.
%   [freqresp: sequence of spectra] = af(input: waveform)
%
% James Salsman, '99-2000; adapted from Malcolm Slaney's mfcc.m [from his 
% "Auditory Toolbox" at www.interval.com.]
%
% Also works when samplingRate is 8000 or fftSize is 256, but I haven't tried 
% both together.

function [freqresp] = af(input)

samplingRate = 16000;                % samples/sec
frameRate = 100;                     % started at 100/sec, maybe use 160 later
windowStep = samplingRate/frameRate; % 16000/100 = 160; shall be even
windowSize = 400;                    % large, but from ETSI standard ES 201108
windowSize = samplingRate/40;        % must be even; use 256 for 11k/s rate 
fftSize = 512;                       % larger better but slower; usually 2^int.

[rows, cols] = size(input);          % sanity-check input dimensions
if (rows > cols)                     % if there are more rows than columns
   input=input';                     % then transpose 
end

% hamming window shapes input signal
%
hamWindow = 0.54 - 0.46*cos(2*pi*(0:windowSize-1)/windowSize);

% how many columns of data we will end up with
%
cols = ceil((length(input) - windowSize)/windowStep);

% Allocate (preextend) all the space we need for the output arrays.
%
freqresp  =  zeros(fftSize/2, cols);
fftData   =  zeros(1, fftSize);
fftMag    =  zeros(1, fftSize);

for c = 1:cols,                             % for each column of data

    first = (c-1)*windowStep + 1;           % determine
    last = first + windowSize - 1;          %   endpoints
    
    % shape the data with a hamming window
    %
    fftData(1:windowSize) = input(first:last).*hamWindow;

    % find the magnitude of the fft
    %
    fftMag = abs(fft(fftData));             % keep only real magnitudes

    freqresp(:, c) = fftMag(1:fftSize/2)';  % emit

end;