function [ ] = envelope( f )
%ENVELOPE Summary of this function goes here
%   Detailed explanation goes here
[~,peaks] = findpeaks(f);
[~,valleys] = findpeaks(-f);
figure();
plot(peaks, f(peaks)); hold on;
plot(valleys, f(valleys));

end

