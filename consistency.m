function [ Dm ] = consistency( y_stft,  ftsize, w, h )
%RECONT Summary of this function goes here
%   Detailed explanation goes here
   z = istft(y_stft, ftsize, w, h);
   z_stft = stft(z, ftsize, w, h);
   %ratio = (abs(z_stft) - abs(y_stft)).^2/abs(y_stft).^2;
   % We ignore a few frames at begininng and end to avoid edge effect.
   P = 35;
   [n_channels, n_frames] = size(y_stft);
   %Dm = sum(sum(ratio(:,(1+P):(n_frames-P))));
   num = 0;
   den = 0;
   for i = 1:n_channels
      for j = (1+P):(n_frames-P)
         num = num + (abs(z_stft(i,j)) - abs(y_stft(i,j))).^2;
         den = den + abs(y_stft(i,j)).^2;
      end
   end
   Dm = num/den;
end

