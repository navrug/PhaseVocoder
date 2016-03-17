%http://www.ee.columbia.edu/ln/rosa/matlab/pvoc/

[y,sr]=audioread('clar.wav');
% 1024 samples is about 60 ms at 16kHz, a good window
z=pvoc(y,0.67,1024);
% Compare original and resynthesis
soundsc(y,sr)
pause(size(y,1)/sr);
soundsc(z,sr) 


%%
[y,sr]=audioread('voix.wav');
[~, channels] = size(y);
y = y(:,1);
e = pvoc(y, 0.5);
f = resample(e,1,2); % NB: 0.8 = 4/5
l = min(length(y), length(f));
soundsc(y(1:l)+f(1:l),sr);

%%
fast = pvoc(d, 1.25, 8192);
deep = resample(fast,5,4); % NB: 0.8 = 4/5

soundsc(fast,sr)
pause(size(fast,1)/sr);
soundsc(deep,sr)

%% Pure sound
sr = 44100;
T = 1.6254;
t = 0:(1/sr):T;
f = 880;
a = 1;
y = a*sin(2*pi*f*t);
soundsc(y,sr);

%% Moving pure sound
sr = 44100;
T = (80*(1024-128)+128)/44100;
t = 0:(1/sr):T;
f_start = 30*44100/1024;
f_end = 40*44100/1024;
f = f_start:((f_end-f_start)/(length(t)-1)):f_end;
a = 1;
y = a*sin(2*pi*f.*t);
soundsc(y,sr);

%%
e = pvoc(y, 2);
%f = resample(e,2,1); % NB: 0.8 = 4/5
soundsc(e,sr);

%% Envelope
envelope(e);

%% Phase
stft_y = stft(y);
[~,max_stft_y] = max(abs(stft_y),[],1);
idx = sub2ind(size(stft_y'), 1:length(max_stft_y), max_stft_y);

[n_channels, n_frames] = size(stft_y);
figure();
plot(1:n_frames,angle(stft_y(idx))); hold on;
plot(max_stft_y)


%% From scratch for alpha = 2

%% Pure average
n=1024;
hop=256;
alpha=1.4;
y_spec = stft(y,n,n,hop);
[n_channels, n_frames] = size(y_spec);
z_spec = zeros(n_channels, 2*n_frames-1);
z_spec(:,2*(1:n_frames)-1) = y_spec;
z_spec(:,2*(1:(n_frames-1))) = 0.5*(y_spec(:,1:(n_frames-1)) + y_spec(:,2:n_frames));
z = istft(z_spec,n,n,hop);
soundsc(z,sr);

%% Average of amplitude, average of phase
y_spec = stft(y,n,n,hop);
[n_channels, n_frames] = size(y_spec);
z_spec = zeros(n_channels, 2*n_frames-1);
z_spec(:,2*(1:n_frames)-1) = y_spec;
z_spec(:,2*(1:(n_frames-1))) = 0.5*(abs(y_spec(:,1:(n_frames-1))) + abs(y_spec(:,2:n_frames)));
angles = 0.5*(angle(y_spec(:,1:(n_frames-1))) + angle(y_spec(:,2:n_frames)));
z_spec(:,2*(1:(n_frames-1))) = z_spec(:,2*(1:(n_frames-1))) .* exp(1i*angles);
z = istft(z_spec,n,n,hop);
soundsc(z,sr);


%% Phase vocoder
n=1024;
hop=256;
alpha=1.4;
Ra = 1/sr;
Rs = alpha*Ra;
y_spec = stft(y,n,n,hop);
y_abs = abs(y_spec);
y_angles = angle(y_spec);
[n_channels, n_frames] = size(y_spec);

% Phase unwrapping
omega = 2*pi/(n/2) * (0:(n/2))';
delta_phi = zeros(size(y_angles));
delta_phi(:,2:n_frames) = y_angles(:,2:n_frames) - y_angles(:,1:n_frames-1);
delta_phi(:,2:n_frames) = bsxfun(@minus, delta_phi(:,2:n_frames), Ra*omega);
delta_p_phi = mod(delta_phi + pi, 2*pi) - pi;

% Basic vodocer phase synthesis
z_angles = alpha * ( bsxfun(@minus, y_angles, y_angles(:,1)) + cumsum(delta_p_phi-delta_phi,2));
z_angles =  bsxfun(@plus, z_angles, alpha*y_angles(:,1));

% Loose phase locking synthesis
z_spec = abs(y_spec) .* exp(1i * z_angles);
loose_lock = zeros(size(z_spec));
loose_lock(1,:) = z_spec(2,:) + z_spec(1,:);
loose_lock(2:n_channels-1,:) = z_spec(3:n_channels,:) + z_spec(2:n_channels-1,:) + z_spec(1:n_channels-2,:);
loose_lock(n_channels,:) = z_spec(n_channels,:) + z_spec(n_channels-1,:);

% Rigid phase locking synthesis
% Find peaks in a frame
peaks = cell(n_frames,1);
for frame = 1:n_frames
   [pks, locs] = findpeaks(abs(y_spec(:,frame)),'MinPeakProminence',2);
   peaks{frame} = locs;
end
peaks = cell(n_frames,1);
for frame = 1:n_frames
   locs = [];
   if y_abs(1,frame) > y_abs(2,frame) && y_abs(1,frame) > y_abs(3,frame)
      locs = [1];
   end
   if y_abs(2,frame) > y_abs(1,frame) && y_abs(2,frame) > y_abs(3,frame) ...
      && y_abs(2,frame) > y_abs(4,frame)
      locs = [locs 2];
   end
   for i=3:n_channels-2
      if y_abs(i,frame) > y_abs(i-1,frame) && y_abs(i,frame) > y_abs(i-2,frame) ...
         && y_abs(i,frame) > y_abs(i+1,frame) && y_abs(i,frame) > y_abs(i+2,frame)
         locs = [locs i];
      end
   end
   if y_abs(n_channels-1,frame) > y_abs(n_channels-3,frame) && y_abs(n_channels-1,frame) > y_abs(n_channels-2,frame) ...
      && y_abs(n_channels-1,frame) > y_abs(n_channels,frame)
      locs = [locs n_channels-1];
   end
   if y_abs(n_channels,frame) > y_abs(n_channels-2,frame) && y_abs(n_channels,frame) > y_abs(n_channels-1,frame)
      locs = [locs n_channels];
   end
   peaks{frame} = locs;
end

% Find closest peak for each channel in each frame
closest = zeros(size(z_spec));
for frame = 1:n_frames
   if length(peaks{frame}) == 1
      closest(:,frame) = peaks{frame}(1);
   else
      closest(1:floor(0.5*(peaks{frame}(1)+peaks{frame}(2))),frame) = peaks{frame}(1);
      for k = 2:length(peaks{frame})-1
         win_start = floor(0.5*(peaks{frame}(k-1)+peaks{frame}(k)))+1;
         win_end = floor(0.5*(peaks{frame}(k)+peaks{frame}(k+1)));
         closest(win_start:win_end , frame) = peaks{frame}(k);
      end
      win_start = floor(0.5*(peaks{frame}(length(peaks{frame})-1)+peaks{frame}(length(peaks{frame}))));
      closest(win_start:n_channels,frame) = peaks{frame}(length(peaks{frame}));
   end
end

% Aligning phase to that of closest peak
idx = sub2ind(size(z_angles), closest,repmat(1:n_frames,n_channels,1));
rigid_lock = z_angles(idx) + y_angles - y_angles(idx);

%%
z_spec = abs(y_spec) .* exp(1i * z_angles);
z = istft(z_spec,n,n,floor(alpha*hop));
soundsc(z,sr);
%%
z_spec = abs(y_spec) .* exp(1i * angle(loose_lock));
z = istft(z_spec,n,n,floor(alpha*hop));
soundsc(z,sr);

%%
z_spec = abs(y_spec) .* exp(1i * rigid_lock);
z = istft(z_spec,n,n,floor(alpha*hop));
soundsc(z,sr);

