%http://www.ee.columbia.edu/ln/rosa/matlab/pvoc/

11[d,sr]=audioread('clar.wav');
sr
% 1024 samples is about 60 ms at 16kHz, a good window
y=pvoc(d,.75,1024);
% Compare original and resynthesis
sound(d,16000)
sound(y,16000) 


%%
[d,sr]=audioread('clar.wav');
[d,sr]=audioread('voix.wav');
[~, channels] = size(d);
d = d(:,1);
e = pvoc(d, 0.5);
f = resample(e,1,2); % NB: 0.8 = 4/5
l = min(length(d), length(f));
soundsc(d(1:l)+f(1:l),sr);

%%
fast = pvoc(d, 1.25, 8192);
deep = resample(fast,5,4); % NB: 0.8 = 4/5

soundsc(fast,sr)
pause(size(fast,1)/sr);
soundsc(deep,sr)

%% Pure sound
sr = 44100;
T = 2;
t = 0:(1/sr):T;
f = 440;
a = 1;
y = a*sin(2*pi*f*t);
sound(y,sr);

%% Moving pure sound
sr = 44100;
T = 0.25;
t = 0:(1/sr):T;
f_start = 30*44100/1024;
f_end = 40*44100/1024;
f = f_start:((f_end-f_start)/(length(t)-1)):f_end;
a = 1;
y = a*sin(2*pi*f.*t);
sound(y,sr);

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


