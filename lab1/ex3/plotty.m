clear all; close all;

N = 8192;

orig = load('am26.dat');
%origReal = orig(1:2:N*2);
%origIm = orig(2:2:N*2);
%orig = origReal.^2.*origIm.^2;
realTrans = load('realTrans');
%imTrans = load('imTrans');

t = linspace(0,1,N) ;
%Plot how it looks initially
subplot(2,1,1)

plot(t,orig);
xlabel('Time (s)');
ylabel('Value');
title('Noisy signal')

subplot(2,1,2)
plot(t,orig.^2);
xlabel('Time (s)');
ylabel('Value');
title('Noisy signal squared')


figure(2)
subplot(2,1,1)
%K vector...? lets see we have n = linspace(-N/2,N/2,N), N = 1024
% Then fn = n/(dt*N);

%Much simpler way to put it is that we go from minimum nyquist freq to cax. with N values all and all. 
%Now the minus frequency is the same as the positive counterpart
trans = (realTrans);
dt = 1/N;
f = linspace(-1/(2*dt),1/(2*dt), N);
%trans goes from 0 to 1/(dt*2), -1/(dt*2) to 0, due to fft implementation.
trans = [trans(N/2+1:N); trans(2:N/2); 0]; %This must be cheating
%trans = trans - min(trans); %Remove offset 
%trans = trans./max(trans); %Normalize

plot(f,trans.*trans)
hold on
%plot([0 0],[0 600], 'r')
xlabel('Frequency (Hz)');
ylabel('Power Magnitude (normalized)');
title('Fourier Transform of Signal Data ')

%Filter
fc = 1152;
bandwidth = 1280-895.6;
filter = exp( (-1/2)* ( (abs(f)-fc)./bandwidth).^2);
subplot(2,1,2)
plot(f,(filter'.*trans).^2)
title('Filtered Fourier Transform');

%Filtered Signal 
figure(3)
subplot(2,1,1)
filtered = load('filteredSignal');
plot(t,filtered);
subplot(2,1,2)
plot(t,filtered.^2)
