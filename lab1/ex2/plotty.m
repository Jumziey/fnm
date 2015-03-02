clear all; close all;

orig = load('init');
realTrans = load('realTrans');
imTrans = load('imTrans');

N = 1024;
t = linspace(0,1,N) ;
%Plot how it looks initially
subplot(2,1,1)

plot(t,orig);
xlabel('Time (s)');
ylabel('Value');
title('0011001100110011 signal')

subplot(2,1,2)
%K vector...? lets see we have n = linspace(-N/2,N/2,N), N = 1024
% Then fn = n/(dt*N);

%Much simpler way to put it is that we go from minimum nyquist freq to cax. with N values all and all. 
%Now the minus frequency is the same as the positive counterpart
trans = (realTrans.*realTrans + imTrans.*imTrans);
dt = 1/N;
f = linspace(-1/(2*dt),1/(2*dt), 1025);
f = f(1:end-1)
%trans goes from 0 to 1/(dt*2), -1/(dt*2) to 0, due to fft implementation.
trans = [trans(N/2+1:N); trans(2:N/2); 0]; %This must be cheating
%trans = trans - min(trans); %Remove offset 
trans = trans./max(trans); %Normalize




plot(f(N/2:N),trans(N/2:N))
hold on
%plot([0 0],[0 600], 'r')
xlabel('Frequency (Hz)');
ylabel('Power Magnitude (normalized)');
title('Fourier Transform of Signal Data u1 = 10')
axis([-1 100 0 1])
%Only real part
