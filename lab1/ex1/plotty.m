clear all; close all;

dOrig = load('init');
dTrans = load('realTrans');

N = 1024;
t = linspace(0,1,N) ;
%Plot how it looks initially
plot(t,dOrig);
xlabel('Time (s)');
ylabel('Value');
figure(2)
%K vector...? lets see we have n = linspace(-N/2,N/2,N), N = 1024
% Then fn = n/(dt*N);

%Much simpler way to put it is that we go from minimum nyquist freq to max. with N values all and all. 
%Now the minus frequency is the same as the positive counterpart
dt = 1/N;
f = linspace(-1/(2*dt),1/(2*dt), 1024);
%dTrans goes from 0 to 1/(dt*2), -1/(dt*2) to 0, due to fft implementation.
dTrans = [dTrans(N/2:N); dTrans(1:N/2-1)];
dTrans = dTrans - min(dTrans); %Remove offset 
dTrans = dTrans./max(dTrans); %Normalize


plot(f,dTrans)
hold on
%plot([0 0],[0 600], 'r')
xlabel('Frequency (Hz)');
ylabel('Magnitude');
alpha = 16*dt;
Ft = @(f) exp(-(1/4) *(2*pi*alpha*f).^2);
plot(f,Ft(f), '.g')
legend('FFT Data', 'Theoretical values');

FreqRes = 1024/1023
%Only real part