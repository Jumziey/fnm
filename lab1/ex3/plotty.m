clear all; close all;

realTrans = load('realTrans');
imagTrans = load('imagTrans');
trans = (realTrans.^2 + imagTrans.^2);
realFilteredTrans = load('realFilteredTrans');
imagFilteredTrans = load('imagFilteredTrans');
fTrans = (realFilteredTrans.^2 + imagFilteredTrans.^2);

N = 8192;
df = 1/(N);

trans = [trans(N/2+1:N-1);0;trans(1:N/2-1);0];
fTrans = [fTrans(N/2+1:N-1);0;fTrans(1:N/2-1);0];
f = linspace(-1/(2*df), 1/(2*df), N+1);
f = f(1:end-1);

figure(1)
subplot(2,1,1)
plot(f,trans)
title('Frequency Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
%axis([0 4500 0 5e6])
subplot(2,1,2)
plot(f,fTrans)
title('Filtered Frequency Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
%axis([0 4500 0 5e6])
matlab2tikz('freqPowSpec.tikz', 'height', '\figureheight', 'width', '\figurewidth');

figure(2)
subplot(2,1,2)
filtered = load('filteredSignal');
plot(0:1/N:1-1/N,filtered.^2)
title('Filtered Signal')
xlabel('time (s)')
ylabel('Value')

subplot(2,1,1)
init = load('am26.dat');
plot(0:1/N:1-1/N,init.^2);
title('Original Signal')
xlabel('time (s)')
ylabel('Value')
matlab2tikz('signal.tikz', 'height', '\figureheight', 'width', '\figurewidth');

