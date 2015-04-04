clear all; close all;

filename = '../wavefunc.dat';
A = importdata(filename,' ',1);
d = A.data;


plot(d(:,1), d(:,4))
hold on
plot(d(:,1), d(:,3))

title('Numerical solution')
%legend('|\psi(x,t)|^2', 'imag(\psi(x,t))')
