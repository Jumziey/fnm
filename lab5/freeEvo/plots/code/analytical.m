
time = 8e-4;
tvec = linspace(0,time,500);
xvec = linspace(0,1,1000);

s0 = 1/40;
k0 = 200;
x0 = 1/4;
phi = @(theta,t) (-theta-k0^2*t)*i

p1 = (s0^2/pi)^(1/4);
p2 = @(theta,t) exp(phi(theta,t))/sqrt((s0^2+2*t*i));
p3 = @(x) exp(k0*i*x);
p4 = @(x,t) exp(-(x-x0-2*k0*t).^2/(2*s0^2+4*i*t)) 

wave = @(x,t,theta) p1*p2(theta,t)*p3(x).*p4(x,t)


for t = time
	theta = atan(2*t/s0)/2;
	figure(1)
	hold off
	plot(xvec,abs(wave(xvec,t,theta)).^2)
	hold on
	plot(xvec, imag(wave(xvec,t,theta)))
	title('Analytical solution')
	legend('|\psi(x,t)|^2', 'imag(\psi(x,t))')
	axis([0 1 -5 25])
	drawnow
end
