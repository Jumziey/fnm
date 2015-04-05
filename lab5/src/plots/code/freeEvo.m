close all; clear all;
prob = {'freeEvo'};
parameters = {'biggestStep', 'periodicInf', 'smallError'};
savePlot = true;

pw = 1;
plotWindows = 2;
for para=parameters
	%Load parameters for analytical solution
	p = importdata(sprintf('../../parameters/%s',char(para)), '=');
	x_min = p.data(1);
	x_max = p.data(2);
	mass = p.data(3);
	nx = p.data(4);
	dt = p.data(5);
	nt = p.data(6);
	cd ../..
	[status, cmdout] = system(sprintf('./%sWave parameters/%s',char(prob),char(para)));
	cd plots/code
	if status ~= 0
		disp(cmdout)
		return;
	end

	%Constants
	s0 = 1/40;
	k0 = 200;
	x0 = 1/4;
	t = dt*nt;
	phi = @(theta,t) (-theta-k0^2*t)*i;
	xvec = linspace(x_min,x_max,nx);
	theta = atan(2*t/s0)/2;
	%Defining analytical solution
	p1 = (s0^2/pi)^(1/4);
	p2 = @(theta,t) exp(phi(theta,t))/sqrt((s0^2+2*t*i));
	p3 = @(x) exp(k0*i*x);
	p4 = @(x,t) exp(-(x-x0-2*k0*t).^2/(2*s0^2+4*i*t));

	wave = @(x,t,theta) p1*p2(theta,t)*p3(x).*p4(x,t);

	%Get ce data
	A = importdata('../../wavefunc.dat',' ',1);
	d = A.data;
	%Analytical over Numerical
	figure(pw)
	plot(d(:,1), d(:,4))
	hold on
	plot(d(:,1), d(:,3))

	%We create 125 dots on the line, seem good
	points = 60;
	jump = ceil(length(d)/points);
	plot(xvec(1:jump:end),abs(wave(xvec(1:jump:end),t,theta)).^2, '.r', 'markersize', 4)
	plot(xvec(1:jump:end), imag(wave(xvec(1:jump:end),t,theta)), '.black', 'markersize', 4)
	legend('|\psi_n(x,t)|^2', 'imag(\psi_n(x,t))','|\psi_a(x,t)|^2', 'imag(\psi_a(x,t))', 'Location', 'northwest')
	title('Analytical Solution over Numerical Solution')
	xlabel('Position')
	
	if savePlot == true
		filename = sprintf('../%sPlot.tikz',char(para));
		matlab2tikz(filename, 'parseStrings', true,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},' 'ticklabel style={font=\tiny}']);
	end
	%Error plot
	figure(pw+1)
	plot(d(:,1), abs(wave(d(:,1),t,theta)).^2 - d(:,4))
	title(sprintf('Error of |\\psi_n(x,t)|^2, dt = %g, dx = %g',dt,(x_max-x_min)/nx))
	ylabel('Magnitude')
	xlabel('Position')
	
	if savePlot == true
		filename = sprintf('../%sError.tikz',char(para));
		matlab2tikz(filename, 'parseStrings', true,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},' 'ticklabel style={font=\tiny}']);	
	end
	pw = pw+plotWindows;
end

