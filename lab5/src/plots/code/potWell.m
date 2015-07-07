close all; clear all;
prob = {'potWell'};
parameters = {'potWell'};
savePlot = true;

pw = 1;
plotWindows = 2;
for para=parameters
	%Need to mark out the potential well in the plot window
	p = importdata(sprintf('../../parameters/%s',char(para)), '=');
	x_min = p.data(2);
	x_max = p.data(3);
	mass = p.data(4);
	nx = p.data(5);
	dt = p.data(6);
	nt = p.data(7);
	
	width = 0.032;
	hwp = ceil((width/((x_max-x_min)/nx))/2);
	potWellStart = ceil(nx/2)-hwp;
	potWellEnd = ceil(nx/2)+hwp;


	%Get ce data
	A = importdata(sprintf('../../wavefunc_%s.dat', char(prob(1))),' ',1);
	d = A.data;
	figure(pw)
	plot(d(:,1), d(:,4))
	hold on
	plot(d(:,1), d(:,3))
	plot([d(potWellStart,1) d(potWellStart,1)], [0 max(d(:,4))], 'color', [1 .5 1]);
	plot([d(potWellEnd,1) d(potWellEnd,1)], [0 max(d(:,4))], 'color', [1 .5 1])
	legend('|\psi_n(x,t)|^2', 'imag(\psi_n(x,t))','Potential Well', 'Location', 'northwest')
	title('Numerical Solution with a Potential Well')
	xlabel('Position')
	
	if savePlot == true
		filename = sprintf('../%sPlot.tikz',char(para));
		matlab2tikz(filename, 'parseStrings', true,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},' 'ticklabel style={font=\tiny}']);
	end
	
	figure(pw+1)
	tr = load(sprintf('../../%s.dat',char(prob)));
	plot(tr(:,1),tr(:,2))
	hold on
	plot(tr(:,1),tr(:,3))
	legend('Reflection', 'Transmission', 'Location','west')
	
	if savePlot == true
		filename = sprintf('../%sTR.tikz',char(para));
		matlab2tikz(filename, 'parseStrings', true,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},' 'ticklabel style={font=\tiny}']);
	end
	
	pw = pw+plotWindows;
end

