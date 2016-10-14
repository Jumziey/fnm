close all; clear all;
prob = {'trvk0'};
parameters = {'trvk0'};
savePlot = true;

pw = 1;
for para=parameters

	%Get ce data
	d = load('../../potWellk0.dat');
	k0 = d(:,1);
	T = d(:,2);
	R = d(:,3);
	figure(pw)
	plot(k0,T)
	hold on
	plot(k0,R)
	legend('Transmission', 'Reflection', 'Location', 'east')
	title('T/R Vs k0')
	xlabel('k0')
	ylabel('Reflection/Transmission Magnitude')

	if savePlot == true
		filename = sprintf('../%sPlot.tikz',char(para));
		matlab2tikz(filename, 'parseStrings', true,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},' 'ticklabel style={font=\tiny}']);
	end
end
