close all; clear all;

gases = {'Ne', 'Ar', 'Kr', 'Xe'};
prop = 'Cv';

srow = length(gases);
j=1;
for gas=gases
		subplot(srow,1,j)
		dataFile = sprintf('../data/%s%s',char(gas),prop);
		data = load(dataFile);
		hold on
		plot(data(:,1),data(:,2), 'colo
		r', [1 .5 0])
		title(sprintf('%s %s Heat Capacity', char(gas),prop))
		ylabel('C_v')
		xlabel('Temp (K)')
		j=j+1;
end

filename = sprintf('../%s.tikz',prop);
matlab2tikz(filename, 'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['ticklabel style={font=\tiny},' 'xlabel style={font=\small}']);

CvAna=[ 3 0.134
		4 0.347
		5 0.749
		6 1.350
		7 2.13
		8 3.10
		9 4.19
		10 5.36
		11 6.66
		12 7.89
		13 9.08
		14 10.10
		15 11.1
		16 12.0
		17 12.9
		18 13.9
		19 15.0
		20 15.9]
