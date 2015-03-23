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

