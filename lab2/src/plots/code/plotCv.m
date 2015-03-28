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
		plot(data(:,1),data(:,2), 'color', [1 .5 0])
		title(sprintf('%s %s Heat Capacity', char(gas),prop))
		ylabel('C_v (j/m^3 K)')
		xlabel('Temp (K)')
		j=j+1;
end

filename = sprintf('../%s.tikz',prop);
matlab2tikz(filename, 'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions', ['ticklabel style={font=\tiny},' 'xlabel style={font=\small}']);

figure(2)

%argon
gas = gases(2);
kb = 1.3806488e-23;
rnn = 3.7477e-10;
CvAna = 3*kb*sqrt(2)^3/(2*rnn^3);

dataFile = sprintf('../data/%s%s',char(gas),prop);
data = load(dataFile);
hold on
plot(data(:,1),data(:,2), 'color', [1 .5 0])
plot([0 max(data(:,1))],[1 1]*CvAna, 'b')
title(sprintf('%s %s Heat Capacity', char(gas),prop))
ylabel('C_v (j/m^3 K)')
xlabel('Temp (K)')
legend('Solid State Physics Value', 'Thermal Physics Value', 'Location', 'southeast')

matlab2tikz('../ArCvResult.tikz', 'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions', ['ticklabel style={font=\tiny},' 'xlabel style={font=\small}']);
