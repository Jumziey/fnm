close all; clear all;

NeFreq200
NeGamma200

gases = {'Ne', 'Ar', 'Kr', 'Xe'};
directions = [200 220 222];
res = 200;
prop = 'Freq';

scol = length(directions);
srow = length(gases);
i=1;j=1;
for gas=gases
	for dir=directions
		subplot(srow,scol,j)
		dataFile = sprintf('../data/%s%s%d',char(gas),prop,dir);
		data = load(dataFile);
		
		%making less markers
		data4 = data(:,4);
		x = linspace(1,res,length(data4));
		plot(x(1:10:end),data4(1:10:end),'xm')
		hold on
		%plus the line
		plot(data4, 'm')
		
		hold on
		plot(data(:,5), 'color', [1 .5 0])
		plot(data(:,6), 'black')
		title(sprintf('%s %s in [%d]', char(gas),prop, dir/2))
		
		hdir=sprintf('[%d]',dir/2);
		ddir=sprintf('[%d]', dir);
		ylabel('$\omega(\mathbf{q})$ (rad/s)')
		set(gca,'XTickLabel',{'[000]' hdir ddir})
		set(gca, 'XTick', [0 res/2 res])

		j=j+1;
	end
	i=i+1;
end




filename = sprintf('../%s.tikz',prop);
matlab2tikz(filename, 'parseStrings', false,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},'...
                       					'ticklabel style={font=\tiny}']);

