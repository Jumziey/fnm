close all; clear all;
dataList = {'forwardData' 'backwardData'};
titles = {'Top Spin [0,4] (s)' 'Top Spin [4,0] (s)'};

i=1;j=1;
for d=dataList
	figure(i)
	ax = gca;
	dataFile = sprintf('../data/%s',char(d));
	data = load(dataFile);
	time = linspace(0,4,length(data));
	
	plot(time,data(:,1)*180/pi)
	hold on
	plot(time,data(:,2)*180/pi)
	plot(time,data(:,3)*180/pi)
	
	if i == 1
		legend('\phi','\psi','\theta', 'Location', 'northwest')
	else
		legend('\phi','\psi','\theta')
	end
	if j ~= 1
		ax.XTickLabel = flip(ax.XTickLabel);
	end
	title(char(titles(j)))
	xlabel('Time (s)')
	ylabel('Degrees')
	
	filename = sprintf('../%s.tikz',char(d));
	matlab2tikz(filename, 'parseStrings', true,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},'...
                       					'ticklabel style={font=\tiny}']);
                       						
	i = i+2;
	j = j+1;
end
