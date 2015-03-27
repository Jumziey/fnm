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
	
	plot(time,data(:,4)*180/pi)
	plot(time,data(:,5)*180/pi)
	plot(time,data(:,6)*180/pi)
	if i == 1
		legend('\phi','\psi','\theta', 'c1', 'c2', 'Energy Diff', 'Location', 'northwest')
	else
		legend('\phi','\psi','\theta', 'c1', 'c2', 'Energy Diff')
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
                       						
	figure(i+1)
	ax = gca;
	plot(time,data(:,4)*180/pi)
	hold on
	plot(time,data(:,5)*180/pi)
	plot(time,data(:,6)*180/pi)
	if j ~= 1
		ax.XTickLabel = flip(ax.XTickLabel);
	end
	title(sprintf('%s Errors',char(titles(j))));
	xlabel('Time (s)')
	ylabel('Degrees')
	if j == 1
		legend('c1', 'c2', 'Energy Diff', 'Location', 'northwest')
	else
		legend('c1', 'c2', 'Energy Diff')
	end
	filename = sprintf('../%sError.tikz',char(d));
	matlab2tikz(filename, 'parseStrings', true,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},'...
                       					'ticklabel style={font=\tiny}']);
	i = i+2;
	j = j+1;
end
