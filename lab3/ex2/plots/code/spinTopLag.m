close all; clear all;
dataList = {'forwardData' 'backwardData'};
titles = {'Top Spin [0,4] (s)' 'Top Spin [4,0] (s)'};
ps = {'bx', 'rd', 'gs', 'mo'};
createImage = false;

i=1;j=1;
for d=dataList
	figure(i)
	ax = gca;
	dataFile = sprintf('../data/%s',char(d));
	data = load(dataFile);
	time = linspace(0,4,length(data));
	
	subplot(2,1,1)
	plot(time,data(:,2)*180/pi)
	subplot(2,1,2)
	plot(time,data(:,1)*180/pi)
	hold on
	plot(time,data(:,3)*180/pi)
	
	subplot(2,1,1)
	if i == 1
		legend('\psi', 'Location', 'northwest')
	else
		legend('\psi')
	end
	if j ~= 1
		ax = gca;
		ax.XTickLabel = flip(ax.XTickLabel);
	end
	title(char(titles(j)))
	xlabel('Time (s)')
	ylabel('Degrees')
	
	subplot(2,1,2)
	if i == 1
		legend('\phi','\theta', 'Location', 'southwest')
	else
		legend('\phi','\theta')
	end
	if j ~= 1
		ax = gca;
		ax.XTickLabel = flip(ax.XTickLabel);
	end
	title(char(titles(j)))
	xlabel('Time (s)')
	ylabel('Degrees')
	
	if(createImage)
		filename = sprintf('../%s.tikz',char(d));
		matlab2tikz(filename, 'parseStrings', true,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},'...
		                     					'ticklabel style={font=\tiny}']);
	end
                       						
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
	%ylabel('Degrees')
	if j == 1
		legend('c1', 'c2', 'Energy Diff', 'Location', 'northwest')
	else
		legend('c1', 'c2', 'Energy Diff')
	end
	
	if(createImage)
		filename = sprintf('../%sError.tikz',char(d));
		matlab2tikz(filename, 'parseStrings', true,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},'...
		                     					'ticklabel style={font=\tiny}']);
	end
	i = i+2;
	j = j+1;
end

fdata = load('../data/forwardData');
bdata = load('../data/backwardData');
figure(i+2)
bdata(:,6) = flip(bdata(:,6));
plot(time,abs(fdata(2:end-1,6) - bdata(:,6)), 'x')
title('Difference in energy error')
ylabel('Energy Error Difference')
xlabel('Time (s)')



