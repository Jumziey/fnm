close all; clear all;
dataList = {'forwardData' 'backwardData'};
titles = {'Top Spin [0,4] (s)' 'Top Spin [4,0] (s)'};

i=1;j=1;
d=dataList(1);
	figure(i)
	ax = gca;
	dataFile = sprintf('../data/%s',char(d));
	data = load(dataFile);
	time = linspace(0,4,length(data));
	
	plot(time,data(:,1))
	hold on
	plot(time,data(:,2))
	plot(time,data(:,3))
	
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
	
                    						
	i = i+2;
	j = j+1;

