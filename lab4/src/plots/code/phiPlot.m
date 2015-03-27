close all; clear all;

phiAn = @(r) 1.0-0.5*(r+2.0).*exp(-r);

dataFiles = {'phi1', 'phi2', 'phi3'};
titles = {'$\varphi(r)$ using analytical for $\varphi_1$', '$\varphi(r)$ using numerical integration for $\varphi_1$', '$\varphi(r)$ using 95\% of analytical for $\varphi_1$'};

%Plot of initial values
i=1;
for dataFile=dataFiles
	figure(i);
	data = load(sprintf('../data/%s',char(dataFile)));
	plot(data(:,1),data(:,2))
	hold on
	plot(data(:,1),data(:,3))
	%axis([0 20 -0.1 1.2]) 
	ylabel('$\varphi(r)$ (statvolt/cm)')
	xlabel('r (cm)')
	title(titles(i))
	legend('$\varphi(r)$', 'error')
	i = i+1;
	filename = sprintf('../%s.tikz',char(dataFile));
	matlab2tikz(filename, 'parseStrings', false,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},' 'ticklabel style={font=\tiny}']);
end

figure(i)
%We now wanna remove the 1st order fit from phi3
p = polyfit(data(end-10:end,1),data(end-10:end,2),1)
plot(data(:,1), data(:,2)-p(1)*data(:,1))
hold on
plot(data(:,1), phiAn(data(:,1)) - (data(:,2)-p(1)*data(:,1)))
	legend('$\varphi(r)$', 'error')

ylabel('$\varphi$ (statvolt/cm)')
xlabel('r (cm)')
title('Fixed $\varphi(r)$ with linear regression')
	matlab2tikz('../phi3reg.tikz', 'parseStrings', false,'height', '\figureheight', 'width', '\figurewidth', 'extraaxisoptions',['title style={font=\small},' 'ticklabel style={font=\tiny}']);
	

xVal = linspace(2,20, 10);
for x=xVal
	ind = find(data(:,1)==x);
	disp(sprintf('%g & %g & %g \\\\ \\hline',x, data(ind,2)-p(1)*x, phiAn(x) - (data(ind,2)-p(1)*x)))
end


