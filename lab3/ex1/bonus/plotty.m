clear all; close all;

load('data');

x = linspace(0,5.1,10000)';
plot(x,data);
hold on

%Reverse engineer

%[-inf,0) = 0

%Fit too [0,2.4] - Straight line
endInd = find(x>2.4,1,'first')-1;
firstX = x(1:endInd);
firstData = data(1:endInd);
p1 = polyfit(firstX,firstData,1);

plot(firstX, p1(1)*firstX + p1(2), 'g')
hold on;
xR = x(endInd+1:end);
dataR = data(endInd+1:end);

%Fit too (2.4, 3.9] - Quadratic or Cubic?
endInd = find(xR>3.9,1,'first')-1;
secX = xR(1:endInd);
secData = dataR(1:endInd);
p2 = polyfit(secX,secData,2);

plot(secX, p2(1)*secX.^2 + p2(2)*secX + p2(3), 'g')
xR = xR(endInd+1:end);
dataR = dataR(endInd+1:end);

%Fit too (3.9,5] - straigt line of value 0.3
endInd = find(xR>5,1,'first')-1;
plot(xR(1:endInd), ones(size(xR(1:endInd)))*0.3, 'g')

% (5, inf] = 0 
