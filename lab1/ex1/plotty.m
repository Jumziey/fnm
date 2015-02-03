clear all; close all;


data = load('afterTransform');
t = 0:1/1024:1-1/1024;
%Plot how it looks initially
plot(t,data);