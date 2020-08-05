close all; clear all; clc
global t

%Double check which parameter you will be working with
%Change the parameters that will be fitted

t = 1:432; %Maximum amount of time - 18 days



filename = 'Data.csv';
delimiterIn = ',';
headerlinesIn = 1;
Testdata = importdata(filename,delimiterIn,headerlinesIn);
% 
%432 is the biggest number on the hours
%Should be set equal to t

