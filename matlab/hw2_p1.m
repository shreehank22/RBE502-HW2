% Homework 2, Problem 1

clear; clc; syms x1 x2 real; % defining the variables

f1 = [-x1 +x1*x2;-x2];

f2 = [-x2-x1*(1 - x1^2 - x2^2);x1 - x2*(1 - x1^2 - x2^2)];

f3 = [ x2*(1-x1^2);-(x1+x2)*(1-x1^2)];

f4 = [-x1-x2;2*x1-x2^3];

