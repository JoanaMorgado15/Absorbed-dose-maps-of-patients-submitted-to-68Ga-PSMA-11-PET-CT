function [N] = calculateIntegral(halfLife,t0,t1)
%Creation of function calculateIntegral where halfLife, t0 and t1 are
%inputs and N is the output

fun = @(t) exp(-(log(2)/halfLife)*t);
%Acumulated activity is defined as the sum of activity in each moment in
%time and equal to A* (exp(-(log(2)/halfLife)*t)dt integral), since the
%decay constant is equal to log(2)/halfLife.

N = integral(fun,t0,t1);
%Integral of the function fun from t0 to t1