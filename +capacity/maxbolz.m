function [px] = maxbolz(C,x,lambda)
%MAX- Summary of this function goes here
%   Detailed explanation goes here
den  = sum(exp(-lambda.*abs(C).^2));
px   = exp(-lambda.*abs(x).^2)./den;
end

