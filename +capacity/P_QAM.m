function[p]=P_QAM(X_string,M)
% INPUT: 
%   X_string string of symbols
%   M order of M QAM constellation

% OUPUT:
%   p: probabilities of X

% Generate the M QAM constellation 
MQAM_stars = qammod((0:M-1).',M,'UnitAveragePower',true,'PlotConstellation',false);
                    
% Estimates the probability of the constellation distribution
% P(x)=(number of occurrency of x)/(overall number of symbols)                    
p  = arrayfun(@(x)length(find(X_string==x)), MQAM_stars) / length(X_string);

end