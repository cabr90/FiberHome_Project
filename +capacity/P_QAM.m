function[p]=P_QAM(X,M)
% INPUT: 
%   X_string string of symbols
%   M order of M QAM constellation

% OUPUT:
%   p: probabilities of X

[n,d] = size(X);

% Generate the M QAM constellation 
MQAM_stars = qammod((0:M-1).'.*ones(M,d),M,'UnitAveragePower',true,'PlotConstellation',false);

if d>1
    Ctmp = [];
    for i = 1:M
        for j = 2:d
            Ctmp = [Ctmp;ones(M,1)*MQAM_stars(i,1:j-1),MQAM_stars(:,j:end)];
        end
    end
    MQAM_stars  = Ctmp;
end


% Estimates the probability of the constellation distribution
% P(x)=(number of occurrency of x)/(overall number of symbols)                    

% pt   = arrayfun( @(x,y) length(find( all( eq( X,[x,y] ),2 ) ) ), MQAM_stars(:,1),MQAM_stars(:,2) )/n;

p  = cellfun( @(x) length(find( all( eq( X,x ),2 ) ) ), num2cell(MQAM_stars,2))/n;


end