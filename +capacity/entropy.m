% Entropy: Returns entropy (in bits) of each column of 'X'
%
% H = Entropy(X)
%
% H = row vector of calculated entropies (in bits)
% X = data to be analyzed
% P = Probablity of each symbol in non-decrescent order


function H = entropy(P)

% Establish size of data
[n m] = size(P);

% Housekeeping
H = zeros(1,m);

for j = 1:m
    % Assemble observed alphabet
%     x = unique(X(:,Column));
%     x = sortrows([real(x) imag(x)])*[1;1i];
		
    % Calculate entropy in bits
    H(j) = -sum(P(:,j) .* log2(P(:,j)));
end


