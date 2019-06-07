classdef capacity_functions
    %CAPACITY_FUNCTIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant = true )
        
        xi  = [-3.436159118837737603327; -2.532731674232789796409; ...
               -1.756683649299881773451; -1.036610829789513654178; ...
               -0.3429013272237046087892; 0.3429013272237046087892;...
               1.036610829789513654178;   1.756683649299881773451; ...
               2.532731674232789796409;   3.436159118837737603327];
        w   = [7.64043285523262062916E-6; 0.001343645746781232692202; ...
               0.0338743944554810631362;  0.2401386110823146864165;   ...
               0.6108626337353257987836;  0.6108626337353257987836;   ...
               0.2401386110823146864165;  0.03387439445548106313617;  ...
               0.001343645746781232692202; 7.64043285523262062916E-6];
           
    end
    
    methods(Static)
        
        function Es = symbol_energy(C,Pk)
            %SYMBOL_ENERGY return the average symbol energy
            %   ES = SYMBOL_ENERGY(C,PK) return the average symbol energy
            %   of the vector C given the vector of the probabilities Pk
            
            Es = sum(Pk.*abs(C).^2);
        end
        
        function b0 = get_mqamsymbs_zero_kth_bit(k,M)
            %GET_MQAMSYMBS_ZERO_KTH_BIT inserts a binary zero inside a
            %positive binary number at a given position.
            %   B0 = GET_MQAMSYMBS_ZERO_KTH_BIT(K,M) return a vector of
            %   numbers in base 10 that has a zero in the left-msb bit when
            %   coded in base 2.
            
            m  = log2(M);
            b0 = de2bi((0:2^(m-1)-1),m,'left-msb');
            b0(:,[1,k]) = b0(:,[k,1]);
            b0 = bi2de(b0,'left-msb');
        end
        
        function mi  = qam_mi(C,M,sg,Pk)
            %QAM_GMI evaluates the Mutual information (MI) for M-QAM
            %using Gauss-Hermite quadrature and assuming an AWGN channel.
            %   MI = QAM_MI(C,M,SG,PK) evaluates the mi for a M-QAM
            %   constallation points C whit propabilities Pk assuming a
            %   noise deviation of SG.
            
            import capacity.*;
            
            xi  = capacity_functions.xi;
            w   = capacity_functions.w;
            
            mi = 0.0;
            sg2 = sg^2;
            %  Cycle through constellation points
            for i = 1:M
                dij   = (C(i) - C).';
                for l=1:length(xi)
                    c_xi    = xi(l)+1i*xi;
                    tmp = sum(exp((-abs(dij).^2 - 2*sg.*real(c_xi.*dij))/sg2).*Pk.',2);
                    mi  = mi - w(l)*sum(w.*log2(tmp).*Pk(i));
                end
            end
            
            mi = mi/pi;
        end
        
        function mi  = qam_montecarlo_mi(z,C,M,sg,Pk)
            %QAM_MONTECARLO_MI evaluates the Mutual information (MI) for M-QAM
            %using Monte-Carlo integration and assuming an AWGN channel.
            %   MI = QAM_MONTECARLO_MI(Z,C,M,SG,PK) given the received noise Z
            %   (transmitted-received symbols), evaluates the mi for M-QAM
            %   constellation points C whit propabilities Pk assuming a
            %   noise deviation of SG.
            
            mi  = 0.0;
            sg2 = sg^2;
            D   = length(z);
            
            %  Cycle through constellation points
            for i = 1:M
                dij   = (C(i) - C).';
                tmp = sum(exp((-abs(dij).^2 - 2*real(z.*dij))/sg2).*Pk.',2);
                mi  = mi - sum(log2(tmp).*Pk(i));
            end
            
            mi = mi/D;
        end
        
        function mi  = pam_mi(C,M,sg,Pk)
            %QAM_GMI evaluates the Mutual information (MI) for M-QAM
            %using Gauss-Hermite quadrature and assuming an AWGN channel.
            %   MI = QAM_MI(C,M,SG,PK) evaluates the mi for a M-QAM
            %   constallation points C whit propabilities Pk assuming a
            %   noise deviation of SG.
            
            import capacity.*;
            
            xi  = capacity_functions.xi;
            w   = capacity_functions.w;
            
            mi = 0.0;
            sg2 = sg^2;
            %  Cycle through constellation points
            for i = 1:M
                dij   = (C-C(i)).';
                for l=1:length(xi)                    
                    tmp = sum(exp((-dij.^2 - sqrt(8)*sg.*xi(l).*dij)/(2*sg2)).*Pk.',2);
                    mi  = mi - w(l)*log2(tmp)*Pk(i);
                end
            end
            
            mi = mi/sqrt(pi);
        end
        
        function gmi = pam_gmi(C,M,sg,Pk)
            %PAM_GMI evaluates the BICM Mutual information (GMI) for PAM
            %using Gauss-Hermite quadrature and assuming an AWGN channel.
            %   GMI = PAM_GMI(C,M,SG,PK) evaluates the gmi for a M-PAM
            %   constellation points C whit propabilities Pk assuming a
            %   noise deviation of SG.
            
            import capacity.*;
            
            xi  = capacity_functions.xi;
            w   = capacity_functions.w;
            
            gmi = 0.0;
            m   = log2(M);
            sg2 = sg^2;
            
            %  Cycle through constellation bit
            for k = 1:m
                % Cycle through binary values
                for b = 0:1
                    bi  = (capacity_functions.get_mqamsymbs_zero_kth_bit(k,M)+b*2^(m-(k-1)-1))+1;
                    bi  = intersect(bi,find(Pk));
                    bj  = bi;
                    Pbk(1+b,k) = sum(Pk(bi));
                    % Cycle through constellation points where k-th bit is
                    % equal to b
                    for i = 1:length(bi)
                        dip   = (C(bi(i)) - C(Pk~=0)).';
                        dij   = (C(bi(i)) - C(bj)).';
                        for l=1:length(xi)
                            tmp_num = sum(exp((-abs(dip).^2-sqrt(8)*sg.*xi(l).*dip)/(2*sg2)).*Pk(Pk~=0).',2);
                            tmp_den = sum(exp((-abs(dij).^2-sqrt(8)*sg.*xi(l).*dij)/(2*sg2)).*(Pk(bj)/Pbk(1+b,k)).',2);
                            gmi = gmi - w(l)*log2(tmp_num./tmp_den).*Pk(bi(i));
                        end
                    end
                end
            end
            
            gmi = gmi/sqrt(pi);
            
            % Add the entropy of the constellation
            gmi = gmi - sum(Pk(Pk~= 0).*log2(Pk(Pk~=0)));
            
            % Remove the entropy of each bit
            gmi = gmi + sum(Pbk(1,:).*log2(Pbk(1,:)));
            gmi = gmi + sum(Pbk(2,:).*log2(Pbk(2,:)));
        end
        
        
        function gmi = qam_gmi(C,M,sg,Pk)
            %QAM_GMI evaluates the BICM Mutual information (GMI) for QAM
            %using Gauss-Hermite quadrature and assuming an AWGN channel.
            %   GMI = QAM_GMI(C,M,SG,PK) evaluates the gmi for a M-QAM
            %   constallation points C whit propabilities Pk assuming a
            %   noise deviation of SG.
            
            import capacity.*;
            
            xi  = capacity_functions.xi;
            w   = capacity_functions.w;
            
            gmi = 0.0;
            m   = log2(M);
            sg2 = sg^2;
            
            %  Cycle through constellation bit
            for k = 1:m
                % Cycle through binary values
                for b = 0:1
                    bi  = (capacity_functions.get_mqamsymbs_zero_kth_bit(k,M)+b*2^(m-(k-1)-1))+1;
                    bi  = intersect(bi,find(Pk));
                    bj  = bi;
                    Pbk(1+b,k) = sum(Pk(bi));
                    % Cycle through constellation points where k-th bit is
                    % equal to b
                    for i = 1:length(bi)
                        dip   = (C(bi(i)) - C(Pk~=0)).';
                        dij   = (C(bi(i)) - C(bj)).';
                        for l=1:length(xi)
                            c_xi    = xi(l)+1i*xi;
                            tmp_num = sum(exp((-abs(dip).^2-2*sg.*real(c_xi.*dip))/sg2).*Pk(Pk~=0).',2);
                            tmp_den = sum(exp((-abs(dij).^2-2*sg.*real(c_xi.*dij))/sg2).*(Pk(bj)/Pbk(1+b,k)).',2);
                            gmi = gmi - w(l)*sum(w.*log2(tmp_num./tmp_den).*Pk(bi(i)));
                        end
                    end
                end
            end
            
            gmi = gmi/pi;
            
            % Add the entropy of the constellation
            gmi = gmi - sum(Pk(Pk~= 0).*log2(Pk(Pk~=0)));
            
            % Remove the entropy of each bit
            gmi = gmi + sum(Pbk(1,:).*log2(Pbk(1,:)));
            gmi = gmi + sum(Pbk(2,:).*log2(Pbk(2,:)));
        end
        
        function gmi = qam_montecarlo_gmi(z,C,M,sg,Pk)
            %QAM_MONTECARLO_GMI evaluates the BICM Mutual information (GMI) for QAM
            %using Monte-Carlo integration and assuming an AWGN channel.
            %   GMI = QAM_MONTECARLO_GMI(Z,C,M,SG,PK) given the received noise Z
            %   (transmitted-received symbols), evaluates the gmi for M-QAM
            %   constallation points C whit propabilities Pk assuming a
            %   noise deviation of SG.
            
            import capacity.*;
            
            gmi = 0.0;
            m   = log2(M);
            sg2 = sg^2;
            D   = length(z);
            %  Cycle through constellation bit
            for k = 1:m
                % Cycle through binary values
                for b = 0:1
                    bi  = (capacity_functions.get_mqamsymbs_zero_kth_bit(k,M)+b*2^(m-(k-1)-1))+1;
                    bi  = intersect(bi,find(Pk));
                    bj  = bi;
                    Pbk(1+b,k) = sum(Pk(bi));
                    % Cycle through constellation points where k-th bit is
                    % equal to b
                    for i = 1:length(bi)
                        dip   = (C(bi(i)) - C(Pk~=0)).';
                        dij   = (C(bi(i)) - C(bj)).';
                        tmp_num = sum(exp((-abs(dip).^2-2.*real(z.*dip))/sg2).*Pk(Pk~= 0).',2);
                        tmp_den = sum(exp((-abs(dij).^2-2.*real(z.*dij))/sg2).*(Pk(bj)/Pbk(1+b,k)).',2);
                        gmi   = gmi - sum(log2(tmp_num./tmp_den))/D.*Pk(bi(i));
                    end
                    
                end
            end
            
            % Add the entropy of the constellation
            gmi = gmi - sum(Pk(Pk~= 0).*log2(Pk(Pk~= 0)));
            
            % Remove the entropy of each bit
            gmi = gmi + sum(Pbk(1,:).*log2(Pbk(1,:)));
            gmi = gmi + sum(Pbk(2,:).*log2(Pbk(2,:)));
            
        end
    end
end

