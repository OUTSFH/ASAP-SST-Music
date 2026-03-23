% function [h, Dh, tt] = hermf(N, M, tm)
%     if mod(N, 2) == 0
%         error('N must be an odd integer');
%     end
%     dt = 2*tm/(N-1);
%     tt = linspace(-tm, tm, N)';
%     g = exp(-tt.^2/2);
%     P = zeros(M + 1, N);
%     P(1, :) = ones(1, N);
%     if M >= 1
%         P(2, :) = 2*tt';
%     end
%     for k = 3:M + 1
%         P(k, :) = 2*tt' .* P(k - 1, :) - 2*(k - 2)*P(k - 2, :);
%     end
%     for k = 1:M + 1
%         Htemp(k, :) = P(k, :).*g' / sqrt(sqrt(pi)*2^(k - 1)*gamma(k)) * sqrt(dt);
%     end
%     h = Htemp(1:M, :)'; 
%     if M >= 1
%         Dh = zeros(N, M); 
%         for k = 1:M
%             Dh(:, k) = (tt .* Htemp(k, :).' - sqrt(2*k)*Htemp(k + 1, :).') * dt;
%         end
%     else
%         Dh = [];
%     end
% end

function [h, Dh] = hermf(N, M, tm)
    dt = 2*tm/(N-1) ; 
    tt = linspace(-tm,tm,N) ; 
    P = [] ; Htemp = [] ; DH = [] ; 
    g = exp(-tt.^2/2) ; 
    P(1,:) = ones(1,N) ; 
    P(2,:) = 2*tt ; 

    for k = 3 : M+1 
        P(k,:) = 2*tt.*P(k-1,:) - 2*(k-2)*P(k-2,:) ; 
    end 

    for k = 1:M+1 
        Htemp(k,:) = P(k,:).*g/sqrt(sqrt(pi)*2^(k-1)*gamma(k))*sqrt(dt) ; 
    end 

    h = Htemp(1:M,:) ; 

    for k = 1:M 
        Dh(k,:) = (tt.*Htemp(k,:) - sqrt(2*(k))*Htemp(k+1,:))*dt ; 
    end
    
end    
