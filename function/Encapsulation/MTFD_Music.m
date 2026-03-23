function [estimated_angle, P_theta] = MTFD_Music(S, theta, num_elements, d, wavelength, P, num_segments)
    
    total_length = size(S,2);
    seg_length = floor(total_length / num_segments);
    total_length_new = num_segments * seg_length;
    S = S(:, 1:total_length_new);
    M = num_elements;
    D_accum = cell(M, M);
    
    for i = 1:M
        for j = 1:M
            D_accum{i,j} = [];
        end
    end

    for seg = 1:num_segments
        idx_start = (seg-1)*seg_length + 1;
        idx_end = seg * seg_length;  
        S_seg = S(:, idx_start:idx_end);  

        for i = 1:M
            for j = 1:M
                s1 = S_seg(i,:);
                s2 = S_seg(j,:);
                N_seg = length(s1);
                if mod(N_seg,2)==0
                    L_seg = N_seg - 1;
                else
                    L_seg = N_seg;
                end
                
                MM = 2^nextpow2(N_seg);
                K = 2 * 2^nextpow2(N_seg);
 
                z1 = fft(s1, K);
                z1(2:K/2) = 2 * z1(2:K/2);
                z1(K/2+1:K) = 0;
                z1 = ifft(z1);
                
                z2 = fft(s2, K);
                z2(2:K/2) = 2 * z2(2:K/2);
                z2(K/2+1:K) = 0;
                z2 = ifft(z2);

                K_TL = zeros(K, N_seg);
                L_half = fix(L_seg/2);
                win = 'blackmanharris';
                g = feval(win, L_seg);

                for n = 1:K
                    for tau = -L_half:L_half
                        G = g(1 + tau + L_half);
                        idx1 = 1 + mod(n-1+tau, K);
                        idx2 = 1 + mod(n-1-tau, K);
                        mm = 1 + mod(tau, MM);
                        K_TL(mm, n) = G * z1(idx1) * conj(z2(idx2));
                    end
                end
                
                TFR = zeros(MM, N_seg);
                temp = fft(K_TL, MM);
                TFR(:, 2:N_seg) = temp(:, 1:N_seg-1);

                if isempty(D_accum{i,j})
                    D_accum{i,j} = TFR(:, 2:N_seg);
                else
                    D_accum{i,j} = D_accum{i,j} + TFR(:, 2:N_seg);
                end
            end
        end
    end

    for i = 1:M
        for j = 1:M
            D_accum{i,j} = D_accum{i,j} / num_segments;
        end
    end

    D_avg = zeros(size(D_accum{1,1}));
    for ii = 1:M
        D_avg = D_avg + D_accum{ii,ii};
    end

    D_avg = D_avg / M;
    thr = 0.9 * max(D_avg(:));
    Tr = abs(D_avg) >= thr;
    [~, n_p] = find(Tr);
    n_p = length(n_p);
    D_s = zeros(M, M);
    
    for m1 = 1:M
        for m2 = 1:M
            D_s(m1, m2) = sum(sum(D_accum{m1, m2} .* Tr)) / n_p;
        end
    end

    theta_rad = deg2rad(theta);
    [V, ~] = eig(D_s);
    Un = V(:, P+1:end);  
    P_theta = zeros(size(theta));
    
    for ii = 1:length(theta)
        a = exp(-1j * 2*pi*d/wavelength * sin(theta_rad(ii)) * (0:M-1));
        P_theta(ii) = 1 / (a * Un * Un' * a');
    end

    P_theta = 10 * log10(abs(P_theta) / max(abs(P_theta)));
    [~, idx] = max(P_theta);
    estimated_angle = theta(idx);

end
