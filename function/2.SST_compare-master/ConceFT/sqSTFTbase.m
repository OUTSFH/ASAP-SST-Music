function [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, tDS, h, Dh, Smooth, Hemi)
    [xrow, xcol] = size(x);
    t = 1:xrow;
    tLen = length(t(1:tDS:end));
    N_freq = length(-0.5 + alpha:alpha:0.5);
    tfrtic = linspace(-0.5, 0.5, N_freq)'; 
    Lidx = round((lowFreq + 0.5) * N_freq) + 1; 
    Hidx = round((highFreq + 0.5) * N_freq);
    fLen = Hidx - Lidx + 1;
    tfrsqtic = linspace(lowFreq, highFreq, fLen)';
    
    if xcol ~= 1
        error('X must have only one column');
    elseif highFreq > 0.5 || lowFreq < -0.5
        error('TopFreq must be a value in [-0.5, 0.5]');
    elseif (tDS < 1) || (rem(tDS, 1))
        error('tDS must be an integer value >= 1');
    end
    
    [hrow, hcol] = size(h);
    Lh = (hrow - 1)/2;
    if (hcol ~= 1) || (rem(hrow, 2) == 0)
        error('H must be a smoothing window with odd length');
    end
    
    tfr = zeros(N_freq, tLen); 
    tfrsq = zeros(fLen, tLen);
    Ex = mean(abs(x).^2);
    Threshold = 1.0e-8 * Ex;
    Mid = round(length(tfrsqtic)/2);
    Delta = 20 * (tfrsqtic(2) - tfrsqtic(1))^2;
    weight = exp(-(tfrsqtic(max(Mid-10, 1):min(Mid+10, length(tfrsqtic))) - tfrsqtic(Mid)).^2 / Delta);
    weight = weight / sum(weight);
    weightIDX = (max(Mid-10, 1):min(Mid+10, length(tfrsqtic))) - Mid;
    
    for tidx = 1:tLen
        ti = t((tidx-1)*tDS + 1);
        tau = -min([round(N_freq/2)-1, Lh, ti-1]):min([round(N_freq/2)-1, Lh, xrow-ti]);
        indices = mod(N_freq + tau, N_freq) + 1;
        norm_h = norm(h(Lh+1+tau));
        
        tf0 = zeros(N_freq, 1);
        tf1 = zeros(N_freq, 1);
        tf0(indices) = x(ti+tau) .* conj(h(Lh+1+tau)) / norm_h;
        tf1(indices) = x(ti+tau) .* conj(Dh(Lh+1+tau)) / norm_h;
        tf0 = fft(tf0); 
        tf0 = fftshift(tf0); 
        tf0 = tf0(1:N_freq); 
        tf1 = fft(tf1);
        tf1 = fftshift(tf1); 
        tf1 = tf1(1:N_freq);
        
        omega = zeros(size(tf1));
        avoid_warn = find(tf0 ~= 0);
        omega(avoid_warn) = round(imag(N_freq * tf1(avoid_warn) ./ tf0(avoid_warn) / (2*pi)));
        
        sst = zeros(fLen, 1);
        for jcol = 1:N_freq
            if abs(tf0(jcol)) > Threshold
                jcolhat = jcol - omega(jcol);
                if (jcolhat <= Hidx) && (jcolhat >= Lidx)
                    if Smooth
                        IDXb = find((jcolhat - Lidx + 1 + weightIDX <= Hidx) & (jcolhat - Lidx + 1 + weightIDX >= Lidx));
                        IDXa = jcolhat - Lidx + 1 + weightIDX(IDXb);
                        if Hemi
                            if real(tf0(jcol)) > 0
                                sst(IDXa) = sst(IDXa) + tf0(jcol) * weight(IDXb);
                            else
                                sst(IDXa) = sst(IDXa) - tf0(jcol) * weight(IDXb);
                            end
                        else
                            sst(IDXa) = sst(IDXa) + tf0(jcol) * weight(IDXb);
                        end
                    else
                        if Hemi
                            if real(tf0(jcol)) > 0
                                sst(jcolhat - Lidx + 1) = sst(jcolhat - Lidx + 1) + tf0(jcol);
                            else
                                sst(jcolhat - Lidx + 1) = sst(jcolhat - Lidx + 1) - tf0(jcol);
                            end
                        else
                            sst(jcolhat - Lidx + 1) = sst(jcolhat - Lidx + 1) + tf0(jcol);
                        end
                    end
                end
            end
        end
        tfr(:, tidx) = tf0;
        tfrsq(:, tidx) = sst * 2 * (tfrtic(2) - tfrtic(1));
    end
end
