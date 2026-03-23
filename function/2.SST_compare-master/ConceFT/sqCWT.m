function [tfr, tfrsq, tfrtic, tfrsqtic] = sqCWT(t, x, freqlow, freqhigh, alpha, opts)
oct = 1;         
scale = 2;        
nvoice = 30;       
Gamma = 1e-8;     
noctave = floor(log2(length(x))) - oct;
dt = t(2) - t(1);
[tfr, tfrtic] = CWT(t, x, oct, scale, nvoice, opts);
Dtfr = (-1i/2/pi/dt)*[tfr(2:end,:) - tfr(1:end-1,:); tfr(end,:)-tfr(end-1,:)] ;
Dtfr((abs(tfr) < Gamma)) = NaN;
omega = Dtfr./tfr;
[tfrsq, tfrsqtic] = SQ(tfr, omega, alpha, scale, nvoice, freqlow, freqhigh);
tfr = transpose(tfr) ;
tfrsq = transpose(tfrsq) ;    
end

function [tfrsq, freq] = SQ(tfd, omega, alpha, scale, nvoice, freqlow, freqhigh)
    omega = abs(omega);  
    [n, nscale] = size(tfd);
    nalpha = floor((freqhigh - freqlow) / alpha) + 1;
    tfrsq = zeros(n, nalpha);
    freq = linspace(freqlow, freqhigh, nalpha);
    nfreq = length(freq);
    
    for b = 1:n            
        for kscale = 1:nscale     
            qscale = scale * (2^(kscale/nvoice));
            if (isfinite(omega(b, kscale)) && (omega(b, kscale) > 0))
                k = floor((omega(b, kscale) - freqlow) / alpha) + 1;               
                if (isfinite(k) && (k < nfreq))
                    ha = freq(k+1) - freq(k);
                    tfrsq(b,k) = tfrsq(b,k) + log(2) * tfd(b,kscale) * sqrt(qscale) / ha / nvoice;
                end
            end
        end
    end
end
