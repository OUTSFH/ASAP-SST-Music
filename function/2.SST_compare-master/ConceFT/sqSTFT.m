function [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, lowFreq, highFreq, alpha, tDS, h, Dh)
    [xrow, ~] = size(x);
    t = 1:length(x);
    tLen = length(t(1:tDS:length(x)));
    N = ceil((highFreq - lowFreq) / alpha) + 1;
    tfrtic = linspace(lowFreq, highFreq, N)';
    
    % 预分配
    tfr = zeros(N, tLen);    % 针对窗函数 h 的结果
    tf3 = zeros(N, tLen);    % 针对窗函数 Dh 的结果
    
    [hrow, ~] = size(h);
    Lidx_h = (hrow-1)/2;
    
    % 计算 STFT
    for tidx = 1:tLen
        ti = t((tidx-1)*tDS+1);
        tau = -min([round(N/2)-1, Lidx_h, ti-1]) : min([round(N/2)-1, Lidx_h, xrow-ti]);
        indices = rem(N + tau, N) + 1;
        norm_h = norm(h(Lidx_h+1+tau));
        tfr(indices, tidx) = x(ti+tau) .* conj(h(Lidx_h+1+tau)) / norm_h;
        tf3(indices, tidx) = x(ti+tau) .* conj(Dh(Lidx_h+1+tau)) / norm_h;
    end
    
    % 对每一列进行 FFT，并对频率轴做 fftshift，使得 0 频出现在中间
    tfr = fft(tfr, N, 1);
    tf3 = fft(tf3, N, 1);
    tfr = fftshift(tfr, 1);
    tf3 = fftshift(tf3, 1);
    
    % 计算重赋值参数 tf3 的值（避免除 0 的情况）
    avoid_warn = find(tfr ~= 0);
    tf3(avoid_warn) = round(imag(N * tf3(avoid_warn) ./ tfr(avoid_warn) / (2*pi)));
    
    % 同步重赋值（synchrosqueezing）
    tfrsq = zeros(N, tLen);
    Ex = mean(abs(x(t(1):t(end))).^2);
    Threshold = 1.0e-8 * Ex;
    
    for icol = 1:tLen
        for jcol = 1:N
            if abs(tfr(jcol,icol)) > Threshold
                jcolhat = jcol - tf3(jcol,icol);
                jcolhat = mod(round(jcolhat-1), N) + 1;
                tfrsq(jcolhat,icol) = tfrsq(jcolhat,icol) + tfr(jcol,icol);
            end
        end
    end
    tfrsqtic = tfrtic;
end
