function [S, F, T] = cus_stft(x, fs, window, overlap, nfft)
 
    N = length(x);
    win_len = length(window);
    hop     = win_len - overlap;

    if N < win_len
        x    = [x; zeros(win_len - N, 1)];
        N    = win_len;
    end

    num_frames = floor((N - win_len)/hop) + 1;

    if mod(nfft,2)~=0
        error('nfft 必须为偶数');
    end

    S = zeros(nfft, num_frames);
    T = (0:num_frames-1) * hop / fs;
    F = (-nfft/2 : nfft/2-1) * (fs / nfft);

    for i = 1:num_frames
        idx1        = (i-1)*hop + 1;
        idx2        = idx1 + win_len - 1;
        frame       = x(idx1:idx2) .* window;
        fft_result  = fft(frame, nfft);
        S(:, i)     = fftshift(fft_result);
    end
end
