function [estimated_angle, P_tf_music] = SQCWT_Music(data_matrix, fs, wavelength, d, theta, lf, hf, sqcwt_alpha, opts, P, r_filter, gamma_filter)

[M, ~] = size(data_matrix);
S_matrix_unfilt = cell(1, M);  
F_matrix_unfilt = cell(1, M);  
T_matrix_unfilt = cell(1, M);  

for ch = 1:M
    signal = data_matrix(ch, :).';
    n = length(signal);
    t = (0:n-1)/fs;  
    [~, tfrsq, ~, tfrsqtic] = sqCWT(t, signal, lf, hf, sqcwt_alpha, opts);
    S_matrix_unfilt{ch} = tfrsq;
    F_matrix_unfilt{ch} = tfrsqtic/1e6; 
    T_matrix_unfilt{ch} = t;
end

filtered_matrix = zeros(size(data_matrix));
for ch = 1:M
    sig = data_matrix(ch, :).';
    [sig_filt, ~] = ASAP_Hankel_1D(sig, r_filter, gamma_filter);
    Lf = length(sig_filt);
    filtered_matrix(ch,1:Lf) = sig_filt.';
end

S_matrix = cell(1, M);
for ch = 1:M
    sig = filtered_matrix(ch, :).';
    n = length(sig);
    t = (0:n-1)/fs;
    [~, tfrsq, ~, ~] = sqCWT(t, sig, lf, hf, sqcwt_alpha, opts);
    S_matrix{ch} = tfrsq;
end

stft_data = [];
for ch = 1:M
    stft_data = [stft_data, reshape(S_matrix{ch}, [], 1)];
end

R = cov(stft_data);
[eigVec, eigVal] = eig(R);
eigVal = diag(eigVal); [~, idx] = sort(eigVal, 'descend');
noiseVec = eigVec(:, idx(P+1:end));
P_tf_music = zeros(1, length(theta));
theta_rad = deg2rad(theta);

for ii = 1:length(theta)
    a = exp(1j * (0:M-1)' * (2*pi/wavelength) * d * sin(theta_rad(ii)));
    P_tf_music(ii) = 1 / (a' * (noiseVec * noiseVec') * a);
end

P_tf_music = 10*log10(abs(P_tf_music) / max(abs(P_tf_music)));
[~, max_idx] = max(P_tf_music);
estimated_angle = theta(max_idx);

end
