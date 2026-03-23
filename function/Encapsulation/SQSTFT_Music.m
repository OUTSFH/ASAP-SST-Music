function [estimated_angle, P_tf_music] = SQSTFT_Music(data_matrix, wavelength, d, theta, lowFreq, highFreq, alpha_param, tDS, h, Dh, P)

[M, ~] = size(data_matrix);
S_unfilt = cell(1, M);
for ch = 1:M
    sig = data_matrix(ch, :).';
    [~, ~, tfrsq, ~] = sqSTFT(sig, lowFreq, highFreq, alpha_param, tDS, h, Dh);
    S_unfilt{ch} = tfrsq;
end

filtered = zeros(size(data_matrix));
for ch = 1:M
    sig = data_matrix(ch, :).';
    [sig_filt, ~] = ASAP_Hankel_1D(sig, 4, 0.5);  % hard-coded r_filter=4, gamma_filter=0.5
    Lf = length(sig_filt);
    filtered(ch,1:Lf) = sig_filt.';
end

S_filt = cell(1, M);
for ch = 1:M
    sig = filtered(ch, :).';
    [~, ~, tfrsq, ~] = sqSTFT(sig, lowFreq, highFreq, alpha_param, tDS, h, Dh);
    S_filt{ch} = tfrsq;
end

stft_data = [];
for ch = 1:M
    stft_data = [stft_data, reshape(S_unfilt{ch}, [], 1)];
end

R = cov(stft_data);
[eigVec, eigVal] = eig(R);
eigVal = diag(eigVal);
[~, idx] = sort(eigVal, 'descend');
noiseVec = eigVec(:, idx(P+1:end));
P_tf_music = zeros(1, length(theta));
theta_rad = deg2rad(theta);

for ii = 1:length(theta)
    a = exp(-1j * (0:M-1)' * (2*pi/wavelength) * d * sin(theta_rad(ii)));
    P_tf_music(ii) = 1 / (a' * (noiseVec * noiseVec') * a);
end

P_tf_music = 10*log10(abs(P_tf_music)/max(abs(P_tf_music)));
[~, max_idx] = max(P_tf_music);
estimated_angle = theta(max_idx);

end
