function [estimated_angle, P_tf_music] = STFT_Music(data_matrix1, fs, window, overlap, nfft, num_elements, wavelength, d, theta, P)

num_channels = size(data_matrix1, 1);
S_matrix = cell(1, num_channels);

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';
    [S_matrix{ch}, ~, ~] = cus_stft(signal, fs, window, overlap, nfft);
end

stft_data = [];
for ch = 1:num_channels
    stft_data = [stft_data, reshape(S_matrix{ch}, [], 1)];
end

R = cov(stft_data);
[eigenVec, eigenVal] = eig(R);
[eigenVal, idx] = sort(diag(eigenVal), 'descend');
eigenVec = eigenVec(:, idx);
noiseVec = eigenVec(:, P+1:end);
P_tf_music = zeros(1, length(theta));

for ii = 1:length(theta)
    a = exp(-1j * (0:num_elements-1)' * (2 * pi / wavelength) * d * sin(deg2rad(theta(ii))));
    noiseMatrix = noiseVec * noiseVec';
    P_tf_music(ii) = 1 / (a' * noiseMatrix * a);
end

P_tf_music = 10 * log10(abs(P_tf_music) / max(abs(P_tf_music)));
[~, max_idx] = max(P_tf_music);
estimated_angle = theta(max_idx);

end
