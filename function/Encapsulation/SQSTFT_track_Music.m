function [estimated_angle, P_tf_music] = SQSTFT_track_Music(data_matrix1, num_elements, d, wavelength, theta, fs, lowFreq, highFreq, alpha_param, tDS, h, Dh, r_filter, gamma_filter, delta_limit, min_length, LOWER_PRCTILE_LIMIT)

num_channels = size(data_matrix1, 1);
M = num_elements;
P = 1; 
c = 3e8;

for ch = 1:num_channels
    signal = data_matrix1(ch,:).';
    [filtered_signal, ~] = ASAP_Hankel_1D(signal, r_filter, gamma_filter);
    data_matrix1(ch, 1:length(filtered_signal)) = filtered_signal.';
end

S_matrix_filt = cell(1, num_channels);
F_matrix_filt = cell(1, num_channels);
T_matrix_filt = cell(1, num_channels);

for ch = 1:num_channels
    signal = data_matrix1(ch,:).';
    [~, tfrtic, tfrsq, tfrsqtic] = sqSTFT(signal, lowFreq, highFreq, alpha_param, tDS, h', Dh');
    S_matrix_filt{ch} = tfrsq;
    F_matrix_filt{ch} = tfrsqtic * fs;
    T_matrix_filt{ch} = (1:length(signal))/fs;
end

delta_limit = 4;      
min_length = 15;      
LOWER_PRCTILE_LIMIT = 95; 
[tf_pts, aligned] = tracks_multiLRmethod(S_matrix_filt,F_matrix_filt,fs,delta_limit,min_length,LOWER_PRCTILE_LIMIT,d,c);

valid_rows = all(aligned~=0,2);
aligned    = aligned(valid_rows,:);
tf_pts     = tf_pts(valid_rows,:);
R  = (aligned' * aligned) / size(aligned, 1); 

% all_time_freq_points = [];
% for ch = 1:num_channels
%     sst_ref = abs(S_matrix_filt{ch}).';
%     [individual_tracks, ~] = tracks_LRmethod(sst_ref, fs, delta_limit, min_length, LOWER_PRCTILE_LIMIT);
%     for i = 1:length(individual_tracks)
%         all_time_freq_points = [all_time_freq_points; individual_tracks{i}];
%     end
% end
% 
% time_freq_points = unique(all_time_freq_points, 'rows');
% stft_data = [];
% 
% for ch = 1:num_channels
%     sst_ch = S_matrix_filt{ch};  
%     points_values = zeros(size(time_freq_points,1), 1);
%     for p = 1:size(time_freq_points,1)
%         t_idx = time_freq_points(p,1); 
%         f_idx = time_freq_points(p,2); 
%         if t_idx <= size(sst_ch,2) && f_idx <= size(sst_ch,1)
%             points_values(p) = sst_ch(f_idx, t_idx);
%         else
%             points_values(p) = 0;
%         end
%     end
%     stft_data = [stft_data, points_values];
% end
% 
% valid_rows = all(stft_data ~= 0, 2);
% stft_data = stft_data(valid_rows, :);
% R = (stft_data' * stft_data) / size(stft_data, 1);

[eigenVec, eigenVal] = eig(R);
[eigenVal, idx] = sort(diag(eigenVal), 'descend');
eigenVec = eigenVec(:, idx);
num_signals = 1;
noiseVec = eigenVec(:, num_signals+1:end);
P_tf_music = zeros(1, length(theta));

for ii = 1:length(theta)
    a = exp(-1j*(0:M-1)'*(2*pi/wavelength)*d*sin(deg2rad(theta(ii))));
    noiseMatrix = noiseVec * noiseVec';
    P_tf_music(ii) = 1 / (a' * noiseMatrix * a);
end

P_tf_music = 10*log10(abs(P_tf_music)/max(abs(P_tf_music)));
[~, max_idx] = max(P_tf_music);
estimated_angle = theta(max_idx);

end
