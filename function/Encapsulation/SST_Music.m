% function [estimated_angle, P_tf_music] = SST_Music(data, fs, wavelength, d, theta, time_win, overlap, nfft, ths, hf, lf, r_filter, gamma_filter)
% 
% num_channels = size(data, 1);
% M = num_channels; 
% S_matrix_unfilt = cell(1, num_channels);    
% 
% for ch = 1:num_channels
%     signal = data(ch, :).';
%     [sst, ~, ~, ~] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
%     S_matrix_unfilt{ch} = sst;
% end
% 
% % for ch = 1:num_channels
% %     signal = data(ch, :).';  
% %     [filtered_signal, ~] = ASAP_Hankel_1D(signal, r_filter, gamma_filter);
% %     data(ch, 1:length(filtered_signal)) = filtered_signal.';
% % end 
% %  
% % S_matrix = cell(1, num_channels);   
% % for ch = 1:num_channels
% %     signal = data(ch, :).';
% %     [sst, ~, ~, ~] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
% %     S_matrix{ch} = sst;
% % end
% 
% stft_data = [];
% for ch = 1:num_channels
%     stft_data = [stft_data, reshape(S_matrix_unfilt{ch}, [], 1)];
% end
% 
% R = cov(stft_data);
% [eigenVec, eigenVal] = eig(R);
% [eigenVal, idx] = sort(diag(eigenVal), 'descend');
% eigenVec = eigenVec(:, idx);
% noiseVec = eigenVec(:, 2:end); 
% theta_rad = deg2rad(theta);
% P_tf_music = zeros(1, length(theta));
% 
% for ii = 1:length(theta)
%     a = exp(-1j * (0:M-1)' * (2*pi/wavelength) * d * sin(theta_rad(ii)));
%     noiseMatrix = noiseVec * noiseVec';
%     P_tf_music(ii) = 1 / (a' * noiseMatrix * a);
% end
% 
% P_tf_music = 10 * log10(abs(P_tf_music) / max(abs(P_tf_music)));
% [~, max_idx] = max(P_tf_music);
% estimated_angle = theta(max_idx);
% 
% end

function [estimated_angles, P_tf_music] = SST_Music(data, fs, wavelength, d, theta, time_win, overlap, nfft, ths, hf, lf, r_filter, gamma_filter)

    num_channels = size(data, 1);
    M = num_channels; 
    S_matrix_unfilt = cell(1, num_channels);    

    for ch = 1:num_channels
        signal = data(ch, :).';
        [sst, ~, ~, ~] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
        S_matrix_unfilt{ch} = sst;
    end

    % 构建复向量矩阵
    stft_data = [];
    for ch = 1:num_channels
        stft_data = [stft_data, reshape(S_matrix_unfilt{ch}, [], 1)];
    end

    % 协方差矩阵及特征值分解
    R = cov(stft_data);
    [eigenVec, eigenVal] = eig(R);
    [eigenVal, idx] = sort(diag(eigenVal), 'descend');
    eigenVec = eigenVec(:, idx);

    % 保留前两个信号，构造噪声子空间
    num_signals = 2;  % 修改为两个信号
    noiseVec = eigenVec(:, num_signals+1:end); 

    % 计算 MUSIC 谱
    theta_rad = deg2rad(theta);
    P_tf_music = zeros(1, length(theta));

    for ii = 1:length(theta)
        a = exp(1j * (0:M-1)' * (2*pi/wavelength) * d * sin(theta_rad(ii)));
        noiseMatrix = noiseVec * noiseVec';
        P_tf_music(ii) = 1 / (a' * noiseMatrix * a);
    end

    % 归一化谱
    P_tf_music = 10 * log10(abs(P_tf_music) / max(abs(P_tf_music)));

    % 找到谱峰对应的两个角度
    [~, sorted_idx] = sort(P_tf_music, 'descend');
    estimated_angles = sort(theta(sorted_idx(1:num_signals)));  % 从小到大排序输出

end
