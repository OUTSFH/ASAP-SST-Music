% function [estimated_angle, P_music] = DESHAPE_Music(data_matrix, num_elements, d, wavelength, theta, fs, time_win, overlap, nfft, ths, hf, lf,r_filter,gamma_filter)
%     
%     num_channels = size(data_matrix, 1);
%     S_cell = cell(1, num_channels);
% 
%     for ch = 1:num_channels
%         x = data_matrix(ch, :).';
%         [deshape, ~, ~, ~, ~, ~] = deShape_2(x, fs, time_win, hf, ths, overlap, nfft, lf, ths);
%         S_cell{ch} = deshape;
%     end
% 
% % for ch = 1:num_channels
% %     signal = data_matrix(ch, :).';  
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
%     stft_data = [];
%     for ch = 1:num_channels
%         stft_data = [stft_data, reshape(S_cell{ch}, [], 1)];
%     end
% 
%     R = cov(stft_data);
%     [eigVec, eigVal] = eig(R);
%     [~, idx] = sort(diag(eigVal), 'descend');
%     noiseVec = eigVec(:, idx(2:end));    
%     P_music = zeros(1, length(theta));
% 
%     for ii = 1:length(theta)
%         a = exp(-1j * (0:num_elements-1)' * (2*pi/wavelength) * d * sin(deg2rad(theta(ii))));
%         P_music(ii) = 1 / (a' * (noiseVec * noiseVec') * a);
%     end
% 
%     P_music = 10*log10(abs(P_music)/max(abs(P_music)));
%     [~, max_idx] = max(P_music);
%     estimated_angle = theta(max_idx);
% 
% end

function [estimated_angles, P_music] = DESHAPE_Music(data_matrix, num_elements, d, wavelength, theta, fs, time_win, overlap, nfft, ths, hf, lf, r_filter, gamma_filter)
    
    num_channels = size(data_matrix, 1);
    S_cell = cell(1, num_channels);

    % deShape 变换
    for ch = 1:num_channels
        x = data_matrix(ch, :).';
        [deshape, ~, ~, ~, ~, ~] = deShape_2(x, fs, time_win, hf, ths, overlap, nfft, lf, ths);
        S_cell{ch} = deshape;
    end

    % 将每个通道的时频结果拼接成矩阵（每列一个通道）
    stft_data = [];
    for ch = 1:num_channels
        stft_data = [stft_data, reshape(S_cell{ch}, [], 1)];
    end

    % 协方差矩阵及特征分解
    R = cov(stft_data);
    [eigVec, eigVal] = eig(R);
    [~, idx] = sort(diag(eigVal), 'descend');

    % 对于两个信号，噪声子空间取从第三大特征向量开始
    noiseVec = eigVec(:, idx(3:end));

    % 计算 MUSIC 谱
    P_music = zeros(1, length(theta));
    for ii = 1:length(theta)
        a = exp(1j * (0:num_elements-1)' * (2*pi/wavelength) * d * sin(deg2rad(theta(ii))));
        P_music(ii) = 1 / (a' * (noiseVec * noiseVec') * a);
    end

    P_music = 10*log10(abs(P_music) / max(abs(P_music)));
    [~, sorted_idx] = sort(P_music, 'descend');
    top2 = sorted_idx(1:2);
    estimated_angles = theta(top2);

end
