% function [estimated_angle, P_tf_music] = SST_track_Music(data, num_elements, d, wavelength, theta, fs, time_win, overlap, nfft, ths, hf, lf, r_filter, gamma_filter, delta_limit, min_length, LOWER_PRCTILE)
% % SST_track_Music 函数：基于平滑稀疏变换（SST）轨迹融合的MUSIC DOA估计算法
% % 输入参数：
% %   data: 多通道信号矩阵（通道×时间，每行一个通道，每列一个采样点）
% %   num_elements: 天线阵元数量
% %   d: 阵元间距（米）
% %   wavelength: 信号波长（米）
% %   theta: 扫描角度范围（度数）
% %   fs: 采样率（Hz）
% %   time_win: SST窗口长度（样本数）
% %   overlap: SST重叠长度（样本数）
% %   nfft: SST的FFT点数
% %   ths: SST重新分配的阈值参数
% %   hf: 输出最高频率（Hz）
% %   lf: 输出最低频率（Hz）
% %   r_filter: ASAP_Hankel_1D低秩参数
% %   gamma_filter: ASAP_Hankel_1D收敛速度参数
% %   delta_limit: 轨迹搜索邻域大小
% %   min_length: 轨迹最小长度
% %   LOWER_PRCTILE: 能量阈值
% % 输出参数：
% %   estimated_angle: 估计的到达角度（度数）
% %   P_tf_music: SST-track-MUSIC空间谱（dB尺度）
% 
% M = num_elements;
% num_channels = size(data, 1);
% c = 3e8;
% 
% % 未滤波 SST 分析 
% S_matrix_unfilt = cell(1, num_channels);
% for ch = 1:num_channels
%     signal = data(ch, :).';
%     [sst, ~, frequency, t_vector] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
%     S_matrix_unfilt{ch} = sst;
% end
% 
% % ASAP_Hankel_1D 滤波处理
% for ch = 1:num_channels
%     signal = data(ch, :).';  
%     [filtered_signal, ~] = ASAP_Hankel_1D(signal, r_filter, gamma_filter);
%     data(ch, 1:length(filtered_signal)) = filtered_signal.';
% end
% 
% % 滤波后 SST 分析 
% S_matrix = cell(1, num_channels);   % 各通道滤波后 SST 谱
% F_matrix = cell(1, num_channels);   % 各通道频率轴
% T_matrix = cell(1, num_channels);   % 各通道时间轴
% 
% for ch = 1:num_channels
%     signal = data(ch, :).';
%     [sst, ~, frequency, t_vector] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
%     S_matrix{ch} = sst;
%     F_matrix{ch} = frequency / 1e6; 
%     T_matrix{ch} = t_vector;
% end
% 
% delta_limit = 4;      
% min_length = 15;      
% LOWER_PRCTILE_LIMIT = 95; 
% 
% [tf_pts, aligned] = tracks_multiLRmethod(S_matrix,F_matrix,fs,delta_limit,min_length,LOWER_PRCTILE_LIMIT,d,c);
% valid_rows = all(aligned~=0,2);
% aligned    = aligned(valid_rows,:);
% tf_pts     = tf_pts(valid_rows,:);
% R  = (aligned' * aligned) / size(aligned, 1); 
% 
% % % 多通道轨迹融合
% % all_time_freq_points = [];
% % for ch = 1:num_channels
% %     sst_ref = abs(S_matrix{ch}).';
% %     [individual_tracks, ~] = tracks_LRmethod(sst_ref, fs, delta_limit, min_length, LOWER_PRCTILE);
% %     for i = 1:length(individual_tracks)
% %         track = individual_tracks{i}; 
% %         all_time_freq_points = [all_time_freq_points; track];
% %     end
% % end
% % time_freq_points = unique(all_time_freq_points, 'rows');
% % 
% % % 从所有通道 SST 数据中提取融合后轨迹点处的复值
% % stft_data = [];
% % for ch = 1:num_channels
% %     sst_ch = S_matrix{ch};  
% %     points_values = zeros(size(time_freq_points,1), 1);
% %     for p = 1:size(time_freq_points,1)
% %         t_idx = time_freq_points(p,1); 
% %         f_idx = time_freq_points(p,2); 
% %         if t_idx <= size(sst_ch,2) && f_idx <= size(sst_ch,1)
% %             points_values(p) = sst_ch(f_idx, t_idx);
% %         else
% %             points_values(p) = 0;
% %         end
% %     end
% %     stft_data = [stft_data, points_values];
% % end
% % valid_rows = all(stft_data ~= 0, 2);
% % stft_data = stft_data(valid_rows, :);
% % 
% % % 基于轨迹点快拍构造协方差矩阵进行 MUSIC DOA 估计
% % R = (stft_data' * stft_data) / size(stft_data, 1);
% 
% [eigenVec, eigenVal] = eig(R);
% [eigenVal, idx] = sort(diag(eigenVal), 'descend');
% eigenVec = eigenVec(:, idx);
% num_signals = 1;
% noiseVec = eigenVec(:, num_signals+1:end);
% P_tf_music = zeros(1, length(theta));
% 
% for ii = 1:length(theta)
%     a = exp(-1j*(0:M-1)'*2*pi*d*sin(deg2rad(theta(ii)))/wavelength);
%     noiseMatrix = noiseVec * noiseVec';
%     P_tf_music(ii) = 1 / (a' * noiseMatrix * a);
% end
% 
% P_tf_music = 10*log10(abs(P_tf_music)/max(abs(P_tf_music)));
% [~, max_idx] = max(P_tf_music);
% estimated_angle = theta(max_idx);
% 
% end
function [estimated_angles, P_tf_music] = SST_track_Music(data, num_elements, d, wavelength, theta, fs, time_win, overlap, nfft, ths, hf, lf, r_filter, gamma_filter, delta_limit, min_length, LOWER_PRCTILE)

M = num_elements;
num_channels = size(data, 1);
c = 3e8;

% === Step 1: 原始 SST 分析 ===
S_matrix_unfilt = cell(1, num_channels);
for ch = 1:num_channels
    signal = data(ch, :).';
    [sst, ~, frequency, t_vector] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
    S_matrix_unfilt{ch} = sst;
end

% === Step 2: ASAP_Hankel_1D 滤波处理 ===
for ch = 1:num_channels
    signal = data(ch, :).';  
    [filtered_signal, ~] = ASAP_Hankel_1D(signal, r_filter, gamma_filter);
    data(ch, 1:length(filtered_signal)) = filtered_signal.';
end

% === Step 3: 滤波后 SST 分析 ===
S_matrix = cell(1, num_channels);   
F_matrix = cell(1, num_channels);   
T_matrix = cell(1, num_channels);   

for ch = 1:num_channels
    signal = data(ch, :).';
    [sst, ~, frequency, t_vector] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
    S_matrix{ch} = sst;
    F_matrix{ch} = frequency / 1e6; 
    T_matrix{ch} = t_vector;
end

% === Step 4: 多通道轨迹提取与融合 ===
[tf_pts, aligned] = tracks_multiLRmethod(S_matrix, F_matrix, fs, delta_limit, min_length, LOWER_PRCTILE, d, c);
valid_rows = all(aligned~=0, 2);
aligned = aligned(valid_rows, :);
tf_pts = tf_pts(valid_rows, :);

% === Step 5: 协方差矩阵与特征分解 ===
R  = (aligned' * aligned) / size(aligned, 1);
[eigenVec, eigenVal] = eig(R);
[eigenVal, idx] = sort(diag(eigenVal), 'descend');
eigenVec = eigenVec(:, idx);

num_signals = 2;  % 由1改为2
noiseVec = eigenVec(:, num_signals+1:end);

% === Step 6: MUSIC 空间谱 ===
P_tf_music = zeros(1, length(theta));
for ii = 1:length(theta)
    a = exp(1j * (0:M-1)' * 2*pi*d * sin(deg2rad(theta(ii))) / wavelength);
    P_tf_music(ii) = 1 / (a' * (noiseVec * noiseVec') * a);
end
P_tf_music = 10 * log10(abs(P_tf_music) / max(abs(P_tf_music)));

% === Step 7: 估计两个角度 ===
[~, sorted_idx] = sort(P_tf_music, 'descend');
estimated_angles = sort(theta(sorted_idx(1:num_signals)));  % 输出为升序

end
