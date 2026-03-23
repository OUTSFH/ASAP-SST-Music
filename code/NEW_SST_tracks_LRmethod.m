clear all
close all

% 添加 PROPACK 工具箱
if exist('.\PROPACK', 'dir')==7
    addpath PROPACK;
else
    fprintf('No PROPACK installed.\n');
    error('Break, PROPACK not installed');
end

% 读取二进制文件，构造 QPSK 信号
fid = fopen('QPSK-1090-2048.bin', 'rb'); 
data = fread(fid, 'float32'); 
fclose(fid);
qpsk_signal = data(1:2:end) + 1j * data(2:2:end);

% 参数设置
f = 1090e6;           % 信号中心频率 1090 MHz
c = 3e8;              % 光速
wavelength = c / f;   % 波长
num_elements = 4;     % 天线阵元数量
d = wavelength * 0.5; % 阵元间距
theta = -90:0.1:90;   % 扫描角度
true_angle = -45;     % 真实来波方向
SNR = 10;             % 信噪比
fs = 30e6;            % 采样率 30 MHz

% SST 参数设置
time_win = 128;         % 窗口长度
overlap = time_win/2;   % 重叠长度
nfft = 256;             % FFT 点数
ths = 0.1;              % SST 重新分配的阈值参数
hf = 15e6;              % 输出最高频率 15 MHz
lf = -15e6;             % 输出最低频率 -15 MHz

% 添加误差参数
amp_sigma_db = 0.5;   % 幅度误差标准差（dB）
phase_sigma_deg = 2;  % 相位误差标准差（度）

% 构建转向矩阵，并添加幅度和相位误差
P = 1;                % 信号源数量
M = num_elements;     % 阵元数量
doa = true_angle;     % 信号源来波方向
lambda = wavelength;  
A = zeros(P, M);

for k = 1:P
    % 理想转向矢量（角度转换为弧度）
    ideal_steering = exp(-1j * 2*pi*d*sin(deg2rad(doa))*(0:M-1)/lambda);
    % 幅度误差（单位 dB）转换为线性刻度
    amp_error = 10.^(randn(1, M)*amp_sigma_db/20);
    % 相位误差（单位度）转换为弧度
    phase_error = deg2rad(randn(1, M)*phase_sigma_deg);
    % 加入误差
    A(k,:) = amp_error .* exp(1j*phase_error) .* ideal_steering;
end
    A = A';

% 构造多通道信号矩阵
qpsk_signal = reshape(qpsk_signal, 1, []);
data_matrix1 = A * qpsk_signal;  

% 添加高斯白噪声
data_matrix1 = awgn(data_matrix1, SNR, 'measured');

% 对各通道添加稀疏干扰
num_channels = size(data_matrix1, 1);
for ch = 1:num_channels
    signal = data_matrix1(ch, :).';  
    n_sig = length(signal);
    alpha = 0.1;       
    c_val = 1;         
    temp = rand(1, n_sig);
    IND = find(temp < alpha);  
    m = length(IND);
    os = zeros(m,1);
    a_val = c_val * mean(abs(real(signal)));
    b_val = c_val * mean(abs(imag(signal)));
    for i = 1:m
        v1 = a_val - 2*a_val*(1 - rand());
        v2 = b_val - 2*b_val*(1 - rand());
        os(i) = v1 + 1i*v2;
    end
    signal(IND) = signal(IND) + os;
    data_matrix1(ch, :) = signal.';
end

% 添加脉冲噪声干扰
pulse_ratio = 0.2;    % 脉冲噪声比例
alfa = 0.75;          % 稳定分布参数 alpha
gamma_param = 1e-2;   % 稳定分布参数 gamma

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';
    n_sig = length(signal);
    IND = find(rand(n_sig,1) < pulse_ratio);
    if ~isempty(IND)
        pulse_noise = zeros(length(IND),1);
        for idx = 1:length(IND)
            r1 = stblrnd(alfa, 0, gamma_param, 0, 1, 1);
            r2 = stblrnd(alfa, 0, gamma_param, 0, 1, 1);
            pulse_noise(idx) = sqrt(1/2) * (r1 + 1i*r2);
        end
        signal(IND) = signal(IND) + pulse_noise;
    end
    data_matrix1(ch, :) = signal.';
end

% for ch = 1:num_channels
%     signal = data_matrix1(ch, :).';
%     n = length(signal);
% 
%     % 计算原始信号功率
%     signal_power = mean(abs(signal).^2);
%    
%     % 根据设定的SNR计算脉冲噪声功率
%     noise_power = signal_power / (10^(SNR / 10));
%     
%     % 生成脉冲噪声（满足稳定分布）
%     r1 = stblrnd(alfa, 0, gamma_param, 0, n, 1);
%     r2 = stblrnd(alfa, 0, gamma_param, 0, n, 1);
%     pulse_noise = sqrt(noise_power / 2) * (r1 + 1i * r2);
% 
%     % 直接将脉冲噪声添加到信号中
%     signal = signal + pulse_noise;
%     
%     % 更新信号矩阵
%     data_matrix1(ch, :) = signal.';
% end

% 保存未滤波数据（含高斯白噪声+稀疏干扰+脉冲噪声）
data_matrix1_unfiltered = data_matrix1;

% 未滤波 SST 分析 
S_matrix_unfilt = cell(1, num_channels);   % 各通道 SST 谱
F_matrix_unfilt = cell(1, num_channels);   % 各通道频率轴
T_matrix_unfilt = cell(1, num_channels);   % 各通道时间轴

for ch = 1:num_channels
    signal = data_matrix1_unfiltered(ch, :).';
    [sst, ~, frequency, t_vector] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
    S_matrix_unfilt{ch} = sst;
    F_matrix_unfilt{ch} = frequency / 1e6; 
    T_matrix_unfilt{ch} = t_vector; 
end

% ASAP_Hankel_1D 滤波处理
r_filter = 4;       % 目标低秩参数
gamma_filter = 0.5; % 收敛速度参数

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';  
    [filtered_signal, ~] = ASAP_Hankel_1D(signal, r_filter, gamma_filter);
    data_matrix1(ch, 1:length(filtered_signal)) = filtered_signal.';
end

% 滤波后 SST 分析 
S_matrix = cell(1, num_channels);   % 各通道滤波后 SST 谱
F_matrix = cell(1, num_channels);   % 各通道频率轴
T_matrix = cell(1, num_channels);   % 各通道时间轴

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';
    [sst, ~, frequency, t_vector] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
    S_matrix{ch} = sst;
    F_matrix{ch} = frequency / 1e6; 
    T_matrix{ch} = t_vector;
end

% 对比未滤波与滤波后的 SST 图
for ch = 1:num_channels
    figure;
    subplot(1,2,1);
    surf(T_matrix_unfilt{ch}, F_matrix_unfilt{ch}, 10*log10(abs(S_matrix_unfilt{ch}).^2), 'EdgeColor', 'none');
    title(['(未滤波) Channel ', num2str(ch)], 'FontSize', 14);
    xlabel('Time (s)');
    ylabel('Frequency (MHz)');
    view(2); colorbar; colormap jet; grid on; caxis([-30 40]);
    
    subplot(1,2,2);
    surf(T_matrix{ch}, F_matrix{ch}, 10*log10(abs(S_matrix{ch}).^2), 'EdgeColor', 'none');
    title(['(滤波后) Channel ', num2str(ch)], 'FontSize', 14);
    xlabel('Time (s)');
    ylabel('Frequency (MHz)');
    view(2); colorbar; colormap jet; grid on;
    set(gcf, 'Position', [100,100,1200,500]);
end

% 多通道轨迹融合
delta_limit = 4;       % 轨迹搜索邻域大小
min_length = 15;       % 轨迹最小长度
LOWER_PRCTILE_LIMIT = 95; % 能量阈值

% 对每个通道均独立提取轨迹，并融合各通道轨迹点
all_time_freq_points = [];
for ch = 1:num_channels
    sst_ref = abs(S_matrix{ch}).';
    [individual_tracks, ~] = tracks_LRmethod(sst_ref, fs, delta_limit, min_length, LOWER_PRCTILE_LIMIT);
    for i = 1:length(individual_tracks)
        track = individual_tracks{i}; 
        all_time_freq_points = [all_time_freq_points; track];
    end
end

% 融合后的轨迹点去重处理
time_freq_points = unique(all_time_freq_points, 'rows');

% 从所有通道 SST 数据中提取融合后轨迹点处的复值
stft_data = [];
for ch = 1:num_channels
    sst_ch = S_matrix{ch};  
    points_values = zeros(size(time_freq_points,1), 1);
    for p = 1:size(time_freq_points,1)
        t_idx = time_freq_points(p,1); 
        f_idx = time_freq_points(p,2); 
        if t_idx <= size(sst_ch,2) && f_idx <= size(sst_ch,1)
            points_values(p) = sst_ch(f_idx, t_idx);
        else
            points_values(p) = 0;
        end
    end
    stft_data = [stft_data, points_values];
end

% 去掉任一通道为 0 的轨迹点
valid_rows = all(stft_data ~= 0, 2);
stft_data = stft_data(valid_rows, :);

% 在各通道滤波后 SST 图上叠加显示融合后的轨迹点
for ch = 1:num_channels
    figure;
    sst_plot = 10*log10(abs(S_matrix{ch}).^2);  
    imagesc(T_matrix{ch}, F_matrix{ch}, sst_plot);
    axis xy;
    colormap jet; colorbar;
    title(['通道 ', num2str(ch), ' SST 图与融合轨迹'], 'FontSize', 14);
    xlabel('Time (s)');
    ylabel('Frequency (MHz)');
    hold on;
   
    for p = 1:size(time_freq_points,1)
        t_idx = time_freq_points(p,1);
        f_idx = time_freq_points(p,2);
        if t_idx <= length(T_matrix{ch}) && f_idx <= length(F_matrix{ch})
            t_val = T_matrix{ch}(t_idx);
            f_val = F_matrix{ch}(f_idx);
            plot(t_val, f_val, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
        end
    end
    hold off;
end

% 基于轨迹点快拍构造协方差矩阵进行 MUSIC DOA 估计
R = (stft_data' * stft_data) / size(stft_data, 1);
[eigenVec, eigenVal] = eig(R);
[eigenVal, idx] = sort(diag(eigenVal), 'descend');
eigenVec = eigenVec(:, idx);

num_signals = 1;
signalVec = eigenVec(:, 1:num_signals);
noiseVec = eigenVec(:, num_signals+1:end);

P_tf_music = zeros(1, length(theta));
for ii = 1:length(theta)
    a = exp(-1j*(0:num_elements-1)'*2*pi*d*sin(deg2rad(theta(ii)))/wavelength);
    noiseMatrix = noiseVec * noiseVec';
    P_tf_music(ii) = 1 / (a' * noiseMatrix * a);
end

P_tf_music = 10*log10(abs(P_tf_music)/max(abs(P_tf_music)));
[~, max_idx] = max(P_tf_music);
estimated_angle = theta(max_idx);
disp(['多通道融合后的 MUSIC 估计到达角度：', num2str(estimated_angle), '°']);

% 绘制基于轨迹点 MUSIC 谱
figure;
plot(theta, P_tf_music, 'LineWidth', 2, 'Color', 'b'); grid on;
title('基于多通道融合轨迹的 MUSIC 谱估计', 'FontSize', 16);
xlabel('角度 (°)', 'FontSize', 14);
ylabel('谱强度 (dB)', 'FontSize', 14);

