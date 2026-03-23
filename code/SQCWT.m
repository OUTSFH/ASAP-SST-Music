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
c = 3e8;              % 光速（m/s）
wavelength = c / f;   % 波长
num_elements = 4;     % 天线阵元数量
d = wavelength * 0.5; % 阵元间距（波长的0.5倍）
theta = -90:0.1:90;   % 扫描角度
true_angle = -45;     % 真实来波方向
SNR = 10;             % 信噪比
fs = 30e6;            % 采样率 30 MHz
sqcwt_freqlow  = -15e6;
sqcwt_freqhigh = 15e6;
sqcwt_alpha    = 0.1e8; % 频率分辨率
opts = [];              % 若有额外参数，可在 opts 中设置

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
    % 幅度误差：服从正态分布（单位 dB），转换为线性刻度
    amp_error = 10.^(randn(1, M)*amp_sigma_db/20);
    % 相位误差：服从正态分布（单位度），转换为弧度
    phase_error = deg2rad(randn(1, M)*phase_sigma_deg);
    % 加入误差
    A(k,:) = amp_error .* exp(1j*phase_error) .* ideal_steering;
end
    A = A';

% 调整 qpsk_signal 维度以匹配矩阵乘法，并生成多通道信号矩阵
qpsk_signal = reshape(qpsk_signal, 1, []);
data_matrix1 = A * qpsk_signal;  

% 添加高斯白噪声
data_matrix1 = awgn(data_matrix1, SNR, 'measured');

% 对每个通道添加稀疏干扰
num_channels = size(data_matrix1, 1);
for ch = 1:num_channels
    signal = data_matrix1(ch, :).';  
    n_sig = length(signal);
    alpha_int = 0.1;       
    c_val = 1;         
    temp = rand(1, n_sig);
    IND = find(temp < alpha_int);  
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
num_channels = size(data_matrix1, 1);

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';
    n = length(signal);
    IND = find(rand(n,1) < pulse_ratio);
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

% 未滤波 sqCWT 时频分析
S_matrix_unfilt = cell(1, num_channels);  % 存储未滤波各通道 sqCWT 谱
F_matrix_unfilt = cell(1, num_channels);  % 存储未滤波各通道频率轴
T_matrix_unfilt = cell(1, num_channels);  % 存储未滤波各通道时间轴

for ch = 1:num_channels
    signal = data_matrix1_unfiltered(ch, :).';
    n = length(signal);
    t = (0:n-1)/fs;  
    [tfr, tfrsq, tfrtic, tfrsqtic] = sqCWT(t, signal, sqcwt_freqlow, sqcwt_freqhigh, sqcwt_alpha, opts);
    S_matrix_unfilt{ch} = tfrsq;
    F_matrix_unfilt{ch} = tfrsqtic/1e6; 
    T_matrix_unfilt{ch} = t;
end

% ASAP_Hankel_1D 滤波处理
r_filter = 4;       % 目标低秩参数
gamma_filter = 0.5; % 收敛速度参数（介于0和1之间）

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';  
    [filtered_signal, s_hat] = ASAP_Hankel_1D(signal, r_filter, gamma_filter);
    data_matrix1(ch, 1:length(filtered_signal)) = filtered_signal.';
end

% 滤波后 sqCWT 时频分析
S_matrix = cell(1, num_channels);    % 存储滤波后各通道 sqCWT 谱
F_matrix = cell(1, num_channels);    % 存储滤波后各通道频率轴
T_matrix = cell(1, num_channels);    % 存储滤波后各通道时间轴

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';
    n = length(signal);
    t = (0:n-1)/fs;
    [tfr, tfrsq, tfrtic, tfrsqtic] = sqCWT(t, signal, sqcwt_freqlow, sqcwt_freqhigh, sqcwt_alpha, opts);
    S_matrix{ch} = tfrsq;
    F_matrix{ch} = tfrsqtic/1e6;  
    T_matrix{ch} = t;
end

% % 对比未滤波与滤波后的 sqCWT 时频图
% for ch = 1:num_channels
%     figure;
%     subplot(1,2,1);
%     surf(T_matrix_unfilt{ch}, F_matrix_unfilt{ch}, 10*log10(abs(S_matrix_unfilt{ch}).^2), 'EdgeColor', 'none');
%     title(['(未滤波) Channel ', num2str(ch)], 'FontSize', 14);
%     xlabel('Time (s)', 'FontSize', 12);
%     ylabel('Frequency (MHz)', 'FontSize', 12);
%     zlabel('Energy (dB)', 'FontSize', 12);
%     view(2); colorbar; colormap jet; grid on;
%     caxis([-30 40]);
% 
%     subplot(1,2,2);
%     surf(T_matrix{ch}, F_matrix{ch}, 10*log10(abs(S_matrix{ch}).^2), 'EdgeColor', 'none');
%     title(['(滤波后) Channel ', num2str(ch)], 'FontSize', 14);
%     xlabel('Time (s)', 'FontSize', 12);
%     ylabel('Frequency (MHz)', 'FontSize', 12);
%     zlabel('Energy (dB)', 'FontSize', 12);
%     view(2); colorbar; colormap jet; grid on;
%     set(gcf, 'Position', [100,100,1200,500]);
% end

% MUSIC DOA 估计
stft_data = [];
for ch = 1:num_channels
    stft_data = [stft_data, reshape(S_matrix{ch}, [], 1)];
end

% 计算协方差矩阵及特征值分解
R = cov(stft_data);
[eigenVec, eigenVal] = eig(R);
[eigenVal, idx] = sort(diag(eigenVal), 'descend');
eigenVec = eigenVec(:, idx);

% 分离信号子空间与噪声子空间
num_signals = 1;
signalVec = eigenVec(:, 1:num_signals);
noiseVec = eigenVec(:, num_signals+1:end);

% 计算 MUSIC 谱
P_tf_music = zeros(1, length(theta));  
for ii = 1:length(theta)
    a = exp(-1j * (0:num_elements-1)' * (2*pi/wavelength) * d * sin(deg2rad(theta(ii))));
    noiseMatrix = noiseVec * noiseVec';
    P_tf_music(ii) = 1 / (a' * noiseMatrix * a);
end

P_tf_music = 10*log10(abs(P_tf_music)/max(abs(P_tf_music)));
[~, max_idx] = max(P_tf_music);
estimated_angle = theta(max_idx);
disp(['估计的到达角度：', num2str(estimated_angle), '°']);

figure;
plot(theta, P_tf_music, 'LineWidth', 2, 'Color', 'b'); grid on;
title('NEW-sqCWT-MUSIC Spectrum Estimation', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Angle (degrees)', 'FontSize', 16);
ylabel('MUSIC Spectrum (dB)', 'FontSize', 16);
xticks(-90:5:90); yticks(-20:5:20);
legend('NEW-sqCWT-MUSIC Spectrum', 'Location', 'Best', 'FontSize', 12);
