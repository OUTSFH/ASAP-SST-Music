% 参数设置
f = 1090e6;           % 信号中心频率 1090 MHz
c = 3e8;              % 光速（米/秒）
wavelength = c / f;   % 波长
num_elements = 4;     % 天线阵元数量
d = wavelength * 0.5; % 阵元间距为波长的0.5倍
theta = -90:0.1:90;   % 扫描角度
true_angle = -40; % 真实来波方向
SNR = 5;              % 信噪比
fs = 30e6;            % 采样率 30 MHz
time_win = 128;       % 窗口长度
overlap = time_win / 2; % 重叠长度
nfft = 256;           % FFT 点数
beta = 6;             % Kaiser窗的 beta 参数
window = kaiser(time_win, beta); % Kaiser窗
ths = 0.1;              % deShape 重新分配的阈值参数
hf = 15e6;              % 输出最高频率 15 MHz
lf = -15e6;             % 输出最低频率 -15 MHz
sqcwt_alpha    = 0.1e8; % 频率分辨率
opts = [];              % 若有额外参数，可在 opts 中设置
WinLen = time_win;
if mod(WinLen,2)==0
    WinLen = WinLen + 1;  % 确保窗长度为奇数
end
num_segments = 20;

% sqSTFT 时频分析
lowFreq = -0.5;    % 归一化频率下限
highFreq = 0.5;    % 归一化频率上限
alpha_param = 1/nfft; % 频率分辨率（归一化）
tDS = 1;           % 时间下采样因子
[h, Dh] = hermf(71, 1, 6); % 生成 Hermite 窗及其导数

% ASAP_Hankel_1D 滤波处理
r_filter = 4;       
gamma_filter = 0.5; 

% 多通道轨迹融合
delta_limit = 4;       % 轨迹搜索邻域大小
min_length = 15;       % 轨迹最小长度
LOWER_PRCTILE_LIMIT = 95; % 能量阈值

% 生成信号矩阵，每个阵元作为一通道
data_matrix1 = load_data('QPSK2.mat');
data_matrix1 = data_matrix1(1:1024,:)';

% 添加高斯白噪声
data_matrix1 = awgn(data_matrix1, SNR, 'measured');

% 在添加完白噪声的基础上，对每个通道再添加稀疏干扰
num_channels = size(data_matrix1, 1);
for ch = 1:num_channels
    signal = data_matrix1(ch, :).';  
    n = length(signal);
    alpha = 0.1;       
    c_val = 1;         
    temp = rand(1, n);
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

% 脉冲噪声干扰
% alpha特征指数，0 < alpha <= 2，α 值越小，脉冲性越强，脉冲流越大；α 值越大，脉冲性越弱，脉冲流越小。
% beta对称参数，-1 <= beta <= 1
% gamma比例参数，gamma > 0，γ 值越大，脉冲噪声的幅度越大，脉冲流越大；γ 值越小，脉冲噪声的幅度越小，脉冲流越小。
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

[estimated_angle1, P_music1] = SQSTFT_trackmultiLRmethod_Music(data_matrix1, num_elements, d, wavelength, theta, fs, lowFreq, highFreq, alpha_param, tDS, h, Dh, r_filter, gamma_filter, delta_limit, min_length, LOWER_PRCTILE_LIMIT);
disp(['估计的到达角度1: ',num2str(estimated_angle1),'°']);
[estimated_angle2, P_music2] = DESHAPE_Music(data_matrix1, num_elements, d, wavelength, theta, fs, time_win, overlap, nfft, ths, hf, lf,r_filter,gamma_filter);
disp(['估计的到达角度2：', num2str(estimated_angle2),'°']);
[estimated_angle3, P_music3] = SST_Music(data_matrix1, fs, wavelength, d, theta, time_win, overlap, nfft, ths, hf, lf, r_filter, gamma_filter);
disp(['估计的到达角度3: ',num2str(estimated_angle3),'°']);
[estimated_angle4, P_music4] = SQCWT_Music(data_matrix1, fs, wavelength, d, theta, lf, hf, sqcwt_alpha, opts, P, r_filter, gamma_filter);
disp(['估计的到达角度4: ',num2str(estimated_angle4),'°']);
[estimated_angle5, P_music5] = SST_track_Music(data_matrix1, num_elements, d, wavelength, theta, fs, time_win, overlap, nfft, ths, hf, lf, r_filter, gamma_filter, delta_limit, min_length, LOWER_PRCTILE_LIMIT);
disp(['估计的到达角度5: ',num2str(estimated_angle5),'°']);
% [estimated_angle6, P_music6] = SQSTFT_Music(data_matrix1, wavelength, d, theta, lowFreq, highFreq, alpha_param, tDS, h, Dh, P);
% disp(['估计的到达角度6: ',num2str(estimated_angle6),'°']);
% [estimated_angle7, P_music7] = SQSTFT_track_Music(data_matrix1, num_elements, d, wavelength, theta, fs, lowFreq, highFreq, alpha_param, tDS, h, Dh, r_filter, gamma_filter, delta_limit, min_length, LOWER_PRCTILE_LIMIT);
% disp(['估计的到达角度7: ',num2str(estimated_angle7),'°']);

figure;
% 绘制 ROOT-TF-Music 谱
plot(theta, P_music1, 'LineWidth', 2, 'Color', 'r' , 'LineStyle', '-.',  'DisplayName', '\fontname{Times New Roman}\fontsize{12}ROOT-TF-Music');
hold on;
% 绘制 NEW-deShape-MUSIC 谱
plot(theta, P_music2, 'LineWidth', 2, 'Color', 'b',   'LineStyle', '-.',  'DisplayName', '\fontname{Times New Roman}\fontsize{12}DESHAPE-MUSIC');
% 绘制 SST-Music 谱
plot(theta, P_music3, 'LineWidth', 2, 'Color', 'g', 'LineStyle', ':', 'DisplayName', '\fontname{Times New Roman}\fontsize{12}SST-Music');
% 绘制 SQCWT-Music 谱
plot(theta, P_music4, 'LineWidth', 2, 'Color', 'm', 'LineStyle', '--', 'DisplayName', '\fontname{Times New Roman}\fontsize{12}SQCWT-Music');
% 绘制 SST-track-Music 谱
plot(theta, P_music5, 'LineWidth', 2, 'Color', 'r',   'LineStyle', '-', 'DisplayName', '\fontname{Times New Roman}\fontsize{12}SST-track-Music');
% 绘制 SQSTFT-Music 谱
% plot(theta, P_music6, 'LineWidth', 2, 'Color', 'k', 'LineStyle', ':', 'DisplayName', 'SQSTFT-Music');
% 绘制 SQSTFT-track-Music 谱
% plot(theta, P_music7, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-.', 'DisplayName', 'SQSTFT-track-Music');
% 绘制真实角度空间谱
yl = get(gca, 'YLim');
hold on;
% 定义样式数组
styles = {'--', '-.'}; % 不同线型
colors = {'k', 'k'};   % 不同颜色

for k = 1:length(true_angle)
    ang = true_angle(k);
    style_idx = mod(k-1, length(styles)) + 1;
    color_idx = mod(k-1, length(colors)) + 1;
    line([ang ang], yl, ...
         'Color', colors{color_idx}, ...
         'LineStyle', styles{style_idx}, ...
         'LineWidth', 2, ...
         'DisplayName', ['真实角度 ' num2str(k)]);
end

legend('show');
hold off;
title('\fontname{宋体}\fontsize{18}空间谱对比图', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('\fontname{宋体}\fontsize{12}角度\fontname{Times New Roman}\fontsize{12}/(°)', 'FontSize', 14);
ylabel('\fontname{宋体}\fontsize{12}功率\fontname{Times New Roman}\fontsize{12}/(dB)', 'FontSize', 14);
grid on; 
legend('Location', 'northwest', 'FontSize', 10); 
xticks(-90:20:90); 
xlim([-90 90]);    
yticks(-30:5:30);  
hold off;
    