% 读取二进制文件
fid = fopen('QPSK-1090-1024.bin', 'rb'); 
data = fread(fid, 'float32'); 
fclose(fid);
qpsk_signal1 = data(1:2:end) + 1j * data(2:2:end);

fid = fopen('xinhao1-1024.bin', 'rb'); 
data = fread(fid, 'float32'); 
fclose(fid);
qpsk_signal2 = data(1:2:end) + 1j * data(2:2:end);

% 参数设置
f = 1090e6;           % 信号中心频率 1090 MHz
c = 3e8;              % 光速（米/秒）
wavelength = c / f;   % 波长
num_elements = 4;     % 天线阵元数量
d = wavelength * 0.5; % 阵元间距为波长的0.5倍
theta = -90:0.1:90;   % 扫描角度
true_angle = [-45, 15]; % 真实来波方向
SNR = 10;             % 信噪比
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

% 添加误差参数
amp_sigma_db = 0.5;   % 幅度误差标准差，单位dB
phase_sigma_deg = 2;  % 相位误差标准差，单位度

% 构建转向矩阵，并添加幅度和相位误差
P = 2;                % 信号源数量
M = num_elements;     % 阵元数量
doa = true_angle;     % 信号源来波方向
lambda = wavelength;  % 波长
A = zeros(M, P);
idx = (0:M-1)';        % 阵元索引列向量

for k = 1:P
    doa_k      = true_angle(k);
    ideal_steer = exp(-1j*2*pi*d*sin(deg2rad(doa_k)).*idx/wavelength);
    amp_err    = (10.^(randn(1,M)*amp_sigma_db/20)).';  % M×1
    phase_err  = deg2rad(randn(1,M)*phase_sigma_deg).';% M×1
    A(:,k)     = amp_err.*exp(1j*phase_err).*ideal_steer;
end

% 调整 qpsk_signal 的维度以匹配矩阵乘法
qpsk_signal1 = reshape(qpsk_signal1, 1, []);
qpsk_signal2 = reshape(qpsk_signal2, 1, []);
S = [qpsk_signal1;
     qpsk_signal2];

% 生成信号矩阵，每个阵元作为一通道
data_matrix1 = A * S;  

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

% for ch = 1:num_channels
%     signal = data_matrix1(ch, :).';
%     n = length(signal);
%     IND = find(rand(n,1) < pulse_ratio);
%     if ~isempty(IND)
%         pulse_noise = zeros(length(IND),1);
%         for idx = 1:length(IND)
%             r1 = stblrnd(alfa, 0, gamma_param, 0, 1, 1);
%             r2 = stblrnd(alfa, 0, gamma_param, 0, 1, 1);
%             pulse_noise(idx) = sqrt(1/2) * (r1 + 1i*r2);
%         end
%         signal(IND) = signal(IND) + pulse_noise;
%     end
%     data_matrix1(ch, :) = signal.';
% end

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';
    n = length(signal);

    % 计算原始信号功率
    signal_power = mean(abs(signal).^2);
   
    % 根据设定的SNR计算脉冲噪声功率
    noise_power = signal_power / (10^(SNR / 10));
    
    % 生成脉冲噪声（满足稳定分布）
    r1 = stblrnd(alfa, 0, gamma_param, 0, n, 1);
    r2 = stblrnd(alfa, 0, gamma_param, 0, n, 1);
    pulse_noise = sqrt(noise_power / 2) * (r1 + 1i * r2);

    % 直接将脉冲噪声添加到信号中
    signal = signal + pulse_noise;
    
    % 更新信号矩阵
    data_matrix1(ch, :) = signal.';
end

% STFT-1
num_channels = size(data_matrix1, 1);
S_matrix = cell(1, num_channels);  
F_matrix = cell(1, num_channels); 
T_matrix = cell(1, num_channels);  

for ch = 1:num_channels
    signal = data_matrix1(ch, :).'; 
    [S_matrix{ch}, F_matrix{ch}, T_matrix{ch}] = cus_stft(signal, fs, window, overlap, nfft);
end

figure('Position', [100, 400, 1000, 300]); 
for ch = 1:num_channels
    left = 0.08 + (ch - 1) * 0.22; 
    bottom = 0.15;
    width = 0.2;
    height = 0.8;
    subplot('Position', [left, bottom, width, height]); 
    F = F_matrix{ch} / 1e6; 
    T = T_matrix{ch};
    S_db = 10 * log10(abs(S_matrix{ch}).^2);
    surf(T, F, S_db, 'EdgeColor', 'none');
    ax = gca; 
    set(gca, 'XTick', [0  1e-5  2e-5  3e-5], 'XTickLabels', {'0','1','2','3'});
    set(gca, 'XLim', [0 3.2e-5]);   
    title(' ', 'FontSize', 12);
    xlabel('\fontname{宋体}\fontsize{10}时间/\fontname{Times New Roman}\fontsize{10}s', 'FontSize', 10);
    if ch == 1
        ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz', 'FontSize', 10);
    else
        h_ylabel = ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz', 'FontSize', 10);
        set(h_ylabel, 'Visible', 'off');
    end
    zlabel('Energy (dB)', 'FontSize', 10);
    view(2); 
    colorbar('vertical', 'Location', 'eastoutside'); 
    colormap jet;
    caxis([-30 40]);
    grid on;
    set(gca, 'FontSize', 10, 'LineWidth', 0.8);
end

% SST-2 
S_matrix_1 = cell(1, num_channels);    
F_matrix_1 = cell(1, num_channels);   
T_matrix_1 = cell(1, num_channels);    

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';
    [sst, ~, frequency, t_vector] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
    S_matrix_1{ch} = sst;
    F_matrix_1{ch} = frequency / 1e6; 
    T_matrix_1{ch} = t_vector; 
end

figure('Position', [100, 400, 1000, 300]); 
for ch = 1:num_channels
    left = 0.08 + (ch - 1) * 0.22; 
    bottom = 0.15;
    width = 0.2;
    height = 0.8;
    subplot('Position', [left, bottom, width, height]); 
    F = F_matrix_1{ch} ; 
    T = T_matrix_1{ch};
    S_db = 10 * log10(abs(S_matrix_1{ch}).^2);
    surf(T, F, S_db, 'EdgeColor', 'none');
    ax = gca; 
    set(gca, 'XTick', [0  1e-5  2e-5  3e-5], 'XTickLabels', {'0','1','2','3'});
    set(gca, 'XLim', [0 3.2e-5]);   
    title(' ', 'FontSize', 12);
    xlabel('\fontname{宋体}\fontsize{10}时间/\fontname{Times New Roman}\fontsize{10}s', 'FontSize', 10);
    if ch == 1
        ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz', 'FontSize', 10);
    else
        h_ylabel = ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz', 'FontSize', 10);
        set(h_ylabel, 'Visible', 'off');
    end
    zlabel('Energy (dB)', 'FontSize', 10);
    view(2); 
    colorbar('vertical', 'Location', 'eastoutside'); 
    colormap jet;
    caxis([-30 40]);
    grid on;
    set(gca, 'FontSize', 10, 'LineWidth', 0.8);
end

% SQSTFT-3 
S_matrix_2 = cell(1, num_channels);
F_matrix_2 = cell(1, num_channels);
T_matrix_2 = cell(1, num_channels);

for ch = 1:num_channels
    signal = data_matrix1(ch,:).';
    [~, tfrtic, tfrsq, tfrsqtic] = sqSTFT(signal, lowFreq, highFreq, alpha_param, tDS, h', Dh');
    S_matrix_2{ch} = tfrsq;
    F_matrix_2{ch} = tfrsqtic * fs;
    T_matrix_2{ch} = (1:length(signal))/fs;
end

figure('Position', [100, 400, 1000, 300]); 
for ch = 1:num_channels
    left = 0.08 + (ch - 1) * 0.22; 
    bottom = 0.15;
    width = 0.2;
    height = 0.8;
    subplot('Position', [left, bottom, width, height]); 
    F_disp = F_matrix_2{ch} / 1e6;
    db_energy_unfilt = 10*log10(abs(S_matrix_2{ch}).^2 + eps);
    surf(T_matrix_2{ch}, F_disp, db_energy_unfilt, 'EdgeColor', 'none');
    axis tight;
    ax = gca; 
    set(gca, 'XTick', [0  1e-5  2e-5  3e-5], 'XTickLabels', {'0','1','2','3'});
    set(gca, 'XLim', [0 3.2e-5]); 
    title(' ', 'FontSize', 14);
    xlabel('\fontname{宋体}\fontsize{10}时间/\fontname{Times New Roman}\fontsize{10}s', 'FontSize', 12);
    if ch == 1
        ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz', 'FontSize', 12);
    else
        h_ylabel = ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz', 'FontSize', 12);
        set(h_ylabel, 'Visible', 'off');
    end
    zlabel('Energy (dB)', 'FontSize', 12);
    view(2);
    clim([max(db_energy_unfilt(:)) - 80, max(db_energy_unfilt(:))]);
    colorbar('vertical', 'Location', 'eastoutside'); 
    colormap jet;
    grid on;
    caxis([-30 40]);
end

% SST_track-4
% ASAP_Hankel_1D 滤波处理
r_filter = 4;       
gamma_filter = 0.5; 

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';  
    [filtered_signal, ~] = ASAP_Hankel_1D(signal, r_filter, gamma_filter);
    data_matrix1(ch, 1:length(filtered_signal)) = filtered_signal.';
end

% 滤波后SST 
S_matrix_3 = cell(1, num_channels);   
F_matrix_3 = cell(1, num_channels);   
T_matrix_3 = cell(1, num_channels);   

for ch = 1:num_channels
    signal = data_matrix1(ch, :).';
    [sst, ~, frequency, t_vector] = SST_1(signal, fs, time_win, overlap, nfft, hf, lf, ths);
    S_matrix_3{ch} = sst;
    F_matrix_3{ch} = frequency / 1e6; 
    T_matrix_3{ch} = t_vector;
end

figure('Position', [100, 400, 1000, 300]); 
for ch = 1:num_channels
    left = 0.08 + (ch - 1) * 0.22; 
    bottom = 0.15;
    width = 0.2;
    height = 0.8;
    subplot('Position', [left, bottom, width, height]); 
    F = F_matrix_3{ch} ; 
    T = T_matrix_3{ch};
    S_db = 10 * log10(abs(S_matrix_3{ch}).^2);
    surf(T, F, S_db, 'EdgeColor', 'none');
    ax = gca; 
    set(gca, 'XTick', [0  1e-5  2e-5  3e-5], 'XTickLabels', {'0','1','2','3'});
    set(gca, 'XLim', [0 3.2e-5]); 
    title(' ', 'FontSize', 12);
    xlabel('\fontname{宋体}\fontsize{10}时间/\fontname{Times New Roman}\fontsize{10}s', 'FontSize', 10);
    if ch == 1
        ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz', 'FontSize', 10);
    else
        h_ylabel = ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz', 'FontSize', 10);
        set(h_ylabel, 'Visible', 'off');
    end
    zlabel('Energy (dB)', 'FontSize', 10);
    view(2); 
    colorbar('vertical', 'Location', 'eastoutside'); 
    colormap jet;
    caxis([-30 40]);
    grid on;
    set(gca, 'FontSize', 10, 'LineWidth', 0.8);
end

% 多通道轨迹融合
delta_limit = 4;       % 轨迹搜索邻域大小
min_length = 15;       % 轨迹最小长度
LOWER_PRCTILE_LIMIT = 95; % 能量阈值

% 对每个通道均独立提取轨迹，并融合各通道轨迹点
all_time_freq_points = [];
for ch = 1:num_channels
    sst_ref = abs(S_matrix_3{ch}).';
    [individual_tracks, ~] = tracks_LRmethod(sst_ref, fs, delta_limit, min_length, LOWER_PRCTILE_LIMIT);
    for i = 1:length(individual_tracks)
        track = individual_tracks{i}; 
        all_time_freq_points = [all_time_freq_points; track];
    end
end

% 融合后的轨迹点去重处理
time_freq_points = unique(all_time_freq_points, 'rows');
stft_data = [];
for ch = 1:num_channels
    sst_ch = S_matrix_3{ch};  
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

figure('Position', [100, 400, 1000, 300]); 
for ch = 1:num_channels
    left = 0.08 + (ch - 1) * 0.22; 
    bottom = 0.15;
    width = 0.2;
    height = 0.8;
    subplot('Position', [left, bottom, width, height]); 
    sst_plot = 10*log10(abs(S_matrix_3{ch}).^2);
    imagesc(T_matrix_3{ch}, F_matrix_3{ch}, sst_plot);
    axis xy;
    colormap jet;
    colorbar;
    ax = gca; 
    set(gca, 'XTick', [0  1e-5  2e-5  3e-5], 'XTickLabels', {'0','1','2','3'});
    set(gca, 'XLim', [0 3.2e-5]);   
    title(' ', 'FontSize', 14);
    xlabel('\fontname{宋体}\fontsize{10}时间/\fontname{Times New Roman}\fontsize{10}s');
    if ch == 1
        ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz');
    else
        h_ylabel = ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz');
        set(h_ylabel, 'Visible', 'off');
    end
    hold on;
    caxis([-30 40]);

    for p = 1:size(time_freq_points, 1)
        t_idx = time_freq_points(p, 1);
        f_idx = time_freq_points(p, 2);
        if t_idx <= length(T_matrix_3{ch}) && f_idx <= length(F_matrix_3{ch})
            t_val = T_matrix_3{ch}(t_idx);
            f_val = F_matrix_3{ch}(f_idx);
            plot(t_val, f_val, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
        end
    end
    hold off;
end

% SQSTFT_filter-5
S_matrix_3 = cell(1, num_channels);
F_matrix_3 = cell(1, num_channels);
T_matrix_3 = cell(1, num_channels);

for ch = 1:num_channels
    signal = data_matrix1(ch,:).';
    [~, tfrtic, tfrsq, tfrsqtic] = sqSTFT(signal, lowFreq, highFreq, alpha_param, tDS, h', Dh');
    S_matrix_3{ch} = tfrsq;
    F_matrix_3{ch} = tfrsqtic * fs/ 1e6;
    T_matrix_3{ch} = (1:length(signal))/fs;
end

figure('Position', [100, 400, 1000, 300]); 
for ch = 1:num_channels
    left = 0.08 + (ch - 1) * 0.22; 
    bottom = 0.15;
    width = 0.2;
    height = 0.8;
    subplot('Position', [left, bottom, width, height]); 
    F_disp = F_matrix_3{ch} ;
    db_energy_unfilt = 10*log10(abs(S_matrix_3{ch}).^2 + eps);
    surf(T_matrix_3{ch}, F_disp, db_energy_unfilt, 'EdgeColor', 'none');
    axis tight;
    ax = gca; 
    set(gca, 'XTick', [0  1e-5  2e-5  3e-5], 'XTickLabels', {'0','1','2','3'});
    set(gca, 'XLim', [0 3.2e-5]);  
    title(' ', 'FontSize', 14);
    xlabel('\fontname{宋体}\fontsize{10}时间/\fontname{Times New Roman}\fontsize{10}s', 'FontSize', 12);
    if ch == 1
        ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz', 'FontSize', 12);
    else
        h_ylabel = ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz', 'FontSize', 12);
        set(h_ylabel, 'Visible', 'off');
    end
    zlabel('Energy (dB)', 'FontSize', 12);
    view(2);
    clim([max(db_energy_unfilt(:)) - 80, max(db_energy_unfilt(:))]);
    colorbar('vertical', 'Location', 'eastoutside'); 
    colormap jet;
    grid on;
    caxis([-30 40]);
end

% % 对每个通道均独立提取轨迹，并融合各通道轨迹点
% all_time_freq_points = [];
% for ch = 1:num_channels
%     sst_ref = abs(S_matrix_4{ch}).';
%     [individual_tracks, ~] = tracks_LRmethod(sst_ref, fs, delta_limit, min_length, LOWER_PRCTILE_LIMIT);
%     for i = 1:length(individual_tracks)
%         track = individual_tracks{i}; 
%         all_time_freq_points = [all_time_freq_points; track];
%     end
% end
% 
% % 融合后的轨迹点去重处理
% time_freq_points = unique(all_time_freq_points, 'rows');
% stft_data = [];
% for ch = 1:num_channels
%     sst_ch = S_matrix_4{ch};  
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
% % 去掉任一通道为 0 的轨迹点
% valid_rows = all(stft_data ~= 0, 2);
% stft_data = stft_data(valid_rows, :);
% 
% % SQSTFT_filter_track-6
% figure('Position', [100, 400, 1000, 300]); 
% for ch = 1:num_channels
%     left = 0.08 + (ch - 1) * 0.22; 
%     bottom = 0.15;
%     width = 0.2;
%     height = 0.8;
%     subplot('Position', [left, bottom, width, height]); 
%     sst_plot = 10*log10(abs(S_matrix_4{ch}).^2);
%     imagesc(T_matrix_4{ch}, F_matrix_4{ch}, sst_plot);
%     axis xy;
%     colormap jet;
%     colorbar;
%     ax = gca; 
%     set(gca, 'XTick', [0 2e-5 4e-5 6e-5], 'XTickLabels', {'0','2','4','6'});
%     set(gca, 'XLim', [0 6e-5]);  
%     title(' ', 'FontSize', 14);
%     xlabel('\fontname{宋体}\fontsize{10}时间/\fontname{Times New Roman}\fontsize{10}s');
%     if ch == 1
%         ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz');
%     else
%         h_ylabel = ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz');
%         set(h_ylabel, 'Visible', 'off');
%     end
%     hold on;
%     caxis([-30 40]);
% 
%     for p = 1:size(time_freq_points, 1)
%         t_idx = time_freq_points(p, 1);
%         f_idx = time_freq_points(p, 2);
%         if t_idx <= length(T_matrix_4{ch}) && f_idx <= length(F_matrix_4{ch})
%             t_val = T_matrix_4{ch}(t_idx);
%             f_val = F_matrix_4{ch}(f_idx);
%             plot(t_val, f_val, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
%         end
%     end
%     hold off;
% end

% % SQSTFT_filter_track1-7
% S_matrix_5 = cell(1, num_channels);
% F_matrix_5 = cell(1, num_channels);
% T_matrix_5 = cell(1, num_channels);
% 
% for ch = 1:num_channels
%     signal = data_matrix1(ch,:).';
%     [~, tfrtic, tfrsq, tfrsqtic] = sqSTFT(signal, lowFreq, highFreq, alpha_param, tDS, h', Dh');
%     S_matrix_5{ch} = tfrsq;
%     F_matrix_5{ch} = tfrsqtic * fs/ 1e6;
%     T_matrix_5{ch} = (1:length(signal))/fs;
% end
% 
% % 提取轨迹
% [tf_pts, aligned] = tracks_multiLRmethod(S_matrix_5,F_matrix_5,fs,delta_limit,min_length,LOWER_PRCTILE_LIMIT,d,c);
% 
% % 清除无效轨迹
% valid_rows = all(aligned~=0,2);
% aligned    = aligned(valid_rows,:);
% tf_pts     = tf_pts(valid_rows,:);
% 
% % 加权因子：能量因子
% E = mean(abs(aligned).^2,2);
% E = normalize(E, 'range');
% 
% % 加权因子：稀疏/平滑加权（用时间跳变差分做平滑度指标）
% dT = abs(diff(tf_pts(:,1)));        
% smoothness = 1 ./ (1 + dT);         
% smoothness = [smoothness; smoothness(end)];  
% SM = normalize(smoothness, 'range');
% 
% % 加权因子：自适应方向加权
% theta_scan = -90:0.5:90;
% R0 = (aligned' * aligned) / size(aligned,1);
% a = @(theta) exp(-1j*2*pi*d*sin(theta*pi/180)*(0:size(aligned,2)-1)'/c);
% P_music = arrayfun(@(th) 1./(a(th)' * (R0 - a(th)*a(th)') * a(th)), theta_scan);
% [~, idx_max] = max(real(P_music));
% theta0 = theta_scan(idx_max); 
% v0 = a(theta0);             
% v0 = v0 / norm(v0);          
% sim = abs(aligned * v0);       
% ADAPT = normalize(sim, 'range');
% 
% % 融合所有加权因子
% gamma_val = 1.5;
% weights = 0.5 * E + 0.25 * SM + 0.25 * ADAPT;
% weights = weights .^ gamma_val;
% weights = normalize(weights, 'range');  
% useAlphaDisplay = true;  
% 
% figure('Position', [100, 400, 1000, 300]);  
% for ch = 1:M
%     left = 0.08 + (ch - 1) * 0.22; 
%     bottom = 0.15;
%     width = 0.2;
%     height = 0.8;
%     subplot('Position', [left, bottom, width, height]); 
%     db_sst = 10*log10(abs(S_matrix_5{ch}).^2 + eps);
%     imagesc(T_matrix_5{ch}, F_matrix_5{ch}, db_sst);
%     axis xy;
%     colormap jet;
%     colorbar;
%     set(gca, 'XTick', [0 2e-5 4e-5 6e-5], 'XTickLabels', {'0','2','4','6'});
%     set(gca, 'XLim', [0 6e-5]);  
%     title(' ', 'FontSize', 14);
%     xlabel('\fontname{宋体}\fontsize{10}时间/\fontname{Times New Roman}\fontsize{10}s');
%     if ch == 1
%         ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz');
%     else
%         h_ylabel = ylabel('\fontname{宋体}\fontsize{10}频率/\fontname{Times New Roman}\fontsize{10}MHz');
%         set(h_ylabel, 'Visible', 'off');
%     end
%     caxis([max(db_sst(:)) - 100, max(db_sst(:))]);
%     hold on;
%     caxis([-30 40]);
%     t_vals = T_matrix_5{ch}(tf_pts(:,1));
%     f_vals = F_matrix_5{ch}(tf_pts(:,2));
% 
%     if useAlphaDisplay
%         scatter(t_vals, f_vals, ...
%             20, 'k', 'filled', ...
%             'MarkerFaceAlpha', 'flat', ...
%             'AlphaData', weights, ...
%             'AlphaDataMapping','scaled');
%     else
% 
%         for p = 1:size(tf_pts,1)
%             t_idx = tf_pts(p,1);
%             f_idx = tf_pts(p,2);
%             if t_idx <= length(T_matrix_5{ch}) && f_idx <= size(S_matrix_5{ch},1)
%                 plot(T_matrix_5{ch}(t_idx), F_matrix_5{ch}(f_idx), ...
%                      'ko','MarkerSize', 3, 'MarkerFaceColor','k');
%             end
%         end
%     end
%     hold off;
% end

