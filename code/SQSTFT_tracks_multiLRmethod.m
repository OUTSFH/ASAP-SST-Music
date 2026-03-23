% 添加 PROPACK 工具箱
if exist('.\PROPACK', 'dir')==7
    addpath PROPACK;
else
    fprintf('No PROPACK installed.\n');
    error('Break, PROPACK not installed');
end

% 读取二进制文件并生成复信号
fid = fopen('QPSK-1090-2048.bin', 'rb'); 
data = fread(fid, 'float32'); 
fclose(fid);
qpsk_signal = data(1:2:end) + 1j * data(2:2:end);

% 参数设置
f = 1090e6;           % 信号中心频率 1090 MHz
c = 3e8;              % 光速（m/s）
wavelength = c / f;   % 波长
num_elements = 4;     % 阵元个数
d = wavelength * 0.5; % 阵元间距（0.5波长）
theta = -90:0.1:90;   % 扫描角度
true_angle = -45;     % 真实来波方向
SNR = 10;             % 信噪比 (dB)
fs = 30e6;            % 采样率 30 MHz
time_win = 256;       % 原始窗口长度
WinLen = time_win;
if mod(WinLen,2)==0
    WinLen = WinLen + 1;  % 确保窗长度为奇数
end
overlap = time_win/2; % 重叠长度
nfft = 256;           % FFT 点数

% 添加误差参数
amp_sigma_db = 0.5;
phase_sigma_deg = 2;

% 构建转向矩阵，并添加幅度与相位误差
P = 1;              % 信号源个数
M = num_elements;   % 阵元个数
doa = true_angle;   % 信号来波方向
lambda = wavelength;
A = zeros(P, M);

for k = 1:P
    ideal_steering = exp(-1j * 2*pi*d*sin(deg2rad(doa))*(0:M-1)/lambda);
    amp_error = 10.^(randn(1,M) * amp_sigma_db/20);
    phase_error = deg2rad(randn(1,M) * phase_sigma_deg);
    A(k,:) = amp_error .* exp(1j*phase_error) .* ideal_steering;
end
    A = A';  

% 构造数据矩阵：每个阵元为一通道
qpsk_signal = reshape(qpsk_signal, 1, []);
data_matrix1 = A * qpsk_signal;

% 添加高斯白噪声
data_matrix1 = awgn(data_matrix1, SNR, 'measured');

% 对每个通道添加稀疏干扰
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
% 低噪声环境：ths = 1 或 ths = 0.9，以保留更多细节
% 一般噪声环境：ths = 0.7 ~ 0.85，减少部分噪声干扰
% 高噪声环境：ths = 0.3 ~ 0.6，只保留主要信号成分
pulse_ratio = 0.2;    % 脉冲噪声比例
alfa = 0.75;          % 稳定分布参数 alpha
gamma_param = 1e-2;   % 稳定分布参数 gamma
num_channels = size(data_matrix1, 1);

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

% sqSTFT 时频分析
lowFreq = -0.5;    % 归一化频率下限
highFreq = 0.5;    % 归一化频率上限
alpha_param = 1/nfft; % 频率分辨率（归一化）
tDS = 1;           % 时间下采样因子
[h, Dh] = hermf(71, 1, 6); % 生成 Hermite 窗及其导数

if mod(length(h),2)==0
    h = h(1:end-1);
    Dh = Dh(1:end-1);
end

% 计算每个通道的 sqSTFT 变换（滤波前）
S_matrix_unfilt = cell(1, num_channels);
F_matrix_unfilt = cell(1, num_channels);
T_matrix_unfilt = cell(1, num_channels);

for ch = 1:num_channels
    signal = data_matrix1(ch,:).';
    [~, tfrtic, tfrsq, tfrsqtic] = sqSTFT(signal, lowFreq, highFreq, alpha_param, tDS, h', Dh');
    S_matrix_unfilt{ch} = tfrsq;
    F_matrix_unfilt{ch} = tfrsqtic * fs;
    T_matrix_unfilt{ch} = (1:length(signal))/fs;
end

% ASAP_Hankel_1D 滤波处理 
r_filter = 4;       % 目标低秩参数
gamma_filter = 0.5; % 收敛速度参数，介于0和1之间

for ch = 1:num_channels
    signal = data_matrix1(ch,:).';
    [filtered_signal, s_hat] = ASAP_Hankel_1D(signal, r_filter, gamma_filter);
    data_matrix1(ch, 1:length(filtered_signal)) = filtered_signal.';
end

% sqSTFT 时频分析（滤波后）
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

% 绘制滤波前后的 sqSTFT 时频图对比 
for ch = 1:num_channels
    figure;
    subplot(1,2,1);
    F_disp = F_matrix_unfilt{ch} / 1e6;
    db_energy_unfilt = 10*log10(abs(S_matrix_unfilt{ch}).^2 + eps);
    surf(T_matrix_unfilt{ch}, F_disp, db_energy_unfilt, 'EdgeColor','none');
    title(['(未滤波) Channel ', num2str(ch)], 'FontSize', 14);
    xlabel('Time (s)','FontSize',12);
    ylabel('Frequency (MHz)','FontSize',12);
    zlabel('Energy (dB)','FontSize',12);
    view(2); 
    clim([max(db_energy_unfilt(:))-80, max(db_energy_unfilt(:))]); 
    colorbar; colormap jet; grid on;
    caxis([-30 40]);
    
    subplot(1,2,2);
    F_disp = F_matrix_filt{ch} / 1e6;
    db_energy_filt = 10*log10(abs(S_matrix_filt{ch}).^2 + eps);
    surf(T_matrix_filt{ch}, F_disp, db_energy_filt, 'EdgeColor','none');
    title(['(滤波后) Channel ', num2str(ch)], 'FontSize', 14);
    xlabel('Time (s)','FontSize',12);
    ylabel('Frequency (MHz)','FontSize',12);
    zlabel('Energy (dB)','FontSize',12);
    view(2); 
    clim([max(db_energy_filt(:))-100, max(db_energy_filt(:))]); 
    colorbar; colormap jet; grid on;
    set(gcf, 'Position', [100,100,1200,500]);
end

% 多通道轨迹融合参数
delta_limit = 4;       % 轨迹搜索邻域大小
min_length = 15;       % 轨迹最小长度
LOWER_PRCTILE_LIMIT = 95; % 能量阈值

% 提取轨迹
[tf_pts, aligned] = tracks_multiLRmethod(S_matrix_filt,F_matrix_filt,fs,delta_limit,min_length,LOWER_PRCTILE_LIMIT,d,c);

% 清除无效轨迹
valid_rows = all(aligned~=0,2);
aligned    = aligned(valid_rows,:);
tf_pts     = tf_pts(valid_rows,:);

% 加权因子：能量因子
E = mean(abs(aligned).^2,2);
E = normalize(E, 'range');

% 加权因子：稀疏/平滑加权（用时间跳变差分做平滑度指标）
dT = abs(diff(tf_pts(:,1)));        
smoothness = 1 ./ (1 + dT);         
smoothness = [smoothness; smoothness(end)];  
SM = normalize(smoothness, 'range');

% 加权因子：自适应方向加权
theta_scan = -90:0.5:90;
R0 = (aligned' * aligned) / size(aligned,1);
a = @(theta) exp(-1j*2*pi*d*sin(theta*pi/180)*(0:size(aligned,2)-1)'/c);
P_music = arrayfun(@(th) 1./(a(th)' * (R0 - a(th)*a(th)') * a(th)), theta_scan);
[~, idx_max] = max(real(P_music));
theta0 = theta_scan(idx_max); 
v0 = a(theta0);             
v0 = v0 / norm(v0);          
sim = abs(aligned * v0);       
ADAPT = normalize(sim, 'range');

% 融合所有加权因子
weights = 0.5*E + 0.25*SM + 0.25*ADAPT;
weights = weights.^1.5;         

% 计算加权协方差矩阵
Xw = aligned .* sqrt(weights);      
R  = (Xw' * Xw) / sum(weights);    

% 在各通道滤波后 sqSTFT 图上叠加显示融合后的轨迹点-true/false
useAlphaDisplay = true; 
for ch = 1:M
    figure;
    db_sst = 10*log10(abs(S_matrix_filt{ch}).^2 + eps);
    imagesc(T_matrix_filt{ch}, F_matrix_filt{ch}/1e6, db_sst);
    axis xy; colormap jet; colorbar;
    title(sprintf('Channel %d', ch));
    xlabel('Time (s)'); ylabel('Freq (MHz)');
    caxis([max(db_sst(:))-100, max(db_sst(:))]);
    hold on;
    t_vals = T_matrix_filt{ch}(tf_pts(:,1));
    f_vals = F_matrix_filt{ch}(tf_pts(:,2)) / 1e6;
    
    if useAlphaDisplay
        scatter(t_vals, f_vals, ...
            40, 'k', 'filled', ...
            'MarkerFaceAlpha', 'flat', ...
            'AlphaData', weights, ...
            'AlphaDataMapping','scaled');
    else

        for p = 1:size(tf_pts,1)
            t_idx = tf_pts(p,1);
            f_idx = tf_pts(p,2);
            if t_idx <= length(T_matrix_filt{ch}) && f_idx <= size(S_matrix_filt{ch},1)
                plot(T_matrix_filt{ch}(t_idx), F_matrix_filt{ch}(f_idx)/1e6, ...
                     'ko','MarkerSize',4,'MarkerFaceColor','k');
            end
        end
    end
    hold off;
end

[eigenVec, eigenVal] = eig(R);
[eigenVal, idx] = sort(diag(eigenVal), 'descend');
eigenVec = eigenVec(:, idx);

% 选择信号子空间
num_signals = 1;
signalVec = eigenVec(:, 1:num_signals);
noiseVec = eigenVec(:, num_signals+1:end);

% 计算 MUSIC 谱
P_tf_music = zeros(1, length(theta));
for ii = 1:length(theta)
    a = exp(-1j*(0:num_elements-1)'*(2*pi/wavelength)*d*sin(deg2rad(theta(ii))));
    noiseMatrix = noiseVec * noiseVec';
    P_tf_music(ii) = 1 / (a' * noiseMatrix * a);
end

P_tf_music = 10*log10(abs(P_tf_music)/max(abs(P_tf_music)));
[~, max_idx] = max(P_tf_music);
estimated_angle = theta(max_idx);
disp(['估计的到达角度：', num2str(estimated_angle), '°']);

figure;
plot(theta, P_tf_music, 'LineWidth', 2, 'Color', 'b');
title('SQSTFT-tracks_multiLRmethod-MUSIC Spectrum Estimation', 'FontSize',16,'FontWeight','bold');
xlabel('Angle (degrees)','FontSize',16);
ylabel('MUSIC Spectrum (dB)','FontSize',16);
xticks(-90:5:90); yticks(-20:5:20); grid on;
legend('SQSTFT-tracks_multiLRmethod-MUSIC Spectrum','Location','Best','FontSize',12);
