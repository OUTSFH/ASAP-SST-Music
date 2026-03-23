function [deshape, ceps, mask, tfr, frequency, quefrency] = deShape_2(x, Fs, hlength, hf, gamma, hop, n, lf, ths)
% Computes the de-shape synchrosqueezing transform of the signal x
% INPUT
%    x      :  信号（列向量）
%    Fs     :  采样率
%    hlength:  窗口长度（样本数）
%    hf     :  输出频率上限（Hz）
%    gamma  :  Cepstral 幂次
%    hop    :  每 hop 个样本计算一次 FFT
%    n      :  FFT 点数
%    lf     :  输出频率下限（Hz），可为负；倒谱处理时取 abs(lf)
%    ths    :  重分配的阈值参数（分位数）
% OUTPUT
%    deshape:  de-shape SST（时频重分布）
%    ceps   :  短时倒谱
%    mask   :  倒谱反变换得到的掩模
%    tfr    :  STFT（经倒谱处理前）的时频表示
%    frequency: 对称频率轴（Hz）
%    quefrency: 倒谱轴（s）

lf_positive = abs(lf);

if nargin < 9
    ths = 1;
end

if nargin < 8
    error('deShape_1: 请设置 FFT 下界 lf.');
end

if nargin < 7
    n = 2^(nextpow2(length(x))) + 1;
end

if nargin < 6
    hop = 1;
end

% window bandwidth
sigma = 0.15;
squeeze_flag = 1;

% 整理输入，确保为列向量
x = x(:);
if any(isnan(x))
    x = interp1(find(~isnan(x)), x(~isnan(x)), 1:length(x), 'pchip', 'extrap');
end

NN = length(x);
t = 1:hop:NN;
num_frames = length(t);

% 设定 FFT 长度（取偶数）
n = n + 1 - rem(n,2);
N = n;

% 保证窗口长度为奇数
hlength = hlength + 1 - rem(hlength,2);
Lh = (hlength - 1) / 2;

% 生成高斯窗及其一阶导数
ex = linspace(-0.5, 0.5, hlength)';
h_win = exp(-ex.^2/(2*sigma^2));
dh = -ex/(sigma^2) .* h_win;

% 初始化 STFT 计算所需矩阵
tfr_temp = zeros(N, num_frames);
tfr2_temp = zeros(N, num_frames);

for icol = 1:num_frames
    ti = t(icol);
    tau = -min([N-1, Lh, ti-1]):min([N-1, Lh, NN-ti]);
    indices = mod(N + tau, N) + 1;
    segment = x(ti+tau);
    tfr_temp(indices, icol) = segment .* h_win(Lh+1+tau);
    tfr2_temp(indices, icol) = segment .* dh(Lh+1+tau);
end

% 采用全 FFT 并 fftshift 得到对称频率轴
tfr_full = fft(tfr_temp, N, 1);
tfr_full = fftshift(tfr_full, 1);
tfr2_full = fft(tfr2_temp, N, 1);
tfr2_full = fftshift(tfr2_full, 1);

% 构造频率向量（根据 N 的奇偶性）
if mod(N,2)==0
    k = (-N/2 : N/2-1).';
else
    k = (-(N-1)/2 : (N-1)/2).';
end
frequency = k * (Fs/N);

% 计算短时倒谱（
ceps_full = ifft(abs(fft(tfr_temp, N, 1)).^gamma, N, 1);
ceps_full = real(ceps_full);

% 倒谱时间轴
quefrency = (0:N-1)'/Fs;

% 构造倒谱掩模 mask 
xorig = 1./(quefrency(2:end)+eps);  
yorig = ceps_full(2:end, :);        
pos_idx = find(frequency >= 0);
M = length(pos_idx);  
xq = abs(frequency(pos_idx)); 
xinterp = flipud(xorig);
yinterp = flipud(yorig);
mask_part = interp1(xinterp, yinterp, xq, 'linear', 0);  

% 构造全掩模 mask
if mod(N,2)==0
    mask = [flipud(mask_part); mask_part];
else
    mask = [flipud(mask_part(2:end,:)); mask_part];
end

% 检查 mask 尺寸是否与 tfr_full 匹配
if ~isequal(size(tfr_full), size(mask))
    error('掩模 mask 尺寸 [%d %d] 与 tfr_full 尺寸 [%d %d] 不匹配', ...
        size(mask,1), size(mask,2), size(tfr_full,1), size(tfr_full,2));
end

% 计算 de-shape STFT
deshape = tfr_full .* mask;

if ~squeeze_flag
    return
end

% 重分配操作：先设定阈值
th_val = quantile(abs(tfr_full(:)), (1-ths));
[N_freq, ~] = size(tfr_full);
omega = -inf(N_freq, num_frames);
valid = abs(tfr_full) > th_val;
omega(valid) = round( (N/hlength) * imag(tfr2_full(valid) ./ (tfr_full(valid))/(2*pi) ) );

% 计算新的频率索引
idx_mat = repmat((1:N_freq).', 1, num_frames);
new_idx = idx_mat - omega;
new_idx(new_idx < 1 | new_idx > N_freq) = idx_mat(new_idx < 1 | new_idx > N_freq);

% 重分配累加
deshape_rs = zeros(size(deshape));
for col = 1:num_frames
    for row = 1:N_freq
        idx = new_idx(row, col);
        deshape_rs(idx, col) = deshape_rs(idx, col) + deshape(row, col);
    end
end
deshape = deshape_rs;

% 按照用户设定的频率范围 [lf, hf] 裁剪时频表示
f_low = lf;  
f_high = hf;
idx_crop = (frequency >= f_low) & (frequency <= f_high);
tfr = tfr_full(idx_crop, :);
mask = mask(idx_crop, :);
deshape = deshape(idx_crop, :);
frequency = frequency(idx_crop);

% 对倒谱也按 lf_positive 进行裁剪
ceps = ceps_full(quefrency <= 2/lf_positive, :);
quefrency = quefrency(quefrency <= 2/lf_positive);

end
