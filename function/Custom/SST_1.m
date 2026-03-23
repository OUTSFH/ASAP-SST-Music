function [sst, tfr, frequency, t_vector] = SST_1(x, Fs, hlength, hop, n, hf, lf, ths) 
% Computes the synchrosqueezing transform (SST) of the signal x
% INPUT
%    x      :  信号（列向量）
%    Fs     :  采样率
%    hlength:  窗口长度（样本数）
%    hop    :  每隔 hop 个样本计算一次 FFT
%    n      :  FFT 点数（频率轴像素数）
%    hf     :  输出显示的最高频率（Hz）
%    lf     :  输出显示的最低频率（Hz）
%    ths    :  重新分配的阈值参数（取值 0~1）
% OUTPUT
%    sst      :  SST 结果
%    tfr      :  对应的 STFT 结果
%    frequency:  频率轴（Hz）
%    t_vector :  时间轴（秒）

switch nargin
    case 7
        ths = 1;
    case 6 
        ths = 1;
        lf = 0;
    case 5
        ths = 1;
        hf = inf;
        lf = 0;
    case 4
        ths = 1;
        n = pow2(nextpow2(length(x))) / 2 + 1;
        hf = inf;
        lf = 0;
    case 3
        ths = 1;
        hop = 1;
        n = pow2(nextpow2(length(x))) / 2 + 1;
        hf = inf;
        lf = 0;
    case 2
        error('请指定窗口长度。')
    case 1
        error('请指定采样率。')
    case 0
        Fs = 200;
        x = 2 * mod(1e-2:1e-2:1e2, 1) - 1;
        hlength = 1001;
        lf = 1;
        hop = 40;
        n = 8000;
        hf = 12;
        ths = 0.5;
        disp('测试代码：2 Hz 锯齿波。')
end

% 窗口带宽参数sigma决定了高斯窗的“宽度”
% sigma 越小，使得指数衰减更快，得到一个窄窗；sigma 越大，则得到一个宽窗
sigma = 0.15;
squeeze_flag = 1;

x = x(:);
if any(isnan(x))
    x = interp1(find(~isnan(x)), x(~isnan(x)), 1:length(x), 'pchip', 'extrap')';
end

NN = length(x);
t_idx = 1:hop:NN;
num_frames = length(t_idx);
t_vector = t_idx / Fs;  

% 保证 n 为偶数
n = n + 1 - rem(n, 2);
N = 2 * (n - 1);

% 保证窗口长度为奇数
hlength = hlength + 1 - rem(hlength, 2);
Lh = (hlength - 1) / 2;

% 构造高斯窗及其导数
ex = linspace(-0.5, 0.5, hlength)';
h_win = exp(-ex.^2/(2*sigma^2));
dh = -ex/(sigma^2) .* h_win;

% 初始化时频矩阵
tfr = zeros(N, num_frames);
tfr2 = zeros(N, num_frames);

for icol = 1:num_frames
    ti = t_idx(icol);
    tau = -min([N-1, Lh, ti-1]):min([N-1, Lh, NN-ti]);
    indices = mod(N + tau, N) + 1;
    rSig = x(ti + tau);
    tfr(indices, icol) = rSig .* h_win(Lh+1+tau);
    tfr2(indices, icol) = rSig .* dh(Lh+1+tau);
end

% FFT 及 fftshift 得到对称频谱
tfr = fft(tfr, N, 1);
tfr = fftshift(tfr, 1);
if squeeze_flag
    tfr2 = fft(tfr2, N, 1);
    tfr2 = fftshift(tfr2, 1);
end

% 构造频率轴（-Fs/2 到 Fs/2）
frequency = linspace(-Fs/2, Fs/2, N)';

% 按指定频率范围裁剪
valid = (frequency >= lf) & (frequency <= hf);
tfr = tfr(valid, :);
if squeeze_flag
    tfr2 = tfr2(valid, :);
end
frequency = frequency(valid);
neta = length(frequency);

if ~squeeze_flag
    sst = tfr;
    return
end

% 计算重新分配的门限（全局统计）
ths_val = quantile(abs(tfr(:)), (1 - ths));

% 计算 omega 算子
df = frequency(2) - frequency(1);
omega = -inf(neta, num_frames);
ups = abs(tfr) > ths_val;
index = repmat((1:neta)', [1, num_frames]);
offset = round( imag(tfr2./tfr) ./ (2*pi*df) );
omega(ups) = index(ups) + offset(ups);
id = omega < 1 | omega > neta | ~ups;
omega(id) = index(id);

% 重分配：利用 accumarray 快速实现
id_linear = omega + neta * (0:num_frames-1);
sst = reshape(accumarray(id_linear(:), tfr(:), [neta*num_frames, 1]), [neta, num_frames]);

end
