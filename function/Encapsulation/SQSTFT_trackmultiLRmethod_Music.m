% function [estimated_angle, P_tf_music] = SQSTFT_trackmultiLRmethod_Music(data_matrix1, num_elements, d, wavelength, theta, fs, lowFreq, highFreq, alpha_param, tDS, h, Dh, r_filter, gamma_filter, delta_limit, min_length, LOWER_PRCTILE_LIMIT)
% 
% num_channels = size(data_matrix1, 1);        
% c = 3e8;         
% for ch = 1:num_channels
%     signal = data_matrix1(ch,:).';
%     [filtered_signal, s_hat] = ASAP_Hankel_1D(signal, r_filter, gamma_filter);
%     data_matrix1(ch, 1:length(filtered_signal)) = filtered_signal.';
% end
% 
% S_matrix_filt = cell(1, num_channels);
% F_matrix_filt = cell(1, num_channels);
% T_matrix_filt = cell(1, num_channels);
% 
% for ch = 1:num_channels
%     signal = data_matrix1(ch,:).';
%     [~, tfrtic, tfrsq, tfrsqtic] = sqSTFT(signal, lowFreq, highFreq, alpha_param, tDS, h', Dh');
%     S_matrix_filt{ch} = tfrsq;
%     F_matrix_filt{ch} = tfrsqtic * fs;
%     T_matrix_filt{ch} = (1:length(signal))/fs;
% end
% 
% delta_limit = 4;      
% min_length = 15;      
% LOWER_PRCTILE_LIMIT = 95; 
% [tf_pts, aligned] = tracks_multiLRmethod(S_matrix_filt,F_matrix_filt,fs,delta_limit,min_length,LOWER_PRCTILE_LIMIT,d,c);
% 
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
% weights = 0.5*E + 0.25*SM + 0.25*ADAPT;
% weights = weights.^1.5;         
% 
% % 计算加权协方差矩阵
% Xw = aligned .* sqrt(weights);      
% R  = (Xw' * Xw) / sum(weights);  
% [eigenVec, eigenVal] = eig(R);
% [eigenVal, idx] = sort(diag(eigenVal), 'descend');
% eigenVec = eigenVec(:, idx);
% num_signals = 1;
% signalVec = eigenVec(:, 1:num_signals);
% noiseVec = eigenVec(:, num_signals+1:end);
% P_tf_music = zeros(1, length(theta));
% 
% for ii = 1:length(theta)
%     a = exp(-1j*(0:num_elements-1)'*(2*pi/wavelength)*d*sin(deg2rad(theta(ii))));
%     noiseMatrix = noiseVec * noiseVec';
%     P_tf_music(ii) = 1 / (a' * noiseMatrix * a);
% end
% 
% P_tf_music = 10*log10(abs(P_tf_music)/max(abs(P_tf_music)));
% [~, max_idx] = max(P_tf_music);
% estimated_angle = theta(max_idx);

function [estimated_angles, P_tf_music] = SQSTFT_trackmultiLRmethod_Music(data_matrix1, num_elements, d, wavelength, theta, fs, lowFreq, highFreq, alpha_param, tDS, h, Dh, r_filter, gamma_filter, delta_limit, min_length, LOWER_PRCTILE_LIMIT)

num_channels = size(data_matrix1, 1);
c = 3e8;
for ch = 1:num_channels
    signal = data_matrix1(ch, :).';
    [filtered_signal, s_hat] = ASAP_Hankel_1D(signal, r_filter, gamma_filter);
    data_matrix1(ch, 1:length(filtered_signal)) = filtered_signal.';
end

% 时频变换
S_matrix_filt = cell(1, num_channels);
F_matrix_filt = cell(1, num_channels);
T_matrix_filt = cell(1, num_channels);
for ch = 1:num_channels
    signal = data_matrix1(ch, :).';
    [~, tfrtic, tfrsq, tfrsqtic] = sqSTFT(signal, lowFreq, highFreq, alpha_param, tDS, h', Dh');
    S_matrix_filt{ch} = tfrsq;
    F_matrix_filt{ch} = tfrsqtic * fs;
    T_matrix_filt{ch} = (1:length(signal)) / fs;
end

% 轨迹提取与融合
[tf_pts, aligned] = tracks_multiLRmethod(S_matrix_filt, F_matrix_filt, fs, delta_limit, min_length, LOWER_PRCTILE_LIMIT, d, c);
valid_rows = all(aligned ~= 0, 2);
aligned    = aligned(valid_rows, :);
tf_pts     = tf_pts(valid_rows, :);

% 加权因子计算
E = mean(abs(aligned).^2, 2);
E = normalize(E, 'range');
dT = abs(diff(tf_pts(:,1)));
smoothness = 1 ./ (1 + dT);
smoothness = [smoothness; smoothness(end)];
SM = normalize(smoothness, 'range');

theta_scan = -90:0.5:90;
R0 = (aligned' * aligned) / size(aligned,1);
a_fun = @(theta) exp(-1j*2*pi*d*sin(theta*pi/180)*(0:size(aligned,2)-1)'/c);
P_music = arrayfun(@(th) 1./(a_fun(th)' * (R0 - a_fun(th)*a_fun(th)') * a_fun(th)), theta_scan);
[~, idx_max] = max(real(P_music));
theta0 = theta_scan(idx_max);
v0 = a_fun(theta0);
v0 = v0 / norm(v0);
sim = abs(aligned * v0);
ADAPT = normalize(sim, 'range');

% 融合权重
weights = 0.5*E + 0.25*SM + 0.25*ADAPT;
weights = weights.^1.5;

% 加权协方差矩阵与特征分解，信号数改为2
Xw = aligned .* sqrt(weights);
R  = (Xw' * Xw) / sum(weights);
[eigenVec, eigenVal] = eig(R);
[eigenVal, idx] = sort(diag(eigenVal), 'descend');
eigenVec = eigenVec(:, idx);
num_signals = 2;  
signalVec = eigenVec(:, 1:num_signals);
noiseVec = eigenVec(:, num_signals+1:end);

% 计算空间谱
P_tf_music = zeros(1, length(theta));
for ii = 1:length(theta)
    steering = exp(1j*(0:num_elements-1)'*(2*pi/wavelength)*d*sin(deg2rad(theta(ii))));
    noiseMatrix = noiseVec * noiseVec';
    P_tf_music(ii) = 1 / (steering' * noiseMatrix * steering);
end

P_tf_music = 10*log10(abs(P_tf_music) / max(abs(P_tf_music)));
% 找到两个峰值对应的角度
[~, sorted_idx] = sort(P_tf_music, 'descend');
top_idx = sorted_idx(1:num_signals);
estimated_angles = theta(top_idx);
end
