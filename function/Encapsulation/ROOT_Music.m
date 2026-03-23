function [estimated_angle, P_root_music] = ROOT_Music(data, theta, num_elements, d, wavelength)

M = num_elements;         
N = size(data, 2);        
R = (data * data') / N;  
[eigenVec, eigenVal] = eig(R);
[eigenVal, idx] = sort(diag(eigenVal), 'descend');
eigenVec = eigenVec(:, idx);
num_signals = 1;
noiseVec = eigenVec(:, num_signals+1:end);  
En = noiseVec * noiseVec';
b = zeros(2*(M-1), 1);

for i = -(M-1):(M-1)
    b(i + M) = sum(diag(En, i));  
end

b = flipud(b);  
rts = roots(b);
distance = 1 - abs(rts);  
[~, sort_idx] = sort(abs(distance));
valid_rts = rts(sort_idx);

DOA = asin((angle(valid_rts(end)) * wavelength) / (2*pi*d)) * 180/pi;
estimated_angle = DOA;
P_root_music = zeros(size(theta));

for ii = 1:length(theta)
    a_scan = exp(1j * (0:M-1)' * (2*pi/wavelength) * d * sin(deg2rad(theta(ii))));
    P_root_music(ii) = 1 / (a_scan' * En * a_scan);
end

P_root_music = 10 * log10(abs(P_root_music) / max(abs(P_root_music)));

end
