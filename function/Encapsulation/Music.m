function [estimated_angle, P_music] = Music(data, theta, num_elements, d, wavelength)

M = num_elements;         
N = size(data, 2);        
R = (data * data') / N;  
[eigenVec, eigenVal] = eig(R);
[eigenVal, idx] = sort(diag(eigenVal), 'descend');
eigenVec = eigenVec(:, idx);
num_signals = 1;
noiseVec = eigenVec(:, num_signals+1:end);  
P_music = zeros(1, length(theta));

for ii = 1:length(theta)
    angle_rad = deg2rad(theta(ii));
    a_scan = exp(1j * (0:M-1)' * (2*pi/wavelength) * d * sin(angle_rad));
    noiseMatrix = noiseVec * noiseVec';
    P_music(ii) = 1 / (a_scan' * noiseMatrix * a_scan);  
end

P_music = 10 * log10(abs(P_music) / max(abs(P_music)));
[~, max_idx] = max(P_music);
estimated_angle = theta(max_idx);

end
