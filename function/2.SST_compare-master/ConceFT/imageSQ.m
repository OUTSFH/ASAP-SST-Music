function imageSQ(t, ytic, M, Qv)
    if nargin < 4, Qv = 1; 
    end
    fz = 20;
    Q = M(:);
    
    fprintf('\n\t\t\t ** Smallest value = %g\n', min(M(:)));
    fprintf('\t\t\t ** Maximal value = %g\n', max(M(:)));
    fprintf('\t\t\t ** 99.8%% value = %g\n', quantile(M(:), 0.998));
    fprintf('\t\t\t ** 0.2%% value = %g\n', quantile(M(:), 0.002));
    
    q = quantile(Q, Qv);
    M(M > q) = q;
    M = M ./ q;
    m = quantile(Q, 0.002);
    M(M < m) = m;
    imagesc(t, ytic, M);
    axis xy;
    set(gca, 'fontsize', fz);
end