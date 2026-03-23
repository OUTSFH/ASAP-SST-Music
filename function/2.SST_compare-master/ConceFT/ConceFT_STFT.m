function [tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_STFT(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT, Smooth, Hemi) 
   
    N = length(x);
    [h, Dh, ~] = hermf(WinLen, dim, supp); 
    fprintf(['Run ordinary STFT-SST (Smooth = ',num2str(Smooth),', Hemi = ',num2str(0),')\n']);
    [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, h(:,1), Dh(:,1), Smooth, 0);
    ConceFT = [];

    if MT > 1
    ConceFT = zeros(size(tfrsq));
    fprintf(['STFT-ConceFT total (Smooth = ',num2str(Smooth),', Hemi = ',num2str(Hemi),'): ',num2str(MT),'; now:     ']);
        for ii = 1:MT
            fprintf('\b\b\b\b');	tmp = sprintf('%4d',ii); fprintf([tmp]);
            rv = randn(dim, 1); 
            rv = rv ./ norm(rv);
            rh = h * rv;   
            rDh = Dh * rv;
            [~, ~, tfrsq_temp, tfrsqtic_temp] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, rh, rDh, Smooth, Hemi);
            ConceFT = ConceFT + tfrsq_temp;
        end

        ConceFT = ConceFT ./ MT;
        fprintf('\n');
        tfrsq = ConceFT;
        
    end
end
