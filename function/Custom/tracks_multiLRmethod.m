% 多阵元协同轨迹提取
function [pts, aligned] = tracks_multiLRmethod(S_cell,F_cell,fs,delta,minlen,lim,d,c)
    
    M = numel(S_cell); 
    all_pts=[];

    for ch=1:M
        tfmag = abs(S_cell{ch}).';
        [trks,~] = tracks_LRmethod(tfmag,fs,delta,minlen,lim);
        for k=1:numel(trks), all_pts=[all_pts;trks{k}]; 
        end
    end

    pts = unique(all_pts,'rows'); 
    P=size(pts,1);
    raw = zeros(P,M);
    
    for ch=1:M
        mat = S_cell{ch}; 
        tmp=zeros(P,1);
        for p=1:P
            t=pts(p,1); 
            f=pts(p,2);
            if t<=size(mat,2)&&f<=size(mat,1), tmp(p)=mat(f,t); 
            end
        end
        raw(:,ch)=tmp;
    end

    pos = (0:M-1)*d; 
    aligned=zeros(size(raw));
    
    for p=1:P
        fk = F_cell{1}(pts(p,2));
        for m=1:M
            tau = pos(m)/c; 
            aligned(p,m)=raw(p,m)*exp(1j*2*pi*fk*tau);
        end
    end
end
