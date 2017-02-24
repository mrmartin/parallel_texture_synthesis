%upsample the coordinates S for a texture of size m at level l
function S_new = upsample_s (S_old, m, l)
    S_new = zeros(2^l,2^l,2);
    S_new(1:2:end,1:2:end,:)=S_old;%copy
    S_new(2:2:end,1:2:end,1)=S_old(:,:,1)+m/2^l;%one down
    S_new(2:2:end,1:2:end,2)=S_old(:,:,2);
    S_new(1:2:end,2:2:end,1)=S_old(:,:,1);%one right
    S_new(1:2:end,2:2:end,2)=S_old(:,:,2)+m/2^l;%one right
    S_new(2:2:end,2:2:end,:)=S_old+m/2^l;%down and right
    S_new = mod(S_new-1,m)+1;
end