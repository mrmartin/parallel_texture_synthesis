%upsample the coordinates S for a texture of size m at level l
function S_new = upsample_s (S_old, m, level_scale, toroidal)
    next_level_size=log2(size(S_old,1))+1;
    S_new = zeros(2^next_level_size,2^next_level_size,2);
    S_new(1:2:end,1:2:end,:)=S_old;%copy
    S_new(2:2:end,1:2:end,1)=S_old(:,:,1)+m/2^level_scale;%one down
    S_new(2:2:end,1:2:end,2)=S_old(:,:,2);
    S_new(1:2:end,2:2:end,1)=S_old(:,:,1);%one right
    S_new(1:2:end,2:2:end,2)=S_old(:,:,2)+m/2^level_scale;%one right
    S_new(2:2:end,2:2:end,:)=S_old+m/2^level_scale;%down and right
    if(toroidal)
        S_new = mod(S_new-1,m)+1;
    else
        S_new = abs(m-1-abs(m-S_new))+1;
    end
end