%jitter the coordinates S for a texture of size m with randomness r, level l,
%and hash parameter h
function S_new = jitter_s(S_old, m, r, l, h)
    rng(h);
    S_new = S_old+floor(r*(rand([size(S_old,1),size(S_old,2),2])*(m/2^(l-1))-(m/2^l)));
    S_new = mod(S_new-1,m)+1;
end