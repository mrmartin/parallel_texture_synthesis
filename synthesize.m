%texture synthesis as per "Parallel Controllable Texture Synthesis"
% by Lefebvre and Hoppe
% Section 3.1 - basic scheme with gaussian stack

%input toroidal texture of 2^l
im_in = imread('input_texture.png');

if(round(log2(size(im_in,1)))~=log2(size(im_in,1)) || round(log2(size(im_in,2)))~=log2(size(im_in,1)))
    disp('error, input texture must be square, of size 2^l')
end

%is the input toroidal?
toroidal = false;

%precompute neighbours into convenient form at each level > 2
neighbourhood = 7;

%precompute similar candidate pixels C
k_candidates = 10;

%perform pixel comparisons in LAB space
colorTransform = makecform('srgb2lab');
im_orig = im_in;
im_in = applycform(im_in, colorTransform);

%[Nexemplar, candidates, pixels_in, m, input_levels] = preprocess_texture(im_in, neighbourhood, k_candidates);

%the input and output need not be the same size
%the scale needs to be the same!
%The last output level will be the finest input level
output_size=128;%power of two
output_levels=log2(output_size);

%S has three dimensions:
% first and second are output x, y
% third is (x, y) in input
random_seed = 2;
rng(random_seed)

%if output_levels is less than input_levels, S_1 will be larger than 1
%if output_levels equals input_levels, S_1 will be [1 x 1 x 2]
%if otput_levels is more than input levels S_1 will be [1 x 1 x 2], but at a higher scale
level_difference=max(1,2^(output_levels-input_levels));
S=reshape(1+floor(rand(level_difference^2,2)*size(im_in,1)),[level_difference level_difference 2])% <--random / deterministic --> S=reshape([16 16],[1 1 2]);
imagesc(reshape(pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3))

corrections = [0 0 9 6 3 1];%tune the number of correction passes per level, for speed and quality
%jitter parameter at each level
r=[1 1 1 0.3 0 0 0];%repmat(0.4,levels,1);

close
hFig = figure(1);
set(hFig, 'Position', [0 205 1150 465])
tic
for l_out=max(1,input_levels-output_levels+1):input_levels
    toc
    disp(['starting to perform synthesis at level ' num2str(l_out)])
    tic
    S=upsample_s(S,m,l_out,toroidal);
    S=jitter_s(S, m, r(l_out), l_out, random_seed, toroidal);
    subplot(1,2,1)
    imagesc(reshape(pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3))
    title('upsampled and jittered')
    %pause(0.5)
    if(l_out>2)
        for c=1:corrections(l_out)
            S=correct_s(S,Nexemplar{l_out},pixels_in,m,candidates{l_out},neighbourhood);
        end
        subplot(1,2,2)
        imagesc(reshape(pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3))
        title('corrected')
        pause(0.5)
    end
end
toc
%close
orig_pixels_in=reshape(im_orig,m^2,3);
imagesc(reshape(orig_pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3))
title('final')
%imwrite(reshape(pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3),'output_texture.png')