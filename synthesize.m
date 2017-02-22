%texture synthesis as per "Parallel Controllable Texture Synthesis"
% by Lefebvre and Hoppe
% Section 3.1 - basic scheme

%input toroidal texture of 2^l
im_in = imread('input_texture.png');

if(round(log2(size(im_in,1)))~=log2(size(im_in,1)) || round(log2(size(im_in,2)))~=log2(size(im_in,1)))
    disp('error, input texture must be square, of size 2^l')
end

%imagesc(repmat(im_in,3,3,1))
m=size(im_in,1);
pixels_in=reshape(im_in,m^2,3);

levels = log2(size(im_in,1));

%precompute neighbours into convenient form !!at each level > 2!!
neighbourhood = 5;

Nexemplar = cell(levels,1);
for l=3:levels
    range=m/2^l*[(neighbourhood-1)/-2:(neighbourhood-1)/2];
    [X,Y]=meshgrid(range, range);
    %remove self
    X=X(setdiff(1:end,neighbourhood^2/2-0.5));
    Y=Y(setdiff(1:end,neighbourhood^2/2-0.5));
    
    im_gaussian_blur = imgaussfilt(repmat(im_in,3,3),2^(levels-l));
    im_gaussian_blur = im_gaussian_blur(m+[1:m],m+[1:m],:);
    
    pixels_gaussian_blur=reshape(im_gaussian_blur,m^2,3);

    Nexemplar{l}=zeros(m^2,neighbourhood^2-1,3);
    for x=1:m
        for y=1:m
            Nexemplar{l}(sub2ind([m,m],x,y),:,:)=reshape(pixels_gaussian_blur(sub2ind([m,m],mod(X+x-1,m)+1,mod(Y+y-1,m)+1),:),[1 neighbourhood^2-1 3]);
        end
    end
end

%S has three dimensions:
% first and second are output x, y
% third is (x, y) in input
S=reshape([1 1],[1 1 2]);
imagesc(reshape(pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3))

corrections = 3;
%jitter parameter at each level
r=[0 0 0 0.05 0 0];%repmat(0.4,levels,1);

for l=1:levels
    S=upsample_s(S,m,l);
    S=jitter_s(S, m, r(l), l, 7);
    imagesc(reshape(pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3))
    title('upsampled and jittered')
    pause(0.5)
    if(l>2)
        for c=1:corrections
            S=correct_s(S,Nexemplar{l},pixels_in,m);
        end
        imagesc(reshape(pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3))
        title('corrected')
        pause(0.5)
    end
end
%close