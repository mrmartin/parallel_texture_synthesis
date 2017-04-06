%texture synthesis as per "Parallel Controllable Texture Synthesis"
% by Lefebvre and Hoppe
% Section 3.1 - basic scheme with gaussian stack

%preprocess_terrain_flow

load('mountains_with_precomputed_flow_vectors')
im_in(:,:,1)=(Z-mean(Z(:)))/max(Z(:));
im_in(:,:,2)=reshape(directions(:,1),128,128)/max(abs(directions(:,1)));
im_in(:,:,3)=reshape(directions(:,2),128,128)/max(abs(directions(:,2)));

neighbourhood = 5;
%add the neighbours of depressions to unusable pixel list
badX=repmat(depressions(:,1),1,neighbourhood*2+1)+repmat(-neighbourhood:neighbourhood,size(depressions,1),1);
badY=repmat(depressions(:,2),1,neighbourhood*2+1)+repmat(-neighbourhood:neighbourhood,size(depressions,1),1);
depression_neighbours=[];
for i=1:size(depressions,1)
    [X Y]=meshgrid(badX(i,:),badY(i,:));
    depression_neighbours=[depression_neighbours;X(:) Y(:)];
end
depressions=unique(depression_neighbours,'rows');

if(round(log2(size(im_in,1)))~=log2(size(im_in,1)) || round(log2(size(im_in,2)))~=log2(size(im_in,1)))
    disp('error, input texture must be square, of size 2^l')
end

%is the input toroidal?
toroidal = false;

%perform pixel comparisons in LAB space
colorTransform = makecform('srgb2lab');
im_orig = im_in;
%im_in = applycform(im_in, colorTransform);

%imagesc(repmat(im_in,3,3,1))
m=size(im_in,1);
pixels_in=reshape(im_in,m^2,3);

input_levels = log2(size(im_in,1));

%precompute neighbours into convenient form at each level > 2
Nexemplar = cell(input_levels,1);
for l=3:input_levels
    range=m/2^l*[(neighbourhood-1)/-2:(neighbourhood-1)/2];
    [X,Y]=meshgrid(range, range);
    %remove self
    X=X(setdiff(1:end,neighbourhood^2/2+0.5));
    Y=Y(setdiff(1:end,neighbourhood^2/2+0.5));
    
    %only perform blur on lower levels
    if(l<input_levels)
        im_gaussian_blur = imgaussfilt(repmat(im_in,3,3),2^(input_levels-l-1));
        im_gaussian_blur = im_gaussian_blur(m+[1:m],m+[1:m],:);
    else
        im_gaussian_blur = im_in;
    end
    
    pixels_gaussian_blur=reshape(im_gaussian_blur,m^2,3);

    %prepare the neighbourhood of each pixel on a line
    Nexemplar{l}=zeros(m^2,neighbourhood^2-1,3);
    for x=1:m
        for y=1:m
            Nexemplar{l}(sub2ind([m,m],x,y),:,:)=reshape(pixels_gaussian_blur(sub2ind([m,m],mod(X+x-1,m)+1,mod(Y+y-1,m)+1),:),[1 neighbourhood^2-1 3]);
            if(max((mod(X+x-1,m)+1)~=(X+x)) || max((mod(Y+y-1,m)+1)~=(Y+y)))
                Nexemplar{l}(sub2ind([m,m],x,y),1)=Inf;%discard those which reach outside the boundary
            end
            if(any(x==depressions(:,1) & y==depressions(:,2)))
                Nexemplar{l}(sub2ind([m,m],x,y),1)=Inf;%discard depressions
            end
        end
    end
end

%the input and output need not be the same size
%the scale needs to be the same!
%The last output level will be the finest input level
output_size=256;%power of two
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
S=reshape(1+floor(rand(level_difference^2,2)*size(im_in,1)),[level_difference level_difference 2]);% <--random / deterministic --> S=reshape([16 16],[1 1 2]);
imagesc(reshape(pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3))

corrections = 4;
%jitter parameter at each level
r=[1 1 1 0.3 0 0 0];%repmat(0.4,levels,1);

pixels_in=reshape(im_orig,m^2,3);
close
tic
for l_out=max(1,input_levels-output_levels+1):input_levels
    S=upsample_s(S,m,l_out,toroidal);
    S=jitter_s(S, m, r(l_out), l_out, random_seed, toroidal);
    %pause(0.5)
    if(l_out>2)
        for c=1:corrections
            S=correct_s(S,Nexemplar{l_out},pixels_in,m);
        end
    end
    output=reshape(pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3);
    subplot(1,2,1)
    imagesc(output(end:-1:1,:,1))
    [X Y]=meshgrid(1:size(output,1),1:size(output,2));
    title(['ouptut at level ' num2str(l_out)])
    subplot(1,2,2)
    quiver(X,Y,output(:,:,2),output(:,:,3))
    axis([1 size(output,1) 1 size(output,2)])
    pause(0.3)%better than drawnow
end
toc
subplot(1,2,1)
title('final')
%imwrite(reshape(pixels_in(sub2ind([m,m],S(:,:,1),S(:,:,2)),:),size(S,1),size(S,2),3),'output_texture.png')