% one channel terrain version!
%texture synthesis as per "Parallel Controllable Texture Synthesis"
% by Lefebvre and Hoppe
% Section 3.1 - basic scheme with gaussian stack
close all

%% load real terrain and prepare gaussian stack
% [im_in, R]=sdtsdemread('~/Datasets/USGS_terrain/1106CEL0.DDF');
%remove all NaN regions
% im_in=im_in(33:1314,43:1102);
clear
from_disk=false;

if(from_disk)
    clear
    load mountains_64_pca.mat
else
    f=fopen('~/Desktop/exported terrains/real_mountains.raw');
    raw=fread(f,1025*1025,'uint16');
    fclose(f);
    extremes=csvread('~/Desktop/exported terrains/canyons.csv');
    low=extremes(1);
    high=extremes(2);
    raw=reshape(raw,1025,1025);
    Z=low+(raw*(high-low)/max(raw(:)));
    
%     [Z, R]=sdtsdemread('exported terrains/DEM_1/1106CEL0.DDF');
%     Z=Z(61:1388,42:1103);%keep non-NaN region only
    im_in=Z(round(linspace(1,size(Z,1)-1,64)),round(linspace(1,size(Z,2)-1,64)));
    
    %get a small subsampled region
    %im_in=Z(1:16:1024,1:16:1024);%im_in(1:16:1024,1:16:1024);
    %remove edges from the rendered texture, because they're unrealistic
    %im_in=Z(round(linspace(50,970,128)),round(linspace(50,970,128)));

    %repeat texture including flipped versions
    %im_in=[im_in im_in(:,end:-1:1); im_in(end:-1:1,:) im_in(end:-1:1,end:-1:1)];

    if(round(log2(size(im_in,1)))~=log2(size(im_in,1)) || round(log2(size(im_in,2)))~=log2(size(im_in,1)))
        disp('error, input texture must be square, of size 2^l')
    end

    %imagesc(repmat(im_in,3,3,1))
    m=size(im_in,1);

    input_levels = log2(m);

    %precompute neighbours into convenient form at each level > 2
    neighbourhood = 5;

    %extend over edges by flipping. Done instead of repeating because it's not toroidal
    im_extended_flipped=[im_in(end:-1:1,end:-1:1) im_in(end:-1:1,:) im_in(end:-1:1,end:-1:1); im_in(:,end:-1:1) im_in im_in(:,end:-1:1); im_in(end:-1:1,end:-1:1) im_in(end:-1:1,:) im_in(end:-1:1,end:-1:1)];

    gaussian_stack = cell(input_levels,1);
    pca_coeffs = cell(input_levels,1);
    means = cell(input_levels,1);
    for l_out=3:input_levels
        range=m/2^l_out*[(neighbourhood-1)/-2:(neighbourhood-1)/2];
        [X,Y]=meshgrid(range, range);
        %remove self
        X=X(setdiff(1:end,neighbourhood^2/2+0.5));
        Y=Y(setdiff(1:end,neighbourhood^2/2+0.5));

        %perform blur on all levels except last
        if(l_out<input_levels)
            im_gaussian_blur = imgaussfilt(im_extended_flipped,2^(input_levels-l_out-3));
            im_gaussian_blur = im_gaussian_blur(m+[1:m],m+[1:m]);
        else
            im_gaussian_blur = im_in;
        end
               
        gaussian_stack{l_out}=zeros(m^2,neighbourhood^2-1);
        for x=1:m
            for y=1:m
                gaussian_stack{l_out}(sub2ind([m,m],x,y),:)=im_gaussian_blur(sub2ind([m,m],mod(X+x-1,m)+1,mod(Y+y-1,m)+1));
                %gaussian_stack{l_out}(sub2ind([m,m],x,y),:)=gaussian_stack{l_out}(sub2ind([m,m],x,y),:)-im_gaussian_blur(x,y);%normalize by own height
                if(max((mod(X+x-1,m)+1)~=(X+x)) || max((mod(Y+y-1,m)+1)~=(Y+y))) % discard those which reach outside the boundary
                    gaussian_stack{l_out}(sub2ind([m,m],x,y),1)=Inf;
%                 else % also remove those which cross the middle, because those are unrealistic flip artefacts
%                     vx=floor((X+x)/(m/2));
%                     vy=floor((Y+y)/(m/2));
%                     if(all(vx==vx(1)) && all(vy==vy(1)))
%                         %ok
%                     else
%                         gaussian_stack{l_out}(sub2ind([m,m],x,y),1)=Inf;
%                     end
                end
            end
        end
        usable=gaussian_stack{l_out}(:,1)~=Inf;
        pca_coeffs{l_out}=pca(gaussian_stack{l_out}(usable,:));
        pca_coeffs{l_out}=pca_coeffs{l_out}(:,1:10);%use first five only
        means{l_out}=mean(gaussian_stack{l_out}(usable,:));
        gaussian_stack{l_out}=(gaussian_stack{l_out}-repmat(means{l_out},size(gaussian_stack{l_out},1),1))*pca_coeffs{l_out};
    end
    
    save mountains_64_pca
end

%% view input

subplot(1,3,1)
[X Y]=meshgrid(1:size(im_in));
surface(X,Y,im_in,'EdgeColor','none','FaceColor','texturemap')
view(3)
title('input')

%%
%the input and output need not be the same size
%the scale needs to be the same!
%The last output level will be the finest input level
output_size=128 ;%power of two
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

corrections = 3;
%jitter parameter at each level
r=repmat(0.2,50,1);%zeros(1,6);%[1 1 1 0.3 0 0];%

for l_out=max(1,input_levels-output_levels+1):input_levels
    S=upsample_s(S,m,l_out,false);%not toroidal
    S=jitter_s(S, m, r(l_out), l_out, random_seed, false);%not toroidal
    if(l_out>2)%perform corrections on higher levels only
        for c=1:corrections
            S=correct_s_pca_strange_results(S,gaussian_stack{l_out},im_in(:),pca_coeffs{l_out},means{l_out},m,1);%requires edit for more than 1 image channels
        end
        [connected, flow_targets]=watershed_basins(output);

        %strategy B
        if(false & l_out>3)
            depression_pixel=find_inflow_region(flow_targets, size(S,1));%each region id is the index of the lowest pixel
            disp(['there are already ' num2str(length(depression_pixel)) ' pixels which need fixing'])
            counter=0;
            while ~isempty(depression_pixel)
                S=correct_s_raise_single_pixel(S,gaussian_stack{l_out},im_in,m,depression_pixel(1));
                output=reshape(im_in(sub2ind([m,m],S(:,:,1),S(:,:,2))),size(S,1),size(S,2));
                smaller_and_equal=true;
                flow_targets=find(find_local_minima(output,smaller_and_equal));
                depression_pixel=find_inflow_region(flow_targets, size(S,1));
                counter=counter+1;
            end
            disp(['it turned out to require ' num2str(counter) ' fixes'])
%             [connected, flow_targets]=watershed_basins(output);
%             imagesc(connected)
%             cmap=polcmap;
%             colormap(h,cmap)
        end
    end
    h=subplot(1,3,2);
    output=reshape(im_in(sub2ind([m,m],S(:,:,1),S(:,:,2))),size(S,1),size(S,2));%index of current output in the input
    %[X Y]=meshgrid(1:size(S,1));
    %surface(X,Y,output,'EdgeColor','none','FaceColor','texturemap')
    %view(3)
    imagesc(output)
    %imagesc(reshape(im_in(sub2ind([m,m],S(:,:,1),S(:,:,2))),size(S,1),size(S,2)))
    title(['output level ' num2str(l_out)])
    drawnow
    pause
end