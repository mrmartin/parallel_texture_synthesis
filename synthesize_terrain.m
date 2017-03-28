% one channel terrain version!
%texture synthesis as per "Parallel Controllable Texture Synthesis"
% by Lefebvre and Hoppe
% Section 3.1 - basic scheme with gaussian stack
close all

%% load real terrain and prepare gaussian stack
[im_in, R]=sdtsdemread('~/Datasets/USGS_terrain/1106CEL0.DDF');
%remove all NaN regions
im_in=im_in(33:1314,43:1102);

% f=fopen('exported terrains/canyons.raw');
% raw=fread(f,1025*1025,'uint16');
% fclose(f);
% extremes=csvread('~/Downloads/exported terrains/canyons.csv');
% low=extremes(1);
% high=extremes(2);
% raw=reshape(raw,1025,1025);
% Z=low+(raw*(high-low)/max(raw));
% min(Z)
% max(Z)
% Z=low+(raw*(high-low)/max(raw(:)));

%get a small subsampled region
im_in=im_in(1:16:1024,1:16:1024);

%repeat texture including flipped versions
%im_in=[im_in im_in(:,end:-1:1); im_in(end:-1:1,:) im_in(end:-1:1,end:-1:1)];

if(round(log2(size(im_in,1)))~=log2(size(im_in,1)) || round(log2(size(im_in,2)))~=log2(size(im_in,1)))
    disp('error, input texture must be square, of size 2^l')
end

%imagesc(repmat(im_in,3,3,1))
m=size(im_in,1);

levels = log2(m);

%precompute neighbours into convenient form at each level > 2
neighbourhood = 5;

%extend over edges by flipping. Done instead of repeating because it's not toroidal
im_extended_flipped=[im_in(end:-1:1,end:-1:1) im_in(end:-1:1,:) im_in(end:-1:1,end:-1:1); im_in(:,end:-1:1) im_in im_in(:,end:-1:1); im_in(end:-1:1,end:-1:1) im_in(end:-1:1,:) im_in(end:-1:1,end:-1:1)];

Nexemplar = cell(levels,1);
for l=2:levels
    range=m/2^l*[(neighbourhood-1)/-2:(neighbourhood-1)/2];
    [X,Y]=meshgrid(range, range);
    %remove self
    X=X(setdiff(1:end,neighbourhood^2/2+0.5));
    Y=Y(setdiff(1:end,neighbourhood^2/2+0.5));
    
    %perform blur on all levels except last
    if(l<levels)
        im_gaussian_blur = imgaussfilt(im_extended_flipped,2^(levels-l-3));
        im_gaussian_blur = im_gaussian_blur(m+[1:m],m+[1:m]);
    else
        im_gaussian_blur = im_in;
    end
    %subplot(1,4,l-2)
    %imagesc(im_gaussian_blur)
    
    Nexemplar{l}=zeros(m^2,neighbourhood^2-1);
    for x=1:m
        for y=1:m
            Nexemplar{l}(sub2ind([m,m],x,y),:)=im_gaussian_blur(sub2ind([m,m],mod(X+x-1,m)+1,mod(Y+y-1,m)+1));
        end
    end
end
%%
% f2=figure
% f1=figure
% pause
for myr = 1%:9
    %S has three dimensions:
    % first and second are output x, y
    % third is (x, y) in input
    S=reshape([1 1],[1 1 2]);

    corrections = 3;
    %jitter parameter at each level
    r=[0 0 0.02 0.9 0.9 0.9];%repmat(0.9,levels,1);%zeros(1,6);%[1 1 1 0.3 0 0];%

    for l=1:levels
%         figure(f1)
        clf
        subplot(1,3,1)
        [X Y]=meshgrid(1:size(im_in));
        surface(X,Y,im_in,'EdgeColor','none','FaceColor','texturemap')
        view(3)
        title('input')
        S=upsample_s(S,m,l,false);%not toroidal
        S=jitter_s(S, m, r(l), l, myr,false);%not toroidal
        %subplot(1,3,1)
        %imagesc(reshape(im_in(sub2ind([m,m],S(:,:,1),S(:,:,2))),size(S,1),size(S,2)))
        %title('upsampled and jittered')
        if(l>1)
            for c=1:corrections
                S=correct_s_one_channel(S,Nexemplar{l},im_in,m);
            end
            subplot(1,3,2)
            output=reshape(im_in(sub2ind([m,m],S(:,:,1),S(:,:,2))),size(S,1),size(S,2));
            [X Y]=meshgrid(1:size(S,1));
            surface(X,Y,output,'EdgeColor','none','FaceColor','texturemap')
            view(3)
            %imagesc(reshape(im_in(sub2ind([m,m],S(:,:,1),S(:,:,2))),size(S,1),size(S,2)))
            title('output')
            drawnow
            %watershed basin computations
            [connected, flow_targets]=watershed_basins(output);

            % strategy A - unfinished
    %         region_ids = find_inflow_region(flow_targets, size(S,1));
    %         while ~isempty(region_ids)
    %             region_ids
    %             S=correct_outflow_region(S,Nexemplar{l},im_in,m,region_ids(1),connected);
    %             output=reshape(im_in(sub2ind([m,m],S(:,:,1),S(:,:,2))),size(S,1),size(S,2));
    %             [connected, flow_targets]=watershed_basins(output);
    %             region_ids = find_inflow_region(flow_targets, size(S,1));
    %         end

            %strategy B
            if(l>3)
                depression_pixel=find_inflow_region(flow_targets, size(S,1));%each region id is the index of the lowest pixel
                disp(['there are already ' num2str(length(depression_pixel)) ' pixels which need fixing'])
                counter=0;
                while ~isempty(depression_pixel)
                    S=correct_s_raise_single_pixel(S,Nexemplar{l},im_in,m,depression_pixel(1));
                    output=reshape(im_in(sub2ind([m,m],S(:,:,1),S(:,:,2))),size(S,1),size(S,2));
                    %[connected, flow_targets]=watershed_basins(output);
                    flow_targets=find(find_local_minima(output));
                    depression_pixel=find_inflow_region(flow_targets, size(S,1));
                    counter=counter+1;
                end
                disp(['it turned out to require ' num2str(counter) ' fixes'])

                [connected, flow_targets]=watershed_basins(output);
                h=subplot(1,3,3);
                imagesc(connected)
                cmap=polcmap;
                colormap(h,cmap)
                title('connected')
                pause(0.5)
            end
        end
    end
%     figure(f2)
%     subplot(3,3,myr)
%     imagesc(output)
%     title(['seed=' num2str(myr)])
end