%% load real terrain, mountains
f=fopen('~/Desktop/exported terrains/real_mountains.raw');
raw=fread(f,1025*1025,'uint16');
fclose(f);
extremes=csvread('~/Desktop/exported terrains/canyons.csv');
low=extremes(1);
high=extremes(2);
raw=reshape(raw,1025,1025);
Z=low+(raw*(high-low)/max(raw(:)));
Z=Z(1:2:1024,1:2:1024);

%% load another terrain, yosemite
% [Z, R]=sdtsdemread('exported terrains/DEM_1/1106CEL0.DDF');
% Z=Z(61:1388,42:1103);
% Z=Z(round(linspace(1,size(Z,1),128)),round(linspace(1,size(Z,2),128)));

%% find best blur level
% for i=linspace(0.5,5,1000)
%     tmp=find_local_minima(imgaussfilt(Z,i),true,0);disp([num2str(i) ' - ' num2str(length(find(tmp(2:end-1,2:end-1))))])
% end

%% smooth the image
Z_orig=Z;
Z=imgaussfilt(Z,3.5);

%% Brute force rain flow simulation
close
%for each point, flow from neighbouring points
rain=ones(size(Z));
flow=rain;
inflow=zeros(size(Z));
for i=1:100
    depressions=[];
    flow_tmp=zeros(size(Z));
    for x=1:size(Z,1)
        for y=1:size(Z,2)
            if(~isnan(Z(x,y)))
                %find lower neighbours
                neighbors=[];
                if(x>1)
                    neighbors=[neighbors;x-1 y];
                    if(y>1)
                        neighbors=[neighbors;x-1 y-1];
                    end
                end
                if(x<size(Z,1))
                    neighbors=[neighbors;x+1 y];
                    if(y<size(Z,2))
                        neighbors=[neighbors;x+1 y+1];
                    end
                end
                if(y>1)
                    neighbors=[neighbors;x y-1];
                    if(x<size(Z,1))
                        neighbors=[neighbors;x+1 y-1];
                    end
                end
                if(y<size(Z,2))
                    neighbors=[neighbors;x y+1];
                    if(x>1)
                        neighbors=[neighbors;x-1 y+1];
                    end
                end
                neighbors=sub2ind(size(Z),neighbors(:,1),neighbors(:,2));
                if(sum(max(0,Z(x,y)-Z(neighbors)))==0)
                    if(x==1 | y==1 | x==size(Z,1) | y==size(Z,2))
                        %ok, on the boundary
                    else
                        depressions=[depressions;x y];
                        %disp(['there are no lower neighbours to [' num2str(x) ', ' num2str(y) ']']);
                        %1)it's a lake
                        % accumulate the inflow over the whole area, and accumulate outflow into single pixel
                        %2)it's an error
                        % do not use this pixel in synthesis
                    end
                end
                distribution=max(0,Z(x,y)-Z(neighbors))/sum(max(0,Z(x,y)-Z(neighbors)));
                distribution(isnan(distribution))=0;
                if(max(isnan(distribution)))
                    x
                    y
                end
                flow_tmp(neighbors)=flow_tmp(neighbors)+flow(x,y).*distribution;
            end
        end
    end
    flow=flow_tmp;
    imagesc(-1*flow)
    pause(0.2)
    inflow=inflow+flow;
end


%% terrain to pixel flow strength and direction
imsize=size(Z,1);
close

subplot(1,3,1)
imagesc(Z(end:-1:1,:))
title('height, showing forbidden pixels')
%plot pixels which don't have a neighbour with more water
% [x y]=ind2sub([imsize imsize],find(find_local_minima(-1*inflow,true,0.03)));

% find pixels from which water does not flow out of the terrain
[colors, outflow]=watershed_basins(Z);%perform with a smoothed image
inflow_regions=find_inflow_region(outflow,imsize);
[x y]=ind2sub([imsize imsize],find(ismember(colors,inflow_regions)));

hold on
plot(y,imsize-x+1,'rx')
depressions=[x y];
pause(1)
% given rainflow on the terrain, we wish to calculate:
%   1)how much water does a pixel receive
inflow=log2(inflow+1);
subplot(1,3,2)
imagesc(-1*inflow(end:-1:1,:))%invert to make water blue
title('inflow amounts')

%   2)where does this water go, using calculated neighbours
% connected=flow_graph(find(flow_graph(:,1)),:);
% [fromx, fromy]=ind2sub([imsize imsize],connected(:,1));
% [tox, toy]=ind2sub([imsize imsize],connected(:,2));
% directions=[fromx-tox fromy-toy];
[Fx, Fy]=gradient(Z);
directions=-1*[Fx(:) Fy(:)];
directions=normr(directions).*repmat(sum(abs(directions),2),1,2);%normalize outflow direction vectors, while keeping those that don't go anywhere 0
directions=directions.*repmat(inflow(:),1,2);
[X Y]=meshgrid(1:imsize,1:imsize);
subplot(1,3,3)
quiver(X,Y,reshape(directions(:,1),imsize,imsize),reshape(directions(:,2),imsize,imsize))
axis([1 imsize 1 imsize])
title('outflow')
%   2)where does this water go, using gradient on terrain
% assume uniformly distributed rainfall, and handle local depressions
% pixels which cannot be used, because they are local depressions
%depressions'
pause(0.2)%drawnow

%% save
Z = Z_orig;
clearvars -except depressions directions Z
save('mountains_with_precomputed_flow_vectors')