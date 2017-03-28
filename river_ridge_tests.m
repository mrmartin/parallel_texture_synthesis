%cd ~/Downloads/DEM_1/
%[Z, R]=sdtsdemread('1106CEL0.DDF');

f=fopen('~/Desktop/exported terrains/canyons.raw');
raw=fread(f,1025*1025,'uint16');
fclose(f);
extremes=csvread('~/Desktop/exported terrains/canyons.csv');
low=extremes(1);
high=extremes(2);
raw=reshape(raw,1025,1025);
Z=low+(raw*(high-low)/max(raw));
min(Z)
max(Z)
Z=low+(raw*(high-low)/max(raw(:)));

%Z=Z(1:5:300,1:5:300);

%% Numerical terrain analysis
% subplot(1,3,1)
% %mapshow(Z,R,'DisplayType','surface')
% view(3)
% axis normal
% 
% subplot(1,3,2)
% [Fx, Fy]=gradient(Z);
% quiver(X,Y,Fx,Fy)
% 
% subplot(1,3,3)
% [X,Y]=meshgrid(1:size(Z,1),1:size(Z,2));
% [K,H,P1,P2] = surfature(X,Y,Z);
% surf(X,Y,Z,H,'facecolor','interp');
% set(gca,'clim',[-1,1])

%% Brute force rain flow simulation
% close
% h=figure
% set(h, 'Position', [0 0 1024 800])
% %for each point, flow from neighbouring points
% rain=ones(size(Z));
% flow_tmp=zeros(size(Z));
% flow=rain;
% for i=1:20
%     for x=1:size(Z,1)
%         for y=1:size(Z,2)
%             if(~isnan(Z(x,y)))
%                 %find lower neighbours
%                 neighbors=[];
%                 if(x>1)
%                     neighbors=[neighbors;x-1 y];
%                     if(y>1)
%                         neighbors=[neighbors;x-1 y-1];
%                     end
%                 end
%                 if(x<size(Z,1))
%                     neighbors=[neighbors;x+1 y];
%                     if(y<size(Z,2))
%                         neighbors=[neighbors;x+1 y+1];
%                     end
%                 end
%                 if(y>1)
%                     neighbors=[neighbors;x y-1];
%                     if(x<size(Z,1))
%                         neighbors=[neighbors;x+1 y-1];
%                     end
%                 end
%                 if(y<size(Z,2))
%                     neighbors=[neighbors;x y+1];
%                     if(x>1)
%                         neighbors=[neighbors;x-1 y+1];
%                     end
%                 end
%                 neighbors=sub2ind(size(Z),neighbors(:,1),neighbors(:,2));
%                 if(sum(max(0,Z(x,y)-Z(neighbors)))==0)
%                     %disp(['there are no lower neighbours to [' num2str(x) ', ' num2str(y) ']']);
%                     %there are three options:
%                         %1) map imperfection, a lower neighbour is within a few surrounding pixels
%                         %2) outflow from the map -> drainage basin sink
%                         %3) local minimum, undefined behavior
%                     %try to find lower pixel/map boundary in the neighborhood
%                     for k=2:4
%                         %get addresses of all surrounding squares
%                         [X,Y]=meshgrid(x-k:x+k,y-k:y+k);
%                         %if any lie outside the map, they are outflow candidates
%                         if(max(X(:)<1 | X(:)>size(Z,1) | Y(:)<1 | Y(:)>size(Z,2)))
%                             disp(['[' num2str(x) ', ' num2str(y) '] is an outflow candidate']);
%                             break
%                         else
%                             neighbors=sub2ind(size(Z),X(:),Y(:));
%                             %remove self from list
%                             neighbors=neighbors(setdiff(1:length(neighbors),(length(neighbors)+1)/2));
%                             if(sum(max(0,Z(x,y)-Z(neighbors)))==0)
%                                 %disp(['there is no lower pixel within neighbourhood k=' num2str(k) ' of [' num2str(x) ', ' num2str(y) ']']);
%                             else
%                                 break
%                             end
%                         end
%                     end
%                 end
%                 distribution=max(0,Z(x,y)-Z(neighbors))/sum(max(0,Z(x,y)-Z(neighbors)));
%                 distribution(isnan(distribution))=0;
%                 if(max(isnan(distribution)))
%                     x
%                     y
%                 end
%                 flow_tmp(neighbors)=flow_tmp(neighbors)+flow(x,y).*distribution;
%             end
%         end
%     end
%     flow=flow_tmp;
%     imagesc(-1*flow)
% %    saveas(h,['high_res_flow_' num2str(i) '.png']) 
%     pause
% end

%% Terrain to flow graph
flow_graph = zeros(numel(Z),2);
outflow=[];
for x=1:size(Z,1)
    for y=1:size(Z,2)
        if(~isnan(Z(x,y)))
            %find lowest neighbour
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
                %disp(['there are no lower neighbours to [' num2str(x) ', ' num2str(y) ']']);
                %there are three options:
                    %1) map imperfection, a lower neighbour is within a few surrounding pixels
                    %2) outflow from the map -> drainage basin sink
                    %3) local minimum, undefined behavior
                %try to find lower pixel/map boundary in the neighborhood
                for k=2:8
                    %get addresses of all surrounding squares
                    [X,Y]=meshgrid(x-k:x+k,y-k:y+k);
                    %if any lie outside the map, they are outflow candidates
                    if(max(X(:)<1 | X(:)>size(Z,1) | Y(:)<1 | Y(:)>size(Z,2)))
                        break
                    %if any of the neighbors are NaN, they are considered
                    %outside the map, and the are also outflow
                    elseif(max(isnan(Z(sub2ind(size(Z),X(:),Y(:))))))
                        break
                    else
                        neighbors=sub2ind(size(Z),X(:),Y(:));
                        %remove self from list
                        neighbors=neighbors(setdiff(1:length(neighbors),(length(neighbors)+1)/2));
                        if(sum(max(0,Z(x,y)-Z(neighbors)))==0)
                            %disp(['there is no lower pixel within neighbourhood k=' num2str(k) ' of [' num2str(x) ', ' num2str(y) ']']);
                        else
                            break
                        end
                    end
                end
                if(max(X(:)<1 | X(:)>size(Z,1) | Y(:)<1 | Y(:)>size(Z,2)) || max(isnan(Z(sub2ind(size(Z),X(:),Y(:))))))
                    disp(['[' num2str(x) ', ' num2str(y) '] is an outflow candidate']);
                    outflow = [outflow;sub2ind(size(Z),x,y)];
                    flow_graph(sub2ind(size(Z),x,y),1)=sub2ind(size(Z),x,y);
                    flow_graph(sub2ind(size(Z),x,y),2)=sub2ind(size(Z),x,y);
                elseif(sum(max(0,Z(x,y)-Z(neighbors)))==0)
                    disp(['k=' num2str(k) ' was not enough at [' num2str(x) ', ' num2str(y) '], so it''s probably is a local depression.'])
                    outflow = [outflow;sub2ind(size(Z),x,y)];
                    flow_graph(sub2ind(size(Z),x,y),1)=sub2ind(size(Z),x,y);
                    flow_graph(sub2ind(size(Z),x,y),2)=sub2ind(size(Z),x,y);
                else
                    distribution=max(0,Z(x,y)-Z(neighbors))/sum(max(0,Z(x,y)-Z(neighbors)));
                    [~,ind]=max(distribution);
                    flow_graph(sub2ind(size(Z),x,y),1)=sub2ind(size(Z),x,y);
                    flow_graph(sub2ind(size(Z),x,y),2)=neighbors(ind);
                end
            end
            if(flow_graph(sub2ind(size(Z),x,y),1)==0)%nothing special has been done
                    distribution=max(0,Z(x,y)-Z(neighbors))/sum(max(0,Z(x,y)-Z(neighbors)));
                    [~,ind]=max(distribution);
                    flow_graph(sub2ind(size(Z),x,y),1)=sub2ind(size(Z),x,y);
                    flow_graph(sub2ind(size(Z),x,y),2)=neighbors(ind);                
            end
        end
    end
end
numel(Z)
length(find(flow_graph(:,1)))+length(find(isnan(Z)))

%go from the root to the leaves, simplifying connections
connected=flow_graph(find(flow_graph(:,1)),:);
disconnected=find(isnan(Z));

%connected(~,1) flows into connected(~,2)
itteration=0;
while ~isempty(setdiff(connected(:,2),outflow))
    itteration=itteration+1;
    disp(['itteration ' num2str(itteration)]);
    for i=1:length(connected)
    %    i
        connected(i,2)=connected(connected(:,1)==connected(i,2),2);
    end
end

colors=Z;
colors(disconnected)=0;
colors(connected(:,1))=connected(:,2);

[X,Y]=meshgrid(1:size(Z,2),1:size(Z,1));
shading interp;
surface(X,Y,Z,colors,'EdgeColor','none','FaceColor','texturemap')
view(3)
%segment colors
polcmap