%% Terrain to watershed basins
function [colors, outflow] = watershed_basins(Z)
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
                    flow_graph(sub2ind(size(Z),x,y),:)=[sub2ind(size(Z),x,y) sub2ind(size(Z),x,y)];
                    outflow = [outflow;sub2ind(size(Z),x,y)];
                else
                    distribution=max(0,Z(x,y)-Z(neighbors))/sum(max(0,Z(x,y)-Z(neighbors)));
                    [~,ind]=max(distribution);
                    flow_graph(sub2ind(size(Z),x,y),1)=sub2ind(size(Z),x,y);
                    flow_graph(sub2ind(size(Z),x,y),2)=neighbors(ind);                
                end
            end
        end
    end

    %go from the root to the leaves, simplifying connections
    connected=flow_graph(find(flow_graph(:,1)),:);
    disconnected=find(isnan(Z));

    %connected(~,1) flows into connected(~,2)
    itteration=0;
    while ~isempty(setdiff(connected(:,2),outflow))
        itteration=itteration+1;
        %disp(['itteration ' num2str(itteration)]);
        for i=1:length(connected)
        %    i
            connected(i,2)=connected(connected(:,1)==connected(i,2),2);
        end
    end
    
    colors=Z;
    colors(disconnected)=0;
    colors(connected(:,1))=connected(:,2);
end