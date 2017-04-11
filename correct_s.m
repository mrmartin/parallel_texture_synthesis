%correct the coordinates according to neighbours
function S_new = correct_s(S_old, Nexemplar, pixels_in, m, candidates, neighbourhood, synthesis_method)
    range=[(neighbourhood-1)/-2:(neighbourhood-1)/2];
    [X,Y]=meshgrid(range, range);
    %remove self
    X=X(setdiff(1:end,neighbourhood^2/2+0.5));
    Y=Y(setdiff(1:end,neighbourhood^2/2+0.5));
    
    candidate_neighbourhood = 3;
    range=[(candidate_neighbourhood-1)/-2:(candidate_neighbourhood-1)/2];
    [cX,cY]=meshgrid(range, range);%do not remove self

    S_x=S_old(:,:,1);
    S_y=S_old(:,:,2);
    
    flow_magnitudes=sqrt(sum(pixels_in.^2,2));

    S_new_i = zeros(size(S_old,1)*size(S_old,2),1);
    %for each pixel in the synthesized texture
    parfor in_i=1:(size(S_old,1)*size(S_old,2))
        [x,y]=ind2sub([size(S_old,1),size(S_old,2)],in_i);
        pixels_in_S_old = sub2ind([size(S_old,1) size(S_old,2)],mod(X+x-1,size(S_old,1))+1,mod(Y+y-1,size(S_old,2))+1);
        neighbours = pixels_in(sub2ind([m,m],S_x(pixels_in_S_old),S_y(pixels_in_S_old)),:);
        if(synthesis_method==1)%match against all pixels
            [~,i]=min(sum(sum(abs(Nexemplar-repmat(reshape(double(neighbours),[1 neighbourhood^2-1 3]),[m^2,1,1])),3),2));
        elseif(synthesis_method==2)%constrain each pixel to be within a range of flow magnitude: at least what's going in and at most what's going out
            if(x>floor(neighbourhood/2) & x<size(S_old,1)-floor(neighbourhood/2)+1 & y>floor(neighbourhood/2) & y<size(S_old,2)-floor(neighbourhood/2)+1)%for non-near-edge pixels
                four_neighbours=sub2ind([size(S_old,1) size(S_old,2)],mod([x-1; x; x+1; x]-1,size(S_old,1))+1,mod([y; y-1; y; y+1]-1,size(S_old,1))+1);
                four_flows=pixels_in(sub2ind([m,m],S_x(four_neighbours),S_y(four_neighbours)),2:3);
                neighbour_flow=four_flows([5 2 7 4]).*[1 1 -1 -1];%get horizontal element of horizontal neighbours and vertical element of vertical neighbours
                out=abs(min([0 neighbour_flow]));%what's going out
                in=max([0 neighbour_flow]);%what's going in
                if(in<=out)
                    flow_candidates=find(flow_magnitudes<out+0.01 & flow_magnitudes>in-0.01);%add threshold
                else
                    flow_candidates=find(flow_magnitudes>in-0.01);%add threshold
                end
                [~,i]=min(sum(sum(abs(Nexemplar(flow_candidates,:,:)-repmat(reshape(double(neighbours),[1 neighbourhood^2-1 3]),[length(flow_candidates),1,1])),3),2));
                i=flow_candidates(i);
            else%pixels on the edge
                relevant=find(x+X>0 & x+X<=size(S_old,1) & y+Y>0 & y+Y<=size(S_old,2));%which neighbours are relevant here? - ignore toroidal element pixels
                [~,i]=min(sum(sum(abs(Nexemplar(:,relevant,:)-repmat(reshape(double(neighbours(relevant,:)),[1 length(relevant) 3]),[m^2,1,1])),3),2));
            end
        else%match against pixels which are neighbours of the candidates
            if(x>floor(neighbourhood/2) & x<size(S_old,1)-floor(neighbourhood/2)+1 & y>floor(neighbourhood/2) & y<size(S_old,2)-floor(neighbourhood/2)+1)%for non-near-edge pixels
                candidate_pixels_in_S_old = sub2ind([size(S_old,1) size(S_old,2)],mod(cX+x-1,size(S_old,1))+1,mod(cY+y-1,size(S_old,2))+1);
                neighbours_of_candidates = candidates(sub2ind([m,m],S_x(candidate_pixels_in_S_old(:)),S_y(candidate_pixels_in_S_old(:))),:);
                neighbours_of_candidates = neighbours_of_candidates(neighbours_of_candidates~=Inf);%remove unusable candidates
                if(isempty(neighbours_of_candidates))%may happen at lowest level, where speed is not an issue
                    %disp('length(neighbours_of_candidates)==0 in correct_s');
                    %match against all pixels
                    [~,i]=min(sum(sum(abs(Nexemplar-repmat(reshape(double(neighbours),[1 neighbourhood^2-1 3]),[m^2,1,1])),3),2));
                else
                    [~,i]=min(sum(sum(abs(Nexemplar(neighbours_of_candidates,:,:)-repmat(reshape(double(neighbours),[1 neighbourhood^2-1 3]),[length(neighbours_of_candidates),1,1])),3),2));
                    i=neighbours_of_candidates(i);
                end
            else
                relevant=find(x+X>0 & x+X<=size(S_old,1) & y+Y>0 & y+Y<=size(S_old,2));%which neighbours are relevant here? - ignore toroidal element pixels
                [~,i]=min(sum(sum(abs(Nexemplar(:,relevant,:)-repmat(reshape(double(neighbours(relevant,:)),[1 length(relevant) 3]),[m^2,1,1])),3),2));
            end
        end
        S_new_i(in_i) = i;
    end
    %convert back to subx suby coords
    [X, Y]=ind2sub([m,m],reshape(S_new_i,size(S_old,1),size(S_old,2)));
    S_new(:,:,1)=X;
    S_new(:,:,2)=Y;
end