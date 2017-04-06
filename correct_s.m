%correct the coordinates according to neighbours
function S_new = correct_s(S_old, Nexemplar, pixels_in, m, candidates, neighbourhood)
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
    function i = correct_one_pixel(x, y)
        pixels_in_S_old = sub2ind([size(S_old,1) size(S_old,2)],mod(X+x-1,size(S_old,1))+1,mod(Y+y-1,size(S_old,2))+1);
        neighbours = pixels_in(sub2ind([m,m],S_x(pixels_in_S_old),S_y(pixels_in_S_old)),:);
        %match against all pixels
        %[~,i]=min(sum(sum(abs(Nexemplar-repmat(reshape(double(neighbours),[1 neighbourhood^2-1 3]),[m^2,1,1])),3),2));
        %match against pixels which are neighbours of the candidates
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
    end

    S_new_i = zeros(size(S_old,1),size(S_old,2));
    %for each pixel in the synthesized texture
    for x=1:size(S_old,1)
        for y=1:size(S_old,2)
            S_new_i(x,y) = correct_one_pixel(x, y);
        end
    end
    %convert back to subx suby coords
    [X, Y]=ind2sub([m,m],S_new_i);
    S_new(:,:,1)=X;
    S_new(:,:,2)=Y;
end