%correct the coordinates according to neighbours
% pca_level is [n_pixels x n_pca_channels] (16384 x 24)
% pca_image is [n_pixels x n_image_channels] (16384 x 1)
function S_new = correct_s_pca(S_old, pca_level, pca_image, m, channels)
    neighbourhood = 5;
    range=[(neighbourhood-1)/-2:(neighbourhood-1)/2];
    [X,Y]=meshgrid(range, range);
    %remove self
    X=X(setdiff(1:end,neighbourhood^2/2+0.5));
    Y=Y(setdiff(1:end,neighbourhood^2/2+0.5));
    
%     %[S_x, S_y]=ind2sub([m,m],S_old);
%     S_x = reshape(S_old(:,:,1),size(S_old,1)^2,1);
%     S_y = reshape(S_old(:,:,2),size(S_old,1)^2,1);
%     S_x=mod(repmat(S_x,1,24)+repmat(X,size(S_old,1)^2,1)-1,m)+1;
%     S_y=mod(repmat(S_y,1,24)+repmat(Y,size(S_old,1)^2,1)-1,m)+1;
%     neighbours=sub2ind([m,m],S_x,S_y);%neighbour indexes
%     neighbourhoods=pca_image(neighbours);%neighbour values
%     
%     %by referencing the pixels in the pca level
%     
%     % a*a'+b*b'-2*a*b' = sum((a-b).^2), so
%     %all_distances = sqrt( bsxfun(@plus,sum(pca_level.^2,2),sum(neighbourhoods.^2,2)') - 2*(pca_level*neighbourhoods') ); %requires too much memory for all_distances, in_size^2 by out_size^2, 16384 x 65536
%     %[~,nearest_neighbours]=min(all_distances);
%     [~,nearest_neighbours] = pdist2(pca_level,neighbourhoods,'euclidean','Smallest',1);
%     [subx,suby]=ind2sub([m,m],reshape(nearest_neighbours,[size(S_old,1),size(S_old,2)]));
%     S_new(:,:,1)=subx;
%     S_new(:,:,2)=suby;

    function [subx, suby] = correct_one_pixel(x, y)
        pixels_in_S_old = sub2ind([size(S_old,1) size(S_old,2)],mod(X+x-1,size(S_old,1))+1,mod(Y+y-1,size(S_old,2))+1);
        S_x=S_old(:,:,1);
        S_y=S_old(:,:,2);
        neighbours = pca_image(sub2ind([m,m],S_x(pixels_in_S_old),S_y(pixels_in_S_old)),:)';
        [~,i]=min(sum(sum(abs(pca_level-repmat(reshape(double(neighbours),[1 neighbourhood^2-1 channels]),[m^2,1,1])),3),2));
        [subx,suby]=ind2sub([m,m],i);
    end

    S_new = zeros(size(S_old));
    %for each pixel in the synthesized texture
    for x=1:size(S_old,1)
        for y=1:size(S_old,2)
            [subx, suby] = correct_one_pixel(x, y);
            S_new(x,y,1)=subx;
            S_new(x,y,2)=suby;
        end
    end
end