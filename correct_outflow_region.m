% correct a region which does not flow out of the synthesized map
function S_new = correct_outflow_region(S_old, Nexemplar, pixels_in, m, region_id, segmentation)
    % there is one of two possible problems with this region:
    %   1)the lowest point is not on the boundary
    %   2)the lowest point on the boundary is not adjacent to a lower point in another region
    %find points on boundary
    segment=segmentation==region_id;
    adjacent_segment_elements = conv2(double(segment), ones(3), 'same')-segment;
    diff_with_max_adjacencies=adjacent_segment_elements-conv2(ones(size(S_old,1)), ones(3), 'same')+ones(size(S_old,1));
    diff_with_max_adjacencies(~segment)=0;
    boundary_pixels=find(diff_with_max_adjacencies~=0);
    %find the lowest point on the boundary
    S_old_vector=reshape(S_old,size(S_old,1)^2,2);
    [lowest_boundary_height,lowest_pixel]=min(pixels_in(sub2ind([m,m],S_old_vector(boundary_pixels,1),S_old_vector(boundary_pixels,2))));
    %find lowest adjacent point outside of the segment
    [a, b]=ind2sub([size(S_old,1) size(S_old,2)],boundary_pixels(lowest_pixel))
    [X, Y]=meshgrid(max(1,a-1):min(size(S_old,1),a+1),max(1,b-1):min(size(S_old,1),b+1));
    outflow_candidates=sub2ind([size(S_old,1) size(S_old,2)],X,Y)';
    outflow_candidates=outflow_candidates(~segment(outflow_candidates));
    outflow_height=min(pixels_in(sub2ind([m,m],S_old_vector(outflow_candidates,1),S_old_vector(outflow_candidates,2))));
    if(lowest_boundary_height>outflow_height)
        disp('the problem is not on the boundary - fix inside')
        S_new = S_old;
    else
        %correct the lowest pixel on the boundary to be higher than its lowest
        %neighbour outside of the region
        neighbourhood = 5;
        range=[(neighbourhood-1)/-2:(neighbourhood-1)/2];
        [X,Y]=meshgrid(range, range);
        %remove self
        X=X(setdiff(1:end,neighbourhood^2/2+0.5));
        Y=Y(setdiff(1:end,neighbourhood^2/2+0.5));

        pixels_in_S_old = sub2ind([size(S_old,1) size(S_old,2)],abs(size(S_old,1)-1-abs(size(S_old,1)-(X+a)))+1,abs(size(S_old,2)-1-abs(size(S_old,2)-(Y+b)))+1);
        S_x=S_old(:,:,1);
        S_y=S_old(:,:,2);
        neighbours = pixels_in(sub2ind([m,m],S_x(pixels_in_S_old),S_y(pixels_in_S_old)));
        %only select out of exemplar pixels higher than outflow_height
        allowable_pixels=find(pixels_in>outflow_height);
        [~,i]=min(sum(sum(abs(Nexemplar(allowable_pixels,:)-repmat(neighbours,[length(allowable_pixels),1])),3),2));
        [subx,suby]=ind2sub([m,m],allowable_pixels(i));

        S_new = S_old;
        S_new(a,b,1)=subx;
        S_new(a,b,2)=suby;
    end
end