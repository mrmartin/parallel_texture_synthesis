%correct the coordinates according to neighbours
function S_new = correct_s_raise_single_pixel(S_old, Nexemplar, pixels_in, m, depression_pixel)
    S_new = S_old;
    %find pixel to raise
    [x, y]=ind2sub([size(S_old,1) size(S_old,2)],depression_pixel);
    
    %find lowest neighbour
    S_old_vector=reshape(S_old,size(S_old,1)^2,2);
    [X, Y]=meshgrid(max(1,x-1):min(size(S_old,1),x+1),max(1,y-1):min(size(S_old,1),y+1));
    %remove self
    tmp=setdiff([X(:) Y(:)],[x y],'rows');
    X=tmp(:,1);
    Y=tmp(:,2);
    outflow_candidates=sub2ind([size(S_old,1) size(S_old,2)],X,Y)';
    min_height=min(pixels_in(sub2ind([m,m],S_old_vector(outflow_candidates,1),S_old_vector(outflow_candidates,2))));
    
    neighbourhood = 5;
    range=[(neighbourhood-1)/-2:(neighbourhood-1)/2];
    [X,Y]=meshgrid(range, range);
    %remove self
    X=X(setdiff(1:end,neighbourhood^2/2+0.5));
    Y=Y(setdiff(1:end,neighbourhood^2/2+0.5));
    
    pixels_in_S_old = sub2ind([size(S_old,1) size(S_old,2)],abs(size(S_old,1)-1-abs(size(S_old,1)-(X+x)))+1,abs(size(S_old,2)-1-abs(size(S_old,2)-(Y+y)))+1);
    S_x=S_old(:,:,1);
    S_y=S_old(:,:,2);
    neighbours = pixels_in(sub2ind([m,m],S_x(pixels_in_S_old),S_y(pixels_in_S_old)));
    allowable_pixels=find(pixels_in>min_height);
    if(isempty(allowable_pixels))
        disp('there is a pixel which needs to be raised, but it can''t');
    end
    [~,i]=min(sum(sum(abs(Nexemplar(allowable_pixels,:)-repmat(neighbours,[length(allowable_pixels),1])),3),2));
    [subx,suby]=ind2sub([m,m],allowable_pixels(i));
    S_new(x,y,1)=subx;
    S_new(x,y,2)=suby;
end