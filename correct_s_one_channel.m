%correct the coordinates according to neighbours
function S_new = correct_s_one_channel(S_old, Nexemplar, pixels_in, m)
    [sx, sy, ~]=size(S_old);
    neighbourhood = 5;
    range=[(neighbourhood-1)/-2:(neighbourhood-1)/2];
    [X,Y]=meshgrid(range, range);
    %remove self
    X=X(setdiff(1:end,neighbourhood^2/2+0.5));
    Y=Y(setdiff(1:end,neighbourhood^2/2+0.5));

    S_new = zeros(size(S_old));
    S_x=S_old(:,:,1);
    S_y=S_old(:,:,2);
    
    function i_new = correct_one_pixel(i_each)
        % abs(63-abs(64-r))+1 converts numbers as expected when
        % flipping the image. So [-5 0 1 35 64 65 70] becomes [7 2 1 35 64 63 58]
        [x,y]=ind2sub([sx sy],i_each);%y=ceil(i_each/sx);x=mod(i_each-1,sx)+1;
        pixels_in_S_old = sub2ind([sx sy],abs(sx-1-abs(sx-(X+x)))+1,abs(sy-1-abs(sy-(Y+y)))+1);
        neighbours = pixels_in(sub2ind([m,m],S_x(pixels_in_S_old),S_y(pixels_in_S_old)));

        [~,i_new]=min(sum(sum(abs(Nexemplar-repmat(neighbours,[m^2,1])),3),2));
    end

    %for each pixel in the synthesized texture, run on all CPUs
    i=1:(sx*sy);
    indexes = arrayfun(@correct_one_pixel, i);

    [subx,suby]=ind2sub([m,m],reshape(indexes,[sx,sy]));
    S_new(:,:,1)=subx;
    S_new(:,:,2)=suby;
end