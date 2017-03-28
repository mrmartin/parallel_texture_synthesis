%correct the coordinates according to neighbours
function S_new = correct_s_one_channel(S_old, Nexemplar, pixels_in, m)
    neighbourhood = 5;
    range=[(neighbourhood-1)/-2:(neighbourhood-1)/2];
    [X,Y]=meshgrid(range, range);
    %remove self
    X=X(setdiff(1:end,neighbourhood^2/2+0.5));
    Y=Y(setdiff(1:end,neighbourhood^2/2+0.5));

    S_new = zeros(size(S_old));
    %for each pixel in the synthesized texture
    for x=1:size(S_old,1)
        for y=1:size(S_old,2)
            % abs(63-abs(64-r))+1 converts numbers as expected when
            % flipping the image. So [-5 0 1 35 64 65 70] becomes [7 2 1 35 64 63 58]
            pixels_in_S_old = sub2ind([size(S_old,1) size(S_old,2)],abs(size(S_old,1)-1-abs(size(S_old,1)-(X+x)))+1,abs(size(S_old,2)-1-abs(size(S_old,2)-(Y+y)))+1);
            S_x=S_old(:,:,1);
            S_y=S_old(:,:,2);
            neighbours = pixels_in(sub2ind([m,m],S_x(pixels_in_S_old),S_y(pixels_in_S_old)));

            [~,i]=min(sum(sum(abs(Nexemplar-repmat(neighbours,[m^2,1])),3),2));
            [subx,suby]=ind2sub([m,m],i);
            S_new(x,y,1)=subx;
            S_new(x,y,2)=suby;
        end
    end
end