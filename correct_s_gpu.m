%correct the coordinates according to neighbours
function S_new = correct_s_gpu(S_old, Nexemplar, pixels_in, m)
    neighbourhood = 5;
    range=[(neighbourhood-1)/-2:(neighbourhood-1)/2];
    [X,Y]=meshgrid(range, range);
    %remove self
    X=X(setdiff(1:end,neighbourhood^2/2-0.5));
    Y=Y(setdiff(1:end,neighbourhood^2/2-0.5));

    size_S_old_1 = size(S_old,1);
    size_S_old_2 = size(S_old,2);

    S_x=S_old(:,:,1);
    S_y=S_old(:,:,2);
    all3=[1 2 3];
    all24=1:24;
    distances_init=zeros(1,m^2);
    function index = correct_one_pixel(x, y)
        pixels_in_S_old = sub2ind_2d(size_S_old_1,mod(X+x-1,size_S_old_1)+1,mod(Y+y-1,size_S_old_2)+1);
        neighbours = pixels_in(sub2ind_2d(m,S_x(pixels_in_S_old),S_y(pixels_in_S_old)),all3);
        % size(abs(Nexemplar-repmat(reshape(double(neighbours),[1 neighbourhood^2-1 3]),[m^2,1,1]))) = 4096 x 24 x 3
        % cur_distances = abs(Nexemplar-repmat(reshape(double(neighbours),[1 neighbourhood^2-1 3]),[m^2,1,1]));
        distances=distances_init;
        for i=1:4096
            cur_distance = abs(Nexemplar(i,all24,all3)-neighbours)
            for j=1:24
                distances(i) = distances(i) + cur_distance(j,1) + cur_distance(j,2) + cur_distance(j,3);
            end
        end
        %distances = sum(sum(abs(Nexemplar-repmat(reshape(double(neighbours),[1 neighbourhood^2-1 3]),[m^2,1,1])),3),2);
        
        minv=-Inf;
        index=0;
        for i=1:m^2
            if distances(i)<minv
                minv = v(i);
                index = i;
            end
        end
    end

    S_new = zeros(size(S_old));
    %for each pixel in the synthesized texture
    x=gpuArray.colon(1, size(S_old,1))'
    y=gpuArray.colon(1, size(S_old,2))
    index = arrayfun(@correct_one_pixel,x, y)
    %[subx,suby]=ind2sub([m,m],index);
    %S_new(x,y,1)=subx;
    %S_new(x,y,2)=suby;
    pause
end