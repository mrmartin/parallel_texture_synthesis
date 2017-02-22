%correct the coordinates according to neighbours
function S_new = correct_s(S_old, Nexemplar, pixels_in, m, l)
    neighbourhood = 5;
    range=m/2^l*[(neighbourhood-1)/-2:(neighbourhood-1)/2];
    [X,Y]=meshgrid(range, range);
    %remove self
    X=X(setdiff(1:end,neighbourhood^2/2-0.5));
    Y=Y(setdiff(1:end,neighbourhood^2/2-0.5));

    S_new = zeros(size(S_old));
    %for each pixel in the synthesized texture
    for x=1:size(S_old,1)
        for y=1:size(S_old,2)
            neighbours = pixels_in(sub2ind([m,m],mod(X+x-1,m)+1,mod(Y+y-1,m)+1),:);
            [~,i]=min(sum(sum(abs(Nexemplar-repmat(reshape(double(neighbours),[1 24 3]),[m^2,1,1])),3),2));
            [subx,suby]=ind2sub([m,m],i);
            S_new(x,y,1)=subx;
            S_new(x,y,2)=suby;
        end
    end
end