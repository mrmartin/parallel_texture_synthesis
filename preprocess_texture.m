function [Nexemplar, candidates, pixels_in, m, input_levels] = preprocess_texture(im_in, neighbourhood, k_candidates);
    %imagesc(repmat(im_in,3,3,1))
    m=size(im_in,1);
    pixels_in=reshape(im_in,m^2,3);

    input_levels = log2(size(im_in,1));

    Nexemplar = cell(input_levels,1);
    candidates = cell(input_levels,1);
    tic
    for l=3:input_levels
        toc
        disp(['starting to perform preprocessing at level ' num2str(l)])
        tic
        range=m/2^l*[(neighbourhood-1)/-2:(neighbourhood-1)/2];
        [X,Y]=meshgrid(range, range);
        %remove self
        X=X(setdiff(1:end,neighbourhood^2/2+0.5));
        Y=Y(setdiff(1:end,neighbourhood^2/2+0.5));

        %only perform blur on lower levels
        if(l<input_levels)
            im_gaussian_blur = imgaussfilt(repmat(im_in,3,3),2^(input_levels-l-1));
            im_gaussian_blur = im_gaussian_blur(m+[1:m],m+[1:m],:);
        else
            im_gaussian_blur = im_in;
        end

        pixels_gaussian_blur=reshape(im_gaussian_blur,m^2,3);

        %prepare the neighbourhood of each pixel on a line
        Nexemplar{l}=zeros(m^2,neighbourhood^2-1,3);
        for x=1:m
            for y=1:m
                Nexemplar{l}(sub2ind([m,m],x,y),:,:)=reshape(pixels_gaussian_blur(sub2ind([m,m],mod(X+x-1,m)+1,mod(Y+y-1,m)+1),:),[1 neighbourhood^2-1 3]);
                if(max((mod(X+x-1,m)+1)~=(X+x)) || max((mod(Y+y-1,m)+1)~=(Y+y)))
                    Nexemplar{l}(sub2ind([m,m],x,y),1)=Inf;%discard those which reach outside the boundary
                end
            end
        end

        candidates{l} = zeros(m^2,k_candidates);
        %for each pixel of the input
        for x=1:m
            for y=1:m
                if(Nexemplar{l}(sub2ind([m,m],x,y),1)==Inf)%not an allowed pixel
                    candidates{l}(sub2ind([m,m],x,y),:)=Inf;
                else%prepare candidates
                    neighbours = reshape(pixels_gaussian_blur(sub2ind([m,m],mod(X+x-1,m)+1,mod(Y+y-1,m)+1),:),[1 neighbourhood^2-1 3]);
                    [~,i]=sort(sum(sum(abs(Nexemplar{l}-repmat(double(neighbours),[m^2,1,1])),3),2));
                    %candidates must be at least 5% of the image apart from each other
                    [dx, dy]=ind2sub([m,m],i);
                    far_enough=i(sqrt(sum([dx-dx(1) dy-dy(1)].^2,2))>m/20);
                    %also allow self
                    candidates{l}(sub2ind([m,m],x,y),:)=[i(1) far_enough(1:k_candidates-1)'];
                end
            end
        end
    end
end