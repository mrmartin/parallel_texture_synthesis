function grid_gpu_test_2
    gridSize = 250;
    grid = rand(gridSize);%gpuArray(rand(gridSize));

    neighbourhoods = zeros((gridSize-2)^2,9);
    [nX, nY]=meshgrid(-1:1,-1:1);
    for x=2:gridSize-1
        for y=2:gridSize-1
            neighbourhoods(sub2ind([gridSize-2 gridSize-2],x-1,y-1),:)=reshape(grid(sub2ind([gridSize gridSize],y+nY,x+nX)),1,9);
        end
    end

    function lowest = min_diff(elem)
        [x, y]=ind2sub([gridSize gridSize],elem);
        tmp=repmat(reshape(grid(sub2ind([gridSize gridSize],y+nY,x+nX)),1,9),(gridSize-2)^2,1);
        [~,lowest] = min(sum(abs(neighbourhoods-tmp),2));
        %account for offset caused by ignoring first rown and column
        [a b]=ind2sub([gridSize-2 gridSize-2],lowest);
        lowest=sub2ind([gridSize gridSize],a+1,b+1);
    end

    [X, Y]=meshgrid(2:gridSize-1,2:gridSize-1);
    elems = sub2ind([gridSize gridSize],X(:),Y(:))'%gpuArray.colon(2, gridSize)';
    arrayfun(@min_diff, elems)
end