% perform local elementwise grid operations on gpu
function grid_gpu_test
    gridSize = 8;
    grid_in = rand(gridSize);
    
    %pad the grid to handle boundary effects
    grid = zeros(gridSize+2);
    grid(2:end-1,2:end-1) = grid_in
    grid = gpuArray(grid);
    
    %find the minimum of a vector of size s
    function index = formin(v,s)
        minv=-Inf;
        index=0;
        for i=1:s
            if v(i)<minv
                minv = v(i);
                index = i;
            end
        end
    end
    
    function X = updateParentGrid(row, col)
        rowU = row-1;  rowD = row+1;
        colL = col-1;  colR = col+1;
        % sum neighbors
        neighbors ...
            = grid(rowU,colL) + grid(row,colL) + grid(rowD,colL) ...
            + grid(rowU,col)                   + grid(rowD,col) ...
            + grid(rowU,colR) + grid(row,colR) + grid(rowD,colR);
        % find lowest neighbor
        lowest = min(min(min(min(min(min(min(grid(rowU,colL),grid(row,colL)),grid(rowD,colL)),grid(rowU,col)),grid(rowD,col)),grid(rowU,colR)),grid(row,colR)),grid(rowD,colR));
        % difference between mean and min
        X = neighbors./8-lowest;
    end

    %exclude boundaries
    rows = gpuArray.colon(2, gridSize+1)';
    cols = gpuArray.colon(2, gridSize+1);
    grid = arrayfun(@updateParentGrid, rows, cols)
end