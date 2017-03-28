function paralleldemo_gpu_stencil()
    gridSize = 100;
    numGenerations = 100;
    initialGrid = (rand(gridSize,gridSize) > .75);
    gpu = gpuDevice();

       function X = updateGrid(X, N)
            p = [1 1:N-1];
            q = [2:N N];
            % Count how many of the eight neighbors are alive.
            neighbors = X(:,p) + X(:,q) + X(p,:) + X(q,:) + ...
                X(p,p) + X(q,q) + X(p,q) + X(q,p);
            % A live cell with two live neighbors, or any cell with
            % three live neighbors, is alive at the next step.
            X = (X & (neighbors == 2)) | (neighbors == 3);
       end

    % CPU
    grid = initialGrid;
    for generation = 1:numGenerations
        grid = updateGrid(grid, gridSize);
    end

    % GPU
    grid = gpuArray(initialGrid);
    for generation = 1:numGenerations
        grid = updateGrid(grid, gridSize);
    end

    % GPU arrayfun
    grid = gpuArray(initialGrid);
        function X = updateParentGrid(row, col, N)
            % Take account of boundary effects
            rowU = max(1,row-1);  rowD = min(N,row+1);
            colL = max(1,col-1);  colR = min(N,col+1);
            % Count neighbors
            neighbors ...
                = grid(rowU,colL) + grid(row,colL) + grid(rowD,colL) ...
                + grid(rowU,col)                   + grid(rowD,col) ...
                + grid(rowU,colR) + grid(row,colR) + grid(rowD,colR);
            % A live cell with two live neighbors, or any cell with
            % three live neighbors, is alive at the next step.
            X = (grid(row,col) & (neighbors == 2)) | (neighbors == 3);
        end

    rows = gpuArray.colon(1, gridSize)';
    cols = gpuArray.colon(1, gridSize);
    for generation = 1:numGenerations
        grid = arrayfun(@updateParentGrid, rows, cols, gridSize);
    end
end