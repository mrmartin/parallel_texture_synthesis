%sub2ind for 2d matrices, accepts x and y scalars or vectors of the same size
function ndx = sub2ind_2d(size_1,x,y)
    ndx = x + (y - 1).*size_1(1);
end
