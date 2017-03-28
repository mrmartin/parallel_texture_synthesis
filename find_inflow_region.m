% using watershed basin segmentations, identify regions which need correction
function region_ids = find_inflow_region(flow_targets, map_size)
    [a, b]=ind2sub([map_size map_size],flow_targets);
    bad = intersect(find(b>1 & b<map_size),find(a>1 & a<map_size));
    region_ids = flow_targets(bad);
end