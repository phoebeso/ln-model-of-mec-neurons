function [hd_grid,dirVec,direction] = hd_map(direction,nbins)

%compute head direction
direction(direction < 0) = direction(direction<0)+2*pi; % go from 0 to 2*pi, without any negative numbers

hd_grid = zeros(length(direction),nbins);
dirVec = 2*pi/nbins/2:2*pi/nbins:2*pi-2*pi/nbins/2;

for i = 1:numel(direction)
    
    % figure out the hd index
    [~, idx] = min(abs(direction(i)-dirVec));
    hd_grid(i,idx) = 1;
  
end

return