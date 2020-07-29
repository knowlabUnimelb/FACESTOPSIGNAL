function olap = overlap(p, q, xi)

%% Compute the amount of olap between 95% HDI's
density = [p(:), q(:)];

cdfs = cumsum(density)./(ones(size(density,1), 1) * sum(density));
hdi = nan(2,size(density,2));
for i = 1:size(density,2)
    hdi(1,i) = xi(find(cdfs(:,i)<=.025, 1, 'last'));
    hdi(2,i) = xi(find(cdfs(:,i)>=.985, 1, 'first'));
end

% XY coordinates of line 1: (hdi(1,1),0)-(hdi(2,1),0)
% XY coordinates of line 2: (hdi(1,2),0)-(hdi(2,2),0)
pairs = allcomb(1:size(density,2), 1:size(density,2)); pairs(pairs(:,1) == pairs(:,2), :) = [];
pairs = unique(sort(pairs, 2), 'rows');
olap = nan(1,size(pairs,1));
for i = 1:size(pairs, 1)
    % Get pairs
    x = hdi(:,pairs(i,:));
    
    % Sort so that the lines are ordered appropriately
    if x(1,1) ~= x(1,2)
        [a, b] = sort(x(1,:), 2);
        x = [a; x(2,b)];
    else
        [a, b] = sort(x(2,:), 2);
        x = [x(1,b); a];
    end
    
    % Check olap
    if x(1,2) < x(2,1) && x(2,1) <= x(2,2) % Then the lines olap
        olap(i) = (x(2,1) - x(1,2))./(x(2,2) - x(1,1)); % Proportion olap
    elseif x(1,2) < x(2,1) && x(2,1) > x(2,2) % Then they completely olap
        olap(i) = 1;
    else
        olap(i) = 0; % If there's no olap, then the proportion olap is 0
    end
end