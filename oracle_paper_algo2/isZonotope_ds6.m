function [isZono, generators] = isZonotope_ds6(vertices)
    % Initialize
    [n, d] = size(vertices);
    isZono = false;
    generators = [];
    
    % Step 1: Get edges via convex hull
    try
        K = convhulln(vertices);
        edges = getUniqueEdges(K, vertices);
    catch
        return
    end
    
    % Step 2: Process edges
    [G_cell, mu] = processEdges(edges, d);
    if isempty(G_cell) || length(G_cell) ~= length(mu)
        return
    end
    
    % Step 3: Check generator counts
    min_count = 2^(d-1);
    valid_gens = {};
    for i = 1:length(G_cell)
        if mu(i) >= min_count
            valid_gens{end+1} = G_cell{i};
        else
            return
        end
    end
    
    % Step 4: Verify via Minkowski sum
    if ~verifyGenerators(vertices, valid_gens)
        return
    end
    
    % Construct proper generator matrix
    isZono = true;
    generators = zeros(d, length(valid_gens));
    for i = 1:length(valid_gens)
        generators(:, i) = valid_gens{i};
    end
    generators = normalizeGenerators(generators);

    generators = 0.5 * generators;
end

%% Helper functions

function edges = getUniqueEdges(K, vertices)
% Extract unique edges from convex hull facets
    edges = {};
    edge_set = containers.Map();
    
    for i = 1:size(K, 1)
        facet = K(i, :);
        n_vertices = length(facet);
        
        for j = 1:n_vertices
            v1 = facet(j);
            v2 = facet(mod(j, n_vertices) + 1);
            
            % Ensure consistent ordering to avoid duplicate edges
            if v1 > v2
                [v1, v2] = deal(v2, v1);
            end
            
            key = sprintf('%d-%d', v1, v2);
            if ~isKey(edge_set, key)
                edge_set(key) = true;
                vec = vertices(v2, :) - vertices(v1, :);
                if norm(vec) > 1e-10 % Ignore zero-length edges
                    edges{end+1} = struct('vector', vec, 'length', norm(vec));
                end
            end
        end
    end
end

function [G_cell, mu] = processEdges(edges, dim)
% Process edges to find candidate generators (Algorithm 3 lines 1-14)
    G_cell = {};
    mu = [];
    directions = zeros(length(edges), dim);
    
    % First pass: collect all directions
    for i = 1:length(edges)
        [dir, vec] = normalizeDirection(edges{i}.vector);
        directions(i,:) = dir;
    end
    
    % Find unique directions
    [unique_dirs, ~, ic] = unique(round(directions*1e8)/1e8, 'rows', 'stable');
    
    % Process each unique direction
    for k = 1:size(unique_dirs, 1)
        idx = find(ic == k);
        lengths = arrayfun(@(x) edges{x}.length, idx);
        
        % Check all edges in this direction have same length
        if max(lengths) - min(lengths) > 1e-6
            G_cell = [];
            mu = [];
            return;
        end
        
        % Store the first vector and count
        G_cell{end+1} = edges{idx(1)}.vector';
        mu(end+1) = length(idx);
    end
end

function valid_gens = checkGeneratorCounts(G_cell, mu, min_count)
% Check each generator appears in enough edges
    valid_gens = {};
    for i = 1:length(G_cell)
        if mu(i) >= min_count
            valid_gens{end+1} = G_cell{i};
        else
            valid_gens = [];
            return;
        end
    end
end

function success = verifyGenerators(vertices, generators)
% Verify generators via Minkowski sum check
    success = false;
    n_vertices = size(vertices, 1);
    
    % Get original vertex count (using convex hull)
    try
        [~, V1] = convhulln(vertices);
        orig_count = size(V1, 1);
    catch
        return;
    end
    
    % Check each generator
    for i = 1:length(generators)
        s = generators{i}';
        
        % Compute P + s (Minkowski sum)
        V_plus = [vertices; vertices + repmat(s, n_vertices, 1)];
        
        try
            [~, V2] = convhulln(V_plus);
            if size(V2, 1) ~= orig_count
                return;
            end
        catch
            return;
        end
    end
    
    success = true;
end

function [dir, vec] = normalizeDirection(v)
% Normalize direction vector: first nonzero coordinate positive
    v = v(:)'; % Ensure row vector
    idx = find(abs(v) > 1e-12, 1);
    if isempty(idx)
        dir = v;
        vec = v;
        return;
    end
    
    if v(idx) < 0
        vec = -v;
    else
        vec = v;
    end
    dir = vec / norm(vec);
end

function generators = normalizeGenerators(generators)
% Ensure all generators have first non-zero element positive
    for i = 1:size(generators, 2)
        g = generators(:, i);
        idx = find(abs(g) > 1e-12, 1);
        if ~isempty(idx) && g(idx) < 0
            generators(:, i) = -g;
        end
    end
end