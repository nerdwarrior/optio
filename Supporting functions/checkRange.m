function output = checkRange(rangeText, maxes, dims)
    if mod(dims,3) ~= 0 || ~prod(isint(maxes)) || length(maxes) ~= 3
        errmsg("checkRange usage error (check inputs)")
        return
    end
    
    output = cell(1,1 + dims);
    output{1} = true;
    
    % Get ranges from inputs
    ranges = split(rangeText, ",");
    for i = 1:dims

        try
            rng = str2num(ranges(i));
        catch ME
            errmsg('Error: please input valid ranges');
            output{1} = false;
            return
        end

        % Abort if not valid range (i - 3 * floor(i/4) lets us loop through
        % maxes with our indices)
        if isempty(rng) || max(rng) > maxes(i - 3 * floor(i/4)) || min(rng) < 1
            errmsg('Error: please input valid ranges');
            output{1} = false;
            return
        end

        % Assign to output{} if OK
        output{i + 1} = rng;
    end
end