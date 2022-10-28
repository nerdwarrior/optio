function dummyOutput = viewSettings(rangeText, modulus, poisson, nelx, nely, nelz, volfrac, penalty, filter, bc, maxit)
    
    % Convert filter input to description
    switch filter
        case 1
            filter = "Density only";
        case 2
            filter = "Density + projection";
        case 3
            filter = "Density + projection + eta optimization";
    end
    
    % Convert bc input to description
    if strcmp(bc,'N')
        bc = "Neumann";
    elseif strcmp(bc,'D')
        bc = "Dirichlect";
    end
    
    % Create message
    msg = msgbox(["Poisson's ratio: " + poisson, "Young's modulus: " + modulus + " GPa",'',...
        "X elements: " + nelx, "Y elements: " + nely, "Z elements: " + nelz,'',...
        "Load indices (x, y, z): (" + rangeText{1} + ", " + rangeText{2} + ", " + rangeText{3} + ")",...
        "Constraint indices (x, y, z): (" + rangeText{4} + ", " + rangeText{5} + ", " + rangeText{6} + ")",'',...
        "Volume fraction: " + volfrac, "Penalty factor: " + penalty, "Filter: " + filter,...
        "Boundary condition: " + bc, "Maximum interations: " + maxit],...
    "Settings");
    uiwait(msg);
    
    dummyOutput = NaN;
end