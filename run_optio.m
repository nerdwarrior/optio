clear

%% Create welcome page
msg = questdlg(["Welcome to opt.io!"; "";...
    "This program will help you solve basic topology optimization problems. The steps are as follows:";...
    "1) Select an existing material, or input custom material settings";...
    "2) Define your design space, load, and constraint"; "3) Define your optimization settings"; "";...
    "For help, please refer to README.md.";...
    "First, let's select your material!"],...
    "opt.io", "Proceed", "Exit", "README", "Proceed");

% Select options based on user input (proceed only if "Proceed" selected)
switch msg
    case {'Exit', ''}
        clear
        return
    case 'README'
        open("README.md")
        clear
        return
end

%% Initialize poisson (to check material loop below)
poisson = 0;

% Ask for material input until valid poisson is provided
while poisson == 0 
    
    % Create list of selections (only one can be selected at once)
    [indx,tf] = listdlg('PromptString',{"Select your material.",''},...
        'SelectionMode','single','ListString',...
        {'Aluminum AlSi10Mg', 'Stainless Steel 316L', 'Titanium 6Al-4V', 'Custom'});
    
    % Catch cases where user exits
    if tf == 0
        clear
        return
    end
    
    % If user didn't exit, assign material properties accordingly    
    switch indx
        case 1
            
            % Source: https://onlinelibrary.wiley.com/doi/abs/10.1002/mawe.201800233
            poisson = 0.35;
            modulus = 67; % GPa
        case 2
            
            % Source: https://super-metals.com/wp-content/uploads/2015/03/SS-316.pdf
            poisson = 0.27;
            modulus = 195; % GPa
        case 3
            
            % Source: https://en.wikipedia.org/wiki/Ti-6Al-4V
            poisson = 0.34;
            modulus = 108; % GPa
        
        % Case for custom material input
        case 4
            
            % Define input dialog values
            prompt = {'Poisson"s ratio:','Young"s modulus (GPa):'};
            dlgtitle = 'Custom material';
            boxdims = [1 35];
            definput = {'0.3','100'};
            
            % Store input dialog inputs in answer
            answer = str2double(inputdlg(prompt,dlgtitle,boxdims,definput));
            
            % Catches empty, NaN, non-positive, and non-numerical values
            if length(answer(answer > 0)) ~= 2
                
                % Return error message before looping
                errmsg('Error: please input a valid Poisson"s ratio and Young"s modulus.');
            else
                
                % Store values
                poisson = answer(1);
                modulus = answer(2);
            end
    end
end

%% Ask for design space dimensinos, load, constraint
dims = 6;
proceed = false;

while proceed == false
    
    % Define input dialog values
    prompt = {'Load x indices:', 'Load y indices:', 'Load z indices:',...
        newline + "Load indices (comma-separated, i.e., x, y, z):",...
        'Constraint indices (comma-separated, i.e., x, y, z):'};
    dlgtitle = 'Design space, load, & constraint';
    boxdims = [1 35];
    definput = ["24","12","12","25, 1:13, 13","1, 1:13, 1:13"];
    
    % Store input dialog inputs in answer{}
    answer = inputdlg(prompt,dlgtitle,boxdims,definput);
    
    % Catch whether user exited
    if length(answer) ~= dims / 3 + 3
        clear
        return
    end
    
    % For convenience and readability, let's extract the answers
    elements = str2double(answer(1:3));
    rangeText = convertCharsToStrings(answer{4}) + "," + convertCharsToStrings(answer{5});
    
    % Catch cases where the elements are not integers >= 10
    if length(elements(elements >= 10 & isint(elements))) ~= 3
        errmsg('Error: please input a valid number of elements (integer, >= 10)');
        continue
    end
    
    % Assign design space dimensions
    nelx = elements(1);
    nely = elements(2);
    nelz = elements(3);
    
    % Build array of range maximums
    maxes = [nelx + 1, nely + 1, nelz + 1];
    
    % Check ranges and store them in a cell array
    ranges = checkRange(rangeText, maxes, dims);
    
    % Change proceed variable
    proceed = ranges{1};
    
    % Cut off proceed variable from ranges{} (now stores 3 sets of load indices,
    % followed by 3 sets of constraint indices)
    ranges = ranges(2:end);
end

%% Ask for optimization settings input until valid inputs provided

while true
    
    % Define input dialog values
    prompt = {'Target volume fraction:', 'Penalty factor:',...
        'Filter setting (1 for density only, 2 for density + projection, 3 for density + projection + eta optimization):',...
        'Boundary condition (N for Neumann, D for Dirichlect',...
        'Maximum number of iterations (min. 20):'};
    dlgtitle = 'Optimization settings';
    boxdims = [1 35];
    definput = {'0.12','3','1','N','100'};
    
    % Store input dialog inputs in answer
    answer = inputdlg(prompt,dlgtitle,boxdims,definput);
    
    % Catch whether user exited
    if length(answer) ~= 5
        clear
        return
    end
        
    % For convenience and readability, let's extract the answers
    volfrac = str2double(answer(1));
    penalty = str2double(answer(2));
    filter = str2double(answer(3));
    bc = upper(answer{4});
    maxit = str2double(answer(5));

    % Catch cases where volfrac is not between 0 and 1 (non-inclusive)
    if volfrac <= 0 || volfrac >= 1
        errmsg('Error: please input a valid volume fraction (decimal between 0 and 1)');
    
    % Catch cases where penalty is not positive integer less than five
    elseif ~ismember(penalty,1:5)
        errmsg('Error: please input a valid penalty factor (integer, 1 through 5)');
    
    % Catch cases where filter is not in [1,2,3]
    elseif ~ismember(filter, [1,2,3])
        errmsg('Error: please input a valid filter (1, 2, or 3)');
    
    % Catch cases where bc is not 'N' or 'D'
    elseif ~strcmp(bc,'N') && ~strcmp(bc,'D')
        errmsg('Error: please input a valid boundary condition (N or D)');
    
    % Catch cases where maxit is not an integer >= 20
    elseif ~isint(maxit) || maxit < 20
        errmsg('Error: please input a valid maximum number of iterations (integer, >= 20)');
 
   % If no errors, break
    else
        break
    end
end

%% Initialize msg to control loop below
msg = 'View settings';

% Keep user on confirmation page until optimization starts or exited
while strcmp(msg,'View settings')
    
    % Confirm that topology optimization is about to start
    msg = questdlg(["You're good to go!", '', 'Select "Optimize" to begin' ...
        'running the algorithm. Exit above to cancel.'],...
        "Confirmation", "Optimize", "View settings", "Exit", "Optimize");
    
    % Select options based on user input (proceed only if "Proceed" selected)
    switch msg
        case {'Exit', ''}
            clear
            return
        case 'View settings'
            rangeText = split(rangeText,",");
            viewSettings(rangeText, modulus, poisson, nelx, nely, nelz, volfrac, penalty, filter, bc, maxit);         
    end
end

%% Try optimization & give options to view/export results
try
    
    % Save mesh output from optimization
    tr = optio(ranges, modulus, poisson, nelx, nely, nelz, volfrac, penalty, filter, bc, maxit);
    
    % View optimization results
    results(tr);
catch ME
    errmsg('Optimization failed. Please check inputs (esp. loads/constraints).');
    clear
    return
end

clear