% Inputs modified to take material parameters from user (E0 is Young's
% modulus of solid, Emin is Young's modulus of void, nu is Poisson's ratio)
function tr = optio(ranges,E0,nu,nelx,nely,nelz,volfrac,penal,ft,ftBC,maxit)
%% -----------------------------------------------------CUSTOM) ASSUMPTIONS
% Modification to reduce necessary function inputs (simplifying usage)

% Young's modulus of void
Emin = 1e-9;

% Allowable amount to move each node per iteration
move = 0.2;

% eta used in projection
eta = 0.5;

% beta (sharpness) used in projection
beta = 1;

% Filter radius
rmin = sqrt(3);
%% -----------------------------------------PRE. 1) CONTINUATION PARAMETERS
% Continuation scheme on some variable var defined as a data structure:
% varCont = {istart, maxVar, isteps, deltaVar};
% Continuation (changing the variable as you go) starts when loop=istart,
% changing by deltaVar each isteps until maxVar is reached

% Define continuation scheme for the penalty factor (steers densities to
% all-solid or all-white by raising "grey" values between densities of 0
% and 1 to this power); standard is 3 (for all topology optimization)
penalCnt = { 1, 1, 25, 0.25 };

% Define continuation scheme for the sharpness factor (smoothing) of the 
% relaxed Heaviside projection used:
% H(˜xe, η, β) = (tanh(βη) + tanh(β(˜xe − η)))/(tanh(βη) + tanh(β(1 − η)))
% Where the projection takes some intermediate density ˜xe (obtained after
% "filtering"; see block starting line 214) and maps it to a physical
% density (0 or 1; material exists or does not exist)
betaCnt  = { 1, 1, 25,    2 };

% Define fitler boundary conditions based on function input ('N' for
% zero-Neumann, 'D' (i.e., the "else" case) for zero-Dirichlect)
% Neumann specifies values of the solution's derivative along the domain
% boundary; Dirichlect specifies the value of the solution itself
if ftBC == 'N', bcF = 'symmetric'; else, bcF = 0; end       


%% -----------------------------------------PRE. 2) DISCRETIZATION FEATURES
% Get the total number of elements by multiplying the number of elements in
% the x, y, z directions (as specified by function inputs)
nEl = nelx * nely * nelz;

% reshape() returns a matrix with the size specified by the last three 
% arguments. The elements are taken column-wise from the first input
% argument.
% The below line of code creates a (1 + nely) x (1 + nelz) x (1 + nelx)
% matrix of "nodes", with the matrix values equalling the node numbers (as
% per the first argument). 
nodeNrs = int32( reshape( 1 : ( 1 + nelx ) * ( 1 + nely ) * ( 1 + nelz ), ...
    1 + nely, 1 + nelz, 1 + nelx ) );                                      

% Build towards a nEl x 24 (24 is the number of degrees of freedom of each
% node) connectivity matrix cMat. For each of the nEl local nodes, cMat
% denotes the global DoF that each of the 24 local DoFs targets (i.e., the
% DoFs of each node affects how the overall structure moves).
% This is done by first defining a column vector cVec, then using vector
% addition to get the full matrix. Row and column vector addition yields
% matrix where (i,j) element is equal to a(j) + b(i).
cVec = reshape( 3 * nodeNrs( 1 : nely, 1 : nelz, 1 : nelx ) + 1, nEl, 1 );
cMat = cVec+int32( [0,1,2,3*(nely+1)*(nelz+1)+[0,1,2,-3,-2,-1],-3,-2,-1,3*(nely+...
   1)+[0,1,2],3*(nely+1)*(nelz+2)+[0,1,2,-3,-2,-1],3*(nely+1)+[-3,-2,-1]]);

% Define the total number of global DoFs across the entire discretization
% (3 per node)
nDof = ( 1 + nely ) * ( 1 + nelz ) * ( 1 + nelx ) * 3;

% The following is an efficient way to construct row and column index
% vectors (iK, jK) that map each value of sK (line 319) to a global
% location using some clever symmetries of matrices (a little beyond
% scope) and cMat.
% deal() sets the left-hand-side variables to the corresponding inputs
% (e.g., [B1, B2] = deal(A) sets B1 = B2 = A, and [B1, B2] = deal(A1, A2)
% sets B1 = A1, B2 = A2). 
[ sI, sII ] = deal( [ ] );

% The function iterates across each of the 24 local DoFs to string together
% 24-item-long "chunks" of matrices together in sI and sII.
% cat() concatenates arrays along the dimension indicated by the first
% argument. repmat() builds arrays by repeating the first argument, using the
% dimensions indicated by the following arguments.
for j = 1 : 24
    sI = cat( 2, sI, j : 24 );
    sII = cat( 2, sII, repmat( j, 1, 24 - j + 1 ) );
end

% We then distribute those values to iK and jK via deal().
[ iK , jK ] = deal( cMat( :,  sI )', cMat( :, sII )' );

% This line simply sorts the indices into one matrix Iar, such that 
% K=sparse(iK, jK, sK) yields a lower-triangular matrix.
Iar = sort( [ iK( : ), jK( : ) ], 2, 'descend' ); clear iK jK     

% Define a fixed matrix based on physics (details of which are beyond the
% scope of these comments) that denotes the stiffness of an element with
% unit Young's modulus (which is specified as a function input, as per my
% modifications). It is based on the Poisson's ratio nu, another input of
% the function (again, as per my modifications). 
Ke = 1/(1+nu)/(2*nu-1)/144 *( [ -32;-6;-6;8;6;6;10;6;3;-4;-6;-3;-4;-3;-6;10;...
    3;6;8;3;3;4;-3;-3; -32;-6;-6;-4;-3;6;10;3;6;8;6;-3;-4;-6;-3;4;-3;3;8;3;...
    3;10;6;-32;-6;-3;-4;-3;-3;4;-3;-6;-4;6;6;8;6;3;10;3;3;8;3;6;10;-32;6;6;...
    -4;6;3;10;-6;-3;10;-3;-6;-4;3;6;4;3;3;8;-3;-3;-32;-6;-6;8;6;-6;10;3;3;4;...
    -3;3;-4;-6;-3;10;6;-3;8;3;-32;3;-6;-4;3;-3;4;-6;3;10;-6;6;8;-3;6;10;-3;...
    3;8;-32;-6;6;8;6;-6;8;3;-3;4;-3;3;-4;-3;6;10;3;-6;-32;6;-6;-4;3;3;8;-3;...
    3;10;-6;-3;-4;6;-3;4;3;-32;6;3;-4;-3;-3;8;-3;-6;10;-6;-6;8;-6;-3;10;-32;...
    6;-6;4;3;-3;8;-3;3;10;-3;6;-4;3;-6;-32;6;-3;10;-6;-3;8;-3;3;4;3;3;-4;6;...
    -32;3;-6;10;3;-3;8;6;-3;10;6;-6;8;-32;-6;6;8;6;-6;10;6;-3;-4;-6;3;-32;6;...
    -6;-4;3;6;10;-3;6;8;-6;-32;6;3;-4;3;3;4;3;6;-4;-32;6;-6;-4;6;-3;10;-6;3;...
    -32;6;-6;8;-6;-6;10;-3;-32;-3;6;-4;-3;3;4;-32;-6;-6;8;6;6;-32;-6;-6;-4;...
    -3;-32;-6;-3;-4;-32;6;6;-32;-6;-32]+nu*[ 48;0;0;0;-24;-24;-12;0;-12;0;...
    24;0;0;0;24;-12;-12;0;-12;0;0;-12;12;12;48;0;24;0;0;0;-12;-12;-24;0;-24;...
    0;0;24;12;-12;12;0;-12;0;-12;-12;0;48;24;0;0;12;12;-12;0;24;0;-24;-24;0;...
    0;-12;-12;0;0;-12;-12;0;-12;48;0;0;0;-24;0;-12;0;12;-12;12;0;0;0;-24;...
    -12;-12;-12;-12;0;0;48;0;24;0;-24;0;-12;-12;-12;-12;12;0;0;24;12;-12;0;...
    0;-12;0;48;0;24;0;-12;12;-12;0;-12;-12;24;-24;0;12;0;-12;0;0;-12;48;0;0;...
    0;-24;24;-12;0;0;-12;12;-12;0;0;-24;-12;-12;0;48;0;24;0;0;0;-12;0;-12;...
    -12;0;0;0;-24;12;-12;-12;48;-24;0;0;0;0;-12;12;0;-12;24;24;0;0;12;-12;...
    48;0;0;-12;-12;12;-12;0;0;-12;12;0;0;0;24;48;0;12;-12;0;0;-12;0;-12;-12;...
    -12;0;0;-24;48;-12;0;-12;0;0;-12;0;12;-12;-24;24;0;48;0;0;0;-24;24;-12;...
    0;12;0;24;0;48;0;24;0;0;0;-12;12;-24;0;24;48;-24;0;0;-12;-12;-12;0;-24;...
    0;48;0;0;0;-24;0;-12;0;-12;48;0;24;0;24;0;-12;12;48;0;-24;0;12;-12;-12;...
    48;0;0;0;-24;-24;48;0;24;0;0;48;24;0;0;48;0;0;48;0;48 ] );

% Set the items at lower-triangular indices (note that Ke0 is a vector, but
% MATLAB can index matrices linearly) to equal the corresponding element of
% the transposed elemental stiffness matrix.
% ones() builds an array of ones, and tril() returns the lower-triangular
% indices. 
Ke0( tril( ones( 24 ) ) == 1 ) = Ke';

% Take Ke0 and reshape it into a 24 x 24 square matrix
Ke0 = reshape( Ke0, 24, 24 );

% Recover the full elemental stiffness operator (matrix)
% diag(n) returns a square diagonal matrix where the elements of the
% n are along the diagonal. 
Ke0 = Ke0 + Ke0' - diag( diag( Ke0 ) );

%% -----------------------------PRE. 3) LOADS, SUPPORTS AND PASSIVE DOMAINS

% Define which DoFs have loads applied (here, the nodeNrs() indices select
% the lowest row of points that are furthest out, modeling the cantilever
% problem). The 3 selects the third global DoF across this edge.
lcDof = 3 * nodeNrs(ranges{2}, ranges{3}, ranges{1});

% This code is an addition to enable more flexible implementation.
% Regardless of the indices inputted above, this line will turn it into a
% column vector (as is required when we use fsparse() later).
lcDof = lcDof(:);

% Define which DoFs are constrained (cannot move). Here, the entire
% front face is selected (square region created by (nely + 1) * (nelz +
% 1)). The 3 selects for all three global DoFs across this face.
% This code is modified from the source code to enable more flexible
% implementation (by allowing users to input x, y, and z indices directly).
% We simply break out each index into its three DoFs, concatenate them
% together with cat(), and sort into ascending order with sort().
fixed = 3 * nodeNrs(ranges{5}, ranges{6}, ranges{4});
fixed = sort(cat(2, fixed(:)', fixed(:)' - 1, fixed(:)' - 2));

% Declare passive solid elements (full density, but not actively being
% manipulated) and passive void elements (0 density, not manipulated)
[ pasS, pasV ] = deal( [], [] );

% Compute the load vector
% S = fsparse(II,JJ,SS,[M N Nzmax]) creates a sparse matrix S from triplet
% data (II,JJ,SS). The inputs are all matrices with either matching or
% singleton dimensions according to rules:
% If II is IM-by-IN, JJ JM-by-JN and SS SM-by-SN, then it is required that
% (1) IN = JN or 1, (2) JM = IM or 1, (3) SM = IM or 1 and (4) SN = JN or
% 1. Last argument specifies dimensions of output (corresponds to global
% DoFS).
% N.B.: A sparse matrix is an efficient way to store/work with matrices that are
% mostly zero/empty. Instead of storing all the values, it stores the
% non-zero/non-empty values and their indices. 
F = fsparse( lcDof, 1, -ones(length(lcDof),1), [ nDof, 1 ] );

% Set free DoFs as the full DoFs minus the fixed DoFs (setdiff() simply
% takes the set difference)
free = setdiff( 1 : nDof, fixed );

% Set active variables as the nodes that are not passive (union() simply
% takes the set union)
act = setdiff( ( 1 : nEl )', union( pasS, pasV ) );                        

%% -------------------------------------- PRE. 4) DEFINE IMPLICIT FUNCTIONS
% Define functions for operations during the main topology optimization
% loop
% @ creates a function handle to let us use these as functions later

% Define the Heaviside projection (see line 32); projection is only enabled
% when ft=2 or ft=3 (ft=1 only uses density filtering). 
prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v(:)-eta)))./...
    (tanh(beta*eta)+tanh(beta*(1-eta)));

% Define projection eta-derivative (for Heaviside)
deta = @(v,eta,beta) - beta * csch( beta ) .* sech( beta * ( v( : ) - eta ) ).^2 .* ...
    sinh( v( : ) * beta ) .* sinh( ( 1 - v( : ) ) * beta );

% Define projection x-derivative (to get intermediate density later; see line 32)
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));

% Apply continuation (based on description at line 22) 
cnt = @(v,vCnt,l) v+(l>=vCnt{1}).*(v<vCnt{2}).*(mod(l,vCnt{3})==0).*vCnt{4};

%% ------------------------------------------------- PRE. 5) PREPARE FILTER

% Define filtering grids (3D matrices) based on the inputted filter radius
% (rmin)
% meshgrid(x,y,z) returns the coordinates for each node in an x*y*z grid
[dy,dz,dx]=meshgrid(-ceil(rmin)+1:ceil(rmin)-1,...
    -ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1 );

% Define non-zero component of filter kernel
h = max( 0, rmin - sqrt( dx.^2 + dy.^2 + dz.^2 ) );

% Define matrix of weights to use for filtering (stays constant throughout
% the optimization loop, so defined first)
% imfilter(A,B,C) filters array A with filter B using option C (e.g.,
% blurring an image A with a motion blur filter B)
Hs = imfilter( ones( nely, nelz, nelx ), h, bcF );

% Copy filter into new variable, which we can update as we go to alter the
% sensitivity during the optimization loop
dHs = Hs;

%% ----------------------- PRE. 6) ALLOCATE AND INITIALIZE OTHER PARAMETERS
% Initialize variables (zeros() builds an array of zeros)
[ x, dsK, dV ] = deal( zeros( nEl, 1 ) );

% Define volume derivative (volfrac is a user input for max volume
% fraction from optimization); this value is constant and encapuslates
% "sensitivity" (see optimization loop)
dV( act, 1 ) = 1/nEl/volfrac;

% Apply volume fraction to the active set
x( act ) = ( volfrac*( nEl - length(pasV) ) - length(pasS) )/length( act );

% Set a density of 1 on solid passive set (as per line 169)
x( pasS ) = 1;

% Create variables to keep track of the physical density (xPhys), the old x
% value (xOld), the change in x (ch), the iteration counter (loop), and
% displacement (U). At each iteration, U is computed by solving Ke0U = F
% (i.e., using the stiffness and load matrices)
[ xPhys, xOld, ch, loop, U ] = deal( x, 1, 1, 0, zeros( nDof, 1 ) );

% Modification to start the preview in a new figure window
figure('Name','Optimization preview','NumberTitle','off','Visible','on');
%% ================================================ START OPTIMIZATION LOOP

% Optimization continues while results have not converged and the maximum
% number of iterations has not been reached
while ch > 1e-6 && loop < maxit
  
  % Increment the iteration counter
  loop = loop + 1;
  
  % ----------- RL. 1) COMPUTE PHYSICAL DENSITY FIELD (AND ETA IF PROJECT.)
  
  % Create the filtered field (corresponds to variable discussed in line
  % 34)
  xTilde = imfilter( reshape( x, nely, nelz, nelx ), h, bcF ) ./ Hs;
  
  % Reshape the filtered field to a column vector and store as physical
  % densities
  xPhys( act ) = xTilde( act );
  
  % Apply projection if specified by the user
  if ft > 1
      
      % Compute the optimal eta in the Heaviside equation (line 34) via the
      % Newton method (stepping closer to the soln using the last-computed
      % value and the derivative) 
      
      % Define a conditional volume fraction for Newton steps (conditional
      % on ft being 3)
      f = ( mean( prj( xPhys, eta, beta ) ) - volfrac )  * (ft == 3);
      
      % Run Newton steps until settled upon optimal eta within error of
      % 1e-6
      while abs( f ) > 1e-6
          
          % Take a Newton step
          eta = eta - f / mean( deta( xPhys, eta, beta ) );
          
          % Update the volume fraction
          f = mean( prj( xPhys, eta, beta ) ) - volfrac;
      end
      
      % Modify the filter weights based on the current projection (uses
      % Newton-stepped eta if ft==3, or takes constant eta otherwise due to
      % ewton steps not being performed)
      dHs = Hs ./ reshape( dprj( xPhys, eta, beta ), nely, nelz, nelx );
      
      % Update the densities to the projected field
      xPhys = prj( xPhys, eta, beta );                                     
  end
  
  % Update the change in x (when unspecified as below, norm() computes the 2-norm of the inputted
  % matrix)
  ch = norm( xPhys - xOld ) ./ nEl;
  
  % Store the densities into the old densities to compute future ch
  xOld = xPhys;
  
  % -------------------------- RL. 2) SETUP AND SOLVE EQUILIBRIUM EQUATIONS
  
  % Create stiffness interpolation based on current densities and the
  % inputted Young's moduli
  sK = ( Emin + xPhys.^penal * ( E0 - Emin ) );
  
  % Compute sK's derivative (on the active nodes)
  dsK( act ) = -penal * ( E0 - Emin ) * xPhys( act ) .^ ( penal - 1 );
  
  % Apply elemental stiffness matrix to the stiffness interpolation
  sK = reshape( Ke( : ) * sK', length( Ke ) * nEl, 1 );
  
  % Obtain the global stiffness matrix using iK and jK values (stored in
  % Iar) and stiffness interpolation
  K = fsparse( Iar( :, 1 ), Iar( :, 2 ), sK, [ nDof, nDof ] );
  
  % Use the Cholesky solver on the free (non-constrained) nodes to obtain
  % the full global stiffness operator
  % chol(A) calculates the Cholesky factor of A using the diagonal and
  % upper triangle of A. chol(A, 'lower') uses the diagonal and lower
  % triangle. A must be positive definite.
  L = chol( K( free, free ), 'lower' );
  
  % Update the displacements of the free (non-constrained) nodes by
  % applying the global stiffness operator
  U( free ) = L' \ ( L \ F( free ) );
  
  % ------------------------------------------ RL. 3) COMPUTE SENSITIVITIES
  
  % Compute derivative (sensitivity) of compliance (objective) by applying elemental
  % stiffness matrix to the displacements, selected using the connectivity
  % matrix. Multiply element-wise with the stiffness interpolation.
  % sum(A,2) returns the sum of the elements of A along the second
  % dimension
  dc = dsK .* sum( ( U( cMat ) * Ke0 ) .* U( cMat ), 2 );
  
  % Apply filter to the compliance sensitivity 
  dc = imfilter( reshape( dc, nely, nelz, nelx ) ./ dHs, h, bcF );
  
  % Apply filter to the volume constraint sensitivity
  dV0 = imfilter( reshape( dV, nely, nelz, nelx ) ./ dHs, h, bcF );
  
  % ----------------- RL. 4) UPDATE DESIGN VARIABLES AND APPLY CONTINUATION
 
  % Uses optimality criteria method to update design variables.
  
  % Get active nodes
  xT = x( act );
  
  % Define current upper (xU) and lower (xL) bounds for moving x, based on
  % user input "move"
  [ xU, xL ] = deal( xT + move, xT - move );
  
  % constant part in resizing rule 
  ocP = xT .* sqrt( - dc( act ) ./ dV0( act ) );
  
  % Estimate the square root of the lagrange multiplier, which is used
  % instead of the design variables themselves for efficiency
  l = [ 0, mean( ocP ) / volfrac ];
  
  % Resize the optimality criteria while lagrange's error is high
  while ( l( 2 ) - l( 1 ) ) / ( l( 2 ) + l( 1 ) ) > 1e-4                   
      
      % Get the lagrange multiplier??
      lmid = 0.5 * ( l( 1 ) + l( 2 ) );
      
      % Update x based on the resizing rules, current lagrange multiplier,
      % and upper/lower limits (specific math beyond scope)
      x( act ) = max( max( min( min( ocP / lmid, xU ), 1 ), xL ), 0 );
      
      % Update the lagrange multiplier
      if mean( x ) > volfrac, l( 1 ) = lmid; else, l( 2 ) = lmid; end
  end
  
  % Apply continuation to the penalty and beta based on definition at line
  % 212
  [penal,beta] = deal(cnt(penal,penalCnt,loop), cnt(beta,betaCnt,loop)); 
  
  % ---------------------------------------------------- RL. 5) PLOT DESIGN
  % Rearrange variables into format for creating isosurfaces
  % shiftdim(A, n) moves the dimension sizes of A by n positions (wrapping
  % around at extremes). For example, if A is 4x2x3, shiftdim(A, 1) yields
  % a 3x4x2 matrix.
  isovals = shiftdim( reshape( xPhys, nely, nelz, nelx ), 2 );
  
  % Apply some smoothing for isosurfaces (smooth3() is for 3D data
  % specifically; 'box' specifies the filter, and 1 specifies the filter size)
  isovals = smooth3( isovals, 'box', 1 );
  
  % Create surface patch over isosurface to display (based on points
  % contained in isovals)
  % isosurface(V, b) builds an isosurface using volume data contained in V,
  % connecting points equal to b with edges (just as contour lines connect
  % points of equal elevation).
  % patch() builds polygonal data by filling in the polygonal regions,
  % returned in this case from isosurface()
  patch(isosurface(isovals, .5),'FaceColor',[0.5,0.5,0.5],'EdgeColor','w');
  
  % Create caps to display (creating a closed surface)
  % isocaps() is very similar to isosurface() but it builds end-cap
  % geometry to create a closed surface.
  patch(isocaps(isovals, .5),'FaceColor',[0.5,0.5,0.5],'EdgeColor','w');
  
  % Display the mesh at the given view angle
  % drawnow updates all figures, view() sets the camera position for the
  % figure, and cla() clears the existing graphics
  drawnow; view( [ 145, 25 ] ); axis equal tight off; cla();
  
  % Modification to reset figure values now, such that even if the user
  % closes the window, it will pop back up with the correct name
  set(gcf,'Visible','on','Name','Optimization preview','NumberTitle','off');
end

%% --------------------------------------------------- Export configuration
% (Modification to output mesh instead of a video)

% Extract faces and vertices from surface and surface caps
[fs, vs] = isosurface(isovals, .5);
[fc, vc] = isocaps(isovals, .5);

% Stitch together two sets into one (for a coherent surface)
f = [fs; fc + length(vs(:,1))];
v = [vs; vc];

% Change the order of points in faces data to flip the mesh normals
f = f(:,[3 2 1]);

% Create mesh
tr = triangulation(f,v);

% Close file preview
close();
end
