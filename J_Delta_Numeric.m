function [J,Delta] = J_Delta_Numeric(s1,s2,lambda1,lambda2,showPlots)
% Based on the MLGWS code by
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Giovanni Cotugno, Tomi Johnson, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------
%
% Outputs calculated J and Delta for the Aubry-Andre Model.  You must
% specify in code the details of the numerics.  In most cases, the default
% should be accurate, but if you go to very deep lattices, you might need
% to set "Gmax" to a higher value to use higher frequency plane waves in
% the approximation of the bloch functions.  At shallow lattices, you may
% need to increase the size of the position space mesh with mbCells.  There
% may be other things that need to be tweaked as well.  See the "User
% Documentation.pdf" file of the MLGWS code for complete details on
% parameters of the calculation.
%
% Note -- outputs of the MLGWS will be dumped into a folder called
% "Wannier_data".  It's slightly too much work to suppress this, and may be
% useful for your reference.
%
%   s1 is the primary lattice depth in units of the recoil energy of a
%   photon from the primary lattice laser
%   s2 is the secondary lattice depth in units of the recoil energy of the 
%   recoil energy of a photon from the secondary lattice laser
%
%   lambda1 (lambda2) is the wavelength of the primary (secondary) lattice.
%   Units used don't matter because only beta = lambda1/lambda2 shows up
%   in calculations
%
%   showPlots is a boolean that determines if plots of the lattice 
%   potential and the wannier states are produced.  If True, plots are
%   made.
%
%   Outputs J and Delta are both given in units of the recoil energy of
%   the primary lattice laser.

    % #########################################################################################
    % INPUT PARAMETERS
    % #########################################################################################

    % ---------------------
    % MAIN INPUT PARAMETERS
    % ---------------------
    % The reciprocal lattice vectors. Lengths are given in units of wavelength, i.e. lambda = 1
    G = 4 * pi;
    % The coordinates of the potential coefficients in the reciprocal lattice, given in units of
    % the reciprocal lattice vectors, G
    hkl = [0 1 -1];
    % The ratios of the corresponding potential coefficients
    vi = [1 -1 -1];
    % Strength of the potential in units of the recoil energy, E_R
    v0 = s1;
    % The composite groups to calculate the generalised Wannier functions for.
    % Example: groups = {1,[2,3]}; Band 1 is isolated, while bands 2 and 3 form a composite group
    groups = {1,2};
    % The lattice sites at which to calculate Wannier functions, given in units of the lattice
    % vectors. Hopping and interaction parameters will be calculated between Wannier functions
    % at each site and at the origin
    mbSites = [0 1];

    % ---------------
    % FILE PARAMETERS
    % ---------------
    filename = '1D-CosineLattice'; % A filename for saving the data (make this parameter dependent for data runs)

    % ----------------
    % BLOCH PARAMETERS
    % ----------------
    recalcBloch = 'true'; % Set true to recalculate the Bloch data
    % The maximum norm of reciprocal lattice vectors included in the mesh for the Fourier
    % decompositions.
    Gmax = 100;
    % No. of real space unit cells / no. of q points in each direction of the 1st BZ.
    % The value must be even
    N = 400;

    % --------------------
    % WANNIER90 PARAMETERS
    % --------------------
    recalcComposite = 'true'; % Set true to recalculate the wannier90 optimised Bloch data (and onwards)
    reloadComposite = 'false'; % Set true to load previous wannier data (recalcComposite must be true)
    parallelTransport = 'true'; % 1D only. Set true to use the Parallel Transport algorithm
    disentangle = 'true'; % Set true to disentangle the bands before running Wannier90
    randomise = 'true'; % Set true to randomise the band indices at each q before disentangling
    iterIso = 1000; % The number of iterations for optimising isolated bands
    iterComp = 1000; % The number of iterations for optimising composite bands
    iterDis = 1000; % The number of iterations for disentangling the bands
    epsilon = 1; % Variable used in Wannier90. Set between 0 and 1
    alphaDis = 1; % Variable used in the disenatngling algorithm. Set between 0 and 1

    % ------------------------------
    % HUBBARD (MANY BODY) PARAMETERS
    % ------------------------------
    recalcManyBody = 'true'; % Set true to recalculate the Many Body data
    calcMode = 'multiprod'; % 'multiprod' is memory intensive. Use 'loop' if this causes issues 
    % mbDensity: density of real space points in each lattice direction
    mbDensity = 100;
    % mbCells: no of cells in each direction to produce the real space grid
    mbCells = [-5 5];

    % --------------
    % Plot paramters
    % --------------
    % plotCells: limits of the super-cell on which to plot the functions
    % e.g. [-1 1] produces a plot including the unit cells centred at -lambda, 0, lambda
    % plotDensity: no. of points in each real space unit cell, in each lattice direction
    plotCells = [-1 1];
    plotDensity = 100;

    % #########################################################################################
    % WANNIER CALCULATION
    % #########################################################################################

    % -------------------
    % PRE-PROCESSING CODE
    % -------------------
    % Add library folder to the Matlab path
    path('./Wannier_library_no_verbose', path);
    % Create data directories if they don't exist
    if exist('./Wannier_data', 'dir') ~= 7
        mkdir('./Wannier_data'); end
    if exist('./Wannier_data/Parameters', 'dir') ~= 7
        mkdir('./Wannier_data/Parameters'); end
    % Calculate total number of bands
    [numBands] = TotalBands(groups);
    % Save the parameters to a file. This is then loaded by RunWannier.m to calculate the
    % Wannier functions
    save(['./Wannier_data/Parameters/' filename]);

    % -----------------------------------------------------------------------------
    % CALCULATE THE BLOCH STATES, THE WANNIER FUNCTIONS, AND THE HUBBARD PARAMETERS
    % -----------------------------------------------------------------------------
    [lattice, recip, potential, bloch, neighbours, wannier90, manyBody] = RunWannier(filename);

    % #########################################################################################
    % PLOTTING
    % #########################################################################################

    % --------------------------------------------------------
    % Produce the plots of the potential and Wannier functions
    % --------------------------------------------------------
    if showPlots
        Plots(lattice, potential, recip, bloch, wannier90, manyBody, plotDensity, plotCells, calcMode)
    end
    
    %% Extracting J from the calculations of MLGWS.
    if imag(manyBody.J(1,1,2))>(1e-10)
        error('The imaginary part of J is larger than 1e-10, which probably means that the numerics are not working.  Try increasing Gmax maybe?')
    end
    assert(mbSites(2)==1,'mbSites(2) must be set to 1 in order to get J correctly. (This can be changed if you also change the element of manyBody.J that is extracted.)')
    J = real(-manyBody.J(1,1,2));
    
    
    %% Calculating Delta
    
    beta = lambda1./lambda2;
    
    % X positions for the real space mesh
    X = transpose(manyBody.SCell.Mesh);
    
    % Wannier state values on the X position mesh
    assert(mbSites(1)==0,'mbSites(1) must be set to 0 in order to get the correct Wannier state. (This can be changed if you also change the elements of manyBody.W that are extracted.)')
    W = manyBody.W(:,1,1);
    
    % Volume element (distance between mesh points in 1D)
    dVol = manyBody.SCell.dVol;
    
    % Integral in definition of Delta
    deltaIntegral = dVol*trapz( cos(4*pi*beta*X).*W.^2);
    
    Delta = 0.5*s2*beta^2*deltaIntegral;
    
end