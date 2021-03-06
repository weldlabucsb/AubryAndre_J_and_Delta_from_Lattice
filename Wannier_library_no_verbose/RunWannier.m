
% Calculates bandstructure, Wannier states, and interaction terms for bosons in an optical lattice.
% 
% REFERENCES (equation numbers are referenced in the code):
% 1. N. Marzari and D. Vanderbilt. Maximally localized generalized Wannier functions for composite energy bands.
% Phys. Rev. B 56(20), 12847-12865 Nov 1997.
% 2. I. Souza, N. Marzari and D. Vanderbilt. Maximally localized Wannier functions for entangled energy bands.
% Phys. Rev. B 65(3), 035109 Dec 2001.
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

function [lattice, recip, potential, bloch, neighbours, wannier90, manyBody] = ...
    RunWannier(filename)

% disp('-----------------------------------------------------------------------------');
% disp('MAXIMALLY LOCALISED GENERALISED WANNIER FUNCTIONS');
% disp('Richard Walters, Giovanni Cotugno, Tomi Johnson, Stephen Clark, Dieter Jaksch');
% disp('Atomic and Laser Physics, Clarendon Laboratory, University of Oxford');
% disp('-----------------------------------------------------------------------------');


% Load the RunWannier parameters
load(['./Wannier_data/Parameters/' filename]);

% --------------
% Fix any dimensional mistakes the user may have made.  
% --------------
dimension  = size(G,1);

if size(mbDensity,2) < dimension
    mbDensity = [mbDensity repmat(mbDensity(:,1),1,dimension-size(mbDensity,2))];
end

if size(plotDensity,2) < dimension
    plotDensity = [plotDensity repmat(plotDensity(:,1),1,dimension-size(plotDensity,2))];
end

if size(mbCells,1) < dimension
    mbCells = [mbCells' repmat(mbCells(1,:)',1,dimension-size(mbCells,1))]';
end
     
if size(plotCells,1) < dimension
    plotCells = [plotCells' repmat(plotCells(1,:)',1,dimension-size(plotCells,1))]';
end
% --------------

% Create an instance of a Lattice
lattice = Lattice(G);
% Redefine the reciprocal lattice basis vectors as those with minimum norm
lattice = lattice.MinimalGSet();
% Create an instance of a Potential
potential = Potential(vi, hkl, v0);

% Create directories if they don't already exist
dataPath = './Wannier_data';
folders = {'Bloch_data', 'Isolated_data', 'Composite_data', 'Many_body_data'};
for n = 1 : length(folders)
    if exist([dataPath '/', char(folders(n))], 'dir') ~= 7
        mkdir([dataPath '/', char(folders(n))]); end
end; clear folders n
% Filename for Bloch data
blochData = [dataPath '/Bloch_data/' filename '.mat'];
% Filename for isolated Wannier data
isolatedData = [dataPath '/Isolated_data/' filename '.mat'];
% Filename for composite Wannier data
compositeData = [dataPath '/Composite_data/' filename '.mat'];
% Filename for many-body data
manyBodyData = [dataPath '/Many_body_data/' filename '.mat'];

if exist(blochData, 'file') ~= 2 || strcmp(recalcBloch, 'true')
    % Create an instance of ReciprocalLattice to obtain a mesh of reciprocal lattice points
    % for calculating the single particle physics. Also determine the difference matrix,
    % and apply an energy cut-off
    recip = ReciprocalLattice(lattice.Dimension, lattice.G, Gmax);
    recip = recip.DifferenceMatrix();
    recip = recip.CutOff();
    % Set the potential coefficients using recip
    potential = potential.Coefficients(lattice, G, recip);
    % Create an instance of Bloch (calculate Bloch states and energies)
    bloch = Bloch(lattice, recip, potential, 'origin', N, numBands);
    % Create an instance of Neighbours using the Bloch quasimomenta mesh
    neighbours = Neighbours(lattice, bloch.Q.Density(1), bloch.Q.MeshInd);
    % Save band structure calculation
    save(blochData, 'lattice', 'potential', 'recip', 'bloch', 'neighbours');
%     disp(['Saved Bloch data in file ' blochData]);
    % Create an instance of Wannier90
    wannier90 = Wannier90(bloch.State, bloch.Q.MeshInd, recip, neighbours);
    % Run the wannier90 algorithm on each band in bloch
    wannier90 = wannier90.Isolated(lattice.Dimension, bloch.Q.Mesh, neighbours, epsilon, iterIso);
    % Save all the objects created thus far
    save(isolatedData, 'wannier90');
%     disp(['Saved isolated band Wannier data in file ' isolatedData]);
else
%     disp('Using Bloch data from file');
    load(blochData);
%     disp('Using isolated Wannier data from file');
    load(isolatedData);
end

% Set which bands belong to composite groups 

wannier90 = wannier90.SetGroups(groups);

if exist(compositeData, 'file') ~= 2 || strcmp(recalcComposite, 'true')
    % Load an existing compositeData file if one exists
    if strcmp(reloadComposite, 'true') && exist(compositeData, 'file') == 2
        disp('Loading composite Wannier data from file');
        load(compositeData);
    elseif strcmp(reloadComposite, 'true') && exist(compositeData, 'file') ~= 2
        disp('WARNING: issued reloadComposite but file does not exist');
    end
    if strcmp(disentangle, 'true')
        % Run the disentangling algorithm on any composite groups
        wannier90 = wannier90.Disentangle(neighbours, alphaDis, iterDis, randomise);
        wannier90 = wannier90.Isolated(lattice.Dimension, bloch.Q.Mesh, neighbours, epsilon, iterIso);
    end
    % Run the wannier90 algorithm for any composite groups
    if strcmp(lattice.Dimension, '1D') && strcmp(parallelTransport, 'true')
        wannier90 = wannier90.ParallelTransport(bloch.Q.Mesh, neighbours); end
    wannier90 = wannier90.Composite(bloch.Q.Mesh, neighbours, epsilon, iterComp);
    % Save the wannier90 instance (now for composite groups)
    save(compositeData, 'wannier90');
%     disp(['Saved composite Wannier data in file ' compositeData]);
else
    disp('Using composite Wannier data from file');
    load(compositeData);
end

if exist(manyBodyData, 'file') ~= 2 || strcmp(recalcManyBody, 'true') ...
        || strcmp(recalcComposite, 'true') || strcmp(recalcBloch, 'true')
    % Create an instance of SuperCell over which the Wannier functions will be calculated
%     mbCells = [-(N-1)/2 (N-1)/2];
    superCell = SuperCell(lattice.Dimension, lattice.R, mbDensity, mbCells, 'true');
    % Create an instance of ManyBody, setting the groups and sites
    manyBody = ManyBody(wannier90.Groups, mbSites, superCell);
    % Calculate the Wannier functions
    manyBody = manyBody.CalculateWannier(lattice, recip, bloch, wannier90, calcMode);
    % Apply the transformation that makes the Wannier functions real to the wannier90 instance
    wannier90 = wannier90.Transform(neighbours, manyBody.UPhase);
    save(compositeData, 'wannier90');
%     disp(['Saved composite Wannier data in file ' compositeData]);
    % Calculate the hopping matrix elements
    manyBody = manyBody.CalculateHopping(lattice, bloch, wannier90);
    % Calculate the interaction matrix elements
    manyBody = manyBody.CalculateInteraction();
    % Save the manyBody instance
    save(manyBodyData, 'manyBody');
%     disp(['Saved many-body data in file ' manyBodyData]);
else
    disp('Using many-body data from file');
    load(manyBodyData);
end

% disp('-------------------------------------------------');
% disp('MAXIMALLY LOCALISED GENERALISED WANNIER FUNCTIONS');
% disp('End of calculations');
% disp('-------------------------------------------------');
