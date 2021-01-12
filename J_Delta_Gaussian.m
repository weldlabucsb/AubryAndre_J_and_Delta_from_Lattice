function [J,Delta] = J_Delta_Gaussian(s1,s2,lambda1,lambda2)
%Calculates J and Delta of the 1D Aubry-Andre Tight Binding Hamiltonian based
%on the assumption that the wannier states are gaussians localized on each
%site (i.e. each well is a harmonic oscillator).
%
%   Based on [Michele Modugno 2009 New J. Phys. 11 033023]
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
%   Outputs J and Delta are both given in units of the recoil energy of
%   the primary lattice laser.

beta = lambda1./lambda2;

J = (4/sqrt(pi))*s1.^(0.75).*exp(-2*sqrt(s1));

Delta = 0.5*(s2.*beta.^2).*exp(-beta.^2./sqrt(s1));
end

