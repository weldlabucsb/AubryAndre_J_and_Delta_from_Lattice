function J = J_Numeric_Approx(s1)
%Calculates J of the 1D Tight Binding Hamiltonian based on an approximate
%analytic form that is good to 1% in the range of 8Er to 30Er where Er is 
%the recoil energy of the lattice.
%
%   Based on [Michele Modugno 2009 New J. Phys. 11 033023] who actually is
%   just referencing [Gerbier F et al 2005 Phys. Rev. A 72 053606  (footnote 2 on page 2)] 
%
%   s1 is the (primary) lattice depth in units of the recoil energy of a
%   photon from the primary lattice laser
J = 1.43*s1.^0.98*exp(-2.07*sqrt(s1));
end

