%% Sierpinski_ball
%
% Function to compute, display, and save a Sierpinski ball
% at any iteration number / depth level.
%
% Author : nicolas.douillet9 (at) gmail.com, 2019-2024.
%
%% Syntax
%
% Sierpinski_ball(nb_it);
%
% Sierpinski_ball(nb_it, option_display);
%
% [V, T] = Sierpinski_ball(nb_it, option_display);
%
%% Description
%
% Sierpinski_ball(nb_it) computes and display the nb_it Sierpinski
% ball of unit radius, based on the regular octahedron.
%
% Sierpinski_ball(nb_it, option_display) displays it when
% option_display is set to logical *true/1 (default), and doesn't
% when it is set to  logical false/0.
%
% [V, T] = Sierpinski_ball(nb_it, option_display) saves the resulting
% vertex coordinates in the array V, and the triangulation in the array T.
%
%% See also
%
% <https://fr.mathworks.com/matlabcentral/fileexchange/73285-sierpinski-sphere-spherpinski Spherpinski> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/79152-sierpinski-octahedron Sierpinski_octahedron> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/82595-sierpinski-cube Sierpinski_cube> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/79439-sierpinski-icosahedron Sierpinski_icosahedron>
% <https://fr.mathworks.com/matlabcentral/fileexchange/79439-sierpinski-cuboctahedron Sierpinski_cuboctahedron> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/158166-pi3mft PI3MFT>
% 
%% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
%
% - option_display : either logical, *true/false or numeric *1/0.
%
%% Output arguments
%
%        [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%        [ |  |  |]
%
%        [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%        [ |  |  |]
%
%% Example #1
% Computes and displays the simple Sierpinski ball at iteration 1

Sierpinski_ball(1);

%% Example #2
% Compute, displays, and saves the Sierpinski ball at iteration 3

[V,T] = Sierpinski_ball(3,true);