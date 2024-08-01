function [V, T] = Sierpinski_ball(nb_it, option_display)
%% Sierpinski_ball : function to compute and display the
% Sierpinski ball at any iteration / depth level.
%
% Author : nicolas.douillet9 (at) gmail.com, 2019-2024.
%
%
% Syntax
%
% Sierpinski_ball(nb_it);
% Sierpinski_ball(nb_it, option_display);
% [V, T] = Sierpinski_ball(nb_it, option_display);
%
%
% Description
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
%
% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
% - option_display : either logical, *true/false or numeric *1/0.
%
%
% Output arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%       [ |  |  |]
%
%
% Example #1
%
% Compute and display the simple Sierpinski ball at iteration 1
%
% Sierpinski_ball(1);
%
%
% Example #2
%
% Compute and display the Sierpinski ball at iteration 3
% Store vertices and triangles.
%
% [V,T] = Sierpinski_ball(3,true);


%% Inputs parsing
assert(nargin < 3,'Too many input arguments.');

if ~nargin
    nb_it = 3;
    option_display = true;
elseif nargin > 0
    assert(isnumeric(nb_it) && nb_it == floor(nb_it) && nb_it >= 0,'nb_it parameter value must be numeric positive or null integer.');
    if nargin > 1
        assert(islogical(option_display) || isnumeric(option_display),'option_display parameter type must be either logical or numeric.');
    else
        option_display = true;
    end
end

warning('on');
if option_display && nb_it > 4    
    warning('%s facets to display ! Make sure your graphic card has enough memory.',num2str(4^(nb_it+1)));    
end
warning('off');


%% Body
% Basis vectors
I = [1 0 0];
J = [0 1 0];
K = [0 0 1];

nb_max_it = 7;
sample_step = 2^(nb_max_it-nb_it);

% Create root / mother meshed tetrahedron
[V,T] = Sierpinski_tetrahedron_iterate(I,J,K,[0 0 0],nb_it, sample_step);

% Project it into the unit ball
sphere_coord_mat = zeros(size(V,1),4); % (r, theta, phi, coeff)
sphere_coord_mat(:,1) = abs(sqrt(sum(V.^2,2))); % r
sphere_coord_mat(:,2) = abs(acos(V(:,3)./sphere_coord_mat(:,1))); % theta
sphere_coord_mat(:,3) = abs(acos(V(:,1)./(sphere_coord_mat(:,1).*sin(sphere_coord_mat(:,2))))); % phi
sphere_coord_mat(:,4) = sin(sphere_coord_mat(:,2)) .* (cos(sphere_coord_mat(:,3)) + sin(sphere_coord_mat(:,3))) + cos(sphere_coord_mat(:,2)); % multiplying coeff

idx = isnan(sphere_coord_mat(:,4));
sphere_coord_mat(idx,4) = 1;
V = V .* repmat(sphere_coord_mat(:,4),[1 3]);

% Perform one Rz 180° rotation, and one Rx 180° rotation such that the
% resulting Sierpinski sphere is based on a regular octahedron
RzV = ([-1 0 0; 0 -1 0; 0 0 1]*V')';
T = cat(1,T,T+size(V,1));
V = cat(1,V,RzV);

RxV = ([1 0 0; 0 -1 0; 0 0 -1]*V')';
T = cat(1,T,T+size(V,1));
V = cat(1,V,RxV);

% Remove duplicated vertices
[V,T] = remove_duplicated_vertices(V,T);

% % Ellipsoid option
% V(:,2) = 0.5*(1+sqrt(5))*V(:,2);
% V(:,1) = 0.25*(1+sqrt(5))^2*V(:,1);

% Remove duplicated triangles
T = unique(sort(T,2), 'rows', 'stable');

%% Display
if option_display
    
    figure;
    set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0]);
    trisurf(T,V(:,1),V(:,2),V(:,3),'EdgeColor',[0 1 0]), shading interp, hold on; % shading interp,
    colormap([0 1 0]);
    axis equal, axis tight, axis off;
    ax = gca;
    ax.Clipping = 'off';
    camlight left;
    
end

end % Sierpinski_ball


%% create_Sierpinski_ball_basis_tetrahedron subfunction
function [M, T, C] = create_Sierpinski_ball_basis_tetrahedron(V1, V2, V3, V4, sample_step)

[M123,T] = sample_triangle(V1',V2',V3',sample_step); % OIJ

% Translation vector
tv = repmat(V1,[size(M123,1),1]);

% Rotations
M134 = ([0 0 -1; 0 1 0; 1 0 0]*(M123-tv)' + tv')';
M142 = ([1 0 0; 0 0 -1; 0 1 0]*(M123-tv)' + tv')';

M234 = sample_triangle(V2',V3',V4',sample_step); % IJK

C = 0.5 * [V1+V2; V1+V3; V1+V4; V2+V3; V2+V4; V3+V4];
 
M = [M123; M134; M142; M234];

% Triplet indices list 
T = [T;
     T +   repmat(size(M123,1),[size(T,1) size(T,2)]);...
     T + 2*repmat(size(M123,1),[size(T,1) size(T,2)]);...
     T + 3*repmat(size(M123,1),[size(T,1) size(T,2)])];

end


%% Sierpinski_tetrahedron_iterate subfunction
function [V, T] = Sierpinski_tetrahedron_iterate(V1, V2, V3, V4, nb_it, sample_step)

Summit_array = [V1; V2; V3; V4];
[Vertex_array, Triangle_array, Middle_edge_array] = create_Sierpinski_ball_basis_tetrahedron(V1,V2,V3,V4,sample_step);

p = 0;

while p ~= nb_it
        
    New_vertex_array      = repmat(Vertex_array,[1 1 4]);
    New_summit_array      = repmat(Summit_array,[1 1 4]);
    New_triangle_array    = repmat(Triangle_array,[1 1 4]);
    New_middle_edge_array = repmat(Middle_edge_array,[1 1 4]);
    
    for j = 1:size(Vertex_array,3) % Loop on current nb tetra               
        
        for i = 1:4
            
            V = Summit_array(i,:,j); % current summit
            
            D = sqrt(sum((Middle_edge_array(:,:,j) - repmat(V, [6 1])).^2,2)); % distance matrix
            [~,idx] = sort(D,1);
            
            New_summit_array(:,:,4*(j-1)+i) = [V; Middle_edge_array(idx(1),:,j); Middle_edge_array(idx(2),:,j); Middle_edge_array(idx(3),:,j)];                        
            
            i_zmax = find(New_summit_array(:,3,4*(j-1)+i) == max(New_summit_array(:,3,4*(j-1)+i)),1);
            i_xmax = find(New_summit_array(:,1,4*(j-1)+i) == max(New_summit_array(:,1,4*(j-1)+i)),1);
            i_ymax = find(New_summit_array(:,2,4*(j-1)+i) == max(New_summit_array(:,2,4*(j-1)+i)),1);
            i_rmin = find(sum(New_summit_array(:,:,4*(j-1)+i).^2, 2) == min(sum(New_summit_array(:,:,4*(j-1)+i).^2, 2)),1); % closest to origin                        
            
            
            New_summit_array(:,:,4*(j-1)+i) = [New_summit_array(i_rmin,:,4*(j-1)+i);...
                                               New_summit_array(i_xmax,:,4*(j-1)+i);...
                                               New_summit_array(i_ymax,:,4*(j-1)+i);...
                                               New_summit_array(i_zmax,:,4*(j-1)+i);...
                                               ];
                     
        
        % Create new Reuleaux : vertices, triangles, middle edge
        [New_vertex_array(:,:,4*(j-1)+i),New_triangle_array(:,:,4*(j-1)+i),New_middle_edge_array(:,:,4*(j-1)+i)] = ...
            create_Sierpinski_ball_basis_tetrahedron(New_summit_array(1,:,4*(j-1)+i),...
                                                     New_summit_array(2,:,4*(j-1)+i),...
                                                     New_summit_array(3,:,4*(j-1)+i),...
                                                     New_summit_array(4,:,4*(j-1)+i), sample_step);                                                                          
        end   
        
    end
    
    Vertex_array      = New_vertex_array;
    Summit_array      = New_summit_array;
    Triangle_array    = New_triangle_array;
    Middle_edge_array = New_middle_edge_array;
    
    p = p+1;
    
end

V = Vertex_array(:,:,1);
T = Triangle_array(:,:,1);

for k = 2: size(Vertex_array,3)
    
    T = cat(1, T, Triangle_array(:,:,k)+size(V,1));
    V = cat(1, V, Vertex_array(:,:,k));    
    
end

end


%% sample_triangle subfunction
function [V, T] = sample_triangle(V1, V2, V3, nbstep)

% Create sampling grid
global Ndim;

Ndim = size(V1,1);

% (V1V2, V1V3) base
u = (V2 - V1);
v = (V3 - V1);

V = zeros(sum(1:nbstep+1),Ndim);

nu = u / norm(u);
nv = v / norm(v);
stepu = norm(u) / nbstep;
stepv = norm(v) / nbstep;
k = 1;

% Sampling & vertices generation
for m = 0:nbstep
    
    for n = 0:nbstep
        
        if m+n <= nbstep % in (V1,V2,V3) triangle conditions ; indices # nb segments
            
            % translation vector
            tv = m*stepu*nu + n*stepv*nv;
            V(k,:) = (V1 + tv)';
            k = k+1;
            
        end
        
    end
    
end

% Index triplets list construction
T = zeros(nbstep^2,3);
row_length = 1 + nbstep;
cum_row_length = row_length;
row_idx = 1;
p = 1;

while p <= nbstep^2 && row_length > 1
    
     i = p;
    
    if p < 2 % "right" triangle serie only
        
        while (i < cum_row_length)
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    else
        
        % Since p >= 2
        while i < cum_row_length % both triangle series
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;            
            T(row_idx,:) = [i i-row_length i+1]; % + upside-down triangles serie
            row_idx = row_idx + 1;            
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    end
    
end

T = sort(T,2);
T = unique(T,'rows','stable');

end % sample_triangle


%% remove_duplicated_vertices subfunction
function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in)

tol = 1e4*eps;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);

end % remove_duplicated_vertices