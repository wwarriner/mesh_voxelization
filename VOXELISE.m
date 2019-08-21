function [ voxels, faces ] = VOXELISE(grid,fv,rays)
% VOXELISE  Voxelise a 3D triangular-polygon mesh.
%==========================================================================
% AUTHOR        Adam H. Aitkenhead
% CONTACT       adam.aitkenhead@christie.nhs.uk
% INSTITUTION   The Christie NHS Foundation Trust
%
% INPUTS
%
%     grid - Mandatory - Grid object
%
%     fv  - Optional  - structure     - Structure containing the faces and vertices
%                                           of the mesh, in the same format as that produced
%                                           by the isosurface command.
%
%     rays - Mandatory - String    - Defines the directions in which ray-tracing
%                                           is performed.  The default is 'xyz', which
%                                           traces in the x,y,z directions and combines
%                                           the results.
%
% OUTPUTS
%
%     gridOUTPUT - Mandatory - PxQxR logical array - Voxelised data (1=>Inside the mesh, 0=>Outside the mesh)
%
% NOTES
%
%   - The mesh must be properly closed (ie. watertight).
%   - Defining rays='xyz' means that the mesh is ray-traced in each
%     of the x,y,z directions, with the overall result being a combination
%     of the result from each direction.  This gives the most reliable
%     result at the expense of computation time.
%   - Tracing in only one direction (eg. rays='z') is faster, but
%     can potentially lead to artefacts where a ray exactly crosses
%     several facet edges.
%
% REFERENCES
%
%   - This code uses a ray intersection method similar to that described by:
%     Patil S and Ravi B.  Voxel-based representation, display and
%     thickness analysis of intricate shapes. Ninth International
%     Conference on Computer Aided Design and Computer Graphics (CAD/CG
%     2005)
%==========================================================================

if nargin < 3
    rays = 'xyz';
end

assert( isa( grid, 'Grid' ) );
assert( grid.dimension_count == 3 );

assert( isstruct( fv ) );
assert( isfield( fv, 'faces' ) );
assert( isfield( fv, 'vertices' ) );

assert( ischar( rays ) );
assert( ismember( numel( rays ), 1 : 3 ) );
assert( all( ismember( rays, 'xyz' ) ) );

mesh = CONVERT_meshformat( fv.faces, fv.vertices );

voxels = zeros( grid.shape );

if contains( rays, 'x' )
    p = [ 2 3 1 ];
    p_inv = [ 3 1 2 ];
    voxels = voxels + permute( VOXELISEinternal( ...
        grid.shape, ...
        grid.data{ p( 1 ) }, ...
        grid.data{ p( 2 ) }, ...
        grid.data{ p( 3 ) }, ...
        mesh( :, p, : ) ...
        ), p_inv );
end

if contains( rays, 'y' )
    p = [ 3 1 2 ];
    p_inv = [ 2 3 1 ];
    voxels = voxels + permute( VOXELISEinternal( ...
        grid.shape, ...
        grid.data{ p( 1 ) }, ...
        grid.data{ p( 2 ) }, ...
        grid.data{ p( 3 ) }, ...
        mesh( :, p, : ) ...
        ), p_inv );
end

if contains( rays, 'z' )
    p = [ 1 2 3 ];
    p_inv = [ 1 2 3 ];
    voxels = voxels + permute( VOXELISEinternal( ...
        grid.shape, ...
        grid.data{ p( 1 ) }, ...
        grid.data{ p( 2 ) }, ...
        grid.data{ p( 3 ) }, ...
        mesh( :, p, : ) ...
        ), p_inv );
end

voxels = voxels >= numel( rays ) ./ 2;

end


function [ v, f ] = VOXELISEinternal(shape,X,Y,Z,mesh)

% prep
v = false( shape );
f = uint32( false( shape ) );

x_index_range = get_index_range( X, mesh( :, 1, : ) );
y_index_range = get_index_range( Y, mesh( :, 2, : ) );
z_range = get_z_range( mesh( :, 3, : ) );

mesh_min = min( mesh, [], 3 );
mesh_max = max( mesh, [], 3 );

% cull faces that are definitely not crossed using square criteria
[ crossed, xp, yp, xi, yi, ray_count ] = ...
    determine_possibly_crossed_faces( ...
    X, ...
    x_index_range, ...
    Y, ...
    y_index_range, ...
    mesh_min, ...
    mesh_max ...
    );

% build data for faces that may be crossed
cross_counts = arrayfun( @(x)numel(x{1}), crossed );
faces = cell2mat( crossed );
cross_count = numel( faces );
rays = zeros( cross_count, 1 );
xpe = zeros( cross_count, 1 );
ype = zeros( cross_count, 1 );
start = 1;
for i = 1 : ray_count
    finish = cross_counts( i ) + start - 1;
    rays( start : finish ) = i;
    xpe( start : finish ) = xp( i );
    ype( start : finish ) = yp( i );
    start = finish + 1;
end

% determine faces that are crossed at a vertex
[ crossed, corrections ] = determine_possibly_crossed_vertices( ...
    mesh, ...
    faces, ...
    rays, ...
    xpe, ...
    ype, ...
    ray_count ...
    );

% cull rays that require corrections
check = ~ismember( rays, find( corrections ) );
faces = faces( check );
rays = rays( check );
xpe = xpe( check );
ype = ype( check );

% determine faces that are actually crossed by rays
f_crossed = determine_actual_crossed_faces( ...
    mesh, ...
    faces, ...
    xpe, ...
    ype, ...
    rays, ...
    ray_count ...
    );

% combine vertex and face crossings
for i = 1 : ray_count
    crossed{ i } = [ crossed{ i }; f_crossed{ i } ];
end

% rebuild data based on actual crossings
cross_counts = arrayfun( @(x)numel(x{1}), crossed );
faces = cell2mat( crossed );
cross_count = numel( faces );
rays = zeros( cross_count, 1 );
xpe = zeros( cross_count, 1 );
xie = zeros( cross_count, 1 );
ype = zeros( cross_count, 1 );
yie = zeros( cross_count, 1 );
start = 1;
for i = 1 : ray_count
    finish = cross_counts( i ) + start - 1;
    rays( start : finish ) = i;
    xpe( start : finish ) = xp( i );
    xie( start : finish ) = xi( i );
    ype( start : finish ) = yp( i );
    yie( start : finish ) = yi( i );
    start = finish + 1;
end

% rays_with_faces = unique( ray );
% rays_without = setdiff( ( 1 : ray_count ).', rays_with_faces );
% corrections( rays_without ) = true;

% compute z-values of face crossings
A = mesh(faces,2,1) .* ( mesh(faces,3,2) - mesh(faces,3,3) ) ...
    + mesh(faces,2,2) .* ( mesh(faces,3,3) - mesh(faces,3,1) ) ...
    + mesh(faces,2,3) .* ( mesh(faces,3,1) - mesh(faces,3,2) );
B = mesh(faces,3,1) .* ( mesh(faces,1,2) - mesh(faces,1,3) ) ...
    + mesh(faces,3,2) .* ( mesh(faces,1,3) - mesh(faces,1,1) ) ...
    + mesh(faces,3,3) .* ( mesh(faces,1,1) - mesh(faces,1,2) );
C = mesh(faces,1,1) .* ( mesh(faces,2,2) - mesh(faces,2,3) ) ...
    + mesh(faces,1,2) .* ( mesh(faces,2,3) - mesh(faces,2,1) ) ...
    + mesh(faces,1,3) .* ( mesh(faces,2,1) - mesh(faces,2,2) );
D = - mesh(faces,1,1) .* ( mesh(faces,2,2).*mesh(faces,3,3) - mesh(faces,2,3).*mesh(faces,3,2) ) ...
    - mesh(faces,1,2) .* ( mesh(faces,2,3).*mesh(faces,3,1) - mesh(faces,2,1).*mesh(faces,3,3) ) ...
    - mesh(faces,1,3) .* ( mesh(faces,2,1).*mesh(faces,3,2) - mesh(faces,2,2).*mesh(faces,3,1) );

TOL = 1e-12;
C( abs( C ) < TOL ) = 0;
z = -( D + A.*xpe + B.*ype ) ./ C;
z_faces = faces;
check = z_range( 1 ) <= z & z <= z_range( 2 );
z = z( check );
z_faces = z_faces( check );
rays = rays( check );
xie = xie( check );
yie = yie( check );

ROUND_TOL = 1e12;
z = round( z * ROUND_TOL ) / ROUND_TOL;
tt = [ rays xie yie z ];
tt = sortrows( unique( tt, 'stable', 'rows' ), 1 : size( tt, 2 ) );
rays = tt( :, 1 );
xie = tt( :, 2 );
yie = tt( :, 3 );
z = tt( :, 4 );
fill_count = numel( z ./ 2 );
ind = 1;
while ind < fill_count
    if rays( ind ) == rays( ind + 1 )
        inside = z( ind ) < Z & Z < z( ind + 1 );
        v( xie( ind ), yie( ind ), inside ) = true;
        ind = ind + 2;
    else
        corrections( rays( ind ) ) = true;
        ind = ind + 1;
    end
end

%======================================================
% USE INTERPOLATION TO FILL IN THE RAYS WHICH COULD NOT BE VOXELISED
%======================================================
%For rays where the voxelisation did not give a clear result, the ray is
%computed by interpolating from the surrounding rays.
countCORRECTIONLIST = sum( corrections );

if countCORRECTIONLIST>0
    xc = xi( corrections );
    yc = yi( corrections );
    %If necessary, add a one-pixel border around the x and y edges of the
    %array.  This prevents an error if the code tries to interpolate a ray at
    %the edge of the x,y grid.
    if min(xc)==x_index_range(1) ...
            || max(xc)==x_index_range(2) ...
            || min(yc)==y_index_range(1) ...
            || max(yc)==y_index_range(2)
        v = [zeros(1,shape(2)+2,shape(3));zeros(shape(1),1,shape(3)),v,zeros(shape(1),1,shape(3));zeros(1,shape(2)+2,shape(3))];
        f = [zeros(1,shape(2)+2,shape(3));zeros(shape(1),1,shape(3)),f,zeros(shape(1),1,shape(3));zeros(1,shape(2)+2,shape(3))];
        xi = xi + 1;
        yi = yi + 1;
    end
    
    for loopC = 1:countCORRECTIONLIST
        c = corrections(loopC);
        xc = xi( c );
        yc = yi( c );
        voxelsforcorrection = squeeze( sum( [ ...
            v(xc-1,yc-1,:),...
            v(xc-1,yc,:),...
            v(xc-1,yc+1,:),...
            v(xc,yc-1,:),...
            v(xc,yc+1,:),...
            v(xc+1,yc-1,:),...
            v(xc+1,yc,:),...
            v(xc+1,yc+1,:),...
            ] ) );
        voxelsforcorrection = (voxelsforcorrection>=4);
        v(xc,yc,voxelsforcorrection) = true;
%         
%         f = [ ...
%             faces(xc-1,yc-1,:),...
%             faces(xc-1,yc,:),...
%             faces(xc-1,yc+1,:),...
%             faces(xc,yc-1,:),...
%             faces(xc,yc+1,:),...
%             faces(xc+1,yc-1,:),...
%             faces(xc+1,yc,:),...
%             faces(xc+1,yc+1,:),...
%             ];
%         f( f==0 ) = nan;
%         f = squeeze( mode( f ) );
%         f( isnan( f ) ) = 0;
%         faces(corrections(loopC,1),corrections(loopC,2),voxelsforcorrection) = f(voxelsforcorrection);
    end
    
    %Remove the one-pixel border surrounding the array, if this was added
    %previously.
    if size(v,1)>numel(X) || size(v,2)>numel(Y)
        v = v(2:end-1,2:end-1,:);
        f = f(2:end-1,2:end-1,:);
    end
    
end %if

%disp([' Ray tracing result: ',num2str(countCORRECTIONLIST),' rays (',num2str(countCORRECTIONLIST/(shape(1)*voxcountY)*100,'%5.1f'),'% of all rays) exactly crossed a facet edge and had to be computed by interpolation.'])

end %function


function [ crossed, xp, yp, xi, yi, ray_count ] = ...
    determine_possibly_crossed_faces( ...
    X, ...
    x_index_range, ...
    Y, ...
    y_index_range, ...
    mesh_min, ...
    mesh_max ...
    )

x_count = diff( x_index_range ) + 1;
y_count = diff( y_index_range ) + 1;
ray_count = x_count * y_count;
crossed = cell( ray_count, 1 );
xp = zeros( ray_count, 1 );
xi = zeros( ray_count, 1 );
yp = zeros( ray_count, 1 );
yi = zeros( ray_count, 1 );
i = 1;
for y_index = y_index_range( 1 ) : y_index_range( 2 )
    y = Y( y_index );
    cross_y = find( y <= mesh_max( :, 2 ) & mesh_min( :, 2 ) <= y );
    for x_index = x_index_range( 1 ) : x_index_range( 2 )
        x = X( x_index );
        check = x <= mesh_max( cross_y, 1 ) & mesh_min( cross_y, 1 ) <= x;
        if any( check )
            crossed{ i } = cross_y( check );
            xp( i ) = x;
            xi( i ) = x_index;
            yp( i ) = y;
            yi( i ) = y_index;
            i = i + 1;
        end
    end
end
crossed( i : end ) = [];
xp( i : end ) = [];
xi( i : end ) = [];
yp( i : end ) = [];
yi( i : end ) = [];
ray_count = i - 1;

end


function [ crossed, corrections ] = determine_possibly_crossed_vertices( ...
    mesh, ...
    faces, ...
    rays, ...
    xpe, ...
    ype, ...
    ray_count ...
    )

corrections = false( ray_count, 1 );
crossed = cell( ray_count, 1 );
check = ( mesh(faces,1,1)==xpe & mesh(faces,2,1)==ype ) ...
    | ( mesh(faces,1,2)==xpe & mesh(faces,2,2)==ype ) ...
    | ( mesh(faces,1,3)==xpe & mesh(faces,2,3)==ype );
check = find( check );
for i = 1 : numel( check )
    vi = check( i );
    [ needs_correction, f ] = find_vertex_crossings( mesh, faces( vi ) );
    ri = rays( vi );
    if needs_correction
        corrections( ri ) = true;
    end
    crossed{ ri } = [ crossed{ ri }; f ];
end
crossed( corrections ) = {[]};

end


function crossed = determine_actual_crossed_faces( ...
    mesh, ...
    faces, ...
    xpe, ...
    ype, ...
    rays, ...
    ray_count ...
    )

inds = ( 1 : numel( faces ) ).';
rem = faces;
xx = xpe;
yy = ype;
fy = mesh(rem,2,2) - ((mesh(rem,2,2)-mesh(rem,2,3)) .* (mesh(rem,1,2)-mesh(rem,1,1))./(mesh(rem,1,2)-mesh(rem,1,3)));
ry = mesh(rem,2,2) - ((mesh(rem,2,2)-mesh(rem,2,3)) .* (mesh(rem,1,2)-xx)./(mesh(rem,1,2)-mesh(rem,1,3)));
check = (fy >= mesh(rem,2,1) & ry >= yy) | (fy <= mesh(rem,2,1) & ry <= yy);

inds = inds( check );
rem = rem( check );
xx = xx( check );
yy = yy( check );
fy = mesh(rem,2,3) - ((mesh(rem,2,3)-mesh(rem,2,1)) .* (mesh(rem,1,3)-mesh(rem,1,2))./(mesh(rem,1,3)-mesh(rem,1,1)));
ry = mesh(rem,2,3) - ((mesh(rem,2,3)-mesh(rem,2,1)) .* (mesh(rem,1,3)-xx)./(mesh(rem,1,3)-mesh(rem,1,1)));
check = (fy >= mesh(rem,2,2) & ry >= yy) | (fy <= mesh(rem,2,2) & ry <= yy);

inds = inds( check );
rem = rem( check );
xx = xx( check );
yy = yy( check );
fy = mesh(rem,2,1) - ((mesh(rem,2,1)-mesh(rem,2,2)) .* (mesh(rem,1,1)-mesh(rem,1,3))./(mesh(rem,1,1)-mesh(rem,1,2)));
ry = mesh(rem,2,1) - ((mesh(rem,2,1)-mesh(rem,2,2)) .* (mesh(rem,1,1)-xx)./(mesh(rem,1,1)-mesh(rem,1,2)));
check = (fy >= mesh(rem,2,3) & ry >= yy) | (fy <= mesh(rem,2,3) & ry <= yy);

inds = inds( check );
rem = rem( check );
crossed = cell( ray_count, 1 );
for i = 1 : numel( inds )
    index = inds( i );
    ri = rays( index );
    crossed{ ri } = [ crossed{ ri }; rem( i ) ];
end

end


function value = get_index_range( x, mesh_x )

r = [ min( mesh_x, [], 'all' ) max( mesh_x, [], 'all' ) ];
r_min = x - r( 1 );
p_min = find( r_min == min( r_min ), 1, 'first' );
r_max = x - r( 2 );
p_max = find( r_max == max( r_max ), 1, 'last' );
value = sort( [ p_min p_max ] );

end


function value = get_z_range( mesh_z )

TOL = 1e-12;
value = [ ...
    min( mesh_z, [], 'all' ) - TOL ...
    max( mesh_z, [], 'all' ) + TOL...
    ];

end


function [ needs_correction, face_cross ] = find_vertex_crossings( mesh, cross )

needs_correction = false;
face_cross = [];
if isempty( cross )
    return;
end

check = zeros( 1, numel( cross ) );
while min( check ) == 0
    index = find( check == 0, 1, 'first' );
    check( index ) = 1;
    
    [ fv.faces, fv.vertices ] = CONVERT_meshformat( mesh( cross, :, : ) );
    adjacent = ismember( fv.faces, fv.faces( index, : ) );
    adjacent = max( adjacent, [], 2 );
    check( adjacent ) = 1;
    
    normals = COMPUTE_mesh_normals( mesh( cross( adjacent ), :, : ) );
    if max( normals( :, 3 ) ) < 0 || min( normals( :, 3 ) ) > 0
        face_cross = [ face_cross cross( index ) ]; %#ok<AGROW>
    else
        needs_correction = true;
        break;
    end
end

end

