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

meshXYZ = CONVERT_meshformat( fv.faces, fv.vertices );

ray_count = numel( rays );

voxels = false( [ grid.shape ray_count ] );
faces = uint32( false( [ grid.shape ray_count ] ) );
direction_count = 0;

if contains( rays, 'x' )
    direction_count = direction_count + 1;
    p = [ 2 3 1 ];
    p_inv = [ 3 1 2 ];
    [ g, f ] = VOXELISEinternal( ...
        grid.shape, ...
        grid.data{ p( 1 ) }, ...
        grid.data{ p( 2 ) }, ...
        grid.data{ p( 3 ) }, ...
        meshXYZ( :, p, : ) ...
        );
    voxels( :, :, :, direction_count ) = permute( g, p_inv );
    faces( :, :, :, direction_count ) = permute( f, p_inv );
    clear( 'g', 'f' );
end

if contains( rays, 'y' )
    direction_count = direction_count + 1;
    p = [ 3 1 2 ];
    p_inv = [ 2 3 1 ];
    [ g, f ] = VOXELISEinternal( ...
        grid.shape, ...
        grid.data{ p( 1 ) }, ...
        grid.data{ p( 2 ) }, ...
        grid.data{ p( 3 ) }, ...
        meshXYZ( :, p, : ) ...
        );
    voxels( :, :, :, direction_count ) = permute( g, p_inv );
    faces( :, :, :, direction_count ) = permute( f, p_inv );
    clear( 'g', 'f' );
end

if contains( rays, 'z' )
    direction_count = direction_count + 1;
    p = [ 1 2 3 ];
    p_inv = [ 1 2 3 ];
    [ g, f ] = VOXELISEinternal( ...
        grid.shape, ...
        grid.data{ p( 1 ) }, ...
        grid.data{ p( 2 ) }, ...
        grid.data{ p( 3 ) }, ...
        meshXYZ( :, p, : ) ...
        );
    voxels( :, :, :, direction_count ) = permute( g, p_inv );
    faces( :, :, :, direction_count ) = permute( f, p_inv );
    clear( 'g', 'f' );
end

if numel(rays)>1
    voxels = sum(voxels,4)>=numel(rays)/2;
    f = sum(faces,4);
    f( sum(faces>0,4)>1 ) = 0;
    faces = f;
end

end


function [voxels,faces] = VOXELISEinternal(shape,X,Y,Z,meshXYZ)

voxels = false( shape );
faces = uint32( false( shape ) );

x_index_range = get_index_range( X, meshXYZ( :, 1, : ) );
y_index_range = get_index_range( Y, meshXYZ( :, 2, : ) );
z_range = get_z_range( meshXYZ( :, 3, : ) );

mesh_min = min( meshXYZ, [], 3 );
mesh_max = max( meshXYZ, [], 3 );

corrections = [];
for y_index = y_index_range( 1 ) : y_index_range( 2 )
    cross_y = find( mesh_min(:,2)<=Y(y_index) & mesh_max(:,2)>=Y(y_index) );
    for x_index = x_index_range( 1 ) : x_index_range( 2 )
        cross = cross_y( mesh_min(cross_y,1)<=X(x_index) & mesh_max(cross_y,1)>=X(x_index) );
        
        if isempty(cross)==0  %Only continue the analysis if some nearby facets were actually identified
            
            % - 2 - For each facet, check if the ray really does cross the facet rather than just passing it close-by:
            
            % GENERAL METHOD:
            % A. Take each edge of the facet in turn.
            % B. Find the position of the opposing vertex to that edge.
            % C. Find the position of the ray relative to that edge.
            % D. Check if ray is on the same side of the edge as the opposing vertex.
            % E. If this is true for all three edges, then the ray definitely passes through the facet.
            %
            % NOTES:
            % A. If a ray crosses exactly on a vertex:
            %    a. If the surrounding facets have normal components pointing in the same (or opposite) direction as the ray then the face IS crossed.
            %    b. Otherwise, add the ray to the correctionlist.
            
            facetCROSSLIST = [];   %Prepare to record all facets which are crossed by the ray.  This array is built on-the-fly, but since
            %it ought to be relatively small (typically a list of <10) should not incur too much of a speed penalty.
            
            %----------
            % - 1 - Check for crossed vertices:
            %----------
            
            % Find which mesh facets contain a vertex which is crossed by the ray:
            vertexCROSSLIST = cross( (meshXYZ(cross,1,1)==X(x_index) & meshXYZ(cross,2,1)==Y(y_index)) ...
                | (meshXYZ(cross,1,2)==X(x_index) & meshXYZ(cross,2,2)==Y(y_index)) ...
                | (meshXYZ(cross,1,3)==X(x_index) & meshXYZ(cross,2,3)==Y(y_index)) ...
                );
            
            if isempty(vertexCROSSLIST)==0  %Only continue the analysis if potential vertices were actually identified
                
                checkindex = zeros(1,numel(vertexCROSSLIST));
                
                while min(checkindex) == 0
                    
                    vertexindex             = find(checkindex==0,1,'first');
                    checkindex(vertexindex) = 1;
                    
                    [temp.faces,temp.vertices] = CONVERT_meshformat(meshXYZ(vertexCROSSLIST,:,:));
                    adjacentindex              = ismember(temp.faces,temp.faces(vertexindex,:));
                    adjacentindex              = max(adjacentindex,[],2);
                    checkindex(adjacentindex)  = 1;
                    
                    coN = COMPUTE_mesh_normals(meshXYZ(vertexCROSSLIST(adjacentindex),:,:));
                    
                    if max(coN(:,3))<0 || min(coN(:,3))>0
                        facetCROSSLIST    = [facetCROSSLIST,vertexCROSSLIST(vertexindex)];
                    else
                        cross = [];
                        corrections    = [ corrections; x_index,y_index ];
                        checkindex(:)     = 1;
                    end
                    
                end
                
            end
            
            %----------
            % - 2 - Check for crossed facets:
            %----------
            
            if isempty(cross)==0  %Only continue the analysis if some nearby facets were actually identified
                
                for loopCHECKFACET = cross'
                    
                    %Check if ray crosses the facet.  This method is much (>>10 times) faster than using the built-in function 'inpolygon'.
                    %Taking each edge of the facet in turn, check if the ray is on the same side as the opposing vertex.
                    
                    Y1predicted = meshXYZ(loopCHECKFACET,2,2) - ((meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,1))/(meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,3)));
                    YRpredicted = meshXYZ(loopCHECKFACET,2,2) - ((meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-X(x_index))/(meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,3)));
                    
                    if (Y1predicted >= meshXYZ(loopCHECKFACET,2,1) && YRpredicted >= Y(y_index)) || (Y1predicted <= meshXYZ(loopCHECKFACET,2,1) && YRpredicted <= Y(y_index))
                        %The ray is on the same side of the 2-3 edge as the 1st vertex.
                        
                        Y2predicted = meshXYZ(loopCHECKFACET,2,3) - ((meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,2))/(meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,1)));
                        YRpredicted = meshXYZ(loopCHECKFACET,2,3) - ((meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-X(x_index))/(meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,1)));
                        
                        if (Y2predicted >= meshXYZ(loopCHECKFACET,2,2) && YRpredicted >= Y(y_index)) || (Y2predicted <= meshXYZ(loopCHECKFACET,2,2) && YRpredicted <= Y(y_index))
                            %The ray is on the same side of the 3-1 edge as the 2nd vertex.
                            
                            Y3predicted = meshXYZ(loopCHECKFACET,2,1) - ((meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,3))/(meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,2)));
                            YRpredicted = meshXYZ(loopCHECKFACET,2,1) - ((meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-X(x_index))/(meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,2)));
                            
                            if (Y3predicted >= meshXYZ(loopCHECKFACET,2,3) && YRpredicted >= Y(y_index)) || (Y3predicted <= meshXYZ(loopCHECKFACET,2,3) && YRpredicted <= Y(y_index))
                                %The ray is on the same side of the 1-2 edge as the 3rd vertex.
                                
                                %The ray passes through the facet since it is on the correct side of all 3 edges
                                facetCROSSLIST = [facetCROSSLIST,loopCHECKFACET];
                                
                            end %if
                        end %if
                    end %if
                    
                end %for
                
                
                %----------
                % - 3 - Find the z coordinate of the locations where the ray crosses each facet or vertex:
                %----------
                
                gridCOzCROSS = zeros(size(facetCROSSLIST));
                gridCOzFACE = zeros(size(facetCROSSLIST));
                for loopFINDZ = facetCROSSLIST
                    
                    % METHOD:
                    % 1. Define the equation describing the plane of the facet.  For a
                    % more detailed outline of the maths, see:
                    % http://local.wasp.uwa.edu.au/~pbourke/geometry/planeeq/
                    %    Ax + By + Cz + D = 0
                    %    where  A = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2)
                    %           B = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2)
                    %           C = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2)
                    %           D = - x1 (y2 z3 - y3 z2) - x2 (y3 z1 - y1 z3) - x3 (y1 z2 - y2 z1)
                    % 2. For the x and y coordinates of the ray, solve these equations to find the z coordinate in this plane.
                    
                    planecoA = meshXYZ(loopFINDZ,2,1)*(meshXYZ(loopFINDZ,3,2)-meshXYZ(loopFINDZ,3,3)) + meshXYZ(loopFINDZ,2,2)*(meshXYZ(loopFINDZ,3,3)-meshXYZ(loopFINDZ,3,1)) + meshXYZ(loopFINDZ,2,3)*(meshXYZ(loopFINDZ,3,1)-meshXYZ(loopFINDZ,3,2));
                    planecoB = meshXYZ(loopFINDZ,3,1)*(meshXYZ(loopFINDZ,1,2)-meshXYZ(loopFINDZ,1,3)) + meshXYZ(loopFINDZ,3,2)*(meshXYZ(loopFINDZ,1,3)-meshXYZ(loopFINDZ,1,1)) + meshXYZ(loopFINDZ,3,3)*(meshXYZ(loopFINDZ,1,1)-meshXYZ(loopFINDZ,1,2));
                    planecoC = meshXYZ(loopFINDZ,1,1)*(meshXYZ(loopFINDZ,2,2)-meshXYZ(loopFINDZ,2,3)) + meshXYZ(loopFINDZ,1,2)*(meshXYZ(loopFINDZ,2,3)-meshXYZ(loopFINDZ,2,1)) + meshXYZ(loopFINDZ,1,3)*(meshXYZ(loopFINDZ,2,1)-meshXYZ(loopFINDZ,2,2));
                    planecoD = - meshXYZ(loopFINDZ,1,1)*(meshXYZ(loopFINDZ,2,2)*meshXYZ(loopFINDZ,3,3)-meshXYZ(loopFINDZ,2,3)*meshXYZ(loopFINDZ,3,2)) - meshXYZ(loopFINDZ,1,2)*(meshXYZ(loopFINDZ,2,3)*meshXYZ(loopFINDZ,3,1)-meshXYZ(loopFINDZ,2,1)*meshXYZ(loopFINDZ,3,3)) - meshXYZ(loopFINDZ,1,3)*(meshXYZ(loopFINDZ,2,1)*meshXYZ(loopFINDZ,3,2)-meshXYZ(loopFINDZ,2,2)*meshXYZ(loopFINDZ,3,1));
                    
                    if abs(planecoC) < 1e-14
                        planecoC=0;
                    end
                    
                    gridCOzCROSS(facetCROSSLIST==loopFINDZ) = (- planecoD - planecoA*X(x_index) - planecoB*Y(y_index)) / planecoC;
                    gridCOzFACE(facetCROSSLIST==loopFINDZ) = loopFINDZ;
                    
                end %for
                
                %Remove values of gridCOzCROSS which are outside of the mesh limits (including a 1e-12 margin for error).
                gridCOzCROSS = gridCOzCROSS( z_range( 1 ) <= gridCOzCROSS & gridCOzCROSS <= z_range( 2 ) );
                
                %Round gridCOzCROSS to remove any rounding errors, and take only the unique values:
                gridCOzCROSS = round(gridCOzCROSS*1e12)/1e12;
                [ gridCOzCROSS, ia ] = unique(gridCOzCROSS);
                gridCOzFACE = gridCOzFACE( ia );
                
                %----------
                % - 4 - Label as being inside the mesh all the voxels that the ray passes through after crossing one facet before crossing another facet:
                %----------
                
                if rem(numel(gridCOzCROSS),2)==0  % Only rays which cross an even number of facets are voxelised
                    
                    for loopASSIGN = 1:(numel(gridCOzCROSS)/2)
                        voxelsINSIDE = (Z>gridCOzCROSS(2*loopASSIGN-1) & Z<gridCOzCROSS(2*loopASSIGN));
                        voxels(x_index,y_index,voxelsINSIDE) = 1;
                        faces(x_index,y_index,find(voxelsINSIDE,1)) = gridCOzFACE(2*loopASSIGN-1);
                        faces(x_index,y_index,find(voxelsINSIDE,1,'last')) = gridCOzFACE(2*loopASSIGN);
                    end %for
                    
                elseif numel(gridCOzCROSS)~=0    % Remaining rays which meet the mesh in some way are not voxelised, but are labelled for correction later.
                    
                    corrections = [ corrections; x_index,y_index ];
                    
                end %if
                
            end
            
        end %if
        
    end %for
end %for


%======================================================
% USE INTERPOLATION TO FILL IN THE RAYS WHICH COULD NOT BE VOXELISED
%======================================================
%For rays where the voxelisation did not give a clear result, the ray is
%computed by interpolating from the surrounding rays.
countCORRECTIONLIST = size(corrections,1);

if countCORRECTIONLIST>0
    
    %If necessary, add a one-pixel border around the x and y edges of the
    %array.  This prevents an error if the code tries to interpolate a ray at
    %the edge of the x,y grid.
    if min(corrections(:,1))==1 || max(corrections(:,1))==numel(X) || min(corrections(:,2))==1 || max(corrections(:,2))==numel(Y)
        voxels     = [zeros(1,voxcountY+2,voxcountZ);zeros(voxcountX,1,voxcountZ),voxels,zeros(voxcountX,1,voxcountZ);zeros(1,voxcountY+2,voxcountZ)];
        faces     = [zeros(1,voxcountY+2,voxcountZ);zeros(voxcountX,1,voxcountZ),faces,zeros(voxcountX,1,voxcountZ);zeros(1,voxcountY+2,voxcountZ)];
        corrections = corrections + 1;
    end
    
    for loopC = 1:countCORRECTIONLIST
        voxelsforcorrection = squeeze( sum( [ voxels(corrections(loopC,1)-1,corrections(loopC,2)-1,:) ,...
            voxels(corrections(loopC,1)-1,corrections(loopC,2),:)   ,...
            voxels(corrections(loopC,1)-1,corrections(loopC,2)+1,:) ,...
            voxels(corrections(loopC,1),corrections(loopC,2)-1,:)   ,...
            voxels(corrections(loopC,1),corrections(loopC,2)+1,:)   ,...
            voxels(corrections(loopC,1)+1,corrections(loopC,2)-1,:) ,...
            voxels(corrections(loopC,1)+1,corrections(loopC,2),:)   ,...
            voxels(corrections(loopC,1)+1,corrections(loopC,2)+1,:) ,...
            ] ) );
        voxelsforcorrection = (voxelsforcorrection>=4);
        voxels(corrections(loopC,1),corrections(loopC,2),voxelsforcorrection) = 1;
        
        f = [ faces(corrections(loopC,1)-1,corrections(loopC,2)-1,:) ,...
            faces(corrections(loopC,1)-1,corrections(loopC,2),:)   ,...
            faces(corrections(loopC,1)-1,corrections(loopC,2)+1,:) ,...
            faces(corrections(loopC,1),corrections(loopC,2)-1,:)   ,...
            faces(corrections(loopC,1),corrections(loopC,2)+1,:)   ,...
            faces(corrections(loopC,1)+1,corrections(loopC,2)-1,:) ,...
            faces(corrections(loopC,1)+1,corrections(loopC,2),:)   ,...
            faces(corrections(loopC,1)+1,corrections(loopC,2)+1,:) ,...
            ];
        f( f==0 ) = nan;
        f = squeeze( mode( f ) );
        f( isnan( f ) ) = 0;
        faces(corrections(loopC,1),corrections(loopC,2),voxelsforcorrection) = f(voxelsforcorrection);
    end %for
    
    %Remove the one-pixel border surrounding the array, if this was added
    %previously.
    if size(voxels,1)>numel(X) || size(voxels,2)>numel(Y)
        voxels = voxels(2:end-1,2:end-1,:);
        faces = faces(2:end-1,2:end-1,:);
    end
    
end %if

%disp([' Ray tracing result: ',num2str(countCORRECTIONLIST),' rays (',num2str(countCORRECTIONLIST/(voxcountX*voxcountY)*100,'%5.1f'),'% of all rays) exactly crossed a facet edge and had to be computed by interpolation.'])

end %function


function value = get_index_range( x, mesh_x )

r = [ min( mesh_x, [], 'all' ) max( mesh_x, [], 'all' ) ];
r_min = x - r( 1 );
p_min = find( r_min == min( r_min ) );
r_max = x - r( 2 );
p_max = find( r_max == max( r_max ) );
value = sort( [ p_min p_max ] );

end


function value = get_z_range( mesh_z )

TOL = 1e-12;
value = [ ...
    min( mesh_z, [], 'all' ) - TOL ...
    max( mesh_z, [], 'all' ) + TOL...
    ];

end

