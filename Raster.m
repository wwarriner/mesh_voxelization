classdef Raster < handle
    
    properties ( SetAccess = private )
        interior(:,:,:) logical
        faces table
    end
    
    properties ( SetAccess = private )
        normals table
    end
    
    methods
        function obj = Raster( grid, fv, rays )
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
            dimensions = find( ismember( 'xyz', rays ) );
            
            obj.grid = grid;
            obj.mesh = mesh;
            obj.dimensions = dimensions;
            
            obj.generate();
        end
        
        function value = get.normals( obj )
            n = COMPUTE_mesh_normals( obj.mesh );
            f = obj.faces;
            value = f( :, : );
            value.faces = [];
            value.dimension = [];
            value.x = n( f.faces, 1 );
            value.y = n( f.faces, 2 );
            value.z = n( f.faces, 3 );
            value = varfun( ...
                @sum, ...
                value, ...
                'groupingvariables', { 'indices' } ...
                );
            value.GroupCount = [];
            value.Properties.VariableNames{'sum_x'} = 'x';
            value.Properties.VariableNames{'sum_y'} = 'y';
            value.Properties.VariableNames{'sum_z'} = 'z';
            t = value{ :, { 'x' 'y' 'z' } };
            v = vecnorm( t, 2, 2 );
            t = t ./ v;
            value{ :, { 'x' 'y' 'z' } } = t;
        end
    end
    
    properties ( Access = private )
        grid Grid
        mesh(:,3,3) double {mustBeReal,mustBeFinite}
        dimensions(1,:) double
    end
    
    methods ( Access = private )
        function generate( obj )
            v = false( obj.grid.shape );
            f = table();
            for i = obj.dimensions
                c = obj.cast( i );
                p_inv = obj.get_inverse_permutation( i );
                v = v + permute( c.interior_array, p_inv );
                face_list = c.get_face_list( p_inv );
                face_list.dimension(:) = i;
                f = [ f; face_list ]; %#ok<AGROW>
            end
            v = v >= numel( obj.dimensions ) ./ 2;
            
            obj.interior = v;
            obj.faces = f;
        end
        
        function rr = cast( obj, dimension )
            p = obj.get_permutation( dimension );
            rr = RayRaster( ...
                obj.grid.points{ p }, ...
                obj.mesh( :, p, : ) ...
                );
        end
    end
    
    properties ( Access = private, Constant )
        PERMUTATIONS = [ ...
            2 3 1; ...
            3 1 2; ...
            1 2 3 ...
            ];
        INVERSE_PERMUTATIONS = [ ...
            3 1 2; ...
            2 3 1; ...
            1 2 3 ...
            ];
    end
    
    methods ( Access = private, Static )
        function p = get_permutation( dimension )
            p = Raster.PERMUTATIONS( dimension, : );
        end
        
        function p = get_inverse_permutation( dimension )
            p = Raster.INVERSE_PERMUTATIONS( dimension, : );
        end
    end
    
end

