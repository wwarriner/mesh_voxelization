classdef Voxels < handle
    
    properties ( SetAccess = private )
        v(:,:,:) logical
    end
    
    methods
        function obj = Voxels( grid, fv, rays )
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
        end
        
        function v = generate( obj )
            v = false( obj.grid.shape );
            for i = obj.dimensions
                v = v + obj.cast( i );
            end
            v = v >= numel( obj.dimensions ) ./ 2;
        end
    end
    
    properties ( Access = private )
        grid Grid
        mesh(:,3,3) double {mustBeReal,mustBeFinite}
        dimensions(1,:) double
    end
    
    methods ( Access = private )
        function v = cast( obj, dimension )
            p = obj.get_permutation( dimension );
            p_inv = obj.get_inverse_permutation( dimension );
            rr = RayRaster( ...
                obj.grid.data{ p }, ...
                obj.mesh( :, p, : ) ...
                );
            v = permute( rr.voxelize(), p_inv );
        end
        
        function v = voxelize( obj )
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
            p = Voxels.PERMUTATIONS( dimension, : );
        end
        
        function p = get_inverse_permutation( dimension )
            p = Voxels.INVERSE_PERMUTATIONS( dimension, : );
        end
    end
    
end

