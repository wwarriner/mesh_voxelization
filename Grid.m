classdef Grid < handle
    
    properties ( SetAccess = private )
        points(1,:) cell
    end
    
    properties ( SetAccess = private, Dependent )
        dimension_count(1,1) uint32
        shape(1,:) uint32
        origin(1,:) double
        envelope Envelope
    end
    
    methods
        function obj = Grid( varargin )
            assert( 1 <= nargin );
            if iscell( varargin{ 1 } )
                in = varargin{ 1 };
            else
                in = varargin;
            end
            
            for i = 1 : nargin
                v = in{ i };
                assert( isa( v, 'double' ) );
                assert( 1 <= numel( v ) );
                assert( isvector( v ) );
                assert( isreal( v ) );
                assert( all( isfinite( v ) ) );
                if 1 < numel( v )
                    assert( all( 0 < diff( v ) ) );
                end
                in{ i } = v( : );
            end
            obj.points = in;
        end
        
        function value = create_zeros( obj, type )
            value = zeros( obj.shape, type );
        end
        
        function value = get.dimension_count( obj )
            value = numel( obj.points );
        end
        
        function value = get.shape( obj )
            value = nan( 1, obj.dimension_count );
            for i = 1 : obj.dimension_count
                value( i ) = numel( obj.points{ i } );
            end
            assert( ~any( isnan( value ) ) );
        end
        
        function value = get.origin( obj )
            value = nan( 1, obj.dimension_count );
            for i = 1 : obj.dimension_count
                d = obj.points{ i };
                value( i ) = d( 1 );
            end
            assert( ~any( isnan( value ) ) );
        end
        
        function value = get.envelope( obj )
            value = Envelope( obj.origin, obj.get_end() );
        end
    end
    
    methods ( Access = private )
        function value = get_end( obj )
            value = nan( 1, obj.dimension_count );
            for i = 1 : obj.dimension_count
                d = obj.points{ i };
                value( i ) = d( end );
            end
            assert( ~any( isnan( value ) ) );
        end
    end
    
end

