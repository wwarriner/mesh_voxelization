
%Plot the original STL mesh:
figure
[stlcoords] = READ_stl('sample.stl');
xco = squeeze( stlcoords(:,1,:) )';
yco = squeeze( stlcoords(:,2,:) )';
zco = squeeze( stlcoords(:,3,:) )';
[hpat] = patch(xco,yco,zco,'b');
axis equal

%Voxelise the STL:
grid = Grid( linspace(-15,15,100), linspace(-15,15,100), linspace(-15,15,100) );
stl = READ_stl( "sample.stl" );
[ fv.faces, fv.vertices ] = CONVERT_meshformat( stl );
%OUTPUTgrid = VOXELISE(grid,fv,'xyz');
VV = Voxels( grid, fv, 'xyz' );
OUTPUTgrid = VV.generate();

%Show the voxelised result:
figure;
subplot(1,3,1);
imagesc(squeeze(sum(OUTPUTgrid,1)));
colormap(gray(256));
xlabel('Z-direction');
ylabel('Y-direction');
axis equal tight

subplot(1,3,2);
imagesc(squeeze(sum(OUTPUTgrid,2)));
colormap(gray(256));
xlabel('Z-direction');
ylabel('X-direction');
axis equal tight

subplot(1,3,3);
imagesc(squeeze(sum(OUTPUTgrid,3)));
colormap(gray(256));
xlabel('Y-direction');
ylabel('X-direction');
axis equal tight
