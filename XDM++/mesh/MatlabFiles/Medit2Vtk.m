function Medit2vtk(  MeditFile, VtkFile )


FidIn= fopen( MeditFile, 'r' );

fgetl( FidIn)
fgetl( FidIn)
fgetl( FidIn)
numPoints = str2num( fgetl( FidIn) )

Points = [];

Points = fscanf( FidIn, '%g %g %g %d ', [ 4, numPoints ])';
fgetl( FidIn ) % triangles
numTriangles = str2num( fgetl( FidIn ) )

Triangles = [];

Triangles = fscanf( FidIn, '%g %g %g %d ', [4, numTriangles ])';
v = randperm( max( Triangles(:, 4)));
Triangles(:, 4) = v(Triangles(:, 4));
WriteMeshVtk( VtkFile, Triangles, Points );

%fscanf( FidIn, '%s %g')
fclose( FidIn );
end