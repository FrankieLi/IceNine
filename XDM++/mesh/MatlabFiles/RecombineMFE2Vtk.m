%%
%%
%%
%%
%%
function RecombineMFE2Vtk( MeditFile, MFEFile, VtkFile )

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

NewPoints = ReadMFEOutput( MFEFile );
WriteMeshVtk( VtkFile, Triangles, NewPoints );
end


function NewPoints = ReadMFEOutput( MFEFile )

fid = fopen( MFEFile, 'r');

nPoints = fscanf( fid, '%d', [1, 1]);

NewPoints = fscanf( fid, '%d %d %g %g %g ', [5, nPoints])';

NewPoints = [NewPoints(:, 3:5), NewPoints(:, 2) ];
fclose( fid );
end