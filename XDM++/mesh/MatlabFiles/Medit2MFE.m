%%
%%
%%  Medit2MFE
%%    Convert Medit files from CGAL to Rollett's group's format
%%    for moving finite element 
%%
%%
function Medit2MFE( MeditFile, TriFilename, PointFilename )
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



%fscanf( FidIn, '%s %g')
fclose( FidIn );

WriteTriangleFile( TriFilename,   Triangles );
WriteNodeFile    ( PointFilename, Points );
% write NodeFile

end


function WriteTriangleFile( filename, Tri )

fid = fopen( filename, 'w');

fprintf(fid, '%d\n', length(Tri) );
IDs = [1:length(Tri)];
fprintf(fid, '%d %d %d %d 0 0 0 -1 %d\n', [IDs', Tri]' );
fclose(fid);
end

function WriteNodeFile( filename, Nodes )

fid = fopen( filename, 'w');
IDs = [1:length(Nodes)];
fprintf(fid, '%d\n', length(Nodes) );
fprintf(fid, '%d 0 %g %g %g \n', [IDs', Nodes(:, 1:3)]' );
fclose(fid);


end