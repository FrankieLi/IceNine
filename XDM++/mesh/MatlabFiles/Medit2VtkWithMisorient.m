function Medit2edit2VtkWithMisorient(  MeditFile, VtkFile, ID2BndFile, OrientationFile)


%%
%%  ReID, Calculate misorientations
%%
ID2BndMap = load( ID2BndFile );
OrientationMap = load( OrientationFile );
ID2BndMap(:, 2:3) = ID2BndMap(:, 2:3) + 1;
findvec = find( ID2BndMap(:, 2) ~= 0);

QOrientList = QuatOfRMat( RMatOfBunge( OrientationMap', 'degrees'));
Misorient = [];
if( size(findvec, 1) > 10000 )
  
  nStep = floor( size(findvec, 1) / 10000 );
  subRegionInd = [1:10000: nStep * 10000 ];
  subRegionInd = [ subRegionInd, size(findvec, 1) + 1];
  for n = 2:length( subRegionInd )
    s1 = subRegionInd( n - 1 );
    s2 = subRegionInd( n ) -1;
    
    fSub = findvec( s1:s2 );
    Misorient = [ Misorient, Misorientation( QOrientList( :, ID2BndMap( fSub, 2 ) ),...
      QOrientList( :, ID2BndMap( fSub, 3 ) ),...
      CubSymmetries() ) ];
  end
else
   Misorient = Misorientation( QOrientList( :, ID2BndMap( fSub, 2 ) ),...
      QOrientList( :, ID2BndMap( fSub, 3 ) ),...
      CubSymmetries() );
end
                         
% ID2MisorientMap = [ID2BndMap(:, 1), ones( length( ID2BndMap ), 1) * 60 ];
% ID2MisorientMap( findvec, 2) = Misorient' * 180 / pi;

ID2MisorientMap = zeros(max(ID2BndMap(:, 1)), 1);

ID2MisorientMap( ID2BndMap(findvec, 1) ) = Misorient' * 180 / pi;

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

nPreStartIndex = min(Triangles(:, 4)) -1;
Triangles(:, 4) = ID2MisorientMap( Triangles(:, 4) - nPreStartIndex  );
WriteMeshVtk( VtkFile, Triangles, Points, 2, 'Misorientation' );

%fscanf( FidIn, '%s %g')
fclose( FidIn );
end