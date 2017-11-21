function SphereDxGenerator( meshDim, radius, filename )


center = [meshDim, meshDim, meshDim]/2;
DxOut = zeros( meshDim, meshDim, meshDim ) + 1;


v = [1:meshDim*meshDim*meshDim];

[Ix, Iy, Iz] = ind2sub(  [meshDim, meshDim, meshDim], v );

test = [Ix' - center(1), Iy' - center(2), Iz' - center(3) ];

r = sqrt( dot(test, test, 2) );


DxOut(v(r < radius )) = 2;
 
 DxOut(v( test(:, 1) < 0 & test(:, 2) < 0 & test(:, 3) < 0 )) = 3; 
 DxOut(v( test(:, 1) < 0 & test(:, 2) < 0 & test(:, 3) > 0 )) = 4; 
 DxOut(v( test(:, 1) < 0 & test(:, 2) > 0 & test(:, 3) > 0 )) = 5; 
DxOut(v(r > radius )) = 1;

size( test ) 
unique(DxOut)
dlmwrite( filename, DxOut, ' ');
end
