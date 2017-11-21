function WriteMeshVtk( filename, Tsurf, TPoints )



TSurfProp  = Tsurf(:, 4);
TPointProp = TPoints(:, 4);

fd = fopen( filename, 'w');
fprintf(fd, '# vtk DataFile Version 2.0\n' );
fprintf(fd, 'Matlab converted vtk file from mesh\n' );
fprintf(fd, 'ASCII\n' );
fprintf(fd, 'DATASET POLYDATA\n' );
fprintf(fd, 'POINTS %i float \n', length( TPoints ) );
fprintf(fd, '%g %g %g\n', TPoints(:, 1:3)' );
fprintf(fd, 'POLYGONS %i %i\n', length(Tsurf), length(Tsurf) * 4 );
fprintf(fd, '3 %i %i %i\n', (Tsurf(:, 1:3) -1 )' );
fprintf(fd, 'CELL_DATA %i\n', length(Tsurf) );
fprintf(fd, 'SCALARS SurfaceID int 1\n');
fprintf(fd, 'LOOKUP_TABLE default\n');
fprintf(fd, '%i\n', Tsurf(:, 4)' );
fclose(fd);
end
