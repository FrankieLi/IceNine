%%
%%
%%  Make boundary points so we can test stuff.
%%
%%
%%
function bnd = MakeBndPoints()

r = 31;

z = [ [-r : 1 : -1 ], 0, [ 1 : 1 : r ] ]';

zeros_ = zeros( length( z ), 1 );
wc = ones( length(z), 1 );
wc(1) = 2.1;
wc(end) = 2.1;
central_line = [ zeros_, zeros_, z, wc ];

theta = [0: 2*pi/90: 2 *pi]'; 
theta = theta(2:end-1);
zeros_ = zeros( length( theta ), 1 );
S      =  r * cos( theta );
z_     =  r * sin( theta );
X = [  S,  zeros_, z_ ];
Y = [  zeros_,  S, z_ ];

bnd = [ X; Y ];

w = ones( length( bnd ), 1 ) * 1;
bnd = [bnd, w];
bnd = [central_line(1, :);central_line(end, :); central_line(2:end-1, :); bnd];
bnd(:, 1) = bnd(:, 1) + r;
bnd(:, 2) = bnd(:, 2) + r;
bnd(:, 3) = bnd(:, 3) + r;
end