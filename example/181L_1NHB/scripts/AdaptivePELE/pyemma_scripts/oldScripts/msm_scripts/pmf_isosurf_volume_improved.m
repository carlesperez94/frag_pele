%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a simple matlab script to plot PMF surfaces from Emma. %
%% Written by D. Lecina                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% T:  Temperature (K) %
%%%%%%%%%%%%%%%%%%%%%%%

function pmf_isosurf(DeltaW, T, rawPoints, pdbFilename, ligandName, printFigures)

if nargin < 6
    printFigures = 0;
end

% Plot concave hull
[R, Rlig] = getPDBCoordinates(pdbFilename, ligandName);


if printFigures == 1
    figure(1)
    scatter3(R(:,1), R(:,2), R(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0.75 0.75])
    hold on
    scatter3(Rlig(:,1), Rlig(:,2), Rlig(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0 1])
    hold on
end

%proteinRadius = input('Protein Region. Radius of concave hull: ');
%explorationRadius = input('Volume definition - Exploration region. Radius of concave hull: ');
proteinRadius = 3;%input('Protein Region. Radius of concave hull: ');
explorationRadius = 3;%input('Volume definition - Exploration region. Radius of concave hull: ');

if printFigures == 1
    close(1)
    %repeated to avoid java exception bug in matlab when working with dual monitors
    figure(1)
    scatter3(R(:,1), R(:,2), R(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0.75 0.75])
    hold on
    scatter3(Rlig(:,1), Rlig(:,2), Rlig(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0 1])
    hold on

    [V,S] = alphavol(R, proteinRadius);

    h = trisurf(S.bnd, R(:,1), R(:,2), R(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.5);
    hold on
    verticesInConcaveHull = get(h, 'vertices');
    facesInConcaveHull = get(h, 'faces');
    facesOut = unifyMeshNormals(facesInConcaveHull, verticesInConcaveHull);
    hold on
end




if printFigures == 1
    figure(2)
    h = trisurf(S.bnd, R(:,1), R(:,2), R(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.5);
    hold on
end

    points = load (rawPoints);
    X = [points(:,1) points(:,2) points(:,3)];
    [V,S] = alphavol(X, explorationRadius);

if printFigures == 1
    h = trisurf(S.bnd, X(:,1), X(:,2), X(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.5);
end

if printFigures == 1
    figure(3)
    scatter3(X(:,1), X(:,2), X(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 1 0])
    hold on
    h = trisurf(S.bnd, X(:,1), X(:,2), X(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.25, 'FaceLighting', 'phong', 'LineWidth', 1.5);
    hold on
    scatter3(Rlig(:,1), Rlig(:,2), Rlig(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0])
    %shading interp;
    %shading faceted
    xlabel('x \AA','interpreter','latex','FontSize',10,'FontName','Times')
    ylabel('y \AA','interpreter','latex','FontSize',10,'FontName','Times')
    zlabel('z \AA','interpreter','latex','FontSize',10,'FontName','Times')
    hold on
end

if printFigures == 1
    figure(4)
    scatter3(Rlig(:,1), Rlig(:,2), Rlig(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0])
    hold on
    h = trisurf(S.bnd, X(:,1), X(:,2), X(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.2);
    hold on
    xlabel('x \AA','interpreter','latex','FontSize',10,'FontName','Times')
    ylabel('y \AA','interpreter','latex','FontSize',10,'FontName','Times')
    zlabel('z \AA','interpreter','latex','FontSize',10,'FontName','Times')
end


disp('In order to speed up the calculations, enter the region to interpolate. The rest will be omitted')
%xmin = input('X min ... ');
%xmax = input('X max ... ');
%ymin = input('Y min ... ');
%ymax = input('Y max ... ');
%zmin = input('Z min ... ');
%zmax = input('Z max ... ');

%if isempty(xmin)
    xmin = floor(min(points(:,1)));
%    xmin = -10;
%end
%if isempty(xmax)
    xmax = ceil(max(points(:,1)));
%    xmax = 10;
%end
%if isempty(ymin)
    ymin = floor(min(points(:,2)));
%    ymin = 10;
%end
%if isempty(ymax)
    ymax = ceil(max(points(:,2)));
%    ymax = 30;
%end
%if isempty(zmin)
    zmin = floor(min(points(:,3)));
%    zmin = 5;
%end
%if isempty(zmax)
    zmax = ceil(max(points(:,3)));
%    zmax = 30;
%end

disp('Domain: ');
fprintf('x:(%f, %f); y:(%f, %f), z:(%f, %f)\n\n', xmin, xmax, ymin, ymax, zmin, zmax);
fprintf('Volume = %f Angstroms**3\n\n', (xmax-xmin)*(ymax-ymin)*(zmax-zmin));


% load 1D data

load 'pmf_xyz.dat'

rawx = pmf_xyz(:,1);
rawy = pmf_xyz(:,2);
rawz = pmf_xyz(:,3);
rawp = pmf_xyz(:,4);

x = [];
y = [];
z = [];
p = [];

for i = 1:size(rawx(:))
    if rawx(i) >= xmin && rawx(i) <= xmax
        if rawy(i) >= ymin && rawy(i) <= ymax
            if rawz(i) >= zmin && rawz(i) <= zmax
                x = [x; rawx(i)];
                y = [y; rawy(i)];
                z = [z; rawz(i)];
                p = [p; rawp(i)];
            end
        end
    end
end

%x = rawx;
%y = rawy;
%z = rawz;
%p = rawp;


if printFigures == 1
    figure(50)
    h = trisurf(S.bnd, X(:,1), X(:,2), X(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.2);
    hold on
    scatter3(Rlig(:,1), Rlig(:,2), Rlig(:,3), 90, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0 0])
    lightangle(-45,30)
    lightangle(45,-30)
    %lightangle(0,0)
    lightangle(178,44)
    hold on
    %scatter3(x_inside, y_inside, z_inside, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 1 1])
    caxis([min(p), max(p)])
    scatter3(x, y, z, 90, p, 'filled')
    colorbar
    xlabel('x \AA','interpreter','latex','FontSize',10,'FontName','Times')
    ylabel('y \AA','interpreter','latex','FontSize',10,'FontName','Times')
    zlabel('z \AA','interpreter','latex','FontSize',10,'FontName','Times')
    hold on
end




%xmin = floor(min(x));
%xmax = ceil(max(x));
%ymin = floor(min(y));
%ymax = ceil(max(y));
%zmin = floor(min(z));
%zmax = ceil(max(z));

% Interpolation
d = input('Enter resolution of meshgrid d  ... ');

dx = xmin: d : xmax;
dy = ymin: d : ymax;
dz = zmin: d : zmax;

disp('number of points:')
disp(size(dx,2)*size(dy,2)*size(dz,2))

%xg, yg and zg are matrices with all the permutations of points
[xg, yg, zg] = meshgrid(dx, dy, dz);



%F = TriScatteredInterp(x, y, z, p, 'natural');
tStart = tic;
%F = scatteredInterpolant(x, y, z, p, 'linear');
F = scatteredInterpolant(x, y, z, p, 'nearest', 'nearest');
hold on

if printFigures == 0
    h = trisurf(S.bnd, X(:,1), X(:,2), X(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.2);
end
verticesInConcaveHull = get(h, 'vertices');
facesInConcaveHull = get(h, 'faces');
facesOut = unifyMeshNormals(facesInConcaveHull, verticesInConcaveHull);
X2 = [xg(:) yg(:) zg(:)];
in = inpolyhedron(facesOut, verticesInConcaveHull, X2, 'FlipNormals', false);

%TODO: Evaluate just in points inside concaveHull
tElapsed = toc(tStart);
pg = F(xg, yg, zg);
fprintf('Interpolation: %d minutes and %f seconds\n',floor(tElapsed/60),rem(tElapsed,60));

min(pg(:))
max(pg(:))


pg_inside = pg;
pg_inside0 = pg;
for i = 1:size(xg(:))
    if in(i) == true
        pg_inside(i) = pg(i);
        pg_inside0(i) = pg(i);
    else
        % := outside
        pg_inside(i) = Inf;
        pg_inside0(i) = 0;
    end
end

dlmwrite ('isosurf2.dat', pg_inside);
save ('isosurf2.mat', 'xg', 'yg', 'zg', 'pg_inside');

x_inside = [];
y_inside = [];
z_inside = [];
g_inside = [];
for i = 1:size(pg_inside, 1)
    for j = 1:size(pg_inside, 2)
        for k = 1:size(pg_inside, 3)
            if pg_inside(i,j,k) ~= Inf 
                x_inside = [x_inside; xmin + d*(j-1)];
                y_inside = [y_inside; ymin + d*(i-1)];
                z_inside = [z_inside; zmin + d*(k-1)];
                g_inside = [g_inside; pg_inside(i,j,k)];
            end
        end
    end
end


%to integrate
min_g_inside = min(g_inside);
for i = 1:size(pg(:))
    pg(i) = pg(i) - min_g_inside;
    pg_inside(i) = pg_inside(i) - min_g_inside;
end
%to plot
for i = 1:size(g_inside(:))
    g_inside(i) = g_inside(i) - min_g_inside;
end

%if printFigures == 1
%    figure(22)
%    %hist(g_inside(:), 20)
%    [counts,centers]  = hist(g_inside(:),20)
%    hist(g_inside(:),20)
%end
min(g_inside(:))
max(g_inside(:))

if printFigures == 1
    figure(5)
    h = trisurf(S.bnd, X(:,1), X(:,2), X(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.2);
    lightangle(-45,30)
    lightangle(45,-30)
    %lightangle(0,0)
    lightangle(178,44)
    hold on
    %scatter3(x_inside, y_inside, z_inside, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 1 1])
    ginsidemin = min(g_inside)
    ginsidemax = max(g_inside)
    caxis([min(g_inside), max(g_inside)])
    scatter3(x_inside, y_inside, z_inside, 45, g_inside, 'filled')
    colorbar
    xlabel('x \AA','interpreter','latex','FontSize',10,'FontName','Times')
    ylabel('y \AA','interpreter','latex','FontSize',10,'FontName','Times')
    zlabel('z \AA','interpreter','latex','FontSize',10,'FontName','Times')
    hold on
end


%disp('Enter boundaries for the integration region')
%xint_min = input('Integration X min ... ');
%xint_max = input('Integration X max ... ');
%yint_min = input('Integration Y min ... ');
%yint_max = input('Integration Y max ... ');
%zint_min = input('Integration Z min ... ');
%zint_max = input('Integragion Z max ... ');
%
%if isempty(xint_min)
%    xint_min = floor(min(points(:,1)));
%%    xint_min = -10;
%end
%if isempty(xint_max)
%    xint_max = ceil(max(points(:,1)));
%    xint_max = -10;
%end
%if isempty(yint_min)
%    yint_min = floor(min(points(:,2)));
%%    yint_min = 10;
%end
%if isempty(yint_max)
%    yint_max = ceil(max(points(:,2)));
%%    yint_max = 30;
%end
%if isempty(zint_min)
%    zint_min = floor(min(points(:,3)));
%%    zint_min = 5;
%end
%if isempty(zint_max)
%    zint_max = ceil(max(points(:,3)));
%%    zint_max = 30;
%end
    xint_min = floor(min(points(:,1)));
    xint_max = ceil(max(points(:,1)));
    yint_min = floor(min(points(:,2)));
    yint_max = ceil(max(points(:,2)));
    zint_min = floor(min(points(:,3)));
    zint_max = ceil(max(points(:,3)));

disp('Integration Domain: ');
fprintf('x:(%f, %f); y:(%f, %f), z:(%f, %f)\n\n', xint_min, xint_max, yint_min, yint_max, zint_min, zint_max);

for i = 1:size(pg_inside(:))
    if xg(i) < xint_min || xg(i) > xint_max
        pg_inside(i) = Inf;
    end
    if yg(i) < yint_min || yg(i) > yint_max
        pg_inside(i) = Inf;
    end
    if zg(i) < zint_min || zg(i) > zint_max
        pg_inside(i) = Inf;
    end
end


disp('*** Volme integration of isosurface ***')
disp('Enter the range of G_pmf values to consider. Remind that all the points with larger isovalue will be ommitted')
minIsovalue = input('Min isovalue ... ');
maxIsovalue = input('Max isovalue ... ');
intervalIsovalue = input('Interval ... ');

arrayOfIsovalues = minIsovalue:intervalIsovalue:maxIsovalue;
solutions = [];

tStart = tic;
minValue = min(p);

disp('  ---------------------------------------------------------------------------------');
disp('        boundary            G     kb*T*log(Vb/Vo)      deltaW      integral ');
disp('  ---------------------------------------------------------------------------------');

integralResults = [];
deltaGEstimations = [];
endmin = min(pg_inside(:))
endmax = max(pg_inside(:))
for isovalue = arrayOfIsovalues
    integralResult = integrateForIsoValue(pg_inside, isovalue, d, 1/(T*0.0019872041));
    integralResults = [integralResults integralResult];
    deltaGEstimation = -0.0019872041*T*log(integralResult/1661) - (DeltaW);
    deltaGEstimations = [deltaGEstimations deltaGEstimation];

    fprintf('\t%f\t%f\t%f\t%f\t%f\n', isovalue, deltaGEstimation, -0.0019872041*T*log(integralResult/1661), DeltaW, integralResult);
end
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));

disp('  ------------------------------------------------------');


figure(10)
scatter(arrayOfIsovalues, integralResults);
figure(11)
scatter(arrayOfIsovalues, deltaGEstimations);


isovalue = input('Choose isovalue to plot ');
plot_isosurface(2, 4, isovalue, pg_inside, xg, yg, zg);

%scatter3(X_in(:,1),X_in(:,2), X_in(:,3), 50, X_in(:,4), 'fill'); %add points to plot
caxis([floor(min(pg_inside(:))), ceil(max(pg_inside(:)))])
title(sprintf('Radius = %.1g',radius), 'Fontsize',10,'FontName','Times')
hold on
%hold off

return


%%%%%%%
% End %
%%%%%%%
