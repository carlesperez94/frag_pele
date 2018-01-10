function correct_kmeans_stationary_distribution(rawPoints, pdbFilename)
    gridSize = 1.;

    load 'matrix/stationarydistribution.dat'
%% store the value of stationary distribution ...

    p = stationarydistribution(:,1);
    load 'discretized/clusterCenters.dat'
%% store the value of x coordinate

    x = clusterCenters(:,1);
    y = clusterCenters(:,2);
    z = clusterCenters(:,3);


    %% create mesh
    points = load (rawPoints);
    xmin = floor(min(points(:,1)));
    xmax = ceil(max(points(:,1)));
    ymin = floor(min(points(:,2)));
    ymax = ceil(max(points(:,2)));
    zmin = floor(min(points(:,3)));
    zmax = ceil(max(points(:,3)));
    dx = xmin: gridSize : xmax;
    dy = ymin: gridSize : ymax;
    dz = zmin: gridSize : zmax;
    [xg, yg, zg] = meshgrid(dx, dy, dz);

    %%% create surface
    explorationRadius = 3; %%% concave hull radii
    X = [points(:,1) points(:,2) points(:,3)];
    [V,S] = alphavol(X, explorationRadius);
    %plots surface
    figure(1)
    h = trisurf(S.bnd, X(:,1), X(:,2), X(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.15);
    hold on
    %%% end create surface
    %%% unify normals
    verticesInConcaveHull = get(h, 'vertices');
    facesInConcaveHull = get(h, 'faces');
    facesOut = unifyMeshNormals(facesInConcaveHull, verticesInConcaveHull);
    %%% end unify normals
    %%% check whether points are in or not
    X2 = [xg(:) yg(:) zg(:)];
    in = inpolyhedron(facesOut, verticesInConcaveHull, X2, 'FlipNormals', false);
    %%% end check
   


    %% There's probably a much clever way to do it
    clusterNumber = [];
    for i = 1:size(x(:))
        clusterNumber = [clusterNumber; i];
    end

    %%%%% CALCULATING VOLUME %%%%%
    volumeF = scatteredInterpolant(x, y, z, clusterNumber, 'nearest', 'nearest');
    interpolatedClusterNumber = volumeF(xg, yg, zg);
    for i = 1:size(xg(:))
        if in(i) == false
            interpolatedClusterNumber(i) = 0;
        end
    end

    x_inside = [];
    y_inside = [];
    z_inside = [];

    interpolatedClusterNumberInside = [];
    for i = 1:size(interpolatedClusterNumber, 1)
        for j = 1:size(interpolatedClusterNumber, 2)
            for k = 1:size(interpolatedClusterNumber, 3)
                if interpolatedClusterNumber(i,j,k) ~= 0 
                    x_inside = [x_inside; xmin + gridSize*(j-1)];
                    y_inside = [y_inside; ymin + gridSize*(i-1)];
                    z_inside = [z_inside; zmin + gridSize*(k-1)];
                    interpolatedClusterNumberInside = [interpolatedClusterNumberInside; interpolatedClusterNumber(i,j,k)];
                end
            end
        end
    end


    volumeOfClusters = zeros(size(x(:)));
    for i = 1:size(interpolatedClusterNumberInside(:))
        cluster = interpolatedClusterNumberInside(i);
        volumeOfClusters(cluster) = volumeOfClusters(cluster) + 1;
    end
    for i = 1:size(volumeOfClusters) 
        volumeOfClusters(i) = volumeOfClusters(i)*gridSize^3;
    end
    %%%%% END CALCULATING VOLUME %%%%%

    %%%print cluster sizes
    caxis([min(volumeOfClusters), max(volumeOfClusters)])
    scatter3(x, y, z, 90, volumeOfClusters, 'filled')
    colorbar

    %%% Normalise with volume
    correctedP = zeros(size(p(:)));
    for i = 1:size(volumeOfClusters) 
        correctedP(i) = p(i)/volumeOfClusters(i);
    end


    %%% Calculate G, and set minG = 0
    kbT = 0.0019872041*300;
    g = - kbT*log(p);
    correctedG = - kbT*log(correctedP);

    g = g - min(g);
    correctedG = correctedG - min(correctedG);
    
    %in order to use it in the pipeline like previous scripts
    %%% Required by pmf.sh
    
    
    save pmf_1d.dat correctedG -ascii
    save x_coord.dat x -ascii
    save y_coord.dat y -ascii
    save z_coord.dat z -ascii

    %%% PLOT PMF FOR THE CORRECTED VALUES
    figure(2)
    title('Normalised kmeans to same volume')
    h1 = plot(x,correctedG,'ko');
    hold on
    h2 = plot(y,correctedG,'r+');
    h3 = plot(z,correctedG,'g*');
    xlabel('x, y, z \AA','interpreter','latex','FontSize',10,'FontName','Times')
    ylabel('G_{pmf} (kcal/mol)','FontSize',10,'FontName','Times')
    
    legend([h1 h2 h3],{'x','y','z'},3)

    %%% PLOT PMF FOR THE ORIGINAL PMF VALUES
    figure(3)
    title('Raw kmeans')
    h1 = plot(x,g,'ko');
    hold on
    h2 = plot(y,g,'r+');
    h3 = plot(z,g,'g*');
    xlabel('x, y, z \AA','interpreter','latex','FontSize',10,'FontName','Times')
    ylabel('G_{pmf} (kcal/mol)','FontSize',10,'FontName','Times')

    legend([h1 h2 h3],{'x','y','z'},3)

    return


    %%% Plot PMF
    figure(4)
    title('Corrected G')
    h = trisurf(S.bnd, X(:,1), X(:,2), X(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.15);
    hold on
    %%% end create surface
    %%% unify normals
    verticesInConcaveHull = get(h, 'vertices');
    facesInConcaveHull = get(h, 'faces');
    facesOut = unifyMeshNormals(facesInConcaveHull, verticesInConcaveHull);
    %%% end unify normals
    %%% check whether points are in or not
    X2 = [xg(:) yg(:) zg(:)];
    in = inpolyhedron(facesOut, verticesInConcaveHull, X2, 'FlipNormals', false);
    %%% end check
    caxis([min(correctedG), max(correctedG)])
    scatter3(x, y, z, 90, correctedG, 'filled')
    colorbar



    figure(5)
    title('Kmeans Gpmf')
    h = trisurf(S.bnd, X(:,1), X(:,2), X(:,3), 'FaceColor', 'blue', 'FaceAlpha', 0.15);
    hold on
    %%% end create surface
    %%% unify normals
    verticesInConcaveHull = get(h, 'vertices');
    facesInConcaveHull = get(h, 'faces');
    facesOut = unifyMeshNormals(facesInConcaveHull, verticesInConcaveHull);
    %%% end unify normals
    %%% check whether points are in or not
    X2 = [xg(:) yg(:) zg(:)];
    in = inpolyhedron(facesOut, verticesInConcaveHull, X2, 'FlipNormals', false);
    %%% end check
    caxis([min(g), max(g)])
    scatter3(x, y, z, 90, g, 'filled')
    colorbar



