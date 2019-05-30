function integral = integrateForIsoValue(matrix, isovalue, d, beta)
%%% matrix is a 3D matrix
for i = 1:size(matrix, 1)
    for j = 1:size(matrix, 2)
        for k = 1:size(matrix, 3)
            if matrix(i,j,k) > isovalue | isnan(matrix(i,j,k)) 
                matrix(i,j,k) = Inf;
            end
        end
    end
end


imin = 1+1;
imax = size(matrix,1) - 1;
jmin = 1+1;
jmax = size(matrix,2) - 1;
kmin = 1+1;
kmax = size(matrix,3) - 1;

integral = 0;
for i = imin:imax
    for j = jmin:jmax
        for k = kmin:kmax
            %average of the 8 closest points
            partialResult = 0;
            allInf = 1;
            for ii = i:i
                for jj = j:j
                    for kk = k:k
                        %exp(-beta * matrix(ii, jj, kk))
                        if matrix(ii,jj,kk) ~= Inf
                            %partialResult = partialResult + matrix(ii,jj,kk);% + exp(- beta * matrix(ii,jj,kk));
                            partialResult = partialResult + exp(- beta * matrix(ii,jj,kk));
                            %partialResult = partialResult + matrix(ii,jj,kk);
                            allInf = 0;
                        end
                    end
                end
            end
            if allInf == 0
                %integral = integral + partialResult/27.;
                integral = integral + partialResult;
                %integral = integral + exp(- beta * partialResult/27.);
            end
        end
    end
end
integral = integral*d^3;


