function plot_isosurface(proteinFigure, volumeFigure, isovalue, matrix3d, x, y, z)

%remove nan's andd smoothen boundaries
for i = 1:size(matrix3d, 1)
    for j = 1:size(matrix3d, 2)
        for k = 1:size(matrix3d, 3)
            if matrix3d(i,j,k) == 0 | isnan(matrix3d(i,j,k)) | matrix3d(i,j,k) > isovalue
                matrix3d(i,j,k) = isovalue + 1;
            end
        end
    end
end

figure(proteinFigure)
matrix3d = smooth3(matrix3d, 'gaussian', 1);
s = patch(isosurface(x, y, z, matrix3d, isovalue));
n = isonormals(x, y, z, matrix3d, s);

grid on
%isocolors(x, y, z, matrix3d, s)
%shading interp;
%light;
%camlight 
%camlight left; camlight right
lightangle(-45,30)
lightangle(45,-30)
%lightangle(0,0)
lightangle(178,44)
%lighting gouraud
lighting phong;
%material shiny
%material metal

%alpha(0.5);

color = [0 0.75 0.75];
set(s,'FaceColor', color,'EdgeColor','none');
title(strcat('Isovalue: ',num2str(isovalue)), 'Fontsize',10,'FontName','Times')
xlabel('x \AA','interpreter','latex','FontSize',10,'FontName','Times')
ylabel('y \AA','interpreter','latex','FontSize',10,'FontName','Times')
zlabel('z \AA','interpreter','latex','FontSize',10,'FontName','Times')
hold on

figure(volumeFigure)
s = patch(isosurface(x, y, z, matrix3d, isovalue));
n = isonormals(x, y, z, matrix3d, s);
set(s,'FaceColor', color,'EdgeColor','none');
title(strcat('Isovalue: ',num2str(isovalue)), 'Fontsize',10,'FontName','Times')
xlabel('x \AA','interpreter','latex','FontSize',10,'FontName','Times')
ylabel('y \AA','interpreter','latex','FontSize',10,'FontName','Times')
zlabel('z \AA','interpreter','latex','FontSize',10,'FontName','Times')
hold on
grid on
lightangle(-45,30)
lightangle(45,-30)
lightangle(178,44)
%lighting gouraud
lighting phong;
