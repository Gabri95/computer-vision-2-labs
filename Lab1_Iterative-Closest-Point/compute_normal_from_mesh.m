function [normals] = compute_normal_from_mesh(cloud)
    
    
DT = delaunayTriangulation(squeeze(cloud(1:2, :))');

TR = triangulation(DT.ConnectivityList, cloud');
normals = vertexNormal(TR);

trisurf(TR,'FaceColor',[0.8 0.8 1.0]);
axis equal
hold on
quiver3(TR.Points(:,1),TR.Points(:,2),TR.Points(:,3), ...
     normals(:,1),normals(:,2),normals(:,3),0.5,'Color','b');
hold off

normals = normals';
end