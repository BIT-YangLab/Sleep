function dist_transform = generate_geodesic_distance_mat(mesh, hemi)
% Copy from CBIG
    % usage: generate_geodesic_distance_mat('fsaverage5/fsaverage6/fsaverage' or 'fs_LR_32k'. , 'lh');
    % example: lh_dist = generate_geodesic_distance_mat('fs_LR_32k', 'lh');
    if(~isempty(strfind(mesh, 'fsaverage')))
        MARS_sbjMesh = CBIG_ReadNCAvgMesh(hemi, mesh, 'inflated','cortex');
    elseif(~isempty(strfind(mesh, 'fs_LR')))
        MARS_sbjMesh = CBIG_read_fslr_surface(hemi, mesh, 'inflated', 'medialwall.annot');
    end
    num_vertices = int32(size(MARS_sbjMesh.vertices, 2));
    maxNeighbors = int32(size(MARS_sbjMesh.vertexNbors, 1));
    dist_transform = zeros(num_vertices, num_vertices);

    vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(MARS_sbjMesh.vertexNbors, 1)),...
        int32(size(MARS_sbjMesh.vertices, 2)), int32(MARS_sbjMesh.vertexNbors), single(MARS_sbjMesh.vertices));
    vertexDist2Nbors = sqrt(vertexDistSq2Nbors);

    for i = 1:num_vertices
        
        dist_transform(i, :) = MARS_DT_Boundary(int32([1:1:num_vertices] == i), num_vertices,...
            maxNeighbors, MARS_sbjMesh.vertexNbors, double(vertexDist2Nbors)); %min_heap assumes double
        
    end
end