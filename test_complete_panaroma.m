    clear all;
    clc;
    
    vl_setup;
    
    path = "/home/kalyan/Desktop/assignment3/west_campus1/";
    all_images = dir(strcat(path, "*.JPG"));
    
    i1 = imread(strcat(path, all_images(1).name));
    out_temp = i1;
    
    h_mat_final = [1 0 0; 0 1 0; 0 0 1];
    
    for main_iterator = 1:(length(all_images) - 1)
      
      i2 = imread(strcat(path, all_images(main_iterator + 1).name));
      
      i1_gray = single(rgb2gray(i1));
      i2_gray = single(rgb2gray(i2));
    
      [f_i1, d_i1] = vl_sift(i1_gray);
      [f_i2, d_i2] = vl_sift(i2_gray);
    
      d_i1 = double(d_i1);
      d_i2 = double(d_i2);
    
      n_d_i1 = diag(sqrt((d_i1')*d_i1));
      n_d_i2 = diag(sqrt((d_i2')*d_i2));
      
      norm_d_i1 = zeros(size(d_i1));
      norm_d_i2 = zeros(size(d_i2));
    
      for n_iterator = 1:size(d_i1, 2)
        norm_d_i1(:, n_iterator) = d_i1(:, n_iterator)./n_d_i1(n_iterator);
      endfor
    
      for n_iterator2 = 1:size(d_i2, 2)
        norm_d_i2(:, n_iterator2) = d_i2(:, n_iterator2)./n_d_i2(n_iterator2);
      endfor
    
      cp_mat = zeros(size(max(size(d_i1, 2), size(d_i2, 2))), 2);
    
      num_match_p = 1;
      for d_iterator = 1:size(d_i1, 2)
      
        sym_test_mat = (norm_d_i1(:, d_iterator)' * norm_d_i2);
        [sorted_mat, indices_mat] = sort(acos(sym_test_mat), 'ascend');
        ratio_cp = sorted_mat(1)/sorted_mat(2);
        if (ratio_cp < 0.8)
          cp_mat(num_match_p, 1) = d_iterator;
          cp_mat(num_match_p, 2) = indices_mat(1);
          num_match_p = num_match_p + 1;
        else
      
        endif
    
      endfor
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %x_prime - first image xcord
      %x - second image xcord
      %y_prime - first image ycord
      %y - secnd image ycord
      %Left Image - I1
      %Right Image - I2
      %%%RANSAC%%%%%%%%%%%%%%%%%%%%%
    
      iter = 200;
      p_to_c = min(size(cp_mat, 1), 10);
      past_error = 100000;
      past_n_inliers = 0;
      h_ransac = zeros(3, 3);
      k = 20000000;
    
      for i = 1:iter
    
        a = 1;
        b = size(cp_mat, 1);
        rand_ind = (b - a).*rand(p_to_c, 1) + a;
        rand_ind = round(rand_ind);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%
      %Computing homography matrix from 1 to 0
      %Destination: given by Right
      %Source: Left
      %%%%%%%%%%%%%%%%%%%%%%%%%%
      
        x_prime_for_H = f_i1(1, cp_mat(rand_ind, 1))';
        y_prime_for_H = f_i1(2, cp_mat(rand_ind, 1))';
      
        x_for_H = f_i2(1, cp_mat(rand_ind, 2))';
        y_for_H = f_i2(2, cp_mat(rand_ind, 2))';
        
        if(numel(x_prime_for_H) >= 4)
        
          h_mat = compute_homography(x_prime_for_H, y_prime_for_H, x_for_H, y_for_H);
          inlier_ind = 1;
      
        inlier_mat = zeros(size(cp_mat, 1), 1);
        counter = 1;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%
      %finding inliers
      %%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:b
      
          x_cord_s = f_i2(1, cp_mat(j, 2));
          y_cord_s = f_i2(2, cp_mat(j, 2));
      
          x_cord_d = f_i1(1, cp_mat(j, 1));
          y_cord_d = f_i1(2, cp_mat(j, 1));
      
      % H * source coordinates = destination coordinates
          l = h_mat * [x_cord_s; y_cord_s; 1];
          x_cord_dash = l(1)/l(3);
          y_cord_dash = l(2)/l(3);
        
      
          if(abs((x_cord_dash - x_cord_d)) + abs((y_cord_dash - y_cord_d)) < 6)
          %if(sqrt((floor(x_cord_d) - floor(x_cord_dash))^2 + (floor(y_cord_d) - floor(y_cord_dash))^2) < 36)
            inlier_mat(j) = 1;
            inlier_ind = inlier_ind + 1;
          endif
      
        endfor
      
      
      if(past_n_inliers < inlier_ind)
        past_n_inliers = inlier_ind;
        h_ransac = h_mat;
        if(inlier_ind > 1)
          inliers_ransac = inlier_mat;
          
        endif
      endif
    endif
    endfor
    
    %%Recomputing the homography for the set of inliers
    non_zero_inliers = find(inliers_ransac);
    x_prime = f_i1(1, cp_mat(non_zero_inliers, 1))';
    y_prime = f_i1(2, cp_mat(non_zero_inliers, 1))';
    x = f_i2(1, cp_mat(non_zero_inliers, 2))';
    y = f_i2(2, cp_mat(non_zero_inliers, 2))';
    
    h_mat_final = compute_homography(x_prime, y_prime, x, y);
    
    transformStructure_Projective = maketform('projective', h_mat_final');
    img_transformed = imtransform(i2,transformStructure_Projective);
    
    [m2, n2 ~] = size(i2);
    [m1, n1, dim] = size(out_temp);
    
    edge_indices = zeros(3, 4);
    edge_indices(:, 1) = h_mat_final * [1; 1; 1];
    edge_indices(:, 2) = h_mat_final * [n2; 1; 1];
    edge_indices(:, 3) = h_mat_final * [n2; m2; 1];
    edge_indices(:, 4) = h_mat_final * [1; m2; 1];
    
    edge_x_coords = edge_indices(1, :)./edge_indices(3, :);
    edge_y_coords = edge_indices(2, :)./edge_indices(3, :);
    
    vert_var = round(min(edge_y_coords));
    vert_off = 0;
    if (vert_var <= 0)
      vert_off = -vert_var + 1;
      vert_var = 1;
    endif
    
    hor_var = round(min(edge_x_coords));
    hor_off = 0;
    if (hor_var <= 0)
      hor_off = -hor_var + 1;
      hor_var = 1;
    endif
    
    [m3 n3 ~] = size(img_transformed);
    stitch_temp(vert_var:vert_var + m3 - 1, hor_var: hor_var + n3 - 1, :) = img_transformed;
    stitch_temp(vert_off + 1: vert_off + m1, hor_off + 1: hor_off + n1, :) = out_temp;
    
    out_temp = stitch_temp;
    i1 = stitch_temp;
    hmat_mat(:, :, main_iterator) = h_mat_final;
    
   endfor