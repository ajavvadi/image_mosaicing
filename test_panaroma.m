clear all;
clc;

vl_setup;

%path = "/home/kalyan/Desktop/assignment3/family_house/";
%all_images = dir("/home/kalyan/Desktop/assignment3/family_house/*.JPG");
%
%for main_iterator = 1:(length(all_images) - 1)
%
%  i1 = imread(all_images(main_iterator).name);
%  i2 = imread(all_images(main_iterator + 1).name);
%  
%  
%
%endfor

i1 = imread('/home/kalyan/Desktop/assignment3/family_house/1.JPG');
i2 = imread('/home/kalyan/Desktop/assignment3/family_house/2.JPG');

%i1 = imread('/home/kalyan/Desktop/assignment3/code/Input_Left.jpg');

%H_input = [1 0 0; 0 10 0; 0 0 1];

%transformStructure_Projective = maketform('projective', H_input');
%i2 = imtransform(i1,transformStructure_Projective);

%i2 = imread('/home/kalyan/Desktop/assignment3/code/Right.jpg');

i1_gray = single(rgb2gray(i1));
i2_gray = single(rgb2gray(i2));

[f_i1, d_i1] = vl_sift(i1_gray);
[f_i2, d_i2] = vl_sift(i2_gray);

%[matches, scores] = vl_ubcmatch(d_i1, d_i2);

d_i1 = double(d_i1);
d_i2 = double(d_i2);

n_d_i1 = diag(sqrt((d_i1')*d_i1));
n_d_i2 = diag(sqrt((d_i2')*d_i2));

for n_iterator = 1:size(d_i1, 2)
  norm_d_i1(:, n_iterator) = d_i1(:, n_iterator)./n_d_i1(n_iterator);
endfor

for n_iterator = 1:size(d_i2, 2)
  norm_d_i2(:, n_iterator) = d_i2(:, n_iterator)./n_d_i2(n_iterator);
endfor

cp_mat = zeros(size(max(size(d_i1, 2), size(d_i2, 2))), 1);

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
  
  end

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
%inliers_ransac = zeros(p_to_c, );
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
  
  x_prime = f_i1(1, cp_mat(rand_ind, 1))';
  y_prime = f_i1(2, cp_mat(rand_ind, 1))';
  
  x = f_i2(1, cp_mat(rand_ind, 2))';
  y = f_i2(2, cp_mat(rand_ind, 2))';

  h_mat = compute_homography(x_prime, y_prime, x, y);
  
  inlier_ind = 1;
  error_term = 0;
  
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
  
  error_term;
  
  if(past_n_inliers < inlier_ind)
    past_n_inliers = inlier_ind;
    h_ransac = h_mat;
    if(inlier_ind > 1)
      inliers_ransac = inlier_mat;
      
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
img21 = imtransform(i2,transformStructure_Projective);

H = h_mat_final;

[M1 N1 dim] = size(i1);
[M2 N2 ~] = size(i2);
[M3 N3 ~] = size(img21);

pt = zeros(3,4);
pt(:,1) = H*[1;1;1];
pt(:,2) = H*[N2;1;1];
pt(:,3) = H*[N2;M2;1];
pt(:,4) = H*[1;M2;1];
x2 = pt(1,:)./pt(3,:);
y2 = pt(2,:)./pt(3,:);

up = round(min(y2));
Yoffset = 0;
if up <= 0
	Yoffset = -up+1;
	up = 1;
end

left = round(min(x2));
Xoffset = 0;
if left<=0
	Xoffset = -left+1;
	left = 1;
end

[M3 N3 ~] = size(img21);
imgout(up:up+M3-1,left:left+N3-1,:) = img21;
	% img1 is above img21
imgout(Yoffset+1:Yoffset+M1,Xoffset+1:Xoffset+N1,:) = i1;








%  x_cord_s_test = f_i2(1, cp_mat(non_zero_inliers, 2));
%  y_cord_s_test = f_i2(2, cp_mat(non_zero_inliers, 2));
%  
%  x_cord_d_test = f_i1(1, cp_mat(non_zero_inliers, 1));
%  y_cord_d_test = f_i1(2, cp_mat(non_zero_inliers, 1));
%  
%  % H * source coordinates = destination coordinates
%  l = h_mat_final * [x_cord_s_test; y_cord_s_test; ones(size(x_cord_s_test))];
%  x_cord_dash = l(1, :)./l(3, :);
%  y_cord_dash = l(2, :)./l(3, :);
%  
%  test_soource = [x_cord_s_test; y_cord_s_test];
%  test_destination = [x_cord_d_test; y_cord_d_test];
%  test_dash = [x_cord_dash; y_cord_dash];
  
  %test_soource(:, 1:5)
  %test_destination(:, 1:5)
  %test_dash(:, 1:5)


%figure;
%imshow(transformedImage);


  %error_term = error_term + sqrt((x_cord_d - x_cord_dash)^2 + (y_cord_d - y_cord_dash)^2);
  
%  if(sqrt((x_cord_d - x_cord_dash)^2 + (y_cord_d - y_cord_dash)^2) < k)
%    k = sqrt((x_cord_d - x_cord_dash)^2 + (y_cord_d - y_cord_dash)^2);
%    index_min = j;
%    coord_s(counter, :) = [x_cord_s y_cord_s];
%    coord_d(counter, :) = [x_cord_d y_cord_d];
%    coord_dash(counter, :) = [x_cord_dash y_cord_dash];
%    counter = counter + 1;
%  endif

%coord_s
%coord_d
%coord_dash
%inliers_ransac'
%h_ransac

%Homography
%  oneColumn = ones(p_to_c,1);
%  zeroColumns = zeros(p_to_c,3);
%  constantColumn = [x_prime; y_prime];
%  topH = [x_prime, y_prime, oneColumn, zeroColumns, (-1.*x.*x_prime), (-1.*y.*x_prime)];
%  bottomH = [zeroColumns, x, y, oneColumn, (-1.*x.*y_prime), (-1.*y.*y_prime)];
%  matrixA = [topH; bottomH];
%
%  hVector = matrixA\constantColumn;
%  
%  hommographyMatrixTranspose = reshape([hVector;1],3,3);
%  
%  h_mat = hommographyMatrixTranspose';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imshow(i1_gray, []);
%figure;
%imshow(i2_gray, []);

%perm_i1 = randperm(size(f_i1, 2));
%sel_i1 = perm_i1(1:50);
%perm_i2 = randperm(size(f_i2, 2));
%sel_i2 = perm_i2(1:50);

%imshow(i1);
%i1_h1 = vl_plotsiftdescriptor(d_i1(:, sel_i1), f_i1(:, sel_i1));
%set(i1_h1, 'color', 'k', 'linewidth', 3);