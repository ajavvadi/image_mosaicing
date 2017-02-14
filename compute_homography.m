## Copyright (C) 2016 Kalyan
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {@var{retval} =} compute_homography (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Kalyan <kalyan@kalyan-XPS-L501X>
## Created: 2016-05-08

function [h_mat] = compute_homography (x_prime_s, y_prime_s, x_s, y_s)

  if(numel(x_prime_s) < 4)
  
    display(1);
    h_mat = [1 0 0; 0 1 0; 0 0 1];
    
  else
      
      [coords_prime coords T_prime T] = compute_dlt(x_prime_s, y_prime_s, x_s, y_s);
    x_prime = (coords_prime(1, :)./coords_prime(3, :))';
    y_prime = (coords_prime(2, :)./coords_prime(3, :))';
    x = (coords(1, :)./coords(3, :))';
    y = (coords(2, :)./coords(3, :))';
    p_to_c = size(x_prime, 1);
    oneColumn = ones(p_to_c,1);
    zeroColumns = zeros(p_to_c,3);
    constantColumn = [x_prime; y_prime];
  
%    topH = [x_prime, y_prime, oneColumn, zeroColumns, (-1.*x.*x_prime), (-1.*y.*x_prime)];
%    bottomH = [zeroColumns, x, y, oneColumn, (-1.*x.*y_prime), (-1.*y.*y_prime)];
%    matrixA = [topH; bottomH];
%
%    hVector = matrixA\constantColumn;
%  
%    hommographyMatrixTranspose = reshape([hVector;1],3,3);
%  
%    h_mat = hommographyMatrixTranspose';
  
  topH = [x_prime, y_prime, oneColumn, zeroColumns, (-1.*x.*x_prime), (-1.*y.*x_prime), -1.*x_prime];
  bottomH = [zeroColumns, x, y, oneColumn, (-1.*x.*y_prime), (-1.*y.*y_prime), -1.*y_prime];
  
  matrixA = [topH; bottomH];
  
  [u d v] = svd(matrixA'*matrixA);
  h_vector = v(:, end);
  h_mat = reshape(h_vector, 3, 3)';
  
  h_mat = inv(T_prime)*h_mat*T;
  
  endif

endfunction
