function [coords_prime, coords, T_prime, T] = compute_dlt (x_prime, y_prime, x, y)

  x_prime_t = sum(x_prime)/numel(x_prime);
  y_prime_t = sum(y_prime)/numel(y_prime);

  x_t = sum(x)/numel(x);
  y_t = sum(y)/numel(y);

  x_prime_prime = x_prime - x_prime_t*ones(size(x_prime));
  y_prime_prime = y_prime - y_prime_t*ones(size(y_prime));

  x_s_prime = x - x_t*ones(size(x));
  y_s_prime = y - y_t*ones(size(y));

  scale_prime = (numel(x_prime) * sqrt(2))./sum((((x_prime_prime.*x_prime_prime) + (y_prime_prime.*y_prime_prime)).^0.5));
  scale_no_prime = (numel(x) * sqrt(2))./sum((((x_s_prime.*x_s_prime) + (y_s_prime.*y_s_prime)).^0.5));

  T_prime = [scale_prime 0 scale_prime*(-1.*x_prime_t); 0 scale_prime scale_prime*(-1.*y_prime_t); 0 0 1];
  T = [scale_no_prime 0 scale_no_prime*(-1.*x_t); 0 scale_no_prime scale_no_prime*(-1.*y_t); 0 0 1];
  
  prime_cords = [x_prime'; y_prime'; ones(1, numel(x_prime))];
  no_prime_cords = [x'; y'; ones(1, numel(x))];
  
  coords_prime = T_prime * prime_cords;
  coords = T * no_prime_cords;
  
endfunction