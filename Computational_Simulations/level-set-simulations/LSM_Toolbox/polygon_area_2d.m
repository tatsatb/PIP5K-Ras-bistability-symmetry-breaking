function area = polygon_area_2d (v)

%% POLYGON_AREA_2D computes the area of a polygon in 2D.
%
%  Discussion:
%
%    AREA = 1/2 * abs ( sum ( 1 <= I <= N ) X(I) * ( Y(I+1) - Y(I-1) ) )
%    where Y(0) should be replaced by Y(N), and Y(N+1) by Y(1).
%
%  Modified:
%
%    30 January 2005
%
%  Modified by: Brett Kutscher
%  Original Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real V(2,N), the vertices.
%
%    Output, real area of the polygon.
%    Output is negative if the contour is listed clockwise
%    Output is positive if the contour is listed counter-clockwise
%
  area = 0.0;
  n = size(v,2);
  
  for i = 1 : n

    im1 = i_wrap ( i-1, 1, n );
    ip1 = i_wrap ( i+1, 1, n );

    area = area + v(1,i) * ( v(2,ip1) - v(2,im1) );

  end

  area = 0.5 * area;
