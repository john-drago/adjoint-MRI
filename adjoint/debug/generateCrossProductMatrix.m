function V = generateCrossProductMatrix( v )
% function will generate cross product matrix for vector v, so that when
% applied to another vector w, i.e., Vw = v x w.

V = [    0, -v(3),  v(2);...
      v(3),     0, -v(1);...
     -v(2),  v(1),     0];

end