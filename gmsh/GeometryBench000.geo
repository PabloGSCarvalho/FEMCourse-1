

IsquadQ = 1;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Definição de parâmtros

a = 5; 
b = 5;
h = 1;
L = 1;
Lf = 1;

n_bc = 2;
nx = 3;
ny = 3;
pr = 1;

// Coordenadas dos pontos

  //Domínio Omega
  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {L, 0, 0, 1e+22};
  Point(3) = {L, h, 0, 1e+22};
  Point(4) = {0, h, 0, 1e+22};

// Fronteiras

  //Domínio Omega  
  Line(1) = {1, 2};
  Line(2) = {2, 3};
  Line(3) = {3, 4};
  Line(4) = {4, 1};  



  Transfinite Line{1,3} = nx Using Progression pr;
  Transfinite Line{2,4} = ny Using Progression pr;
 // Transfinite Line{7} = nx Using Progression pr;


// Definição da superfície 

  Line Loop(5) = {1, 2, 3, 4};
//   Line Loop(6) = {7, 3, 4, 5};

  Plane Surface(6) = {5};
//  Plane Surface(7) = {6};

 // Line{7} in Surface{6};

  Transfinite Surface {6} = {1,2,3,4};
//  Transfinite Surface {7} = {5,6,3,4};

  If(IsquadQ)

  Recombine Surface {6};
//  Recombine Surface {7};

  EndIf

 // Physical Volume("internal") = {1};
 // Extrude {0, 0, 10} {
 //  //Surface{110};
 //  Layers{1};
 //  Recombine;
 // }

  Physical Surface("Omega") = {6};
  Physical Line("bottom") = {1};
  Physical Line("top") = {3};
  Physical Line("right") = {2};
  Physical Line("left") = {4};
//  Physical Point("PointLeft") = {5};
//  Physical Point("PointRight") = {6};
  
  //Physical Line("holes") = {holes[]};  
  //Physical Surface("interface") = {23};

  
  Coherence Mesh;


