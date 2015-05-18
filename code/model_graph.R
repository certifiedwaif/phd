library(DiagrammeR)
grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true]

  # several 'node' statements
  node [shape = box,
        fontname = Helvetica]
  A; B; C; D; E; F

  node [shape = circle,
        fixedsize = true,
        width = 0.9] // sets as circles
  1; 2; 3; 4; 5; 6; 7; 8

  # several 'edge' statements
  A->1; B->2; B->3; B->4; C->A
  1->D; E->A; 2->4; 1->5; 1->F
  E->6; 4->6; 5->7; 6->7; 3->8
}
")

grViz('
digraph zip {
  layout = dot
  node [shape = circle]
  y; x; r; beta ; rho; sigmabeta;sigmau

  edge [color = blue]
  y->x; y->r;
  r->rho
  x->beta; x->u
  beta->sigmabeta
  u->sigmau
}')