function [ order ] = graphtopoorder( dg )
%Sort a directed acyclic graph from a sparse input matrix
  %  obj = biograph(dg);
   % order_cells = sort3(obj);
   % order = cell2mat(order_cells);
   G = digraph(dg);
   order = toposort(G);
end
