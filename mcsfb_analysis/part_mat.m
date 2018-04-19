function [ row_part_ids ] = part_mat( A , col_part_ids , param )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    param=struct;
end

dim=size(A);
N=dim(1);
if N~=dim(2)
    error('A should be a square matrix');
end

parts=unique(col_part_ids);
num_part=length(parts);
row_part_ids=zeros(N,1);
completed_columns=zeros(N,1);

for i=1:(num_part-1)
   unassigned_columns=(completed_columns==0);
   E = A(row_part_ids==0,unassigned_columns);   

   sel_col= ( col_part_ids == parts(i));
   rem_col= col_part_ids(unassigned_columns);
   col1_inds =find(rem_col==parts(i));
   [ p1 , p2 ] = part_mat_2_by_2( E , col1_inds,param);

   unassigned_rows=find(row_part_ids==0);
   row_part_ids(unassigned_rows(p1))=i;
   
   completed_columns(sel_col)=1;

end
row_part_ids(unassigned_rows(p2))=num_part;

end

