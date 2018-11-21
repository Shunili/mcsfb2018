function [ p1 , p2 ] = part_mat_2_by_2( E , col1_inds , param )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% add a rank check on E? or better yet a condition number check on the
% output matrix piece

dim=size(E);
N=dim(1);
if N~=dim(2)
    error('E should be a square matrix');
end

K=length(col1_inds);
if K>N
    error('Too many column indices specified');
end

max_its=N;

if nargin < 3
    param=struct;
end

if ~isfield(param,'init_method')
    init_method='greedy_init';
else
    init_method=param.init_method;
end

if isfield(param,'row_guess1')
    row_guess1=unique(param.row_guess1);
else
    switch init_method
        case 'greedy_init' % basically an implementation of the Steinitz exchange lemma
% else
%     if isfield(param,'greedy_init')
%         greedy_init=param.greedy_init;
%     else
%         greedy_init=1;
%     end
%     if greedy_init % basically an implementation of the Steinitz exchange lemma
%tic
            row_guess1=zeros(1,K);
            B=eye(N);
            chosen=zeros(N,1); 
            for j=1:K % speed up by not refactoring the entire B matrix each loop?
                if mod(j,100)==0
                    iter=j
                end
                x=B\E(:,col1_inds(j));
                x=(ones(N,1)-chosen).*x; %set the already chosen vertices to zeros
                [~,b]=max(abs(x)); %b is the index
                chosen(b)=1;
                row_guess1(j)=b;
                B(:,b)=E(:,col1_inds(j));
            end
%initial_guess_time=toc
        case 'greedy_init2'
            row_guess1=zeros(1,K);
            chosen=zeros(N,1);
            for j=1:K
                max_sigma_min=0;
                selection=-1;
                for i=1:N
                    if ~chosen(i)
                        new_set=[row_guess1(1:j-1),i];
                        sv= svds(E(new_set,col1_inds),1,'smallest');
                        if sv > max_sigma_min
                            selection=i;
                        end
                    end      
                end
                chosen(selection)=1;
                row_guess1(j)=selection;
            end
        case 'rand_init'
            norm_Uk=sum(E.^2, 2);
            weights=norm_Uk/sum(norm_Uk);
            row_guess1 = datasample(1:N, K, 'Replace', false, 'Weights', weights); 
        case 'hybrid_init'    
            if N<3000
                param.init_method='greedy_init';
                [p1,p2]=part_mat_2_by_2(E,col1_inds,param);
                return;
            else
                 param.init_method='rand_init';
                [p1,p2]=part_mat_2_by_2(E,col1_inds,param);
                return;
            end
                
        otherwise
%    elseif 
%    else
            row_guess1=1:K;
    end
end
%tic
if (sum(row_guess1==0)>0)
    error('not all row guesses assigned');
end
if length(unique(row_guess1))~=K
    error('not enough entries in row guess')
end

row_guess2=setdiff(1:N,row_guess1);
idx2=[row_guess2,row_guess1];

if ~strcmp(init_method,'greedy_init') && ~strcmp(init_method,'greedy_init2')
    idx1=[row_guess1,row_guess2];
    E1=E(:,col1_inds);
    E1=E1';
    [R1,piv1] = rref(E1(:,idx1)); 
    %piv1=sort(idx1(piv1),'ascend');
    piv1=idx1(piv1);
else
    piv1=row_guess1;
end

E2=E(:,setdiff(1:N,col1_inds));
E2=E2';
[R2,piv2] = rref(E2(:,idx2));
%piv2=sort(idx2(piv2),'ascend'); 
piv2=idx2(piv2);

% % Check full rank
% rank_check1=K-rank(E1(:,piv1))
% rank_check2=N-K-rank(E2(:,piv2))

%assigned=sort([piv1,piv2],'ascend');
assigned=[piv1,piv2];
unassigned=setdiff(1:N,assigned);

%middle_step_time=toc
%tic
if ~isempty(unassigned)
    idx1inv=zeros(N,1);
    idx1inv(idx1)=1:N;
    R1=R1(:,idx1inv);

    idx2inv=zeros(N,1);
    idx2inv(idx2)=1:N;
    R2=R2(:,idx2inv); 
end

while ~isempty(unassigned)        
    dup=intersect(piv1,piv2);
    to_check=unassigned(1); %pick y to start
    checked=[];
    nodes=to_check;
    edges=[];
    found=0;
    for j=1:max_its
        %chain_len=j
        new_nodes=[];
        new_edges=[];
        checked=[checked,to_check'];
        unchecked=setdiff(1:N,checked);
        for k=1:length(to_check)
            x1=(R1(:,to_check(k))~=0);
            y1=piv1(x1);
            x2=(R2(:,to_check(k))~=0);
            y2=piv2(x2);
            exch1=intersect(y1,dup); 
            exch2=intersect(y2,dup);
            %randomize which is checked first? put these two if
            %else into a single group?
            if ~isempty(exch1)
                exch1=exch1(randperm(length(exch1)));
                sel_ind=find(nodes(:,j)==to_check(k));
                sel_chain_nodes=[nodes(sel_ind,:),exch1];
                if isempty(edges)
                    sel_chain_edges=1;
                else
                    sel_chain_edges=[edges(sel_ind,:),1];
                end
                c1=find(sel_chain_edges==1);
                c2=find(sel_chain_edges==2);
                %piv1=sort([setdiff(piv1,sel_chain_nodes(c1+1)),sel_chain_nodes(c1)],'ascend'); 
                %piv2=sort([setdiff(piv2,sel_chain_nodes(c2+1)),sel_chain_nodes(c2)],'ascend'); 
                piv1=[setdiff(piv1,sel_chain_nodes(c1+1)),sel_chain_nodes(c1)]; 
                piv2=[setdiff(piv2,sel_chain_nodes(c2+1)),sel_chain_nodes(c2)];
                found=1;
                break;
            elseif ~isempty(exch2)
                exch2=exch2(randperm(length(exch2)));
                sel_ind=find(nodes(:,j)==to_check(k));
                sel_chain_nodes=[nodes(sel_ind,:),exch2];
                if isempty(edges)
                    sel_chain_edges=2;
                else
                    sel_chain_edges=[edges(sel_ind,:),2];
                end
                c1=find(sel_chain_edges==1);
                c2=find(sel_chain_edges==2);
                %piv1=sort([setdiff(piv1,sel_chain_nodes(c1+1)),sel_chain_nodes(c1)],'ascend'); 
                piv1=[setdiff(piv1,sel_chain_nodes(c1+1)),sel_chain_nodes(c1)];
                %piv2=sort([setdiff(piv2,sel_chain_nodes(c2+1)),sel_chain_nodes(c2)],'ascend'); 
                piv2=[setdiff(piv2,sel_chain_nodes(c2+1)),sel_chain_nodes(c2)]; 
                found=1;
                break;
            elseif ~isempty(setdiff(unchecked,dup))
                % randomize the order here?
                next1=intersect(y1,unchecked);
                checked=unique([checked,y1]);
                unchecked=setdiff(1:N,checked);
                next2=intersect(y2,unchecked);
                checked=unique([checked,y2]);
                unchecked=setdiff(1:N,checked);
                next=[next1,next2];
                if ~isempty(next)
                    old_ind=find(nodes(:,j)==to_check(k));
                    new_nodes=[new_nodes;[repmat(nodes(old_ind,:),length(next),1),next']];
                    if isempty(edges)
                        new_edges=[new_edges;[ones(length(next1),1);2*ones(length(next2),1)]];
                    else
                        new_edges=[new_edges;[repmat(edges(old_ind,:),length(next),1),[ones(length(next1),1);2*ones(length(next2),1)]]];
                    end
                end
            end
        end  
        if found
            break;
        else
            to_check=unique(new_nodes(:,j+1));
            to_check=to_check(randperm(length(to_check)));
            nodes=new_nodes;
            edges=new_edges;
        end
    end
    
    %assigned=sort([piv1,piv2],'ascend');
    assigned=[piv1,piv2];
    unassigned=setdiff(1:N,assigned);

    if isempty(unassigned) 
        break;
    end

    idx1=[piv1,setdiff(1:N,piv1)];
    [R1,piv1] = rref(E1(:,idx1));
    idx1inv=zeros(N,1);
    idx1inv(idx1)=1:N;
    R1=R1(:,idx1inv); 
    %piv1=sort(idx1(piv1),'ascend');
    piv1=idx1(piv1);
    
    idx2=[piv2,setdiff(1:N,piv2)];
    [R2,piv2] = rref(E2(:,idx2));
    idx2inv=zeros(N,1);
    idx2inv(idx2)=1:N;
    R2=R2(:,idx2inv); 
    %piv2=sort(idx2(piv2),'ascend'); 
    piv2=idx2(piv2);
end

%adjust_time=toc

p1=piv1;
p2=piv2;

% % double check full rank matrices
% end_check1=K-rank(E(piv1,1:K))
% end_check2=N-K-rank(E(piv2,K+1:N))
% 
% % double check partition
% proper_num=N-length(piv1)-length(piv2)
% leftover=setdiff(1:N,[piv1,piv2])


end

