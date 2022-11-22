% Routine implementing a simple preferential attachment (B-A) model for network growth
% The probability that a new vertex attaches to a given old vertex is proportional to the (total) vertex degree
% Vertices arrive one at a time
% INPUTs: n - final (desired) number of vertices, m - # edges to attach at every step
% OUTPUTs: edge list, [number of edges x 3]
% NOTE: Assume undirected simple graph
% Source: "The Structure and Function of Complex Networks", M.E.J. Newman;  "Emergence of Scaling in Random Networks" B-A.
% GB, March 18, 2006

function [mat_final, el_final, Table] = preferential_attachment_ms(n,m)
%figure,imagesc(zeros(59)); %ms

vertices = 2; % default: 2
if not(vertices<=n); fprintf('the number of final nodes is smaller than the initial\n');  return; end
el=[1 2 1; 2 1 1];  % start with an edge

while vertices < n
    vertices=vertices+1;  % add new vertex

    if m>=vertices
      for node=1:vertices-1
        el = [el; node vertices 1];
        el = [el; vertices node 1];
      end
      continue
    end
    
    deg=[];               % compute nodal degrees for this iteration
    for v=1:vertices; 
        deg=[deg; v numel(find(el(:,1)==v))];
    end
    deg=sortrows(deg);
        
    % add m edges
    r = randsample(deg(:,1),m,'true',deg(:,2)/max(deg(:,2)));
    while not(length(unique(r))==length(r))
      r = randsample(deg(:,1),m,'true',deg(:,2)/max(deg(:,2)));
    end

    for node=1:length(r)
      el = [el; r(node) vertices 1];
      el = [el; vertices r(node) 1];  
    end
    
    %% Mini GTA (MS) continued
    mat = (edgeL2adj(el)); 
    %num_nodes=vertices 
    Table.deg_growth(:,vertices) = [sum(edgeL2adj(el)) zeros(1,n-vertices) ]';
    Table.net(vertices)           = {edgeL2adj(el)};
    Table.edges(vertices,:)      = sum(sum(edgeL2adj(el)))./(n*(n-1));
end
mat_final = mat;
el_final = el;

%% Prep plot (get fr coordinates)
%[coordinates, disps] = fr_graph_viz(mat_final);
 
%% Plot dev
figure; 
% subplot(1,4,1),imagesc(Table.deg_growth);
% subplot(1,4,2),plot(Table.edges);
% subplot(1,4,3),imagesc(mat);
% subplot(1,4,4),hist(reshape(sum(mat),1,[]));
% figure,gplot(mat_final,coordinates); hold on; plot(coordinates(1:n,1),coordinates(1:n,2),'r.'); %box off; axis square;
