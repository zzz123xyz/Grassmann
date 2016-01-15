% Projection metric distance measurement on Grassmannian manifolds
% Written by Arnold Wiliem *2011*
% manifold_toolbox
function [cost] = Gm_dist_projection(X,Y)

[U,D,V] = svd(X' * Y);
cost = sqrt(length(diag(D)) - sum(diag(D).^2));