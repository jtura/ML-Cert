%% Exploring bounds
clear all
% number particles
n = 5;

% one/two body (always needed!)
one_two = ones(2,n);

% other marginals
rest = zeros(n-2, n);

% all
all = [one_two;rest];

% changing
all(3,1) = 1;

%% SDP

tic
[aa, E] = go1D(all);
toc