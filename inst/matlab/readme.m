% example of running hcp:

% (1) load data

load example_data

% input data: Y (matrix of gene expression, nxg where n is number of subjects and g is number of genes);
%             F (matrix of known covariates, nxd where d is number of covariates)


% (2) take log of read counts

Y = log2(2+Y);

% (3) standardize the data

pkg load statistics
source standardize.m

Yn = standardize(Y')';
Fn = standardize(F')';

Yn(1:5, 1:5)
Fn(1:5)

% (4) set the model parameters

k = 10
lambda = 20;
lambda2 = 1;
lambda3 = 1;
iter = 100;

% (5) run HCP
source hidden_covariates_model.m

[Z,B,U] = hidden_covariates_model(Fn,Yn,k,lambda,lambda2,lambda3,iter);

Z(1:5, 1:5)
B(1:5, 1:5)
U(1:5)

% (6) get the residual data:

Res = Yn-Z*B;
dlmwrite("Res.txt", Res, "\t")

