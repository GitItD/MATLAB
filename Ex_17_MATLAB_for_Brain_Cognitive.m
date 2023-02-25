%% 17.6.1
%{
Why does the following line of code not produce a 2-by-1,000 matrix?
Add one character to make it produce a 2-by-1,000 matrix, and then
add one more character (without changing the previous character) to
produce a 1,000-by-2 matrix.
x = [ 1*randn(1,1000) .5*randn(1,1000) ];
%}
x = [ 1*randn(2, 500), .5*randn(2, 500) ];
x2 = [ 1*randn(2, 500), .5*randn(2, 500) ]';
whos
%% 17.6.2
%{
Relative eigenvalue magnitudes encode the proportion of variance
explained. Divide all eigenvalues by their sum and multiply by 100 to
transform the eigenvalues to percent variance explained. How much
variance is explained by the first 50 eigenfaces? Write code to compute
how many eigenfaces you would need to explain 73% (or any other
percent) of the variance.
%}
A = randi(100, 100, 100);
[eigvec, eigval] = eig(A);
eigval =  bsxfun(@rdivide, diag(eigval), sum(diag(eigval)));
disp(['Eigenvalue multiplied by 100 is ' (num2str(sum(eigval*100)))])
disp(['Variance explained by first 50 eigenfaces is equal to ' (num2str(sum(eigval(end-50:end))))])
disp(['Variance explained by first 73% of eigenfaces is equal to ' (num2str(sum(eigval(end-73:end))))])

% ChatGPT
var_exp = sum(eigval(1:50)) / sum(eigval);
fprintf('Variance explained by first 50 eigenfaces: %.2f%%\n', var_exp*100);
target_var = 0.73; % or any other desired proportion of variance
prop_var = 0;
num_eigenfaces = 0;
while prop_var < target_var
    num_eigenfaces = num_eigenfaces + 1;
    prop_var = prop_var + eigval(num_eigenfaces) / sum(eigval);
end
fprintf('Number of eigenfaces to explain %.2f%% of variance: %d\n', target_var*100, num_eigenfaces);

%% 17.6.5
%{
Build on figure 17.2 by repeating with 3D data. Youfll need to use
plot3 instead of plot, and you might need to do a bit of research to
figure out how to construct a 3D rotation matrix. When plotting, use
the view function to make sure the data and eigenvectors are easily
visible.
%}
% data
x = [ 1*randn(1000,1), .4*randn(1000,1), .6*randn(1000,1) ];

% Define the rotation angles in degrees
theta_x = 30; % Rotation angle around X axis
theta_y = 45; % Rotation angle around Y axis
theta_z = 60; % Rotation angle around Z axis

% Convert the angles to radians
theta_x = deg2rad(theta_x);
theta_y = deg2rad(theta_y);
theta_z = deg2rad(theta_z);

% Create the rotation matrices around each axis
Rx = [1 0 0; 0 cos(theta_x) -sin(theta_x); 0 sin(theta_x) cos(theta_x)];
Ry = [cos(theta_y) 0 sin(theta_y); 0 1 0; -sin(theta_y) 0 cos(theta_y)];
Rz = [cos(theta_z) -sin(theta_z) 0; sin(theta_z) cos(theta_z) 0; 0 0 1];

% Compute the 3D rotation matrix
R = Rz*Ry*Rx;

% rotate data
y = x*R;
z = y*R;

% PCA of x (original data)
x = bsxfun(@minus, x, mean(x,1));
covmat = (x'*x) / (length(x) - 1);
[evecsX, evalsX] = eig(covmat);

% PCA of y (correlated data)
y = bsxfun(@minus,y,mean(y,1));
covmat = (y'*y) / (length(y)-1);
[evecsY,evalsY] = eig(covmat);

%PCA of z (correlated data)
z = bsxfun(@minus, z, mean(z, 1));
covmat = (z'*z)/(length(z) - 1);
[evecsZ, evalsZ] = eig(covmat);

figure(2), clf
% plot original data
subplot(121)
plot3(x(:,1), y(:,2), z(:,3), 'm.','markersize', 5)
set(gca,'xlim', [-5 5], 'ylim', [-5 5], 'zlim', [-5 5])
hold on
plot(evalsX(1, 1)*[0, evecsX(1,1)], evalsX(1, 1)*[0, evecsX(2, 1)],'k','linew',4)
plot(evalsX(2, 2)*[0, evecsX(1,2)], evalsX(2, 2)*[0, evecsX(2, 2)],'k','linew',4)
plot(evalsY(1,1)*[0 evecsY(1,1)], evalsY(1,1)*[0 evecsY(2,1)], 'k', 'linew', 4)
plot(evalsY(2,2)*[0 evecsY(1,2)], evalsY(2,2)*[0 evecsY(2,2)], 'k', 'linew', 4)
plot(evalsZ(1,1)*[0 evecsY(1,1)], evalsZ(1,1)*[0 evecsZ(2,1)], 'k', 'linew', 4)
plot(evalsZ(2,2)*[0 evecsZ(1,2)], evalsZ(2,2)*[0 evecsZ(2,2)], 'k', 'linew', 4)
xlabel('x-axis'), ylabel('y-axis'), zlabel('z-axis')
axis square


subplot(122)
plot3(y(:,1), y(:,2), z(:,3), 'm.', 'markersize',5)
set(gca,'xlim', [-5 5], 'ylim', [-5 5], 'zlim', [-5 5])
hold on
plot(evalsY(1,1)*[0 evecsY(1,1)], evalsY(1,1)*[0 evecsY(2,1)], 'k', 'linew', 4)
plot(evalsY(2,2)*[0 evecsY(1,2)], evalsY(2,2)*[0 evecsY(2,2)], 'k', 'linew', 4)
plot(evalsZ(1,1)*[0 evecsY(1,1)], evalsZ(1,1)*[0 evecsZ(2,1)], 'k', 'linew', 4)
plot(evalsZ(2,2)*[0 evecsZ(1,2)], evalsZ(2,2)*[0 evecsZ(2,2)], 'k', 'linew', 4)
xlabel('x-axis'), ylabel('y-axis'), zlabel('z-axis')
axis square

% compute component scores
pc1 = y*evecsY(:,1);
pc2 = y*evecsY(:,2);
pc3 = z*evecsZ(:,3);

figure(3), clf
plot3(pc3, pc2, pc1, 'm.')
set(gca,'xlim',[-5 5],'ylim',[-5 5], 'zlim', [-5 5])
xlabel('PC1 axis'), ylabel('PC2 axis'), zlabel('PC3 axis')
hold on
plot([0 pc1(1,1)], [0 pc1(2,1)], 'k', 'linew', 4)
plot([0 pc2(1,2)], [0 pc2(2,2)], 'k', 'linew', 4)
plot([0 pc3(1,1)], [0 pc3(2,1)], 'k', 'linew', 4)
axis square

%% 17.6.6
%{
I claimed that MATLAB returns eigenvectors with unit length. Test this
in MATLAB by writing a script that will generate random square symmetric
matrices of any size, and compute the lengths of all eigenvectors.
If you get complex eigenvectors, youfll need to compute their
magnitude. Try this again with square nonsymmetric matrices.
%}

%% 17.6.9
%{
We can define the reconstruction error of each eigenface as the sum of
squared errors of the difference between the PCA-reconstructed face
versus the original face. Compute this for varying numbers of components
(1 to 200) used to reconstruct each of 500 faces. Plot the reconstruction
error as a function of the number of components, averaged
over all 500 faces. Show the variance over faces using a patch.
%}

%% 17.6.11
%{
A PCA of a covariance matrix can be used to create multivariate random
correlated data. Write code to implement this based on the following
description. Then test your code by computing the correlation
matrix. Create an N-by-N covariance matrix (N is the number of channels)
that contains only positive values. Then take its eigendecomposition,
then multiply the eigenvector matrix by the square root of the
eigenvalue matrix (vectors on the left, values on the right). Finally,
right-multiply the transpose of that matrix by an K-by-N matrix of randomly
generated numbers (K time points).
%}

%% 17.6.14
%{
Redo figure 17.2 using ICA instead of PCA. Make a figure that shows
both methods, similar to figure 17.7. You can see that PCA rotated the
axes to show the directions of maximum variance, whereas ICA
unmixed the data into uncoupled sources. Neither approach is
gbetterh; they are different methods, and in some cases they can be
used to obtain comparable outcomes.
%}

%% 17.6.15
%{
One of the many functions of PCA is to variance-normalize a matrix,
which is also called sphering a matrix. This means reconstructing the
data to have equal variance in all principal directions. The formula is
y? = ySL?1/2, where y is the original data, S is the matrix of eigenvectors,
ƒ© is the matrix of eigenvalues, and ƒ©?1/2 indicates the square root of the
matrix inverse of ƒ©. Implement this formula for the data used in figures
17.2 and 17.7. Then plot the new data ( ?y) on top of the original data
(y). What is the effect of sphering the data? What are the eigenvalues
and eigenvectors of the sphered data?
%}