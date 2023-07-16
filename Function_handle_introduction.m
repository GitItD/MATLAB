% Function handle introduction:

sin([0 pi/2 pi 3/2*pi 2*pi])
h = @sin
feval(h,[0 pi/2 pi 3/2*pi 2*pi])
% Feval help: https://www.mathworks.com/help/matlab/ref/feval.html
fplot(h,[0 2*pi])

%{
Comparing both results illustrates that passing the function handle
worked as expected.
You might wonder what the big deal is. It is arguably as easy?if not easier?to just
type the values directly into the sin function than to formally declare a function handle.
Of course, you would be right to be skeptical. However?at the very least?you will save
time typing when you use the same function over and over again?given that you use function
handles that are shorter than the function itself. Moreover, you can create more succinct
code, which is always a concern as your programs get longer and more intricate.
More importantly, there are functions that actually do useful stuff with function handles.
For example, fplot plots a given function over a specified range.
%}

AA = integral(h,0,pi)
BB = integral(h,0,2*pi)
CC = integral(h,0,pi/2)

%{ 
You'll need to turn your input into a function handle in order to use
integrals. Just try using >> DD = integral(T,0,pi)
%}

T = cos([0 pi/2 pi 3/2*pi 2*pi]);

%{
In addition, you can not only tag pre-existing MATLAB functions, but also declare your
own functions and tag them with a function handle, as follows:
%}
QQ = @(x) x.^5-9.*x^4+8.*x^3-2.*x.^2+x+500;
figure
fplot(QQ, [0 10])

%{
The function handle works when declaring an unknown variable, such as X,
but you wouldn't be able to assign it without the function handle. For
instance, try assigning: AA = X^2 + X
%}

% Ex 2.30:
woa = @(w) w.^3 + w.^2 + 100;
figure
fplot(woa, [1 5])
wow = integral(woa, 1, 5)
easy = @(ez) 2*ez
wow_easy = integral(easy, 1, 2) %#ok<*NOPTS>

%{ 
As you can see, function handles are mainly used to call other functions by
their names. Other functions that takes function handles include:
%}

I = @ones;
II = I(4)
III = I

XIN = @sin
CUZ = @cos
YAN = @tan
XIN(1) + CUZ(1) + YAN(0)