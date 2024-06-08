clear;clc;

Dp = 25;
Dm = 25;
s = 1;

tmax = 5;
dt = 0.01;
xspan = -10:0.1:10;
yspan = -10:0.1:10;

[X,Y] = meshgrid(xspan,yspan);
t = 0;
p = rand(size(X));
m = rand(size(Y));
dpdtMatrix = zeros(size(p));
dmdtMatrix = zeros(size(m));

while t < tmax
    diffusionMatrixP = Dp * laplace(p);
    diffusionMatrixM = Dm * laplace(m);
    for ix = 1:length(xspan)
        for iy = 1:length(yspan)
            pos = [xspan(ix),yspan(iy)];
            pLoc = p(ix,iy);
            mLoc = m(ix,iy);

            dpdtMatrix(ix,iy) = diffusionMatrixP(ix,iy) + f(pLoc,mLoc,s,pos);
            dmdtMatrix(ix,iy) = diffusionMatrixM(ix,iy) + g(pLoc,mLoc,pos);
        end
    end
    p = max(p + dpdtMatrix * dt,0);
    m = max(m + dmdtMatrix * dt,0);
    t = t + dt;

    figure(1)
    surf(X,Y,p)
    title("P53 at t="+t)
    xlabel("x")
    ylabel("y")
    xlim([-10,10])
    ylim([-10,10])

    figure(2)
    surf(X,Y,m)
    title("Mdm2 at t="+t)
    xlabel("x")
    ylabel("y")
    xlim([-10,10])
    ylim([-10,10])
end


function out = laplace(A)
    kernel = [0,1,0;
             1,-4,1;
              0,1,0];
    out = conv2(A,kernel,"same");
end

function out = indicator(pos)
    x = pos(1);
    y = pos(2);
    r = sqrt(x^2+y^2);

    if r <= 1
        out = 1;
    else
        out = 0;
    end
end

function out = f(p,m,s,pos) %f and g are missing indicator function stuff
    k1 = 1;
    alph = 2;
    kcat = 1;
    km = 10;
    q = 10;

    I = indicator(pos);

    out = I * k1 * (1 + alph*s^2) / (1 + s^2) - kcat*p*m / (km + p) * (1-I) - q*p;
end

function out = g(p,m,pos)
    k2 = 1;
    alph = 2;
    q = 10;

    I = indicator(pos);

    out = I * k2 * (1 + alph*p^2) / (1 + p^2) - q*m;
end

% function u = eulerStep(u0,pos,odeHandle,dt)
%     dxdt = odeHandle(u0,pos);
%     u = u0 + dxdt*dt;
% end