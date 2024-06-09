clear;clc;

Dp = 8172;
Dm = 8172;
s = 10;

tmax = 5/60;
dt = 0.00001;
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

    p(1,:) = zeros(size(p(1,:)));
    p(end,:) = zeros(size(p(end,:)));
    p(:,1) = zeros(size(p(:,1)));
    p(:,end) = zeros(size(p(end,:)));

    m(1,:) = zeros(size(p(1,:)));
    m(end,:) = zeros(size(p(end,:)));
    m(:,1) = zeros(size(p(:,1)));
    m(:,end) = zeros(size(p(end,:)));

    figure(1)
    surf(X,Y,p)
    zlabel("P (nM)")
    xlabel("x (\mu m)")
    ylabel("y (\mu m)")
    title("P at t="+3600*t+"s")
    xlim([-10,10])
    ylim([-10,10])

    figure(2)
    surf(X,Y,m)
    title("Mdm2 at t="+3600*t+"s")
    zlabel("M (nM)")
    xlabel("x (\mu m)")
    ylabel("y (\mu m)")
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

    out = r <= 2.5;
end

function out = f(p,m,s,pos) %f and g are missing indicator function stuff
    k1 = 1000;
    alph = 2;
    kcat = 11;
    km = 10^5;
    q = 0.1;

    I = indicator(pos);

    out = I * k1 * (1 + alph*s^2) / (1 + s^2) - kcat*p*m / (km + p) * (1-I) - q*p;
end

function out = g(p,m,pos)
    k2 = 0.47;
    alph = 2;
    q = 0.2;

    I = indicator(pos);

    out = I * k2 * (1 + alph*p^2) / (1 + p^2) - q*m;
end

% function u = eulerStep(u0,pos,odeHandle,dt)
%     dxdt = odeHandle(u0,pos);
%     u = u0 + dxdt*dt;
% end