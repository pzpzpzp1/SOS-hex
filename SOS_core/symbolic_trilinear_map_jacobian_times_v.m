%{
The trilinear map for the hex with given vertices, evaluated at the given points
	V	vertices of hex, 8x3
	p	points to evaluate at, Nx3
    vec   vector to multiply jacobian by
returns:
	Dphi jacobian values of the trilinear map, 3x3xN
%}
function [coeffmat, inputs] = symbolic_trilinear_map_jacobian_times_v(V, p, vec)
    if nargin == 0
        syms p1 p2 p3 real;
        syms v1 v2 v3 real;
        syms V11 V12 V13 ...
             V21 V22 V23 ...
             V31 V32 V33 ...
             V41 V42 V43 ...
             V51 V52 V53 ...
             V61 V62 V63 ...
             V71 V72 V73 ...
             V81 V82 V83 real;
         p = [p1 p2 p3]; vec=[v1 v2 v3]';
         V = [V11 V12 V13 ;...
                V21 V22 V23 ;...
                V31 V32 V33 ;...
                V41 V42 V43 ;...
                V51 V52 V53 ;...
                V61 V62 V63 ;...
                V71 V72 V73 ;...
                V81 V82 V83];
    end
    
    inputs.V = V; inputs.p = p; inputs.vec = vec;
    
    u = reshape(p(:,1),1,1,[]);
    v = reshape(p(:,2),1,1,[]);
    w = reshape(p(:,3),1,1,[]);
    V11 = V(1,1); V12 = V(1,2); V13 = V(1,3);
    V21 = V(2,1); V22 = V(2,2); V23 = V(2,3);
    V31 = V(3,1); V32 = V(3,2); V33 = V(3,3);
    V41 = V(4,1); V42 = V(4,2); V43 = V(4,3);
    V51 = V(5,1); V52 = V(5,2); V53 = V(5,3);
    V61 = V(6,1); V62 = V(6,2); V63 = V(6,3);
    V71 = V(7,1); V72 = V(7,2); V73 = V(7,3);
    V81 = V(8,1); V82 = V(8,2); V83 = V(8,3);
    
	Dphi = [[V41.*v.*(w - 1) - V31.*v.*(w - 1) + V51.*w.*(v - 1) - V61.*w.*(v - 1) - V11.*(v - 1).*(w - 1) + V21.*(v - 1).*(w - 1) + V71.*v.*w - V81.*v.*w, V21.*u.*(w - 1) - V31.*u.*(w - 1) + V51.*w.*(u - 1) - V81.*w.*(u - 1) - V11.*(u - 1).*(w - 1) + V41.*(u - 1).*(w - 1) - V61.*u.*w + V71.*u.*w, V51.*(u - 1).*(v - 1) - V11.*(u - 1).*(v - 1) - V31.*u.*v + V71.*u.*v + V21.*u.*(v - 1) + V41.*v.*(u - 1) - V61.*u.*(v - 1) - V81.*v.*(u - 1)];...
        [V42.*v.*(w - 1) - V32.*v.*(w - 1) + V52.*w.*(v - 1) - V62.*w.*(v - 1) - V12.*(v - 1).*(w - 1) + V22.*(v - 1).*(w - 1) + V72.*v.*w - V82.*v.*w, V22.*u.*(w - 1) - V32.*u.*(w - 1) + V52.*w.*(u - 1) - V82.*w.*(u - 1) - V12.*(u - 1).*(w - 1) + V42.*(u - 1).*(w - 1) - V62.*u.*w + V72.*u.*w, V52.*(u - 1).*(v - 1) - V12.*(u - 1).*(v - 1) - V32.*u.*v + V72.*u.*v + V22.*u.*(v - 1) + V42.*v.*(u - 1) - V62.*u.*(v - 1) - V82.*v.*(u - 1)];...
        [V43.*v.*(w - 1) - V33.*v.*(w - 1) + V53.*w.*(v - 1) - V63.*w.*(v - 1) - V13.*(v - 1).*(w - 1) + V23.*(v - 1).*(w - 1) + V73.*v.*w - V83.*v.*w, V23.*u.*(w - 1) - V33.*u.*(w - 1) + V53.*w.*(u - 1) - V83.*w.*(u - 1) - V13.*(u - 1).*(w - 1) + V43.*(u - 1).*(w - 1) - V63.*u.*w + V73.*u.*w, V53.*(u - 1).*(v - 1) - V13.*(u - 1).*(v - 1) - V33.*u.*v + V73.*u.*v + V23.*u.*(v - 1) + V43.*v.*(u - 1) - V63.*u.*(v - 1) - V83.*v.*(u - 1)]];
    
    Dphivec = Dphi*vec;
    simplify(Dphivec);
%     for i=1:numel(V)
%         coeffmat(:,i) = diff(Dphivec,V(i));
%     end
    coeffmat = jacobian(Dphivec,V(:));
    
end

%{
syms u v w real;
syms V11 V12 V13 real;
syms V21 V22 V23 real;
syms V31 V32 V33 real;
syms V41 V42 V43 real;
syms V51 V52 V53 real;
syms V61 V62 V63 real;
syms V71 V72 V73 real;
syms V81 V82 V83 real;

symp = [u v w];
symV = [V11 V12 V13;...
V21 V22 V23 ;...
V31 V32 V33 ;...
V41 V42 V43 ;...
V51 V52 V53 ;...
V61 V62 V63 ;...
V71 V72 V73 ;...
V81 V82 V83 ;];
phi = symbolic_trilinear_map(V, p);
J11 = diff(phi(1),u);
J12 = diff(phi(1),v);
J13 = diff(phi(1),w);
J21 = diff(phi(2),u);
J22 = diff(phi(2),v);
J23 = diff(phi(2),w);
J31 = diff(phi(3),u);
J32 = diff(phi(3),v);
J33 = diff(phi(3),w);
Jsym = [J11 J12 J13; J21 J22 J23; J31 J32 J33;];
%}

