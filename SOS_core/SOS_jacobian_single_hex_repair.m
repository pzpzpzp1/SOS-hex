%% init
clear all; close all;
Vfit = rand_hex(1,.5);

model = SOS_jacobian_batch_model(4);
minjacdet1 = model(Vfit); 
plot_hex_surface(Vfit); colorbar; title(sprintf('start hex: %f',minjacdet1)); rotate3d on;

%% build problem
[coeffmat, inputs] = symbolic_trilinear_map_jacobian_times_v;
V = inputs.V;
u = inputs.p(1); v = inputs.p(2); w = inputs.p(3);
v1 = inputs.vec(1); v2 = inputs.vec(2); v3 = inputs.vec(3);
VV = sdpvar(25, 25);
domain = [u; 1 - u; v; 1 - v; w; 1 - w; 1 - (v1^2+v2^2+v3^2);];
vars = [u v w v1 v2 v3];

epsilon = 2;
AA = coeffmat'*coeffmat;
mineigval = dot(AA(:),reshape(VV(2:end,2:end),[],1))-epsilon; 

k=4; 
[s1, s1_c] = polynomial(vars, k);
[s2, s2_c] = polynomial(vars, k);
[s3, s3_c] = polynomial(vars, k);
[s4, s4_c] = polynomial(vars, k);
[s5, s5_c] = polynomial(vars, k);
[s6, s6_c] = polynomial(vars, k);
[s7, s7_c] = polynomial(vars, k);

F_sos = [sos(s1), sos(s2), sos(s3), sos(s4), sos(s5), sos(s6)]
F_positiv = sos(mineigval - [s1,s2,s3,s4,s5,s6,s7]*domain);

actualobj = trace(VV) -1 - 2*dot(VV(2:end,1), Vfit(:));

[F, obj, m] = sosmodel([F_sos, F_positiv, VV>=0; VV(1)==1], actualobj, [], [s1_c; s2_c; s3_c; s4_c; s5_c; s6_c; s7_c;]);
diagnostics = optimize(F, obj, []);

%% extract result
VVval = value(VV);
eigs(VVval,5,'lm')'
Vnew = reshape(VVval(2:end,1),8,3);
minjacdet2 = model(Vnew); 
plot_hex_surface(Vnew); colorbar; title(sprintf('repaired: %f',minjacdet2)); rotate3d on;

%% compute bounds
vv = VVval(2:end,2:end);
lowerbound = trace(vv) - 2*dot(Vnew(:),Vfit(:)) + dot(Vfit(:),Vfit(:));
% NOTE: upperbound is only valid if rounding still satisfies the constraints, which it may not (demonstrated possible).
upperbound = norm(Vfit(:)-Vnew(:))^2; 

% result is not rank 1 but the eigenvalue ratio is ~100x. Rounding based on largest eigval fails to provide feasible solution (though it is non-degen for high enough epsilon)






