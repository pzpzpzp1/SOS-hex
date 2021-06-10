%% check that cube and hollow sphere constraints are archimedean. Success seems bound=4 is good enough.
clear all; close all;
k=4;
sdpvar u v w bound v1 v2 v3;

% set up model
vars1 = [u v w];
vars2 = [v1 v2 v3];
vars = [vars1 vars2];
vnorm_g = 1 - (v1^2+v2^2+v3^2); % |v|^2 <= 1

[s1, s1_c] = polynomial(vars, k);
[s2, s2_c] = polynomial(vars, k);
[s3, s3_c] = polynomial(vars, k);
[s4, s4_c] = polynomial(vars, k);
[s5, s5_c] = polynomial(vars, k);
[s6, s6_c] = polynomial(vars, k);
[s7, s7_c] = polynomial(vars, k);

domain = [u; 1 - u; v; 1 - v; w; 1 - w; ... % u in 0-1 cube
    vnorm_g;]; % hypercube

F_sos = [sos(s1), sos(s2), sos(s3), sos(s4), sos(s5), sos(s6)]; 

vobj = bound - (u^2 + v^2 + w^2 + v1^2 + v2^2 + v3^2);
F_positiv = sos(vobj - [s1, s2, s3, s4, s5, s6, s7]*domain);

[F, obj, m] = sosmodel([F_sos, F_positiv], bound, [], [s1_c; s2_c; s3_c; s4_c; s5_c; s6_c; bound;...
s7_c;]);    

diagnostics = optimize(F, obj, []);
diagnostics.problem
value(bound)



