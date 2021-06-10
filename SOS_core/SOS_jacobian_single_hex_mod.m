%{
if norm([v1 v2 v3])==0, then cube is not degenerate. means the only element of the nullspace is 0. 

note: This turns out to be a higher degree because the nullspace is a constraint instead of objective.

Computes the minimum eigenvalue of (gradX^T gradX) 
	V			the vertices of the hex
	k			optional, the degree of the SOS relaxation (default 4)
returns:
	--- tbd
%}
function [J, argmin, second_eig, diagnostics] = SOS_jacobian_single_hex_mod(V, k)
    if nargin == 0;
        V = rand_hex(1,.3);
    end
    
	if nargin < 2; k = 4; end

	sdpvar u v w bound v1 v2 v3;
    % syms u v w real;

	% set up model
	vars1 = [u v w];
    vars2 = [v1 v2 v3];
    vars = [u, v, w, v1, v2, v3];
%     vars = [u, v, w];

	Dphi = symbolic_trilinear_map_jacobian(V, [u,v,w]);
    Dphiv = Dphi*vars(4:6)'; % Dphi*v == 0
    vnorm_g = 1 - (v1^2+v2^2+v3^2); % |v|^2 <= 1
    d = randn(3,1); d = d/norm(d);
    vobj = dot([v1 v2 v3], d); % maximize v1
%     J = dot(cross(Dphi(:,1),Dphi(:,2)),Dphi(:,3));
%     J = symbolic_jacobian_det(V, [u, v, w]);

	[s1, s1_c] = polynomial(vars, k);
	[s2, s2_c] = polynomial(vars, k);
	[s3, s3_c] = polynomial(vars, k);
	[s4, s4_c] = polynomial(vars, k);
	[s5, s5_c] = polynomial(vars, k);
	[s6, s6_c] = polynomial(vars, k);
    [s7, s7_c] = polynomial(vars, k);
% 	[s8, s8_c] = polynomial(vars2, k);
% 	[s9, s9_c] = polynomial(vars2, k);
% 	[s10, s10_c] = polynomial(vars2, k);
% 	[s11, s11_c] = polynomial(vars2, k);
% 	[s12, s12_c] = polynomial(vars2, k);
    [r1, r1_c] = polynomial(vars, k);
    [r2, r2_c] = polynomial(vars, k);
    [r3, r3_c] = polynomial(vars, k);

	domain = [u; 1 - u; v; 1 - v; w; 1 - w; ... % u in 0-1 cube
        vnorm_g;];
%            1-v1; 1+v1; 1-v2; 1+v2; 1-v3; 1+v3;]; % v1,v2,v3 in [-1,1] cube
%     domain = [u; 1 - u; v; 1 - v; w; 1 - w; ] % u in 0-1 cube
    
	F_sos = [sos(s1), sos(s2), sos(s3), sos(s4), sos(s5), sos(s6),...
       sos(s7)]; %, sos(s8), sos(s9), sos(s10), sos(s11), sos(s12)];
%     F_sos = [sos(s1), sos(s2), sos(s3), sos(s4), sos(s5), sos(s6)];
% 	F_positiv = sos(vobj - bound - [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12]*domain - [r1,r2,r3]*Dphiv);
    F_positiv = sos(vobj - bound - [s1, s2, s3, s4, s5, s6, s7]*domain - [r1,r2,r3]*Dphiv);
%     F_positiv = sos(J - bound - [s1, s2, s3, s4, s5, s6]*domain);

	[F, obj, m] = sosmodel([F_sos, F_positiv], -bound, [], [s1_c; s2_c; s3_c; s4_c; s5_c; s6_c; bound;...
    s7_c; r1_c; r2_c; r3_c;]);    
    % s7_c; s8_c; s9_c; s10_c; s11_c; s12_c; r1_c; r2_c; r3_c;]);
%     [F, obj, m] = sosmodel([F_sos, F_positiv], -bound, [], [s1_c; s2_c; s3_c; s4_c; s5_c; s6_c; bound;]);
	
	% solve the problem
	diagnostics = optimize(F, obj, []);

	% extract the answers
	J = value(bound);

	moment_mat = dual(F(8));
	moment_mat = moment_mat/moment_mat(1,1);

	eig_vals = eig(moment_mat);
	second_eig = eig_vals(2);

	point = moment_mat(2:4, 1);
	argmin = point(end:-1:1);
end