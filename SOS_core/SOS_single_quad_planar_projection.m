function SOS_single_quad_planar_projection(Vfit, k)
    yalmip('clear'); close all;
    if nargin == 0
        Vfit = rand_hex(1,.3); 
        Vfit = Vfit(1:4,:);
        k=8;
    end
    
    plot_hex_surface(repmat(Vfit,2,1)); colorbar; title('start hex'); rotate3d on;
    for i=1:4; text(Vfit(i,1),Vfit(i,2),Vfit(i,3),['<-- ' num2str(i)]); end

    %% build problem
    sdpvar bound;
    % polys == 0 defines V's making a planar hex mesh.
    [polys, V] = buildpolys();

%     domain = [polys; N - V(:)'*V(:);];

%     yalmip('clear')
%     [s1, s1_c] = polynomial(Vf(1,:), k);
%     [s2, s2_c] = polynomial(Vf(2,:), k);
%     [s3, s3_c] = polynomial(Vf(3,:), k);
%     [s4, s4_c] = polynomial(Vf(4,:), k);
%     [s5, s5_c] = polynomial(Vf(5,:), k);
%     [s6, s6_c] = polynomial(Vf(6,:), k);
    
    vars = V(:);
    [s1, s1_c] = polynomial(vars, k);
    %[s7, s7_c] = polynomial(vars, k);

    mainobj = (V(:)-Vfit(:))'*(V(:)-Vfit(:)); % strictly positive
    
    F_positiv = sos(mainobj - bound - s1*polys);

%     [F, obj, m] = sosmodel([F_sos, F_positiv], -bound, [], [s1_c; s2_c; s3_c; s4_c; s5_c; s6_c; s7_c; bound]);
    [F, obj, m] = sosmodel([F_positiv], -bound, [], [s1_c; bound]);
    diagnostics = optimize(F, obj, []);

    %% extract result
	J = value(bound);

	moment_mat = dual(F(2));
	moment_mat = moment_mat/moment_mat(1,1);

	eig_vals = sort(eig(moment_mat));
	second_eig = eig_vals(end-1);

	point = moment_mat(2:13, 1);
	Vargmin = reshape(point(end:-1:1),3,4)';
    diff = norm(Vargmin - Vfit,'fro');
    plot_hex_surface(repmat(Vargmin,2,1)); colorbar; title(['out hex: ' num2str(diff)]); rotate3d on;
    
    figure; plot(value(s1_c)); title('monomial coeffs');
end

function [polys, V] = buildpolys()
    sdpvar V11 V12 V13 ...
                 V21 V22 V23 ...
                 V31 V32 V33 ...
                 V41 V42 V43 ...
                 real;
    V = [V11 V12 V13 ;...
        V21 V22 V23 ;...
        V31 V32 V33 ;...
        V41 V42 V43 ];
    
    p1 = V(1,:);
    p2 = V(2,:);
    p3 = V(3,:);
    p4 = V(4,:);
    polys = dot(cross(p2-p1,p4-p1),p3-p1);
end


