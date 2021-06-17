function SOS_single_hex_planar_projection(Vfit, k)
    yalmip('clear')
    if nargin == 0
        Vfit = rand_hex(1,.3);
        k=5;
    end
    
    % recenter and get bounding box
%     Vfit = Vfit - mean(Vfit);
%     N = max(vecnorm(Vfit,2,2)*2)^2; % bounding radius for compactness
    
    plot_hex_surface(Vfit); colorbar; title('start hex'); rotate3d on;
    for i=1:8
        text(Vfit(i,1),Vfit(i,2),Vfit(i,3),['<-- ' num2str(i)])
    end

    %% build problem
    sdpvar bound;
    % polys == 0 defines V's making a planar hex mesh.
    [polys, V, Vf] = buildpolys();

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
    [s2, s2_c] = polynomial(vars, k);
    [s3, s3_c] = polynomial(vars, k);
    [s4, s4_c] = polynomial(vars, k);
    [s5, s5_c] = polynomial(vars, k);
    [s6, s6_c] = polynomial(vars, k);
    %[s7, s7_c] = polynomial(vars, k);

    mainobj = (V(:)-Vfit(:))'*(V(:)-Vfit(:)) + .1; % strictly positive
    %F_sos = [sos(s7)];
%     F_positiv = sos(mainobj - bound - [s1,s2,s3,s4,s5,s6,s7]*domain);
    F_positiv = sos(mainobj - bound - [s1,s2,s3,s4,s5,s6]*polys);

%     [F, obj, m] = sosmodel([F_sos, F_positiv], -bound, [], [s1_c; s2_c; s3_c; s4_c; s5_c; s6_c; s7_c; bound]);
    [F, obj, m] = sosmodel([F_positiv], -bound, [], [s1_c; s2_c; s3_c; s4_c; s5_c; s6_c; bound]);
    diagnostics = optimize(F, obj, []);

    %% extract result
	J = value(bound);

	moment_mat = dual(F(2));
	moment_mat = moment_mat/moment_mat(1,1);

	eig_vals = sort(eig(moment_mat));
	second_eig = eig_vals(end-1);

	point = moment_mat(2:25, 1);
	Vargmin = reshape(point(end:-1:1),3,8)';
    diff = norm(Vargmin - Vfit,'fro');
    plot_hex_surface(Vargmin); colorbar; title(['out hex: ' num2str(diff)]); rotate3d on;
end

function [polys, V, Vf] = buildpolys()
    sdpvar V11 V12 V13 ...
                 V21 V22 V23 ...
                 V31 V32 V33 ...
                 V41 V42 V43 ...
                 V51 V52 V53 ...
                 V61 V62 V63 ...
                 V71 V72 V73 ...
                 V81 V82 V83 real;
    V = [V11 V12 V13 ;...
        V21 V22 V23 ;...
        V31 V32 V33 ;...
        V41 V42 V43 ;...
        V51 V52 V53 ;...
        V61 V62 V63 ;...
        V71 V72 V73 ;...
        V81 V82 V83];
    faceIndexes = [1 2 3 4; 5 6 7 8; 1 5 6 2; 2 6 7 3; 3 7 8 4; 4 8 5 1];
    for i=1:6
        V1234 = V(faceIndexes(i,:),:);
        Vf(i,:) = V1234(:)';
        p1 = V1234(1,:);
        p2 = V1234(2,:);
        p3 = V1234(3,:);
        p4 = V1234(4,:);
        polys(i,:) = dot(cross(p2-p1,p4-p1),p3-p1);
    end
end


